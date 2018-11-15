#snippet header
# Author: Gholamali (Ali) rahnavard
# Email: gholamali.rahnavard@gmail.com
# This script includes functions for visualizing overall output of MaAsLin2 and
# individual associations as scatterplot and boxplot
# Date: `r paste(date())` 
# --------------

# Load libararies 
for( lib in c('ggplot2',"grid",'pheatmap')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}

# MaAsLin2 theme based on Nature journal requirements
nature_theme <- theme_bw() + theme(axis.text.x = element_text(size = 8, vjust = 1),
                                   axis.text.y = element_text(size = 8, hjust = 1),
                                   axis.title=element_text(size = 10  ),
                                   plot.title =element_text(size=7, face='bold'),
                                   legend.title=element_text(size=6, face='bold'),
                                   legend.text=element_text(size=6),
                                   axis.line = element_line(colour = 'black', size = .25),
                                   axis.line.x = element_line(colour = 'black', size = .25), 
                                   axis.line.y = element_line(colour = 'black', size = .25),
                                   panel.border = element_blank(), 
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank())

## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}


# MaAsLin2 heatmap function for overall view of associations 
maaslin2_heatmap <- function(output_results, title = NULL, cell_value = "qval", data_label = 'data', metadata_label = 'metadata',
                             border_color = "grey93", format =  NA, color = colorRampPalette(c("darkblue","grey90", "darkred"))(50),
                             col_rotate = 45) {#)
  # read MaAsLin output
  
  df <- read.table( output_results,
                    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE)
  
  if (dim(df)[1] < 1){
    print('There is no association to plot!')
    return (NULL)
  }
  
  # rotate colnames
  if (col_rotate == 45) {
    ## 'Overwrite' default draw_colnames with your own version 
    assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                      ns=asNamespace("pheatmap"))
  }
  metadata <- df$metadata
  data <- df$feature
  value <- NA
  # values to use for coloring the heatmap
  # and set the colorbar boundaries
  if (cell_value == "pval"){
    value <- -log(df$pval)*sign(df$coef)
    value <- pmax(-2, pmin(2, value))
    if (is.null(title)) title<-"Significant associations (-log(pval)*sign(coeff))"
  }else if(cell_value == "qval"){
    value <- -log(df$qval)*sign(df$coef)
    value <- pmax(-2, pmin(2, value))
    if (is.null(title)) title<-"Significant associations (-log(qval)*sign(coeff))"
  }else if(cell_value == "coef"){
    value <- df$coef
    if (is.null(title)) title<-"Significant associations (coeff)"
  }
  n <- length(unique(data))
  m <- length(unique(metadata))
  a = matrix(0, nrow=n, ncol=m)
  for (i in 1:length(metadata)){
    #if (abs(a[as.numeric(metadata[i]), as.numeric(metadata[i])]) > abs(value[i]))
    #  next
    a[as.numeric(data)[i], as.numeric(metadata)[i]] <- value[i]
  }
  rownames(a) <- levels(data)
  colnames(a) <- levels(metadata)
  p <- NULL
  tryCatch({ 
    p <- pheatmap::pheatmap(
      a, cellwidth = 5, cellheight = 5,   # changed to 3
      main = title,
      fontsize = 6,
      kmeans_k = NA,
      border=TRUE,
      show_rownames = T, show_colnames = T,
      scale="none",
      cluster_rows = FALSE, cluster_cols = TRUE,
      clustering_distance_rows = "euclidean", 
      clustering_distance_cols = "euclidean",
      legend=TRUE,
      border_color = border_color,
      color = color,
      treeheight_row=0,
      treeheight_col=0,
      display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "-", "")),  nrow(a)))
    }, error=function(err){ 
      logging::logerror("Unable to plot heatmap") 
    })
  return(p)
}

save_heatmap <- function(results_file, heatmap_file, title = NULL, cell_value = "qval", data_label = 'data', metadata_label = 'metadata',
                         border_color = "grey93", format =  NA, color = colorRampPalette(c("blue","grey90", "red"))(50)) {
  # generate a heatmap and save it to a pdf
  heatmap <- maaslin2_heatmap(results_file, title, cell_value, data_label, metadata_label, border_color, color)

  if (!is.null(heatmap)) {
    pdf(heatmap_file)
    print(heatmap)
    dev.off()  
  }
}

maaslin2_association_plots <- function(metadata, features, output_results, write_to='./')
{ 
  #MaAslin2 scatter plot function and theme
  
  # combine the data and metadata to one datframe using common rows
  # read MaAsLin output
  if (is.character(metadata)){
    metadata <- read.table( metadata, header = TRUE,
                          row.names = 1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  }
  if (is.character(features)){
    features <- read.table( features, header = TRUE,
                          row.names = 1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  }
  
  common_rows <- intersect(rownames(features), rownames(metadata))
  input_df_all <- cbind(features[common_rows, ,drop=FALSE], metadata[common_rows, ,drop=FALSE])
 
  # read MaAsLin output
  if (is.character(output_results)){
    output_df_all <- read.table( output_results, header = TRUE,
                               row.names = NULL, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  }else{
    output_df_all <- output_results
  }
  
  if (dim(output_df_all)[1] < 1){
    print('There is no association to plot!')
    return (NULL)
  }
 
  logging::loginfo("Plotting associations from most to least significant, grouped by metadata")
  metadata_types <- unlist(output_df_all[, 'metadata'])
  metadata_labels <- unlist(metadata_types[!duplicated(metadata_types)])
  metadata_number <- 1
  for (label in metadata_labels) {

    # for file name replace any non alphanumeric with underscore
    plot_file <- paste(write_to, "/", gsub("[^[:alnum:]_]", "_", label), ".pdf", sep="")
    data_index <- which(label==metadata_types)
    logging::loginfo("Plotting data for metadata number %s, %s", metadata_number, label) 
    pdf(plot_file, onefile=TRUE)

    for (i in data_index){
      x_label <- as.character(output_df_all[i, 'metadata'])
      y_label <- as.character(output_df_all[i, 'feature'])
      input_df <- input_df_all[c(x_label,y_label)]
      colnames(input_df) <- c("x", "y")
      # if Metadata is continuous generate a scatter plot
      # Continuous is defined as numerical with more than 2 values (to exclude binary data)
      temp_plot <- NULL
      if (is.numeric(input_df[1,'x']) & length(unique(input_df[['x']])) > 2){
        logging::loginfo("Creating scatter plot for continuous data, %s vs %s", x_label, y_label)
        temp_plot <- ggplot2::ggplot(data=input_df, 
          ggplot2::aes(as.numeric(as.character(x)), as.numeric(as.character(y)))) +
          ggplot2::geom_point( fill = 'darkolivegreen4', color = 'darkolivegreen4', alpha = .5, shape = 21, size = 1.5, stroke = 0.05) + 
          ggplot2::scale_x_continuous(limits=c(min(input_df['x']), max(input_df['x'])))+
          ggplot2::scale_y_continuous(limits=c(min(input_df['y']), max(input_df['y'])))+
          ggplot2::stat_smooth(method = "glm", color ='blue', na.rm = T)+ 
          ggplot2::guides(alpha='none')+ggplot2::labs("")+
          ggplot2::xlab(x_label) +  ggplot2::ylab(y_label) + nature_theme
      }else{
        # if Metadata is categorical generate a boxplot
        ### check if the variable is categorical
        logging::loginfo("Creating boxplot for catgorical data, %s vs %s", x_label, y_label)
        input_df['x'] <- sapply(input_df['x'], as.character) 
        temp_plot <- ggplot2::ggplot(data=input_df, aes(x, y)) + 
          ggplot2::geom_boxplot(aes(fill =input_df$x), notch = T,
                                outlier.alpha = 0.25, na.rm = T,
                                show.legend = F) + 
          ggplot2::scale_fill_brewer(palette="Spectral")
      
        # format the figure to default nature format, remove legend, add x/y labels
        temp_plot <- temp_plot + nature_theme +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          ggplot2::xlab(x_label) +  ggplot2::ylab(y_label) +
          theme(legend.position="none")
    }
    stdout <- capture.output(print(temp_plot),type="message")
    logging::logdebug(stdout)
  }
  dev.off()
  metadata_number <- metadata_number + 1
  }
}
