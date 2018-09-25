#snippet header
# Author: Gholamali (Ali) rahnavard
# Email: gholamali.rahnavard@gmail.com
# This script includes functions for visualizing overall output of MaAsLin2 and
# individual associations as scatterplot and boxplot
# Date: `r paste(date())` 
# --------------

# Load libararies 
library(ggplot2)
library(pheatmap)
#library(ggcorrplot)

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

# MaAsLin2 heatmap function for overall view of associations 
maaslin2_heatmap <- function(output_results, title = "", cell_value = "Q.value", data_label = 'Data', metadata_label = 'Metadata',
                             border_color = "grey93", format =  NA, color = colorRampPalette(c("darkblue","grey90", "darkred"))(50)) {#)
  # read MaAsLin output
  
  df <- read.table( output_results,
                    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE)
  
  if (dim(df)[1] < 1){
    print('There is no association to plot!')
    return (NULL)
  }
  metadata <- df$Metadata
  data <- df$Feature
  value <- NA
  # values to use for coloring the heatmap
  # and set the colorbar boundaries
  if (cell_value == "P.value"){
    value <- -log(df$P.value)*sign(df$Coefficient)
    value <- pmax(-2, pmin(2, value))
  }else if(cell_value == "Q.value"){
    value <- -log(df$Q.value)*sign(df$Coefficient)
    value <- pmax(-2, pmin(2, value))
  }else if(cell_value == "Coefficient"){
    value <- df$Coefficient
  }
  n <- length(unique(metadata))
  m <- length(unique(data))
  print(n, m)
  a = matrix(0, nrow=n, ncol=m)
  print(dim(a))
  for (i in 1:length(metadata)){
    #if (abs(a[as.numeric(metadata[i]), as.numeric(metadata[i])]) > abs(value[i]))
    #  next
    a[as.numeric(metadata)[i], as.numeric(data)[i]] <- value[i]
  }
  rownames(a) <- levels(metadata)
  colnames(a) <- levels(data)
  p <- pheatmap::pheatmap(
    a, cellwidth = 7, cellheight = 7,   # changed to 3
    main = title,
    fontsize = 6,
    kmeans_k = NA,
    border=TRUE,
    show_rownames = T, show_colnames = T,
    scale="none",
    #clustering_method = "complete",
    cluster_rows = FALSE, cluster_cols = TRUE,
    clustering_distance_rows = "euclidean", 
    clustering_distance_cols = "euclidean",
    legend=TRUE,
    border_color = border_color,
    color = color,
    treeheight_row=0,
    treeheight_col=0,
    display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "-", "")),  nrow(a))
  )
  return(p)
}

save_heatmap <- function(results_file, heatmap_file, title = "", cell_value = "Q.value", data_label = 'Data', metadata_label = 'Metadata',
                         border_color = "grey93", format =  NA, color = colorRampPalette(c("blue","grey90", "red"))(50)) {
  # generate a heatmap and save it to a pdf
  heatmap <- maaslin2_heatmap(results_file, title, cell_value, data_label, metadata_label, border_color, color)
  
  ggplot2::ggsave(filename=heatmap_file, plot=heatmap$gtable, width = 135, height = 100, units = "mm", dpi = 350)
}

maaslin2_association_plots <- function(metadata_path, features_path,
                                       output_results, write_to_file = F, write_to='./')
{ 
  #MaAslin2 scatter plot function and theme
  
  # combine the data and metadata to one datframe using common rows
  # read MaAsLin output
  metadata <- read.table( metadata_path, header = TRUE,
                          row.names = 1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  features <- read.table( features_path, header = TRUE,
                          row.names = 1, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  
  common_rows <- intersect(rownames(features), rownames(metadata))
  input_df_all <- cbind(features[common_rows,], metadata[common_rows,])
  
  # read MaAsLin output
  output_df_all <- read.table( output_results, header = TRUE,
                               row.names = NULL, sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
  
  if (dim(output_df_all)[1] < 1){
    print('There is no association to plot!')
    return (NULL)
  }
  # a list to store association(scatter or boxplot) plot of all associations 
  association_plot <- vector(mode="list", length=dim(output_df_all)[1])
  
  for (i in 1:dim(output_df_all)[1]){
    #print(i)
    #i <- 1
    x_label <- as.character(output_df_all[i, 'Metadata'])
    y_label <- as.character(output_df_all[i, 'Feature'])
    input_df <- input_df_all[c(x_label,y_label)]
    colnames(input_df) <- c("x", "y")
    # if Metadata is continuous generate a scatter plot
    temp_plot <- NULL
    if (is.numeric(input_df[1,'x'])){
      temp_plot <- ggplot2::ggplot(data=input_df, 
        ggplot2::aes(as.numeric(as.character(x)), as.numeric(as.character(y)))) +
        ggplot2::geom_point( fill = 'darkolivegreen4', color = 'darkolivegreen4', alpha = .5, shape = 21, size = 1.5, stroke = 0.05) + 
        ggplot2::scale_x_continuous(limits=c(min(input_df['x']), max(input_df['x'])))+
        ggplot2::scale_y_continuous(limits=c(min(input_df['y']), max(input_df['y'])))+
        ggplot2::stat_smooth(method = "glm", color ='blue', na.rm = T)+ 
        ggplot2::guides(alpha='none')+ggplot2::labs("")+
        ggplot2::xlab(x_label) +  ggplot2::ylab(y_label) + nature_theme
    }else{
      # if Metadata is categorical generate a Jitter plot with boxplot
      ### check if the variable is categorical
      temp_plot <- ggplot2::ggplot(data=input_df, aes(x, y)) + 
        ggplot2::geom_boxplot(notch = TRUE) +
        ggplot2::geom_jitter(position = position_jitter(0.5), ggplot2::aes(colour = x))
      
      # fomrat the figure to defualt nature format and add some 
      # statistics (Q.value and Coefficient) to the plot 
      temp_plot <- temp_plot + nature_theme +
        annotate(geom="text", x=(max(input_df['x'],na.rm = T) - 
                                   .5*(max(input_df['y'], na.rm = T) - 
                                         min(input_df['x'], na.rm = T))),
                 y=max(input_df['y'], na.rm = T) -0.012, 
                 label=paste("Spearman correlation: ", str(output_df[i, 'Coefficient']),"q-value = ",
                             str(output_df[i, 'Q.value'])), color="black", size=rel(2), fontface="italic")+
        ggplot2::guides(legend.position=NULL)+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
    association_plot[[i]] <- temp_plot 
    if (write_to_file){
      ggplot2::ggsave(filename=paste(write_to,'/', i, '.pdf', sep = ''), plot=temp_plot,
                      width = 6.5, height = 5, units = "cm", dpi = 300)
      #dev.off()
    }
  }
  return (association_plot)
}
