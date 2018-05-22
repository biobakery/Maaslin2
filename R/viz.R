snippet header
# Author: Gholamali (Ali) rahnavard
# Email: gholamali.rahnavard@gmail.com
# This script includes functions for visualzing overall output of MaAsLin2 and
# inditual associations as scatterplot or boxplot
# Date: `r paste(date())` 
# --------------

# Load libararies 
library(ggplot2)
library(pheatmap)

# MaAsLin2 theme based on Nutur journal rquiremnts
nature_theme <- theme(axis.text.x = element_text(size = 8, vjust = 1),
                      axis.text.y = element_text(size = 8, hjust = 1),
                      axis.title=element_text(size = 10  ),
                      plot.title =element_text(size=7, face='bold'),
                      legend.title=element_text(size=6, face='bold'),
                      legend.text=element_text(size=6),
                      axis.line = element_line(colour = 'black', size = .25),
                      axis.line.x = element_line(colour = 'black', size = .25), axis.line.y = element_line(colour = 'black', size = .25)) 

# MaAsLin2 heatmap function for overall view of associations 
maaslin2_heatmap <- function(maaslin_output, title = "", cell_value = "Q.value", data_label = 'Data', metadata_label = 'Metadata',
                            border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(50)) {#)
  # read MaAsLin output
  df <- read.table( maaslin_output,
                    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE)
  metadata <- df$Variable
  data <- df$Feature
  value <- NA
  # values to use for coloring the heatmap
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
  a = matrix(0, nrow=n, ncol=m)
  for (i in 1:length(metadata)){
    #if (abs(a[as.numeric(metadata[i]), as.numeric(metadata[i])]) > abs(value[i]))
    #  next
    a[as.numeric(metadata)[i], as.numeric(data)[i]] <- value[i]
  }
  rownames(a) <- levels(metadata)
  colnames(a) <- levels(data)
  
  #colnames(a) <-  sapply(colnames(a), pcl.sub )
  
  p <- pheatmap(a, cellwidth = 7, cellheight = 7,   # changed to 3
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
                display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "|", "")),  nrow(a)))
  return(p)
}

maaslin2_scatter <- function(input_df, output_df, write_to_file = F, output_path='./'){ 
  #MaAslin2 scatter plot function and theme
  
  input_df = as.data.frame(t(rbind(bugs, metabolites)))
  for (i in (1:dim(output_df)[1])){
    scatter_plot[i] <- ggplot(data=bugs_metabolites,aes(output_df[i, 'variable'], output_df[i, 'feature'])) +
      geom_point( aes(), fill = 'darkolivegreen4', color = 'darkolivegreen4', alpha = .5, shape = 21, size = 1.5, stroke = 0.05) + 
      scale_x_continuous(limits=c(min(input_df[output_df[i, 'variable']]), max(input_df[output_df[i, 'variable']])))+
      scale_y_continuous(limits=c(min(input_df[output_df[i, 'feature']]), max(input_df[output_df[i, 'feature']])))+
      stat_smooth(method = "glm", color ='blue')+ 
      guides(alpha='none')+labs("")+
      xlab(output_df[i, 'variable']) +  ylab(output_df[i, 'feature']) + nature_theme+
      annotate(geom="text", x=(max(input_df[output_df[i, 'variable']],na.rm = T) - 
                              .5*(max(input_df[output_df[i, 'variable']], na.rm = T) - 
                              min(input_df[output_df[i, 'variable']], na.rm = T))),
               y=max(input_df[output_df[i, 'feature']], na.rm = T) -0.012, 
               label=paste("Spearman correlation: ", str(output_df[i, 'correllation']),"q-value = ",
                           str(output_df[i, 'fdr'])), color="black", size=rel(2), fontface="italic")+
      guides(legend.position=NULL)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    if (write_to_file)
      ggsave(filename=paste(output_path,'/' + str(i), '.pdf'), plot=scatter_plot[i],
             width = 6, height = 6, units = "cm", dpi = 300)
  }
  return (scatter_plot)
}