# load Maaslin2
library(Maaslin2)

# test case for th MaAsLin2 heatmap
output_path = '~/Documents/mgx/maaslin2/tests/'
metadata_path <- '~/Documents/mgx/maaslin2/tests/example1_metadata.txt'
#
features_path <- '~/Documents/mgx/maaslin2/tests/example1_features.txt'
#
Maaslin2(features_path, metadata_path, output_path)
output_results = '~/Documents/mgx/maaslin2/tests/significant_results.tsv'
maaslin2_heatmap_overview_plot <- maaslin2_heatmap(output_results, title = "", cell_value = "Q.value", data_label = 'Data', 
                                                    metadata_label = 'Metadata', border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(500))

ggsave(filename='~/Documents/Hutlab/maaslin2/tests/maaslin2_heatmap_overview_plot.pdf',  
       plot=maaslin2_heatmap_overview_plot$gtable, width = 135, height = 100, units = "mm", dpi = 350)

# test association plots


write_to_file = T
plots <- maaslin2_association_plots(metadata_path, features_path, output_results, write_to_file, write_to = '~/Documents/Hutlab/maaslin2/tests/' )
plots[[1]]
plots[[2]]

