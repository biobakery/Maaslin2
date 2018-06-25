

# test case for th MaAsLin2 heatmap
output_df = '~/Documents/Hutlab/maaslin2/tests/example/maaslin2_output.txt'

maaslin2_heatmap_overview_plot <- maaslin2_heatmap (output_df, title = "", cell_value = "Q.value", data_label = 'Data', 
                             metadata_label = 'Metadata', border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(500))

ggsave(filename='~/Documents/Hutlab/maaslin2/tests/example/maaslin2_heatmap_overview_plot.pdf',  
       plot=maaslin2_heatmap_overview_plot$gtable, width = 135, height = 100, units = "mm", dpi = 350)

#test scatter
input_df = '~/Documents/Hutlab/maaslin2/tests/example/Bug_Stool.tsv'


write_to_file = F
output_path='~/Documents/Hutlab/maaslin2/tests/example/'
plots <- maaslin2_association_plots(input_df, output_df, write_to_file, output_path)
plots[[1]]
