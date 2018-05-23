

# test case for th MaAsLin2 heatmap

maaslin2_heatmap_overview_plot <- maaslin2_heatmap (output_df, title = "", cell_value = "Q.value", data_label = 'Data', 
                             metadata_label = 'Metadata', border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(500))

ggsave(filename='~/Documents/maaslin2_heatmap_overview_plot',  
       plot=maaslin2_heatmap_overview_plot$gtable, width = 135, height = 100, units = "mm", dpi = 350)

#test scatter
input_df = '~/Documents/input.txt'
output_df = '~/Documents/output.txt'

write_to_file = F
output_path='~/Documents/'
plots <- maaslin2_association_plots(input_df, output_df, write_to_file, output_path)
plots[[1]]
