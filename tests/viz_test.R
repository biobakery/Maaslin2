# test case for th MaAsLin2 heatmap
maaslin2_output <- '~/Downloads/ali.txt'
srb_ddr <- maaslin2_heatmap ('~/Downloads/ali.txt', title = "", cell_value = "Q.value", data_label = 'Data', 
                             metadata_label = 'Metadata', border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(500))

ggsave(filename='~/Documents/01_srb_ddr.pdf', #'~/Dropbox (Huttenhower Lab)/hutlab/long/sulfur/figures/bugs/01_srb_ddr.pdf', 
       plot=srb_ddr$gtable, width = 135, height = 100, units = "mm", dpi = 350)

#test scatter
bugs <- read.table('~/Documents/HMP2/bugs-metabolites_bugs_residual_subject.tsv',  header = TRUE, row.names = 1,sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
metabolites <- read.table('~/Documents/HMP2/bugs-metabolites_metabolites_residual_subject.tsv',  header = TRUE, row.names = 1,sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
