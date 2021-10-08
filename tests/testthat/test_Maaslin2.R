library(testthat)
library(Maaslin2)

expected_results_run1 <- read.table("expected_results_run1.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")
expected_results_run2 <- read.table("expected_results_run2.tsv", header = TRUE, stringsAsFactors=FALSE, sep="\t")

features <- read.table(system.file(package="Maaslin2","extdata","HMP2_taxonomy.tsv"), header = TRUE, row.names = 1, sep="\t")
metadata <- read.table(system.file(package="Maaslin2","extdata","HMP2_metadata.tsv"), header = TRUE, row.names = 1, sep="\t")

# run with dysbiosis as a nested categorical variable
fit_data <- Maaslin2(features, metadata, 'output', transform = "AST",
    fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
    random_effects = c('site', 'subject'),
    reference = 'diagnosis,nonIBD',
    standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)
maaslin_results = read.table(file.path("output","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run1$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run1$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run1$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(expected_results_run1$N.not.0[1:50],equals(maaslin_results$N.not.0[1:50]))
expect_that(round(as.numeric(expected_results_run1$pval.value[1:50]),10),equals(round(as.numeric(maaslin_results$pval.value[1:50]),10)))
expect_that(round(as.numeric(expected_results_run1$qval.value[1:50]),10),equals(round(as.numeric(maaslin_results$qval.value[1:50]),10)))

# run with with dysbiosis as a continous variable
fit_data <- Maaslin2(features, metadata, 'output2', normalization = "NONE", transform = "AST",
    fixed_effects = c('diagnosis', 'dysbiosis', 'antibiotics', 'age'),
    random_effects = c('site', 'subject'),
    reference = 'diagnosis,nonIBD',
    standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)
maaslin_results = read.table(file.path("output2","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run2$metadata[1:50],equals(maaslin_results$metadata[1:50]))
expect_that(expected_results_run2$feature[1:50],equals(maaslin_results$feature[1:50]))
expect_that(round(expected_results_run2$N[1:50],10),equals(round(maaslin_results$N[1:50],10)))
expect_that(expected_results_run2$N.not.0[1:50],equals(maaslin_results$N.not.0[1:50]))
expect_that(round(as.numeric(expected_results_run2$pval.value[1:50]),10),equals(round(as.numeric(maaslin_results$pval.value[1:50]),10)))
expect_that(round(as.numeric(expected_results_run2$qval.value[1:50]),10),equals(round(as.numeric(maaslin_results$qval.value[1:50]),10)))


# test the prevalence and abundance filtering
logging::basicConfig(level = 'FINEST')
logging::addHandler(logging::writeToFile,
                    file = "tmp.log", level = "DEBUG")
logging::setLevel(20, logging::getHandler('basic.stdout'))

relab_matrix  <- head(features)
counts_matrix <- ceiling(relab_matrix * 100000 )

# the data should be the same whether or not the input data is relative or absolute.
relab_filtered  <- do_prevalence_abundance_filtering(unfiltered_data = relab_matrix, min_abundance = .1, min_prevalence = .1, min_variance = 0)
counts_filtered <- do_prevalence_abundance_filtering(unfiltered_data = counts_matrix, min_abundance = .1, min_prevalence = .1, min_variance = 0)
expect_that(ceiling(relab_filtered*100000), equals(counts_filtered))


# the prevalence filtering should be the same.
relab_noab_filtered  <- do_prevalence_abundance_filtering(unfiltered_data = relab_matrix, min_abundance = 0, min_prevalence = .1, min_variance = 0)
counts_noab_filtered <- do_prevalence_abundance_filtering(unfiltered_data = counts_matrix, min_abundance = 0, min_prevalence = .1, min_variance = 0)
expect_that(ceiling(relab_noab_filtered*100000), equals(counts_noab_filtered))

# the abund filtering should be the same.
relab_noprev_filtered  <- do_prevalence_abundance_filtering(unfiltered_data = relab_matrix, min_abundance = .1, min_prevalence = 0, min_variance = 0)
counts_noprev_filtered <- do_prevalence_abundance_filtering(unfiltered_data = counts_matrix, min_abundance = .1, min_prevalence = 0, min_variance = 0)
expect_that(ceiling(relab_noprev_filtered*100000), equals(counts_noprev_filtered))


# this shows the desired behavior 
# -- for relative abundances, abundances are treated as a percent
# -- for counts it is treated as a percentage if < 0 
# -- fir counts it is treated as a count if > 0
relab_percent_filtered  <- do_prevalence_abundance_filtering(unfiltered_data = relab_matrix, min_abundance = 20, min_prevalence = .1, min_variance = 0)
counts_percent_filtered <- do_prevalence_abundance_filtering(unfiltered_data = counts_matrix, min_abundance = .2, min_prevalence = .1, min_variance = 0)
counts_ammount_filtered <- do_prevalence_abundance_filtering(unfiltered_data = counts_matrix, min_abundance = 20, min_prevalence = .1, min_variance = 0)
expect_that(ceiling(relab_percent_filtered*100000), equals(counts_percent_filtered))
expect_false(ncol(relab_percent_filtered) == ncol(counts_ammount_filtered))
