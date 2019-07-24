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
    standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)
maaslin_results = read.table(file.path("output","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run1$metadata,equals(maaslin_results$metadata))
expect_that(expected_results_run1$feature,equals(maaslin_results$feature))
expect_that(expected_results_run1$N,equals(maaslin_results$N))
expect_that(expected_results_run1$N.not.0,equals(maaslin_results$N.not.0))
expect_that(expected_results_run1$pval.value,equals(maaslin_results$pval.value))
expect_that(expected_results_run1$qval.value,equals(maaslin_results$qval.value))

# run with with dysbiosis as a continous variable
fit_data <- Maaslin2(features, metadata, 'output2', normalization = "NONE", transform = "AST",
    fixed_effects = c('diagnosis', 'dysbiosis', 'antibiotics', 'age'),
    random_effects = c('site', 'subject'),
    standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)
maaslin_results = read.table(file.path("output2","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results_run2$metadata,equals(maaslin_results$metadata))
expect_that(expected_results_run2$feature,equals(maaslin_results$feature))
expect_that(expected_results_run2$N,equals(maaslin_results$N))
expect_that(expected_results_run2$N.not.0,equals(maaslin_results$N.not.0))
expect_that(expected_results_run2$pval.value,equals(maaslin_results$pval.value))
expect_that(expected_results_run2$qval.value,equals(maaslin_results$qval.value))
