library(testthat)
library(Maaslin2)

expected_results <- read.table("expected_results.tsv", header = TRUE, stringsAsFactors=FALSE)

features <- read.table(system.file(package="Maaslin2","extdata","example1_features.txt"), header = TRUE, row.names = 1)
metadata <- read.table(system.file(package="Maaslin2","extdata","example1_metadata.txt"), header = TRUE, row.names = 1)

fit_data <- Maaslin2(features, metadata, 'output')
maaslin_results = read.table(file.path("output","significant_results.tsv"), header = TRUE, stringsAsFactors=FALSE)

expect_that(expected_results$metadata,equals(maaslin_results$metadata))
expect_that(expected_results$feature,equals(maaslin_results$feature))
expect_that(expected_results$coef,equals(maaslin_results$coef))
expect_that(expected_results$stderr,equals(maaslin_results$stderr))
expect_that(expected_results$N,equals(maaslin_results$N))
expect_that(expected_results$N.not.0,equals(maaslin_results$N.not.0))
expect_that(expected_results$pval.value,equals(maaslin_results$pval.value))
expect_that(expected_results$qval.value,equals(maaslin_results$qval.value))

# run again with different options to check turning off normalization/transform
# and that all results are provided even if significance is reduced
fit_data <- Maaslin2(features, metadata, 'output_nondefault_settings', normalization='NONE', transform='NONE', min_prevalence=0, max_significance=0.001, plot_heatmap=FALSE, plot_scatter=FALSE)
expect_that(length(fit_data$results$pval[fit_data$results$pval < 0.001]),equals(3))

