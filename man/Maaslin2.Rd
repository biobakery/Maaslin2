\name{Maaslin2}
\alias{Maaslin2}
\title{
MaAsLin2 is the next generation of MaAsLin. MaAsLin is a multivariate
    statistical framework that finds associations between clinical 
    metadata and potentially high-dimensional experimental data.
}
\description{
MaAsLin2 was developed to find associations between microbiome 
    multi'omics features and complex metadata in population-scale 
    epidemiological studies. The software includes multiple 
    analysis methods, normalization, and transform options to 
    customize analysis for your specific study.
}
\usage{
Maaslin2(
    input_data,
    input_metadata,
    output,
    min_abundance = 0.0,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 0.25,
    random_effects = NULL,
    fixed_effects = NULL,
    correction = "BH",
    standardize = TRUE,
    cores = 1,
    plot_heatmap = TRUE,
    plot_scatter = TRUE
)
}
\arguments{
    \item{input_data}{
    The tab-delimited input file of features.
}
    \item{input_metadata}{
    The tab-delimited input file of metadata.
}
    \item{output}{
    The output folder to write results.
}
    \item{min_abundance}{
    The minimum abundance for each feature.
}
    \item{min_prevalence}{
    The minimum percent of samples for which a feature 
        is detected at minimum abundance.
}
    \item{max_significance}{
    The q-value threshold for significance.
}
    \item{normalization}{
    The normalization method to apply.
}
    \item{transform}{
    The transform to apply.
}
    \item{analysis_method}{
    The analysis method to apply.
}
    \item{random_effects}{
    The random effects for the model, comma-delimited for multiple effects.
}
    \item{fixed_effects}{
    The fixed effects for the model, comma-delimited for multiple effects.
}
    \item{correction}{
    The correction method for computing the q-value.
}
    \item{standardize}{
    Apply z-score so continuous metadata are on the same scale.
}
    \item{plot_heatmap}{
    Generate a heatmap for the significant associations.
}
    \item{plot_scatter}{
    Generate scatter plots for the significant associations.
}
    \item{cores}{
    The number of R processes to run in parallel.
}
}
\value{
    Data.frame containing the results from applying the model.
}
\author{
    Himel Mallick<hmallick@broadinstitute.org>,\cr
    Ali Rahnavard<rah@broadinstitute.org>,\cr
    Maintainers: Lauren McIver<lauren.j.mciver@gmail.com>,\cr
}
\examples{
    input_data <- system.file(
        'extdata','HMP2_taxonomy.tsv', package="Maaslin2")
    input_metadata <-system.file(
        'extdata','HMP2_metadata.tsv', package="Maaslin2")
    fit_data <- Maaslin2(
        input_data, input_metadata,'demo_output', transform = "AST",
        fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'),
        random_effects = c('site', 'subject'),
        standardize = FALSE)
}
