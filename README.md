# MaAsLin2 User Manual #

MaAsLin2 is the next generation of MaAsLin.

[MaAsLin](https://huttenhower.sph.harvard.edu/maaslin) is a multivariate statistical framework that finds associations between clinical metadata and potentially high-dimensional experimental data.

If you use the MaAsLin2 software, please cite our manuscript: Himel Mallick, Timothy L. Tickle, Lauren J. McIver, Gholamali Rahnavard, George Weingart, Joseph N. Paulson, Siyuan Ma, Boyu Ren, Emma Schwager, Ayshwarya Subramanian, Eric A. Franzosa, Hector Corrada Bravo, Curtis Huttenhower. "Multivariable Association in Population-scale Meta'omic Surveys" (In Preparation).

If you have questions, please email the [MaAsLin Users Google Group](https://groups.google.com/forum/#!forum/maaslin-users).

--------------------------------------------

## Contents ##
* [Description](#markdown-header-description)
* [Requirements](#markdown-header-requirements)
* [Installation](#markdown-header-installation)
* [How to Run](#markdown-header-how-to-run)
    * [Input Files](#markdown-header-input-files)
    * [Output Files](#markdown-header-output-files)
    * [Run a Demo](#markdown-header-run-a-demo)
    * [Options](#markdown-header-options)
* [Troubleshooting](#markdown-header-troubleshooting)

## Description ##

MaAsLin2 was developed to find associations between microbiome multi'omics features and complex metadata in population-scale epidemiological studies. The software includes multiple analysis methods, normalization, and transform options to customize analysis for your specific study. 

## Requirements ##

MaAsLin2 is an R package that can be run on the command line or as an R function. It requires the following R packages included in Biocondutor and CRAN (Comprehensive R Archive Network). Please install these packages before running MaAsLin2.

* Bioconductor packages
    * [edgeR: Empirical Analysis of Digital Gene Expression Data in R](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
    * [metagenomeSeq: Statistical analysis for sparse high-throughput sequencing](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)
    * These packages can be installed through Bioconductor by first sourcing biocLite with ``source("https://bioconductor.org/biocLite.R")`` and then installing each package with ``biocLite("edgeR")``.
* CRAN packages
    * [pscl: Political Science Computational Laboratory](https://cran.r-project.org/web/packages/pscl/pscl.pdf)
    * [pbapply: Adding Progress Bar to '*apply' Functions](https://cran.rstudio.com/web/packages/pbapply/index.html)
    * [car: Companion to Applied Regression](https://cran.rstudio.com/web/packages/car/index.html)
    * [nlme: Linear and Nonlinear Mixed Effects Models](https://cran.rstudio.com/web/packages/nlme/index.html)
    * [dplyr: A Grammer of Data Manipulation](https://cran.rstudio.com/web/packages/dplyr/index.html)
    * [vegan: Community Ecology Package](https://cran.rstudio.com/web/packages/vegan/index.html)
    * [chemometrics: Multivariate Statistical Analysis in Chemometrics](https://cran.rstudio.com/web/packages/chemometrics/index.html) 
    * [ggplot2: Create Elegant Data Visualizations Using the Grammer of Graphics](https://cran.rstudio.com/web/packages/ggplot2/index.html)
    * [pheatmap: Pretty Heatmaps](https://cran.rstudio.com/web/packages/pheatmap/index.html)
    * [cplm: Compound Poisson Linear Models](https://cran.rstudio.com/web/packages/cplm/index.html)
    * [hash: Full feature implementation of hash/associated arrays/dictionaries](https://cran.rstudio.com/web/packages/hash/index.html)
    * [logging: R logging package](https://cran.rstudio.com/web/packages/logging/index.html)
    * These packages can be installed in R with ``install.packages('pscl')`` or from the command line ``$ R -q -e "install.packages('pscl', repos='http://cran.r-project.org')"`` individually (for those packages which you do not yet have installed) or as a set by providing the complete list as a vector.

## Installation ##

MaAsLin2 can be run from the command line or as an R function. If only running from the command line, you do not need to install the MaAsLin2 package.

1. Download the source: [MaAsLin2.tar.gz](https://bitbucket.org/biobakery/maaslin2/get/tip.tar.gz)
2. Decompress the download: ``$ tar xzvf maaslin2.tar.gz``
3. Install the package (only required if running as an R function): ``$ R CMD INSTALL maaslin2``

## How to Run ##

MaAsLin2 can be run from the command line or as an R function. Both methods require the same arguments, have the same options, and use the same default settings.

To run from the command line: ``$ Maaslin2.R $DATA $METADATA $OUTPUT``

* Provide the full path to the MaAsLin2 executable (ie ./R/Maaslin2.R if you are in the source folder).
* Replace ``$DATA`` with the path to your data (or features) file. 
* Replace ``$METADATA`` with the path to your metadata file.
* Replace ``$OUTPUT`` with the path to the folder to write the output.

To run from R as a function: 

```
$ R
> library(Maaslin2)
> Maaslin2(data, metadata, output)
```

### Input Files ###

MaAsLin2 requires two input files.

1. Data (or features) file
    * This file is tab-delimited formatted with features as columns and samples as rows.
    * Possible features in this file include data like taxonomic or gene abundances.
2. Metadata file
    * This file is tab-delimited formatted with metadata as columns and samples as rows.
    * Possible metadata in this file include gender or age.

Please note the same samples must be included in both files. Also these samples must be in the same order in both files.

### Output Files ###

MaAsLin2 generates two types of output files: data and visualization.

1. Data output files
    * all_results.tsv : This file contains all of the association results ordered by increasing q-value.
    * significant_results.tsv : This file is a subset of the data in the first file. It only includes those associations with q-values less than or equal to the significance threshold.
    * maaslin2.log : This file contains all of the debug information for the run. It includes all settings, warnings, errors, and steps run.
2. Visualization output files
    * heatmap.pdf : This file contains a heatmap of the significant associations.
    * 1*.pdf : These files are scatter plots with one generated for each significant association.

### Run a Demo ###

Example input files can be found in the tests folder of the MaAsLin2 source. 

To run: ``$ Maaslin2.R maaslin2/tests/example1_data.txt maaslin2/tests/example1_metadata.txt demo_output``

When running this command, all output files will be written to a folder named ``demo_output``.

### Options ###

Run MaAsLin2 help to print a list of the options and the default settings.


```
$ Maaslin2.R --help
Usage: ./R/Maaslin2.R [options] <data.tsv> <metadata.tsv> <output_folder>


Options:
	-h, --help
		Show this help message and exit

	-a MIN_ABUNDANCE, --min_abundance=MIN_ABUNDANCE
		The minimum abundance for each feature [ Default: 0 ]

	-p MIN_PREVALENCE, --min_prevalence=MIN_PREVALENCE
		The minimum percent of samples for which a feature is detected at minimum abundance [ Default: 0.1 ]

	-s MAX_SIGNIFICANCE, --max_significance=MAX_SIGNIFICANCE
		The q-value threshold for significance [ Default: 0.25 ]

	-n NORMALIZATION, --normalization=NORMALIZATION
		The normalization method to apply [ Default: TSS ] [ Choices: TSS, CLR, CSS, NONE, TMM ]

	-t TRANSFORM, --transform=TRANSFORM
		The transform to apply [ Default: LOG ] [ Choices: LOG, LOGIT, AST, NONE ]

	-m ANALYSIS_METHOD, --analysis_method=ANALYSIS_METHOD
		The analysis method to apply [ Default: CPLM ] [ Choices: CPLM, LM, NEGBIN, ZICP, ZINB ]
```

## Troubleshooting ##

1. Question: When I run from the command line I see the error ``Maaslin2.R: command not found``. How do I fix this? 
    * Answer: Provide the full path to the executable when running Maaslin2.R.
2. Question: When I run as a function I see the error ``Error in library(Maaslin2): there is no package called 'Maaslin2'``. How do I fix this? 
    * Answer: Install the R package and then try loading the library again.
3. Question: When I try to install the R package I see errors about dependencies not being installed. Why is this?
    * Answer: Installing the R package will not automatically install the packages MaAsLin2 requires. Please install the dependencies and then install the MaAsLin2 R package.

