#### Remove all variables from the workspace ####
rm(list = ls())

# First set directory to the repo and source all functions
# setwd('/Users/hmallick/Dropbox (Personal)/Repos/Maaslin2')

# Install and load Packages (Specific to this Example)
if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('pkgmaker', 'data.table', 'tibble')

# Load Utility Functions
pkgmaker::source_files('./R', '*.R')

# Load Datasets
features<-data.table::fread('./tests/example1_features.txt', header=TRUE)
metadata<-data.table::fread('./tests/example1_metadata.txt', header=TRUE)

# Convert to Proper Data Frame with Rownames
features<-tibble::column_to_rownames(as.data.frame(features), '#')
metadata<-tibble::column_to_rownames(as.data.frame(metadata), '#')

# Example Run 1 (LM + TSS + AST)
LM_example_1<-fit.LM(features, metadata,normalization='TSS', transformation='AST')

# Example Run 2 (LM + CLE + NONE)
LM_example_2<-fit.LM(features, metadata, normalization='CLR', transformation='NONE')

# Example Run 3 (LM + CLE + NONE)
LM_example_3<-fit.LM(features, metadata, normalization='TSS', transformation='LOG')
