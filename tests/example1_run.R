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

# Example Run 1 (LM + CLR + NONE)
Example_1<-fit.LM(features, metadata,normalization='CLR', transformation='NONE')

# Example Run 2 (LM + CSS + LOG)
Example_2<-fit.LM(features, metadata,normalization='CSS', transformation='LOG')

# Example Run 3 (LM + TSS + LOGIT)
Example_3<-fit.LM(features, metadata,normalization='TSS', transformation='AST')

# Example Run 4 (CPLM + TSS + AST)
Example_4<-fit.CPLM(features, metadata, normalization='TSS', transformation='LOG')

# Example Run 5 (ZICP + NONE + NONE)
Example_5<-fit.ZICP(features, metadata, normalization='NONE', transformation='NONE')

# Example Run 6 (negbin + TMM + NONE)
Example_6<-fit.negbin(features, metadata, normalization='TMM', transformation='NONE')

# Example Run 7 (ZINB + NONE + NONE)
Example_7<-fit.ZINB(features, metadata, normalization='NONE', transformation='NONE')



