#!/usr/bin/env Rscript

##########
# MaAsLin2

# Copyright (c) 2018 Harvard School of Public Health

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########

# load in the required libraries, report an error if they are not installed
library("optparse")

# set the default choices
normalization_choices <- c("TSS","CLR","CSS","TSS","NONE","TMM")
analysis_method_choices <- c("LM","CPLM","ZICP","NEGBIN","ZINB")
transform_choices <- c("LOG","LOGIT","AST","NONE")

# set the default run options
args <- list()
args$input_data <- NULL
args$input_metadata <- NULL
args$output <- NULL
args$min_abundance <- 0.0
args$min_prevalence <- 0.1
args$normalization <- normalization_choices[1]
args$transform <- transform_choices[1]
args$analysis_method <- analysis_method_choices[1]

# add command line arguments
options <- OptionParser(usage = "%prog [options] <data.tsv> <metadata.tsv> <output_folder>")
options <- add_option(options, c("-a","--min_abundance"), type="double", dest="min_abundance", 
    default=args$min_abundance, help="The minimum abundance for each feature [ Default: %default ]")
options <- add_option(options, c("-p","--min_prevalence"), type="double", dest="min_prevalence", 
    default=args$min_prevalence, help="The minimum percent of samples for which a feature is detected [ Default: %default ]")
options <- add_option(options, c("-n","--normalization"), type="character", dest="normalization", 
    default=args$normalization, help=paste("The normalization method to apply [ Default: %default ] [ Choices:",toString(normalization_choices),"]"))
options <- add_option(options, c("-t","--transform"), type="character", dest="transform", 
    default=args$transform, help=paste("The transform to apply [ Default: %default ] [ Choices:",toString(transform_choices),"]"))
options <- add_option(options, c("-m","--analysis_method"), type="character", dest="analysis_method", 
    default=args$analysis_method, help=paste("The analysis method to apply [ Default: %default ] [ Choices:",toString(analysis_method_choices),"]"))

# main maaslin2 function with defaults set to the same as those used on the command line
Maaslin2 <- function(input_data, input_metadata, output, min_abundance=args$min_abundance, 
    min_prevalence=args$min_prevalence, normalization=args$normalization, transform=args$transform, 
    analysis_method=args$analysis_method)
{
    # read in the data and metadata
    data <- read.table(input_data, header=TRUE, sep = "\t", row.names = 1)
    metadata <- read.table(input_metadata, header=TRUE, sep = "\t", row.names = 1)

    # create an output folder if it does not exist
    if (!file.exists(output)) {
        print("Creating output folder")
        dir.create(output)
    }

    # check valid normalization option selected
    if (! normalization %in% normalization_choices) {
        stop(paste("Please select a normalization from the list of available options:", toString(normalization_choices)))
    }

    # check valid transform option selected
    if (! transform %in% transform_choices) {
        stop(paste("Please select a transform from the list of available options:", toString(transform_choices)))
    }

    if (analysis_method == "LM") {
        results <- fit.LM(data, metadata, normalization=normalization, transform=transform)
    } else if (analysis_method == "CPLM") {
        results <- fit.CPLM(data, metadata, normalization=normalization, transform=transform)
    } else if (analysis_method == "ZICP") {
        results <- fit.ZICP(data, metadata, normalization=normalization, transform=transform)
    } else if (analysis_method == "NEGBIN") {
        results <- fit.negbin(data, metadata, normalization=normalization, transform=transform)
    } else if (analysis_method == "ZINB") {
        results <- fit.ZINB(data, metadata, normalization=normalization, transform=transform)
    } else {
        stop(paste("Please select an analysis method from the list of available options:", toString(analysis_method_choices)))
    }

    # count the total values (non NA) for each feature
    results$N <- apply(results, 1, FUN = function(x) length(data[,x[1]]))
    results$N.not.zero <- apply(results, 1, FUN = function(x) length(which(data[,x[1]] > 0)))

    # write the results to a file
    write.table(results[c("metadata","feature","metadata","coef","N","N.not.zero","pval","qval")], file=file.path(output,"results.tsv"), sep="\t", quote=FALSE, col.names=c("Variable","Feature","Value","Coefficient","N","N.not.0","P.value","Q.value"))

}

# this evaluates to true if script is being called directly as an executable
if(identical(environment(), globalenv()) &&
    !length( grep( "^source\\(", sys.calls()))) 
{

    # source all R in Maaslin2 package, relative to this folder (same method as original maaslin)
    script_options <-commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=","",script_options[grep("--file=", script_options)])
    script_dir <- dirname(script_path)
    script_name <- basename(script_path)
    
    for(R_file in dir(script_dir, pattern = "*.R"))
    {
        if(! ( R_file == script_name )) source( file.path(script_dir, R_file) )
    }

    # get command line options and positional arguments
    parsed_arguments = parse_args(options,positional_arguments=TRUE)
    current_args <- parsed_arguments[["options"]]
    positional_args <- parsed_arguments[["args"]]
    # check three positional arguments are provided
    if(length(positional_args)!= 3) {
        print_help(options)
        stop("Please provide the required positional arguments <data.tsv> <metadata.tsv> <output_folder>")
    }
    
    # call maaslin with the command line options
    Maaslin2(positional_args[1], positional_args[2], positional_args[3],
        current_args$min_abundance, current_args$min_prevalence, 
        current_args$normalization, current_args$transform,
        current_args$analysis_method) 
}
