#!/usr/bin/env Rscript

###############################################################################
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
###############################################################################

# load in the required libraries, report an error if they are not installed

for (lib in c('optparse', 'logging', 'data.table', 'dplyr')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

###############################################################
# If running on the command line, load other Maaslin2 modules #
###############################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
        !length(grep("^source\\(", sys.calls()))) {
    # source all R in Maaslin2 package, relative to this folder
    # same method as original maaslin
    script_options <- commandArgs(trailingOnly = FALSE)
    script_path <-
        sub("--file=", "", script_options[grep("--file=", script_options)])
    script_dir <- dirname(script_path)
    script_name <- basename(script_path)
    
    for (R_file in dir(script_dir, pattern = "*.R$"))
    {
        if (!(R_file == script_name))
            source(file.path(script_dir, R_file))
    }
    library("stats", "utils", "grDevices", quietly = TRUE)
}

###########################
# Set the default options #
###########################

normalization_choices <- c("TSS", "CLR", "CSS", "NONE", "TMM")
analysis_method_choices_names <-
    c("LM", "CPLM", "NEGBIN", "ZINB")
transform_choices <- c("LOG", "LOGIT", "AST", "NONE")
valid_choice_method_norm <- hash::hash()
valid_choice_method_norm[[analysis_method_choices_names[3]]] <-
    normalization_choices[3:5]
valid_choice_method_norm[[analysis_method_choices_names[4]]] <-
    normalization_choices[3:5]
valid_choice_method_transform <- analysis_method_choices_names[1:2]
valid_choice_transform_norm <- hash::hash()
valid_choice_transform_norm[[transform_choices[2]]] <-
    normalization_choices[c(1, 4)]
valid_choice_transform_norm[[transform_choices[3]]] <-
    normalization_choices[c(1, 4)]
correction_choices <-
    c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY")

# set the default run options
args <- list()
args$input_data <- NULL
args$input_metadata <- NULL
args$output <- NULL
args$min_abundance <- 0.0
args$min_prevalence <- 0.1
args$min_variance <- 0.0
args$max_significance <- 0.25
args$normalization <- normalization_choices[1]
args$transform <- transform_choices[1]
args$analysis_method <- analysis_method_choices_names[1]
args$random_effects <- NULL
args$fixed_effects <- NULL
args$correction <- correction_choices[1]
args$standardize <- TRUE
args$plot_heatmap <- TRUE
args$heatmap_first_n <- 50
args$plot_scatter <- TRUE
args$max_pngs <- 10
args$save_scatter <- FALSE
args$cores <- 1
args$save_models <- FALSE
args$reference <- NULL

##############################
# Add command line arguments #
##############################

options <-
    optparse::OptionParser(usage = paste(
        "%prog [options]",
        " <data.tsv> ",
        "<metadata.tsv> ",
        "<output_folder>" 
        )
    )
options <-
    optparse::add_option(
        options,
        c("-a", "--min_abundance"),
        type = "double",
        dest = "min_abundance",
        default = args$min_abundance,
        help = paste0("The minimum abundance for each feature",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-p", "--min_prevalence"),
        type = "double",
        dest = "min_prevalence",
        default = args$min_prevalence,
        help = paste0("The minimum percent of samples for which",
            "a feature is detected at minimum abundance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-b", "--min_variance"),
        type = "double",
        dest = "min_variance",
        default = args$min_variance,
        help = paste0("Keep features with variances",
            "greater than value",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-s", "--max_significance"),
        type = "double",
        dest = "max_significance",
        default = args$max_significance,
        help = paste0("The q-value threshold for significance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-n", "--normalization"),
        type = "character",
        dest = "normalization",
        default = args$normalization,
        help = paste(
            "The normalization method to apply",
            "[ Default: %default ] [ Choices:",
            toString(normalization_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-t", "--transform"),
        type = "character",
        dest = "transform",
        default = args$transform,
        help = paste(
            "The transform to apply [ Default: %default ] [ Choices:",
            toString(transform_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-m", "--analysis_method"),
        type = "character",
        dest = "analysis_method",
        default = args$analysis_method,
        help = paste(
            "The analysis method to apply [ Default: %default ] [ Choices:",
            toString(analysis_method_choices_names),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-r", "--random_effects"),
        type = "character",
        dest = "random_effects",
        default = args$random_effects,
        help = paste("The random effects for the model, ",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-f", "--fixed_effects"),
        type = "character",
        dest = "fixed_effects",
        default = args$fixed_effects,
        help = paste("The fixed effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: all ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-c", "--correction"),
        type = "character",
        dest = "correction",
        default = args$correction,
        help = paste("The correction method for computing",
            "the q-value [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-z", "--standardize"),
        type = "logical",
        dest = "standardize",
        default = args$standardize,
        help = paste("Apply z-score so continuous metadata are on",
            "the same scale [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-l", "--plot_heatmap"),
        type = "logical",
        dest = "plot_heatmap",
        default = args$plot_heatmap,
        help = paste("Generate a heatmap for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-i", "--heatmap_first_n"),
        type = "double",
        dest = "heatmap_first_n",
        default = args$heatmap_first_n,
        help = paste("In heatmap, plot top N features with significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-o", "--plot_scatter"),
        type = "logical",
        dest = "plot_scatter",
        default = args$plot_scatter,
        help = paste("Generate scatter plots for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-g", "--max_pngs"),
        type = "double",
        dest = "max_pngs",
        default = args$max_pngs,
        help = paste("The maximum number of scatterplots for signficant",
                     "associations to save as png files [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-O", "--save_scatter"),
        type = "logical",
        dest = "save_scatter",
        default = args$save_scatter,
        help = paste("Save all scatter plot ggplot objects",
                     "to an RData file [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-e", "--cores"),
        type = "double",
        dest = "cores",
        default = args$cores,
        help = paste("The number of R processes to ",
            "run in parallel [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-j", "--save_models"),
        type = "logical",
        dest = "save_models",
        default = args$save_models,
        help = paste("Return the full model outputs ",
                     "and save to an RData file [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-d", "--reference"),
        type = "character",
        dest = "reference",
        default = args$reference,
        help = paste("The factor to use as a reference for",
            "a variable with more than two levels",
            "provided as a string of 'variable,reference'",
            "semi-colon delimited for multiple variables [ Default: NA ]"
        )
    )

option_not_valid_error <- function(message, valid_options) {
    logging::logerror(paste(message, ": %s"), toString(valid_options))
    stop("Option not valid", call. = FALSE)
}

#######################################################
# Main maaslin2 function (defaults same command line) #
#######################################################

Maaslin2 <-
    function(
        input_data,
        input_metadata,
        output,
        min_abundance = 0.0,
        min_prevalence = 0.1,
        min_variance = 0.0,
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
        heatmap_first_n = 50,
        plot_scatter = TRUE,
        max_pngs = 10,
        save_scatter = FALSE,
        save_models = FALSE,
        reference = NULL)
    {
        # Allow for lower case variables
        normalization <- toupper(normalization)
        transform <- toupper(transform)
        analysis_method <- toupper(analysis_method)
        
        # Match variable ignoring case then set correctly as required for p.adjust
        correction <- correction_choices[match(toupper(correction), toupper(correction_choices))]

        #################################################################
        # Read in the data and metadata, create output folder, init log #
        #################################################################
        # if a character string then this is a file name, else it 
        # is a data frame
        if (is.character(input_data) && file.exists(input_data)) {
            data <-
                data.frame(data.table::fread(
                    input_data, header = TRUE, sep = "\t"),
                        row.names = 1)
            if (nrow(data) == 1) {
                # read again to get row name
                data <- read.table(input_data, header = TRUE, row.names = 1)
            }
        } else if (is.data.frame(input_data)) {
            if (!tibble::has_rownames(input_data)) {
              stop("If supplying input_data as a data frame, it must have appropriate rownames!")
            }
            data <- as.data.frame(input_data) # in case it's a tibble or something
        } else if (is.matrix(input_data)) {
            logging::logwarn("Input is a matrix, passing through as.data.frame() .")
            data <- as.data.frame(input_data)
        } else {
            stop("input_data is neither a file nor a data frame!")
        }

        if (is.character(input_metadata) && file.exists(input_metadata)) {
            metadata <-
                data.frame(data.table::fread(
                    input_metadata, header = TRUE, sep = "\t"),
                        row.names = 1)
            if (nrow(metadata) == 1) {
                metadata <- read.table(input_metadata,
                    header = TRUE,
                    row.names = 1)
            }
        } else if (is.data.frame(input_metadata)) {
            if (!tibble::has_rownames(input_metadata)) {
              stop("If supplying input_metadata as a data frame, it must have appropriate rownames!")
            }
            metadata <- as.data.frame(input_metadata) # in case it's a tibble or something
        } else {
          stop("input_metadata is neither a file nor a data frame!")
        } 
        # create an output folder and figures folder if it does not exist
        if (!file.exists(output)) {
            print("Creating output folder")
            dir.create(output)
        }
        
        features_folder <- file.path(output, "features")
        if (!file.exists(features_folder)) {
            print("Creating output feature tables folder")
            dir.create(features_folder)
        }
        
        fits_folder <- file.path(output, "fits")
        if (!file.exists(fits_folder)) {
            print("Creating output fits folder")
            dir.create(fits_folder)
        }

        if (plot_heatmap || plot_scatter) {
            figures_folder <- file.path(output, "figures")
            if (!file.exists(figures_folder)) {
                print("Creating output figures folder")
                dir.create(figures_folder)
            }
        }
    
        
        # create log file (write info to stdout and debug level to log file)
        # set level to finest so all log levels are reviewed
        log_file <- file.path(output, "maaslin2.log")
        # remove log file if already exists (to avoid append)
        if (file.exists(log_file)) {
            print(paste("Warning: Deleting existing log file:", log_file))
            unlink(log_file)
        }
        logging::basicConfig(level = 'FINEST')
        logging::addHandler(logging::writeToFile, 
            file = log_file, level = "DEBUG")
        logging::setLevel(20, logging::getHandler('basic.stdout'))
        
        #####################
        # Log the arguments #
        #####################
        
        logging::loginfo("Writing function arguments to log file")
        logging::logdebug("Function arguments")
        if (is.character(input_data)) {
            logging::logdebug("Input data file: %s", input_data)
        }
        if (is.character(input_metadata)) {
            logging::logdebug("Input metadata file: %s", input_metadata)
        }
        logging::logdebug("Output folder: %s", output)
        logging::logdebug("Min Abundance: %f", min_abundance)
        logging::logdebug("Min Prevalence: %f", min_prevalence)
        logging::logdebug("Normalization: %s", normalization)
        logging::logdebug("Transform: %s", transform)
        logging::logdebug("Analysis method: %s", analysis_method)
        logging::logdebug("Max significance: %f", max_significance)
        logging::logdebug("Random effects: %s", random_effects)
        logging::logdebug("Fixed effects: %s", fixed_effects)
        logging::logdebug("Correction method: %s", correction)
        logging::logdebug("Standardize: %s", standardize)
        logging::logdebug("Cores: %d", cores)
        
        ####################################
        # Check valid options are selected #
        ####################################
        
        # Check valid normalization option selected
        logging::loginfo("Verifying options selected are valid")
        if (!normalization %in% normalization_choices) {
            option_not_valid_error(
                paste(
                    "Please select a normalization",
                    "from the list of available options"),
                toString(normalization_choices)
            )
        }
        
        # check valid transform option selected
        if (!transform %in% transform_choices) {
            option_not_valid_error(
                "Please select a transform from the list of available options",
                toString(transform_choices)
            )
        }
        
        # check valid method option selected
        if (!analysis_method %in% analysis_method_choices_names) {
            option_not_valid_error(
                paste(
                    "Please select an analysis method",
                    "from the list of available options"),
                toString(analysis_method_choices_names)
            )
        }
        
        # check valid correction method selected
        if (!correction %in% correction_choices) {
            option_not_valid_error(
                paste("Please select a correction method",
                    "from the list of available options"),
                toString(correction_choices)
            )
        }
        
        # check a valid choice combination is selected
        for (limited_method in hash::keys(
            valid_choice_method_norm)) {
            if (analysis_method == limited_method) {
                if (!normalization %in% 
                    valid_choice_method_norm[[limited_method]]) {
                    option_not_valid_error(
                        paste0("This method can only be used ",
                            "with a subset of normalizations. ",
                            "Please select from the following list"
                        ),
                        toString(valid_choice_method_norm[[limited_method]])
                    )
                }
            }
        }
        for (limited_transform in hash::keys(valid_choice_transform_norm)) {
            if (transform == limited_transform) {
                if (!normalization %in% 
                    valid_choice_transform_norm[[limited_transform]]) {
                    option_not_valid_error(
                        paste0("This transform can only be used",
                            " with a subset of normalizations. ",
                            "Please select from the following list"
                        ),
                        toString(
                            valid_choice_transform_norm[[limited_transform]])
                    )
                }
            }
        }
        
        # check that the transform can be applied to the method selected
        if (transform != "NONE")
        {
            if (!analysis_method %in% valid_choice_method_transform) {
                option_not_valid_error(
                    paste0("The transform selected can only be used",
                        " with some methods. ",
                        "Please select from the following list"
                    ),
                    toString(valid_choice_method_transform)
                )
            }
        }
        # check that plots are generated if to be saved
        if (!plot_scatter && save_scatter) {
            logging::logerror("Scatter plots cannot be saved if they are not plotted")
            stop("Option not valid", call. = FALSE)
        }
        
        ###############################################################
        # Determine orientation of data in input and reorder to match #
        ###############################################################
        
        logging::loginfo("Determining format of input files")
        samples_row_row <- intersect(rownames(data), rownames(metadata))
        if (length(samples_row_row) > 0) {
            # this is the expected formatting so do not modify data frames
            logging::loginfo(
                paste(
                    "Input format is data samples",
                    "as rows and metadata samples as rows"))
        } else {
            samples_column_row <- intersect(colnames(data), rownames(metadata))

            if (length(samples_column_row) == 0) {
                # modify possibly included special chars in sample names in metadata
                rownames(metadata) <- make.names(rownames(metadata))
            
                samples_column_row <- intersect(colnames(data), rownames(metadata))
            }

            if (length(samples_column_row) > 0) {
                logging::loginfo(
                    paste(
                        "Input format is data samples",
                        "as columns and metadata samples as rows"))
                # transpose data frame so samples are rows
                data <- as.data.frame(t(data))
                logging::logdebug(
                    "Transformed data so samples are rows")
            } else {
                samples_column_column <- 
                    intersect(colnames(data), colnames(metadata))
                if (length(samples_column_column) > 0) {
                    logging::loginfo(
                        paste(
                            "Input format is data samples",
                            "as columns and metadata samples as columns"))
                    data <- as.data.frame(t(data))
                    metadata <- type.convert(as.data.frame(t(metadata)))
                    logging::logdebug(
                        "Transformed data and metadata so samples are rows")
                } else {
                    samples_row_column <- 
                        intersect(rownames(data), colnames(metadata))

                    if (length(samples_row_column) == 0) {
                        # modify possibly included special chars in sample names in data
                        rownames(data) <- make.names(rownames(data))
            
                        samples_row_column <- intersect(rownames(data), colnames(metadata))
                    }

                    if (length(samples_row_column) > 0) {
                        logging::loginfo(
                            paste(
                                "Input format is data samples",
                                "as rows and metadata samples as columns"))
                        metadata <- type.convert(as.data.frame(t(metadata)))
                        logging::logdebug(
                            "Transformed metadata so samples are rows")
                    } else {
                        logging::logerror(
                            paste("Unable to find samples in data and",
                                "metadata files.",
                                "Rows/columns do not match."))
                        logging::logdebug(
                            "Data rows: %s", 
                            paste(rownames(data), collapse = ","))
                        logging::logdebug(
                            "Data columns: %s", 
                            paste(colnames(data), collapse = ","))
                        logging::logdebug(
                            "Metadata rows: %s", 
                            paste(rownames(metadata), collapse = ","))
                        logging::logdebug(
                            "Metadata columns: %s",
                            paste(colnames(data), collapse = ","))
                        stop()
                    }
                }
            }
        }
       
        # replace unexpected characters in feature names
        colnames(data) <- make.names(colnames(data))
 
        # check for samples without metadata
        extra_feature_samples <-
            setdiff(rownames(data), rownames(metadata))
        if (length(extra_feature_samples) > 0)
            logging::logdebug(
                paste("The following samples were found",
                    "to have features but no metadata.",
                    "They will be removed. %s"),
                paste(extra_feature_samples, collapse = ",")
            )
        
        # check for metadata samples without features
        extra_metadata_samples <-
            setdiff(rownames(metadata), rownames(data))
        if (length(extra_metadata_samples) > 0)
            logging::logdebug(
                paste("The following samples were found",
                    "to have metadata but no features.",
                    "They will be removed. %s"),
                paste(extra_metadata_samples, collapse = ",")
            )
        
        # get a set of the samples with both metadata and features
        intersect_samples <- intersect(rownames(data), rownames(metadata))
        logging::logdebug(
            "A total of %s samples were found in both the data and metadata",
            length(intersect_samples)
        )
        
        # now order both data and metadata with the same sample ordering
        logging::logdebug(
            "Reordering data/metadata to use same sample ordering")
        data <- data[intersect_samples, , drop = FALSE]
        metadata <- metadata[intersect_samples, , drop = FALSE]
        
        ###########################################
        # Compute the formula based on user input #
        ###########################################
        
        random_effects_formula <- NULL
        # use all metadata if no fixed effects are provided
        if (is.null(fixed_effects)) {
            fixed_effects <- colnames(metadata)
        } else {
            fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
            # remove any fixed effects not found in metadata names
            to_remove <- setdiff(fixed_effects, colnames(metadata))
            if (length(to_remove) > 0)
                logging::logwarn(
                    paste("Feature name not found in metadata",
                        "so not applied to formula as fixed effect: %s"),
                    paste(to_remove, collapse = " , ")
                )
            fixed_effects <- setdiff(fixed_effects, to_remove)
            if (length(fixed_effects) == 0) {
                logging::logerror("No fixed effects included in formula.")
                stop()
            }
        }
        
        if (!is.null(random_effects)) {
            random_effects <-
                unlist(strsplit(random_effects, ",", fixed = TRUE))
            # subtract random effects from fixed effects
            common_variables <- intersect(fixed_effects, random_effects)
            if (length(common_variables) > 0) {
                logging::logwarn(
                    paste("Feature name included as fixed and random effect,",
                          "check that this is intended: %s"),
                    paste(common_variables, collapse = " , ")
                )
            }
            # remove any random effects not found in metadata
            to_remove <- setdiff(random_effects, colnames(metadata))
            if (length(to_remove) > 0)
                logging::logwarn(
                    paste("Feature name not found in metadata",
                        "so not applied to formula as random effect: %s"),
                    paste(to_remove, collapse = " , ")
                )
            random_effects <- setdiff(random_effects, to_remove)
            
            # create formula
            if (length(random_effects) > 0) {
                random_effects_formula_text <-
                    paste(
                        "expr ~ (1 | ",
                        paste(
                            random_effects,
                            ")",
                            sep = '',
                            collapse = " + (1 | "
                        ),
                        sep = '')
                logging::loginfo("Formula for random effects: %s",
                    random_effects_formula_text)
                random_effects_formula <-
                    tryCatch(
                        as.formula(random_effects_formula_text),
                        error = function(e)
                            stop(
                                paste(
                                    "Invalid formula for random effects: ",
                                    random_effects_formula_text
                                )
                            )
                    )
            }
        }
        
        # reduce metadata to only include fixed/random effects in formula
        effects_names <- union(fixed_effects, random_effects)
        metadata <- metadata[, effects_names, drop = FALSE]
        
        # create the fixed effects formula text
        formula_text <-
            paste("expr ~ ", paste(fixed_effects, collapse = " + "))
        logging::loginfo("Formula for fixed effects: %s", formula_text)
        formula <-
            tryCatch(
                as.formula(formula_text),
                error = function(e)
                    stop(
                        paste(
                            "Invalid formula.",
                            "Please provide a different formula: ",
                            formula_text
                        )
                    )
            )
        
        #########################################################
        # Filter data based on min abundance and min prevalence #
        #########################################################

        if (is.null(reference)) {
            reference <- ","
        }
        split_reference <- unlist(strsplit(reference, "[,;]"))
        
        # for each fixed effect, check that a reference level has been set if necessary: number of levels > 2 and metadata isn't already an ordered factor
        for (i in fixed_effects) {
            # don't check for or require reference levels for numeric metadata
            if (is.numeric(metadata[,i])) {
                next
            }
            # respect ordering if a factor is explicitly passed in with no reference set
            if (is.factor(metadata[,i]) && !(i %in% split_reference)) {
                logging::loginfo(paste("Factor detected for categorial metadata '", 
                                       i, "'. Provide a reference argument or manually set factor ordering to change reference level.", sep=""))
                next
            }
            
            # set metadata as a factor (ordered alphabetically)
            metadata[,i] <- as.factor(metadata[,i])
            mlevels <- levels(metadata[,i])
            
            # get reference level for variable being considered, returns NA if not found
            ref <- split_reference[match(i, split_reference)+1]
            
            # if metadata has 2 levels, allow but don't require setting reference level, otherwise require it
            if ((length(mlevels) == 2)) {
                if(!is.na(ref)) {
                    metadata[,i] = relevel(metadata[,i], ref = ref)
                }
            } else if (length(mlevels) > 2) {
                if (!is.na(ref)) {
                    metadata[,i] = relevel(metadata[,i], ref = ref)
                } else {
                    stop(paste("Please provide the reference for the variable '",
                               i, "' which includes more than 2 levels: ",
                               paste(as.character(mlevels), collapse=", "), ".", sep=""))   
                } 
            } else {
                stop("Provided categorical metadata has fewer than 2 unique, non-NA values.")
            }
        }
 
        unfiltered_data <- data
        unfiltered_metadata <- metadata
        
        # require at least total samples * min prevalence values 
        # for each feature to be greater than min abundance
        logging::loginfo(
            "Filter data based on min abundance and min prevalence")
        total_samples <- nrow(unfiltered_data)
        logging::loginfo("Total samples in data: %d", total_samples)
        min_samples <- total_samples * min_prevalence
        logging::loginfo(
            paste("Min samples required with min abundance",
                "for a feature not to be filtered: %f"),
            min_samples
        )
        
        # Filter by abundance using zero as value for NAs
        data_zeros <- unfiltered_data
        data_zeros[is.na(data_zeros)] <- 0
        filtered_data <-
            unfiltered_data[, 
                colSums(data_zeros > min_abundance) > min_samples,
                drop = FALSE]
        total_filtered_features <-
            ncol(unfiltered_data) - ncol(filtered_data)
        logging::loginfo("Total filtered features: %d", total_filtered_features)
        filtered_feature_names <-
            setdiff(names(unfiltered_data), names(filtered_data))
        logging::loginfo("Filtered feature names from abundance and prevalence filtering: %s",
            toString(filtered_feature_names))
        
        #################################
        # Filter data based on variance #
        #################################
        
        sds <- apply(filtered_data, 2, sd)
        variance_filtered_data <- filtered_data[, which(sds > min_variance), drop = FALSE]
        variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
        logging::loginfo("Total filtered features with variance filtering: %d", variance_filtered_features)
        variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
        logging::loginfo("Filtered feature names from variance filtering: %s",
                         toString(variance_filtered_feature_names))
        filtered_data <- variance_filtered_data
       
        ######################
        # Normalize features #
        ######################
        
        logging::loginfo(
            "Running selected normalization method: %s", normalization)
        filtered_data_norm <-
            normalizeFeatures(filtered_data, normalization = normalization)
        
        ################################
        # Standardize metadata, if set #
        ################################
        
        if (standardize) {
            logging::loginfo(
                "Applying z-score to standardize continuous metadata")
            metadata <- dplyr::mutate_if(metadata, is.numeric, scale)
        } else {
            logging::loginfo("Bypass z-score application to metadata")
        }
        
        ############################
        # Transform and run method #
        ############################
       
        # transform features
        logging::loginfo("Running selected transform method: %s", transform)
        filtered_data_norm_transformed <-
            transformFeatures(filtered_data_norm, transformation = transform)
        
        # apply the method to the data with the correction
        logging::loginfo(
            "Running selected analysis method: %s", analysis_method)
        fit_data <-
            fit.data(
                filtered_data_norm_transformed,
                metadata,
                analysis_method,
                formula = formula,
                random_effects_formula = random_effects_formula,
                correction = correction,
                save_models = save_models,
                cores = cores
            )
        
        #################################################################
        # Count the total values for each feature (untransformed space) #
        #################################################################
        
        logging::loginfo("Counting total values for each feature")
        fit_data$results$N <-
            apply(
                fit_data$results,
                1,
                FUN = function(x)
                    length(filtered_data_norm[, x[1]])
            )
        fit_data$results$N.not.zero <-
            apply(
                fit_data$results,
                1,
                FUN = function(x)
                    length(which(filtered_data_norm[, x[1]] > 0))
            )
        
        ################################
        # Write out the raw model fits #
        ################################
        
        if (save_models) {
            model_file = file.path(fits_folder, "models.rds")
            # remove models file if already exists (since models append)
            if (file.exists(model_file)) {
                logging::logwarn(
                    "Deleting existing model objects file: %s", model_file)
                unlink(model_file)
            }
            logging::loginfo("Writing model objects to file %s", model_file)
            saveRDS(fit_data$fits, file = model_file)   
        }
        
        ##########################################
        # Write processed feature tables to file #
        ##########################################
        
        filtered_file = file.path(features_folder, "filtered_data.tsv")
        logging::loginfo("Writing filtered data to file %s", filtered_file)
        write.table(
            data.frame("feature" = rownames(filtered_data), filtered_data), 
            file = filtered_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
            )
        
        filtered_data_norm_file = file.path(features_folder, "filtered_data_norm.tsv")
        logging::loginfo("Writing filtered, normalized data to file %s", filtered_data_norm_file)
        write.table(
            data.frame("feature" = rownames(filtered_data_norm), filtered_data_norm), 
            file = filtered_data_norm_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
        )
        
        filtered_data_norm_transformed_file = file.path(features_folder, "filtered_data_norm_transformed.tsv")
        logging::loginfo("Writing filtered, normalized, transformed data to file %s", filtered_data_norm_transformed_file)
        write.table(
            data.frame("feature" = rownames(filtered_data_norm_transformed), filtered_data_norm_transformed), 
            file = filtered_data_norm_transformed_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
        )
        
        ###########################
        # Write residuals to file #
        ###########################
        
        residuals_file = file.path(fits_folder, "residuals.rds")
        # remove residuals file if already exists (since residuals append)
        if (file.exists(residuals_file)) {
            logging::logwarn(
                "Deleting existing residuals file: %s", residuals_file)
            unlink(residuals_file)
        }
        logging::loginfo("Writing residuals to file %s", residuals_file)
        saveRDS(fit_data$residuals, file = residuals_file)
        
        ###############################
        # Write fitted values to file #
        ###############################
        
        fitted_file = file.path(fits_folder, "fitted.rds")
        # remove fitted file if already exists (since fitted append)
        if (file.exists(fitted_file)) {
          logging::logwarn(
            "Deleting existing fitted file: %s", fitted_file)
          unlink(fitted_file)
        }
        logging::loginfo("Writing fitted values to file %s", fitted_file)
        saveRDS(fit_data$fitted, file = fitted_file)
        
        #########################################################
        # Write extracted random effects to file (if specified) #
        #########################################################
        
        if (!is.null(random_effects)) {
          ranef_file = file.path(fits_folder, "ranef.rds")
          # remove ranef file if already exists (since ranef append)
          if (file.exists(ranef_file)) {
            logging::logwarn(
              "Deleting existing ranef file: %s", ranef_file)
            unlink(ranef_file)
          }
          logging::loginfo("Writing extracted random effects to file %s", ranef_file)
          saveRDS(fit_data$ranef, file = ranef_file)
        }
        
        #############################
        # Write all results to file #
        #############################
        
        results_file <- file.path(output, "all_results.tsv")
        logging::loginfo(
            "Writing all results to file (ordered by increasing q-values): %s",
            results_file)
        ordered_results <- fit_data$results[order(fit_data$results$qval), ]
        # Remove any that are NA for the q-value
        ordered_results <-
            ordered_results[!is.na(ordered_results$qval), ]
        write.table(
            ordered_results[c(
                "feature",
                "metadata",
                "value",
                "coef",
                "stderr",
                "N",
                "N.not.zero",
                "pval",
                "qval")],
            file = results_file,
            sep = "\t",
            quote = FALSE,
            col.names = c(
                "feature",
                "metadata",
                "value",
                "coef",
                "stderr",
                "N",
                "N.not.0",
                "pval",
                "qval"
            ),
            row.names = FALSE
        )
        
        ###########################################
        # Write results passing threshold to file #
        ###########################################
        
        significant_results <-
            ordered_results[ordered_results$qval <= max_significance, ]
        significant_results_file <-
            file.path(output, "significant_results.tsv")
        logging::loginfo(
            paste("Writing the significant results",
                "(those which are less than or equal to the threshold",
                "of %f ) to file (ordered by increasing q-values): %s"),
            max_significance,
            significant_results_file
        )
        write.table(
            significant_results[c(
                "feature",
                "metadata",
                "value",
                "coef",
                "stderr",
                "N",
                "N.not.zero",
                "pval",
                "qval")],
            file = significant_results_file,
            sep = "\t",
            quote = FALSE,
            col.names = c(
                "feature",
                "metadata",
                "value",
                "coef",
                "stderr",
                "N",
                "N.not.0",
                "pval",
                "qval"
            ),
            row.names = FALSE
        )
        
        #######################################################
        # Create visualizations for results passing threshold #
        #######################################################
        
        if (plot_heatmap) {
            heatmap_file <- file.path(output, "heatmap.pdf")
            logging::loginfo(
                "Writing heatmap of significant results to file: %s",
                heatmap_file)
            save_heatmap(significant_results_file, heatmap_file, figures_folder,
                first_n = heatmap_first_n)
        }
        
        if (plot_scatter) {
            logging::loginfo(
                paste("Writing association plots",
                    "(one for each significant association)",
                    "to output folder: %s"),
                output
            )
            saved_plots <- maaslin2_association_plots(
                unfiltered_metadata,
                filtered_data,
                significant_results_file,
                output,
                figures_folder,
                max_pngs,
                save_scatter)
            if (save_scatter) {
                scatter_file <- file.path(figures_folder, "scatter_plots.rds")
                # remove plots file if already exists
                if (file.exists(scatter_file)) {
                    logging::logwarn(
                        "Deleting existing scatter plot objects file: %s", scatter_file)
                    unlink(scatter_file)
                }
                logging::loginfo("Writing scatter plot objects to file %s", scatter_file)
                saveRDS(saved_plots, file = scatter_file)   
            }
        }
        
        return(fit_data)
    }

###########################################################################
# If running on the command line, get arguments and call maaslin function #
###########################################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
        !length(grep("^source\\(", sys.calls()))) {
    # get command line options and positional arguments
    parsed_arguments = optparse::parse_args(options, 
        positional_arguments = TRUE)
    current_args <- parsed_arguments[["options"]]
    positional_args <- parsed_arguments[["args"]]
    # check three positional arguments are provided
    if (length(positional_args) != 3) {
        optparse::print_help(options)
        stop(
            paste("Please provide the required",
                "positional arguments",
                "<data.tsv> <metadata.tsv> <output_folder>")
        )
    }
    
    # call maaslin with the command line options
    fit_data <-
        Maaslin2(
            positional_args[1],
            positional_args[2],
            positional_args[3],
            current_args$min_abundance,
            current_args$min_prevalence,
            current_args$min_variance,
            current_args$normalization,
            current_args$transform,
            current_args$analysis_method,
            current_args$max_significance,
            current_args$random_effects,
            current_args$fixed_effects,
            current_args$correction,
            current_args$standardize,
            current_args$cores,
            current_args$plot_heatmap,
            current_args$heatmap_first_n,
            current_args$plot_scatter,
            current_args$max_pngs,
            current_args$save_scatter,
            current_args$save_models,
            current_args$reference
        )
}
