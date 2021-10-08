# Load Required Packages
for (lib in c('vegan', 'chemometrics', 'car', 'metagenomeSeq', 'edgeR')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

###################
##   Filtering   ##
###################
do_prevalence_abundance_filtering <- function(
    unfiltered_data,
    min_abundance = 0.0,
    min_prevalence = 0.1,
    min_variance = 0.0){
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
    use_raw_abund_treshold = TRUE
    # the example hmp dataset has relative abundance sample sums up to 1.0000004
    if (all(rowSums(data_zeros) <=1.01)){
        logging::loginfo("relative abundances detected")
        if(min_abundance >1){
            logging::loginfo("Min_abundance > 1 not informative for relative abundance data; converting min_abundance to percentage")
            min_abundance = min_abundance/100
            logging::loginfo("Min_abundance threshold: %s, or %s%%", min_abundance, 100*min_abundance)
            use_raw_abund_treshold = FALSE
        }
    } else{
        logging::loginfo("Count (or normalized) data detected")
        if(min_abundance < 1 & all(data_zeros[data_zeros !=0 ] > min_abundance) ){
            if(min_abundance != 0){
                logging::loginfo("Min_abundance is less than 1  and no non-zero values in the data are less than 1. Treating minimum abundance as a percentage. To avoid this behavior, set min_abundance above the lowest non-zero value in data")
                use_raw_abund_treshold = FALSE
            }
        }
    }
    if(use_raw_abund_treshold){
        filtered_data <-
            unfiltered_data[,
                            colSums(data_zeros > min_abundance) > min_samples,
                            drop = FALSE]
        perc_lost_abundance <- rowSums(data_zeros <= min_abundance)/ncol(data_zeros)
    } else{
        # do filtering on prop table, carry names back to counts table
        # margin 1 means do proportions rowwise
        data_zeros_prop <- prop.table(data.matrix(data_zeros), margin = 1)
        filtered_data_prop <-
            unfiltered_data[,
                            colSums(data_zeros_prop > min_abundance ) > min_samples,
                            drop = FALSE]
        filtered_data <- unfiltered_data[
            rownames(filtered_data_prop),
            colnames(filtered_data_prop)]
        perc_lost_abundance <- rowSums(data_zeros_prop <= min_abundance)/ncol(data_zeros_prop)
    }
    total_filtered_features <-
        ncol(unfiltered_data) - ncol(filtered_data)
    logging::loginfo("Total filtered features: %d", total_filtered_features)
    filtered_feature_names <-
        setdiff(names(unfiltered_data), names(filtered_data))
    logging::loginfo("Mean percent of features removed from abundance filtering: %s",
                     toString(mean(perc_lost_abundance)))
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
    
    return(filtered_data)
}
###################
## Transformation #
###################

transformFeatures = function(features, transformation) {
    if (transformation == 'LOG')     {
        features <- apply(features, 2, LOG)
    }
    
    if (transformation == 'LOGIT')     {
        features <- apply(features, 2, LOGIT)
    }
    
    if (transformation == 'AST')     {
        features <- apply(features, 2, AST)
    }
    
    return(features)
}


##################
## Normalization #
##################

normalizeFeatures = function(features, normalization) {
    if (normalization == 'TSS')
    {
        features <- TSSnorm(features)
    }
    
    if (normalization == 'CLR')
    {
        features <- CLRnorm(features)
    }
    
    if (normalization == 'CSS')
    {
        features <- CSSnorm(features)
    }
    
    if (normalization == 'TMM')
    {
        features <- TMMnorm(features)
    }
    
    if (normalization == 'NONE')
    {
        features <- features
    }
    
    return(features)
}

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # TSS Normalizing the Data
    features_TSS <-
        vegan::decostand(
            features_norm,
            method = "total",
            MARGIN = 1,
            na.rm = TRUE)
    
    # Convert back to data frame
    features_TSS <- as.data.frame(features_TSS)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_TSS) <- dd
    
    
    # Return
    return(features_TSS)
}


######################
## CLR Normalization #
######################

# Apply CLR Normalization To A Dataset

CLRnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # CLR Normalizing the Data
    features_CLR <- chemometrics::clr(features_norm + 1)
    
    # Convert back to data frame
    features_CLR <- as.data.frame(features_CLR)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_CLR) <- dd
    
    
    # Return
    return(features_CLR)
}

######################
## CSS Normalization #
######################

# Apply CSS Normalization To A Dataset

CSSnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # CSS Normalizing the Data
    # Create the metagenomeSeq object
    MGS = metagenomeSeq::newMRexperiment(
        t(features_norm),
        featureData = NULL,
        libSize = NULL,
        normFactors = NULL
    )
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = metagenomeSeq::cumNorm(MGS, p = metagenomeSeq::cumNormStat(MGS))
    # Save the normalized data as data.frame
    features_CSS = as.data.frame(t(
        metagenomeSeq::MRcounts(MGS, norm = TRUE, log = FALSE)))
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_CSS) <- dd
    
    
    # Return as list
    return(features_CSS)
}

######################
## TMM Normalization #
######################

# Apply TMM Normalization To A Dataset

TMMnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # TMM Normalizing the Data
    X <- t(features_norm)
    
    libSize = edgeR::calcNormFactors(X, method = "TMM")
    eff.lib.size = colSums(X) * libSize
    
    ref.lib.size = mean(eff.lib.size)
    #Use the mean of the effective library sizes as a reference library size
    X.output = sweep(X, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
    #Normalized read counts
    
    # Convert back to data frame
    features_TMM <- as.data.frame(t(X.output))
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_TMM) <- dd
    
    
    # Return as list
    return(features_TMM)
}

#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST <- function(x) {
    return(sign(x) * asin(sqrt(abs(x))))
}

########################
# Logit Transformation #
########################

# Zero-inflated Logit Transformation (Does not work well for microbiome data)
LOGIT <- function(x) {
    y <- car::logit(x, adjust = 0)
    y[!is.finite(y)] <- 0
    return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#     y<-0.5*log(x/(1-x)) + 10
#     y[!is.finite(y)]<-0
#     return(y)
# }

######################
# Log Transformation #
######################

# Log Transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}
