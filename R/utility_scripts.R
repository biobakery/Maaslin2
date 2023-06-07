# Load Required Packages
for (lib in c('vegan', 'chemometrics', 'car', 'metagenomeSeq', 'edgeR')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
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
    y <- sign(x) * asin(sqrt(abs(x)))
    if(any(is.na(y))) {
        logging::logerror(
            paste0("AST transform is only valid for values between -1 and 1. ",
                   "Please select an appropriate normalization option or ",
                   "normalize your data prior to running."))
        stop()
    }
    return(y)
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
