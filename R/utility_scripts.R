# Install or Load Required Packages

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('vegan', 'chemometrics', 'car')

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(features) {
  
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd<-colnames(features_norm)
  
  # TSS Normalizing the Data
  features_TSS <- vegan::decostand(features_norm, method="total", MARGIN=1)
  
  # Convert back to data frame
  features_TSS<-as.data.frame(features_TSS)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TSS) <- dd;
  
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
  dd<-colnames(features_norm)
  
  # CLR Normalizing the Data
  features_CLR<-chemometrics::clr(features_norm+1)
  
  # Convert back to data frame
  features_CLR<-as.data.frame(features_CLR)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_CLR) <- dd;
  
  # Return
  return(features_CLR)
}

#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST<-function(x){
  return(sign(x)*asin(sqrt(abs(x))))
}

########################
# Logit Transformation #
########################

# Zero-inflated Logit Transformation (Does not work well for microbiome data)
LOGIT<-function(x){
  y<-car::logit(x, adjust=0)
  y[!is.finite(y)]<-0
  return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#   y<-0.5*log(x/(1-x)) + 10
#   y[!is.finite(y)]<-0
#   return(y)
# }

######################
# Log Transformation #
######################

# Log Transformation
LOG<-function(x){
  return(log(x+1))
}
