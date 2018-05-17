# Install or Load Required Packages

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('pbapply', 'car', 'nlme', 'dplyr')

# Fit Linear Model To A Dataset

fit.LM <- function(features, 
                   metadata, 
                   normalization='TSS', 
                   transformation='LOG', 
                   randomEffect=FALSE){
  
  
  # Basic Error Messages for Illegal Combinations 
  if (!transformation %in% c('AST', 'LOG', 'LOGIT', 'NONE')) {
    stop ('Transformation should be one of AST/LOG/LOGIT/NONE.')
  }
  
  if (!normalization %in% c('TSS', 'CSS', 'TMM', 'CLR')) {
    stop ('Normalization should be one of TSS/CSS/TMM/CLR.')
  }
  
  if (normalization!='TSS' && transformation %in% c('AST', 'LOG', 'LOGIT')) {
    stop ('This combination of normalization and transformation is invalid.')
  }

  #######################################
  # Apply Normalization to the Features #
  #######################################
  
  if (normalization=='TSS')
  {
    features<-TSSnorm(features)
  }
  
  if (normalization=='CLR')
  {
    features<-CLRnorm(features)
  }
  
  ######################################################################
  # Apply Per-Feature Modeling Followed by User-defined Transformation #
  ######################################################################

  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    featuresVector <- features[, x]
    
    # Transform
    if (transformation =='LOG') featuresVector<-LOG(featuresVector);
    if (transformation =='LOGIT') featuresVector<-LOGIT(featuresVector);
    if (transformation =='AST') featuresVector<-AST(featuresVector);
    
    # Fit Model
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    fit <- tryCatch({
          fit1 <- glm(formula, data = dat_sub, family='gaussian')
        }, error=function(err){
          fit1 <- try({glm(formula, data = dat_sub, family='gaussian')}) 
          return(fit1)
        })
    
    # Gather Output
    if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
          } 
    else{
          print(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          colnames(para)<-c('coef', 'pval')
          para$metadata<-colnames(metadata)
          para$features<-colnames(features)[x]
          rownames(para)<-NULL
        }
    return(para)
  })    
   
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = "fdr"))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  return(paras)  
}

