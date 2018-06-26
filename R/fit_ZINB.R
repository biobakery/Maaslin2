# Install or Load Required Packages

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org') 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'pbapply', 'pscl')

# Fit Zero-inflated Negative Binomial (ZINB) To A Dataset

fit.ZINB <- function(features, 
                   metadata, 
                   normalization ='NONE', 
                   transformation ='NONE', 
                   randomEffect = FALSE){
  
  #######################################
  # Apply Normalization to the Features #
  #######################################
  
  features<-normalizeFeatures(features, normalization = normalization)
  
  ######################################################################
  # Apply Per-Feature Modeling Followed by User-defined Transformation #
  ######################################################################
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    featuresVector <- features[, x]
    
    # Fit Model
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    fit <- tryCatch({
      fit1 <- pscl::zeroinfl(formula, data = dat_sub, dist = "negbin")
    }, error=function(err){
      fit1 <- try({pscl::zeroinfl(formula, data = dat_sub)}) 
      return(fit1)
    })
    
    # Gather Output
    if (class(fit) != "try-error"){
          para<-as.data.frame(summary(fit)$coefficients$count)[-c(1, (ncol(metadata)+2)),-c(2:3)]
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
          para$feature<-colnames(features)[x]
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

