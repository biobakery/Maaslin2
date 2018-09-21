# Load Required Packages
for( lib in c('dplyr', 'pbapply', 'pscl')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}

# Fit Zero-inflated Negative Binomial (ZINB) To A Dataset

fit.ZINB <- function(features, 
                   metadata, 
                   normalization ='NONE', 
                   transformation ='NONE', 
                   formula = NULL,
                   correction = "BH",
                   residuals_file = NULL){
  
  ##########################################################
  # Apply Normalization and Transformation to the Features #
  ##########################################################
  
  features<-normalizeFeatures(features, normalization = normalization)
  
  ######################################################################
  # Apply Per-Feature Modeling Followed by User-defined Transformation #
  ######################################################################
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    # Extract Features One by One
    featuresVector <- features[, x]
    
    # Fit Model
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    if (is.null(formula)) formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    fit <- tryCatch({
      fit1 <- pscl::zeroinfl(formula, data = dat_sub, dist = "negbin")
    }, error=function(err){
      fit1 <- try({pscl::zeroinfl(formula, data = dat_sub)}) 
      return(fit1)
    })
    
    # Gather Output
    if (all(class(fit) != "try-error")){
          pscl_summary <- summary(fit)$coefficients$count
          para<-as.data.frame(pscl_summary)[-c(1, (ncol(metadata)+2)),-c(2:3)]
          para$name<-rownames(pscl_summary)[c(2:11)]
          if (!is.null(residuals_file)) write(paste("Residuals for feature",x,paste(residuals(fit), collapse=",")),file=residuals_file,append=TRUE)
          colnames(para)<-c('coef', 'pval', 'name')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
          } 
    else{
          logging::logwarn(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          para$name<-colnames(metadata)
          colnames(para)<-c('coef', 'pval', 'name')
          para$metadata<-colnames(metadata)
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
    return(para)
  })    
   
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = correction))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata', 'name'), dplyr::everything())
  return(paras)  
}

