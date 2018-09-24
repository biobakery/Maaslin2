# Load Required Packages
for( lib in c('car', 'dplyr', 'pbapply', 'cplm')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}


# Fit Compound Poisson Linear Model (CPLM) To A Dataset

fit.CPLM <- function(features, 
                   metadata, 
                   normalization ='TSS', 
                   transformation ='LOG', 
                   formula = NULL,
                   correction = "BH",
                   residuals_file = NONE){
  
  ##########################################################
  # Apply Normalization and Transformation to the Features #
  ##########################################################
  
  features<-normalizeFeatures(features, normalization = normalization)
  features<-transformFeatures(features, transformation = transformation)  
  
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
      fit1 <- cplm::cpglm(formula, data = dat_sub)
    }, error=function(err){
      fit1 <- try({cplm::cpglm(formula, data = dat_sub)}) 
      return(fit1)
    })
    
    # Gather Output
    if (all(class(fit) != "try-error")){
          cplm_out<-capture.output( cplm_summary <- cplm::summary(fit)$coefficients )
          para<-as.data.frame(cplm_summary)[-1,-c(2:3)]
          para$name<-rownames(cplm_summary)[-1]
          logging::logdebug("Summary output\n%s", paste(cplm_out, collapse="\n"))
          if (!is.null(residuals_file)) write(paste("Residuals for feature",x,paste(residuals(fit), collapse=",")),file=residuals_file,append=TRUE)
          colnames(para)<-c('coef', 'pval', 'name')
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
          } 
    else{
          logging::logwarn(paste("Fitting problem for feature", x, "returning NA"))
          para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          para$name<-colnames(metadata)
          colnames(para)<-c('coef', 'pval', 'name')
          para$feature<-colnames(features)[x]
          rownames(para)<-NULL
        }
    return(para)
  })    
   
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = correction))
  # determine the metadata names from the model names
  metadata_names<-colnames(metadata)
  # order the metadata names by decreasing length
  metadata_names_ordered<-metadata_names[order(nchar(metadata_names),decreasing=TRUE)]
  # find the metadata name based on the match to the beginning of the string
  extract_metadata_name <- function(name) {
    return(metadata_names_ordered[mapply(startsWith,name,metadata_names_ordered)][1])
  }
  paras$metadata<-unlist(lapply(paras$name,extract_metadata_name))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata', 'name'), dplyr::everything())
  return(paras)  
}

