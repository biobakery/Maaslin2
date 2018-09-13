# Load Required Packages
for( lib in c('pbapply', 'car', 'nlme', 'dplyr')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}

# Fit Linear Model To A Dataset

fit.LM <- function(features, 
                   metadata, 
                   normalization ='TSS', 
                   transformation ='LOG', 
                   randomEffect = NULL,
                   formula = NULL){
  
  #######################################
  # Apply Normalization to the Features #
  #######################################
  
  features<-normalizeFeatures(features, normalization = normalization)
  
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
    if (is.null(formula)) formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    if (!is.null(randomEffect)) {
        fit <- tryCatch({
              fit1 <- nlme::lme(fixed=formula, random=randomEffect[["formula"]], data = dat_sub)
            }, error=function(err){
              fit1 <- try({nlme::lme(fixed=formula, random=randomEffect[["formula"]], data = dat_sub)}) 
              return(fit1)
            })
        # Gather Output
        if (all(class(fit) != "try-error")){
              para<-as.data.frame(coef(summary(fit)))[-1,-c(2:4)]
              } 
        else{
              print(paste("Fitting problem for feature", x, "returning NA"))
              para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
            }
        colnames(para)<-c('coef', 'pval')
        para$metadata<-setdiff(colnames(metadata),randomEffect[["names"]])
    } 
    else {
        fit <- tryCatch({
              fit1 <- glm(formula, data = dat_sub, family='gaussian')
            }, error=function(err){
              fit1 <- try({glm(formula, data = dat_sub, family='gaussian')}) 
              return(fit1)
            })
        # Gather Output
        if (all(class(fit) != "try-error")){
              para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
              } 
        else{
              print(paste("Fitting problem for feature", x, "returning NA"))
              para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
            }
        colnames(para)<-c('coef', 'pval')
        para$metadata<-colnames(metadata)
    }

    para$feature<-colnames(features)[x]
    rownames(para)<-NULL

    return(para)
  })    
   
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = "fdr"))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), dplyr::everything())
  rownames(paras)<-NULL
  return(paras)  
}

