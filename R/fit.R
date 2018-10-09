# Load Required Packages
for( lib in c('dplyr', 'pbapply', 'MASS', 'lme4', 'car', 'cplm', 'nlme', 'pscl', 'parallel')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}

# fit the data using the model selected and applying the correction
fit.data <- function(features, metadata, model, formula = NULL, random_effects_formula = NULL, correction = "BH", cores = 1){

  # set the formula default to all fixed effects if not provided  
  if (is.null(formula)) formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))

  #############################################################
  # Determine the function and summary for the model selected #
  #############################################################

  if (model=="LM") {
    if (is.null(random_effects_formula)) {
      model_function <- function(formula, data, na.action) { return(glm(formula, data = data, family='gaussian', na.action = na.action)) }
      summary_function <- function(fit) {
        lm_summary <- summary(fit)$coefficients
        para<-as.data.frame(lm_summary)[-1,-c(2:3)]
        para$name<-rownames(lm_summary)[-1]
        return(para)     
      }
    } else {
      model_function <- function(formula, data, na.action) { return(nlme::lme(fixed=formula, random=random_effects_formula, data = data, na.action = na.action)) }
      summary_function <- function(fit) {
        lm_summary<-coef(summary(fit))
        para<-as.data.frame(lm_summary)[-1,-c(2:4)]
        para$name<-rownames(lm_summary)[-1]
        return(para)
      }
    }
  }

  if (model=="CPLM") {
    model_function <- cplm::cpglm
    summary_function <- function(fit) {
      cplm_out<-capture.output( cplm_summary <- cplm::summary(fit)$coefficients )
      para<-as.data.frame(cplm_summary)[-1,-c(2:3)]
      para$name<-rownames(cplm_summary)[-1]
      logging::logdebug("Summary output\n%s", paste(cplm_out, collapse="\n"))
      return(para)
    }
  }

  if (model=="NEGBIN") {
    model_function <- MASS::glm.nb
    summary_function <- function(fit) {
      glm_summary <- summary(fit)$coefficients
      para<-as.data.frame(glm_summary)[-1,-c(2:3)]
      para$name<-rownames(glm_summary)[-1]
      return(para)
    }
  }

  if (model=="ZICP") {
    model_function <- cplm::zcpglm
    summary_function <- function(fit) {
      cplm_out<-capture.output( cplm_summary <- summary(fit)$coefficients$tweedie )
      para<-as.data.frame(cplm_summary)[-1,-c(2:3)]
      para$name<-rownames(cplm_summary)[-1]
      logging::logdebug("Summary output\n%s", paste(cplm_out, collapse="\n"))
      return(para)
    }
  }

  if (model=="ZINB") {
    model_function <- pscl::zeroinfl
    summary_function <- function(fit) {
      pscl_summary <- summary(fit)$coefficients$count
      para<-as.data.frame(pscl_summary)[-c(1, (ncol(metadata)+2)),-c(2:3)]
      para$name<-rownames(pscl_summary)[c(2:11)]
      return(para)
    }
  }
  
  #######################################
  # Init cluster for parallel computing #
  #######################################

  cluster <- NULL
  if (cores > 1)
  {  
      logging::loginfo("Creating cluster of %s R processes", cores)
      cluster <- parallel::makeCluster(cores)
  }

  ############################## 
  # Apply per-feature modeling #
  ##############################
  outputs <- pbapply::pblapply(1:ncol(features), cl=cluster, function(x){
    
    # Extract Features One by One
    featuresVector <- features[, x]
    
    # Fit Model
    logging::loginfo("Fitting model to feature number %d, %s", x, colnames(features)[x])
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
    fit <- tryCatch({
      fit1 <- model_function(formula, data = dat_sub, na.action = na.exclude)
    }, error=function(err){
      fit1 <- try({model_function(formula, data = dat_sub, na.action = na.exclude)}) 
      return(fit1)
    })
    
    # Gather Output
    output<-list()
    if (all(class(fit) != "try-error")){
          output$para<-summary_function(fit)
          output$residuals<-residuals(fit)
          } 
    else{
          logging::logwarn(paste("Fitting problem for feature", x, "returning NA"))
          output$para<-as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
          output$para$name<-colnames(metadata)
          output$residuals<-NA
        }
    colnames(output$para)<-c('coef', 'pval', 'name')
    output$para$feature<-colnames(features)[x]
    return(output)
  })    

  # stop the cluster
  if (! is.null(cluster) ) parallel::stopCluster(cluster)

  # bind the results for each feature
  paras<-do.call(rbind, lapply(outputs, function(x) { return(x$para) }))
  residuals<-do.call(rbind, lapply(outputs, function(x) { return(x$residuals) }))
  row.names(residuals)<-colnames(features)

  ################################
  # Apply correction to p-values #
  ################################

  paras$qval<-as.numeric(p.adjust(paras$pval, method = correction))

  #####################################################
  # Determine the metadata names from the model names #
  #####################################################

  metadata_names<-colnames(metadata)
  # order the metadata names by decreasing length
  metadata_names_ordered<-metadata_names[order(nchar(metadata_names),decreasing=TRUE)]
  # find the metadata name based on the match to the beginning of the string
  extract_metadata_name <- function(name) {
    return(metadata_names_ordered[mapply(startsWith,name,metadata_names_ordered)][1])
  }
  paras$metadata<-unlist(lapply(paras$name,extract_metadata_name))
  # compute the value as the model contrast minus metadata
  paras$value<-mapply(function(x,y) { if (x==y) x else gsub(x,"",y) }, paras$metadata, paras$name)

  ##############################
  # Sort by decreasing q-value #
  ##############################

  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata', 'value'), dplyr::everything())
  return(list("results"=paras,"residuals"=residuals))  
}

