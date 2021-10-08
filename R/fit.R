# Load Required Packages
for (lib in c(
    'dplyr',
    'pbapply',
    'lmerTest',
    'car',
    'parallel',
    'MuMIn',
    'glmmTMB',
    'MASS',
    'cplm',
    'pscl'
)) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# fit the data using the model selected and applying the correction
fit.data <-
    function(
        features,
        metadata,
        model,
        formula = NULL,
        random_effects_formula = NULL,
        correction = "BH",
        cores = 1,
        save_fits_to=NA) {

        # set the formula default to all fixed effects if not provided
        if (is.null(formula))
            formula <-
                as.formula(paste(
                    "expr ~ ", 
                    paste(colnames(metadata), 
                    collapse = "+")))

        if (!(is.null(random_effects_formula))) {
            formula <-
                paste(
                    '. ~', 
                    #paste(all.vars(formula)[-1], collapse = ' + '), 
                    paste(labels(terms(tmp)), collapse = ' + '), 
                    '.', 
                    sep = ' + ')
	    formula <- update(random_effects_formula, formula)
	}

        #############################################################
        # Determine the function and summary for the model selected #
        #############################################################
        
        ################
        # Linear Model #
        ################
        
        if (model == "LM") {
            if (is.null(random_effects_formula)) {
                model_function <-
                    function(formula, data, na.action) {
                        return(glm(
                            formula,
                            data = data,
                            family = 'gaussian',
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                    lm_summary <- summary(fit)$coefficients
                    para <- as.data.frame(lm_summary)[-1, -3]
                    para$name <- rownames(lm_summary)[-1]
                    return(para)
                }
            } else {
              ranef_function <- lme4::ranef
              model_function <-
                    function(formula, data, na.action) {
                        return(lmerTest::lmer(
                            formula, 
                            data = data, 
                            na.action = na.action))
                    }
                summary_function <- function(fit) {
                    lm_summary <- coef(summary(fit))
                    para <- as.data.frame(lm_summary)[-1, -c(3:4)]
                    para$name <- rownames(lm_summary)[-1]
                    return(para)
                }
            }
        }
        
        ####################
        # Compound Poisson #
        ####################
      
        if (model == "CPLM") {
            if (is.null(random_effects_formula)) {
              model_function <- cplm::cpglm
              summary_function <- function(fit) {
                cplm_out <-
                  capture.output(
                    cplm_summary <- cplm::summary(fit)$coefficients)
                    para <- as.data.frame(cplm_summary)[-1, -3]
                    para$name <- rownames(cplm_summary)[-1]
                    logging::logdebug(
                        "Summary output\n%s", 
                        paste(cplm_out, collapse = "\n"))
                    return(para)
                    }
              } else {
                ranef_function <- glmmTMB::ranef
                model_function <-
                    function(formula, data, na.action) {
                        return(glmmTMB::glmmTMB(
                            formula,
                            data = data,
                            family=glmmTMB::tweedie(link = "log"),
                            ziformula = ~0,
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                  glmmTMB_summary <- coef(summary(fit))
                  para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                  para$name <- rownames(glmmTMB_summary$cond)[-1]
                  return(para)
                }
              }
          }
        
        #####################
        # Negative Binomial #
        #####################
        
        if (model == "NEGBIN") {
            if (is.null(random_effects_formula)) {
                model_function <- MASS::glm.nb
                summary_function <- function(fit) {
                  glm_summary <- summary(fit)$coefficients
                  para <- as.data.frame(glm_summary)[-1, -3]
                  para$name <- rownames(glm_summary)[-1]
                  return(para)
                  }
            } else {
              ranef_function <- glmmTMB::ranef
              model_function <-
                    function(formula, data, na.action) {
                        return(glmmTMB::glmmTMB(
                            formula,
                            data = data,
                            family=glmmTMB::nbinom2(link = "log"),
                            ziformula = ~0,
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                  glmmTMB_summary <- coef(summary(fit))
                  para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                  para$name <- rownames(glmmTMB_summary$cond)[-1]
                  return(para)
                }
            }
          }
        
        ###################################
        # Zero-inflated Negative Binomial #
        ###################################
        
        if (model == "ZINB") {
          if (is.null(random_effects_formula)) {
            model_function <-
              function(formula, data, na.action) {
                return(pscl::zeroinfl(
                  formula,
                  data = data,
                  dist = "negbin",
                  na.action = na.action))
                }
            summary_function <- function(fit) {
              pscl_summary <- summary(fit)$coefficients$count
              para <-as.data.frame(pscl_summary)[-c(1, (ncol(metadata) + 2)), -3]
	            para$name <- rownames(pscl_summary)[-c(1, (ncol(metadata) + 2))]
	            return(para)
	            }
            } else {
              ranef_function <- glmmTMB::ranef
              model_function <-
                function(formula, data, na.action) {
                  return(glmmTMB::glmmTMB(
                    formula,
                    data = data,
                    family=glmmTMB::nbinom2(link = "log"),
                    ziformula = ~1,
                    na.action = na.action))
                  }
              summary_function <- function(fit) {
                glmmTMB_summary <- coef(summary(fit))
                para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                para$name <- rownames(glmmTMB_summary$cond)[-1]
                return(para)
              }
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
        outputs <-
            pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
                # Extract Features One by One
                featuresVector <- features[, x]
                
                # Fit Model
                logging::loginfo(
                    "Fitting model to feature number %d, %s",
                    x,
                    colnames(features)[x])
                dat_sub <-
                    data.frame(expr = as.numeric(featuresVector), metadata)
                fit <- tryCatch({
                    fit1 <-
                        model_function(
                            formula, 
                            data = dat_sub, 
                            na.action = na.exclude)
                }, error = function(err) {
                    fit1 <-
                        try({
                            model_function(
                                formula, 
                                data = dat_sub, 
                                na.action = na.exclude)
                        })
                    return(fit1)
                })
                # Gather Output
                output <- list()
                if (all(!inherits(fit, "try-error"))) {
                    output$para <- summary_function(fit)
                    output$residuals <- residuals(fit)
                    output$fitted <- fitted(fit)
                    output$fit <- fit
                    output$name <- colnames(features)[x]
                    if (!(is.null(random_effects_formula))) {
                      l <- ranef_function(fit)
                      d<-as.vector(unlist(l))
                      names(d)<-unlist(lapply(l, row.names))
                      output$ranef<-d
                      }
                    }
                else
                  {
                    logging::logwarn(paste(
                        "Fitting problem for feature", 
                        x, 
                        "returning NA"))
                    output$para <-
                        as.data.frame(matrix(NA, 
                            nrow = ncol(metadata), ncol = 3))
                    output$para$name <- colnames(metadata)
                    output$residuals <- NA
                    output$fitted <- NA
                    output$fit <- NA
                    if (!(is.null(random_effects_formula))) output$ranef <- NA
                  }
                colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
                output$para$feature <- colnames(features)[x]
                return(output)
                })
        
        # stop the cluster
        if (!is.null(cluster))
            parallel::stopCluster(cluster)
        
        # bind the results for each feature
        paras <-
            do.call(rbind, lapply(outputs, function(x) {
                return(x$para)
            }))
        residuals <-
            do.call(rbind, lapply(outputs, function(x) {
                return(x$residuals)
            }))
        row.names(residuals) <- colnames(features)
        
        fitted <-
          do.call(rbind, lapply(outputs, function(x) {
            return(x$fitted)
          }))
        row.names(fitted) <- colnames(features)   
        
        if (!(is.null(random_effects_formula))) {
          ranef <-
            do.call(rbind, lapply(outputs, function(x) {
              return(x$ranef)
            }))
          row.names(ranef) <- colnames(features) 
        }
 
        #################################
        # Save Fit objects for analysis #
        #################################
        if (!is.na(save_fits_to)){
          logging::loginfo("Saving %s fit objects", length(outputs))
          fits <- lapply(outputs, function(x){x$fit})
          names(fits) <- sapply(outputs, function(x){x$name})
          saveRDS(fits, file = save_fits_to)
          
        }
        ################################
        # Apply correction to p-values #
        ################################
        
        paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
        
        #####################################################
        # Determine the metadata names from the model names #
        #####################################################
        
        metadata_names <- colnames(metadata)
        # order the metadata names by decreasing length
        metadata_names_ordered <-
            metadata_names[order(
                nchar(metadata_names), decreasing = TRUE)]
        # find the metadata name based on the match 
        # to the beginning of the string
        extract_metadata_name <- function(name) {
            return(metadata_names_ordered[mapply(
                startsWith, 
                name, 
                metadata_names_ordered)][1])
        }
        paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
        # compute the value as the model contrast minus metadata
        paras$value <-
            mapply(function(x, y) {
                if (x == y)
                    x
                else
                    gsub(x, "", y)
            }, paras$metadata, paras$name)
        
        ##############################
        # Sort by decreasing q-value #
        ##############################
        
        paras <- paras[order(paras$qval, decreasing = FALSE), ]
        paras <-
            dplyr::select(
                paras,
                c('feature', 'metadata', 'value'),
                dplyr::everything())
        rownames(paras)<-NULL
        
        if (!(is.null(random_effects_formula))) {
          return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = ranef))
        } else {
          return(list("results" = paras, "residuals" = residuals, "fitted" = fitted))
        }
    }        
          
