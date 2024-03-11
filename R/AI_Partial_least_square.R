

#' Title
#'
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param gen_name
#' @param message
#' @param center
#' @param para_tunning
#' @param ncomp
#' @param pls_paras_tunning
#' @param resample_method_tune
#' @param N_feature_impo
#' @param core
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
AI_pls <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   message = TRUE,
                   scale = TRUE,
                   para_tunning = FALSE,
                   ncomp = 3,
                   pls_paras_tunning= c(ncomp = 10), # number of components
                   resample_method_tune = "cv", # c("cv","boot")
                   N_feature_impo = 10,
                   core = NULL,
                   ...

){
  msg <- sprintf("==================================================\n")

  # Define cross-validation control
  cv <- caret::trainControl(method = resample_method_tune, number = 10)

  yy = as.numeric(pheno_object[, response])
  # Perform cross-validated PLS regression
  if(isTRUE(para_tunning)){

    ncomp <- pls_paras_tunning$ncomp
  cv_results <- caret::train(x = geno_omic_object,
                             y = yy,
                             method = "pls",
                             trControl = cv, tuneGrid = data.frame(ncomp = 1:ncomp),
                             metric = "RMSE")

  # Extract the cross-validated RMSE values
  cv_rmse <- cv_results$results$RMSE

  # Find the optimal number of components with minimum RMSE
  optimal_components <- which.min(cv_rmse)

  } else {

    pls_model <- pls::plsr(yy~ geno_omic_object,
                           scale = FALSE,
                           center = FALSE,
                           ncomp = optimal_components,
                           validation = "none")
    cumulative_explained_variance <- cumsum(pls::explvar(pls_model))
    # Find the number of components explaining at least 90% of the variance
    num_components <- which(cumulative_explained_variance >= 90)[1]
    if(ncomp< num_components){
      optimal_components <- num_components
  message(paste(msg, "The number of component provided explain less than 90% of the variance. We make adjustment as this might affect final result."))
    }else {
    optimal_components <-  ncomp

    }
  }

  # Fit the final PLS model with the optimal number of components
  #if(!is.null(geno_omic_object)){
  GID <- rownames(geno_omic_object)
  pls_model <- pls::plsr(yy~ geno_omic_object,
                         scale = FALSE,
                         center = FALSE,
                         ncomp = optimal_components,
                         validation = "none")

  if(!is.null(geno_omic_test_object)){
    GID <- rownames(geno_omic_test_object)
  # Fit PLS model with optimal number of components

  # Predict the response variable using the final PLS model
    AI_preds <- stats::predict(pls_model,
                                newdata =geno_omic_test_object,
                                ncomp = optimal_components)

  } else {
    if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){
      AI_preds <- stats::predict(pls_model,
                                  newdata =geno_omic_object,
                                  ncomp = optimal_components)

    }

  }

  AI_preds <- data.frame(name = GID,
                         Predicted_value = as.data.frame(AI_preds),
                         Standard_error = NA,
                         PEV = NA,
                         Reliability = NA,
                         stringsAsFactors = FALSE)


  names(AI_preds)[1:2] = c(gen_name, "Predicted_value")

  row.names(AI_preds) <- NULL
  #}
  model_para <-  data.frame(stat = "component",
                             summary = optimal_components)
  colnames(model_para)[1:2] <- c("stat", "summary")
  output = list(model_para,
                AI_preds,
                pls_model)


  names(output) = c("model_parameters", "predicted_values",
                    "trained_model")

  return(output)

}
