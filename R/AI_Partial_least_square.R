

pls_predict_boost <- function(data, indices, best_ncomp, geno_test) {
  x_train <- data[indices, -1]  # ensure 'drop = FALSE' to keep the data frame structure if one column
  y_train <- data[indices, 1]

  AI_fit <- pls::plsr(y_train ~ as.matrix(x_train), ncomp = best_ncomp)
  return(stats::predict(AI_fit, newdata = as.matrix(geno_test), ncomp = best_ncomp,
                 reshape= TRUE))
}

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
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   ncomp = 3,
                   pls_paras_tunning= list(ncomp = 10), # number of components
                   resample_method_tune = "cv", # c("cv","boot")
                   N_feature_impo = 10,
                   CI_width_thresholds = c(0.33, 0.66),
                   high_reliability_thres = 0.9,
                   low_reliability_thres = 0.5,
                   n_components = 20,
                   threshold = 100,
                   target = "test_set",
                   iqr_multiplier = 1.5,
                   interval_width_high_threshold = NULL,
                   interval_width_low_threshold = NULL,
                   interval_width_moderate_threshold = NULL,
                   n_bootstrap = 100,
                   system_database = FALSE,
                   ...

){
  msg <- sprintf("==================================================\n")
## Standardizing features is beneficial in PLS to
## ensure that variables with larger scales do not dominate the model
  # if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
  #   if(isTRUE(scaling) || isFALSE(scaling)){
  #     geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
  #   }
  #
  # }
  # if(!is.null(geno_omic_test_object)){
  #   if(isTRUE(scaling) || isFALSE(scaling)){
  #     geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
  #   }
  #
  # }
  # Define cross-validation control
  cv <- caret::trainControl(method = resample_method_tune, number = 10)
################
  scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
  geno_omic_object <- stats::predict(scaler, geno_omic_object)

  if(!is.null(geno_omic_test_object)){
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
  }

################
  y_train <-  pheno_object[, response]
  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

  data <- cbind(y=y_train_scaled, geno_omic_object)


  # Perform cross-validated PLS regression
  if(isTRUE(para_tunning)){

    ncomp <- pls_paras_tunning$ncomp
  cv_results <- caret::train(x = geno_omic_object,
                             y = y_train_scaled,
                             method = "pls",
                             trControl = cv,
                             tuneGrid = data.frame(ncomp = 1:ncomp),
                             metric = "RMSE")

  # Extract the cross-validated RMSE values
  cv_rmse <- cv_results$results$RMSE

  # Find the optimal number of components with minimum RMSE
  optimal_components <- which.min(cv_rmse)


  if(!is.null(geno_omic_test_object)){
    GID <- rownames(geno_omic_test_object)
    boot_results <- boot::boot(data = data,
                               statistic=pls_predict_boost,
                               R=n_bootstrap,
                               best_ncomp = optimal_components,
                               geno_test = geno_omic_test_object)

  } else {
    GID <- rownames(geno_omic_object)
    boot_results <- boot::boot(data = data,
                               statistic=pls_predict_boost,
                               R=n_bootstrap,
                               best_ncomp = optimal_components,
                               geno_test = geno_omic_object)
  }


  } else {



  #   pls_model <- pls::plsr(yy~ geno_omic_object,
  #                          scale = FALSE,
  #                          center = FALSE,
  #                          #ncomp = optimal_components,
  #                          validation = "none")
  #   cumulative_explained_variance <- cumsum(pls::explvar(pls_model))
  #   # Find the number of components explaining at least 90% of the variance
  #   num_components <- which(cumulative_explained_variance >= 90)[1]
  #   # Handle case where no components meet the criterion
  #   if (is.na(num_components) | length(num_components)==0) {
  #     num_components <- which(cumulative_explained_variance >= 50)[1]
  #     if (is.na(num_components) | length(num_components)==0) {
  #       num_components <- length(cumulative_explained_variance)
  #     }
  #     #num_components <- length(cumulative_explained_variance)  # Use max number of components or some default
  #     message("No components explain at least 90% of the variance.")
  #   }
  #
  #   if (is.null(ncomp) | !is.numeric(ncomp)) {
  #     ncomp <- num_components  # Default to using 'num_components' if 'ncomp' is not defined
  #   }
  #
  #   if(ncomp< num_components){
  #     optimal_components <- num_components
  # message(paste(msg, "The number of component provided explain less than 90% of the variance. We make adjustment as this might affect final result."))
  #   }else {
  #   optimal_components <-  ncomp
  #
  #   }


    if(is.null(ncomp)){
    pls_model <- pls::plsr(y_train_scaled ~ as.matrix(geno_omic_object), validation = "CV", segments = 5,
                           center = FALSE)

    # Get the cross-validated RMSEP values
    rmsep_values <- pls::RMSEP(pls_model)

    # Extract the RMSEP values for cross-validation
    rmsep_cv <- rmsep_values$val["CV", , ]

    # Find the optimal number of components
    optimal_components <- which.min(rmsep_cv)

    } else {
      optimal_components <- ncomp
    }

    if(!is.null(geno_omic_test_object)){
      GID <- rownames(geno_omic_test_object)
      boot_results <- boot::boot(data = data,
                                 statistic=pls_predict_boost,
                                 R=n_bootstrap,
                                 best_ncomp = optimal_components,
                                 geno_test = geno_omic_test_object)

    } else {
      GID <- rownames(geno_omic_object)
      boot_results <- boot::boot(data = data,
                                 statistic=pls_predict_boost,
                                 R=n_bootstrap,
                                 best_ncomp = optimal_components,
                                 geno_test = geno_omic_object)
    }



  }

# Construct output
pred_variances <- apply(boot_results$t, 2, var)
pred_SE <- apply(boot_results$t, 2, sd)
AI_pred <- apply(boot_results$t, 2, mean)
genetic_var <- var(AI_pred)
AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)

result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(boot_results = boot_results,
                                                       CI_width_thresholds = CI_width_thresholds)
result_rel <-  reliability_thresholds(prediction_error_var = pred_variances,
                                      genetic_var = genetic_var,
                                      high_reliability_thres = high_reliability_thres,
                                      low_reliability_thres = low_reliability_thres)

composite_reliability <- composite_reliability_tst(geno_trn = geno_omic_object,
                                                   geno_tst = geno_omic_test_object,
                                                   geno_tst_trn = if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) geno_omic_object else NULL,
                                                   names_tst = NULL,
                                                   names_trn = NULL,
                                                   n_components = n_components,
                                                   threshold = threshold,
                                                   target = target,
                                                   interval_width = result_rel_MPIW$Uncertainty,
                                                   CI_width_thresholds = CI_width_thresholds,
                                                   interval_width_high_threshold = interval_width_high_threshold,
                                                   interval_width_low_threshold = interval_width_low_threshold,
                                                   apply_pca = TRUE)
AI_preds <- data.frame(name = GID,
                       Predicted_value = AI_pred_reverted,
                       Standard_error = pred_SE,
                       PEV = pred_variances,
                       lower_bound = result_rel_MPIW$lower_bound,
                       upper_bound = result_rel_MPIW$upper_bound,
                       Uncertainty = result_rel_MPIW$Uncertainty,
                       Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
                       Reliability = result_rel$reliability,
                       Reliability_remarks = result_rel$remarks,
                       Reliability_percentage = result_rel$reliability_percentage,
                       Composite_reliability = composite_reliability$trustworthiness,
                       Composite_reliability_percentage = composite_reliability$reliability_percentage,
                       stringsAsFactors = FALSE)

names(AI_preds)[1] <- c(gen_name)

  row.names(AI_preds) <- NULL

  diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = boot_results,
                                                      GID_names = GID,
                                                      CI_width_thresholds = CI_width_thresholds,
                                                      predictions = AI_pred_reverted,
                                                      standard_errors = pred_SE,
                                                      prediction_error_var = pred_variances,
                                                      genetic_var = genetic_var,
                                                      confidence_level = 0.95,
                                                      model_for_CI_cal = "ML",
                                                      composite_reliability_score = composite_reliability$reliability_score,
                                                      composite_reliability = composite_reliability$trustworthiness,
                                                      composite_reliability_percentage = composite_reliability$reliability_percentage,
                                                      #threshold = NULL,
                                                      high_reliability_thres = high_reliability_thres,
                                                      low_reliability_thres = low_reliability_thres,
                                                      system_database = system_database)
  #}
  model_para <-  data.frame(stat = "component",
                             summary = optimal_components)

  colnames(model_para)[1:2] <- c("stat", "summary")
  output = list(model_parameters = model_para,
                predicted_values = AI_preds,
                diagnostic_plots = diagnostic_plots)


  return(output)

}
