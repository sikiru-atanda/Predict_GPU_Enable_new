
knn_predict_boost <- function(data, indices, best_k, geno_test) {
  x_train <- data[indices, -1]  # ensure 'drop = FALSE' to keep the data frame structure if one column
  y_train <- data[indices, 1]

  AI_fit <-  caret::knnreg(x = x_train,
                         y = y_train,
                         k = best_k)
  return(stats::predict(AI_fit, geno_test, reshape= TRUE))
}

#' Title
#'
#' @param gen_name
#' @param core
#' @param message
#' @param center
#' @param para_tunning
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param knn_paras_tunning
#' @param k
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_knn <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   knn_paras_tunning= list(k = seq(3, 21, by = 2)),
                   k = 5,
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


  ## KNN relies on distance metrics to find the nearest neighbors,
  ## so scaling the features is critical to ensure that
  ## all dimensions contribute equally to the distance calculations.
  if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
    scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)

    # if(isTRUE(scaling) || isFALSE(scaling)){
    #   geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
    # }

  }
  if(!is.null(geno_omic_test_object)){
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
    # if(isTRUE(scaling) || isFALSE(scaling)){
    #   geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
    # }

  }

  y_train <-  pheno_object[, response]
  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

  data <- cbind(y=y_train_scaled, geno_omic_object)
  #msg <- "\n==================================================\n"

     if(isTRUE(para_tunning)){
       # create hyperparameter grid
if(is.null(knn_paras_tunning) | length(knn_paras_tunning)!=0){
       AI_grid <- expand.grid(k = knn_paras_tunning$k)

} else {
  AI_grid <- expand.grid(
    #k = data.frame(k = seq(3, 21, by = 2)),
    k = seq(3, 21, by = 2)
  )

}

       AI_trcontrol = caret::trainControl(method = "cv",
                                          number = AI_cv_nfolds,
                                          verboseIter = TRUE,
                                          returnData = FALSE,
                                          returnResamp = "all",
                                          allowParallel = TRUE)

       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {

         GID <- rownames(geno_omic_object)
         AI_fit = caret::train(x = geno_omic_object,
                               y = pheno_object[, response],
                               trControl = AI_trcontrol,
                               tuneGrid = AI_grid,
                               method = "knn")

         best_k <- AI_fit$bestTune$k

if(!is.null(geno_omic_test_object)){

  GID <- rownames(geno_omic_test_object)
  boot_results <- boot::boot(data,
                             statistic=knn_predict_boost,
                             best_k = best_k,
                             geno_test = geno_omic_test_object,
                             R=n_bootstrap)

} else {

  if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){

    GID <- rownames(geno_omic_object)
    boot_results <- boot::boot(data,
                               statistic=knn_predict_boost,
                               best_k = best_k,
                               geno_test = geno_omic_object,
                               R=n_bootstrap)
}

}



         model_para <- data.frame(stat ="k",
                                  summary = AI_fit$bestTune$k)


       } else {
         stop("Training set missing")
       }

       ####### Start when their is no need for tunning


     } else {


       if(!is.null(geno_omic_test_object)){

         GID <- rownames(geno_omic_test_object)
         boot_results <- boot::boot(data,
                                    statistic=knn_predict_boost,
                                    best_k = k,
                                    geno_test = geno_omic_test_object,
                                    R=n_bootstrap)

       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){

           GID <- rownames(geno_omic_object)
           boot_results <- boot::boot(data,
                                      statistic=knn_predict_boost,
                                      best_k = k,
                                      geno_test = geno_omic_object,
                                      R=n_bootstrap)
         }

       }

       model_para <- data.frame(stat ="k",
                                summary = k)

     } ## End of when no need for tunning.

  # Construct output
  #boot_results$t <- revert_scaling(boot_results$t, y_scaler)

  revert_scaling_ml <- function(scaled_values, scaler_mean, scaler_sd) {
    scaled_values * scaler_sd + scaler_mean
  }

  boot_results$t <- apply(
    boot_results$t,
    2,
    function(col_vec) revert_scaling_ml(col_vec, y_scaler$mean, y_scaler$std)
  )

  pred_variances <- apply(boot_results$t, 2, var)
  pred_SE <- apply(boot_results$t, 2, sd)
  AI_pred <- apply(boot_results$t, 2, mean)
  genetic_var <- var(AI_pred)
  AI_pred_reverted <-  AI_pred
  #AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)

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
                         #Composite_reliability = composite_reliability$trustworthiness,
                         #Composite_reliability_percentage = composite_reliability$reliability_percentage,
                         stringsAsFactors = FALSE)

  names(AI_preds)[1] <- c(gen_name)

  diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = boot_results,
                                                      GID_names = GID,
                                                      CI_width_thresholds = CI_width_thresholds,
                                                      predictions = AI_pred_reverted,
                                                      standard_errors = pred_SE,
                                                      prediction_error_var = pred_variances,
                                                      genetic_var = genetic_var,
                                                      confidence_level = 0.95,
                                                      model_for_CI_cal = "ML",
                                                      #composite_reliability_score = composite_reliability$reliability_score,
                                                      #composite_reliability = composite_reliability$trustworthiness,
                                                      #composite_reliability_percentage = composite_reliability$reliability_percentage,
                                                      #threshold = NULL,
                                                      high_reliability_thres = high_reliability_thres,
                                                      low_reliability_thres = low_reliability_thres,
                                                      system_database = system_database)





     output <- list(model_parameters = model_para,
                   predicted_values = AI_preds,
                   diagnostic_plots = diagnostic_plots
                   )




return(output)



 }
