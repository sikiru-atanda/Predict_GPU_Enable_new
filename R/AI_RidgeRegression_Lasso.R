
#' Title
#'
#' @param gen_name
#' @param core
#' @param message
#' @param para_tunning
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param ...
#' @param lasso_paras_tunning
#' @param lambda
#' @param GS_model
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_RidgeRegression_Lasso <- function(pheno_object=NULL,
                                     geno_omic_object = NULL,
                                     geno_omic_test_object = NULL,
                                     response=NULL,
                                     gen_name=NULL,
                                     message = TRUE,
                                     scaling = TRUE,
                                     centering = FALSE,
                                     omic_count = NULL,
                                     AI_cv_nfolds = 5,
                                     para_tunning = FALSE,
                                     lasso_paras_tunning= list(lambda_tune=seq(0.000001,0.9,length.out=100)^4),
                                     lambda_rr = NULL,
                                     GS_model = c("Lasso",
                                                "Ridge_Regression"),
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

  msg <- "\n==================================================\n"

  if(is.null(geno_omic_object) & is.null(pheno_object)) {

    stop(print(paste(msg,"provide matrix of the predictors and the data.frame of the Y variable.")), call. = FALSE)
  }

  if(!is.null(geno_omic_object)){
    GID <- rownames(geno_omic_object)
    scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)

  }

  cols_with_na <- which(colSums(is.na(geno_omic_object)) > 0)
  if(length(cols_with_na)!=0){
    geno_omic_object <- geno_omic_object[, -cols_with_na]
    if(!is.null(geno_omic_test_object)){
      geno_omic_test_object <- geno_omic_test_object[, -cols_with_na]
    }
  }

  if(!is.null(geno_omic_test_object)){
    test_label <- rownames(geno_omic_test_object)
    geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
    GID <- rownames(geno_omic_test_object)
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)

  }

  y_train <-  pheno_object[, response]
  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

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
  #msg <- "\n==================================================\n"
  #####################################################################
  # if (is.null(lambda_rr) & !inherits(lambda_rr, 'numeric')){
  #
  #   lambda <-  seq(0.000001,0.9,length.out=100)^4
  # }

  if(GS_model=="Ridge_Regression"){

    alpha <-  0
  }

  if(GS_model=="Lasso"){

    alpha <-  1
  }


## when length of response variable is 1
#########################
#
#      if(isTRUE(para_tunning)){
#        # create hyperparameter grid
#        AI_grid <- expand.grid(
#          #k = data.frame(k = seq(11,85,by = 2)))
#          lambda = data.frame(alpha = alpha, lambda = lambda_tune)
#
#        )
#
#        AI_trcontrol <- caret::trainControl(method = "cv",
#                                           number = AI_cv_nfolds,
#                                           verboseIter = TRUE,
#                                           returnData = FALSE,
#                                           returnResamp = "all",
#                                           allowParallel = TRUE)
#
#        if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
#          GID <- rownames(geno_omic_object)
#          AI_fit <-  caret::train(x = geno_omic_object,
#                                y = pheno_object[, response],
#                                trControl = AI_trcontrol,
#                                tuneGrid = AI_grid,
#                                GS_model = "glment")
#
#          if(!is.null(geno_omic_test_object)){
#            GID <- rownames(geno_omic_test_object)
#          AI_preds <- stats::predict(AI_fit,
#                                     geno_omic_test_object,
#                                     reshape = TRUE)
#
#          } else {
#
#            if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)){
#              AI_preds <- stats::predict(AI_fit,
#                                         geno_omic_object,
#                                         reshape = TRUE)
#            }
#          }
#
#          AI_preds <- data.frame(name = GID,
#                                 Predicted_value = as.data.frame(AI_preds),
#                                 Standard_error = NA,
#                                 PEV = NA,
#                                 Reliability = NA,
#                                 stringsAsFactors = FALSE)
#          names(AI_preds)[1:2] <-  c(gen_name, "Predicted_value")
#
#          bestTune <- c(AI_fit$bestTune, AI_fit$GS_model)
#
#
#        } else {
#          stop("Training set missing")
#        }
#        ####### Start when their is no need for tunning
#      } else {

  lasso_predict_boost <- function(data, indices, alpha, best_lambda, geno_test) {
    x_train <- data[indices, -1]  # ensure 'drop = FALSE' to keep the data frame structure if one column
    y_train <- data[indices, 1]

    fit <- glmnet::glmnet(as.matrix(x_train), y_train,
                          alpha=alpha, lambda=best_lambda,
                          standardize = FALSE)
    return(stats::predict(object = fit, s=best_lambda, newx=geno_test))
  }

  data = cbind(y = y_train_scaled, geno_omic_object)

  cv_lasso <- glmnet::cv.glmnet(as.matrix(geno_omic_object),
                                y_train_scaled,
                                alpha = alpha)

  # Optimal λ value
  best_lambda <- cv_lasso$lambda.min

  if(!is.null(geno_omic_test_object)){
    #GID <- rownames(geno_omic_test_object)
  boot_results <- boot::boot(data = data,
                             statistic=lasso_predict_boost,
                             R=n_bootstrap,
                             alpha = alpha,
                             best_lambda = best_lambda,
                             geno_test = geno_omic_test_object)

  } else {
    #GID <- rownames(geno_omic_object)
    boot_results <- boot::boot(data = data,
                               statistic=lasso_predict_boost,
                               R=n_bootstrap,
                               alpha = alpha,
                               best_lambda = best_lambda,
                               geno_test = geno_omic_object)
  }


  # Construct output
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
  #AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)
  AI_pred_reverted <-  AI_pred


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


  if(!is.null(geno_omic_test_object)){

    train_test_label <- ifelse(rownames(geno_omic_test_object)%in%test_label, "Test", "Train")
  } else {
    train_test_label <- rep("Train", nrow(geno_omic_object))
  }

  AI_preds <- data.frame(name = GID,
                         Predicted_value = AI_pred_reverted,
                         Train_Test_Label = train_test_label,
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


     #} ## End of when no need for tunning.


     model_para <- data.frame(stat = c("lambda", "alpha"),
                              summary = best_lambda,
                              stringsAsFactors = FALSE)

     colnames(model_para)[1:2] <- c("stat", "summary")

     output <-  list(model_parameters = model_para,
                     predicted_values = AI_preds,
                     diagnostic_plots = diagnostic_plots)


return(output)


 }
