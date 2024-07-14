

# Function to train and predict using random for bootstrapping
train_predict_randomForest <- function(data_label_geno, indices, test_geno, ntree = NULL,
                                       mtry=NULL, maxnodes = NULL) {
  train_data <- data_label_geno[, -1]
  y_train <- data_label_geno[, 1]
  # Subset the data
  if(((is.null(mtry) & is.null(maxnodes)) & is.null(ntree))){
    model <- randomForest::randomForest(x = train_data,
                                        y = y_train
                                        )

  }else if(is.null(mtry) & is.null(maxnodes)){
    model <- randomForest::randomForest(x = train_data,
                                        y = y_train,
                                        ntree = ntree)

  } else if (!is.null(mtry) & is.null(maxnodes)){
    model <- randomForest::randomForest(x = train_data,
                                        y = y_train,
                                        ntree = ntree,
                                        mtry = mtry)


  } else if (!is.null(mtry) & !is.null(maxnodes)){
    model <- randomForest::randomForest(x = train_data,
                                        y = y_train,
                                        ntree = ntree,
                                        mtry = mtry,
                                        maxnodes = maxnodes)


  } else {

    if (is.null(mtry) & !is.null(maxnodes)){
      model <- randomForest::randomForest(x = train_data,
                                          y = y_train,
                                          ntree = ntree,
                                          maxnodes = maxnodes)

    }

  }
  # Predict on the original data
  if(!is.null(test_geno)){
    pred <- stats::predict(model, as.matrix(test_geno), reshape = TRUE)
  } else{
    pred <- stats::predict(model, as.matrix(train_data), reshape = TRUE)
  }
  return(pred)
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
#' @param rf_paras_tunning
#' @param ntree
#' @param mtry
#' @param maxnodes
#' @param importance
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_randomForest <- function(pheno_object=NULL,
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
                            rf_paras_tunning= list(mtry = TRUE,
                                                   ntree = c(500, 1000, 1500),
                                                   nodesize = c(1, 5, 10),
                                                   maxnodes = c(30, 50, NULL)),  # NULL means no limit),
                            ntree=500,
                            mtry = NULL,
                            maxnodes = NULL,
                            importance=TRUE,
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

  if(!is.null(geno_omic_object)){
    scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)

  }

  if(!is.null(geno_omic_test_object)){
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)

  }

  y_train <-  pheno_object[, response]
  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

  # if(!is.null(geno_omic_object)){
  #   if(isTRUE(scaling)){
  #     geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
  #   } else{
  #     if(isTRUE(centering) && !is.null(omic_count)){
  #       geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = FALSE)
  #     }
  #   }
  # }

  # if(!is.null(geno_omic_test_object)){
  #   if(isTRUE(scaling)){
  #     geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
  #   } else{
  #     if(isTRUE(centering) && !is.null(omic_count)){
  #       geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = FALSE)
  #     }
  #   }
  # }
  #########################

  if(isTRUE(para_tunning)){
    # create hyperparameter grid

    #n_features <- ncol(geno_omic_object)

    if(isTRUE(mtry)){
      mtry <- c(sqrt(ncol(geno_omic_object)), sqrt(ncol(geno_omic_object))/2, ncol(geno_omic_object)/3)
    } else {
      stop(print(paste(msg,"set mtry equal TRUE: mtry = TRUE")), call. = FALSE)

    }

    AI_grid <- expand.grid(
      mtry = mtry,
      ntree = para_tunning$ntree,
      nodesize = para_tunning$nodesize,
      maxnodes = para_tunning$maxnodes  # NULL means no limit
    )

    AI_trcontrol = caret::trainControl(method = "cv",
                                       number = AI_cv_nfolds,
                                       verboseIter = TRUE,
                                       returnData = FALSE,
                                       returnResamp = "all",
                                       allowParallel = TRUE)

    if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
      GID <- rownames(geno_omic_object)
      AI_fit <-  caret::train(x = geno_omic_object,
                            y = y_train_scaled,
                            trControl = AI_trcontrol,
                            tuneGrid = AI_grid,
                            method = "rf")

    }
    ####### Start when their is no need for tunning
    ntree <-  AI_fit$bestTune$ntree
    maxnodes <-  AI_fit$bestTune$maxnodes
    mtry <-  AI_fit$bestTune$mtry

  }


  if(!is.null(geno_omic_test_object)){
       GID <- rownames(geno_omic_test_object)
  } else {
    if(!is.null(geno_omic_object)){
      GID <- rownames(geno_omic_object)
    }
  }

  data_label_geno <- cbind(y_train_scaled,geno_omic_object)


    if(!is.null(geno_omic_object) & !is.null(pheno_object)) {

      boot_results <- boot::boot(
                                data = data_label_geno,
                                statistic = train_predict_randomForest,
                                ntree = ntree,
                                maxnodes = maxnodes,
                                mtry = mtry,
                                R = n_bootstrap,  # Number of bootstrap samples
                                #sim = "ordinary",
                                test_geno = geno_omic_test_object
                              )

      # Extract bootstrap predictions
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
                                            low_reliability_thres = low_reliability_thres
      )

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


      names(AI_preds)[1] <-  c(gen_name)

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



    } else {

      stop(print(paste(msg,"Training set missing.")), call. = FALSE)


    }


  #} ## End of when no need for tunning.


  model_para <- data.frame(stat = c("ntree","mtry", "maxnodes"),
                           summary = c(ntree, mtry, maxnodes),
                           stringsAsFactors = FALSE)

  colnames(model_para)[1:2] <- c("stat", "summary")

  output <-  list(model_parameters = model_para,
                  predicted_values = AI_preds,
                  diagnostic_plots = diagnostic_plots)


  return(output)

}
