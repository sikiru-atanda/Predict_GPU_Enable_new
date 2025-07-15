
# Cross-validation function
early_stop_for_nround_xgb <- function(X, y, n_folds=5, max_rounds=500, early_stop_rounds=50, params = NULL) {

  if(isFALSE(is.matrix(X))) X <- as.matrix(X)

  # If the dataset size is less than 100, use k-fold cross-validation
  if (nrow(X) < 100) {
    folds <- caret::createFolds(y, k = n_folds, list = TRUE, returnTrain = TRUE)
    best_rounds <- c()

    for (i in seq_along(folds)) {
      train_indices <- folds[[i]]
      val_indices <- setdiff(seq_along(y), train_indices)

      dtrain <- xgboost::xgb.DMatrix(data = X[train_indices, ], label = y[train_indices])
      dval <- xgboost::xgb.DMatrix(data = X[val_indices, ], label = y[val_indices])

      watchlist <- list(train = dtrain, eval = dval)

      xgb_model <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = max_rounds,
        watchlist = watchlist,
        early_stopping_rounds = early_stop_rounds,
        print_every_n = 10,
        maximize = FALSE,
        verbose = 0
      )

      best_rounds <- c(best_rounds, xgb_model$best_iteration)
    }

    best_num_rounds <- ceiling(mean(best_rounds))

  } else {
    # Otherwise, split the data into training and validation sets
    set.seed(42)
    train_indices <- caret::createDataPartition(y, p = 0.8, list = FALSE)
    val_indices <- setdiff(seq_along(y), train_indices)

    dtrain <- xgboost::xgb.DMatrix(data = X[train_indices, ], label = y[train_indices])
    dval <- xgboost::xgb.DMatrix(data = X[val_indices, ], label = y[val_indices])

    watchlist <- list(train = dtrain, eval = dval)

    xgb_model <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = max_rounds,
      watchlist = watchlist,
      early_stopping_rounds = early_stop_rounds,
      print_every_n = 10,
      maximize = FALSE,
      verbose = 0
    )

    best_num_rounds <- xgb_model$best_iteration
  }

  return(best_num_rounds)
}
# Function to train and predict using xgboost for bootstrapping
train_predict_xgboost <- function(data_label_geno, indices, test_geno, params, nrounds) {
  train_data <- data_label_geno[, -1]
  y_train <- data_label_geno[, 1]
  # Subset the data
  dtrain_boot <- xgboost::xgb.DMatrix(data = as.matrix(train_data[indices,]), label = y_train[indices])

  # Train the model
  model <- xgboost::xgboost(data = dtrain_boot, params = params, nrounds = nrounds, verbose = 0)

  # Predict on the original data
  if(!is.null(test_geno)){
    pred <- stats::predict(model, as.matrix(test_geno), reshape = TRUE)
  } else{
    pred <- stats::predict(model, as.matrix(train_data), reshape = TRUE)
  }
  return(pred)
}


#' Title Extreme Gradient Boosting Machine Learning Genomic Selection Pipeline
#' Hyper-parameter tunning  of the parameters is allowed if desired by user
#' The parameters are:
#' Iter_tune : number of boosting iterations
#' learnining_rate_tune:  simply means how fast the model learns.
#' Each tree added modifies the overall model.
#' The magnitude of the modification is controlled by learning rate.
#' The lower the learning rate, the slower the model learns.
#' The advantage of slower learning rate is that the model becomes more robust
#' and efficient.
#'
#'
#' @param gen_name column name containing individuals/genotypes
#' @param message
#' @param xgb_paras_tunning parameter to for tunning
#' @param learning_rate how slow/fast the model learn
#' @param max_depth max_depth refers to the number of leaves of each tree
#' @param subsample  This help to reduce the correlation between results from individual learners.
#' @param booster  to determine if the model is for regression or classification problem
#' @param iteration Number of iteration
#' @param para_tunning  if user required parameter tunning
#' @param pheno_object phenotypic object NA is allowed
#' @param geno_omic_object multi-omic data, NA not allowed
#' @param geno_omic_test_object multi-omic data for testing set if not present in geno_omic_object
#' @param response  y variables/lables
#' @param core number of ram for paralllel job
#' @param resample_method_tune
#' @param number_of_fold_tune
#' @param N_feature_impo number of feature/ x variables to extract based on the importance/weight
#' @param ...
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

AI_Xgb <- function(pheno_object=NULL,
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
                   xgb_paras_tunning= list(Iter_tune = seq(100, 500, 100), # number of boosting iterations
                                        learning_rate_tune = seq(0.01, 0.1, 0.01), # learning rate, low value means model is more robust to overfitting
                                        max_depth = seq(3, 15, 2),
                                        rate_drop = seq(0.05, 0.5, 0.05),
                                        skip_drop = seq(0.05, 1, 0.1),
                                        xgb_gamma = seq(0, 1, 0.01),
                                        colsample_bytree = seq(0.1, 1, 0.1),
                                        min_child_weight = seq(1, 10, 2),
                                        subsample = seq(0.2, 1, 0.1),
                                        L2_tune = seq(0, 1, 0.01), #  for linear gbL2 Regularization (Ridge Regression)
                                        L1_tune = seq(0, 1, 0.01)), # for linear gb
                   resample_method_tune = "cv", # c("cv","boot")
                   number_of_fold_tune = 5,
                   learning_rate = 0.01,
                   max_depth = 6,
                   subsample = 0.7,
                   xgb_booster =  "dart", #"gbtree", # "gblinear",
                   xgb_alpha = 0.001, ## xgboost linear
                   xgb_lambda = 1.0,  # xgboost linear #  dart L2 regularization term on weights (default: 1)
                   xgb_gamma = 0.01, #it acts as a regularization parameter for controlling tree complexity.  dart L1 regularization term on weights (default: 0)
                   iteration = 100,
                   xgb_rate_drop = 0.1,
                   xgb_skip_drop = 0.5,
                   xgb_objective = "reg:squarederror",
                   xgb_sample_type = "uniform",
                   xgb_normalize_type = "tree",
                   #early_stopping_rounds_xgb = TRUE,
                   N_feature_impo = 10,
                   colsample_bytree = 0.7,
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
                   early_stop_for_iteration_xgb = FALSE,
                   system_database = FALSE,
                   ...){
#browser()
  # Auto-adjust nthread based on environment

  #nthread <- if (future::nbrOfWorkers() > 1) 1 else parallel::detectCores(logical = FALSE)
  nthread <- 1 #set_per_worker_threads()

  msg <- "\n==================================================\n"

  if(is.null(geno_omic_object) & is.null(pheno_object)) {

    stop(print(paste(msg,"provide matrix of the predictors and the data.frame of the Y variable.")), call. = FALSE)
  }

  if(!is.null(geno_omic_object)){
    GID <- rownames(geno_omic_object)
    scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)


    cols_with_na <- which(colSums(is.na(geno_omic_object)) > 0)
    if(length(cols_with_na)!=0){
      geno_omic_object <- geno_omic_object[, -cols_with_na]
      if(!is.null(geno_omic_test_object)){
        geno_omic_test_object <- geno_omic_test_object[, -cols_with_na]
      }
    }
    #geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
    # if(isTRUE(scaling)){
    #   geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
    # } else{
    #   if(isTRUE(centering) && !is.null(omic_count)){
    #     geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = FALSE)
    #   }
    # }
  }



  if(!is.null(geno_omic_test_object)){

    test_label <- rownames(geno_omic_test_object)
    geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
    GID <- rownames(geno_omic_test_object)

    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
    # if(isTRUE(scaling)){
    #   geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
    # } else{
    #   if(isTRUE(centering) && !is.null(omic_count)){
    #     geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = FALSE)
    #   }
    # }
  }
  #msg <- "\n==================================================\n"
## when length of response variable is 1
#########################
  y_train <-  pheno_object[, response]
  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]


  #cat("Rows in geno_omic_object:", nrow(geno_omic_object), "\n")
  #cat("Length of y_train_scaled:", length(y_train_scaled), "\n")

#data <- cbind(y=y_train_scaled, geno_omic_object)

if(isTRUE(para_tunning)){

  xgb_grid_tree <- NULL

  xgb_grid_linear <- NULL

  xgb_dart_tree <- NULL

  xgb_trcontrol <-  caret::trainControl(method = "cv",
                                        number = AI_cv_nfolds,
                                        verboseIter = FALSE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)

  if (xgb_booster == "dart") {
    # Define the parameter grid for a tree-based model
    xgb_dart_tree <-  expand.grid(
      nrounds = xgb_paras_tunning$Iter_tune,   # Number of boosting rounds
      max_depth = xgb_paras_tunning$max_depth,  # Varying tree depths
      eta = xgb_paras_tunning$learning_rate_tune,  # Learning rate
      rate_drop = xgb_paras_tunning$rate_drop,
      skip_drop = xgb_paras_tunning$rate_drop,
      colsample_bytree = xgb_paras_tunning$colsample_bytree,  # Subsample ratio of columns when constructing each tree
      min_child_weight = xgb_paras_tunning$min_child_weight,  # Minimum sum of instance weight needed in a child
      subsample = xgb_paras_tunning$subsample
    )

  }

  if (xgb_booster == "gbtree") {
    # Define the parameter grid for a tree-based model
    xgb_grid_tree <-  expand.grid(
      nrounds = unique(xgb_paras_tunning$Iter_tune),   # Number of boosting rounds
      eta = unique(xgb_paras_tunning$learning_rate_tune),  # Learning rate
      max_depth = unique(xgb_paras_tunning$max_depth),  # Varying tree depths
      gamma = unique(xgb_paras_tunning$xgb_gamma),  # Minimum loss reduction required for further partition
      colsample_bytree = unique(xgb_paras_tunning$colsample_bytree),  # Subsample ratio of columns when constructing each tree
      min_child_weight = unique(xgb_paras_tunning$min_child_weight),  # Minimum sum of instance weight needed in a child
      subsample = unique(xgb_paras_tunning$subsample)
    )

  }

if(xgb_booster=="gblinear"){
  # Define the parameter grid for a linear model
  xgb_grid_linear <-  expand.grid(
    nrounds = unique(xgb_paras_tunning$Iter_tune),  # Number of boosting rounds
    eta = unique(xgb_paras_tunning$learning_rate_tune),  # Learning rate
    lambda = unique(xgb_paras_tunning$L2_tune),  # L2 Regularization
    alpha = unique(xgb_paras_tunning$L1_tune)   # L1 Regularization
  )
}

######

    if(!is.null(geno_omic_object)){

      #GID <- rownames(geno_omic_object)
      if(!is.null(xgb_grid_linear)){
   xgb_fit <-  caret::train(x = geno_omic_object,
                   y = y_train_scaled,
                   trControl = xgb_trcontrol,
                   tuneGrid = xgb_grid_linear,
                   method = "xgbLinear")

   best_params <- list(
     booster = xgb_booster,
     objective = "reg:squarederror",
     eta = xgb_fit$bestTune$eta,
     lambda = xgb_fit$bestTune$lambda,  # L2 Regularization
     #max_depth = xgb_fit$bestTune$max_depth,
     alpha = xgb_fit$bestTune$alpha,  # L1 Regularization
     nthread = nthread
   )
   nrounds <- xgb_fit$bestTune$nrounds
      } else {
        if(!is.null(xgb_grid_tree)){
          xgb_fit <-  caret::train(x = geno_omic_object,
                                   y = y_train_scaled,
                                   trControl = xgb_trcontrol,
                                   tuneGrid = xgb_grid_tree,
                                   method = "xgbTree")

          best_params <- list(
            booster = xgb_booster,
            objective = "reg:squarederror",
            eta = xgb_fit$bestTune$eta,
            gamma = xgb_fit$bestTune$gamma,  # Minimum loss reduction required for further partition
            max_depth = xgb_fit$bestTune$max_depth,
            colsample_bytree = xgb_fit$bestTune$colsample_bytree,  # Subsample ratio of columns when constructing each tree
            min_child_weight = xgb_fit$bestTune$min_child_weight,  # Minimum sum of instance weight needed in a child
            subsample = xgb_fit$bestTune$subsample,
            nthread = nthread

          )
          nrounds <- xgb_fit$bestTune$nrounds
          }

        if(!is.null(xgb_dart_tree)){
          # Custom DART method for caret with early stopping
          dart_model <- list(
            type = "Regression",
            library = "xgboost",
            loop = NULL,
            parameters = data.frame(
              parameter = c("nrounds", "max_depth", "eta", "rate_drop", "skip_drop"),
              class = rep("numeric", 5), ## 5 is number of parameters to tune
              label = c("nrounds", "max_depth", "eta", "rate_drop", "skip_drop")
            ),
            grid = function(x, y, len = NULL, search = "grid") {
              expand.grid(
                nrounds = unique(xgb_grid_tree$nrounds),
                max_depth = unique(xgb_grid_tree$max_depth),
                eta = unique(xgb_grid_tree$eta),
                rate_drop = unique(xgb_grid_tree$rate_drop),
                skip_drop = unique(xgb_grid_tree$skip_drop)
              )
            },
            fit = function(x, y, wts, param, lev, last, classProbs, ...) {
              dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x), label = y)
              watchlist <- list(train = dtrain)
              params <- list(
                booster = "dart",
                objective = "reg:squarederror",
                eta = param$eta,
                max_depth = param$max_depth,
                rate_drop = param$rate_drop,
                skip_drop = param$skip_drop
              )
              xgb_model <- xgboost::xgb.train(params = params, data = dtrain, nrounds = param$nrounds, watchlist = watchlist,
                                              early_stopping_rounds = 10, verbose = 0, ...)
              return(xgb_model)
            },
            predict = function(modelFit, newdata, submodels = NULL) {
              stats::predict(modelFit, newdata = xgboost::xgb.DMatrix(data = as.matrix(newdata)))
            },
            prob = NULL
          )

          # Train the model using caret
          xgb_fit <- caret::train(
            x = geno_omic_object,
            y = y_train_scaled,
            trControl = xgb_trcontrol,
            tuneGrid = xgb_dart_tree,
            method = dart_model
          )

          best_params <- list(
            booster = xgb_booster,
            objective = "reg:squarederror",
            eta = xgb_fit$bestTune$eta,
            max_depth = xgb_fit$bestTune$max_depth,
            rate_drop = xgb_fit$bestTune$rate_drop,
            skip_drop = xgb_fit$bestTune$skip_drop,
            sample_type = "uniform",
            normalize_type = "tree",
            nthread = nthread

          )
          nrounds <- xgb_fit$bestTune$nrounds

        }

      }

   if(!is.null(geno_omic_test_object)) {

     data_label_geno = cbind(y_train_scaled,geno_omic_object)

     boot_results <- boot::boot(
       data = data_label_geno,
       statistic = train_predict_xgboost,
       R = n_bootstrap,  # Number of bootstrap samples
       #sim = "ordinary",
       test_geno = geno_omic_test_object,
       params = best_params,
       nrounds = nrounds
     )


     # AI_preds <- stats::predict(xgb_fit,
     #                           geno_omic_test_object,
     #                           reshape = TRUE)

   } else {

     #GID <- rownames(geno_omic_object)
     if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

       data_label_geno = cbind(y_train_scaled, geno_omic_object)

       boot_results <- boot::boot(
         data = data_label_geno,
         statistic = train_predict_xgboost,
         R = n_bootstrap,  # Number of bootstrap samples
         #sim = "ordinary",
         test_geno = geno_omic_object,
         params = best_params,
         nrounds = nrounds
       )

       # AI_preds <- stats::predict(xgb_fit,
       #                             geno_omic_object,
       #                             reshape = TRUE)

     }
   }

      # Extract bootstrap predictions

   bestTune <- c(xgb_fit$bestTune, xgb_fit$method)


   # res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
   #                                 X_train = geno_omic_object,
   #                                 N_feature_impo = N_feature_impo,
   #                                 xgb_booster = xgb_booster)

    } else {

      stop(paste(msg,'Missing training set.'), call. = FALSE)

    }

    ####### Start when their is no need for tunning
    } else {

      if (xgb_booster == "gbtree") {
          xgb_params <- list(
          booster = xgb_booster,
          eta = learning_rate,
          max_depth = max_depth,
          #gamma = xgb_gamma,
          subsample = subsample,
          colsample_bytree = colsample_bytree,
          objective = "reg:squarederror",
          eval_metric = c("rmse", "rmsle", "mape"),
          nthread = nthread
        )
      } else if (xgb_booster == "gblinear") {
        xgb_params <- list(
          booster = xgb_booster,
          max_depth = max_depth,
          subsample = subsample,
          colsample_bytree = colsample_bytree,
          #alpha = xgb_alpha,
          #lambda = xgb_lambda,
          eta = learning_rate,
          objective = "reg:squarederror",
          eval_metric = c("rmse", "rmsle", "mape"),
          nthread = nthread
        )
      } else {
        if (xgb_booster == "dart") {
          xgb_params <- list(
            booster = xgb_booster,
            max_depth = max_depth,
            subsample = subsample,
            colsample_bytree = colsample_bytree,
            #sample_type = "uniform",
            #normalize_type = "tree",
            #rate_drop = xgb_rate_drop,
            #skip_drop = xgb_skip_drop,
            #alpha = xgb_alpha,
            #lambda = xgb_lambda,
            eta = learning_rate,
            objective = "reg:squarederror",
            eval_metric = c("rmse", "rmsle", "mape"),
            nthread = nthread
          )
        }
      }

      if(isTRUE(early_stop_for_iteration_xgb)){
      best_num_rounds <- early_stop_for_nround_xgb(X=geno_omic_object,
                                                   y=y_train_scaled,
                                                   n_folds=5,
                                                   max_rounds= iteration,
                                                   early_stop_rounds=10,
                                                   params=xgb_params)

      nrounds <-  best_num_rounds

      } else {
        nrounds <- iteration
      }

      if(!is.null(geno_omic_test_object)) {
        #GID <- rownames(geno_omic_test_object)
        data_label_geno = cbind(y_train_scaled,geno_omic_object)

        boot_results <- boot::boot(
          data = data_label_geno,
          statistic = train_predict_xgboost,
          R = n_bootstrap,  # Number of bootstrap samples
          #sim = "ordinary",
          test_geno = geno_omic_test_object,
          params = xgb_params,
          nrounds = nrounds
        )


      } else {

        if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
          #GID <- rownames(geno_omic_object)
          data_label_geno = cbind(y_train_scaled,geno_omic_object)

          boot_results <- boot::boot(
            data = data_label_geno,
            statistic = train_predict_xgboost,
            R = n_bootstrap,  # Number of bootstrap samples
            #sim = "ordinary",
            test_geno = geno_omic_object,
            params = xgb_params,
            nrounds = nrounds
          )



          # AI_preds <- stats::predict(xgb_fit,
          #                             geno_omic_object,
          #                             reshape = TRUE)

        }
      }

      # res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
      #                                 X_train = geno_omic_object,
      #                                 N_feature_impo = N_feature_impo,
      #                                 xgb_booster = xgb_booster)



    } ## End of when no need for tunning.

# Extract bootstrap predictions
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


names(AI_preds)[1] <-  c(gen_name)

#####

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

if(isFALSE(para_tunning)){
  #model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)
   model_para <-  data.frame(stat = "nrounds",
                            summary = nrounds,
                            stringsAsFactors = FALSE)
} else{
  model_para <- unlist(xgb_fit$bestTune)
  model_para <- as.data.frame(unlist(model_para))

  model_para <-  data.frame(stat = rownames(model_para),
                            summary = model_para,
                            stringsAsFactors = FALSE)

}

  #model_para <- as.data.frame(unlist(model_para))
  # model_para <-  data.frame(stat = rownames(model_para),
  #                           summary = model_para,
  #                           stringsAsFactors = FALSE)

  colnames(model_para)[1:2] <- c("stat", "summary")

  rownames(model_para) <- NULL

output <-  list(model_parameters = model_para,
                predicted_values = AI_preds,
                diagnostic_plots = diagnostic_plots
              #feature_weight = res_feature
              #trained_model = xgb_fit
              )


return(output)


 }


# Create DMatrix
# geno_omic_object_train <- xgboost::xgb.DMatrix(data = as.matrix(geno_omic_object),
#                                                label = pheno_object[, response])
#
#
# if(isTRUE(early_stopping_rounds_xgb)){
#   # Create 10 bins based on quantiles
#   # pheno_object <- pheno_object |>
#   #   dplyr::mutate(bin = dplyr::ntile(!!dplyr::sym(response), 10))  # Creating 10 bins based on quantiles
#   #
#   #
#   # # Set number of folds for cross-validation
#   # n_folds <- 5
#   #
#   # # Create DMatrix
#   # #dtrain <- xgboost::xgb.DMatrix(data = as.matrix(geno_omic_object), label = pheno_object[[response]])
#   #
#   #
#   # # Define early stopping rounds
#   # early_stopping_fraction <- 0.2  # 20% of iterations
#   # if(is.null(iteration)) iteration <- 100  # Define total number of iterations, adjust as needed
#   # early_stopping_rounds <- max(50, floor(iteration * early_stopping_fraction))  # At least 50 rounds
#   #
#   # # Perform k-fold cross-validation with early stopping
#   # cv_results <- xgboost::xgb.cv(
#   #   params = xgb_params,
#   #   data = geno_omic_object_train,
#   #   nrounds = iteration,
#   #   nfold = n_folds,
#   #   early_stopping_rounds = early_stopping_rounds,
#   #   maximize = FALSE, # or TRUE depending on your metric
#   #   verbose = 0
#   # )
#   #
#   # # Get the best number of rounds
#   # best_iteration <- cv_results$best_iteration
#   #
#   # # Train the final model using the best number of rounds
#   # xgb_fit <- xgboost::xgb.train(
#   #   params = xgb_params,
#   #   data = geno_omic_object_train,
#   #   nrounds = best_iteration,
#   #   watchlist = list(train = geno_omic_object_train),
#   #   verbose = 0
#   # )
#   pheno_object <- pheno_object |>
#     dplyr::mutate(bin = dplyr::ntile(!!dplyr::sym(response), 10))  # Creating 10 bins based on quantiles
#
#   #set.seed(123)
#   train_indices <- caret::createDataPartition(pheno_object$bin, p = 0.8, list = FALSE)
#   training <- pheno_object[train_indices, ]
#   validation <- pheno_object[-train_indices, ]
#
#   train_geno <- geno_omic_object[rownames(geno_omic_object)%in%training[[gen_name]], ]
#   val_geno <- geno_omic_object[rownames(geno_omic_object)%in%validation[[gen_name]], ]
#
#   early_stopping_fraction <- 0.2  # 10% of iterations
#   early_stopping_rounds <- max(50, floor(iteration * early_stopping_fraction))  # At least 50 rounds
#
#   # Creating DMatrix objects
#   train_dmatrix <- xgboost::xgb.DMatrix(data = train_geno, label = training[[response]])
#   eval_dmatrix <- xgboost::xgb.DMatrix(data = val_geno, label = validation[[response]])
#
#
#   watchlist <- list(train = train_dmatrix, eval = eval_dmatrix)
#
#   xgb_fit <- xgboost::xgb.train(
#     params = xgb_params,
#     data = geno_omic_object_train,
#     nrounds = iteration,
#     early_stopping_rounds = early_stopping_rounds,
#     watchlist = watchlist,
#     maximize = FALSE, ##since these are eval_metric = c("rmse", "rmsle", "mape")
#     verbose = 0
#   )
#
# } else {
#
# xgb_fit <- xgboost::xgb.train(
#   params = xgb_params,
#   data = geno_omic_object_train,
#   nrounds = iteration,
#   #early_stopping_rounds = early_stopping_rounds,
#   #watchlist = watchlist,
#   #maximize = FALSE, ##since these areeval_metric = c("rmse", "rmsle", "mape")
#   verbose = 0
# )
# }


# if(!is.null(geno_omic_test_object)){
#   GID <- rownames(geno_omic_test_object)
#   AI_preds <- stats::predict(xgb_fit,
#                            as.matrix(geno_omic_test_object),
#                            reshape = TRUE)
#
# } else {
#   if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){
#
#     AI_preds <- stats::predict(xgb_fit,
#                                 as.matrix(geno_omic_object),
#                                 reshape = TRUE)
#
#
#   }
#
# }
