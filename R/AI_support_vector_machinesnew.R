select_best_model <- function(model_list,
                              preferred_metrics=c("RMSE", "Accuracy", "MAE")) {

  scores <- lapply(model_list, function(model) {
    available_metrics <- intersect(names(model$results), preferred_metrics)
    if (length(available_metrics) > 0) {
      # Assuming lower values are better for selected metrics, flip if needed
      metric_scores <- sapply(available_metrics, function(metric) {
        if (metric == "Accuracy") {
          max(model$results[[metric]])  # Higher is better
        } else if(metric == "MAE"){
          min(model$results[[metric]])
        }else {
          if(metric == "RMSE"){
            min(model$results[[metric]])
          }
        }
        # Lower is better for RMSE
      })
      c(BestScore=min(metric_scores), BestMetric=names(metric_scores)[which.min(metric_scores)])
    } else {
      c(BestScore=NA, BestMetric=NA)
    }
  })

  # Remove models without any of the preferred metrics
  scores <- scores[!sapply(scores, function(x) is.na(x["BestScore"]))]

  # Identify the model with the best overall score
  best_model_index <- which.min(sapply(scores, function(x) x["BestScore"]))
  best_model_name <- names(scores)[best_model_index]

  list(BestModel=model_list[[best_model_name]], BestModelName=best_model_name, Metrics=scores[[best_model_index]])
}


AI_svm <- function(pheno_object=NULL,
                    geno_omic_object=NULL,
                    geno_omic_test_object=NULL,
                    response=NULL,
                    gen_name=NULL,
                    message=TRUE,
                    scaling=TRUE,
                    centering=FALSE,
                    omic_count=NULL,
                    para_tunning=FALSE,
                    AI_cv_nfolds=5,
                    svm_paras_tunning=list(
                    kernel = c("Gaussian", "Linear",
                               "Polynomial", "Hyperbolic_tangent"),
                     #cost = 10^seq(-2, 2, by = 1),
                     offset_value = seq(-2, 2, length.out = 5),
                     sigma = c(0.01, 0.05, 0.1),
                     C = c(1, 10, 100), ## for radial kernel
                     gamma_value = 10^seq(-4, -1, length.out = 4),
                     degree = c(3, 4),  # Default values, used only for polynomial
                     scale = c(0.1, 1) # used only for polynomial
                   ),
                   svm_kernel="Gaussian",
                   svm_type = "eps-regression",
                   sigma_value=0.1,
                   C_value=1,
                   degree_value=3,
                   scale_value=1,
                   offset_value=0,
                   gamma_value = NULL,
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
                   ...) {

  msg <- sprintf("==================================================\n")

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


  data <- cbind(y=y_train_scaled, geno_omic_object)

  # Map user-friendly kernel names to e1071 kernel types
  kernel_type <- switch(svm_kernel,
                        Gaussian="radial",
                        Polynomial="polynomial",
                        Linear="linear",
                        Hyperbolic_tangent="sigmoid")

  # Define the tuning grid and kernel parameters
  #gamma_value <- 1 / (2 * sigma_value^2)  # Convert sigma to gamma

  kernel_params <- list(degree=degree_value,
                        #gamma=gamma_value,
                        coef0=offset_value,
                        cost=C_value
                        #scale=scale_value
                        )

  # Perform scaling if required
  # if (isTRUE(scaling)) {
  #   if (!is.null(geno_omic_object)) {
  #     geno_omic_object <- scale(geno_omic_object, center=TRUE, scale=TRUE)
  #   }
  #   if (!is.null(geno_omic_test_object)) {
  #     geno_omic_test_object <- scale(geno_omic_test_object, center=TRUE, scale=TRUE)
  #   }
  # }

  # Parameter tuning and model training
  if (isTRUE(para_tunning)) {
    # Setup for caret to handle SVM training and cross-validation
    # train_control <- caret::trainControl(method="cv", number=AI_cv_nfolds,
    #                               verboseIter=TRUE, returnData=FALSE,
    #                               returnResamp="all", allowParallel=TRUE)
    #
    # # Define the tuning grid based on the svm_paras_tunning structure
    # tuning_grid <- expand.grid(kernel=svm_paras_tunning$kernel,
    #                            degree=svm_paras_tunning$degree,
    #                            #scale=svm_paras_tunning$scale,
    #                            C=svm_paras_tunning$C
    #                            #gamma=1 / (2 * svm_paras_tunning$sigma^2)
    #                            )  # Convert sigma to gamma

    # # Train the SVM model using caret for parameter tuning
    # svm_model <- caret::train(x=geno_omic_object, y=pheno_object[, response],
    #                           method="svmRadial",  # Change depending on kernel choice
    #                           trControl=train_control, tuneGrid=tuning_grid)
    # create hyperparameter grid

    # Include gamma for kernels that need it (not needed for linear unless specified)
    if ("Gaussian" %in% svm_paras_tunning$kernel) {
      svm_paras_tunning$kernel <- "radial"
      tuning_grid_radial <- expand.grid(sigma = svm_paras_tunning$sigma,
                                        C = svm_paras_tunning$C)
    }
    ##
    if ("Polynomial" %in% svm_paras_tunning$kernel) {
      svm_paras_tunning$kernel <- "polynomial"
      tuning_grid_poly <- expand.grid(degree = svm_paras_tunning$degree,
                                      #scale = svm_paras_tunning$scale,
                                      C = svm_paras_tunning$C)
    }
    ##
    if ("Linear" %in% svm_paras_tunning$kernel) {
      svm_paras_tunning$kernel <- "linear"
      tuning_grid_linear <- expand.grid(C = data.frame(C =svm_paras_tunning$C))
    }

    if ("Hyperbolic_tangent" %in% svm_paras_tunning$kernel) {
      svm_paras_tunning$kernel <- "sigmoid"
      tuning_grid_linear <- expand.grid(C = data.frame(C =svm_paras_tunning$C),
                                        gamma = data.frame(C =svm_paras_tunning$gamma_value),
                                        coef0 = data.frame(C =svm_paras_tunning$offset_value))
    }
    ###
    # List of SVM models to evaluate
    svm_models <- list(
      svmRadial = list(method = "svmRadial", tuneGrid = tuning_grid_radial),
      svmPoly = list(method = "svmPoly", tuneGrid = tuning_grid_poly),
      svmLinear = list(method = "svmLinear", tuneGrid = tuning_grid_linear)
      #svmSigmoid = list(method = "svmSigmoid", tuneGrid = tuning_grid_linear)
    )

    AI_trcontrol = caret::trainControl(method = "cv",
                                       number = AI_cv_nfolds,
                                       verboseIter = TRUE,
                                       returnData = FALSE,
                                       returnResamp = "all",
                                       allowParallel = TRUE)

    # Train the SVM model using caret for parameter tuning
    # Function to train and evaluate models
    train_and_evaluate <- function(model_info, data, labels, control) {
      caret::train(
        x = data,
        y = labels,
        method = model_info$method,
        trControl = control,
        tuneGrid = model_info$tuneGrid
        #preProcess = "scale"  # scale and center predictors

      )
    }

      results <- lapply(svm_models,
                        train_and_evaluate,
                        data = geno_omic_object,
                        labels = pheno_object[, response],
                        control = AI_trcontrol)

      # Apply the function to the results
      best_model_info <- select_best_model(results)

      best_kernel <- best_model_info$BestModelName
      ##
      best_C <- NULL
      best_sigma <- NULL
      best_degree <- NULL
      #best_scale <- NULL
      best_gamma <- NULL
      best_coef0 <- NULL



      if("svmLinear" %in% best_model_info$BestModelName){
        best_C = best_model_info$BestModel$bestTune$C
        best_kernel <- "linear"
      } else if("svmRadial" %in% best_model_info$BestModelName){
        best_kernel <-"radial"
        best_C = best_model_info$BestModel$bestTune$C
        best_sigma = best_model_info$BestModel$bestTune$sigma
      } else if ("svmPoly" %in% best_model_info$BestModelName){
        best_kernel <-"polynomial"
        best_C = best_model_info$BestModel$bestTune$C
        best_degree = best_model_info$BestModel$bestTune$degree
        #best_scale = best_model_info$BestModel$bestTune$scale

      } else {
        if ("sigmoid" %in% best_model_info$BestModelName){
          best_C = best_model_info$BestModel$bestTune$C
          best_gamma = best_model_info$BestModel$bestTune$gamma
          best_coef0 = best_model_info$BestModel$bestTune$coef0

        }
      }
      #svm_model <- best_model_info$BestModel

      # Function to fit SVM and make predictions
      svm_predict_boost <- function(data, indices,
                                    geno_test ,
                                    best_C,
                                    #best_sigma,
                                    best_degree,
                                    #best_scale,
                                    best_gamma,
                                    best_coef0,
                                    best_kernel,
                                    svm_type
                                    ) {
        # geno <- data$geno
        # pheno <- data$pheno
        # geno_test <- data$geno_test
        # best_C <- data$best_C
        # best_sigma <- data$best_sigma
        # best_degree <- data$best_degree
        # best_scale <- data$best_scale
        # best_gamma <- data$best_gamma
        # best_coef0 <- data$best_coef0
        # best_kernel <- data$best_kernel
        # response <- data$response
        x_train <- data[indices, -1]  # ensure 'drop = FALSE' to keep the data frame structure if one column
        y_train <- data[indices, 1]

        if("linear"%in%best_kernel){
          svm_model <- e1071::svm(x = x_train, y = y_train,  kernel = best_kernel, cost = best_C, type = svm_type, scale = FALSE)
        } else if("radial" %in% best_kernel){
          svm_model <- e1071::svm(x = x_train, y = y_train,  kernel = best_kernel, cost = best_C, type = svm_type, scale = FALSE
                                  #sigma =  best_sigma
                                  )
        } else if("polynomial" %in% best_kernel){
          svm_model <- e1071::svm(x = x_train, y = y_train,  kernel = best_kernel, cost = best_C,  degree = best_degree,
                                  #scale = best_scale,
                                  scale = FALSE,
                                  type = svm_type)
        } else if("sigmoid" %in% best_kernel){
          svm_model <- e1071::svm(x = x_train, y = y_train,  kernel = best_kernel,
                                  cost = best_C,  gamma = best_gamma,
                                  coef0 = best_coef0, type = svm_type, scale = FALSE)
        }
        if(!is.null(geno_test)){
          predictions <- stats::predict(svm_model, as.matrix(geno_test))  # predicting on the training data itself for bootstrap

        } else {
          predictions <- stats::predict(svm_model, as.matrix(data[, -1]))  # predicting on the training data itself for bootstrap
        }
        return(predictions)
      }

      # data_list <- list(
      #   geno = geno_omic_object,
      #   pheno = pheno_object,
      #   geno_test = geno_omic_test_object,
      #   best_C = best_C,
      #   best_sigma = best_sigma,
      #   best_degree = best_degree,
      #   best_scale = best_scale,
      #   best_gamma = best_gamma,
      #   best_coef0 = best_coef0,
      #   best_kernel = best_kernel,
      #   response = response  # Make sure to define which column is your response
      # )
      # Bootstrapping
      boot_results <- boot::boot(data = data, statistic = svm_predict_boost,
                                 R = n_bootstrap,
                                 geno_test = geno_omic_test_object,
                                 best_C = best_C,
                                 #best_sigma = best_sigma,
                                 best_degree = best_degree,
                                 #best_scale = best_scale,
                                 best_gamma = best_gamma,
                                 best_coef0 = best_coef0,
                                 best_kernel = best_kernel,
                                 #response = response,
                                 svm_type = svm_type)
      # Calculate standard error of the predictions
      # prediction_se <- apply(boot_results$t, 2, sd)  # Assuming predictions are in rows
      # PEV <- apply(boot_results$t, 2, var)
      #
      # # Standard errors of predictions
      # predictions <- apply(boot_results$t, 2, mean)

      # Train SVM without parameter tuning
  } else {

    # Train SVM without parameter tuning
    para_index <- which(!sapply(kernel_params, is.null))
    if(length(para_index)!=0){
      kernel_params <-  kernel_params[para_index]
      # para_names <- names(kernel_params)
      # for(name in names(kernel_params)) {
      #   assign(name, kernel_params[[name]], envir = environment())
      # }

    } else {
      kernel_params <- list()
      kernel_params[c("degree", "coef0", "cost")] <- c(3, 0, 1)

      # cost <- 1
      # degree <- 3
      # coef0 <- 0
    }

    para_names <- names(kernel_params)

    # Ensure default cost is set if not already specified
    if (!"cost" %in% para_names) {
      kernel_params$cost <- 1
    }


    svm_predict_boost <- function(data, indices,
                                  geno_test ,
                                  kernel_params,
                                  kernel_type,
                                  svm_type,
                                  para_names
                                  ) {

      x_train <- data[indices, -1]  # ensure 'drop = FALSE' to keep the data frame structure if one column
      y_train <- data[indices, 1]

      if ("linear" %in% kernel_type) {
        svm_model <- e1071::svm(x = x_train, y = y_train, kernel = kernel_type, cost = kernel_params$cost, type =  svm_type, scale = FALSE)
      } else if ("radial" %in% kernel_type) {
        if (!"gamma" %in% para_names) {
          svm_model <- e1071::svm(x = x_train, y = y_train, cost = kernel_params$cost, type =  svm_type, scale = FALSE)
        } else {
          svm_model <- e1071::svm(x = x_train, y = y_train, cost = kernel_params$cost,type =  svm_type,
                                  gamma = gamma, scale = FALSE)
        }
      } else if ("polynomial" %in% kernel_type) {
        if (all(c("scale", "degree") %in% para_names)) {
          svm_model <- e1071::svm(x = x_train, y = y_train, cost = kernel_params$cost, degree = kernel_params$degree,
                                  #scale = scale,
                                  scale = FALSE,
                                  type =  svm_type)
        } else {
          svm_model <- e1071::svm(x = x_train, y = y_train, type =  svm_type, scale = FALSE)
        }
      } else if ("sigmoid" %in% kernel_type) {
        if ("coef0" %in% para_names) {
          svm_model <- e1071::svm(x = x_train, y = y_train, cost = kernel_params$cost,
                                  coef0 = kernel_params$coef0, type =  svm_type, scale = FALSE)
        } else {
          svm_model <- e1071::svm(x = x_train, y = y_train, type =  svm_type, scale = FALSE)
        }
      }
      if(!is.null(geno_test)){
        predictions <- stats::predict(svm_model, geno_test)  # predicting on the testing data for bootstrap

      } else {
        predictions <- stats::predict(svm_model, data[, -1])  # predicting on the training data itself for bootstrap
      }
      return(predictions)
    }

    if(!is.null(geno_omic_test_object)){

    boot_results <- boot::boot(data = data, statistic = svm_predict_boost,
                               R = n_bootstrap,
                               geno_test = geno_omic_test_object,
                               kernel_params = kernel_params,
                               kernel_type = kernel_type,
                               svm_type = svm_type,
                               para_names = para_names
                               )

    } else {

      boot_results <- boot::boot(data = data, statistic = svm_predict_boost,
                                 R = n_bootstrap,
                                 geno_test = geno_omic_object,
                                 kernel_params = kernel_params,
                                 kernel_type = kernel_type,
                                 svm_type = svm_type,
                                 para_names = para_names
      )

    }
    # Calculate standard error of the predictions
    # prediction_se <- apply(boot_results$t, 2, sd)  # Assuming predictions are in rows
    # PEV <- apply(boot_results$t, 2, var)
    #
    # # Standard errors of predictions
    # predictions <- apply(boot_results$t, 2, mean)

    # Define the SVM model based on the kernel type
    # if ("linear" %in% kernel_type) {
    #   svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response], kernel = kernel_type, cost = kernel_params$cost)
    # } else if ("radial" %in% kernel_type) {
    #   if (!"gamma" %in% para_names) {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response], cost = kernel_params$cost)
    #   } else {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response], cost = kernel_params$cost,
    #                             gamma = kernel_params$gamma)
    #   }
    # } else if ("polynomial" %in% kernel_type) {
    #   if (all(c("scale", "degree") %in% para_names)) {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response], cost = kernel_params$cost, degree = kernel_params$degree, scale = kernel_params$scale)
    #   } else {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response])
    #   }
    # } else if ("sigmoid" %in% kernel_type) {
    #   if ("coef0" %in% para_names) {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response], cost = kernel_params$cost, coef0 = kernel_params$coef0)
    #   } else {
    #     svm_model <- e1071::svm(x = geno_omic_object, y = pheno_object[, response])
    #   }
    # }

    # if (!is.null(geno_omic_test_object)) {
    #   predictions <- stats::predict(svm_model, geno_omic_test_object)
    #   name_source <- rownames(geno_omic_test_object)
    # } else {
    #   predictions <- stats::predict(svm_model, geno_omic_object)
    #   name_source <- rownames(geno_omic_object)
    # }
  }

  # Predict using the trained model
  if (!is.null(geno_omic_test_object)) {
    #predictions <- stats::predict(svm_model, geno_omic_test_object)
    name_source <- rownames(geno_omic_test_object)
  } else {
    #predictions <- stats::predict(svm_model, geno_omic_object)
    name_source <- rownames(geno_omic_object)
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

  AI_preds <- data.frame(name = name_source,
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

  #####

  diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = boot_results,
                                                      GID_names = name_source,
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

  return(list(model_parameters= kernel_params,
              predicted_values=AI_preds,
              diagnostic_plots = diagnostic_plots))
}
