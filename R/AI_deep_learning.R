

generate_dynamic_layers <- function(input_size, num_hidden_layers, scaling_factor = 0.5, max_neurons = 1000) {

  msg <- "\n==================================================\n"

  if (length(num_hidden_layers) > 1 || length(num_hidden_layers) == 0) {
    stop(paste(msg, "num_hidden_layers should be a vector of length 1", call. = FALSE))
  }

  layers <- numeric(num_hidden_layers)

  # Cap the input size to prevent large layers
  capped_input_size <- ifelse(input_size >= 1000, max_neurons, input_size)

  # The first hidden layer could start as a fraction of the (capped) input size
  layers[1] <- min(floor(capped_input_size * scaling_factor), max_neurons)

  # Subsequent layers reduce in size, following the scaling factor and capped by max_neurons
  for (i in 2:num_hidden_layers) {
    layers[i] <- min(floor(layers[i - 1] * scaling_factor), max_neurons)
    if (layers[i] < 1) break  # Stop adding layers if neurons fall below 1
  }

  return(layers[layers > 0])  # Return only valid layers
}


# Function to train and predict using dp for bootstrapping
train_predict_deeplearning <- function(data_label_geno, indices, test_geno,
                                       num_hidden_layers = NULL,
                                       neurons_per_layer = NULL,
                                       learning_rate = 0.001,
                                       epochs = 10,
                                       batch_size = 32,
                                       validation_split = 0.2,
                                       l2_regularizer_dp = 0.001,
                                       dropout_rate = 0.5,
                                       n_blocks = 2,
                                       dense_layers_cnn = c(128, 64),
                                       kernel_size = 3,
                                       n_neurons_per_block = NULL,
                                       deep_learning_model = "mlp_with_attention",
                                       attention_on_final_layer = TRUE,
                                       attention_across_multiple_layers = FALSE,
                                       batch_normalization = TRUE
                                       ) {
  train_data <- data_label_geno[, -1]
  y_train <- data_label_geno[, 1]
  # Subset the data
  model <- deep_learning_model_utilityy(X_train = train_data,
                                        y_train = y_train,
                                        num_hidden_layers = num_hidden_layers,
                                        neurons_per_layer = neurons_per_layer,
                                        learning_rate = learning_rate,
                                        epochs = epochs,
                                        batch_size = batch_size,
                                        dropout_rate = dropout_rate,
                                        l2_regularizer_dp = l2_regularizer_dp,
                                        validation_split = validation_split,
                                        n_blocks = n_blocks,
                                        kernel_size = kernel_size,
                                        dense_layers_cnn = dense_layers_cnn,
                                        n_neurons_per_block = n_neurons_per_block,
                                        deep_learning_model = deep_learning_model,
                                        attention_on_final_layer = attention_on_final_layer,
                                        attention_across_multiple_layers = attention_across_multiple_layers,
                                        batch_normalization = batch_normalization,
                                        para_tunning = FALSE)
  # Predict on the original data
  if(!is.null(test_geno)){
    pred <-  model$predict(test_geno)
  } else{
    pred <-  model$predict(train_data)
  }
  return(pred)
}

get_best_model <- function(results, hyperparam_combinations) {
  best_metric_value <- Inf  # Start with the highest possible error for regression or lowest possible accuracy for classification
  best_model <- NULL
  best_hyperparameters <- NULL

  for (i in seq_along(results)) {
    key <- names(results)[i]
    history <- results[[key]]$history

    # Determine which metric to use based on the available keys in history
    if ("val_accuracy" %in% names(history$history)) {
      metric_name <- "val_accuracy"
      comparison_function <- max
      best_metric_value <- 0  # Initialize to the lowest possible accuracy
    } else if ("val_mean_absolute_error" %in% names(history$history)) {
      metric_name <- "val_mean_absolute_error"
      comparison_function <- min
      best_metric_value <- Inf  # Initialize to the highest possible error
    } else {
      warning("No recognized metric found in the model history.")
      next
    }

    # Extract the metric value
    current_metric_value <- comparison_function(history$history[[metric_name]])

    # Compare and store the best model and hyperparameters
    if ((metric_name == "val_accuracy" && current_metric_value > best_metric_value) ||
        (metric_name == "val_mean_absolute_error" && current_metric_value < best_metric_value)) {
      best_metric_value <- current_metric_value
      best_model <- results[[key]]$model
      best_hyperparameters <- hyperparam_combinations[i, ]  # Extract the corresponding hyperparameters
    }
  }

  return(list(best_model = best_model, best_hyperparameters = best_hyperparameters))
}


#' Title
#'
#' @param num_hidden_layers
#' @param neurons_per_layer
#' @param learning_rate
#' @param epochs
#' @param param_grid
#' @param validation_split
#' @param early_stop
#' @param batch_size
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param gen_name
#' @param message
#' @param scale
#' @param para_tunning
#'
#' @return
#' @export
#'
#' @examples
deep_learning_model <- function(pheno_object=NULL,
                                y = NULL,
                                omics_data = NULL,
                                crossval = FALSE,
                                tst = NULL,
                                geno_omic_object = NULL,
                                geno_omic_test_object = NULL,
                                response=NULL,
                                gen_name=NULL,
                                message = TRUE,
                                scaling = TRUE,
                                centering = FALSE,
                                omic_count = NULL,
                                num_hidden_layers = 1,
                                neurons_per_layer = NULL,
                                learning_rate_dp = 0.001,
                                epochs = 10,
                                batch_size = 32 ,
                                para_tunning = FALSE,
                                param_grid = NULL,
                                validation_split = 0.2,
                                early_stop = TRUE,
                                l2_regularizer_dp = 0.001,
                                dropout_rate = 0.5,
                                deep_learning_model = "mlp_with_attention",
                                n_blocks = 2,
                                dense_layers_cnn = c(128, 64),
                                kernel_size = 3,
                                n_neurons_per_block = NULL,
                                attention_on_final_layer = TRUE,
                                attention_across_multiple_layers = FALSE,
                                batch_normalization = TRUE,
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
#browser()
  msg <- "\n==================================================\n"
  if(is.null(dense_layers_cnn)) dense_layers_cnn <-  64
  if(!is.null(geno_omic_object)){
    scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)

  }

  if(!is.null(geno_omic_test_object)){
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)

  }

  if(!is.null(omics_data)){
    scaler <- caret::preProcess(omics_data, method = c("center", "scale"))
    omics_data <- stats::predict(scaler, omics_data)
    omics_data <- stats::predict(scaler, omics_data)

  }

  #"ResNet",
  if(any(deep_learning_model%in%c("mlp_with_attention", "mlp", "cnn"))){
    #if(is.null(neurons_per_layer)&& is.null(num_hidden_layers)){
    if(is.null(num_hidden_layers)){
      stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
    }
  }

  if(deep_learning_model=="ResNet"){
    #if(is.null(n_neurons_per_block) || is.null(n_blocks)){
    if(is.null(n_blocks)){
      stop(paste(msg, "n_blocks can't be NULL."), call. = FALSE)
    }
  }

  if(is.null(neurons_per_layer)&& !is.null(num_hidden_layers)){
    if(!is.null(geno_omic_object) && isFALSE(crossval)){
    neurons_per_layer <- generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
    } else if(!is.null(omics_data) && isTRUE(crossval)){
      neurons_per_layer <- generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)

    } else {
      stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
}
  }



  if(deep_learning_model == "ResNet"){
    if(is.null(n_blocks)) n_blocks <- 1
    #if(is.null(dense_layers_cnn)) dense_layers_cnn <-  64
      if(is.null(n_neurons_per_block)&& !is.null(n_blocks)){
        if(!is.null(geno_omic_object) && isFALSE(crossval)){
        n_neurons_per_block <-  generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = n_blocks, scaling_factor = 0.5)

        } else if(!is.null(omics_data) && isTRUE(crossval)){
          n_neurons_per_block <-  generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = n_blocks, scaling_factor = 0.5)
        } else {
          stop(paste(msg, "n_neurons_per_block and n_blocks can't be NULL."), call. = FALSE)
        }
      }

    if(n_blocks!=length(n_neurons_per_block)){
      stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
    }

    neurons_per_layer <- n_neurons_per_block
    num_hidden_layers <- n_blocks
  } else {
    if(is.null(num_hidden_layers)) num_hidden_layers <- 1
    if(is.null(neurons_per_layer)&& !is.null(num_hidden_layers)){
      if(!is.null(geno_omic_object) && isFALSE(crossval)){
        neurons_per_layer <-  generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)

      } else if(!is.null(omics_data) && isTRUE(crossval)){
        neurons_per_layer <-  generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
      } else {
        stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
      }
    }
  }

  if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers) else num_hidden_layers <- as.integer(1)
  if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.numeric(neurons_per_layer) else neurons_per_layer <- as.numeric(64)
  if(!is.null(batch_size))   batch_size <- as.integer(batch_size) else  batch_size <- as.integer(30)
  if(!is.null(epochs)) epochs <- as.integer(epochs) else epochs <- as.integer(10)
  if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate) else dropout_rate <- as.numeric(0.5)
  if(!is.null(learning_rate_dp)) learning_rate <- as.numeric(learning_rate_dp) else learning_rate <- as.numeric(0.01)
  if(!is.null(l2_regularizer_dp)) l2_regularizer_dp <- as.integer(l2_regularizer_dp)

  if(!is.null(geno_omic_test_object)){
    GID <- rownames(geno_omic_test_object)
  } else {
    if(!is.null(geno_omic_object)){
      GID <- rownames(geno_omic_object)
    }
  }

   np <- reticulate::import("numpy")

   if(isFALSE(crossval)){
  y_train <- as.numeric(pheno_object[, response])
  # y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
  #
  # # Predict on the training data and get the scaled values
  # y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

   } else{
     tryCatch({
       if(isTRUE(crossval)){
         y_train <- y
         y_train = y_train[-tst]
         y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

         # Predict on the training data and get the scaled values
         y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]

         model_dp <- deep_learning_model_utilityy(X_train = omics_data[-tst, ],
                                               #y_train = y_train[-tst],
                                               y_train = y_train,
                                               num_hidden_layers = num_hidden_layers,
                                               neurons_per_layer = neurons_per_layer,
                                               learning_rate = learning_rate,
                                               epochs = epochs,
                                               batch_size = batch_size,
                                               dropout_rate = dropout_rate,
                                               l2_regularizer_dp = l2_regularizer_dp,
                                               validation_split = validation_split,
                                               n_blocks = num_hidden_layers,
                                               n_neurons_per_block = neurons_per_layer,
                                               kernel_size = kernel_size,
                                               dense_layers_cnn = dense_layers_cnn,
                                               deep_learning_model = deep_learning_model,
                                               attention_on_final_layer = attention_on_final_layer,
                                               attention_across_multiple_layers = attention_across_multiple_layers,
                                               batch_normalization = batch_normalization,
                                               para_tunning = FALSE)

         preds <- model_dp$predict(omics_data[tst, ])
         preds <- as.data.frame(preds)
         rm(model_dp)
         return(preds[, 1])
       }
     }, error = function(e) {
       message("An error occurred: ", e$message)
       return(NULL) # or handle the error as needed
     })
   }

  # Scale the training labels
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]


  data_label_geno <- cbind(y_train_scaled,geno_omic_object)

  if (isTRUE(para_tunning) && is.null(param_grid)) {
    stop(paste(msg, "param_grid must be provided when tuning is enabled."), call. = FALSE)
  }

  if(is.null(neurons_per_layer)) neurons_per_layer <- list(ncol(geno_omic_object)/2)

  # Create and compile the model


  if (isTRUE(para_tunning)) {
    # Perform grid search for hyperparameter tuning
    # Generate all combinations of hyperparameters
    hyperparam_combinations <- expand.grid(param_grid)

    # Function to check if the number of hidden layers matches the length of neurons_per_layer
    check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
      neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
      return(length(neuron_list) == num_hidden_layers)
    }

    # Loop through the data frame and remove rows with mismatched hidden layers and neurons
    valid_rows <- apply(hyperparam_combinations, 1, function(row) {
      num_hidden_layers <- as.integer(row["num_hidden_layers"])
      neurons_per_layer <- row["neurons_per_layer"]
      check_hidden_layers(num_hidden_layers, neurons_per_layer)
    })

    # Filter the data frame to keep only valid rows
    hyperparam_combinations <- hyperparam_combinations[valid_rows, ]

    grid_results <- grid_search_deep_learning(X_train = geno_omic_object,
                                             y_train = y_train_scaled,
                                             param_grid = param_grid,
                                             epochs = epochs,
                                             batch_size = batch_size,
                                             validation_split = validation_split,
                                             early_stop = early_stop,
                                             n_blocks = n_blocks,
                                             n_neurons_per_block = n_neurons_per_block,
                                             dense_layers_cnn = dense_layers_cnn,
                                             kernel_size = kernel_size,
                                             learning_rate = learning_rate_dp,
                                             dropout_rate = dropout_rate,
                                             l2_regularizer_dp = l2_regularizer_dp,
                                             para_tunning = para_tunning)



    best_model <- get_best_model(grid_results, hyperparam_combinations)
    best_hyperparameter <- best_model$best_hyperparameters

    # Initialize variables only if they exist in the hyperparams
    if ("num_hidden_layers" %in% names(best_hyperparameter )) {
      num_hidden_layers <- as.numeric(best_hyperparameter $num_hidden_layers)
    }

    if ("neurons_per_layer" %in% names(best_hyperparameter )) {
      if (inherits(best_hyperparameter $neurons_per_layer, "list")) {
        neurons_per_layer <- unlist(best_hyperparameter $neurons_per_layer)
      } else {
        neurons_per_layer <- best_hyperparameter$neurons_per_layer
      }
    }

    if ("learning_rate" %in% names(best_hyperparameter )) {
      learning_rate <- best_hyperparameter $learning_rate_dp
    }

    if ("dropout_rate" %in% names(best_hyperparameter )) {
      dropout_rate <- best_hyperparameter $dropout_rate
    }

    if ("epochs" %in% names(best_hyperparameter)) {
      epochs <- best_hyperparameter$epochs
    }

    if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers) else num_hidden_layers <- as.integer(1)
    if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.numeric(neurons_per_layer) else neurons_per_layer <- as.numeric(64)
    if(!is.null(batch_size))   batch_size <- as.integer(batch_size) else  batch_size <- as.integer(30)
    if(!is.null(epochs)) epochs <- as.integer(epochs) else epochs <- as.integer(10)
    if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate) else dropout_rate <- as.numeric(0.5)
    if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate) else learning_rate <- as.numeric(0.01)

}

    if (!is.null(geno_omic_test_object)) {
      X_test <- np$array(as.matrix(geno_omic_test_object))
      #predictions <- best_model$best_model$predict(X_test)

      boot_results <- boot::boot(
                                data = data_label_geno,
                                statistic = train_predict_deeplearning,
                                num_hidden_layers = num_hidden_layers,
                                neurons_per_layer = neurons_per_layer,
                                learning_rate = learning_rate,
                                epochs = epochs,
                                batch_size = batch_size,
                                l2_regularizer_dp = l2_regularizer_dp,
                                dropout_rate = dropout_rate,
                                validation_split = validation_split,
                                R = n_bootstrap,  # Number of bootstrap samples
                                #sim = "ordinary",
                                test_geno =X_test
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
                                                          #composite_reliability_score = composite_reliability$reliability_score,
                                                          #composite_reliability = composite_reliability$trustworthiness,
                                                          #composite_reliability_percentage = composite_reliability$reliability_percentage,
                                                          #threshold = NULL,
                                                          high_reliability_thres = high_reliability_thres,
                                                          low_reliability_thres = low_reliability_thres,
                                                          system_database = system_database)


     # return(list(predicted_values = predictions, trained_model = best_model))
    } else {

      boot_results <- boot::boot(
        data = data_label_geno,
        statistic = train_predict_deeplearning,
        num_hidden_layers = num_hidden_layers,
        neurons_per_layer = neurons_per_layer,
        learning_rate = learning_rate,
        epochs = epochs,
        batch_size = batch_size,
        l2_regularizer_dp = l2_regularizer_dp,
        dropout_rate = dropout_rate,
        R = n_bootstrap,  # Number of bootstrap samples
        #sim = "ordinary",
        test_geno =NULL
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
                                                          #composite_reliability_score = composite_reliability$reliability_score,
                                                          #composite_reliability = composite_reliability$trustworthiness,
                                                          #composite_reliability_percentage = composite_reliability$reliability_percentage,
                                                          #threshold = NULL,
                                                          high_reliability_thres = high_reliability_thres,
                                                          low_reliability_thres = low_reliability_thres,
                                                          system_database = system_database)
    }
# browser()
# print(num_hidden_layers)
# print(learning_rate)
# print(neurons_per_layer)
  # model_para <- data.frame(stat = c("num_hidden_layers","learning_rate", "neurons_per_layer"),
  #                          summary = c(num_hidden_layers, learning_rate, neurons_per_layer),
  #                          stringsAsFactors = FALSE)

model_para <- data.frame(
                      stat = c("num_hidden_layers", "learning_rate", "neurons_per_layer"),
                      summary = I(list(num_hidden_layers, learning_rate, neurons_per_layer)),
                      stringsAsFactors = FALSE
                    )

if(length(neurons_per_layer)>1){
  # Flatten the list in the 'summary' column and replicate 'stat' values accordingly
  expanded_stat <- unlist(lapply(1:nrow(model_para), function(i) {
    rep(model_para$stat[i], length(model_para$summary[[i]]))
  }))

  expanded_summary <- unlist(model_para$summary)

  # Create a new data frame with the expanded values
  model_para <- data.frame(stat = expanded_stat, summary = expanded_summary, stringsAsFactors = FALSE)

}

  colnames(model_para)[1:2] <- c("stat", "summary")

  output <-  list(model_parameters = model_para,
                  predicted_values = AI_preds,
                  diagnostic_plots = diagnostic_plots)


  return(output)

}


