
# Utility function to get the best model from grid search results
#' Title
#'
#' @param results
#' @param param_grid
#'
#' @return
#' @export
#'
#' @examples
get_best_model <- function(results, param_grid) {
  best_accuracy <- 0
  best_model <- NULL

  for (key in names(results)) {
    history <- results[[key]]$history
    validation_accuracy <- max(results[[key]]$history$metrics[[2]])

    if (validation_accuracy > best_accuracy) {
      best_accuracy <- validation_accuracy
      best_model <- results[[key]]$model
      best_hyperparameters <- unlist(strsplit(key, "_")[[1]])

    }
  }

  # Ensure that the number of hyperparameters extracted matches the number specified in param_grid
  if (all(!c("epochs", "batch_size") %in% names(param_grid))) {
    names_grid <- c(names(param_grid), "epochs", "batch_size")
  } else {
    names_grid <- names(param_grid)
  }

  # if (length(best_hyperparameters) != length(unlist(param_grid))) {
  #   stop("Number of hyperparameters obtained from grid search results does not match the number specified in param_grid.")
  # }

  # Create a data frame with hyperparameter names and the best value for each parameter
  # hyperparameters_df <- data.frame(
  #   parameter = names_grid,
  #   value = sapply(strsplit(best_hyperparameters, " ")[[1]], tail, 1),
  #   stringsAsFactors = FALSE,
  #   row.names = NULL
  # )

  res <-  list(best_model = best_model,
               best_hyperparameters = best_hyperparameters)
  return(res)
}
######
# get_best_model <- function(results, param_grid) {
#   if (length(results) == 0) {
#     stop("No results provided.")
#   }
#
#   best_accuracy <- 0
#   best_model <- NULL
#   best_hyperparameters <- NULL
#
#   for (key in names(results)) {
#     if (!is.null(results[[key]]$history)) {
#       validation_accuracy <- max(results[[key]]$history$metrics[[2]])
#
#       if (validation_accuracy > best_accuracy) {
#         best_accuracy <- validation_accuracy
#         best_model <- results[[key]]$model
#         best_hyperparameters <- key
#       }
#     }
#   }
#
#   if (is.null(best_model)) {
#     stop("No valid models found in the results.")
#   }
#
#   # Extract the hyperparameters dynamically based on the param_grid provided by the user
#   hyperparameters <- strsplit(best_hyperparameters, "_")[[1]]
#   hyperparameters_names <- names(param_grid)
#   hyperparameters_df <- data.frame(matrix(ncol = length(param_grid), nrow = 1))
#   colnames(hyperparameters_df) <- hyperparameters_names
#
#   for (i in seq_along(hyperparameters_names)) {
#     hyperparameters_df[, i] <- hyperparameters[i]
#   }
#
#   return(list(best_model = best_model, best_hyperparameters = hyperparameters_df))
# }



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
                                geno_omic_object = NULL,
                                geno_omic_test_object = NULL,
                                response=NULL,
                                gen_name=NULL,
                                message = TRUE,
                                scale = TRUE,
                                num_hidden_layers = 1,
                                neurons_per_layer = NULL,
                                learning_rate = 0.001,
                                epochs = 10,
                                batch_size = 32 ,
                                para_tunning = FALSE,
                                param_grid = NULL,
                                validation_split = 0.2,
                                early_stop = TRUE
) {
  # Validate parameters
  # if (is.null(X_train) || is.null(y_train)) {
  #   stop("X_train and y_train must be provided.")
  # }

  y_train <- as.numeric(pheno_object[, response])
  if (para_tunning && is.null(param_grid)) {
    stop("param_grid must be provided when tuning is enabled.")
  }

  if(is.null(neurons_per_layer)) neurons_per_layer <- list(ncol(geno_omic_object)/2)

  # Create and compile the model


  if (isTRUE(para_tunning)) {
    # Perform grid search for hyperparameter tuning
    grid_results <- grid_search_deep_learning(geno_omic_object,
                                             y_train,
                                             param_grid,
                                             epochs,
                                             batch_size,
                                             validation_split,
                                             early_stop)

    best_model <- get_best_model(grid_results, param_grid)

    if (!is.null(geno_omic_test_object)) {
      predictions <- stats::predict(best_model$best_model, geno_omic_test_object)
      predictions <- data.frame(name = rownames(geno_omic_object),
                                Predicted_value = as.data.frame(predictions),
                                Standard_error = NA,
                                PEV = NA,
                                Reliability = NA,
                                stringsAsFactors = FALSE)
      names(predictions)[1:2] <- c(gen_name,"Predicted_value")
      # return(list(predictions = predictions,
      #             best_model = best_model,
      #             grid_results = best_model$best_hyperparameters))
      return(list(predicted_values = predictions, trained_model = best_model))
    } else {

      if (!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        predictions <- stats::predict(best_model$best_model, geno_omic_object)
        predictions <- data.frame(name = rownames(geno_omic_object),
                                  Predicted_value = as.data.frame(predictions),
                                  Standard_error = NA,
                                  PEV = NA,
                                  Reliability = NA,
                                  stringsAsFactors = FALSE)
        names(predictions)[1:2] <- c(gen_name,"Predicted_value")
      }
      return(list(predicted_values = predictions, trained_model = best_model))
    }
  } else {
    model <- deep_learning_model_utility(geno_omic_object,
                                         y_train,
                                         num_hidden_layers,
                                         neurons_per_layer,
                                         learning_rate,
                                         epochs,
                                         batch_size)
    # Use the provided model for prediction
    if (!is.null(geno_omic_test_object)) {
      predictions <- stats::predict(model, geno_omic_test_object)
      predictions <- data.frame(name = rownames(geno_omic_object),
                                Predicted_value = as.data.frame(predictions),
                                Standard_error = NA,
                                PEV = NA,
                                Reliability = NA,
                                stringsAsFactors = FALSE)
      names(predictions)[1:2] <- c(gen_name,"Predicted_value")
      return(list(predicted_values = predictions, trained_model = model))
    } else {
      if (!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        predictions <- stats::predict(model, geno_omic_object)
        predictions <- data.frame(name = rownames(geno_omic_object),
                                  Predicted_value = as.data.frame(predictions),
                                  Standard_error = NA,
                                  PEV = NA,
                                  Reliability = NA,
                                  stringsAsFactors = FALSE)
        names(predictions)[1:2] <- c(gen_name,"Predicted_value")
      }
      return(list(predicted_values = predictions, trained_model = model))
    }
  }
}


