

# Function to validate the number of hidden layers against the neurons_per_layer list
validate_layers <- function(num_hidden_layers, neurons_per_layer) {
  msg <- "\n==================================================\n"
  if (length(num_hidden_layers) != length(neurons_per_layer)) {
    stop(paste(msg, "The length of num_hidden_layers must match the length of neurons_per_layer"), call. = FALSE)
  }

  for (i in seq_along(num_hidden_layers)) {
    if (num_hidden_layers[i] != length(neurons_per_layer[[i]])) {
      stop(paste(msg, paste("Mismatch in hidden layers and neurons per layer at", num_hidden_layers[i])), call. = FALSE)
    }
  }

  return(TRUE)
}

#' Title
#'
#' @param X_train
#' @param y_train
#' @param param_grid
#' @param epochs
#' @param batch_size
#' @param validation_split
#' @param early_stop
#'
#' @return
#' @export
#'
#' @examples
grid_search_deep_learning <- function(X_train,
                                      y_train,
                                      param_grid,
                                      epochs,
                                      batch_size,
                                      validation_split = 0.2,
                                      early_stop = TRUE) {
  results <- list()

  msg <- "\n==================================================\n"

  # Check if validation_split is within valid range
  if (early_stop && (validation_split < 0 || validation_split >= 1)) {
    stop(paste(msg, "Validation split must be between 0 and 1."), call. = FALSE)
  }

  # Call the validation function
 #test_neuron_hidden_layer <-  tryCatch({
    validate_layers(param_grid$num_hidden_layers, param_grid$neurons_per_layer)
  #   cat("Validation successful: The number of hidden layers matches the neurons per layer configuration.\n")
  # }, error = function(e) {
  #   cat("Validation error:", e$message, "\n")
  #   return(NULL)
  # })

 # if(is.null(test_neuron_hidden_layer))
  # Generate all combinations of hyperparameters
  hyperparam_combinations <- expand.grid(param_grid)

  for (i in 1:nrow(hyperparam_combinations)) {
    # Extract hyperparameter values for the current combination
    hyperparams <- hyperparam_combinations[i, ]

    num_hidden_layers <- hyperparams$num_hidden_layers
    #neurons_per_layer <- hyperparams$neurons_per_layer
    neurons_per_layer <- if(inherits(hyperparams$neurons_per_layer, "list")) unlist(hyperparams$neurons_per_layer) else hyperparams$neurons_per_layer
    learning_rate <- hyperparams$learning_rate

    # Ensure that the number of neurons_per_layer matches the number of num_hidden_layers
    if (length(neurons_per_layer) != num_hidden_layers) {
      warning(paste(msg, "Length of neurons_per_layer does not match num_hidden_layers. Skipping this combination."), call. = FALSE)
      next
    }

    # Create and compile the model with current hyperparameters
    model <- deep_learning_model_utility(X_train,
                                         y_train,
                                         num_hidden_layers,
                                         neurons_per_layer,
                                         learning_rate,
                                         epochs,
                                         batch_size)

    # Train the model with or without early stopping
    if (early_stop) {
      # Create a validation set from the training data
      indices <- sample.int(nrow(X_train), size = floor(validation_split * nrow(X_train)))
      X_val <- X_train[-indices, ]
      y_val <- y_train[-indices]

      # Define early stopping callback
      early_stopping <- keras::callback_early_stopping(monitor = "val_loss", patience = 3)

      # Train the model with early stopping
      history <- model |> keras::fit(
        x = X_train,
        y = y_train,
        epochs = epochs,
        batch_size = batch_size,
        verbose = 0,
        validation_data = list(X_val, y_val),
        callbacks = list(early_stopping)
      )
    } else {
      # Train the model without early stopping
      history <- model |> keras::fit(
        x = X_train,
        y = y_train,
        epochs = epochs,
        batch_size = batch_size,
        verbose = 0
      )
    }

    # Store the model and training history with a systematic name
    result_name <- paste0("Model_", i)  # You can adjust the naming scheme as needed
    results[[result_name]] <- list(
      model = model,
      history = history
    )
  }

  return(results)
}
