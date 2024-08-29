

# Function to validate the number of hidden layers against the neurons_per_layer list
validate_layers <- function(num_hidden_layers, neurons_per_layer) {
  msg <- "\n==================================================\n"
  if(is.null(num_hidden_layers) || is.null(neurons_per_layer)){
  if (length(num_hidden_layers) != length(neurons_per_layer)) {
    stop(paste(msg, "The length of num_hidden_layers must match the length of neurons_per_layer"), call. = FALSE)
  }
  }

if(length(num_hidden_layers)>1){
  for (i in seq_along(num_hidden_layers)) {
    if (num_hidden_layers[i] != length(neurons_per_layer[[i]])) {
      stop(paste(msg, paste("Mismatch in hidden layers and neurons per layer at", num_hidden_layers[i])), call. = FALSE)
    }
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
                                      learning_rate,
                                      dropout_rate = 0.5,
                                      kernel_size = 3,
                                      dense_layers_cnn,
                                      l2_regularizer_dp,
                                      n_blocks = 2,
                                      n_neurons_per_block,
                                      deep_learning_model = "mlp_with_attention",
                                      attention_on_final_layer = TRUE,
                                      attention_across_multiple_layers = FALSE,
                                      batch_normalization = TRUE,
                                      para_tunning = TRUE,
                                      early_stop = TRUE) {
  results <- list()
  np <- reticulate::import("numpy")
  msg <- "\n==================================================\n"

  if(deep_learning_model == "ResNet"){
    if(!is.null(n_blocks) && !is.null(n_neurons_per_block)){
    if(n_blocks==length(n_neurons_per_block)){
      stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
    }
    } else {
      stop(paste(msg,"n_neurons_per_block and n_block can't be NULL."), call. = FALSE)
}
    neurons_per_layer <- n_neurons_per_block
    num_hidden_layers <- n_blocks
  }

  X_train <- np$array(as.matrix(X_train))
  y_train <- np$array(y_train)
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

  for (i in 1:nrow(hyperparam_combinations)) {
    # Extract hyperparameter values for the current combination
    hyperparams <- hyperparam_combinations[i, ]

    # Initialize variables only if they exist in the hyperparams
    if ("num_hidden_layers" %in% names(hyperparams)) {
      num_hidden_layers <- as.integer(hyperparams$num_hidden_layers)
    }

    if ("neurons_per_layer" %in% names(hyperparams)) {
      if (inherits(hyperparams$neurons_per_layer, "list")) {
        neurons_per_layer <- unlist(hyperparams$neurons_per_layer)
      } else {
        neurons_per_layer <- hyperparams$neurons_per_layer
      }
    }

    if ("learning_rate" %in% names(hyperparams)) {
      learning_rate <- hyperparams$learning_rate
    }

    if ("dropout_rate" %in% names(hyperparams)) {
      dropout_rate <- hyperparams$dropout_rate
    }

    if ("epochs" %in% names(hyperparams)) {
      epochs <- hyperparams$epochs
    }

    if ("l2_regularizer_dp" %in% names(hyperparams)) {
      l2_regularizer_dp <- hyperparams$l2_regularizer_dp
    }

    if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers)
    if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.integer(neurons_per_layer)
    if(!is.null(batch_size))   batch_size <- as.integer(batch_size)
    if(!is.null(epochs)) epochs <- as.integer(epochs)
    if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate)
    if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate)

       # Ensure that the number of neurons_per_layer matches the number of num_hidden_layers
    if (length(neurons_per_layer) != num_hidden_layers) {
      warning(paste(msg, "Length of neurons_per_layer does not match num_hidden_layers. Skipping this combination."), call. = FALSE)
      next
    }

    # Create and compile the model with current hyperparameters
    model <- deep_learning_model_utilityy(X_train = X_train,
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
                                         n_neurons_per_block = n_neurons_per_block,
                                         deep_learning_model = deep_learning_model,
                                         attention_on_final_layer = attention_on_final_layer,
                                         attention_across_multiple_layers = attention_across_multiple_layers,
                                         batch_normalization = batch_normalization,
                                         para_tunning = FALSE)

       # Train the model with or without early stopping
    if (early_stop) {
      # Create a validation set from the training data
      set.seed(123)
      indices <- sample.int(nrow(X_train), size = floor(validation_split * nrow(X_train)))
      X_val <- as.matrix(X_train[indices, ])
      y_val <- y_train[indices]
      X_trainn <- X_train[-indices, ]
      y_trainn <- y_train[-indices]

      # Define early stopping callback
      early_stopping <- keras::callback_early_stopping(monitor = "val_loss", patience = 3)


      # Train the model with early stopping
      history <- model$fit(
        x = X_trainn,
        y = y_trainn,
        epochs = epochs,
        batch_size = batch_size,
        validation_data = list(X_val, y_val),
        callbacks = list(early_stopping),
        verbose = 0
      )
    } else {
      # Train the model without early stopping

      history <- model$fit(
        x = X_trainn,
        y = y_trainn,
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

    rm(model)
  }

  return(results)
}
