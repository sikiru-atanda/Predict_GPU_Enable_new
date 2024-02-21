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
#' @importFrom magrittr |>
#' @examples
grid_search_deeplearningOLD <- function(X_train,
                        y_train,
                        param_grid,
                        epochs,
                        batch_size,
                        validation_split = 0.2,
                        early_stop = TRUE) {
  results <- list()

  # Default values for hyperparameters
  default_values <- list(
    num_hidden_layers = 1,
    neurons_per_layer = ncol(X_train)/2,
    learning_rate = 0.001,
    epochs = 10,
    batch_size = 32
  )

  # Replace default values with user-provided values if available
  for (param_name in names(param_grid)) {
    default_values[[param_name]] <- param_grid[[param_name]]
  }

  if (length(param_grid[["num_hidden_layers"]]) != length(param_grid[["neurons_per_layer"]])) {
    stop("Length of num_hidden_layers must match length of neurons_per_layer.")
  }

  for (num_hidden_layers in default_values$num_hidden_layers) {
    for (neurons_per_layer in default_values$neurons_per_layer) {
      for (learning_rate in default_values$learning_rate) {
        for (epoch in default_values$epochs) {
          for (batch_size_val in default_values$batch_size) {
            # Create and compile the model with current hyperparameters
            model <- deep_learning_model_utility(X_train,
                                                 y_train,
                                                 num_hidden_layers,
                                                 neurons_per_layer,
                                                 learning_rate,
                                                 epoch,
                                                 batch_size_val)
            if(isFALSE(early_stop)){
            # Train the model
            history <- model |> keras::fit(
              x = X_train, y = y_train,
              epochs = epoch,
              batch_size = batch_size_val,
              verbose = 0
            )

            } else {
              if(is.null(validation_split)| is.na(validation_split)){
                stop("Provide validation_split to implement early_stop functionality")
              }
              # Create a validation set from the training data
              validation_split <- validation_split  # 20% of the data for validation
              indices <- sample.int(nrow(X_train), size = floor(validation_split * nrow(X_train)))
              X_val <- X_train[-indices, ]
              y_val <- y_train[-indices]

              # Define early stopping callback
              early_stopping <- keras::callback_early_stopping(monitor = "val_loss", patience = 3)

              # Train the model with early stopping
              history <- model |> keras::fit(
                x = X_train,
                y = y_train,
                epochs = epoch,
                batch_size = batch_size_val,
                verbose = 0,
                validation_data = list(X_val, y_val),
                callbacks = list(early_stopping)
              )
            }

            # Save the model and training history
            results[[paste(num_hidden_layers,
                           neurons_per_layer,
                           learning_rate,
                           epoch,
                           batch_size_val, collapse = "_")]] <- list(
                             model = model,
                             history = history
                           )
          }
        }
      }
    }
  }

  return(results)
}
