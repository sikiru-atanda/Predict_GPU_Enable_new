#' Build and train a simple feed-forward MLP for genomic prediction
#'
#' A thin utility that constructs and fits a multi-layer-perceptron regressor
#' on the supplied training features and response, using the requested depth /
#' width, learning rate and batch size. Used by the deep-learning wrappers.
#'
#' @param X_train Numeric matrix of training features (rows = individuals,
#'   columns = features).
#' @param y_train Numeric vector of training responses aligned with `X_train`.
#' @param num_hidden_layers Integer; number of hidden layers in the MLP.
#' @param neurons_per_layer Integer; neurons per hidden layer.
#' @param learning_rate Numeric; optimiser learning rate.
#' @param epochs Integer; number of training epochs.
#' @param batch_size Integer; mini-batch size used during training.
#'
#' @return The trained model object.
#' @export
#' @examples
deep_learning_model_utility <- function(X_train,
                                       y_train,
                                       num_hidden_layers,
                                       neurons_per_layer,
                                       learning_rate,
                                       epochs,
                                       batch_size) {
  input_dim <- ncol(X_train)

  model <- keras::keras_model_sequential()

  model |>
    keras::layer_dense(units = neurons_per_layer[[1]],
                       activation = 'relu', input_shape = c(input_dim))
  # Adding the input layer
  if(length(neurons_per_layer)>1){

    # Adding the hidden layers
    for (i in 2:num_hidden_layers) {
      if (i <= length(neurons_per_layer)) {
        model |>
          keras::layer_dense(units = neurons_per_layer[[i]], activation = 'relu')
      } else {
        # Handle case where there are fewer elements in neurons_per_layer than num_hidden_layers
        # You might want to print a warning or take appropriate action
        # For example, you could choose a default number of neurons or stop adding more layers
      }
    }

  }


  # Determine the appropriate loss function, activation function, and metric based on the nature of the response variable
  if (length(unique(y_train)) == 2) {
    loss_function <- 'binary_crossentropy'
    output_activation <- 'sigmoid'  # For binary classification
    metric <- 'accuracy'            # Accuracy is suitable
  } else if (length(unique(y_train)) > 2 && length(unique(y_train)) < 10) {
    loss_function <- 'categorical_crossentropy'
    output_activation <- 'softmax'  # For multi-class classification
    metric <- 'accuracy'            # Accuracy might be suitable but consider other metrics
  } else {
    loss_function <- 'mean_squared_error'
    output_activation <- 'linear'    # For regression
    metric <- 'mean_absolute_error' # Use MAE for regression
  }

  # Adding the output layer with appropriate activation function
  model |>
    keras::layer_dense(units = 1, activation = output_activation)

  # Compile the model
  model |> keras::compile(
    loss = loss_function,
    optimizer = keras::optimizer_adam(learning_rate = learning_rate),
    metrics = metric
  )

  return(model)
}
