#' Title
#'
#' @param X_train
#' @param y_train
#' @param num_hidden_layers
#' @param neurons_per_layer
#' @param learning_rate
#' @param epochs
#' @param batch_size
#' @param output_optimizer
#'
#' @return
#' @export
#'
#' @examples
deep_learning_model_utility <- function(X_train,
                                        y_train,
                                        num_hidden_layers,
                                        neurons_per_layer,
                                        learning_rate,
                                        epochs,
                                        batch_size,
                                        dropout_rate,
                                        output_optimizer = "adam") {

  msg <- "\n==================================================\n"
  # Error handling for output optimizer
  valid_optimizers <- c("adam", "adamax", "sgd", "rmsprop", "adadelta", "nadam")
  if (!(output_optimizer %in% valid_optimizers)) {
    stop(paste(msg, "Invalid output optimizer. Choose from: ", paste(valid_optimizers, collapse = ", ")), call. = FALSE)
  }

  validate_layers(num_hidden_layers, neurons_per_layer)
  input_dim <- ncol(X_train)

  np <- reticulate::import("numpy")
  X_train <- np$array(as.matrix(X_train))
  y_train <- np$array(y_train)
  # Determines the input dimensionality based on the number of columns in X_train

  # Initializes a sequential Keras model
  model <- keras::keras_model_sequential()

  # Adds the input layer to the model with the specified number of neurons and activation function
  model |> keras::layer_dense(units = neurons_per_layer[[1]],
                              activation = 'relu', input_shape = c(input_dim))

  # Adding the hidden layers
  for (i in 2:num_hidden_layers) {
    if (i <= length(neurons_per_layer)) {
      # Adds each hidden layer to the model with the specified number of neurons and activation function
      model |> keras::layer_dense(units = neurons_per_layer[[i]], activation = 'relu')
      # Adds dropout to the layer if dropout_rate is provided
      if (!is.null(dropout_rate)) {
        model |> keras::layer_dropout(rate = dropout_rate)
      }
    } else {
      # Handle case where there are fewer elements in neurons_per_layer than num_hidden_layers
      warning(paste(msg, "Fewer elements in neurons_per_layer than num_hidden_layers."), call. = FALSE)
      break
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
  model |> keras::layer_dense(units = 1, activation = output_activation)

  # Initialize output optimizer based on user input
  switch(output_optimizer,
         adam = output_optimizer <- keras::optimizer_adam(learning_rate = learning_rate),
         adamax = output_optimizer <- keras::optimizer_adamax(learning_rate = learning_rate),
         sgd = output_optimizer <- keras::optimizer_sgd(learning_rate = learning_rate),
         rmsprop = output_optimizer <- keras::optimizer_rmsprop(learning_rate = learning_rate),
         adadelta = output_optimizer <- keras::optimizer_adadelta(learning_rate = learning_rate),
         nadam = output_optimizer <- keras::optimizer_nadam(learning_rate = learning_rate)
  )

  # Compile the model
  model |> keras::compile(
    loss = loss_function,
    #optimizer = keras::optimizer_adam(learning_rate = learning_rate),
    optimizer = output_optimizer,
    metrics = metric
  )

  return(model)
}
