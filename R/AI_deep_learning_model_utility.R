

# Function to create a residual block
residual_block <- function(input_tensor,
                           units,
                           l2_regularizer_dp = 0.001,
                           dropout_rate = 0.2,
                           batch_normalization = TRUE) {


  # Import the necessary Keras regularizer
  tf <- reticulate::import("tensorflow")
  keras <- reticulate::import("keras")

  if (!is.null(l2_regularizer_dp)) {
    l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
  }

  # Apply the first dense layer with specified number of units and L2 regularization
  # The ReLU activation function is used to introduce non-linearity

  if(isFALSE(batch_normalization)){
    if(!is.null(l2_regularizer_dp)){
      x <- keras::layer_dense(input_tensor, units = units,
                              kernel_regularizer = l2,
                              activation = "relu") # ReLU activation for intermediate layers
    } else{
      x <- keras::layer_dense(input_tensor, units = units,
                              activation = "relu") # ReLU activation for intermediate layers
    }
  } else {
    # Apply the first dense layer with L2 regularization
    if(!is.null(l2_regularizer_dp)){
      x <- keras::layer_dense(input_tensor, units = units,
                              kernel_regularizer = l2)
    } else {
      x <- keras::layer_dense(input_tensor, units = units)
    }

    # Apply batch normalization after the dense layer
    x <- keras::layer_batch_normalization(x)
    # Apply ReLU activation after batch normalization
    x <- keras::layer_activation(x, activation = "relu")
  }
  # Apply dropout to the layer's output to prevent overfitting
  if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)

  # Apply a second dense layer with the same number of units and L2 regularization
  # This layer does not have an activation function, as the ReLU will be applied after the skip connection
  if(!is.null(l2_regularizer_dp)){
    x <- keras::layer_dense(x, units = units,
                            kernel_regularizer = l2)
  } else {
    x <- keras::layer_dense(x, units = units)
  }

  # Apply batch normalization again before the skip connection
  if(isTRUE(batch_normalization)){
    x <- keras::layer_batch_normalization(x)

  }

  # Ensure the input tensor has the same number of units as the output tensor of the dense layers
  # If the dimensions do not match, apply a linear transformation to adjust the input tensor
  # Match the shape of input_tensor with x, if necessary
  if (input_tensor$shape[[2]] != units) {
    if(!is.null(l2_regularizer_dp)){
      input_tensor <- keras::layer_dense(input_tensor, units = units,
                                         kernel_regularizer = l2,
                                         activation = "linear")
    } else {
      input_tensor <- keras::layer_dense(input_tensor, units = units,
                                         activation = "linear")
    }
  }

  # Skip connection
  # Perform the skip connection by adding the original input tensor to the output tensor
  # This helps preserve the identity information and enables the network to learn residuals
  x <- keras::layer_add(list(x, input_tensor))

  # Apply ReLU activation to the combined tensor (after the skip connection)
  # This introduces non-linearity after the residual addition
  x <- keras::layer_activation(x, activation = "relu")

  # Return the final tensor output of the residual block
  return(x)
}

#####
# Function to build a ResNet model
build_resnet_model <- function(input_shape,
                               n_blocks,
                               n_neurons_per_block,
                               l2_regularizer_dp,
                               dropout_rate,
                               loss_function,
                               optimizer,
                               metric,
                               output_activation) {


  if(is.null(n_blocks)) n_blocks <- 2

  # Define the input layer with the specified input shape
  inputs <- keras::layer_input(shape = input_shape)
  x <- inputs

  # Build the ResNet by stacking the specified number of residual blocks
  for (i in 1:n_blocks) {
    # Add a residual block with the specified number of neurons, L2 regularization, and dropout rate
    x <- residual_block(x, units = n_neurons_per_block[i],
                        l2_regularizer_dp = l2_regularizer_dp,
                        dropout_rate = dropout_rate)

  }

  # # Add the output layer with a single unit (neuron) and the dynamically determined activation function
  # # The output layer's activation depends on whether the task is regression, binary classification, or multi-class classification
  x <- keras::layer_dense(x, units = 1, activation = output_activation)  # Output layer for regression

  # # Create the Keras model object, defining inputs and outputs
  model <- keras::keras_model(inputs = inputs, outputs = x)

  # Compile the model with the determined loss function, optimizer, and evaluation metrics
  model$compile(
    loss = loss_function,
    optimizer = optimizer,
    metrics = list(metric)
  )

  return(model)
 #return(x)
}

# Function to build a CNN model
build_cnn_model <- function(input_shape,
                            num_hidden_layers = 3,
                            neurons_per_layer,
                            l2_regularizer_dp,
                            dropout_rate,
                            dense_layers_cnn = c(256, 128, 64),
                            kernel_size = 3,
                            batch_normalization = TRUE,
                            optimizer,
                            metric,
                            loss_function,
                            output_activation,
                            validation_split = 0.2) {

  # Import the necessary Keras regularizer
  tf <- reticulate::import("tensorflow")
  keras <- reticulate::import("keras")

  if (!is.null(l2_regularizer_dp)) {
    l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
  }

  if (!is.null(l2_regularizer_dp)) {
    l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
  }

  inputs <- keras::layer_input(shape = input_shape)
  x <- inputs

  # Loop through each layer to add convolutional layers
  for (i in 1:num_hidden_layers) {
    # Add a 1D convolutional layer
    # - filters = n_neurons_per_layer[i]: Number of output filters (neurons) in the convolution.
    # - kernel_size = 3: The size of the convolution window.
    # - padding = "same": Ensures the output size matches the input size.
    # - activation = NULL: Activation is applied after batch normalization.
    # - kernel_regularizer = l2: Applies L2 regularization if l2_reg is specified.
    if(isTRUE(batch_normalization)){
      # Apply batch normalization to stabilize and accelerate training
      if(!is.null(l2_regularizer_dp)){
        x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
                                  padding = "same", activation = NULL,
                                  kernel_regularizer = l2)
      } else {
        x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
                                  padding = "same", activation = NULL)
      }
    } else{
      if(!is.null(l2_regularizer_dp)){
        x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
                                  padding = "same", activation = "relu",
                                  kernel_regularizer = l2)
      } else {
        x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
                                  padding = "same", activation = "relu")
      }
    }
    # Apply dropout to prevent overfitting
    if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)
  }

  if(length(dense_layers_cnn)>1) {
    # Sort the vector in descending order
    dense_layers_sorted <-  sort(dense_layers_cnn, decreasing = TRUE)

  }

  # Flatten the output from the convolutional layers to prepare it for the dense layers
  x <- keras::layer_flatten(x)

  # Dynamically add multiple dense layers
  for (units in dense_layers_cnn) {
    if(isTRUE(batch_normalization)){
      x <- keras::layer_dense(x, units = units, activation = NULL)
      x <- keras::layer_batch_normalization(x)
    }
    x <- keras::layer_activation(x, activation = "relu")
    # Apply dropout to further prevent overfitting
    if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)
  }

  # Final output layer
  x <- keras::layer_dense(x, units = 1, activation = output_activation)

  model <- keras::keras_model(inputs = inputs, outputs = x)

  model$compile(
    loss = loss_function,
    optimizer = optimizer,
    metrics = list(metric)
  )

  return(model)
  #return(x)
}

# Helper function to build MLP layers
build_mlp_layers <- function(input,
                             neurons_per_layer,
                             l2_regularizer_dp,
                             dropout_rate,
                             batch_normalization) {
  output <- input
  tf <- reticulate::import("tensorflow")
  keras <- reticulate::import("keras")
  if (!is.null(l2_regularizer_dp)) {
    l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
  }

  for (i in seq_along(neurons_per_layer)) {
    if (isTRUE(batch_normalization)) {
      output <- output |> keras::layer_dense(units = neurons_per_layer[i],
                                             kernel_regularizer = l2,
                                             activation = NULL) |>
        keras::layer_batch_normalization() |>
        keras::layer_activation('relu')
    } else {
      output <- output |> keras::layer_dense(units = neurons_per_layer[i],
                                             kernel_regularizer = l2,
                                             activation = 'relu')
    }

    if (!is.null(dropout_rate)) {
      output <- output |> keras::layer_dropout(rate = dropout_rate)
    }
  }
  return(output)
}


# Helper function to build attention layers across multiple layers
build_attention_layers <- function(input,
                                   neurons_per_layer,
                                   l2_regularizer_dp,
                                   dropout_rate,
                                   batch_normalization) {

  tf <- reticulate::import("tensorflow")
  keras <- reticulate::import("keras")
  if (!is.null(l2_regularizer_dp)) {
    l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
  }

  hidden_layers <- list()
  output <- input

  for (i in seq_along(neurons_per_layer)) {
    # Add a dense layer with optional L2 regularization
    if (!is.null(l2_regularizer_dp)) {
      output <- output |> keras::layer_dense(units = neurons_per_layer[i],
                                             kernel_regularizer = l2,
                                             activation = NULL)
    } else {
      output <- output |> keras::layer_dense(units = neurons_per_layer[i],
                                             activation = NULL)
    }
    # Apply batch normalization if enabled
    if (isTRUE(batch_normalization)) {
      output <- output |> keras::layer_batch_normalization() |>
        keras::layer_activation('relu')
    } else {
      output <- output |> keras::layer_activation('relu')
    }

    # Apply dropout if specified
    if (!is.null(dropout_rate)) {
      output <- output |> keras::layer_dropout(rate = dropout_rate)
    }

    hidden_layers[[i]] <- output
  }

  # Concatenate all hidden layers' outputs
  concatenated_output <- keras::layer_concatenate(hidden_layers)

  # Apply the attention mechanism
  attention_probs <- concatenated_output |> keras::layer_dense(units = sum(neurons_per_layer), activation = 'softmax')
  attention_output <- keras::layer_multiply(list(concatenated_output, attention_probs))

  return(attention_output)
}

# Helper function to get the optimizer
get_optimizer <- function(optimizer_name, learning_rate) {
  switch(optimizer_name,
         adam = keras::optimizer_adam(learning_rate = learning_rate),
         adamax = keras::optimizer_adamax(learning_rate = learning_rate),
         sgd = keras::optimizer_sgd(learning_rate = learning_rate),
         rmsprop = keras::optimizer_rmsprop(learning_rate = learning_rate),
         adadelta = keras::optimizer_adadelta(learning_rate = learning_rate),
         nadam = keras::optimizer_nadam(learning_rate = learning_rate))
}



deep_learning_model_utilityy <- function(X_train,
                                         y_train,
                                         num_hidden_layers = 2,
                                         neurons_per_layer = NULL,
                                         learning_rate = 0.001,
                                         epochs = 32,
                                         batch_size = 5,
                                         n_blocks  = 2,
                                         n_neurons_per_block = NULL,
                                         dense_layers_cnn = c(256, 128, 64),
                                         kernel_size = 3,
                                         l2_regularizer_dp = 0.001,
                                         dropout_rate = 0.5,
                                         para_tunning = FALSE,
                                         output_optimizer = "adam",
                                         deep_learning_model = "mlp_with_attention",
                                         attention_on_final_layer = TRUE,
                                         attention_across_multiple_layers = FALSE,
                                         batch_normalization = TRUE,
                                         validation_split = 0.2) {
#browser()
  msg <- "\n==================================================\n"
  # Convert input parameters to appropriate types
  if(!is.null(batch_size))   batch_size <- as.integer(batch_size)
  if(!is.null(epochs)) epochs <- as.integer(epochs)
  if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate)

  if(isTRUE(attention_on_final_layer) && isTRUE(attention_across_multiple_layers)) attention_across_multiple_layers <- FALSE

  if(deep_learning_model == "ResNet"){
    if(!is.null(n_blocks) && !is.null(n_neurons_per_block)){
    if(n_blocks!=length(n_neurons_per_block)){
      stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
    }
    } else {
      stop(paste(msg, "n_block and n_neurons_per_block, can't be NULL."), call. = FALSE)
}
    neurons_per_layer <- n_neurons_per_block
    num_hidden_layers <- n_blocks
  }


  # Validate optimizer
  valid_optimizers <- c("adam", "adamax", "sgd", "rmsprop", "adadelta", "nadam")
  # Error handling for output optimizer
  if (!(output_optimizer %in% valid_optimizers)) {
    stop(paste(msg, "Invalid output optimizer. Choose from: ", paste(valid_optimizers, collapse = ", ")), call. = FALSE)
  }

  # Validate hidden layers and neurons per layer
  if (isTRUE(para_tunning)) {
    validate_layers(num_hidden_layers, neurons_per_layer)
  }

  if(inherits(neurons_per_layer, "list")) neurons_per_layer <- unlist(neurons_per_layer)
  if(inherits(num_hidden_layers, "list")) num_hidden_layers <- unlist(num_hidden_layers)
  if(length(num_hidden_layers)>1) stop(paste(msg, "num_hidden_layers should be vector of length 1"), call. = FALSE)
  if (num_hidden_layers!= length(neurons_per_layer)) {
    stop(paste(msg, paste("Mismatch in hidden layers and neurons per layer at", num_hidden_layers)), call. = FALSE)
  }

  # Determine the appropriate loss function, activation function, and metric based on the nature of the response variable
  if (length(unique(y_train))!= length(y_train) & length(unique(y_train)) == 2) {
    loss_function <- 'binary_crossentropy'
    output_activation <- 'sigmoid'  # For binary classification
    metric <- 'accuracy'            # Accuracy is suitable
  } else if (length(unique(y_train))!= length(y_train) && length(unique(y_train)) > 2 && length(unique(y_train)) < 10) {
    loss_function <- 'categorical_crossentropy'
    output_activation <- 'softmax'  # For multi-class classification
    metric <- 'accuracy'            # Accuracy might be suitable but consider other metrics
  } else {
    loss_function <- 'mean_squared_error'
    output_activation <- 'linear'    # For regression
    metric <- 'mean_absolute_error' # Use MAE for regression
  }

  ### Determine the optimizer
  optimizer <- get_optimizer(output_optimizer, learning_rate)

  input_shape <- ncol(X_train)

  if(deep_learning_model == "cnn") input_shape <- c(ncol(X_train), 1)

  # Convert input data to numpy arrays
  np <- reticulate::import("numpy")
  X_train <- np$array(as.matrix(X_train), dtype = "float32")
  y_train <- np$array(y_train, dtype = "float32")

  # Define the input layer
  input <- keras::layer_input(shape = c(ncol(X_train)))

  # Build the model depending on the method
  if (deep_learning_model == "mlp_with_attention" || deep_learning_model == "mlp" && !attention_across_multiple_layers) {
    output <- build_mlp_layers(input, neurons_per_layer, l2_regularizer_dp, dropout_rate, batch_normalization)
  } else if (deep_learning_model == "mlp_with_attention" && attention_across_multiple_layers) {
    output <- build_attention_layers(input, neurons_per_layer, l2_regularizer_dp, dropout_rate, batch_normalization)
  }


  # Add final layers based on attention configuration
  if (attention_on_final_layer && deep_learning_model == "mlp_with_attention") {
    attention_probs <- output |> keras::layer_dense(units = neurons_per_layer[num_hidden_layers],
                                                    activation = 'softmax')
    output <- keras::layer_multiply(list(output, attention_probs)) |>
      keras::layer_dense(units = 1, activation = output_activation)

    # Compile the model
    model_dp <- keras::keras_model(inputs = input, outputs = output)
    model_dp$compile(loss = loss_function,
                     optimizer = optimizer,
                     metrics = list(metric))


  } else {
    if(deep_learning_model == "mlp" || (deep_learning_model == "mlp_with_attention" && isTRUE(attention_across_multiple_layers))) {
      # Standard MLP or MLP with attention across multiple layers
      output <- output |> keras::layer_dense(units = 1, activation = output_activation)

    # Compile the model
    model_dp <- keras::keras_model(inputs = input,
                                   outputs = output)
    model_dp$compile(loss = loss_function,
                     optimizer = optimizer,
                     metrics = list(metric))

    }
  }

  if(deep_learning_model == "ResNet"){

    model_dp <- build_resnet_model(input_shape = input_shape,
                                   n_blocks = num_hidden_layers,
                                   n_neurons_per_block = neurons_per_layer,
                                   l2_regularizer_dp = l2_regularizer_dp,
                                   dropout_rate = dropout_rate,
                                   loss_function = loss_function,
                                   optimizer = optimizer,
                                   metric = metric,
                                   output_activation = output_activation)
  }

  if(deep_learning_model =="cnn"){
    model_dp <- build_cnn_model(input_shape = input_shape,
                                num_hidden_layers = num_hidden_layers,
                                neurons_per_layer = neurons_per_layer,
                                l2_regularizer_dp = l2_regularizer_dp,
                                dropout_rate = dropout_rate,
                                dense_layers_cnn = dense_layers_cnn,
                                kernel_size = kernel_size,
                                batch_normalization = batch_normalization,
                                optimizer = optimizer,
                                metric = metric,
                                loss_function = loss_function,
                                output_activation = output_activation,
                                validation_split = validation_split)
  }
  # Compile the model
  #model_dp <- keras::keras_model(inputs = input, outputs = output)
  #optimizer <- get_optimizer(output_optimizer, learning_rate)

  # model_dp$compile(loss = loss_function,
  #               optimizer = optimizer,
  #               metrics = list(metric))

  # Callbacks for early stopping and visual feedback
  callbacks <- list(
    keras::callback_early_stopping(monitor = "val_loss", mode = 'min', patience = 50),
    keras::callback_lambda(on_epoch_end = function(epoch, logs) {
      if (epoch %% 20 == 0) cat("\n")
      cat(".")
    })
  )

  # Fit the model
  model_fit <- model_dp$fit(x = X_train,
                           y = y_train,
                           epochs = epochs,
                           batch_size = batch_size,
                           validation_split = validation_split,
                           verbose = 0, callbacks = callbacks)

  return(model_dp)
}



