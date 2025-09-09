# keep your existing get_dl_module()

# torch_predict <- function(model, X, device = NULL) {
#   mod <- get_dl_module()
#   mod$predict(model, X, device = device)
# }

# torch_fit_model <- function(X, y, ...) {
#   mod <- get_dl_module()
#   mod$fit_model(X, y, ...)
# }
# torch_fit_model <- function(X, y, ...) {
#   mod  <- get_dl_module()
#   dots <- list(...)
#
#   # --- Legacy → current arg renames ---
#   # model type
#   if (!is.null(dots$deep_learning_model) && is.null(dots$model_type)) {
#     mt <- tolower(as.character(dots$deep_learning_model))
#     dots$model_type <- switch(
#       mt,
#       "ft" = "ft_transformer",
#       "dcn" = "dcnv2",
#       mt
#     )
#     dots$deep_learning_model <- NULL
#   }
#
#   # attention / batch-norm / regularization / dropout
#   if (!is.null(dots$attention_on_final_layer) && is.null(dots$final_attention)) {
#     dots$final_attention <- isTRUE(dots$attention_on_final_layer)
#     dots$attention_on_final_layer <- NULL
#   }
#   if (!is.null(dots$batch_normalization) && is.null(dots$batch_norm)) {
#     dots$batch_norm <- isTRUE(dots$batch_normalization)
#     dots$batch_normalization <- NULL
#   }
#   if (!is.null(dots$l2_regularizer_dp) && is.null(dots$l2_weight_decay)) {
#     dots$l2_weight_decay <- as.numeric(dots$l2_regularizer_dp)
#     dots$l2_regularizer_dp <- NULL
#   }
#   if (!is.null(dots$dropout_rate) && is.null(dots$dropout)) {
#     dots$dropout <- as.numeric(dots$dropout_rate)
#     dots$dropout_rate <- NULL
#   }
#   if (!is.null(dots$compile) && is.null(dots$compile_model)) {
#     dots$compile_model <- isTRUE(dots$compile)
#     dots$compile <- NULL
#   }
#
#   # CNN pool aliases (in case old names slipped through)
#   if (!is.null(dots$pool_kernel) && is.null(dots$cnn_pool_kernel)) {
#     dots$cnn_pool_kernel <- as.integer(dots$pool_kernel); dots$pool_kernel <- NULL
#   }
#   if (!is.null(dots$pool_stride) && is.null(dots$cnn_pool_stride)) {
#     dots$cnn_pool_stride <- as.integer(dots$pool_stride); dots$pool_stride <- NULL
#   }
#   if (!is.null(dots$pool_padding) && is.null(dots$cnn_pool_padding)) {
#     dots$cnn_pool_padding <- as.integer(dots$pool_padding); dots$pool_padding <- NULL
#   }
#
#   # Call Python
#   args <- c(list(X = X, y = y), dots)
#   do.call(mod$fit_model, args)
# }


torch_fit_model <- function(X, y, ...) {
  mod  <- get_dl_module()
  dots <- list(...)

  # ---- Catch duplicate kwargs early (e.g., dense_layers_cnn twice) ----
  dups <- names(dots)[duplicated(names(dots))]
  if (length(dups)) {
    # keep it loud & helpful
    stop(
      "Duplicate kwargs passed to torch_fit_model(): ",
      paste(unique(dups), collapse = ", "),
      "\nValues seen:\n",
      paste(sprintf("  - %s: %s", dups, vapply(dots[dups], function(v) paste(capture.output(str(v)), collapse=" "), "")),
            collapse = "\n"),
      call. = FALSE
    )
  }

  # --- Legacy → current arg renames ---
  if (!is.null(dots$deep_learning_model) && is.null(dots$model_type)) {
    mt <- tolower(as.character(dots$deep_learning_model))
    dots$model_type <- switch(mt, "ft" = "ft_transformer", "dcn" = "dcnv2", mt)
    dots$deep_learning_model <- NULL
  }
  if (!is.null(dots$attention_on_final_layer) && is.null(dots$final_attention)) {
    dots$final_attention <- isTRUE(dots$attention_on_final_layer)
    dots$attention_on_final_layer <- NULL
  }
  if (!is.null(dots$batch_normalization) && is.null(dots$batch_norm)) {
    dots$batch_norm <- isTRUE(dots$batch_normalization)
    dots$batch_normalization <- NULL
  }
  if (!is.null(dots$l2_regularizer_dp) && is.null(dots$l2_weight_decay)) {
    dots$l2_weight_decay <- as.numeric(dots$l2_regularizer_dp)
    dots$l2_regularizer_dp <- NULL
  }
  if (!is.null(dots$dropout_rate) && is.null(dots$dropout)) {
    dots$dropout <- as.numeric(dots$dropout_rate)
    dots$dropout_rate <- NULL
  }
  if (!is.null(dots$compile) && is.null(dots$compile_model)) {
    dots$compile_model <- isTRUE(dots$compile)
    dots$compile <- NULL
  }

  # CNN pool aliases (old → new)
  if (!is.null(dots$pool_kernel) && is.null(dots$cnn_pool_kernel)) {
    dots$cnn_pool_kernel <- as.integer(dots$pool_kernel);  dots$pool_kernel <- NULL
  }
  if (!is.null(dots$pool_stride) && is.null(dots$cnn_pool_stride)) {
    dots$cnn_pool_stride <- as.integer(dots$pool_stride);  dots$pool_stride <- NULL
  }
  if (!is.null(dots$pool_padding) && is.null(dots$cnn_pool_padding)) {
    dots$cnn_pool_padding <- as.integer(dots$pool_padding); dots$pool_padding <- NULL
  }

  # OPTIONAL: auto-dedupe (keep last) if you’d rather not hard-stop:
  # dots <- dots[!duplicated(names(dots), fromLast = TRUE)]

  args <- c(list(X = X, y = y), dots)
  do.call(mod$fit_model, args)
}

# type = "raw" (logits / regression output), "prob", or "class"
# "auto" -> probs for classification, raw for regression.
torch_predict <- function(model, X, device = NULL,
                          type = c("auto", "raw", "prob", "class"),
                          threshold = 0.5) {
  type <- match.arg(type)
  mod <- get_dl_module()

  # Get model outputs from Python (logits for classification, mean for hetero reg)
  logits <- mod$predict(model, X, device = device)
  Z <- as.matrix(logits)
  n <- nrow(Z); k <- ncol(Z)

  # Try to read the model's declared task; fall back to shape
  task <- tryCatch(tolower(as.character(model$task)), error = function(e) NULL)
  if (is.null(task)) task <- if (k > 1) "multiclass" else "regression"

  if (type == "auto") {
    type <- if (task %in% c("binary", "multiclass")) "prob" else "raw"
  }
  if (type == "raw") {
    return(if (k == 1) as.numeric(Z[, 1]) else Z)
  }

  if (task == "binary") {
    # logits -> probability -> class
    v <- if (k == 1) as.numeric(Z[, 1]) else as.numeric(Z[, 1])
    p <- 1 / (1 + exp(-v))  # sigmoid
    if (type == "prob") return(p)
    return(as.integer(p >= threshold))  # 0/1
  } else if (task == "multiclass") {
    # stable softmax row-wise
    Zc <- sweep(Z, 1, apply(Z, 1, max), FUN = "-")
    E  <- exp(Zc)
    P  <- sweep(E, 1, rowSums(E), FUN = "/")
    if (type == "prob") return(P)
    cls <- max.col(P, ties.method = "first")
    return(as.integer(cls - 1L))  # 0-based to match Python remap
  } else {
    # regression (or GP-DKL mean, or heteroscedastic mean already handled in Python)
    return(if (k == 1) as.numeric(Z[, 1]) else Z)
  }
}

# Convenience: return mean & (optional) variance.
# - For GP-DKL: returns Bayesian mean/variance.
# - For heteroscedastic reg: returns mean/variance from the model’s second channel.
# - For classification / plain models: variance may be NULL.
torch_predict_uncertainty <- function(model, X, device = NULL, prefer_bayesian = TRUE) {
  mod <- get_dl_module()
  out <- mod$predict_with_uncertainty(model, X, device = device, prefer_bayesian = prefer_bayesian)
  # reticulate already converts the Python dict -> R list (mean, variance, kind)
  out
}

# Small helper aliases if you like:
torch_predict_proba <- function(model, X, device = NULL) {
  torch_predict(model, X, device = device, type = "prob")
}
torch_predict_class <- function(model, X, device = NULL, threshold = 0.5) {
  torch_predict(model, X, device = device, type = "class", threshold = threshold)
}


# New PyTorch-backed utility (replaces the old Keras version)
# deep_learning_model_utilityy <- function(
#     X_train, y_train,
#     # Common MLP/CNN/ResNet knobs
#     num_hidden_layers = NULL,
#     neurons_per_layer = NULL,
#     learning_rate = 1e-3,
#     epochs = 32,
#     batch_size = 64,
#     l2_regularizer_dp = 1e-3,
#     dropout_rate = 0.5,
#     validation_split = 0.2,
#     n_blocks = 2,
#     n_neurons_per_block = NULL,
#     deep_learning_model = "mlp_with_attention",   # "mlp","mlp_with_attention","resnet","cnn","ft_transformer","saint","tabnet","node","deepfm","dcnv2","nam","moe","gp_dkl"
#     attention_on_final_layer = TRUE,
#     attention_across_multiple_layers = FALSE,
#     batch_normalization = TRUE,
#     dense_layers_cnn = c(256, 128, 64),
#     kernel_size = 3,
#     # CNN max-pool controls
#     cnn_use_max_pool = FALSE,
#     cnn_pool_kernel  = 2,
#     cnn_pool_stride  = 2,
#     cnn_pool_padding = 0,
#     separable=TRUE,
#     dilations=c(1,2,4),
#     use_se=TRUE,
#     norm_type="group",
#     use_max_pool=TRUE,
#     pool_type="conv",
#     use_global_pool=TRUE,
#     # Compile / reproducibility / device
#     compile_model = TRUE,
#     deterministic = FALSE,
#     random_seed = NULL,
#     device = NULL,
#     # AMP / grad-clip
#     use_amp = TRUE,
#     max_grad_norm = 1.0,
#     # Heteroscedastic regression (regression tasks only)
#     heteroscedastic = FALSE,
#     # Optimizer & class weights
#     optimizer_name = "adam",
#     auto_class_weights = FALSE,
#
#     # FT-Transformer
#     ft_d_model = 192, ft_heads = 8, ft_layers = 4, ft_ff_mult = 4,
#     ft_dropout = 0.1, ft_token_dropout = 0.1, ft_use_cls = TRUE,
#     ft_scalar_tokenizer = TRUE,   # keep default behavior in Python
#
#     # SAINT
#     saint_d_model = 128, saint_heads = 8, saint_layers = 4, saint_ff_mult = 4,
#     saint_dropout = 0.1, saint_token_dropout = 0.1, saint_use_cls = TRUE,
#
#     # Grouping controls for FT/SAINT and NAM/MoE
#     use_grouping = TRUE, group_trigger = 2048, group_method = "auto",
#     init_group_size = 64, max_tokens = 1024, kmeans_batch = 4096, kmeans_iter = 100,
#
#     # TabNet
#     tabnet_steps = 5, tabnet_feature_dim = 64, tabnet_output_dim = 64,
#     tabnet_gamma = 1.5, tabnet_lambda_sparse = 1e-4,
#
#     # NODE
#     node_trees = 8, node_depth = 3,
#
#     # DeepFM
#     deepfm_k = 16, deepfm_hidden = c(128, 64),
#
#     # DCN (v2-lite)
#     dcn_layers = 3, dcn_hidden = c(256, 128),
#
#     # NAM
#     nam_hidden = c(32, 16), nam_activation = "relu",
#     nam_add_linear = TRUE, nam_l1 = 1e-4,
#
#     # MoE
#     moe_n_experts = 4, moe_expert_hidden = c(128, 64),
#     moe_gate_hidden = 128, moe_temperature = 1.0,
#     moe_sparse_topk = NULL, moe_entropy_reg = 0.0,
#
#     # GP-DKL (+ RFF fallback)
#     gp_use_variational = TRUE, gp_num_inducing = 512, gp_feature_dim = 64,
#     gp_kernel = "rbf", gp_ard = TRUE, gp_lr_mult = 0.5,
#     rff_features = 1024, rff_lengthscale = 1.0, rff_deep_hidden = c(128)
# ) {
#   # Normalize model type string (add nam, moe, gp_dkl)
#   model_type <- switch(tolower(deep_learning_model),
#                        "mlp"                = "mlp",
#                        "mlp_with_attention" = "mlp_with_attention",
#                        "resnet"             = "resnet",   "resnet50" = "resnet", "resnet34" = "resnet",
#                        "cnn"                = "cnn",
#                        "ft_transformer"     = "ft_transformer", "ft" = "ft_transformer",
#                        "saint"              = "saint",
#                        "tabnet"             = "tabnet",
#                        "node"               = "node",
#                        "deepfm"             = "deepfm",
#                        "dcn"                = "dcnv2",    "dcnv2" = "dcnv2",
#                        "nam"                = "nam",
#                        "moe"                = "moe",
#                        "gp_dkl"             = "gp_dkl",
#                        stop("Unsupported deep_learning_model: ", deep_learning_model)
#   )
#
#   # ---- ResNet-specific mapping (blocks <-> neurons_per_layer) ----
#   # if (identical(model_type, "resnet")) {
#   #   if (!is.null(n_neurons_per_block)) neurons_per_layer <- n_neurons_per_block
#   #   if (is.null(num_hidden_layers) && !is.null(n_blocks)) num_hidden_layers <- as.integer(n_blocks)
#   #   if (!is.null(num_hidden_layers) && !is.null(neurons_per_layer) &&
#   #       length(neurons_per_layer) != as.integer(num_hidden_layers)) {
#   #     stop("Mismatch between n_blocks (num_hidden_layers) and n_neurons_per_block (neurons_per_layer).")
#   #   }
#   #   if (is.null(neurons_per_layer)) {
#   #     stop("Provide n_neurons_per_block (mapped to neurons_per_layer) for ResNet.")
#   #   }
#   # }
#
#   # ---- Provide default neurons_per_layer for models that need it ----
#   uses_hidden <- model_type %in% c("mlp", "mlp_with_attention", "resnet", "cnn")
#   if (is.null(neurons_per_layer) && uses_hidden) {
#     if (identical(model_type, "cnn")) neurons_per_layer <- c(64, 64, 64) else neurons_per_layer <- c(256, 128)
#   }
#
#   # For models that don't use hidden counts, avoid sending num_hidden_layers
#   send_num_hidden <- if (uses_hidden) num_hidden_layers else NULL
#   send_neurons    <- if (uses_hidden) neurons_per_layer else NULL
#
#   # ---- y: make sure it’s numeric-coded even if character ----
#   y_vec <-
#     if (is.character(y_train)) as.numeric(factor(y_train)) else
#       if (is.factor(y_train))    as.numeric(y_train) else
#         as.numeric(y_train)
#
#   # ---- Call into Python (torch_fit_model -> fit_model) ----
#   fit <- torch_fit_model(
#     X = as.matrix(X_train),
#     y = y_vec,
#
#     model_type = model_type,
#     num_hidden_layers = if (!is.null(send_num_hidden)) as.integer(send_num_hidden) else NULL,
#     neurons_per_layer = if (!is.null(send_neurons)) as.integer(send_neurons) else NULL,
#
#     learning_rate = as.numeric(learning_rate),
#     epochs        = as.integer(epochs),
#     batch_size    = as.integer(batch_size),
#     l2_weight_decay = as.numeric(l2_regularizer_dp),
#     dropout       = as.numeric(dropout_rate),
#     optimizer_name = as.character(optimizer_name),
#
#     final_attention = isTRUE(attention_on_final_layer),
#     attention_across_multiple_layers = isTRUE(attention_across_multiple_layers),
#     batch_norm     = isTRUE(batch_normalization),
#
#     validation_split = as.numeric(validation_split),
#     compile_model    = isTRUE(compile_model),
#     device           = if (is.null(device)) NULL else as.character(device),
#
#     # CNN core + max-pool
#     kernel_size    = as.integer(kernel_size),
#     dense_layers_cnn = as.integer(dense_layers_cnn),
#     cnn_use_max_pool = isTRUE(cnn_use_max_pool),
#     cnn_pool_kernel  = as.integer(cnn_pool_kernel),
#     cnn_pool_stride  = as.integer(cnn_pool_stride),
#     cnn_pool_padding = as.integer(cnn_pool_padding),
#     separable=isTRUE(separable),
#     dilations=as.integer(dilations),
#     use_se=isTRUE(use_se),
#     norm_type=norm_type,
#     use_max_pool=isTRUE(use_max_pool),
#     pool_type=pool_type,
#     use_global_pool=isTRUE(use_global_pool),
#
#     # Repro / AMP / grad-clip
#     deterministic = isTRUE(deterministic),
#     random_seed   = if (is.null(random_seed)) NULL else as.integer(random_seed),
#     use_amp       = isTRUE(use_amp),
#     max_grad_norm = if (is.null(max_grad_norm)) NULL else as.numeric(max_grad_norm),
#
#     # Heteroscedastic regression
#     heteroscedastic = isTRUE(heteroscedastic),
#
#     # FT-Transformer
#     ft_d_model = as.integer(ft_d_model),
#     ft_heads   = as.integer(ft_heads),
#     ft_layers  = as.integer(ft_layers),
#     ft_ff_mult = as.integer(ft_ff_mult),
#     ft_dropout = as.numeric(ft_dropout),
#     ft_token_dropout = as.numeric(ft_token_dropout),
#     ft_use_cls = isTRUE(ft_use_cls),
#     ft_scalar_tokenizer = isTRUE(ft_scalar_tokenizer),
#
#     # SAINT
#     saint_d_model = as.integer(saint_d_model),
#     saint_heads   = as.integer(saint_heads),
#     saint_layers  = as.integer(saint_layers),
#     saint_ff_mult = as.integer(saint_ff_mult),
#     saint_dropout = as.numeric(saint_dropout),
#     saint_token_dropout = as.numeric(saint_token_dropout),
#     saint_use_cls = isTRUE(saint_use_cls),
#
#     # Grouping controls (FT/SAINT + NAM/MoE)
#     use_grouping  = isTRUE(use_grouping),
#     group_trigger = as.integer(group_trigger),
#     group_method  = as.character(group_method),
#     init_group_size = as.integer(init_group_size),
#     max_tokens    = as.integer(max_tokens),
#     kmeans_batch  = as.integer(kmeans_batch),
#     kmeans_iter   = as.integer(kmeans_iter),
#
#     # TabNet
#     tabnet_steps        = as.integer(tabnet_steps),
#     tabnet_feature_dim  = as.integer(tabnet_feature_dim),
#     tabnet_output_dim   = as.integer(tabnet_output_dim),
#     tabnet_gamma        = as.numeric(tabnet_gamma),
#     tabnet_lambda_sparse= as.numeric(tabnet_lambda_sparse),
#
#     # NODE
#     node_trees = as.integer(node_trees),
#     node_depth = as.integer(node_depth),
#
#     # DeepFM
#     deepfm_k      = as.integer(deepfm_k),
#     deepfm_hidden = as.integer(deepfm_hidden),
#
#     # DCN
#     dcn_layers = as.integer(dcn_layers),
#     dcn_hidden = as.integer(dcn_hidden),
#
#     # NAM
#     nam_hidden      = as.integer(nam_hidden),
#     nam_activation  = as.character(nam_activation),
#     nam_add_linear  = isTRUE(nam_add_linear),
#     nam_l1          = as.numeric(nam_l1),
#
#     # MoE
#     moe_n_experts    = as.integer(moe_n_experts),
#     moe_expert_hidden= as.integer(moe_expert_hidden),
#     moe_gate_hidden  = as.integer(moe_gate_hidden),
#     moe_temperature  = as.numeric(moe_temperature),
#     moe_sparse_topk  = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#     moe_entropy_reg  = as.numeric(moe_entropy_reg),
#
#     # GP-DKL (+ RFF fallback)
#     gp_use_variational = isTRUE(gp_use_variational),
#     gp_num_inducing    = as.integer(gp_num_inducing),
#     gp_feature_dim     = as.integer(gp_feature_dim),
#     gp_kernel          = as.character(gp_kernel),
#     gp_ard             = isTRUE(gp_ard),
#     gp_lr_mult         = as.numeric(gp_lr_mult),
#     rff_features       = as.integer(rff_features),
#     rff_lengthscale    = as.numeric(rff_lengthscale),
#     rff_deep_hidden    = as.integer(rff_deep_hidden),
#
#     # Class weighting (classification only)
#     auto_class_weights = isTRUE(auto_class_weights)
#   )
#
#   # Python returns (model, history)
#   return(fit)
# }





#
#
# # Function to create a residual block
# residual_block <- function(input_tensor,
#                            units,
#                            l2_regularizer_dp = 0.001,
#                            dropout_rate = 0.2,
#                            batch_normalization = TRUE) {
#
#
#   # Import the necessary Keras regularizer
#   tf <- reticulate::import("tensorflow")
#   keras <- reticulate::import("keras")
#
#   if (!is.null(l2_regularizer_dp)) {
#     l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
#   }
#
#   # Apply the first dense layer with specified number of units and L2 regularization
#   # The ReLU activation function is used to introduce non-linearity
#
#   if(isFALSE(batch_normalization)){
#     if(!is.null(l2_regularizer_dp)){
#       x <- keras::layer_dense(input_tensor, units = units,
#                               kernel_regularizer = l2,
#                               activation = "relu") # ReLU activation for intermediate layers
#     } else{
#       x <- keras::layer_dense(input_tensor, units = units,
#                               activation = "relu") # ReLU activation for intermediate layers
#     }
#   } else {
#     # Apply the first dense layer with L2 regularization
#     if(!is.null(l2_regularizer_dp)){
#       x <- keras::layer_dense(input_tensor, units = units,
#                               kernel_regularizer = l2)
#     } else {
#       x <- keras::layer_dense(input_tensor, units = units)
#     }
#
#     # Apply batch normalization after the dense layer
#     x <- keras::layer_batch_normalization(x)
#     # Apply ReLU activation after batch normalization
#     x <- keras::layer_activation(x, activation = "relu")
#   }
#   # Apply dropout to the layer's output to prevent overfitting
#   if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)
#
#   # Apply a second dense layer with the same number of units and L2 regularization
#   # This layer does not have an activation function, as the ReLU will be applied after the skip connection
#   if(!is.null(l2_regularizer_dp)){
#     x <- keras::layer_dense(x, units = units,
#                             kernel_regularizer = l2)
#   } else {
#     x <- keras::layer_dense(x, units = units)
#   }
#
#   # Apply batch normalization again before the skip connection
#   if(isTRUE(batch_normalization)){
#     x <- keras::layer_batch_normalization(x)
#
#   }
#
#   # Ensure the input tensor has the same number of units as the output tensor of the dense layers
#   # If the dimensions do not match, apply a linear transformation to adjust the input tensor
#   # Match the shape of input_tensor with x, if necessary
#   if (input_tensor$shape[[2]] != units) {
#     if(!is.null(l2_regularizer_dp)){
#       input_tensor <- keras::layer_dense(input_tensor, units = units,
#                                          kernel_regularizer = l2,
#                                          activation = "linear")
#     } else {
#       input_tensor <- keras::layer_dense(input_tensor, units = units,
#                                          activation = "linear")
#     }
#   }
#
#   # Skip connection
#   # Perform the skip connection by adding the original input tensor to the output tensor
#   # This helps preserve the identity information and enables the network to learn residuals
#   x <- keras::layer_add(list(x, input_tensor))
#
#   # Apply ReLU activation to the combined tensor (after the skip connection)
#   # This introduces non-linearity after the residual addition
#   x <- keras::layer_activation(x, activation = "relu")
#
#   # Return the final tensor output of the residual block
#   return(x)
# }
#
# #####
# # Function to build a ResNet model
# build_resnet_model <- function(input_shape,
#                                n_blocks,
#                                n_neurons_per_block,
#                                l2_regularizer_dp,
#                                dropout_rate,
#                                loss_function,
#                                optimizer,
#                                metric,
#                                output_activation) {
#
#
#   if(is.null(n_blocks)) n_blocks <- 2
#
#   # Define the input layer with the specified input shape
#   inputs <- keras::layer_input(shape = input_shape)
#   x <- inputs
#
#   # Build the ResNet by stacking the specified number of residual blocks
#   for (i in 1:n_blocks) {
#     # Add a residual block with the specified number of neurons, L2 regularization, and dropout rate
#     x <- residual_block(x, units = n_neurons_per_block[i],
#                         l2_regularizer_dp = l2_regularizer_dp,
#                         dropout_rate = dropout_rate)
#
#   }
#
#   # # Add the output layer with a single unit (neuron) and the dynamically determined activation function
#   # # The output layer's activation depends on whether the task is regression, binary classification, or multi-class classification
#   x <- keras::layer_dense(x, units = 1, activation = output_activation)  # Output layer for regression
#
#   # # Create the Keras model object, defining inputs and outputs
#   model <- keras::keras_model(inputs = inputs, outputs = x)
#
#   # Compile the model with the determined loss function, optimizer, and evaluation metrics
#   model$compile(
#     loss = loss_function,
#     optimizer = optimizer,
#     metrics = list(metric)
#   )
#
#   return(model)
#  #return(x)
# }
#
# # Function to build a CNN model
# build_cnn_model <- function(input_shape,
#                             num_hidden_layers = 3,
#                             neurons_per_layer,
#                             l2_regularizer_dp,
#                             dropout_rate,
#                             dense_layers_cnn = c(256, 128, 64),
#                             kernel_size = 3,
#                             batch_normalization = TRUE,
#                             optimizer,
#                             metric,
#                             loss_function,
#                             output_activation,
#                             validation_split = 0.2) {
#
#   # Import the necessary Keras regularizer
#   tf <- reticulate::import("tensorflow")
#   keras <- reticulate::import("keras")
#
#   if (!is.null(l2_regularizer_dp)) {
#     l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
#   }
#
#   if (!is.null(l2_regularizer_dp)) {
#     l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
#   }
#
#   inputs <- keras::layer_input(shape = input_shape)
#   x <- inputs
#
#   # Loop through each layer to add convolutional layers
#   for (i in 1:num_hidden_layers) {
#     # Add a 1D convolutional layer
#     # - filters = n_neurons_per_layer[i]: Number of output filters (neurons) in the convolution.
#     # - kernel_size = 3: The size of the convolution window.
#     # - padding = "same": Ensures the output size matches the input size.
#     # - activation = NULL: Activation is applied after batch normalization.
#     # - kernel_regularizer = l2: Applies L2 regularization if l2_reg is specified.
#     if(isTRUE(batch_normalization)){
#       # Apply batch normalization to stabilize and accelerate training
#       if(!is.null(l2_regularizer_dp)){
#         x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
#                                   padding = "same", activation = NULL,
#                                   kernel_regularizer = l2)
#       } else {
#         x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
#                                   padding = "same", activation = NULL)
#       }
#     } else{
#       if(!is.null(l2_regularizer_dp)){
#         x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
#                                   padding = "same", activation = "relu",
#                                   kernel_regularizer = l2)
#       } else {
#         x <- keras::layer_conv_1d(x, filters = neurons_per_layer[i], kernel_size = kernel_size,
#                                   padding = "same", activation = "relu")
#       }
#     }
#     # Apply dropout to prevent overfitting
#     if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)
#   }
#
#   if(length(dense_layers_cnn)>1) {
#     # Sort the vector in descending order
#     dense_layers_sorted <-  sort(dense_layers_cnn, decreasing = TRUE)
#
#   }
#
#   # Flatten the output from the convolutional layers to prepare it for the dense layers
#   x <- keras::layer_flatten(x)
#
#   # Dynamically add multiple dense layers
#   for (units in dense_layers_cnn) {
#     if(isTRUE(batch_normalization)){
#       x <- keras::layer_dense(x, units = units, activation = NULL)
#       x <- keras::layer_batch_normalization(x)
#     }
#     x <- keras::layer_activation(x, activation = "relu")
#     # Apply dropout to further prevent overfitting
#     if(!is.null(dropout_rate)) x <- keras::layer_dropout(x, rate = dropout_rate)
#   }
#
#   # Final output layer
#   x <- keras::layer_dense(x, units = 1, activation = output_activation)
#
#   model <- keras::keras_model(inputs = inputs, outputs = x)
#
#   model$compile(
#     loss = loss_function,
#     optimizer = optimizer,
#     metrics = list(metric)
#   )
#
#   return(model)
#   #return(x)
# }
#
# # Helper function to build MLP layers
# build_mlp_layers <- function(input,
#                              neurons_per_layer,
#                              l2_regularizer_dp,
#                              dropout_rate,
#                              batch_normalization) {
#   output <- input
#   tf <- reticulate::import("tensorflow")
#   keras <- reticulate::import("keras")
#   if (!is.null(l2_regularizer_dp)) {
#     l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
#   }
#
#   for (i in seq_along(neurons_per_layer)) {
#     if (isTRUE(batch_normalization)) {
#       output <- output |> keras::layer_dense(units = neurons_per_layer[i],
#                                              kernel_regularizer = l2,
#                                              activation = NULL) |>
#         keras::layer_batch_normalization() |>
#         keras::layer_activation('relu')
#     } else {
#       output <- output |> keras::layer_dense(units = neurons_per_layer[i],
#                                              kernel_regularizer = l2,
#                                              activation = 'relu')
#     }
#
#     if (!is.null(dropout_rate)) {
#       output <- output |> keras::layer_dropout(rate = dropout_rate)
#     }
#   }
#   return(output)
# }
#
#
# # Helper function to build attention layers across multiple layers
# build_attention_layers <- function(input,
#                                    neurons_per_layer,
#                                    l2_regularizer_dp,
#                                    dropout_rate,
#                                    batch_normalization) {
#
#   tf <- reticulate::import("tensorflow")
#   keras <- reticulate::import("keras")
#   if (!is.null(l2_regularizer_dp)) {
#     l2 <- keras$regularizers$l2(as.numeric(l2_regularizer_dp))
#   }
#
#   hidden_layers <- list()
#   output <- input
#
#   for (i in seq_along(neurons_per_layer)) {
#     # Add a dense layer with optional L2 regularization
#     if (!is.null(l2_regularizer_dp)) {
#       output <- output |> keras::layer_dense(units = neurons_per_layer[i],
#                                              kernel_regularizer = l2,
#                                              activation = NULL)
#     } else {
#       output <- output |> keras::layer_dense(units = neurons_per_layer[i],
#                                              activation = NULL)
#     }
#     # Apply batch normalization if enabled
#     if (isTRUE(batch_normalization)) {
#       output <- output |> keras::layer_batch_normalization() |>
#         keras::layer_activation('relu')
#     } else {
#       output <- output |> keras::layer_activation('relu')
#     }
#
#     # Apply dropout if specified
#     if (!is.null(dropout_rate)) {
#       output <- output |> keras::layer_dropout(rate = dropout_rate)
#     }
#
#     hidden_layers[[i]] <- output
#   }
#
#   # Concatenate all hidden layers' outputs
#   concatenated_output <- keras::layer_concatenate(hidden_layers)
#
#   # Apply the attention mechanism
#   attention_probs <- concatenated_output |> keras::layer_dense(units = sum(neurons_per_layer), activation = 'softmax')
#   attention_output <- keras::layer_multiply(list(concatenated_output, attention_probs))
#
#   return(attention_output)
# }
#
# # Helper function to get the optimizer
# get_optimizer <- function(optimizer_name, learning_rate) {
#   switch(optimizer_name,
#          adam = keras::optimizer_adam(learning_rate = learning_rate),
#          adamax = keras::optimizer_adamax(learning_rate = learning_rate),
#          sgd = keras::optimizer_sgd(learning_rate = learning_rate),
#          rmsprop = keras::optimizer_rmsprop(learning_rate = learning_rate),
#          adadelta = keras::optimizer_adadelta(learning_rate = learning_rate),
#          nadam = keras::optimizer_nadam(learning_rate = learning_rate))
# }
#
#
#
# deep_learning_model_utilityy <- function(X_train,
#                                          y_train,
#                                          num_hidden_layers = 2,
#                                          neurons_per_layer = NULL,
#                                          learning_rate = 0.001,
#                                          epochs = 32,
#                                          batch_size = 5,
#                                          n_blocks  = 2,
#                                          n_neurons_per_block = NULL,
#                                          dense_layers_cnn = c(256, 128, 64),
#                                          kernel_size = 3,
#                                          l2_regularizer_dp = 0.001,
#                                          dropout_rate = 0.5,
#                                          para_tunning = FALSE,
#                                          output_optimizer = "adam",
#                                          deep_learning_model = "mlp_with_attention",
#                                          attention_on_final_layer = TRUE,
#                                          attention_across_multiple_layers = FALSE,
#                                          batch_normalization = TRUE,
#                                          validation_split = 0.2) {
# #browser()
#   msg <- "\n==================================================\n"
#   # Convert input parameters to appropriate types
#   if(!is.null(batch_size))   batch_size <- as.integer(batch_size)
#   if(!is.null(epochs)) epochs <- as.integer(epochs)
#   if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate)
#
#   if(isTRUE(attention_on_final_layer) && isTRUE(attention_across_multiple_layers)) attention_across_multiple_layers <- FALSE
#
#   if(deep_learning_model == "ResNet"){
#     if(!is.null(n_blocks) && !is.null(n_neurons_per_block)){
#     if(n_blocks!=length(n_neurons_per_block)){
#       stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
#     }
#     } else {
#       stop(paste(msg, "n_block and n_neurons_per_block, can't be NULL."), call. = FALSE)
# }
#     neurons_per_layer <- n_neurons_per_block
#     num_hidden_layers <- n_blocks
#   }
#
#
#   # Validate optimizer
#   valid_optimizers <- c("adam", "adamax", "sgd", "rmsprop", "adadelta", "nadam")
#   # Error handling for output optimizer
#   if (!(output_optimizer %in% valid_optimizers)) {
#     stop(paste(msg, "Invalid output optimizer. Choose from: ", paste(valid_optimizers, collapse = ", ")), call. = FALSE)
#   }
#
#   # Validate hidden layers and neurons per layer
#   if (isTRUE(para_tunning)) {
#     validate_layers(num_hidden_layers, neurons_per_layer)
#   }
#
#   if(inherits(neurons_per_layer, "list")) neurons_per_layer <- unlist(neurons_per_layer)
#   if(inherits(num_hidden_layers, "list")) num_hidden_layers <- unlist(num_hidden_layers)
#   if(length(num_hidden_layers)>1) stop(paste(msg, "num_hidden_layers should be vector of length 1"), call. = FALSE)
#   if (num_hidden_layers!= length(neurons_per_layer)) {
#     stop(paste(msg, paste("Mismatch in hidden layers and neurons per layer at", num_hidden_layers)), call. = FALSE)
#   }
#
#   # Determine the appropriate loss function, activation function, and metric based on the nature of the response variable
#   if (length(unique(y_train))!= length(y_train) & length(unique(y_train)) == 2) {
#     loss_function <- 'binary_crossentropy'
#     output_activation <- 'sigmoid'  # For binary classification
#     metric <- 'accuracy'            # Accuracy is suitable
#   } else if (length(unique(y_train))!= length(y_train) && length(unique(y_train)) > 2 && length(unique(y_train)) < 10) {
#     loss_function <- 'categorical_crossentropy'
#     output_activation <- 'softmax'  # For multi-class classification
#     metric <- 'accuracy'            # Accuracy might be suitable but consider other metrics
#   } else {
#     loss_function <- 'mean_squared_error'
#     output_activation <- 'linear'    # For regression
#     metric <- 'mean_absolute_error' # Use MAE for regression
#   }
#
#   ### Determine the optimizer
#   optimizer <- get_optimizer(output_optimizer, learning_rate)
#
#   input_shape <- ncol(X_train)
#
#   if(deep_learning_model == "cnn") input_shape <- c(ncol(X_train), 1)
#
#   # Convert input data to numpy arrays
#   np <- reticulate::import("numpy")
#   X_train <- np$array(as.matrix(X_train), dtype = "float32")
#   y_train <- np$array(y_train, dtype = "float32")
#
#   # Define the input layer
#   input <- keras::layer_input(shape = c(ncol(X_train)))
#
#   # Build the model depending on the method
#   if (deep_learning_model == "mlp_with_attention" || deep_learning_model == "mlp" && !attention_across_multiple_layers) {
#     output <- build_mlp_layers(input, neurons_per_layer, l2_regularizer_dp, dropout_rate, batch_normalization)
#   } else if (deep_learning_model == "mlp_with_attention" && attention_across_multiple_layers) {
#     output <- build_attention_layers(input, neurons_per_layer, l2_regularizer_dp, dropout_rate, batch_normalization)
#   }
#
#
#   # Add final layers based on attention configuration
#   if (attention_on_final_layer && deep_learning_model == "mlp_with_attention") {
#     attention_probs <- output |> keras::layer_dense(units = neurons_per_layer[num_hidden_layers],
#                                                     activation = 'softmax')
#     output <- keras::layer_multiply(list(output, attention_probs)) |>
#       keras::layer_dense(units = 1, activation = output_activation)
#
#     # Compile the model
#     model_dp <- keras::keras_model(inputs = input, outputs = output)
#     model_dp$compile(loss = loss_function,
#                      optimizer = optimizer,
#                      metrics = list(metric))
#
#
#   } else {
#     if(deep_learning_model == "mlp" || (deep_learning_model == "mlp_with_attention" && isTRUE(attention_across_multiple_layers))) {
#       # Standard MLP or MLP with attention across multiple layers
#       output <- output |> keras::layer_dense(units = 1, activation = output_activation)
#
#     # Compile the model
#     model_dp <- keras::keras_model(inputs = input,
#                                    outputs = output)
#     model_dp$compile(loss = loss_function,
#                      optimizer = optimizer,
#                      metrics = list(metric))
#
#     }
#   }
#
#   if(deep_learning_model == "ResNet"){
#
#     model_dp <- build_resnet_model(input_shape = input_shape,
#                                    n_blocks = num_hidden_layers,
#                                    n_neurons_per_block = neurons_per_layer,
#                                    l2_regularizer_dp = l2_regularizer_dp,
#                                    dropout_rate = dropout_rate,
#                                    loss_function = loss_function,
#                                    optimizer = optimizer,
#                                    metric = metric,
#                                    output_activation = output_activation)
#   }
#
#   if(deep_learning_model =="cnn"){
#     model_dp <- build_cnn_model(input_shape = input_shape,
#                                 num_hidden_layers = num_hidden_layers,
#                                 neurons_per_layer = neurons_per_layer,
#                                 l2_regularizer_dp = l2_regularizer_dp,
#                                 dropout_rate = dropout_rate,
#                                 dense_layers_cnn = dense_layers_cnn,
#                                 kernel_size = kernel_size,
#                                 batch_normalization = batch_normalization,
#                                 optimizer = optimizer,
#                                 metric = metric,
#                                 loss_function = loss_function,
#                                 output_activation = output_activation,
#                                 validation_split = validation_split)
#   }
#   # Compile the model
#   #model_dp <- keras::keras_model(inputs = input, outputs = output)
#   #optimizer <- get_optimizer(output_optimizer, learning_rate)
#
#   # model_dp$compile(loss = loss_function,
#   #               optimizer = optimizer,
#   #               metrics = list(metric))
#
#   # Callbacks for early stopping and visual feedback
#   callbacks <- list(
#     keras::callback_early_stopping(monitor = "val_loss", mode = 'min', patience = 50),
#     keras::callback_lambda(on_epoch_end = function(epoch, logs) {
#       if (epoch %% 20 == 0) cat("\n")
#       cat(".")
#     })
#   )
#
#   # Fit the model
#   model_fit <- model_dp$fit(x = X_train,
#                            y = y_train,
#                            epochs = epochs,
#                            batch_size = batch_size,
#                            validation_split = validation_split,
#                            verbose = 0, callbacks = callbacks)
#
#   return(model_dp)
# }
#
#
#
