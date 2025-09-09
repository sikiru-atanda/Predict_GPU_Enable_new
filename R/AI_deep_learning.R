
####Mapping of models

model_param_prefix <- list(
  Conv1DNet           = "cnn",
  ResNet              = "resnet",
  TabTransformer      = "ft_transformer",
  TabAttention        = "saint",
  TabNet              = "tabnet",
  LightTreeNet        = "node",
  FactorNet           = "deepfm",
  CrossNet            = "dcnv2",
  NeuralAdditive      = "nam",
  MixtureOfExperts    = "moe",
  GPNet               = "gp_dkl",
  DenseAttentionNet   = "mlp_with_attention",
  DenseNeuralNet      = "mlp"
)

### Expected hyperparameters per prefix

expected_params <- list(
  cnn = c("cnn_neurons_per_layer", "cnn_kernel_size", "cnn_dense_layers",
          "cnn_use_max_pool", "cnn_pool_kernel", "cnn_pool_stride",
          "cnn_pool_padding", "cnn_learning_rate", "cnn_separable",
          "cnn_dilations", "cnn_use_se", "cnn_norm_type", "cnn_pool_type",
          "cnn_use_global_pool"),
  resnet = c("resnet_neurons_per_layer", "n_blocks", "resnet_learning_rate"),
  ft_transformer = c("ft_d_model", "ft_heads", "ft_layers", "ft_ff_mult",
                     "ft_dropout", "ft_token_dropout", "ft_use_cls"),
  saint = c("saint_d_model", "saint_heads", "saint_layers", "saint_ff_mult",
            "saint_dropout", "saint_token_dropout", "saint_use_cls"),
  tabnet = c("tabnet_steps", "tabnet_feature_dim", "tabnet_output_dim",
             "tabnet_gamma", "tabnet_lambda_sparse"),
  node = c("node_trees", "node_depth"),
  deepfm = c("deepfm_k", "deepfm_hidden"),
  dcnv2 = c("dcn_layers", "dcn_hidden"),
  nam = c("nam_hidden", "nam_activation", "nam_add_linear", "nam_l1"),
  moe = c("moe_n_experts", "moe_expert_hidden", "moe_gate_hidden",
          "moe_temperature", "moe_sparse_topk", "moe_entropy_reg"),
  gp_dkl = c("gp_use_variational", "gp_num_inducing", "gp_feature_dim",
             "gp_kernel", "gp_ard", "gp_lr_mult",
             "rff_features", "rff_lengthscale", "rff_deep_hidden"),
  mlp = c("mlp_neurons_per_layer", "mlp_learning_rate", "heteroscedastic"),
  mlp_with_attention = c("mlp_neurons_per_layer", "mlp_learning_rate",
                         "final_attention", "attention_across_multiple_layers", "heteroscedastic")
)

#
#### Validator function

validate_model_params <- function(models, params = NULL, para_tunning = FALSE, params_env = parent.frame()) {

  if (!is.null(params) && !isTRUE(para_tunning)) {
    stop("Hyperparameter list `params` must be provided when para_tunning = FALSE.
         Or set para_tunning = TRUE to enable hyperparameter tuning.", call. = FALSE)
  }

  for (m in models) {
    prefix <- model_param_prefix[[m]]
    if (is.null(prefix)) next

    expected <- expected_params[[prefix]]
    if (is.null(expected)) next

    # Source: list for tuning or environment for normal usage
    if (isTRUE(para_tunning)) {
      if (is.null(params)) stop("params list must be provided for hyperparameter tuning", call. = FALSE)
      source <- params
      get_val <- function(x) source[[x]]
      has_val <- function(x) !is.null(source[[x]])
    } else {
      source <- NULL
      get_val <- function(x) get(x, envir = params_env, inherits = TRUE)
      has_val <- function(x) exists(x, envir = params_env, inherits = TRUE)
    }

    # ----------------------
    # Check presence
    # ----------------------
    missing <- expected[!vapply(expected, has_val, logical(1))]
    if (length(missing) > 0) {
      stop(sprintf("%s: Missing required hyperparameters: %s",
                   m, paste(missing, collapse = ", ")), call. = FALSE)
    }

    # ----------------------
    # Semantic / sanity checks
    # ----------------------
    if (prefix == "cnn" && get_val("cnn_kernel_size") <= 0) {
      stop("cnn_kernel_size must be > 0", call. = FALSE)
    }
    if (prefix == "resnet" && length(get_val("resnet_neurons_per_layer")) != get_val("n_blocks")) {
      stop("ResNet: length of resnet_neurons_per_layer must match n_blocks", call. = FALSE)
    }
    if (prefix == "tabnet" && get_val("tabnet_steps") < 1) {
      stop("TabNet: tabnet_steps must be >= 1", call. = FALSE)
    }
    if (prefix == "moe" && get_val("moe_n_experts") < 1) {
      stop("MixtureNet: moe_n_experts must be >= 1", call. = FALSE)
    }
    if (prefix == "gp_dkl" && get_val("gp_num_inducing") < 1) {
      stop("GPNet: gp_num_inducing must be >= 1", call. = FALSE)
    }

    if (prefix %in% c("mlp", "mlp_with_attention")) {
      neurons <- get_val(paste0(prefix, "_neurons_per_layer"))
      nh <- if (has_val(paste0(prefix, "_num_hidden_layers"))) get_val(paste0(prefix, "_num_hidden_layers")) else NULL
      if (!is.null(nh) && length(neurons) != nh) {
        stop(sprintf("%s: num_hidden_layers must match length of neurons_per_layer", m), call. = FALSE)
      }
    }
  }
}

dl_env <- new.env(parent = emptyenv())

get_dl_module <- function() {
  if (!is.null(dl_env$mod)) return(dl_env$mod)

  pkg <- utils::packageName()
  pyfile <- system.file("python/dl_models.py", package = pkg)

  # dev fallbacks if running via load_all()
  if (!nzchar(pyfile)) {
    candidates <- c("inst/python/dl_models.py", "python/dl_models.py", "../inst/python/dl_models.py")
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) pyfile <- normalizePath(hit[1], winslash = "/")
  }
  if (!nzchar(pyfile)) stop("dl_models.py not found under inst/python/.", call. = FALSE)

  pydir <- dirname(pyfile)
  # import as a proper Python module
  mod <- reticulate::import_from_path("dl_models", path = pydir, delay_load = FALSE, convert = TRUE)
  dl_env$mod <- mod
  mod
}

# init_dp_module <- function(prefer_gpu = TRUE) {
#   ba_set_env()
#   ensure_pydeps(engine = "dl", prefer_gpu = prefer_gpu, cuda_version = "auto")
#   m <- get_dl_module()
#   if (is.null(m)) stop("dl_models.py not found.", call. = FALSE)
#   m
# }
# usage:



# Ensure NumPy + PyTorch are available in the current reticulate interpreter
ensure_pydeps <- function(prefer_gpu = TRUE, cuda_version = c("auto","cpu","cu121","cu118")) {
  cuda_version <- match.arg(cuda_version)
  mod <- get_dl_module()
  info <- mod$setup_deps(prefer_gpu = isTRUE(prefer_gpu), cuda = cuda_version)
  invisible(info)
}
#####

# keep your existing get_dl_module()


generate_dynamic_layers <- function(input_size, num_hidden_layers, scaling_factor = 0.5, max_neurons = 1000) {
  msg <- "\n==================================================\n"
  if (length(num_hidden_layers) != 1) {
    stop(paste(msg, "num_hidden_layers should be a vector of length 1", call. = FALSE))
  }
  layers <- numeric(num_hidden_layers)
  capped_input_size <- ifelse(input_size >= 1000, max_neurons, input_size)
  layers[1] <- min(floor(capped_input_size * scaling_factor), max_neurons)
  for (i in 2:num_hidden_layers) {
    layers[i] <- min(floor(layers[i - 1] * scaling_factor), max_neurons)
    if (layers[i] < 1) break
  }
  layers[layers > 0]
}




# train_predict_deeplearning <- function(
#     data_label_geno, indices, test_geno,
#     # --- core ---
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
#
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
#     ft_scalar_tokenizer = TRUE,
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
#   # 1) Split y / X from the combined matrix
#   train_data_full <- data_label_geno[, -1, drop = FALSE]
#   y_train_full    <- data_label_geno[,  1]
#
#   # 2) Subset by bootstrap indices for training
#   X_boot <- train_data_full[indices, , drop = FALSE]
#   y_boot <- y_train_full[indices]
#
#   # 3) Per-replicate seed for reproducibility
#   seed_base <- if (is.null(random_seed)) 0L else as.integer(random_seed)
#   seed_this <- as.integer((seed_base + sum(as.integer(indices))) %% .Machine$integer.max)
#
#   # 4) Train one model on the bootstrap sample
#   fit <- try(
#     deep_learning_model_utilityy(
#       X_train = X_boot,
#       y_train = y_boot,
#
#       num_hidden_layers     = num_hidden_layers,
#       neurons_per_layer     = neurons_per_layer,
#       learning_rate         = learning_rate,
#       epochs                = epochs,
#       batch_size            = batch_size,
#       l2_regularizer_dp     = l2_regularizer_dp,
#       dropout_rate          = dropout_rate,
#       validation_split      = validation_split,
#       n_blocks              = n_blocks,
#       n_neurons_per_block   = n_neurons_per_block,
#       deep_learning_model   = deep_learning_model,
#       attention_on_final_layer       = attention_on_final_layer,
#       attention_across_multiple_layers = attention_across_multiple_layers,
#       batch_normalization   = batch_normalization,
#       dense_layers_cnn      = dense_layers_cnn,
#       kernel_size           = kernel_size,
#
#       # Compile / device / seeds
#       compile_model = compile_model,
#       deterministic = deterministic,
#       random_seed   = seed_this,
#       device        = device,
#
#       # --- NEW knobs passed through to Python ---
#       # CNN max-pool
#       cnn_use_max_pool   = cnn_use_max_pool,
#       cnn_pool_kernel    = as.integer(cnn_pool_kernel),
#       cnn_pool_stride    = as.integer(cnn_pool_stride),
#       cnn_pool_padding   = as.integer(cnn_pool_padding),
#       separable=isTRUE(separable),
#       dilations=as.integer(dilations),
#       use_se=isTRUE(use_se),
#       norm_type=norm_type,
#       use_max_pool=isTRUE(use_max_pool),
#       pool_type=pool_type,
#       use_global_pool=isTRUE(use_global_pool),
#
#       # Optim / AMP / grad-clip / class weights
#       optimizer_name     = as.character(optimizer_name),
#       auto_class_weights = isTRUE(auto_class_weights),
#       use_amp            = isTRUE(use_amp),
#       max_grad_norm      = if (is.null(max_grad_norm)) NULL else as.numeric(max_grad_norm),
#
#       # Heteroscedastic (regression only)
#       heteroscedastic    = isTRUE(heteroscedastic),
#
#       # FT-Transformer
#       ft_d_model       = as.integer(ft_d_model),
#       ft_heads         = as.integer(ft_heads),
#       ft_layers        = as.integer(ft_layers),
#       ft_ff_mult       = as.integer(ft_ff_mult),
#       ft_dropout       = as.numeric(ft_dropout),
#       ft_token_dropout = as.numeric(ft_token_dropout),
#       ft_use_cls       = isTRUE(ft_use_cls),
#       ft_scalar_tokenizer = isTRUE(ft_scalar_tokenizer),
#
#       # SAINT
#       saint_d_model       = as.integer(saint_d_model),
#       saint_heads         = as.integer(saint_heads),
#       saint_layers        = as.integer(saint_layers),
#       saint_ff_mult       = as.integer(saint_ff_mult),
#       saint_dropout       = as.numeric(saint_dropout),
#       saint_token_dropout = as.numeric(saint_token_dropout),
#       saint_use_cls       = isTRUE(saint_use_cls),
#
#       # Grouping (FT/SAINT/NAM/MoE)
#       use_grouping   = isTRUE(use_grouping),
#       group_trigger  = as.integer(group_trigger),
#       group_method   = as.character(group_method),
#       init_group_size= as.integer(init_group_size),
#       max_tokens     = as.integer(max_tokens),
#       kmeans_batch   = as.integer(kmeans_batch),
#       kmeans_iter    = as.integer(kmeans_iter),
#
#       # TabNet
#       tabnet_steps         = as.integer(tabnet_steps),
#       tabnet_feature_dim   = as.integer(tabnet_feature_dim),
#       tabnet_output_dim    = as.integer(tabnet_output_dim),
#       tabnet_gamma         = as.numeric(tabnet_gamma),
#       tabnet_lambda_sparse = as.numeric(tabnet_lambda_sparse),
#
#       # NODE
#       node_trees = as.integer(node_trees),
#       node_depth = as.integer(node_depth),
#
#       # DeepFM
#       deepfm_k      = as.integer(deepfm_k),
#       deepfm_hidden = as.integer(deepfm_hidden),
#
#       # DCN
#       dcn_layers = as.integer(dcn_layers),
#       dcn_hidden = as.integer(dcn_hidden),
#
#       # NAM
#       nam_hidden     = as.integer(nam_hidden),
#       nam_activation = as.character(nam_activation),
#       nam_add_linear = isTRUE(nam_add_linear),
#       nam_l1         = as.numeric(nam_l1),
#
#       # MoE
#       moe_n_experts     = as.integer(moe_n_experts),
#       moe_expert_hidden = as.integer(moe_expert_hidden),
#       moe_gate_hidden   = as.integer(moe_gate_hidden),
#       moe_temperature   = as.numeric(moe_temperature),
#       moe_sparse_topk   = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#       moe_entropy_reg   = as.numeric(moe_entropy_reg),
#
#       # GP-DKL (+ RFF fallback)
#       gp_use_variational = isTRUE(gp_use_variational),
#       gp_num_inducing    = as.integer(gp_num_inducing),
#       gp_feature_dim     = as.integer(gp_feature_dim),
#       gp_kernel          = as.character(gp_kernel),
#       gp_ard             = isTRUE(gp_ard),
#       gp_lr_mult         = as.numeric(gp_lr_mult),
#       rff_features       = as.integer(rff_features),
#       rff_lengthscale    = as.numeric(rff_lengthscale),
#       rff_deep_hidden    = as.integer(rff_deep_hidden)
#     ),
#     silent = TRUE
#   )
#
#   # 5) On failure, return NAs of correct length so boot() keeps going
#   if (inherits(fit, "try-error") || is.null(fit) || length(fit) < 1L) {
#     n_out <- if (!is.null(test_geno)) nrow(test_geno) else nrow(train_data_full)
#     return(rep(NA_real_, n_out))
#   }
#
#   model <- fit[[1L]]  # (model, history)
#
#   # 6) Predict on requested target (force vector output)
#   preds <- try(
#     {
#       if (!is.null(test_geno)) {
#         torch_predict(model, as.matrix(test_geno), device = device, type = "class")
#       } else {
#         torch_predict(model, as.matrix(train_data_full), device = device, type = "class")
#       }
#     },
#     silent = TRUE
#   )
#
#   if (inherits(preds, "try-error") || is.null(preds)) {
#     n_out <- if (!is.null(test_geno)) nrow(test_geno) else nrow(train_data_full)
#     return(rep(NA_real_, n_out))
#   }
#
#   as.numeric(preds)
# }




# get_best_model <- function(results, hyperparam_combinations) {
#   best_val <- Inf
#   best_model <- NULL
#   best_hparams <- NULL
#
#   for (i in seq_along(results)) {
#     key <- names(results)[i]
#     hist <- results[[key]]$history
#     # Py side returns history with train_loss / val_loss
#     if (is.null(hist$val_loss)) {
#       warning("History missing 'val_loss' for key: ", key)
#       next
#     }
#     current <- min(unlist(hist$val_loss))
#     if (is.finite(current) && current < best_val) {
#       best_val <- current
#       best_model <- results[[key]]$model
#       best_hparams <- hyperparam_combinations[i, , drop = FALSE]
#     }
#   }
#
#   list(best_model = best_model, best_hyperparameters = best_hparams)
# }
get_best_model <- function(
    results,
    hyperparam_combinations = NULL,
    metric     = c("val_loss", "train_loss"),
    direction  = c("min", "max"),
    agg        = c("min", "last", "mean", "median"),
    na_rm      = TRUE
) {
  metric    <- match.arg(metric)
  direction <- match.arg(direction)
  agg       <- match.arg(agg)

  if (!length(results)) stop("Empty results list.")

  # --- helpers -------------------------------------------------------
  # Convert reticulate objects or R lists -> numeric vector
  .to_numeric_vec <- function(x) {
    if (is.null(x)) return(numeric(0))
    # If it's a Python object, try py_to_r
    if (inherits(x, "python.builtin.object") || inherits(x, "numpy.ndarray")) {
      if (requireNamespace("reticulate", quietly = TRUE)) {
        x <- reticulate::py_to_r(x)
      } else {
        return(numeric(0))
      }
    }
    v <- suppressWarnings(as.numeric(unlist(x, recursive = TRUE, use.names = FALSE)))
    if (length(v) == 0L) return(numeric(0))
    if (na_rm) v <- v[is.finite(v)]
    v
  }

  # Safely grab a metric vector from history (fallback to train_loss if needed)
  .hist_metric <- function(hist, key) {
    if (is.null(hist)) return(numeric(0))
    # reticulate dict/list both support [[
    v <- try(hist[[key]], silent = TRUE)
    if (!inherits(v, "try-error") && !is.null(v)) {
      out <- .to_numeric_vec(v)
      if (length(out)) return(out)
    }
    if (identical(key, "val_loss")) {
      v2 <- try(hist[["train_loss"]], silent = TRUE)
      if (!inherits(v2, "try-error") && !is.null(v2)) {
        out2 <- .to_numeric_vec(v2)
        if (length(out2)) return(out2)
      }
    }
    numeric(0)
  }

  # Score a single history
  score_of <- function(hist) {
    vec <- .hist_metric(hist, metric)
    if (!length(vec)) return(list(score = NA_real_, epoch = NA_integer_))
    sc <- switch(
      agg,
      min    = min(vec),
      last   = tail(vec, 1),
      mean   = mean(vec),
      median = median(vec)
    )
    ep <- switch(
      agg,
      min    = which.min(vec),
      last   = length(vec),
      mean   = which.min(abs(vec - sc)),
      median = which.min(abs(vec - sc))
    )
    list(score = as.numeric(sc), epoch = as.integer(ep))
  }

  # Sanitize and flatten a config: drop huge/heavy fields and keep scalars/small vectors
  .sanitize_config <- function(cfg) {
    if (is.null(cfg)) return(list())
    drop_keys <- c("X","y","ft_group_index","ft_group_counts","feature_groups")
    keep <- setdiff(names(cfg), drop_keys)
    out <- list()
    for (k in keep) {
      v <- cfg[[k]]
      if (is.null(v)) next
      # unwrap python objects if possible
      if (inherits(v, "python.builtin.object") || inherits(v, "numpy.ndarray")) {
        if (requireNamespace("reticulate", quietly = TRUE)) v <- reticulate::py_to_r(v) else next
      }
      # skip matrices/arrays/data.frames/lists of complex stuff
      if (is.matrix(v) || is.array(v) || is.data.frame(v)) next
      if (is.list(v)) {
        # keep only small simple lists
        vv <- unlist(v, recursive = TRUE, use.names = FALSE)
        if (!length(vv)) next
        v <- vv
      }
      # coerce vectors to a compact character
      if (length(v) > 20) {
        v <- paste0(paste(head(as.character(v), 8), collapse = ","), ",…")
      } else if (length(v) > 1) {
        v <- paste(as.character(v), collapse = ",")
      }
      # allow logical/numeric/character scalars/strings
      if (length(v) == 1) v <- as.character(v)
      out[[k]] <- v
    }
    out
  }

  # Bind list of named lists (possibly different keys) to a data.frame of characters
  rows_to_df <- function(lst) {
    if (!length(lst)) return(data.frame())
    keys <- unique(unlist(lapply(lst, names)))
    if (!length(keys)) return(data.frame())
    df <- data.frame(matrix(NA_character_, nrow = length(lst), ncol = length(keys)),
                     stringsAsFactors = FALSE)
    names(df) <- keys
    for (i in seq_along(lst)) {
      if (is.null(lst[[i]]) || !length(lst[[i]])) next
      nm <- names(lst[[i]])
      df[i, nm] <- as.character(unlist(lst[[i]], use.names = FALSE))
    }
    df
  }

  # --- compute scores -------------------------------------------------
  n <- length(results)
  scores <- rep(NA_real_, n)
  epochs <- rep(NA_integer_, n)
  cfg_rows <- vector("list", n)

  for (i in seq_len(n)) {
    r  <- results[[i]]
    sc <- try(score_of(r$history), silent = TRUE)
    if (!inherits(sc, "try-error") && length(sc)) {
      scores[i] <- sc$score
      epochs[i] <- sc$epoch
    }
    if (!is.null(hyperparam_combinations)) {
      cfg_rows[[i]] <- as.list(hyperparam_combinations[i, , drop = FALSE])
    } else if (!is.null(r$config)) {
      cfg_rows[[i]] <- .sanitize_config(r$config)
    } else {
      cfg_rows[[i]] <- list()
    }
  }

  ok <- which(is.finite(scores))
  if (!length(ok)) stop("No valid histories with metric '", metric, "' were found.")

  best_i <- if (direction == "min") ok[which.min(scores[ok])] else ok[which.max(scores[ok])]

  # --- leaderboard ----------------------------------------------------
  leaderboard <- data.frame(
    index = ok,
    score = scores[ok],
    epoch = epochs[ok],
    stringsAsFactors = FALSE
  )
  ord <- if (direction == "min") order(leaderboard$score, leaderboard$epoch) else order(-leaderboard$score, leaderboard$epoch)
  leaderboard <- leaderboard[ord, , drop = FALSE]

  # Attach hyperparams
  if (!is.null(hyperparam_combinations)) {
    hp_df <- hyperparam_combinations[leaderboard$index, , drop = FALSE]
    rownames(hp_df) <- NULL
    leaderboard <- cbind(leaderboard, hp_df)
  } else {
    cfg_ok  <- cfg_rows[leaderboard$index]
    hp_rows <- lapply(cfg_ok, identity)
    hp_df   <- rows_to_df(hp_rows)
    if (nrow(hp_df)) leaderboard <- cbind(leaderboard, hp_df)
  }

  # --- output ---------------------------------------------------------
  list(
    best_model           = results[[best_i]]$model,
    best_history         = results[[best_i]]$history,
    best_index           = best_i,
    best_score           = scores[best_i],
    best_epoch           = epochs[best_i],
    best_hyperparameters = if (!is.null(hyperparam_combinations)) {
      hyperparam_combinations[best_i, , drop = FALSE]
    } else {
      as.data.frame(t(unlist(.sanitize_config(results[[best_i]]$config))), stringsAsFactors = FALSE)
    },
    leaderboard          = leaderboard
  )
}

#####
# deep_learning_model <- function(pheno_object = NULL,
#                                 y = NULL,
#                                 omics_data = NULL,
#                                 crossval = FALSE,
#                                 tst = NULL,
#                                 geno_omic_object = NULL,
#                                 geno_omic_test_object = NULL,
#                                 response = NULL,
#                                 gen_name = NULL,
#                                 message = TRUE,
#                                 scaling = TRUE,
#                                 centering = FALSE,
#                                 omic_count = NULL,
#                                 #num_hidden_layers = 1,
#                                 #neurons_per_layer = NULL,
#                                 #learning_rate_dp = 0.001,
#                                 #epochs = 10,
#                                 #batch_size = 32,
#                                 para_tunning = FALSE,
#                                 param_grid = NULL,
#                                 #validation_split = 0.2,
#                                 early_stop = TRUE,
#                                 #l2_regularizer_dp = 0.001,
#                                 #dropout_rate = 0.5,
#                                 deep_learning_model = "mlp_with_attention",
#                                 #n_blocks = 2,
#                                 #dense_layers_cnn = c(128, 64),
#                                 #kernel_size = 3,
#                                 #n_neurons_per_block = NULL,
#                                 #attention_on_final_layer = TRUE,
#                                 #attention_across_multiple_layers = FALSE,
#                                 #batch_normalization = TRUE,
#                                 CI_width_thresholds = c(0.33, 0.66),
#                                 high_reliability_thres = 0.9,
#                                 low_reliability_thres = 0.5,
#                                 n_components = 20,
#                                 threshold = 100,
#                                 target = "test_set",
#                                 iqr_multiplier = 1.5,
#                                 interval_width_high_threshold = NULL,
#                                 interval_width_low_threshold = NULL,
#                                 interval_width_moderate_threshold = NULL,
#                                 n_bootstrap = 100,
#                                 system_database = FALSE,
#                                 compile_model = TRUE,
#                                 deterministic = FALSE,
#                                 random_seed = NULL,
#                                 device = NULL,
#
#                                 optimizer_name = "adam",
#                                 use_amp        = TRUE,
#                                 max_grad_norm  = 1.0,
#                                 auto_class_weights = FALSE,
#                                 ##### cnn
#                                 cnn_neurons_per_layer = as.integer(c(64, 64, 64)),
#                                 cnn_kernel_size = 3L,
#                                 cnn_dense_layers = as.integer(c(256, 128, 64)),
#                                 cnn_use_max_pool = FALSE,
#                                 cnn_pool_kernel = 2L,
#                                 cnn_pool_stride = 2L,
#                                 cnn_pool_padding = 0L,
#                                 cnn_learning_rate = 1e-3,
#                                 cnn_separable=TRUE,
#                                 cnn_dilations=c(1,2,4),
#                                 cnn_use_se=TRUE,
#                                 cnn_norm_type="group",
#                                 cnn_pool_type="conv",
#                                 cnn_use_global_pool=FALSE,
#                                 ##### resnet
#                                 resnet_neurons_per_block = as.integer(c(256, 128, 64)),
#                                 resnet_blocks = 3,
#                                 resnet_learning_rate = 1e-3,
#                                 #### ft_transformer
#                                 ft_d_model = 192L,
#                                 ft_heads = 8L,
#                                 ft_layers = 3L,
#                                 ft_ff_mult = 4L,
#                                 ft_dropout = 0.1,
#                                 ft_token_dropout = 0.0,
#                                 ft_use_cls = TRUE,
#                                 #### saint
#                                 saint_d_model = 128L,
#                                 saint_heads = 8L,
#                                 saint_layers = 3L,
#                                 saint_ff_mult = 4L,
#                                 saint_dropout = 0.1,
#                                 saint_token_dropout = 0.0,
#                                 saint_use_cls       = TRUE,
#                                 ###### Grouping controls (FT/SAINT and NAM/MoE)
#                                 use_grouping   = FALSE,
#                                 group_trigger  = 2048,
#                                 group_method   = "auto",
#                                 init_group_size = 64,
#                                 max_tokens      = 1024,
#                                 kmeans_batch    = 4096,
#                                 kmeans_iter     = 100,
#                                 #### tabnet
#                                 tabnet_steps = 5L,
#                                 tabnet_feature_dim = 64L,
#                                 tabnet_output_dim = 64L,
#                                 tabnet_gamma = 1.5,
#                                 tabnet_lambda_sparse = 1e-4,
#                                 #### node
#                                 node_trees = 8L,
#                                 node_depth = 3L,
#                                 #### deepfm
#                                 deepfm_k = 16L,
#                                 deepfm_hidden = as.integer(c(128, 64)),
#                                 #### dcnv2
#                                 dcn_layers = 3L,
#                                 dcn_hidden = as.integer(c(256, 128, 64)),
#                                 #### nam
#                                 nam_hidden = as.integer(c(32, 16)),
#                                 nam_activation = "relu",
#                                 nam_add_linear = TRUE,
#                                 nam_l1 = 1e-4,
#                                 ### moe
#                                 moe_n_experts = 4L,
#                                 moe_expert_hidden = as.integer(c(128, 64)),
#                                 moe_gate_hidden = 128L,
#                                 moe_temperature = 1.0,
#                                 moe_sparse_topk = NA,
#                                 moe_entropy_reg = 0.0,
#                                 ### gp_dkl/RFF knobs
#                                 gp_use_variational = TRUE,
#                                 gp_num_inducing = 256L,
#                                 gp_feature_dim = 64L,
#                                 gp_kernel = "rbf",
#                                 gp_ard = TRUE,
#                                 gp_lr_mult = 0.5,
#                                 rff_features       = 1024,
#                                 rff_lengthscale    = 1.0,
#                                 rff_deep_hidden    = c(128),
#                                 ##### General dp
#                                 model_type = "resnet",
#                                 epochs = 10,
#                                 batch_size = 64 ,
#                                 dropout = 0.2,
#                                 l2_weight_decay = 1e-4,
#                                 l2_regularizer_dp = 0.001,
#                                 dropout_rate = 0.5,
#                                 batch_norm = TRUE,
#                                 validation_split = 0.2,
#                                 compile_model = FALSE,
#                                 deterministic = TRUE,
#                                 random_seed = 123,
#                                 device = NULL,
#                                 #### mlp and attention
#                                 mlp_neurons_per_layer = as.integer(c(128, 64)),
#                                 mlp_learning_rate = 1e-3,
#                                 final_attention = TRUE,
#                                 attention_across_multiple_layers = TRUE,
#                                 heteroscedastic = TRUE,
#                                 ...) {
#
#
#   msg <- "\n==================================================\n"
#
#   # 0) Ensure deps + load module (prefer GPU unless device=='cpu')
#   prefer_gpu <- !identical(device, "cpu")
#   ensure_pydeps(prefer_gpu = prefer_gpu, cuda_version = "auto")
#   mod <- get_dl_module()
#
#   # 1) Inputs and scaling
#   if (is.null(geno_omic_object) && is.null(pheno_object) && isFALSE(crossval)) {
#     stop(paste(msg, "provide matrix of the predictors and the data.frame of the Y variable."), call. = FALSE)
#   }
#
#   if (!is.null(geno_omic_object)) {
#     GID <- rownames(geno_omic_object)
#     scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
#     geno_omic_object <- stats::predict(scaler, geno_omic_object)
#     cols_with_na <- which(colSums(is.na(geno_omic_object)) > 0)
#     if (length(cols_with_na) != 0) {
#       geno_omic_object <- geno_omic_object[, -cols_with_na, drop = FALSE]
#       if (!is.null(geno_omic_test_object)) {
#         geno_omic_test_object <- geno_omic_test_object[, -cols_with_na, drop = FALSE]
#       }
#     }
#   }
#
#   if (!is.null(geno_omic_test_object)) {
#     test_label <- rownames(geno_omic_test_object)
#     geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
#     GID <- rownames(geno_omic_test_object)
#     geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
#   }
#
#   if (!is.null(omics_data) && isTRUE(crossval)) {
#     scaler <- caret::preProcess(omics_data, method = c("center", "scale"))
#     omics_data <- stats::predict(scaler, omics_data)
#     cols_with_na <- which(colSums(is.na(omics_data)) > 0)
#     if (length(cols_with_na) != 0) {
#       omics_data <- omics_data[, -cols_with_na, drop = FALSE]
#     }
#   }
#
#   if (isTRUE(crossval)) {
#     para_tunning <- FALSE
#     param_grid <- NULL
#   }
#
#   # Helper: canonicalize model type -> Python
#   # model_type <- switch(tolower(deep_learning_model),
#   #                      "mlp"               = "mlp",
#   #                      "mlp_with_attention"= "mlp_with_attention",
#   #                      "resnet"            = "resnet",
#   #                      "cnn"               = "cnn",
#   #                      "ft_transformer"    = "ft_transformer", "ft" = "ft_transformer",
#   #                      "saint"             = "saint",
#   #                      "tabnet"            = "tabnet",
#   #                      "node"              = "node",
#   #                      "deepfm"            = "deepfm",
#   #                      "dcn"               = "dcnv2", "dcnv2" = "dcnv2",
#   #                      "nam"               = "nam",
#   #                      "moe"               = "moe",
#   #                      "gp_dkl"            = "gp_dkl",
#   #                      stop("Unsupported deep_learning_model: ", deep_learning_model)
#   # )
#   #
#   # # 2) Hidden-layer requirements (classic dense/CNN; ResNet handled below)
#   # wants_dense <- model_type %in% c("mlp_with_attention", "mlp", "cnn")  # ResNet separate
#   # if (wants_dense && is.null(num_hidden_layers)) {
#   #   stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#   # }
#   # if (identical(model_type, "resnet") && is.null(n_blocks)) {
#   #   stop(paste(msg, "n_blocks can't be NULL."), call. = FALSE)
#   # }
#   # if (is.null(neurons_per_layer) && !is.null(num_hidden_layers) && (wants_dense || identical(model_type, "resnet"))) {
#   #   if (!is.null(geno_omic_object) && isFALSE(crossval)) {
#   #     neurons_per_layer <- generate_dynamic_layers(ncol(geno_omic_object), num_hidden_layers, 0.5)
#   #   } else if (!is.null(omics_data) && isTRUE(crossval)) {
#   #     neurons_per_layer <- generate_dynamic_layers(ncol(omics_data), num_hidden_layers, 0.5)
#   #   } else {
#   #     stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#   #   }
#   # }
#   # if (identical(model_type, "resnet")) {
#   #   if (is.null(n_blocks)) n_blocks <- 1L
#   #       # n_neurons_per_block <-  neurons_per_layer
#   #       # num_hidden_layers <- n_blocks
#   #   if (is.null(n_neurons_per_block)) {
#   #     if (!is.null(geno_omic_object) && isFALSE(crossval)) {
#   #       n_neurons_per_block <- generate_dynamic_layers(ncol(geno_omic_object), n_blocks, 0.5)
#   #     } else if (!is.null(omics_data) && isTRUE(crossval)) {
#   #       n_neurons_per_block <- generate_dynamic_layers(ncol(omics_data), n_blocks, 0.5)
#   #     } else {
#   #       stop(paste(msg, "n_neurons_per_block and n_blocks can't be NULL."), call. = FALSE)
#   #     }
#   #   }
#   #   if (n_blocks != length(n_neurons_per_block)) {
#   #     stop(paste(msg, "Mismatch in hidden number of block and neurons per layer"), call. = FALSE)
#   #   }
#   #   neurons_per_layer <- n_neurons_per_block
#   #   num_hidden_layers <- n_blocks
#   # } else if (wants_dense) {
#   #   if (is.null(num_hidden_layers)) num_hidden_layers <- 1L
#   #   if (is.null(neurons_per_layer)) {
#   #     if (!is.null(geno_omic_object) && isFALSE(crossval)) {
#   #       neurons_per_layer <- generate_dynamic_layers(ncol(geno_omic_object), num_hidden_layers, 0.5)
#   #     } else if (!is.null(omics_data) && isTRUE(crossval)) {
#   #       neurons_per_layer <- generate_dynamic_layers(ncol(omics_data), num_hidden_layers, 0.5)
#   #     } else {
#   #       stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#   #     }
#   #   }
#   # }
#
#   # Normalize scalars
#   num_hidden_layers <- if (!is.null(num_hidden_layers)) as.integer(num_hidden_layers) else NULL
#   neurons_per_layer <- if (!is.null(neurons_per_layer)) as.integer(neurons_per_layer) else NULL
#   batch_size        <- as.integer(if (!is.null(batch_size)) batch_size else 30L)
#   epochs            <- as.integer(if (!is.null(epochs)) epochs else 10L)
#   dropout_rate      <- as.numeric(if (!is.null(dropout_rate)) dropout_rate else 0.5)
#   learning_rate     <- as.numeric(if (!is.null(learning_rate_dp)) learning_rate_dp else 0.01)
#   l2_regularizer_dp <- as.numeric(if (!is.null(l2_regularizer_dp)) l2_regularizer_dp else 0)
#
#   separable <- isTRUE(cnn_separable)
#   dilations <- as.integer(cnn_dilations)
#   use_se <- isTRUE(cnn_use_se)
#   norm_type <- cnn_norm_type
#   #use_max_pool <- isTRUE(use_max_pool)
#   pool_type <- cnn_pool_type
#   use_global_pool <- isTRUE(cnn_use_global_pool)
#   # 3) Response & y scaling
#   if (isFALSE(crossval)) {
#     y_train <- as.numeric(pheno_object[, response])
#   }
#   y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
#   y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
#   data_label_geno <- cbind(y_train_scaled, geno_omic_object)
#
#   # ---------- HYPERPARAMETER TUNING (kept as in your original, but now respects extras) ----------
#   if (isTRUE(para_tunning) && isFALSE(crossval)) {
#     if (is.null(param_grid)) stop(paste(msg, "param_grid must be provided when tuning is enabled."), call. = FALSE)
#
#     # Allow ResNet aliases in param_grid
#     if (identical(model_type, "resnet")) {
#       if (!is.null(param_grid$n_blocks)) names(param_grid)[names(param_grid) == "n_blocks"] <- "num_hidden_layers"
#       if (!is.null(param_grid$n_neurons_per_block)) names(param_grid)[names(param_grid) == "n_neurons_per_block"] <- "neurons_per_layer"
#     }
#
#     combos <- expand.grid(param_grid, stringsAsFactors = FALSE)
#     parse_neurons <- function(x) {
#       if (is.null(x)) return(NULL)
#       if (is.numeric(x)) return(as.integer(x))
#       if (is.list(x))    return(as.integer(unlist(x)))
#       if (is.character(x)) return(as.integer(strsplit(x, ",")[[1]]))
#       NULL
#     }
#
#     # Validate (only for models that use layered widths)
#     if ("num_hidden_layers" %in% names(combos) && "neurons_per_layer" %in% names(combos) &&
#         (wants_dense || identical(model_type, "resnet"))) {
#       ok <- vapply(seq_len(nrow(combos)), function(i) {
#         nl <- suppressWarnings(as.integer(combos$num_hidden_layers[i]))
#         np <- parse_neurons(combos$neurons_per_layer[i][[1]])
#         if (is.null(nl) || is.na(nl) || is.null(np)) return(FALSE)
#         length(np) == nl
#       }, logical(1))
#       combos <- combos[ok, , drop = FALSE]
#       if (nrow(combos) == 0) stop(paste(msg, "No valid hyperparameter combinations after validation."), call. = FALSE)
#     }
#
#     # Simple internal tuner on val_loss
#     best_idx <- NA_integer_; best_score <- Inf; tried <- 0L
#     for (i in seq_len(nrow(combos))) {
#       hp <- combos[i, , drop = FALSE]
#       nl  <- if ("num_hidden_layers" %in% names(hp)) as.integer(hp$num_hidden_layers) else num_hidden_layers
#       npl <- if ("neurons_per_layer" %in% names(hp)) parse_neurons(hp$neurons_per_layer[[1]]) else neurons_per_layer
#       if ((wants_dense || identical(model_type, "resnet")) && !is.null(nl) && !is.null(npl) && length(npl) != nl) next
#
#       # Core overrides
#       lr  <- if ("learning_rate" %in% names(hp)) as.numeric(hp$learning_rate) else learning_rate
#       dr  <- if ("dropout_rate"  %in% names(hp)) as.numeric(hp$dropout_rate)  else dropout_rate
#       ep  <- if ("epochs"        %in% names(hp)) as.integer(hp$epochs)        else epochs
#       l2w <- if ("l2_regularizer_dp" %in% names(hp)) as.numeric(hp$l2_regularizer_dp) else l2_regularizer_dp
#       ks  <- if ("kernel_size"   %in% names(hp)) as.integer(hp$kernel_size)   else kernel_size
#       dlc <- if ("dense_layers_cnn" %in% names(hp)) parse_neurons(hp$dense_layers_cnn[[1]]) else as.integer(dense_layers_cnn)
#
#       # Fit one combo (deterministic for fairness)
#       res <- try({
#         mod$fit_model(
#           X = as.matrix(if (!is.null(geno_omic_object)) geno_omic_object else omics_data),
#           y = as.numeric(if (!is.null(pheno_object)) pheno_object[, response] else y),
#
#           model_type = model_type,
#           num_hidden_layers = nl,
#           neurons_per_layer = if (!is.null(npl)) as.integer(npl) else NULL,
#
#           learning_rate = lr,
#           epochs = ep,
#           batch_size = as.integer(batch_size),
#           l2_weight_decay = l2w,
#           dropout = dr,
#           optimizer_name = optimizer_name,
#
#           final_attention = isTRUE(attention_on_final_layer),
#           attention_across_multiple_layers = isTRUE(attention_across_multiple_layers),
#           batch_norm = isTRUE(batch_normalization),
#
#           validation_split = as.numeric(validation_split),
#           compile_model = isTRUE(compile_model),
#           device = if (is.null(device)) NULL else as.character(device),
#
#           deterministic = TRUE,
#           random_seed = as.integer(if (is.null(random_seed)) 123L else random_seed),
#
#           # AMP / grad-clip / class weights
#           use_amp = isTRUE(use_amp),
#           max_grad_norm = as.numeric(max_grad_norm),
#           auto_class_weights = isTRUE(auto_class_weights),
#
#           # CNN
#           kernel_size = ks,
#           dense_layers_cnn = as.integer(dlc),
#           cnn_use_max_pool = isTRUE(cnn_use_max_pool),
#           cnn_pool_kernel  = as.integer(cnn_pool_kernel),
#           cnn_pool_stride  = as.integer(cnn_pool_stride),
#           cnn_pool_padding = as.integer(cnn_pool_padding),
#           separable=isTRUE(separable),
#           dilations=as.integer(dilations),
#           use_se=isTRUE(use_se),
#           norm_type=norm_type,
#           #use_max_pool=isTRUE(use_max_pool),
#           pool_type=pool_type,
#           use_global_pool=isTRUE(use_global_pool),
#
#           # Heteroscedastic
#           heteroscedastic = isTRUE(heteroscedastic),
#
#           # FT / SAINT (+ grouping)
#           ft_d_model = as.integer(ft_d_model),
#           ft_heads   = as.integer(ft_heads),
#           ft_layers  = as.integer(ft_layers),
#           ft_ff_mult = as.integer(ft_ff_mult),
#           ft_dropout = as.numeric(ft_dropout),
#           ft_token_dropout = as.numeric(ft_token_dropout),
#           ft_use_cls = isTRUE(ft_use_cls),
#
#           saint_d_model = as.integer(saint_d_model),
#           saint_heads   = as.integer(saint_heads),
#           saint_layers  = as.integer(saint_layers),
#           saint_ff_mult = as.integer(saint_ff_mult),
#           saint_dropout = as.numeric(saint_dropout),
#           saint_token_dropout = as.numeric(saint_token_dropout),
#           saint_use_cls = isTRUE(saint_use_cls),
#
#           use_grouping = isTRUE(use_grouping),
#           group_trigger = as.integer(group_trigger),
#           group_method  = as.character(group_method),
#           init_group_size = as.integer(init_group_size),
#           max_tokens      = as.integer(max_tokens),
#           kmeans_batch    = as.integer(kmeans_batch),
#           kmeans_iter     = as.integer(kmeans_iter),
#
#           # TabNet
#           tabnet_steps         = as.integer(tabnet_steps),
#           tabnet_feature_dim   = as.integer(tabnet_feature_dim),
#           tabnet_output_dim    = as.integer(tabnet_output_dim),
#           tabnet_gamma         = as.numeric(tabnet_gamma),
#           tabnet_lambda_sparse = as.numeric(tabnet_lambda_sparse),
#
#           # NODE
#           node_trees = as.integer(node_trees),
#           node_depth = as.integer(node_depth),
#
#           # DeepFM / DCN
#           deepfm_k      = as.integer(deepfm_k),
#           deepfm_hidden = as.integer(deepfm_hidden),
#
#           dcn_layers = as.integer(dcn_layers),
#           dcn_hidden = as.integer(dcn_hidden),
#
#           # NAM
#           nam_hidden     = as.integer(nam_hidden),
#           nam_activation = as.character(nam_activation),
#           nam_add_linear = isTRUE(nam_add_linear),
#           nam_l1         = as.numeric(nam_l1),
#
#           # MoE
#           moe_n_experts     = as.integer(moe_n_experts),
#           moe_expert_hidden = as.integer(moe_expert_hidden),
#           moe_gate_hidden   = as.integer(moe_gate_hidden),
#           moe_temperature   = as.numeric(moe_temperature),
#           moe_sparse_topk   = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#           moe_entropy_reg   = as.numeric(moe_entropy_reg),
#
#           # GP-DKL / RFF
#           gp_use_variational = isTRUE(gp_use_variational),
#           gp_num_inducing    = as.integer(gp_num_inducing),
#           gp_feature_dim     = as.integer(gp_feature_dim),
#           gp_kernel          = as.character(gp_kernel),
#           gp_ard             = isTRUE(gp_ard),
#           gp_lr_mult         = as.numeric(gp_lr_mult),
#           rff_features       = as.integer(rff_features),
#           rff_lengthscale    = as.numeric(rff_lengthscale),
#           rff_deep_hidden    = as.integer(rff_deep_hidden)
#         )
#       }, silent = TRUE)
#
#       if (inherits(res, "try-error")) next
#       tried <- tried + 1L
#       hist <- res[[2L]]
#       vloss <- try(suppressWarnings(as.numeric(hist$val_loss)), silent = TRUE)
#       if (inherits(vloss, "try-error") || length(vloss) == 0 || all(!is.finite(vloss))) next
#       score <- min(vloss, na.rm = TRUE)
#       if (is.finite(score) && score < best_score) { best_score <- score; best_idx <- i }
#     }
#
#     if (!is.na(best_idx) && message) {
#       lock_best_hp(combos, best_idx, env = environment())
#       message("Tuning finished: best val_loss = ", signif(best_score, 6))
#     }else {
#       warning(paste(msg, "No valid hyperparameter combination could be trained; falling back to provided defaults."))
#     }
#   }
#   # ---------- END TUNING ----------
#
#   # 4) CROSS-VALIDATION BRANCH
#   if (isTRUE(crossval)) {
#     out <- tryCatch({
#       y_train_cv <- as.numeric(y)
#       y_train_cv <- y_train_cv[-tst]
#
#       py_fit <- mod$fit_model(
#         X = as.matrix(omics_data[-tst, , drop = FALSE]),
#         y = y_train_cv,
#
#         model_type = model_type,
#         num_hidden_layers = num_hidden_layers,
#         neurons_per_layer = if (!is.null(neurons_per_layer)) as.integer(neurons_per_layer) else NULL,
#
#         learning_rate = learning_rate_dp,
#         epochs = epochs,
#         batch_size = batch_size,
#         l2_weight_decay = l2_regularizer_dp,
#         dropout = dropout_rate,
#         optimizer_name = optimizer_name,
#
#         final_attention = isTRUE(attention_on_final_layer),
#         attention_across_multiple_layers = isTRUE(attention_across_multiple_layers),
#         batch_norm = isTRUE(batch_normalization),
#
#         validation_split = as.numeric(validation_split),
#         compile_model = isTRUE(compile_model),
#         device = if (is.null(device)) NULL else as.character(device),
#
#         deterministic = isTRUE(deterministic),
#         random_seed   = if (is.null(random_seed)) NULL else as.integer(random_seed),
#
#         # AMP / grad-clip / class weights
#         use_amp = isTRUE(use_amp),
#         max_grad_norm = as.numeric(max_grad_norm),
#         auto_class_weights = isTRUE(auto_class_weights),
#
#         # CNN
#         kernel_size = as.integer(kernel_size),
#         dense_layers_cnn = as.integer(dense_layers_cnn),
#         cnn_use_max_pool = isTRUE(cnn_use_max_pool),
#         cnn_pool_kernel  = as.integer(cnn_pool_kernel),
#         cnn_pool_stride  = as.integer(cnn_pool_stride),
#         cnn_pool_padding = as.integer(cnn_pool_padding),
#         separable=isTRUE(separable),
#         dilations=as.integer(dilations),
#         use_se=isTRUE(use_se),
#         norm_type=norm_type,
#         #use_max_pool=isTRUE(use_max_pool),
#         pool_type=pool_type,
#         use_global_pool=isTRUE(use_global_pool),
#
#         # Heteroscedastic
#         heteroscedastic = isTRUE(heteroscedastic),
#
#         # FT / SAINT (+ grouping)
#         ft_d_model = as.integer(ft_d_model),
#         ft_heads   = as.integer(ft_heads),
#         ft_layers  = as.integer(ft_layers),
#         ft_ff_mult = as.integer(ft_ff_mult),
#         ft_dropout = as.numeric(ft_dropout),
#         ft_token_dropout = as.numeric(ft_token_dropout),
#         ft_use_cls = isTRUE(ft_use_cls),
#
#         saint_d_model = as.integer(saint_d_model),
#         saint_heads   = as.integer(saint_heads),
#         saint_layers  = as.integer(saint_layers),
#         saint_ff_mult = as.integer(saint_ff_mult),
#         saint_dropout = as.numeric(saint_dropout),
#         saint_token_dropout = as.numeric(saint_token_dropout),
#         saint_use_cls = isTRUE(saint_use_cls),
#
#         use_grouping   = isTRUE(use_grouping),
#         group_trigger  = as.integer(group_trigger),
#         group_method   = as.character(group_method),
#         init_group_size = as.integer(init_group_size),
#         max_tokens      = as.integer(max_tokens),
#         kmeans_batch    = as.integer(kmeans_batch),
#         kmeans_iter     = as.integer(kmeans_iter),
#
#         # TabNet
#         tabnet_steps         = as.integer(tabnet_steps),
#         tabnet_feature_dim   = as.integer(tabnet_feature_dim),
#         tabnet_output_dim    = as.integer(tabnet_output_dim),
#         tabnet_gamma         = as.numeric(tabnet_gamma),
#         tabnet_lambda_sparse = as.numeric(tabnet_lambda_sparse),
#
#         # NODE
#         node_trees = as.integer(node_trees),
#         node_depth = as.integer(node_depth),
#
#         # DeepFM / DCN
#         deepfm_k      = as.integer(deepfm_k),
#         deepfm_hidden = as.integer(deepfm_hidden),
#         dcn_layers    = as.integer(dcn_layers),
#         dcn_hidden    = as.integer(dcn_hidden),
#
#         # NAM
#         nam_hidden     = as.integer(nam_hidden),
#         nam_activation = as.character(nam_activation),
#         nam_add_linear = isTRUE(nam_add_linear),
#         nam_l1         = as.numeric(nam_l1),
#
#         # MoE
#         moe_n_experts     = as.integer(moe_n_experts),
#         moe_expert_hidden = as.integer(moe_expert_hidden),
#         moe_gate_hidden   = as.integer(moe_gate_hidden),
#         moe_temperature   = as.numeric(moe_temperature),
#         moe_sparse_topk   = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#         moe_entropy_reg   = as.numeric(moe_entropy_reg),
#
#         # GP-DKL / RFF
#         gp_use_variational = isTRUE(gp_use_variational),
#         gp_num_inducing    = as.integer(gp_num_inducing),
#         gp_feature_dim     = as.integer(gp_feature_dim),
#         gp_kernel          = as.character(gp_kernel),
#         gp_ard             = isTRUE(gp_ard),
#         gp_lr_mult         = as.numeric(gp_lr_mult),
#         rff_features       = as.integer(rff_features),
#         rff_lengthscale    = as.numeric(rff_lengthscale),
#         rff_deep_hidden    = as.integer(rff_deep_hidden)
#       )
#
#       preds <- mod$predict(py_fit[[1L]], as.matrix(omics_data[tst, , drop = FALSE]),
#                            device = if (is.null(device)) NULL else as.character(device))
#       as.numeric(preds)
#     }, error = function(e) {
#       message("An error occurred: ", e$message)
#       NULL
#     })
#     return(out)
#   }
#
#   # 5) BOOTSTRAP TRAIN/PREDICT  (forward all NEW knobs into statistic)
#   if (!is.null(geno_omic_test_object)) {
#     X_test <- as.matrix(geno_omic_test_object)
#
#     boot_results <- boot::boot(
#       data = data_label_geno,
#       statistic = train_predict_deeplearning,
#       num_hidden_layers = as.integer(num_hidden_layers),
#       neurons_per_layer = if (!is.null(neurons_per_layer)) as.integer(neurons_per_layer) else NULL,
#       learning_rate = as.numeric(learning_rate_dp),
#       epochs = as.integer(epochs),
#       batch_size = as.integer(batch_size),
#       l2_regularizer_dp = as.numeric(l2_regularizer_dp),
#       dropout_rate = as.numeric(dropout_rate),
#       validation_split = as.numeric(validation_split),
#       R = n_bootstrap,
#       test_geno = X_test,
#       deep_learning_model = deep_learning_model,
#       attention_on_final_layer = attention_on_final_layer,
#       attention_across_multiple_layers = attention_across_multiple_layers,
#       batch_normalization = batch_normalization,
#       dense_layers_cnn = as.integer(dense_layers_cnn),
#       kernel_size = as.integer(kernel_size),
#       n_blocks = as.integer(n_blocks),
#       n_neurons_per_block = if (!is.null(n_neurons_per_block)) as.integer(n_neurons_per_block) else NULL,
#       compile_model = compile_model,
#       deterministic = deterministic,
#       random_seed = random_seed,
#       device = device,
#
#       # Optim/AMP/weights/clip
#       optimizer_name = optimizer_name,
#       use_amp = isTRUE(use_amp),
#       max_grad_norm = as.numeric(max_grad_norm),
#       auto_class_weights = isTRUE(auto_class_weights),
#
#       # CNN max-pool
#       cnn_use_max_pool = isTRUE(cnn_use_max_pool),
#       cnn_pool_kernel  = as.integer(cnn_pool_kernel),
#       cnn_pool_stride  = as.integer(cnn_pool_stride),
#       cnn_pool_padding = as.integer(cnn_pool_padding),
#       separable=isTRUE(separable),
#       dilations=as.integer(dilations),
#       use_se=isTRUE(use_se),
#       norm_type=norm_type,
#       #use_max_pool=isTRUE(use_max_pool),
#       pool_type=pool_type,
#       use_global_pool=isTRUE(use_global_pool),
#
#       # Heteroscedastic
#       heteroscedastic  = isTRUE(heteroscedastic),
#
#       # FT / SAINT (+ grouping)
#       ft_d_model       = as.integer(ft_d_model),
#       ft_heads         = as.integer(ft_heads),
#       ft_layers        = as.integer(ft_layers),
#       ft_ff_mult       = as.integer(ft_ff_mult),
#       ft_dropout       = as.numeric(ft_dropout),
#       ft_token_dropout = as.numeric(ft_token_dropout),
#       ft_use_cls       = isTRUE(ft_use_cls),
#
#       saint_d_model       = as.integer(saint_d_model),
#       saint_heads         = as.integer(saint_heads),
#       saint_layers        = as.integer(saint_layers),
#       saint_ff_mult       = as.integer(saint_ff_mult),
#       saint_dropout       = as.numeric(saint_dropout),
#       saint_token_dropout = as.numeric(saint_token_dropout),
#       saint_use_cls       = isTRUE(saint_use_cls),
#
#       use_grouping   = isTRUE(use_grouping),
#       group_trigger  = as.integer(group_trigger),
#       group_method   = as.character(group_method),
#       init_group_size = as.integer(init_group_size),
#       max_tokens      = as.integer(max_tokens),
#       kmeans_batch    = as.integer(kmeans_batch),
#       kmeans_iter     = as.integer(kmeans_iter),
#
#       # TabNet
#       tabnet_steps         = as.integer(tabnet_steps),
#       tabnet_feature_dim   = as.integer(tabnet_feature_dim),
#       tabnet_output_dim    = as.integer(tabnet_output_dim),
#       tabnet_gamma         = as.numeric(tabnet_gamma),
#       tabnet_lambda_sparse = as.numeric(tabnet_lambda_sparse),
#
#       # NODE
#       node_trees          = as.integer(node_trees),
#       node_depth          = as.integer(node_depth),
#
#       # DeepFM / DCN
#       deepfm_k            = as.integer(deepfm_k),
#       deepfm_hidden       = as.integer(deepfm_hidden),
#       dcn_layers          = as.integer(dcn_layers),
#       dcn_hidden          = as.integer(dcn_hidden),
#
#       # NAM
#       nam_hidden     = as.integer(nam_hidden),
#       nam_activation = as.character(nam_activation),
#       nam_add_linear = isTRUE(nam_add_linear),
#       nam_l1         = as.numeric(nam_l1),
#
#       # MoE
#       moe_n_experts     = as.integer(moe_n_experts),
#       moe_expert_hidden = as.integer(moe_expert_hidden),
#       moe_gate_hidden   = as.integer(moe_gate_hidden),
#       moe_temperature   = as.numeric(moe_temperature),
#       moe_sparse_topk   = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#       moe_entropy_reg   = as.numeric(moe_entropy_reg),
#
#       # GP-DKL / RFF
#       gp_use_variational = isTRUE(gp_use_variational),
#       gp_num_inducing    = as.integer(gp_num_inducing),
#       gp_feature_dim     = as.integer(gp_feature_dim),
#       gp_kernel          = as.character(gp_kernel),
#       gp_ard             = isTRUE(gp_ard),
#       gp_lr_mult         = as.numeric(gp_lr_mult),
#       rff_features       = as.integer(rff_features),
#       rff_lengthscale    = as.numeric(rff_lengthscale),
#       rff_deep_hidden    = as.integer(rff_deep_hidden)
#     )
#
#     revert_scaling_ml <- function(scaled_values, scaler_mean, scaler_sd) {
#       scaled_values * scaler_sd + scaler_mean
#     }
#     boot_results$t <- apply(boot_results$t, 2, function(col_vec) {
#       revert_scaling_ml(col_vec, y_scaler$mean, y_scaler$std)
#     })
#
#     pred_variances <- apply(boot_results$t, 2, var)
#     pred_SE <- apply(boot_results$t, 2, sd)
#     AI_pred <- apply(boot_results$t, 2, mean)
#     genetic_var <- var(AI_pred)
#     AI_pred_reverted <- AI_pred
#
#     result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
#       boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
#     )
#     result_rel <- reliability_thresholds(
#       prediction_error_var = pred_variances,
#       genetic_var = genetic_var,
#       high_reliability_thres = high_reliability_thres,
#       low_reliability_thres = low_reliability_thres
#     )
#
#     train_test_label <- if (!is.null(geno_omic_test_object)) {
#       ifelse(rownames(geno_omic_test_object) %in% test_label, "Test", "Train")
#     } else {
#       rep("Train", nrow(geno_omic_object))
#     }
#
#     AI_preds <- data.frame(
#       name = GID,
#       Predicted_value = AI_pred_reverted,
#       Train_Test_Label = train_test_label,
#       Standard_error = pred_SE,
#       PEV = pred_variances,
#       lower_bound = result_rel_MPIW$lower_bound,
#       upper_bound = result_rel_MPIW$upper_bound,
#       Uncertainty = result_rel_MPIW$Uncertainty,
#       Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#       Reliability = result_rel$reliability,
#       Reliability_remarks = result_rel$remarks,
#       Reliability_percentage = result_rel$reliability_percentage,
#       stringsAsFactors = FALSE
#     )
#     names(AI_preds)[1] <- c(gen_name)
#
#     diagnostic_plots <- diagnostic_plot_true_prediction(
#       boot_results = boot_results,
#       GID_names = GID,
#       CI_width_thresholds = CI_width_thresholds,
#       predictions = AI_pred_reverted,
#       standard_errors = pred_SE,
#       prediction_error_var = pred_variances,
#       genetic_var = genetic_var,
#       confidence_level = 0.95,
#       model_for_CI_cal = "ML",
#       high_reliability_thres = high_reliability_thres,
#       low_reliability_thres = low_reliability_thres,
#       system_database = system_database
#     )
#
#   } else {
#     boot_results <- boot::boot(
#       data = data_label_geno,
#       statistic = train_predict_deeplearning,
#       num_hidden_layers = as.integer(num_hidden_layers),
#       neurons_per_layer = if (!is.null(neurons_per_layer)) as.integer(neurons_per_layer) else NULL,
#       learning_rate = as.numeric(learning_rate_dp),
#       epochs = as.integer(epochs),
#       batch_size = as.integer(batch_size),
#       l2_regularizer_dp = as.numeric(l2_regularizer_dp),
#       dropout_rate = as.numeric(dropout_rate),
#       R = n_bootstrap,
#       test_geno = NULL,
#       deep_learning_model = deep_learning_model,
#       attention_on_final_layer = attention_on_final_layer,
#       attention_across_multiple_layers = attention_across_multiple_layers,
#       batch_normalization = batch_normalization,
#       dense_layers_cnn = as.integer(dense_layers_cnn),
#       kernel_size = as.integer(kernel_size),
#       n_blocks = as.integer(n_blocks),
#       n_neurons_per_block = if (!is.null(n_neurons_per_block)) as.integer(n_neurons_per_block) else NULL,
#       compile_model = compile_model,
#       deterministic = deterministic,
#       random_seed = random_seed,
#       device = device,
#
#       # Optim/AMP/weights/clip
#       optimizer_name = optimizer_name,
#       use_amp = isTRUE(use_amp),
#       max_grad_norm = as.numeric(max_grad_norm),
#       auto_class_weights = isTRUE(auto_class_weights),
#
#       # CNN max-pool
#       cnn_use_max_pool = isTRUE(cnn_use_max_pool),
#       cnn_pool_kernel  = as.integer(cnn_pool_kernel),
#       cnn_pool_stride  = as.integer(cnn_pool_stride),
#       cnn_pool_padding = as.integer(cnn_pool_padding),
#       separable=isTRUE(separable),
#       dilations=as.integer(dilations),
#       use_se=isTRUE(use_se),
#       norm_type=norm_type,
#       #use_max_pool=isTRUE(use_max_pool),
#       pool_type=pool_type,
#       use_global_pool=isTRUE(use_global_pool),
#
#       # Heteroscedastic
#       heteroscedastic  = isTRUE(heteroscedastic),
#
#       # FT / SAINT (+ grouping)
#       ft_d_model       = as.integer(ft_d_model),
#       ft_heads         = as.integer(ft_heads),
#       ft_layers        = as.integer(ft_layers),
#       ft_ff_mult       = as.integer(ft_ff_mult),
#       ft_dropout       = as.numeric(ft_dropout),
#       ft_token_dropout = as.numeric(ft_token_dropout),
#       ft_use_cls       = isTRUE(ft_use_cls),
#
#       saint_d_model       = as.integer(saint_d_model),
#       saint_heads         = as.integer(saint_heads),
#       saint_layers        = as.integer(saint_layers),
#       saint_ff_mult       = as.integer(saint_ff_mult),
#       saint_dropout       = as.numeric(saint_dropout),
#       saint_token_dropout = as.numeric(saint_token_dropout),
#       saint_use_cls       = isTRUE(saint_use_cls),
#
#       use_grouping   = isTRUE(use_grouping),
#       group_trigger  = as.integer(group_trigger),
#       group_method   = as.character(group_method),
#       init_group_size = as.integer(init_group_size),
#       max_tokens      = as.integer(max_tokens),
#       kmeans_batch    = as.integer(kmeans_batch),
#       kmeans_iter     = as.integer(kmeans_iter),
#
#       # TabNet
#       tabnet_steps         = as.integer(tabnet_steps),
#       tabnet_feature_dim   = as.integer(tabnet_feature_dim),
#       tabnet_output_dim    = as.integer(tabnet_output_dim),
#       tabnet_gamma         = as.numeric(tabnet_gamma),
#       tabnet_lambda_sparse = as.numeric(tabnet_lambda_sparse),
#
#       # NODE
#       node_trees          = as.integer(node_trees),
#       node_depth          = as.integer(node_depth),
#
#       # DeepFM / DCN
#       deepfm_k            = as.integer(deepfm_k),
#       deepfm_hidden       = as.integer(deepfm_hidden),
#       dcn_layers          = as.integer(dcn_layers),
#       dcn_hidden          = as.integer(dcn_hidden),
#
#       # NAM
#       nam_hidden     = as.integer(nam_hidden),
#       nam_activation = as.character(nam_activation),
#       nam_add_linear = isTRUE(nam_add_linear),
#       nam_l1         = as.numeric(nam_l1),
#
#       # MoE
#       moe_n_experts     = as.integer(moe_n_experts),
#       moe_expert_hidden = as.integer(moe_expert_hidden),
#       moe_gate_hidden   = as.integer(moe_gate_hidden),
#       moe_temperature   = as.numeric(moe_temperature),
#       moe_sparse_topk   = if (is.null(moe_sparse_topk)) NULL else as.integer(moe_sparse_topk),
#       moe_entropy_reg   = as.numeric(moe_entropy_reg),
#
#       # GP-DKL / RFF
#       gp_use_variational = isTRUE(gp_use_variational),
#       gp_num_inducing    = as.integer(gp_num_inducing),
#       gp_feature_dim     = as.integer(gp_feature_dim),
#       gp_kernel          = as.character(gp_kernel),
#       gp_ard             = isTRUE(gp_ard),
#       gp_lr_mult         = as.numeric(gp_lr_mult),
#       rff_features       = as.integer(rff_features),
#       rff_lengthscale    = as.numeric(rff_lengthscale),
#       rff_deep_hidden    = as.integer(rff_deep_hidden)
#     )
#
#     pred_variances <- apply(boot_results$t, 2, var)
#     pred_SE <- apply(boot_results$t, 2, sd)
#     AI_pred <- apply(boot_results$t, 2, mean)
#     genetic_var <- var(AI_pred)
#     AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)
#
#     result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
#       boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
#     )
#     result_rel <- reliability_thresholds(
#       prediction_error_var = pred_variances,
#       genetic_var = genetic_var,
#       high_reliability_thres = high_reliability_thres,
#       low_reliability_thres = low_reliability_thres
#     )
#
#     train_test_label <- if (!is.null(geno_omic_test_object)) {
#       ifelse(rownames(geno_omic_test_object) %in% test_label, "Test", "Train")
#     } else {
#       rep("Train", nrow(geno_omic_object))
#     }
#
#     AI_preds <- data.frame(
#       name = GID,
#       Predicted_value = AI_pred_reverted,
#       Standard_error = pred_SE,
#       Train_Test_Label = train_test_label,
#       PEV = pred_variances,
#       lower_bound = result_rel_MPIW$lower_bound,
#       upper_bound = result_rel_MPIW$upper_bound,
#       Uncertainty = result_rel_MPIW$Uncertainty,
#       Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#       Reliability = result_rel$reliability,
#       Reliability_remarks = result_rel$remarks,
#       Reliability_percentage = result_rel$reliability_percentage,
#       stringsAsFactors = FALSE
#     )
#     names(AI_preds)[1] <- c(gen_name)
#
#     diagnostic_plots <- diagnostic_plot_true_prediction(
#       boot_results = boot_results,
#       GID_names = GID,
#       CI_width_thresholds = CI_width_thresholds,
#       predictions = AI_pred_reverted,
#       standard_errors = pred_SE,
#       prediction_error_var = pred_variances,
#       genetic_var = genetic_var,
#       confidence_level = 0.95,
#       model_for_CI_cal = "ML",
#       high_reliability_thres = high_reliability_thres,
#       low_reliability_thres = low_reliability_thres,
#       system_database = system_database
#     )
#   }
#
#   # Model parameter summary
#   model_para <- data.frame(
#     stat = c("num_hidden_layers", "learning_rate", "neurons_per_layer"),
#     summary = I(list(num_hidden_layers, learning_rate_dp, neurons_per_layer)),
#     stringsAsFactors = FALSE
#   )
#   if (!is.null(neurons_per_layer) && length(neurons_per_layer) > 1) {
#     expanded_stat <- unlist(lapply(1:nrow(model_para), function(i) {
#       rep(model_para$stat[i], length(model_para$summary[[i]]))
#     }))
#     expanded_summary <- unlist(model_para$summary)
#     model_para <- data.frame(stat = expanded_stat, summary = expanded_summary, stringsAsFactors = FALSE)
#   }
#   colnames(model_para)[1:2] <- c("stat", "summary")
#
#   list(
#     model_parameters = model_para,
#     predicted_values = AI_preds,
#     diagnostic_plots = diagnostic_plots
#   )
# }







#'
#'
#' generate_dynamic_layers <- function(input_size, num_hidden_layers, scaling_factor = 0.5, max_neurons = 1000) {
#'
#'   msg <- "\n==================================================\n"
#'
#'   if (length(num_hidden_layers) > 1 || length(num_hidden_layers) == 0) {
#'     stop(paste(msg, "num_hidden_layers should be a vector of length 1", call. = FALSE))
#'   }
#'
#'   layers <- numeric(num_hidden_layers)
#'
#'   # Cap the input size to prevent large layers
#'   capped_input_size <- ifelse(input_size >= 1000, max_neurons, input_size)
#'
#'   # The first hidden layer could start as a fraction of the (capped) input size
#'   layers[1] <- min(floor(capped_input_size * scaling_factor), max_neurons)
#'
#'   # Subsequent layers reduce in size, following the scaling factor and capped by max_neurons
#'   for (i in 2:num_hidden_layers) {
#'     layers[i] <- min(floor(layers[i - 1] * scaling_factor), max_neurons)
#'     if (layers[i] < 1) break  # Stop adding layers if neurons fall below 1
#'   }
#'
#'   return(layers[layers > 0])  # Return only valid layers
#' }
#'
#'
#' # Function to train and predict using dp for bootstrapping
#' train_predict_deeplearning <- function(data_label_geno, indices, test_geno,
#'                                        num_hidden_layers = NULL,
#'                                        neurons_per_layer = NULL,
#'                                        learning_rate = 0.001,
#'                                        epochs = 10,
#'                                        batch_size = 32,
#'                                        validation_split = 0.2,
#'                                        l2_regularizer_dp = 0.001,
#'                                        dropout_rate = 0.5,
#'                                        n_blocks = 2,
#'                                        dense_layers_cnn = c(128, 64),
#'                                        kernel_size = 3,
#'                                        n_neurons_per_block = NULL,
#'                                        deep_learning_model = "mlp_with_attention",
#'                                        attention_on_final_layer = TRUE,
#'                                        attention_across_multiple_layers = FALSE,
#'                                        batch_normalization = TRUE
#'                                        ) {
#'
#'   # Check and configure TensorFlow GPU memory growth before training
#'   if (reticulate::py_module_available("tensorflow")) {
#'     tf <- reticulate::import("tensorflow", delay_load = TRUE)
#'     physical_devices <- tf$config$list_physical_devices("GPU")
#'
#'     if (length(physical_devices) > 0) {
#'       for (dev in physical_devices) {
#'         tf$config$experimental$set_memory_growth(dev, TRUE)
#'       }
#'     }
#'   }
#'
#'   train_data <- data_label_geno[, -1]
#'   y_train <- data_label_geno[, 1]
#'   # Subset the data
#'   model <- deep_learning_model_utilityy(X_train = train_data,
#'                                         y_train = y_train,
#'                                         num_hidden_layers = num_hidden_layers,
#'                                         neurons_per_layer = neurons_per_layer,
#'                                         learning_rate = learning_rate,
#'                                         epochs = epochs,
#'                                         batch_size = batch_size,
#'                                         dropout_rate = dropout_rate,
#'                                         l2_regularizer_dp = l2_regularizer_dp,
#'                                         validation_split = validation_split,
#'                                         n_blocks = n_blocks,
#'                                         kernel_size = kernel_size,
#'                                         dense_layers_cnn = dense_layers_cnn,
#'                                         n_neurons_per_block = n_neurons_per_block,
#'                                         deep_learning_model = deep_learning_model,
#'                                         attention_on_final_layer = attention_on_final_layer,
#'                                         attention_across_multiple_layers = attention_across_multiple_layers,
#'                                         batch_normalization = batch_normalization,
#'                                         para_tunning = FALSE)
#'   # Predict on the original data
#'   if(!is.null(test_geno)){
#'     pred <-  model$predict(test_geno)
#'   } else{
#'     pred <-  model$predict(train_data)
#'   }
#'
#'   if (reticulate::py_module_available("tensorflow")) {
#'     keras::k_clear_session()
#'   }
#'
#'   rm(model)
#'   return(pred)
#' }
#'
#' get_best_model <- function(results, hyperparam_combinations) {
#'   best_metric_value <- Inf  # Start with the highest possible error for regression or lowest possible accuracy for classification
#'   best_model <- NULL
#'   best_hyperparameters <- NULL
#'
#'   for (i in seq_along(results)) {
#'     key <- names(results)[i]
#'     history <- results[[key]]$history
#'
#'     # Determine which metric to use based on the available keys in history
#'     if ("val_accuracy" %in% names(history$history)) {
#'       metric_name <- "val_accuracy"
#'       comparison_function <- max
#'       best_metric_value <- 0  # Initialize to the lowest possible accuracy
#'     } else if ("val_mean_absolute_error" %in% names(history$history)) {
#'       metric_name <- "val_mean_absolute_error"
#'       comparison_function <- min
#'       best_metric_value <- Inf  # Initialize to the highest possible error
#'     } else {
#'       warning("No recognized metric found in the model history.")
#'       next
#'     }
#'
#'     # Extract the metric value
#'     current_metric_value <- comparison_function(history$history[[metric_name]])
#'
#'     # Compare and store the best model and hyperparameters
#'     if ((metric_name == "val_accuracy" && current_metric_value > best_metric_value) ||
#'         (metric_name == "val_mean_absolute_error" && current_metric_value < best_metric_value)) {
#'       best_metric_value <- current_metric_value
#'       best_model <- results[[key]]$model
#'       best_hyperparameters <- hyperparam_combinations[i, ]  # Extract the corresponding hyperparameters
#'     }
#'   }
#'
#'   return(list(best_model = best_model, best_hyperparameters = best_hyperparameters))
#' }
#'
#'
#' #' Title
#' #'
#' #' @param num_hidden_layers
#' #' @param neurons_per_layer
#' #' @param learning_rate
#' #' @param epochs
#' #' @param param_grid
#' #' @param validation_split
#' #' @param early_stop
#' #' @param batch_size
#' #' @param pheno_object
#' #' @param geno_omic_object
#' #' @param geno_omic_test_object
#' #' @param response
#' #' @param gen_name
#' #' @param message
#' #' @param scale
#' #' @param para_tunning
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' deep_learning_model <- function(pheno_object=NULL,
#'                                 y = NULL,
#'                                 omics_data = NULL,
#'                                 crossval = FALSE,
#'                                 tst = NULL,
#'                                 geno_omic_object = NULL,
#'                                 geno_omic_test_object = NULL,
#'                                 response=NULL,
#'                                 gen_name=NULL,
#'                                 message = TRUE,
#'                                 scaling = TRUE,
#'                                 centering = FALSE,
#'                                 omic_count = NULL,
#'                                 num_hidden_layers = 1,
#'                                 neurons_per_layer = NULL,
#'                                 learning_rate_dp = 0.001,
#'                                 epochs = 10,
#'                                 batch_size = 32 ,
#'                                 para_tunning = FALSE,
#'                                 param_grid = NULL,
#'                                 validation_split = 0.2,
#'                                 early_stop = TRUE,
#'                                 l2_regularizer_dp = 0.001,
#'                                 dropout_rate = 0.5,
#'                                 deep_learning_model = "mlp_with_attention",
#'                                 n_blocks = 2,
#'                                 dense_layers_cnn = c(128, 64),
#'                                 kernel_size = 3,
#'                                 n_neurons_per_block = NULL,
#'                                 attention_on_final_layer = TRUE,
#'                                 attention_across_multiple_layers = FALSE,
#'                                 batch_normalization = TRUE,
#'                                 CI_width_thresholds = c(0.33, 0.66),
#'                                 high_reliability_thres = 0.9,
#'                                 low_reliability_thres = 0.5,
#'                                 n_components = 20,
#'                                 threshold = 100,
#'                                 target = "test_set",
#'                                 iqr_multiplier = 1.5,
#'                                 interval_width_high_threshold = NULL,
#'                                 interval_width_low_threshold = NULL,
#'                                 interval_width_moderate_threshold = NULL,
#'                                 n_bootstrap = 100,
#'                                 system_database = FALSE,
#'                                 ...) {
#' #browser()
#'
#'   # Check and configure TensorFlow GPU memory growth before training
#'   if (reticulate::py_module_available("tensorflow")) {
#'     tf <- reticulate::import("tensorflow", delay_load = TRUE)
#'     physical_devices <- tf$config$list_physical_devices("GPU")
#'
#'     if (length(physical_devices) > 0) {
#'       for (dev in physical_devices) {
#'         tf$config$experimental$set_memory_growth(dev, TRUE)
#'       }
#'     }
#'   }
#'
#'   msg <- "\n==================================================\n"
#'   if(is.null(dense_layers_cnn)) dense_layers_cnn <-  64
#'   if (is.null(geno_omic_object) && is.null(pheno_object) && isFALSE(crossval)) {
#'     stop(paste(msg, "provide matrix of the predictors and the data.frame of the Y variable."), call. = FALSE)
#'   }
#'
#'   if(!is.null(geno_omic_object)){
#'     GID <- rownames(geno_omic_object)
#'     scaler <- caret::preProcess(geno_omic_object, method = c("center", "scale"))
#'     geno_omic_object <- stats::predict(scaler, geno_omic_object)
#'
#'     cols_with_na <- which(colSums(is.na(geno_omic_object)) > 0)
#'     if(length(cols_with_na)!=0){
#'       geno_omic_object <- geno_omic_object[, -cols_with_na]
#'       if(!is.null(geno_omic_test_object)){
#'         geno_omic_test_object <- geno_omic_test_object[, -cols_with_na]
#'       }
#'     }
#'   }
#'
#'
#'
#'   if(!is.null(geno_omic_test_object)){
#'
#'     test_label <- rownames(geno_omic_test_object)
#'     geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
#'     GID <- rownames(geno_omic_test_object)
#'
#'     geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
#'
#'   }
#'
#' #### This is for cross-validation. Not typical omic data
#'   if(!is.null(omics_data) && isTRUE(crossval)){
#'     scaler <- caret::preProcess(omics_data, method = c("center", "scale"))
#'     omics_data <- stats::predict(scaler, omics_data)
#'     omics_data <- stats::predict(scaler, omics_data)
#'     #####
#'     cols_with_na <- which(colSums(is.na(omics_data)) > 0)
#'     if(length(cols_with_na)!=0){
#'       omics_data <- omics_data[, -cols_with_na]
#'
#'     }
#'
#'   }
#'
#'   if(isTRUE(crossval)){
#'     para_tunning <-  FALSE
#'     param_grid <-  NULL
#'   }
#'
#'   #"ResNet",
#'   if(any(deep_learning_model%in%c("mlp_with_attention", "mlp", "cnn"))){
#'     #if(is.null(neurons_per_layer)&& is.null(num_hidden_layers)){
#'     if(is.null(num_hidden_layers)){
#'       stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#'     }
#'   }
#'
#'   if(deep_learning_model=="ResNet"){
#'     #if(is.null(n_neurons_per_block) || is.null(n_blocks)){
#'     if(is.null(n_blocks)){
#'       stop(paste(msg, "n_blocks can't be NULL."), call. = FALSE)
#'     }
#'   }
#'
#'   if(is.null(neurons_per_layer)&& !is.null(num_hidden_layers)){
#'     if(!is.null(geno_omic_object) && isFALSE(crossval)){
#'     neurons_per_layer <- generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
#'     } else if(!is.null(omics_data) && isTRUE(crossval)){
#'       neurons_per_layer <- generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
#'
#'     } else {
#'       stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#' }
#'   }
#'
#'
#'
#'   if(deep_learning_model == "ResNet"){
#'     if(is.null(n_blocks)) n_blocks <- 1
#'     #if(is.null(dense_layers_cnn)) dense_layers_cnn <-  64
#'       if(is.null(n_neurons_per_block)&& !is.null(n_blocks)){
#'         if(!is.null(geno_omic_object) && isFALSE(crossval)){
#'         n_neurons_per_block <-  generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = n_blocks, scaling_factor = 0.5)
#'
#'         } else if(!is.null(omics_data) && isTRUE(crossval)){
#'           n_neurons_per_block <-  generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = n_blocks, scaling_factor = 0.5)
#'         } else {
#'           stop(paste(msg, "n_neurons_per_block and n_blocks can't be NULL."), call. = FALSE)
#'         }
#'       }
#'
#'     if(n_blocks!=length(n_neurons_per_block)){
#'       stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
#'     }
#'
#'     neurons_per_layer <- n_neurons_per_block
#'     num_hidden_layers <- n_blocks
#'   } else {
#'     if(is.null(num_hidden_layers)) num_hidden_layers <- 1
#'     if(is.null(neurons_per_layer)&& !is.null(num_hidden_layers)){
#'       if(!is.null(geno_omic_object) && isFALSE(crossval)){
#'         neurons_per_layer <-  generate_dynamic_layers(input_size = ncol(geno_omic_object), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
#'
#'       } else if(!is.null(omics_data) && isTRUE(crossval)){
#'         neurons_per_layer <-  generate_dynamic_layers(input_size = ncol(omics_data), num_hidden_layers = num_hidden_layers, scaling_factor = 0.5)
#'       } else {
#'         stop(paste(msg, "neurons_per_layer and num_hidden_layers can't be NULL."), call. = FALSE)
#'       }
#'     }
#'   }
#'
#'   if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers) else num_hidden_layers <- as.integer(1)
#'   if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.numeric(neurons_per_layer) else neurons_per_layer <- as.numeric(64)
#'   if(!is.null(batch_size))   batch_size <- as.integer(batch_size) else  batch_size <- as.integer(30)
#'   if(!is.null(epochs)) epochs <- as.integer(epochs) else epochs <- as.integer(10)
#'   if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate) else dropout_rate <- as.numeric(0.5)
#'   if(!is.null(learning_rate_dp)) learning_rate <- as.numeric(learning_rate_dp) else learning_rate <- as.numeric(0.01)
#'   if(!is.null(l2_regularizer_dp)) l2_regularizer_dp <- as.integer(l2_regularizer_dp)
#'
#'   # if(!is.null(geno_omic_test_object)){
#'   #   GID <- rownames(geno_omic_test_object)
#'   # } else {
#'   #   if(!is.null(geno_omic_object)){
#'   #     GID <- rownames(geno_omic_object)
#'   #   }
#'   # }
#'
#'    np <- reticulate::import("numpy")
#'
#'    if(isFALSE(crossval)){
#'   y_train <- as.numeric(pheno_object[, response])
#'   # y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
#'   #
#'   # # Predict on the training data and get the scaled values
#'   # y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
#'
#'    } else{
#'      tryCatch({
#'        if(isTRUE(crossval)){
#'          y_train <- y
#'          y_train = y_train[-tst]
#'          y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
#'
#'          # Predict on the training data and get the scaled values
#'          y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
#'
#'          model_dp <- deep_learning_model_utilityy(X_train = omics_data[-tst, ],
#'                                                #y_train = y_train[-tst],
#'                                                y_train = y_train,
#'                                                num_hidden_layers = num_hidden_layers,
#'                                                neurons_per_layer = neurons_per_layer,
#'                                                learning_rate = learning_rate,
#'                                                epochs = epochs,
#'                                                batch_size = batch_size,
#'                                                dropout_rate = dropout_rate,
#'                                                l2_regularizer_dp = l2_regularizer_dp,
#'                                                validation_split = validation_split,
#'                                                n_blocks = num_hidden_layers,
#'                                                n_neurons_per_block = neurons_per_layer,
#'                                                kernel_size = kernel_size,
#'                                                dense_layers_cnn = dense_layers_cnn,
#'                                                deep_learning_model = deep_learning_model,
#'                                                attention_on_final_layer = attention_on_final_layer,
#'                                                attention_across_multiple_layers = attention_across_multiple_layers,
#'                                                batch_normalization = batch_normalization,
#'                                                para_tunning = FALSE)
#'
#'          preds <- model_dp$predict(omics_data[tst, ])
#'          preds <- as.data.frame(preds)
#'          ###
#'          if (reticulate::py_module_available("tensorflow")) {
#'            keras::k_clear_session()
#'          }
#'
#'          rm(model_dp)
#'          return(preds[, 1])
#'        }
#'      }, error = function(e) {
#'        message("An error occurred: ", e$message)
#'        return(NULL) # or handle the error as needed
#'      })
#'    }
#'
#'   # Scale the training labels
#'   y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
#'
#'   # Predict on the training data and get the scaled values
#'   y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
#'
#'
#'   data_label_geno <- cbind(y_train_scaled,geno_omic_object)
#'
#'   if (isTRUE(para_tunning) && is.null(param_grid) && isFALSE(crossval)) {
#'     stop(paste(msg, "param_grid must be provided when tuning is enabled."), call. = FALSE)
#'   }
#'
#'
#'   if(is.null(neurons_per_layer)) neurons_per_layer <- list(ncol(geno_omic_object)/2)
#'
#'   # Create and compile the model
#'
#'
#'   if (isTRUE(para_tunning)) {
#'     # Perform grid search for hyperparameter tuning
#'     # Generate all combinations of hyperparameters
#'     hyperparam_combinations <- expand.grid(param_grid)
#'
#'     # Function to check if the number of hidden layers matches the length of neurons_per_layer
#'     check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
#'       neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
#'       return(length(neuron_list) == num_hidden_layers)
#'     }
#'
#'     # Loop through the data frame and remove rows with mismatched hidden layers and neurons
#'     valid_rows <- apply(hyperparam_combinations, 1, function(row) {
#'       num_hidden_layers <- as.integer(row["num_hidden_layers"])
#'       neurons_per_layer <- row["neurons_per_layer"]
#'       check_hidden_layers(num_hidden_layers, neurons_per_layer)
#'     })
#'
#'     # Filter the data frame to keep only valid rows
#'     hyperparam_combinations <- hyperparam_combinations[valid_rows, ]
#'
#'     grid_results <- grid_search_deep_learning(X_train = geno_omic_object,
#'                                              y_train = y_train_scaled,
#'                                              param_grid = param_grid,
#'                                              epochs = epochs,
#'                                              batch_size = batch_size,
#'                                              validation_split = validation_split,
#'                                              early_stop = early_stop,
#'                                              n_blocks = n_blocks,
#'                                              n_neurons_per_block = n_neurons_per_block,
#'                                              dense_layers_cnn = dense_layers_cnn,
#'                                              kernel_size = kernel_size,
#'                                              learning_rate = learning_rate_dp,
#'                                              dropout_rate = dropout_rate,
#'                                              l2_regularizer_dp = l2_regularizer_dp,
#'                                              para_tunning = para_tunning)
#'
#'
#'
#'     best_model <- get_best_model(grid_results, hyperparam_combinations)
#'     best_hyperparameter <- best_model$best_hyperparameters
#'
#'     # Initialize variables only if they exist in the hyperparams
#'     if ("num_hidden_layers" %in% names(best_hyperparameter )) {
#'       num_hidden_layers <- as.numeric(best_hyperparameter $num_hidden_layers)
#'     }
#'
#'     if ("neurons_per_layer" %in% names(best_hyperparameter )) {
#'       if (inherits(best_hyperparameter $neurons_per_layer, "list")) {
#'         neurons_per_layer <- unlist(best_hyperparameter $neurons_per_layer)
#'       } else {
#'         neurons_per_layer <- best_hyperparameter$neurons_per_layer
#'       }
#'     }
#'
#'     if ("learning_rate" %in% names(best_hyperparameter )) {
#'       learning_rate <- best_hyperparameter $learning_rate_dp
#'     }
#'
#'     if ("dropout_rate" %in% names(best_hyperparameter )) {
#'       dropout_rate <- best_hyperparameter $dropout_rate
#'     }
#'
#'     if ("epochs" %in% names(best_hyperparameter)) {
#'       epochs <- best_hyperparameter$epochs
#'     }
#'
#'     if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers) else num_hidden_layers <- as.integer(1)
#'     if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.numeric(neurons_per_layer) else neurons_per_layer <- as.numeric(64)
#'     if(!is.null(batch_size))   batch_size <- as.integer(batch_size) else  batch_size <- as.integer(30)
#'     if(!is.null(epochs)) epochs <- as.integer(epochs) else epochs <- as.integer(10)
#'     if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate) else dropout_rate <- as.numeric(0.5)
#'     if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate) else learning_rate <- as.numeric(0.01)
#'
#' }
#'
#'     if (!is.null(geno_omic_test_object)) {
#'       X_test <- np$array(as.matrix(geno_omic_test_object))
#'       #predictions <- best_model$best_model$predict(X_test)
#'
#'       boot_results <- boot::boot(
#'                                 data = data_label_geno,
#'                                 statistic = train_predict_deeplearning,
#'                                 num_hidden_layers = num_hidden_layers,
#'                                 neurons_per_layer = neurons_per_layer,
#'                                 learning_rate = learning_rate,
#'                                 epochs = epochs,
#'                                 batch_size = batch_size,
#'                                 l2_regularizer_dp = l2_regularizer_dp,
#'                                 dropout_rate = dropout_rate,
#'                                 validation_split = validation_split,
#'                                 R = n_bootstrap,  # Number of bootstrap samples
#'                                 #sim = "ordinary",
#'                                 test_geno =X_test
#'                               )
#'
#'       # Extract bootstrap predictions
#'       revert_scaling_ml <- function(scaled_values, scaler_mean, scaler_sd) {
#'         scaled_values * scaler_sd + scaler_mean
#'       }
#'
#'       boot_results$t <- apply(
#'         boot_results$t,
#'         2,
#'         function(col_vec) revert_scaling_ml(col_vec, y_scaler$mean, y_scaler$std)
#'       )
#'       pred_variances <- apply(boot_results$t, 2, var)
#'       pred_SE <- apply(boot_results$t, 2, sd)
#'       AI_pred <- apply(boot_results$t, 2, mean)
#'       genetic_var <- var(AI_pred)
#'       #AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)
#'       AI_pred_reverted <-  AI_pred
#'
#'
#'       result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(boot_results = boot_results,
#'                                                              CI_width_thresholds = CI_width_thresholds)
#'       result_rel <-  reliability_thresholds(prediction_error_var = pred_variances,
#'                                             genetic_var = genetic_var,
#'                                             high_reliability_thres = high_reliability_thres,
#'                                             low_reliability_thres = low_reliability_thres
#'       )
#'
#'       composite_reliability <- composite_reliability_tst(geno_trn = geno_omic_object,
#'                                                          geno_tst = geno_omic_test_object,
#'                                                          geno_tst_trn = if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) geno_omic_object else NULL,
#'                                                          names_tst = NULL,
#'                                                          names_trn = NULL,
#'                                                          n_components = n_components,
#'                                                          threshold = threshold,
#'                                                          target = target,
#'                                                          interval_width = result_rel_MPIW$Uncertainty,
#'                                                          CI_width_thresholds = CI_width_thresholds,
#'                                                          interval_width_high_threshold = interval_width_high_threshold,
#'                                                          interval_width_low_threshold = interval_width_low_threshold,
#'                                                          apply_pca = TRUE)
#'
#'       if(!is.null(geno_omic_test_object)){
#'
#'         train_test_label <- ifelse(rownames(geno_omic_test_object)%in%test_label, "Test", "Train")
#'       } else {
#'         train_test_label <- rep("Train", nrow(geno_omic_object))
#'       }
#'
#'       AI_preds <- data.frame(name = GID,
#'                              Predicted_value = AI_pred_reverted,
#'                              Train_Test_Label = train_test_label,
#'                              Standard_error = pred_SE,
#'                              PEV = pred_variances,
#'                              lower_bound = result_rel_MPIW$lower_bound,
#'                              upper_bound = result_rel_MPIW$upper_bound,
#'                              Uncertainty = result_rel_MPIW$Uncertainty,
#'                              Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#'                              Reliability = result_rel$reliability,
#'                              Reliability_remarks = result_rel$remarks,
#'                              Reliability_percentage = result_rel$reliability_percentage,
#'                              #Composite_reliability = composite_reliability$trustworthiness,
#'                              #Composite_reliability_percentage = composite_reliability$reliability_percentage,
#'                              stringsAsFactors = FALSE)
#'
#'
#'       names(AI_preds)[1] <-  c(gen_name)
#'
#'       diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = boot_results,
#'                                                           GID_names = GID,
#'                                                           CI_width_thresholds = CI_width_thresholds,
#'                                                           predictions = AI_pred_reverted,
#'                                                           standard_errors = pred_SE,
#'                                                           prediction_error_var = pred_variances,
#'                                                           genetic_var = genetic_var,
#'                                                           confidence_level = 0.95,
#'                                                           model_for_CI_cal = "ML",
#'                                                           #composite_reliability_score = composite_reliability$reliability_score,
#'                                                           #composite_reliability = composite_reliability$trustworthiness,
#'                                                           #composite_reliability_percentage = composite_reliability$reliability_percentage,
#'                                                           #threshold = NULL,
#'                                                           high_reliability_thres = high_reliability_thres,
#'                                                           low_reliability_thres = low_reliability_thres,
#'                                                           system_database = system_database)
#'
#'
#'      # return(list(predicted_values = predictions, trained_model = best_model))
#'     } else {
#'
#'       boot_results <- boot::boot(
#'         data = data_label_geno,
#'         statistic = train_predict_deeplearning,
#'         num_hidden_layers = num_hidden_layers,
#'         neurons_per_layer = neurons_per_layer,
#'         learning_rate = learning_rate,
#'         epochs = epochs,
#'         batch_size = batch_size,
#'         l2_regularizer_dp = l2_regularizer_dp,
#'         dropout_rate = dropout_rate,
#'         R = n_bootstrap,  # Number of bootstrap samples
#'         #sim = "ordinary",
#'         test_geno =NULL
#'       )
#'
#'       # Extract bootstrap predictions
#'       pred_variances <- apply(boot_results$t, 2, var)
#'       pred_SE <- apply(boot_results$t, 2, sd)
#'       AI_pred <- apply(boot_results$t, 2, mean)
#'       genetic_var <- var(AI_pred)
#'       AI_pred_reverted <- revert_scaling(AI_pred, y_scaler)
#'
#'
#'       result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(boot_results = boot_results,
#'                                                              CI_width_thresholds = CI_width_thresholds)
#'       result_rel <-  reliability_thresholds(prediction_error_var = pred_variances,
#'                                             genetic_var = genetic_var,
#'                                             high_reliability_thres = high_reliability_thres,
#'                                             low_reliability_thres = low_reliability_thres
#'       )
#'
#'       composite_reliability <- composite_reliability_tst(geno_trn = geno_omic_object,
#'                                                          geno_tst = geno_omic_test_object,
#'                                                          geno_tst_trn = if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) geno_omic_object else NULL,
#'                                                          names_tst = NULL,
#'                                                          names_trn = NULL,
#'                                                          n_components = n_components,
#'                                                          threshold = threshold,
#'                                                          target = target,
#'                                                          interval_width = result_rel_MPIW$Uncertainty,
#'                                                          CI_width_thresholds = CI_width_thresholds,
#'                                                          interval_width_high_threshold = interval_width_high_threshold,
#'                                                          interval_width_low_threshold = interval_width_low_threshold,
#'                                                          apply_pca = TRUE)
#'
#'       if(!is.null(geno_omic_test_object)){
#'
#'         train_test_label <- ifelse(rownames(geno_omic_test_object)%in%test_label, "Test", "Train")
#'       } else {
#'         train_test_label <- rep("Train", nrow(geno_omic_object))
#'       }
#'       AI_preds <- data.frame(name = GID,
#'                              Predicted_value = AI_pred_reverted,
#'                              Standard_error = pred_SE,
#'                              Train_Test_Label = train_test_label,
#'                              PEV = pred_variances,
#'                              lower_bound = result_rel_MPIW$lower_bound,
#'                              upper_bound = result_rel_MPIW$upper_bound,
#'                              Uncertainty = result_rel_MPIW$Uncertainty,
#'                              Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#'                              Reliability = result_rel$reliability,
#'                              Reliability_remarks = result_rel$remarks,
#'                              Reliability_percentage = result_rel$reliability_percentage,
#'                              #Composite_reliability = composite_reliability$trustworthiness,
#'                              #Composite_reliability_percentage = composite_reliability$reliability_percentage,
#'                              stringsAsFactors = FALSE)
#'
#'
#'       names(AI_preds)[1] <-  c(gen_name)
#'
#'       diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = boot_results,
#'                                                           GID_names = GID,
#'                                                           CI_width_thresholds = CI_width_thresholds,
#'                                                           predictions = AI_pred_reverted,
#'                                                           standard_errors = pred_SE,
#'                                                           prediction_error_var = pred_variances,
#'                                                           genetic_var = genetic_var,
#'                                                           confidence_level = 0.95,
#'                                                           model_for_CI_cal = "ML",
#'                                                           #composite_reliability_score = composite_reliability$reliability_score,
#'                                                           #composite_reliability = composite_reliability$trustworthiness,
#'                                                           #composite_reliability_percentage = composite_reliability$reliability_percentage,
#'                                                           #threshold = NULL,
#'                                                           high_reliability_thres = high_reliability_thres,
#'                                                           low_reliability_thres = low_reliability_thres,
#'                                                           system_database = system_database)
#'     }
#' # browser()
#' # print(num_hidden_layers)
#' # print(learning_rate)
#' # print(neurons_per_layer)
#'   # model_para <- data.frame(stat = c("num_hidden_layers","learning_rate", "neurons_per_layer"),
#'   #                          summary = c(num_hidden_layers, learning_rate, neurons_per_layer),
#'   #                          stringsAsFactors = FALSE)
#'
#' model_para <- data.frame(
#'                       stat = c("num_hidden_layers", "learning_rate", "neurons_per_layer"),
#'                       summary = I(list(num_hidden_layers, learning_rate, neurons_per_layer)),
#'                       stringsAsFactors = FALSE
#'                     )
#'
#' if(length(neurons_per_layer)>1){
#'   # Flatten the list in the 'summary' column and replicate 'stat' values accordingly
#'   expanded_stat <- unlist(lapply(1:nrow(model_para), function(i) {
#'     rep(model_para$stat[i], length(model_para$summary[[i]]))
#'   }))
#'
#'   expanded_summary <- unlist(model_para$summary)
#'
#'   # Create a new data frame with the expanded values
#'   model_para <- data.frame(stat = expanded_stat, summary = expanded_summary, stringsAsFactors = FALSE)
#'
#' }
#'
#'   colnames(model_para)[1:2] <- c("stat", "summary")
#'
#'   output <-  list(model_parameters = model_para,
#'                   predicted_values = AI_preds,
#'                   diagnostic_plots = diagnostic_plots)
#'
#'
#'   return(output)
#'
#' }
#'
#'
