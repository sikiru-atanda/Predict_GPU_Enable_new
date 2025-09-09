
gget <- function(x, key, default = NULL) {
  if (is.null(x)) return(default)
  if (is.list(x)) {
    v <- x[[key]]
    return(if (!is.null(v) && length(v)) v else default)
  }
  if (is.environment(x)) {
    if (exists(key, envir = x, inherits = FALSE)) {
      v <- get(key, envir = x, inherits = FALSE)
      return(if (!is.null(v) && length(v)) v else default)
    }
    return(default)
  }
  default
}

gget_chr <- function(x, key, default) {
  v <- gget(x, key, default)
  if (is.function(v)) return(NA_character_)
  tryCatch(as.character(v), error = function(e) NA_character_)
}

lock_best_hp <- function(combos, best_idx, env = parent.frame()) {
  if (is.null(combos) || is.na(best_idx) || best_idx < 1L) return(invisible(NULL))
  hp <- combos[best_idx, , drop = FALSE]

  #gget <- function(df, key) if (key %in% names(df)) df[[key]][[1]] else NULL
  parse_int_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.numeric(x)) return(as.integer(x))
    if (is.list(x))    return(as.integer(unlist(x, use.names = FALSE)))
    if (is.character(x)) return(as.integer(strsplit(gsub("\\s+", "", x), ",")[[1]]))
    NULL
  }
  parse_bool <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) return(tolower(x) %in% c("true","t","1","yes","y"))
    NULL
  }
  set_if <- function(name, val) if (!is.null(val) && !is.na(val)) assign(name, val, envir = env)

  ## Core training knobs
  set_if("num_hidden_layers",   as.integer(gget(hp, "num_hidden_layers")))
  # neurons_per_layer may be a list/char; parse to int vector
  npl <- parse_int_vec(gget(hp, "neurons_per_layer"))
  set_if("neurons_per_layer", npl)

  set_if("learning_rate_dp",  as.numeric(gget(hp, "learning_rate")))
  set_if("dropout_rate",      as.numeric(gget(hp, "dropout_rate")))
  set_if("epochs",            as.integer(gget(hp, "epochs")))
  set_if("batch_size",        as.integer(gget(hp, "batch_size")))
  set_if("l2_regularizer_dp", as.numeric(gget(hp, "l2_regularizer_dp")))
  set_if("compile_model",     isTRUE(parse_bool(gget(hp, "compile_model"))))
  set_if("deterministic",     isTRUE(parse_bool(gget(hp, "deterministic"))))
  set_if("random_seed",       as.integer(gget(hp, "random_seed")))
  set_if("device",            as.character(gget(hp, "device")))

  ## Attention / BN
  set_if("attention_on_final_layer",      isTRUE(parse_bool(gget(hp, "attention_on_final_layer"))))
  set_if("attention_across_multiple_layers", isTRUE(parse_bool(gget(hp, "attention_across_multiple_layers"))))
  set_if("batch_normalization",           isTRUE(parse_bool(gget(hp, "batch_normalization"))))

  ## CNN (+ max-pool)
  set_if("kernel_size",       as.integer(gget(hp, "kernel_size")))
  set_if("dense_layers_cnn",  parse_int_vec(gget(hp, "dense_layers_cnn")))
  set_if("cnn_use_max_pool",  isTRUE(parse_bool(gget(hp, "cnn_use_max_pool"))))
  set_if("cnn_pool_kernel",   as.integer(gget(hp, "cnn_pool_kernel")))
  set_if("cnn_pool_stride",   as.integer(gget(hp, "cnn_pool_stride")))
  set_if("cnn_pool_padding",  as.integer(gget(hp, "cnn_pool_padding")))



  ## Optimizer / AMP / grad-clip / class weights
  set_if("optimizer_name",     as.character(gget(hp, "optimizer_name")))
  set_if("use_amp",            isTRUE(parse_bool(gget(hp, "use_amp"))))
  set_if("max_grad_norm",      as.numeric(gget(hp, "max_grad_norm")))
  set_if("auto_class_weights", isTRUE(parse_bool(gget(hp, "auto_class_weights"))))

  ## Regression variant
  set_if("heteroscedastic",  isTRUE(parse_bool(gget(hp, "heteroscedastic"))))

  ## FT-Transformer
  set_if("ft_d_model",        as.integer(gget(hp, "ft_d_model")))
  set_if("ft_heads",          as.integer(gget(hp, "ft_heads")))
  set_if("ft_layers",         as.integer(gget(hp, "ft_layers")))
  set_if("ft_ff_mult",        as.integer(gget(hp, "ft_ff_mult")))
  set_if("ft_dropout",        as.numeric(gget(hp, "ft_dropout")))
  set_if("ft_token_dropdown", as.numeric(gget(hp, "ft_token_dropout")))  # harmless if name differs
  set_if("ft_token_dropout",  as.numeric(gget(hp, "ft_token_dropout")))
  set_if("ft_use_cls",        isTRUE(parse_bool(gget(hp, "ft_use_cls"))))

  ## SAINT
  set_if("saint_d_model",       as.integer(gget(hp, "saint_d_model")))
  set_if("saint_heads",         as.integer(gget(hp, "saint_heads")))
  set_if("saint_layers",        as.integer(gget(hp, "saint_layers")))
  set_if("saint_ff_mult",       as.integer(gget(hp, "saint_ff_mult")))
  set_if("saint_dropout",       as.numeric(gget(hp, "saint_dropout")))
  set_if("saint_token_dropout", as.numeric(gget(hp, "saint_token_dropout")))
  set_if("saint_use_cls",       isTRUE(parse_bool(gget(hp, "saint_use_cls"))))

  ## Grouping (FT/SAINT + NAM/MoE)
  set_if("use_grouping",    isTRUE(parse_bool(gget(hp, "use_grouping"))))
  set_if("group_trigger",   as.integer(gget(hp, "group_trigger")))
  set_if("group_method",    as.character(gget(hp, "group_method")))
  set_if("init_group_size", as.integer(gget(hp, "init_group_size")))
  set_if("max_tokens",      as.integer(gget(hp, "max_tokens")))
  set_if("kmeans_batch",    as.integer(gget(hp, "kmeans_batch")))
  set_if("kmeans_iter",     as.integer(gget(hp, "kmeans_iter")))

  ## TabNet
  set_if("tabnet_steps",         as.integer(gget(hp, "tabnet_steps")))
  set_if("tabnet_feature_dim",   as.integer(gget(hp, "tabnet_feature_dim")))
  set_if("tabnet_output_dim",    as.integer(gget(hp, "tabnet_output_dim")))
  set_if("tabnet_gamma",         as.numeric(gget(hp, "tabnet_gamma")))
  set_if("tabnet_lambda_sparse", as.numeric(gget(hp, "tabnet_lambda_sparse")))

  ## NODE
  set_if("node_trees", as.integer(gget(hp, "node_trees")))
  set_if("node_depth", as.integer(gget(hp, "node_depth")))

  ## DeepFM / DCN
  set_if("deepfm_k",      as.integer(gget(hp, "deepfm_k")))
  set_if("deepfm_hidden", parse_int_vec(gget(hp, "deepfm_hidden")))
  set_if("dcn_layers",    as.integer(gget(hp, "dcn_layers")))
  set_if("dcn_hidden",    parse_int_vec(gget(hp, "dcn_hidden")))

  ## NAM
  set_if("nam_hidden",     parse_int_vec(gget(hp, "nam_hidden")))
  set_if("nam_activation", as.character(gget(hp, "nam_activation")))
  set_if("nam_add_linear", isTRUE(parse_bool(gget(hp, "nam_add_linear"))))
  set_if("nam_l1",         as.numeric(gget(hp, "nam_l1")))

  ## MoE
  set_if("moe_n_experts",     as.integer(gget(hp, "moe_n_experts")))
  set_if("moe_expert_hidden", parse_int_vec(gget(hp, "moe_expert_hidden")))
  set_if("moe_gate_hidden",   as.integer(gget(hp, "moe_gate_hidden")))
  set_if("moe_temperature",   as.numeric(gget(hp, "moe_temperature")))
  set_if("moe_sparse_topk",   { v <- gget(hp, "moe_sparse_topk"); if (is.null(v)) NULL else as.integer(v) })
  set_if("moe_entropy_reg",   as.numeric(gget(hp, "moe_entropy_reg")))

  ## GP-DKL / RFF
  set_if("gp_use_variational", isTRUE(parse_bool(gget(hp, "gp_use_variational"))))
  set_if("gp_num_inducing",    as.integer(gget(hp, "gp_num_inducing")))
  set_if("gp_feature_dim",     as.integer(gget(hp, "gp_feature_dim")))
  set_if("gp_kernel",          as.character(gget(hp, "gp_kernel")))
  set_if("gp_ard",             isTRUE(parse_bool(gget(hp, "gp_ard"))))
  set_if("gp_lr_mult",         as.numeric(gget(hp, "gp_lr_mult")))
  set_if("rff_features",       as.integer(gget(hp, "rff_features")))
  set_if("rff_lengthscale",    as.numeric(gget(hp, "rff_lengthscale")))
  set_if("rff_deep_hidden",    parse_int_vec(gget(hp, "rff_deep_hidden")))

  invisible(NULL)
}

# ------------------------------------------------------------------
# Neuron parsing & validation helpers
# ------------------------------------------------------------------
parse_neurons <- function(x, input_size, num_hidden_layers, scaling = 0.5) {
  # Accept NULL / "" / "auto" / "dynamic" -> generate dynamically
  if (is.null(x)) return(generate_dynamic_layers(input_size, num_hidden_layers, scaling))
  if (is.character(x) && length(x) == 1L && nchar(trimws(x)) == 0L) {
    return(generate_dynamic_layers(input_size, num_hidden_layers, scaling))
  }
  if (is.list(x)) x <- unlist(x, use.names = FALSE)

  if (is.character(x)) {
    key <- tolower(trimws(x[1L]))
    if (key %in% c("auto", "dynamic")) {
      return(generate_dynamic_layers(input_size, num_hidden_layers, scaling))
    }
    x <- as.integer(strsplit(x, "\\s*,\\s*")[[1]])
  }
  as.integer(x)
}

# Robust validator; only enforces when BOTH are provided and meaningful.
validate_layers <- function(num_hidden_layers, neurons_per_layer) {
  msg <- "\n==================================================\n"
  # Nothing to validate
  if (is.null(num_hidden_layers) || is.null(neurons_per_layer)) return(TRUE)

  # Scalar case
  if (is.numeric(num_hidden_layers) && !is.list(neurons_per_layer)) {
    # Allow "auto"/"dynamic" and empty strings
    if (is.character(neurons_per_layer)) {
      s <- tolower(trimws(neurons_per_layer[1]))
      if (s %in% c("auto", "dynamic", "")) return(TRUE)
      neurons_per_layer <- as.integer(strsplit(s, "\\s*,\\s*")[[1]])
    }
    if (length(neurons_per_layer) != as.integer(num_hidden_layers)) {
      stop(paste(msg, "The length of neurons_per_layer must equal num_hidden_layers"), call. = FALSE)
    }
    return(TRUE)
  }

  # Vector/list case
  if (length(num_hidden_layers) != length(neurons_per_layer)) {
    stop(paste(msg, "num_hidden_layers and neurons_per_layer vectors must have equal length"), call. = FALSE)
  }
  for (i in seq_along(num_hidden_layers)) {
    k <- suppressWarnings(as.integer(num_hidden_layers[[i]]))
    v <- neurons_per_layer[[i]]
    if (is.null(k) || is.na(k)) next
    if (is.character(v)) {
      s <- tolower(trimws(v[1]))
      if (s %in% c("auto", "dynamic", "")) next
      v <- strsplit(v, "\\s*,\\s*")[[1]]
    }
    if (is.list(v)) v <- unlist(v, use.names = FALSE)
    if (length(v) != k) {
      stop(paste(msg, "Mismatch at index", i, ": expected", k, "layers, got", length(v)), call. = FALSE)
    }
  }
  TRUE
}

# Internal: normalize a param_grid (named list) to a data.frame of combinations
norm_param_grid <- function(param_grid) {
  if (is.null(param_grid) || !length(param_grid)) {
    stop("param_grid must be a non-empty named list.", call. = FALSE)
  }
  as.data.frame(do.call(expand.grid, c(param_grid, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))))
}

# Map user-facing model names to Python codes (extended)
# .map_model_type <- function(x) {
#   switch(tolower(x),
#          "mlp"                = "mlp",
#          "mlp_with_attention" = "mlp_with_attention",
#          "resnet"             = "resnet",
#          "cnn"                = "cnn",
#          "ft"                 = "ft_transformer",
#          "ft_transformer"     = "ft_transformer",
#          "saint"              = "saint",
#          "tabnet"             = "tabnet",
#          "node"               = "node",
#          "deepfm"             = "deepfm",
#          "dcn"                = "dcnv2",
#          "dcnv2"              = "dcnv2",
#          "nam"                = "nam",
#          "moe"                = "moe",
#          "gp_dkl"             = "gp_dkl",
#          x  # passthrough; Python will error if unknown
#   )
# }

# ------------------------------------------------------------------
# Grid search (now includes FT/SAINT grouping, NAM, MoE, GP-DKL, AMP, etc.)
# ------------------------------------------------------------------
grid_search_deep_learning <- function(
    X_train, y_train, param_grid,
    # Common training defaults
    epochs, batch_size, validation_split = 0.2,
    learning_rate, dropout_rate = 0.5,
    l2_regularizer_dp,
    compile_model = TRUE,
    deterministic = TRUE,
    random_seed = 123,
    device = NULL,
    early_stop = TRUE,  # signature compatibility; Python does early-stop internally

    # Optimizer / AMP / clip / weights
    optimizer_name = "adam",
    use_amp = TRUE,
    max_grad_norm = 1.0,
    auto_class_weights = FALSE,

    # Classic CNN defaults
    kernel_size = 3,
    dense_layers_cnn = c(256, 128, 64),
    # CNN max-pool
    cnn_use_max_pool = FALSE,
    cnn_pool_kernel  = 2,
    cnn_pool_stride  = 2,
    cnn_pool_padding = 0,

    # ResNet defaults
    n_blocks = 2,
    n_neurons_per_block = NULL,

    # Model choice (can be overridden per-row in param_grid)
    deep_learning_model = "mlp_with_attention",
    attention_on_final_layer = TRUE,
    attention_across_multiple_layers = FALSE,
    batch_normalization = TRUE,

    # Heteroscedastic regression
    heteroscedastic = FALSE,

    # FT-Transformer
    ft_d_model = 192, ft_heads = 8, ft_layers = 4, ft_ff_mult = 4,
    ft_dropout = 0.1, ft_token_dropout = 0.1, ft_use_cls = TRUE,

    # SAINT
    saint_d_model = 128, saint_heads = 8, saint_layers = 4, saint_ff_mult = 4,
    saint_dropout = 0.1, saint_token_dropout = 0.1, saint_use_cls = TRUE,

    # Grouping controls (FT/SAINT + NAM/MoE)
    use_grouping = TRUE, group_trigger = 2048, group_method = "auto",
    init_group_size = 64, max_tokens = 1024, kmeans_batch = 4096, kmeans_iter = 100,

    # TabNet
    tabnet_steps = 5, tabnet_feature_dim = 64, tabnet_output_dim = 64,
    tabnet_gamma = 1.5, tabnet_lambda_sparse = 1e-4,

    # NODE
    node_trees = 8, node_depth = 3,

    # DeepFM
    deepfm_k = 16, deepfm_hidden = c(128, 64),

    # DCN / DCNv2
    dcn_layers = 3, dcn_hidden = c(256, 128),

    # NAM
    nam_hidden = c(32, 16), nam_activation = "relu",
    nam_add_linear = TRUE, nam_l1 = 1e-4,

    # MoE
    moe_n_experts = 4, moe_expert_hidden = c(128, 64),
    moe_gate_hidden = 128, moe_temperature = 1.0,
    moe_sparse_topk = NULL, moe_entropy_reg = 0.0,

    # GP-DKL (+ RFF fallback)
    gp_use_variational = TRUE, gp_num_inducing = 512, gp_feature_dim = 64,
    gp_kernel = "rbf", gp_ard = TRUE, gp_lr_mult = 0.5,
    rff_features = 1024, rff_lengthscale = 1.0, rff_deep_hidden = c(128)
) {
  X <- as.matrix(X_train)
  y <- as.numeric(y_train)

  # helpers
  parse_int_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.numeric(x)) return(as.integer(x))
    if (is.list(x))    return(as.integer(unlist(x)))
    if (is.character(x)) return(as.integer(strsplit(gsub("\\s+", "", x), ",")[[1]]))
    stop("Expected numeric/list/character vector for an integer vector hyperparameter.")
  }
  parse_bool <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) return(tolower(x) %in% c("true","t","1","yes","y"))
    FALSE
  }
  # gget <- function(lst, name, default) {
  #   if (name %in% names(lst) && !is.null(lst[[name]]) && !is.na(lst[[name]])) lst[[name]] else default
  # }
  map_model <- function(s) {
    switch(tolower(s),
           "mlp"                = "mlp",
           "mlp_with_attention" = "mlp_with_attention",
           "resnet"             = "resnet",
           "cnn"                = "cnn",
           "ft"                 = "ft_transformer",
           "ft_transformer"     = "ft_transformer",
           "saint"              = "saint",
           "tabnet"             = "tabnet",
           "node"               = "node",
           "deepfm"             = "deepfm",
           "dcn"                = "dcnv2",
           "dcnv2"              = "dcnv2",
           "nam"                = "nam",
           "moe"                = "moe",
           "gp_dkl"             = "gp_dkl",
           stop("Unsupported deep_learning_model: ", s)
    )
  }

  check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
    neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
    return(length(neuron_list) == num_hidden_layers)
  }
  combos <- expand.grid(param_grid, stringsAsFactors = FALSE)
  if (nrow(combos) == 0L) stop("param_grid produced 0 combinations.")

  valid_rows <- apply(combos, 1, function(row) {
    num_hidden_layers <- as.integer(row["num_hidden_layers"])
    neurons_per_layer <- row["neurons_per_layer"]
    check_hidden_layers(num_hidden_layers, neurons_per_layer)
  })

  # Filter the data frame to keep only valid rows
  combos <- combos[valid_rows, ]
  # validate hidden shapes when relevant
  check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
    neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
    return(length(neuron_list) == num_hidden_layers)
  }
  # if ("num_hidden_layers" %in% names(combos) && "neurons_per_layer" %in% names(combos)) {
  #   uses_hidden_row <- function(row_model) {
  #     mt <- map_model(row_model)
  #     mt %in% c("mlp","mlp_with_attention","resnet","cnn")
  #   }
  #   keep <- logical(nrow(combos))
  #   for (i in seq_len(nrow(combos))) {
  #     if (!uses_hidden_row(gget(combos[i,], "deep_learning_model", deep_learning_model))) {
  #       keep[i] <- TRUE
  #     } else {
  #       nl <- suppressWarnings(as.integer(combos$num_hidden_layers[i]))
  #       np <- combos$neurons_per_layer[i]
  #       keep[i] <- !is.na(nl) && !is.null(np) && check_hidden_layers(nl, np)
  #     }
  #   }
  #   combos <- combos[keep, , drop = FALSE]
  #   if (nrow(combos) == 0L) stop("All combinations filtered out by hidden-layer checks.")
  # }
  if ("num_hidden_layers" %in% names(combos) && "neurons_per_layer" %in% names(combos)) {

    uses_hidden_row <- function(row_model) {
      rm_chr <- tryCatch(
        if (is.function(row_model)) NA_character_ else as.character(row_model),
        error = function(e) NA_character_
      )
      mt <- tryCatch(map_model(rm_chr), error = function(e) NA_character_)
      isTRUE(mt %in% c("mlp", "mlp_with_attention", "resnet", "cnn"))
    }

    keep <- logical(nrow(combos))
    for (i in seq_len(nrow(combos))) {
      mdl <- gget_chr(combos[i, , drop = FALSE], "deep_learning_model",
                      default = tryCatch(as.character(deep_learning_model), error = function(e) NA_character_))

      if (!uses_hidden_row(mdl)) {
        keep[i] <- TRUE
      } else {
        nl <- suppressWarnings(as.integer(combos$num_hidden_layers[i]))
        np <- combos$neurons_per_layer[i]
        # If np is a factor/list, coerce to character first
        if (is.factor(np)) np <- as.character(np)
        keep[i] <- !is.na(nl) && !is.null(np) && check_hidden_layers(nl, np)
      }
    }

    combos <- combos[keep, , drop = FALSE]
    if (nrow(combos) == 0L) stop("All combinations filtered out by hidden-layer checks.")
  }

  results <- vector("list", length = nrow(combos))

  for (i in seq_len(nrow(combos))) {
    hp <- combos[i, , drop = FALSE]

    model_type <- map_model(gget(hp, "deep_learning_model", deep_learning_model))
    uses_hidden <- model_type %in% c("mlp","mlp_with_attention","resnet","cnn")

    # hidden
    nlayers <- gget(hp, "num_hidden_layers", NULL); if (!is.null(nlayers)) nlayers <- as.integer(nlayers)
    neurons <- gget(hp, "neurons_per_layer", NULL); if (!is.null(neurons)) neurons <- parse_int_vec(neurons)

    # resnet mapping
    if (identical(model_type, "resnet")) {
      nb <- gget(hp, "n_blocks", n_blocks)
      nnpb <- gget(hp, "n_neurons_per_block", n_neurons_per_block)
      if (!is.null(nnpb)) neurons <- parse_int_vec(nnpb)
      if (!is.null(nb))   nlayers <- as.integer(nb)
      if (!is.null(nlayers) && !is.null(neurons) && length(neurons) != nlayers) {
        stop(sprintf("Row %d: Mismatch between n_blocks (%s) and n_neurons_per_block length (%s).", i, nlayers, length(neurons)))
      }
      if (is.null(neurons)) stop(sprintf("Row %d: For ResNet, provide n_neurons_per_block.", i))
    }
    if (uses_hidden && !is.null(nlayers) && !is.null(neurons) && length(neurons) != nlayers) {
      stop(sprintf("Row %d: length(neurons_per_layer) must equal num_hidden_layers.", i))
    }

    # core scalars
    lr  <- as.numeric(gget(hp, "learning_rate",       learning_rate))
    dr  <- as.numeric(gget(hp, "dropout_rate",        dropout_rate))
    ep  <- as.integer(gget(hp, "epochs",              epochs))
    bs  <- as.integer(gget(hp, "batch_size",          batch_size))
    l2w <- as.numeric(gget(hp, "l2_regularizer_dp",   l2_regularizer_dp))

    cmp <- isTRUE(parse_bool(gget(hp, "compile_model", compile_model)))
    det <- isTRUE(parse_bool(gget(hp, "deterministic", deterministic)))
    rse <- as.integer(gget(hp, "random_seed", random_seed))

    fin_att  <- isTRUE(parse_bool(gget(hp, "attention_on_final_layer", attention_on_final_layer)))
    across   <- isTRUE(parse_bool(gget(hp, "attention_across_multiple_layers", attention_across_multiple_layers)))
    bn       <- isTRUE(parse_bool(gget(hp, "batch_normalization", batch_normalization)))

    hetero <- isTRUE(parse_bool(gget(hp, "heteroscedastic", heteroscedastic)))

    ks   <- as.integer(gget(hp, "kernel_size", kernel_size))
    dlcn <- parse_int_vec(gget(hp, "dense_layers_cnn", dense_layers_cnn))
    mp_use     <- isTRUE(parse_bool(gget(hp, "cnn_use_max_pool", cnn_use_max_pool)))
    mp_kernel  <- as.integer(gget(hp, "cnn_pool_kernel",  cnn_pool_kernel))
    mp_stride  <- as.integer(gget(hp, "cnn_pool_stride",  cnn_pool_stride))
    mp_padding <- as.integer(gget(hp, "cnn_pool_padding", cnn_pool_padding))

    # optim/AMP/clip/weights
    opt_name  <- as.character(gget(hp, "optimizer_name", optimizer_name))
    amp_use   <- isTRUE(parse_bool(gget(hp, "use_amp", use_amp)))
    grad_clip <- as.numeric(gget(hp, "max_grad_norm", max_grad_norm))
    cls_wt    <- isTRUE(parse_bool(gget(hp, "auto_class_weights", auto_class_weights)))

    # FT / SAINT
    ft_par <- list(
      ft_d_model = as.integer(gget(hp, "ft_d_model", ft_d_model)),
      ft_heads   = as.integer(gget(hp, "ft_heads",   ft_heads)),
      ft_layers  = as.integer(gget(hp, "ft_layers",  ft_layers)),
      ft_ff_mult = as.integer(gget(hp, "ft_ff_mult", ft_ff_mult)),
      ft_dropout = as.numeric(gget(hp, "ft_dropout", ft_dropout)),
      ft_token_dropout = as.numeric(gget(hp, "ft_token_dropout", ft_token_dropout)),
      ft_use_cls = isTRUE(parse_bool(gget(hp, "ft_use_cls", ft_use_cls)))
    )
    saint_par <- list(
      saint_d_model = as.integer(gget(hp, "saint_d_model", saint_d_model)),
      saint_heads   = as.integer(gget(hp, "saint_heads",   saint_heads)),
      saint_layers  = as.integer(gget(hp, "saint_layers",  saint_layers)),
      saint_ff_mult = as.integer(gget(hp, "saint_ff_mult", saint_ff_mult)),
      saint_dropout = as.numeric(gget(hp, "saint_dropout", saint_dropout)),
      saint_token_dropout = as.numeric(gget(hp, "saint_token_dropout", saint_token_dropout)),
      saint_use_cls = isTRUE(parse_bool(gget(hp, "saint_use_cls", saint_use_cls)))
    )
    group_par <- list(
      use_grouping = isTRUE(parse_bool(gget(hp, "use_grouping", use_grouping))),
      group_trigger = as.integer(gget(hp, "group_trigger", group_trigger)),
      group_method  = as.character(gget(hp, "group_method", group_method)),
      init_group_size = as.integer(gget(hp, "init_group_size", init_group_size)),
      max_tokens      = as.integer(gget(hp, "max_tokens", max_tokens)),
      kmeans_batch    = as.integer(gget(hp, "kmeans_batch", kmeans_batch)),
      kmeans_iter     = as.integer(gget(hp, "kmeans_iter", kmeans_iter))
    )

    tabnet_par <- list(
      tabnet_steps       = as.integer(gget(hp, "tabnet_steps",       tabnet_steps)),
      tabnet_feature_dim = as.integer(gget(hp, "tabnet_feature_dim", tabnet_feature_dim)),
      tabnet_output_dim  = as.integer(gget(hp, "tabnet_output_dim",  tabnet_output_dim)),
      tabnet_gamma       = as.numeric(gget(hp, "tabnet_gamma",       tabnet_gamma)),
      tabnet_lambda_sparse = as.numeric(gget(hp, "tabnet_lambda_sparse", tabnet_lambda_sparse))
    )
    node_par <- list(
      node_trees = as.integer(gget(hp, "node_trees", node_trees)),
      node_depth = as.integer(gget(hp, "node_depth", node_depth))
    )
    deepfm_par <- list(
      deepfm_k      = as.integer(gget(hp, "deepfm_k", deepfm_k)),
      deepfm_hidden = parse_int_vec(gget(hp, "deepfm_hidden", deepfm_hidden))
    )
    dcn_par <- list(
      dcn_layers = as.integer(gget(hp, "dcn_layers", dcn_layers)),
      dcn_hidden = parse_int_vec(gget(hp, "dcn_hidden", dcn_hidden))
    )
    nam_par <- list(
      nam_hidden     = parse_int_vec(gget(hp, "nam_hidden", nam_hidden)),
      nam_activation = as.character(gget(hp, "nam_activation", nam_activation)),
      nam_add_linear = isTRUE(parse_bool(gget(hp, "nam_add_linear", nam_add_linear))),
      nam_l1         = as.numeric(gget(hp, "nam_l1", nam_l1))
    )
    moe_par <- list(
      moe_n_experts     = as.integer(gget(hp, "moe_n_experts", moe_n_experts)),
      moe_expert_hidden = parse_int_vec(gget(hp, "moe_expert_hidden", moe_expert_hidden)),
      moe_gate_hidden   = as.integer(gget(hp, "moe_gate_hidden", moe_gate_hidden)),
      moe_temperature   = as.numeric(gget(hp, "moe_temperature", moe_temperature)),
      moe_sparse_topk   = { val <- gget(hp, "moe_sparse_topk", moe_sparse_topk); if (is.null(val)) NULL else as.integer(val) },
      moe_entropy_reg   = as.numeric(gget(hp, "moe_entropy_reg", moe_entropy_reg))
    )
    gpdkl_par <- list(
      gp_use_variational = isTRUE(parse_bool(gget(hp, "gp_use_variational", gp_use_variational))),
      gp_num_inducing    = as.integer(gget(hp, "gp_num_inducing", gp_num_inducing)),
      gp_feature_dim     = as.integer(gget(hp, "gp_feature_dim", gp_feature_dim)),
      gp_kernel          = as.character(gget(hp, "gp_kernel", gp_kernel)),
      gp_ard             = isTRUE(parse_bool(gget(hp, "gp_ard", gp_ard))),
      gp_lr_mult         = as.numeric(gget(hp, "gp_lr_mult", gp_lr_mult)),
      rff_features       = as.integer(gget(hp, "rff_features", rff_features)),
      rff_lengthscale    = as.numeric(gget(hp, "rff_lengthscale", rff_lengthscale)),
      rff_deep_hidden    = parse_int_vec(gget(hp, "rff_deep_hidden", rff_deep_hidden))
    )

    args <- list(
      X = X, y = y,
      model_type = model_type,
      learning_rate = lr, epochs = ep, batch_size = bs,
      l2_weight_decay = l2w, dropout = dr,
      optimizer_name = opt_name,
      final_attention = fin_att, attention_across_multiple_layers = across,
      batch_norm = bn, validation_split = as.numeric(validation_split),
      compile_model = cmp, device = if (is.null(device)) NULL else as.character(device),
      deterministic = det, random_seed = rse,
      heteroscedastic = hetero,
      # AMP/clip/weights
      use_amp = amp_use, max_grad_norm = grad_clip, auto_class_weights = cls_wt,
      # CNN
      kernel_size = ks, dense_layers_cnn = dlcn,
      cnn_use_max_pool = mp_use, cnn_pool_kernel = mp_kernel,
      cnn_pool_stride = mp_stride, cnn_pool_padding = mp_padding
    )
    if (uses_hidden) {
      if (!is.null(nlayers)) args$num_hidden_layers <- nlayers
      if (!is.null(neurons)) args$neurons_per_layer <- as.integer(neurons)
    }
    args <- c(args, ft_par, saint_par, group_par, tabnet_par, node_par, deepfm_par, dcn_par, nam_par, moe_par, gpdkl_par)

    fit <- do.call(torch_fit_model, args)

    results[[paste0("Model_", i)]] <- list(
      model   = fit[[1L]],
      history = fit[[2L]],
      config  = args
    )
  }

  results
}


# --- helpers ---------------------------------------------------------------

# Parse neurons from a variety of inputs:
#  - integer vector: c(256,128)
#  - list element:   list(c(256,128))
#  - character:      "256,128"
#  - "auto":         generate from input size + num_hidden_layers
#' parse_neurons <- function(x, input_size, num_hidden_layers, scaling = 0.5) {
#'   if (is.null(x)) return(generate_dynamic_layers(input_size, num_hidden_layers, scaling))
#'   if (is.list(x)) x <- unlist(x, use.names = FALSE)
#'   if (is.character(x)) {
#'     if (tolower(x) %in% c("auto", "dynamic")) {
#'       return(generate_dynamic_layers(input_size, num_hidden_layers, scaling))
#'     }
#'     x <- as.integer(strsplit(x, "\\s*,\\s*")[[1]])
#'   }
#'   as.integer(x)
#' }
#'
#' # Robust validator (vectorized or per-row)
#' validate_layers <- function(num_hidden_layers, neurons_per_layer) {
#'   msg <- "\n==================================================\n"
#'   # vectorized case: both are length-1 scalars for a single model
#'   if (is.numeric(num_hidden_layers) && !is.list(neurons_per_layer)) {
#'     if (length(neurons_per_layer) != as.integer(num_hidden_layers)) {
#'       stop(paste(msg, "The length of neurons_per_layer must equal num_hidden_layers"), call. = FALSE)
#'     }
#'     return(TRUE)
#'   }
#'   # list/data.frame case: check each pair
#'   if (length(num_hidden_layers) != length(neurons_per_layer)) {
#'     stop(paste(msg, "num_hidden_layers and neurons_per_layer vectors must have equal length"), call. = FALSE)
#'   }
#'   for (i in seq_along(num_hidden_layers)) {
#'     k <- as.integer(num_hidden_layers[[i]])
#'     v <- neurons_per_layer[[i]]
#'     if (is.character(v)) v <- strsplit(v, "\\s*,\\s*")[[1]]
#'     if (is.list(v)) v <- unlist(v, use.names = FALSE)
#'     if (length(v) != k) {
#'       stop(paste(msg, "Mismatch at index", i, ": expected", k, "layers, got", length(v)), call. = FALSE)
#'     }
#'   }
#'   TRUE
#' }
#'
#' # Internal: normalize a param_grid (named list) to a data.frame of combinations
#' norm_param_grid <- function(param_grid) {
#'   if (is.null(param_grid) || !length(param_grid)) {
#'     stop("param_grid must be a non-empty named list.", call. = FALSE)
#'   }
#'   as.data.frame(do.call(expand.grid, c(param_grid, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))))
#' }
#'
#' # Map user-facing model names to Python codes
#' .map_model_type <- function(x) {
#'   switch(tolower(x),
#'          "resnet" = "resnet",
#'          "cnn" = "cnn",
#'          "mlp_with_attention" = "mlp_with_attention",
#'          "mlp" = "mlp",
#'          x
#'   )
#' }
#'
#' # --- main: grid search -----------------------------------------------------
#'
#' #' Grid search for deep learning (PyTorch backend)
#' #' @param X_train matrix/data.frame
#' #' @param y_train numeric/labels
#' #' @param param_grid named list of candidate values; supports:
#' #'   num_hidden_layers, neurons_per_layer, learning_rate, dropout_rate, epochs,
#' #'   batch_size, l2_regularizer_dp, deep_learning_model, n_blocks, n_neurons_per_block,
#' #'   kernel_size, dense_layers_cnn, attention_on_final_layer, attention_across_multiple_layers,
#' #'   batch_normalization, validation_split
#' #' @param parallel logical; use future.apply if TRUE
#' #' @param workers optional integer for future::plan()
#' #' @return named list: each item has $model (python object) and $history (list with train_loss/val_loss)
#' #' @export
#' grid_search_deep_learning <- function(
#'     X_train, y_train, param_grid,
#'     # Common training defaults
#'     epochs, batch_size, validation_split = 0.2,
#'     learning_rate, dropout_rate = 0.5,
#'     l2_regularizer_dp,
#'     compile_model = TRUE,
#'     deterministic = TRUE,
#'     random_seed = 123,
#'     device = NULL,
#'     early_stop = TRUE,  # kept for signature compatibility (Python handles early stopping internally)
#'
#'     # Classic CNN defaults
#'     kernel_size = 3,
#'     dense_layers_cnn = c(256, 128, 64),
#'     # CNN max-pool knobs (NEW)
#'     cnn_use_max_pool = FALSE,
#'     cnn_pool_kernel  = 2,
#'     cnn_pool_stride  = 2,
#'     cnn_pool_padding = 0,
#'
#'     # ResNet defaults
#'     n_blocks = 2,
#'     n_neurons_per_block = NULL,
#'
#'     # Model choice (can be overridden per-row in param_grid)
#'     deep_learning_model = "mlp_with_attention",
#'     attention_on_final_layer = TRUE,
#'     attention_across_multiple_layers = FALSE,
#'     batch_normalization = TRUE,
#'
#'     # Heteroscedastic regression (NEW; only effective for regression tasks)
#'     heteroscedastic = FALSE,
#'
#'     # FT-Transformer (NEW)
#'     ft_d_model = 192, ft_heads = 8, ft_layers = 4, ft_ff_mult = 4,
#'     ft_dropout = 0.1, ft_token_dropout = 0.1, ft_use_cls = TRUE,
#'
#'     # SAINT (NEW)
#'     saint_d_model = 128, saint_heads = 8, saint_layers = 4, saint_ff_mult = 4,
#'     saint_dropout = 0.1, saint_token_dropout = 0.1,
#'
#'     # TabNet (NEW)
#'     tabnet_steps = 5, tabnet_feature_dim = 64, tabnet_output_dim = 64,
#'     tabnet_gamma = 1.5, tabnet_lambda_sparse = 1e-4,
#'
#'     # NODE (NEW)
#'     node_trees = 8, node_depth = 3,
#'
#'     # DeepFM (NEW)
#'     deepfm_k = 16, deepfm_hidden = c(128, 64),
#'
#'     # DCN / DCNv2 (NEW)
#'     dcn_layers = 3, dcn_hidden = c(256, 128)
#' ) {
#'   X <- as.matrix(X_train)
#'   y <- as.numeric(y_train)
#'
#'   # --- helpers ---------------------------------------------------------
#'   parse_int_vec <- function(x) {
#'     if (is.null(x)) return(NULL)
#'     if (is.numeric(x)) return(as.integer(x))
#'     if (is.list(x))    return(as.integer(unlist(x)))
#'     if (is.character(x)) {
#'       s <- gsub("\\s+", "", x)
#'       return(as.integer(strsplit(s, ",")[[1]]))
#'     }
#'     stop("Expected numeric/list/character vector for an integer vector hyperparameter.")
#'   }
#'   parse_bool <- function(x) {
#'     if (is.null(x)) return(NULL)
#'     if (is.logical(x)) return(x)
#'     if (is.numeric(x)) return(x != 0)
#'     if (is.character(x)) return(tolower(x) %in% c("true", "t", "1", "yes", "y"))
#'     FALSE
#'   }
#'   gget <- function(lst, name, default) {
#'     if (name %in% names(lst) && !is.null(lst[[name]]) && !is.na(lst[[name]])) lst[[name]] else default
#'   }
#'   map_model <- function(s) {
#'     switch(tolower(s),
#'            "mlp"                = "mlp",
#'            "mlp_with_attention" = "mlp_with_attention",
#'            "resnet"             = "resnet",
#'            "cnn"                = "cnn",
#'            "ft"                 = "ft_transformer",
#'            "ft_transformer"     = "ft_transformer",
#'            "saint"              = "saint",
#'            "tabnet"             = "tabnet",
#'            "node"               = "node",
#'            "deepfm"             = "deepfm",
#'            "dcn"                = "dcnv2",
#'            "dcnv2"              = "dcnv2",
#'            stop("Unsupported deep_learning_model: ", s)
#'     )
#'   }
#'
#'   # Expand the grid
#'   combos <- expand.grid(param_grid, stringsAsFactors = FALSE)
#'
#'   check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
#'     neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
#'     return(length(neuron_list) == num_hidden_layers)
#'   }
#'
#'   # Loop through the data frame and remove rows with mismatched hidden layers and neurons
#'   valid_rows <- apply(combos, 1, function(row) {
#'     num_hidden_layers <- as.integer(row["num_hidden_layers"])
#'     neurons_per_layer <- row["neurons_per_layer"]
#'     check_hidden_layers(num_hidden_layers, neurons_per_layer)
#'   })
#'
#'   # Filter the data frame to keep only valid rows
#'   combos <- combos[valid_rows, ]
#'   if (nrow(combos) == 0L) stop("param_grid produced 0 combinations.")
#'
#'   results <- vector("list", length = nrow(combos))
#'
#'   for (i in seq_len(nrow(combos))) {
#'     hp <- combos[i, , drop = FALSE]
#'
#'     # Per-combination model choice (falls back to function default)
#'     model_type <- map_model(gget(hp, "deep_learning_model", deep_learning_model))
#'
#'     # These four models use explicit hidden layer specs
#'     uses_hidden <- model_type %in% c("mlp", "mlp_with_attention", "resnet", "cnn")
#'
#'     # --- Hidden-layer handling ----------------------------------------
#'     # num_hidden_layers
#'     nlayers <- gget(hp, "num_hidden_layers", NULL)
#'     if (!is.null(nlayers)) nlayers <- as.integer(nlayers)
#'
#'     # neurons_per_layer
#'     neurons <- gget(hp, "neurons_per_layer", NULL)
#'     if (!is.null(neurons)) neurons <- parse_int_vec(neurons)
#'
#'     # ResNet block mapping
#'     if (identical(model_type, "resnet")) {
#'       # Allow n_blocks / n_neurons_per_block in grid to override function defaults
#'       nb <- gget(hp, "n_blocks", n_blocks)
#'       nnpb <- gget(hp, "n_neurons_per_block", n_neurons_per_block)
#'       if (!is.null(nnpb)) neurons <- parse_int_vec(nnpb)
#'       if (!is.null(nb))   nlayers <- as.integer(nb)
#'       if (!is.null(nlayers) && !is.null(neurons) && length(neurons) != nlayers) {
#'         stop(sprintf("Row %d: Mismatch between n_blocks (%s) and n_neurons_per_block length (%s).",
#'                      i, nlayers, length(neurons)))
#'       }
#'       if (is.null(neurons)) stop(sprintf("Row %d: For ResNet, provide n_neurons_per_block.", i))
#'     }
#'
#'     # Provide safe defaults for models that need hidden layers (if not given)
#'     if (uses_hidden && is.null(neurons)) {
#'       neurons <- if (identical(model_type, "cnn")) c(64, 64, 64) else c(256, 128)
#'     }
#'
#'     # If num_hidden_layers was provided, validate for uses_hidden models
#'     if (uses_hidden && !is.null(nlayers) && !is.null(neurons) && length(neurons) != nlayers) {
#'       stop(sprintf("Row %d: length(neurons_per_layer) must equal num_hidden_layers.", i))
#'     }
#'
#'     # --- Scalars with per-row overrides -------------------------------
#'     lr  <- as.numeric(gget(hp, "learning_rate",       learning_rate))
#'     dr  <- as.numeric(gget(hp, "dropout_rate",        dropout_rate))
#'     ep  <- as.integer(gget(hp, "epochs",              epochs))
#'     bs  <- as.integer(gget(hp, "batch_size",          batch_size))
#'     l2w <- as.numeric(gget(hp, "l2_regularizer_dp",   l2_regularizer_dp))
#'
#'     # Compilation & reproducibility
#'     cmp <- isTRUE(parse_bool(gget(hp, "compile_model", compile_model)))
#'     det <- isTRUE(parse_bool(gget(hp, "deterministic", deterministic)))
#'     rse <- as.integer(gget(hp, "random_seed", random_seed))
#'
#'     # Attention/BN
#'     fin_att  <- isTRUE(parse_bool(gget(hp, "attention_on_final_layer", attention_on_final_layer)))
#'     across   <- isTRUE(parse_bool(gget(hp, "attention_across_multiple_layers", attention_across_multiple_layers)))
#'     bn       <- isTRUE(parse_bool(gget(hp, "batch_normalization", batch_normalization)))
#'
#'     # Heteroscedastic regression
#'     hetero <- isTRUE(parse_bool(gget(hp, "heteroscedastic", heteroscedastic)))
#'
#'     # CNN extras
#'     ks   <- as.integer(gget(hp, "kernel_size", kernel_size))
#'     dlcn <- parse_int_vec(gget(hp, "dense_layers_cnn", dense_layers_cnn))
#'     mp_use     <- isTRUE(parse_bool(gget(hp, "cnn_use_max_pool", cnn_use_max_pool)))
#'     mp_kernel  <- as.integer(gget(hp, "cnn_pool_kernel",  cnn_pool_kernel))
#'     mp_stride  <- as.integer(gget(hp, "cnn_pool_stride",  cnn_pool_stride))
#'     mp_padding <- as.integer(gget(hp, "cnn_pool_padding", cnn_pool_padding))
#'
#'     # FT-Transformer
#'     ft_par <- list(
#'       ft_d_model = as.integer(gget(hp, "ft_d_model", ft_d_model)),
#'       ft_heads   = as.integer(gget(hp, "ft_heads",   ft_heads)),
#'       ft_layers  = as.integer(gget(hp, "ft_layers",  ft_layers)),
#'       ft_ff_mult = as.integer(gget(hp, "ft_ff_mult", ft_ff_mult)),
#'       ft_dropout = as.numeric(gget(hp, "ft_dropout", ft_dropout)),
#'       ft_token_dropout = as.numeric(gget(hp, "ft_token_dropout", ft_token_dropout)),
#'       ft_use_cls = isTRUE(parse_bool(gget(hp, "ft_use_cls", ft_use_cls)))
#'     )
#'
#'     # SAINT
#'     saint_par <- list(
#'       saint_d_model = as.integer(gget(hp, "saint_d_model", saint_d_model)),
#'       saint_heads   = as.integer(gget(hp, "saint_heads",   saint_heads)),
#'       saint_layers  = as.integer(gget(hp, "saint_layers",  saint_layers)),
#'       saint_ff_mult = as.integer(gget(hp, "saint_ff_mult", saint_ff_mult)),
#'       saint_dropout = as.numeric(gget(hp, "saint_dropout", saint_dropout)),
#'       saint_token_dropout = as.numeric(gget(hp, "saint_token_dropout", saint_token_dropout))
#'     )
#'
#'     # TabNet
#'     tabnet_par <- list(
#'       tabnet_steps       = as.integer(gget(hp, "tabnet_steps",       tabnet_steps)),
#'       tabnet_feature_dim = as.integer(gget(hp, "tabnet_feature_dim", tabnet_feature_dim)),
#'       tabnet_output_dim  = as.integer(gget(hp, "tabnet_output_dim",  tabnet_output_dim)),
#'       tabnet_gamma       = as.numeric(gget(hp, "tabnet_gamma",       tabnet_gamma)),
#'       tabnet_lambda_sparse = as.numeric(gget(hp, "tabnet_lambda_sparse", tabnet_lambda_sparse))
#'     )
#'
#'     # NODE
#'     node_par <- list(
#'       node_trees = as.integer(gget(hp, "node_trees", node_trees)),
#'       node_depth = as.integer(gget(hp, "node_depth", node_depth))
#'     )
#'
#'     # DeepFM
#'     deepfm_par <- list(
#'       deepfm_k      = as.integer(gget(hp, "deepfm_k", deepfm_k)),
#'       deepfm_hidden = parse_int_vec(gget(hp, "deepfm_hidden", deepfm_hidden))
#'     )
#'
#'     # DCN
#'     dcn_par <- list(
#'       dcn_layers = as.integer(gget(hp, "dcn_layers", dcn_layers)),
#'       dcn_hidden = parse_int_vec(gget(hp, "dcn_hidden", dcn_hidden))
#'     )
#'
#'     # Build argument list for Python
#'     args <- list(
#'       X = X, y = y,
#'       model_type = model_type,
#'       learning_rate = lr,
#'       epochs = ep,
#'       batch_size = bs,
#'       l2_weight_decay = l2w,
#'       dropout = dr,
#'       optimizer_name = "adam",
#'       final_attention = fin_att,
#'       attention_across_multiple_layers = across,
#'       batch_norm = bn,
#'       validation_split = as.numeric(validation_split),
#'       compile_model = cmp,
#'       device = if (is.null(device)) NULL else as.character(device),
#'       deterministic = det,
#'       random_seed = rse,
#'       heteroscedastic = hetero,
#'
#'       # CNN
#'       kernel_size = ks,
#'       dense_layers_cnn = dlcn,
#'       cnn_use_max_pool = mp_use,
#'       cnn_pool_kernel  = mp_kernel,
#'       cnn_pool_stride  = mp_stride,
#'       cnn_pool_padding = mp_padding
#'     )
#'
#'     # Only send hidden-layer specs for models that use them
#'     if (uses_hidden) {
#'       if (!is.null(nlayers)) args$num_hidden_layers <- nlayers
#'       if (!is.null(neurons)) args$neurons_per_layer <- as.integer(neurons)
#'     }
#'
#'     # Attach family-specific params
#'     args <- c(args, ft_par, saint_par, tabnet_par, node_par, deepfm_par, dcn_par)
#'
#'     # ---- call Python ----
#'     fit <- do.call(torch_fit_model, args)
#'
#'     results[[paste0("Model_", i)]] <- list(
#'       model   = fit[[1L]],
#'       history = fit[[2L]],
#'       config  = args
#'     )
#'   }
#'
#'   results
#' }








#'
#'
#' # Function to validate the number of hidden layers against the neurons_per_layer list
#' validate_layers <- function(num_hidden_layers, neurons_per_layer) {
#'   msg <- "\n==================================================\n"
#'   if(is.null(num_hidden_layers) || is.null(neurons_per_layer)){
#'     if(length(num_hidden_layers)==1 && num_hidden_layers>1) num_hidden_layers <- seq(1:num_hidden_layers)
#'   if (length(num_hidden_layers) != length(neurons_per_layer)) {
#'     stop(paste(msg, "The length of num_hidden_layers must match the length of neurons_per_layer"), call. = FALSE)
#'   }
#'   }
#'
#' if(length(num_hidden_layers)>1){
#'   for (i in seq_along(num_hidden_layers)) {
#'     if (num_hidden_layers[i] != length(neurons_per_layer[[i]])) {
#'       stop(paste(msg, paste("Mismatch in hidden layers and neurons per layer at", num_hidden_layers[i])), call. = FALSE)
#'     }
#'   }
#' }
#'   return(TRUE)
#' }
#'
#' #' Title
#' #'
#' #' @param X_train
#' #' @param y_train
#' #' @param param_grid
#' #' @param epochs
#' #' @param batch_size
#' #' @param validation_split
#' #' @param early_stop
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' grid_search_deep_learning <- function(X_train,
#'                                       y_train,
#'                                       param_grid,
#'                                       epochs,
#'                                       batch_size,
#'                                       validation_split = 0.2,
#'                                       learning_rate,
#'                                       dropout_rate = 0.5,
#'                                       kernel_size = 3,
#'                                       dense_layers_cnn,
#'                                       l2_regularizer_dp,
#'                                       n_blocks = 2,
#'                                       n_neurons_per_block,
#'                                       deep_learning_model = "mlp_with_attention",
#'                                       attention_on_final_layer = TRUE,
#'                                       attention_across_multiple_layers = FALSE,
#'                                       batch_normalization = TRUE,
#'                                       para_tunning = TRUE,
#'                                       early_stop = TRUE) {
#'   results <- list()
#'   np <- reticulate::import("numpy")
#'   msg <- "\n==================================================\n"
#'
#'   if(deep_learning_model == "ResNet"){
#'     if(!is.null(n_blocks) && !is.null(n_neurons_per_block)){
#'     if(n_blocks==length(n_neurons_per_block)){
#'       stop(paste(msg, paste("Mismatch in hidden number of block and neurons per layer at", n_neurons_per_block)), call. = FALSE)
#'     }
#'     } else {
#'       stop(paste(msg,"n_neurons_per_block and n_block can't be NULL."), call. = FALSE)
#' }
#'     neurons_per_layer <- n_neurons_per_block
#'     num_hidden_layers <- n_blocks
#'     if(length(num_hidden_layers)==1 && num_hidden_layers>1) num_hidden_layers <- seq(1:num_hidden_layers)
#'
#'     if(!is.null(param_grid$n_blocks)){
#'       if(length(param_grid$n_blocks)==1 && param_grid$n_blocks>1) param_grid$n_blocks <- seq(1:param_grid$n_blocks)
#'       names(param_grid)[names(param_grid)%in%"n_blocks"] <- "num_hidden_layers"
#'     }
#'     if(!is.null(param_grid$n_neurons_per_block)){
#'       names(param_grid)[names(param_grid)%in%"n_neurons_per_block"] <- "neurons_per_layer"
#'     }
#'   }
#'
#'   X_train <- np$array(as.matrix(X_train))
#'   y_train <- np$array(y_train)
#'   # Check if validation_split is within valid range
#'   if (early_stop && (validation_split < 0 || validation_split >= 1)) {
#'     stop(paste(msg, "Validation split must be between 0 and 1."), call. = FALSE)
#'   }
#'
#'   # Call the validation function
#'  #test_neuron_hidden_layer <-  tryCatch({
#'     validate_layers(param_grid$num_hidden_layers, param_grid$neurons_per_layer)
#'     if(!is.null(param_grid$num_hidden_layers)){
#'       if(length(param_grid$num_hidden_layers)==1 && param_grid$num_hidden_layers>1) param_grid$num_hidden_layers <- seq(1:param_grid$num_hidden_layers)
#'     }
#'   #   cat("Validation successful: The number of hidden layers matches the neurons per layer configuration.\n")
#'   # }, error = function(e) {
#'   #   cat("Validation error:", e$message, "\n")
#'   #   return(NULL)
#'   # })
#'
#'  # if(is.null(test_neuron_hidden_layer))
#'   # Generate all combinations of hyperparameters
#'   hyperparam_combinations <- expand.grid(param_grid)
#'
#'   # Function to check if the number of hidden layers matches the length of neurons_per_layer
#'   check_hidden_layers <- function(num_hidden_layers, neurons_per_layer) {
#'     neuron_list <- strsplit(as.character(neurons_per_layer), ",")[[1]]
#'     return(length(neuron_list) == num_hidden_layers)
#'   }
#'
#'   # Loop through the data frame and remove rows with mismatched hidden layers and neurons
#'   valid_rows <- apply(hyperparam_combinations, 1, function(row) {
#'     num_hidden_layers <- as.integer(row["num_hidden_layers"])
#'     neurons_per_layer <- row["neurons_per_layer"]
#'     check_hidden_layers(num_hidden_layers, neurons_per_layer)
#'   })
#'
#'   # Filter the data frame to keep only valid rows
#'   hyperparam_combinations <- hyperparam_combinations[valid_rows, ]
#'
#'   for (i in 1:nrow(hyperparam_combinations)) {
#'     # Extract hyperparameter values for the current combination
#'     hyperparams <- hyperparam_combinations[i, ]
#'
#'     # Initialize variables only if they exist in the hyperparams
#'     if ("num_hidden_layers" %in% names(hyperparams)) {
#'       num_hidden_layers <- as.integer(hyperparams$num_hidden_layers)
#'     }
#'
#'     if ("neurons_per_layer" %in% names(hyperparams)) {
#'       if (inherits(hyperparams$neurons_per_layer, "list")) {
#'         neurons_per_layer <- unlist(hyperparams$neurons_per_layer)
#'       } else {
#'         neurons_per_layer <- hyperparams$neurons_per_layer
#'       }
#'     }
#'
#'     if ("learning_rate" %in% names(hyperparams)) {
#'       learning_rate <- hyperparams$learning_rate
#'     }
#'
#'     if ("dropout_rate" %in% names(hyperparams)) {
#'       dropout_rate <- hyperparams$dropout_rate
#'     }
#'
#'     if ("epochs" %in% names(hyperparams)) {
#'       epochs <- hyperparams$epochs
#'     }
#'
#'     if ("l2_regularizer_dp" %in% names(hyperparams)) {
#'       l2_regularizer_dp <- hyperparams$l2_regularizer_dp
#'     }
#'
#'     if(!is.null(num_hidden_layers) && !is.list(num_hidden_layers)) num_hidden_layers <- as.integer(num_hidden_layers)
#'     if(!is.null(neurons_per_layer) && !is.list(neurons_per_layer)) neurons_per_layer <- as.integer(neurons_per_layer)
#'     if(!is.null(batch_size))   batch_size <- as.integer(batch_size)
#'     if(!is.null(epochs)) epochs <- as.integer(epochs)
#'     if(!is.null(dropout_rate)) dropout_rate <- as.numeric(dropout_rate)
#'     if(!is.null(learning_rate)) learning_rate <- as.numeric(learning_rate)
#'
#'        # Ensure that the number of neurons_per_layer matches the number of num_hidden_layers
#'     if (length(neurons_per_layer) != num_hidden_layers) {
#'       warning(paste(msg, "Length of neurons_per_layer does not match num_hidden_layers. Skipping this combination."), call. = FALSE)
#'       next
#'     }
#'
#'     # Create and compile the model with current hyperparameters
#'     model <- deep_learning_model_utilityy(X_train = X_train,
#'                                          y_train = y_train,
#'                                          num_hidden_layers = num_hidden_layers,
#'                                          neurons_per_layer = neurons_per_layer,
#'                                          learning_rate = learning_rate,
#'                                          epochs = epochs,
#'                                          batch_size = batch_size,
#'                                          dropout_rate = dropout_rate,
#'                                          l2_regularizer_dp = l2_regularizer_dp,
#'                                          validation_split = validation_split,
#'                                          n_blocks = n_blocks,
#'                                          n_neurons_per_block = n_neurons_per_block,
#'                                          deep_learning_model = deep_learning_model,
#'                                          attention_on_final_layer = attention_on_final_layer,
#'                                          attention_across_multiple_layers = attention_across_multiple_layers,
#'                                          batch_normalization = batch_normalization,
#'                                          para_tunning = FALSE)
#'
#'        # Train the model with or without early stopping
#'     if (early_stop) {
#'       # Create a validation set from the training data
#'       set.seed(123)
#'       indices <- sample.int(nrow(X_train), size = floor(validation_split * nrow(X_train)))
#'       X_val <- as.matrix(X_train[indices, ])
#'       y_val <- y_train[indices]
#'       X_trainn <- X_train[-indices, ]
#'       y_trainn <- y_train[-indices]
#'
#'       # Define early stopping callback
#'       early_stopping <- keras::callback_early_stopping(monitor = "val_loss", patience = 3)
#'
#'
#'       # Train the model with early stopping
#'       history <- model$fit(
#'         x = X_trainn,
#'         y = y_trainn,
#'         epochs = epochs,
#'         batch_size = batch_size,
#'         validation_data = list(X_val, y_val),
#'         callbacks = list(early_stopping),
#'         verbose = 0
#'       )
#'     } else {
#'       # Train the model without early stopping
#'
#'       history <- model$fit(
#'         x = X_trainn,
#'         y = y_trainn,
#'         epochs = epochs,
#'         batch_size = batch_size,
#'         verbose = 0
#'       )
#'
#'     }
#'
#'     # Store the model and training history with a systematic name
#'     result_name <- paste0("Model_", i)  # You can adjust the naming scheme as needed
#'     results[[result_name]] <- list(
#'       model = model,
#'       history = history
#'     )
#'
#'     rm(model)
#'   }
#'
#'   return(results)
#' }
