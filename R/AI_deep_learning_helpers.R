# =========================================================
# Utilities
# =========================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
compact <- function(x) x[!vapply(x, is.null, logical(1))]

coerce_int_vec <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
  if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
  as.integer(x)
}
coerce_num <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
  if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
  as.numeric(x)
}
coerce_logi <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
  if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
  isTRUE(x)
}
assert_matrix <- function(X, name = "X") {
  if (!is.matrix(X)) X <- try(as.matrix(X))
  if (!is.matrix(X)) stop(name, " must be coercible to a numeric matrix")
  storage.mode(X) <- "double"
  X
}

# Infer task from y
infer_task <- function(y) {
  yv <- if (is.factor(y)) as.integer(y) else y
  if (is.numeric(yv) && length(unique(na.omit(yv))) > 10) "regression" else "classification"
}
# Encode y for classification
encode_y <- function(y, task) {
  if (task == "classification") {
    if (is.character(y)) y <- factor(y)
    if (is.factor(y)) return(as.integer(y)) # 1..K
    lv <- sort(unique(y))
    match(y, lv)
  } else {
    as.numeric(y)
  }
}

# Map old->canonical commons and coerce
normalize_common <- function(args) {
  canon <- list(
    learning_rate   = NULL, epochs = NULL, batch_size = NULL, l2_weight_decay = NULL,
    dropout = NULL, batch_norm = NULL, validation_split = NULL, compile_model = NULL,
    deterministic = NULL, random_seed = NULL, device = NULL, use_amp = NULL,
    max_grad_norm = NULL, optimizer_name = NULL, heteroscedastic = NULL,
    auto_class_weights = NULL
  )
  if (!is.null(args$l2_regularizer_dp) && is.null(args$l2_weight_decay))
    args$l2_weight_decay <- args$l2_regularizer_dp
  if (!is.null(args$dropout_rate) && is.null(args$dropout))
    args$dropout <- args$dropout_rate
  if (!is.null(args$batch_normalization) && is.null(args$batch_norm))
    args$batch_norm <- args$batch_normalization

  canon$learning_rate    <- coerce_num(args$learning_rate %||% 1e-3)
  canon$epochs           <- coerce_int_vec(args$epochs %||% 30L)
  canon$batch_size       <- coerce_int_vec(args$batch_size %||% 64L)
  canon$l2_weight_decay  <- coerce_num(args$l2_weight_decay %||% 1e-3)
  canon$dropout          <- coerce_num(args$dropout %||% 0.5)
  canon$batch_norm       <- coerce_logi(args$batch_norm %||% TRUE)
  canon$validation_split <- coerce_num(args$validation_split %||% 0.2)
  canon$compile_model    <- coerce_logi(args$compile_model %||% TRUE)
  canon$deterministic    <- coerce_logi(args$deterministic %||% FALSE)
  canon$random_seed      <- if (is.null(args$random_seed)) NULL else coerce_int_vec(args$random_seed)
  canon$device           <- if (is.null(args$device)) NULL else as.character(args$device)
  canon$use_amp          <- coerce_logi(args$use_amp %||% TRUE)
  canon$max_grad_norm    <- if (is.null(args$max_grad_norm)) NULL else coerce_num(args$max_grad_norm)
  canon$optimizer_name   <- as.character(args$optimizer_name %||% "adam")
  canon$heteroscedastic  <- coerce_logi(args$heteroscedastic %||% FALSE)
  canon$auto_class_weights <- coerce_logi(args$auto_class_weights %||% FALSE)
  canon
}

# Aliases
apply_aliases <- function(dots, aliases) {
  for (src in names(aliases)) {
    dst <- aliases[[src]]
    if (!is.null(dots[[src]]) && is.null(dots[[dst]])) dots[[dst]] <- dots[[src]]
    dots[[src]] <- NULL
  }
  dots
}

# =========================================================
# GP kernel canonicalizer
# =========================================================
GP_KERNEL_CHOICES <- c(
  "rbf","matern12","matern32","matern52","rq","linear","polynomial","cosine","periodic","spectral_mixture"
)
GP_KERNEL_ALIASES <- c("se"="rbf","sqexp"="rbf","exp_quad"="rbf","cos"="cosine","matern"="matern32")
canonicalize_gp_kernel <- function(k) {
  if (is.null(k) || length(k) == 0) return("rbf")
  k <- tolower(trimws(as.character(k)))
  if (k %in% names(GP_KERNEL_ALIASES)) k <- unname(GP_KERNEL_ALIASES[k])
  if (!(k %in% GP_KERNEL_CHOICES)) {
    stop(sprintf("Unknown gp_kernel='%s'. Allowed: %s", k, paste(GP_KERNEL_CHOICES, collapse=", ")), call. = FALSE)
  }
  k
}

# =========================================================
# Alias maps (prefixed -> canonical)
# =========================================================
aliases_common_prefixed <- function(prefix) {
  stats::setNames(
    c("learning_rate","epochs","batch_size","l2_weight_decay","l2_weight_decay",
      "dropout","dropout","batch_norm","validation_split","compile_model",
      "deterministic","random_seed","device","use_amp","max_grad_norm",
      "optimizer_name","heteroscedastic","auto_class_weights"),
    c(paste0(prefix,"_learning_rate"),
      paste0(prefix,"_epochs"),
      paste0(prefix,"_batch_size"),
      paste0(prefix,"_l2_weight_decay"),
      paste0(prefix,"_l2_regularizer_dp"),
      paste0(prefix,"_dropout"),
      paste0(prefix,"_dropout_rate"),
      paste0(prefix,"_batch_norm"),
      paste0(prefix,"_validation_split"),
      paste0(prefix,"_compile_model"),
      paste0(prefix,"_deterministic"),
      paste0(prefix,"_random_seed"),
      paste0(prefix,"_device"),
      paste0(prefix,"_use_amp"),
      paste0(prefix,"_max_grad_norm"),
      paste0(prefix,"_optimizer_name"),
      paste0(prefix,"_heteroscedastic"),
      paste0(prefix,"_auto_class_weights"))
  )
}
aliases_cnn <- function() c(
  "cnn_neurons_per_layer"="neurons_per_layer","cnn_kernel_size"="kernel_size","cnn_dense_layers"="dense_layers_cnn",
  "cnn_use_max_pool"="cnn_use_max_pool","use_max_pool"="cnn_use_max_pool",
  "cnn_pool_kernel"="cnn_pool_kernel","cnn_pool_stride"="cnn_pool_stride","cnn_pool_padding"="cnn_pool_padding",
  "cnn_separable"="separable","cnn_dilations"="dilations","cnn_use_se"="use_se","cnn_norm_type"="norm_type",
  "cnn_pool_type"="pool_type","cnn_use_global_pool"="use_global_pool",
  aliases_common_prefixed("cnn")
)
aliases_mlp <- function() c(
  "mlp_neurons_per_layer"="neurons_per_layer","mlp_final_attention"="final_attention",
  "mlp_attention_across_layers"="attention_across_multiple_layers",
  aliases_common_prefixed("mlp")
)
aliases_ft <- function() c(
  "ft_d_model"="ft_d_model","ft_heads"="ft_heads","ft_layers"="ft_layers","ft_ff_mult"="ft_ff_mult",
  "ft_dropout"="ft_dropout","ft_token_dropout"="ft_token_dropout","ft_use_cls"="ft_use_cls",
  "ft_scalar_tokenizer"="ft_scalar_tokenizer",
  "ft_use_grouping"="use_grouping","ft_group_trigger"="group_trigger","ft_group_method"="group_method",
  "ft_init_group_size"="init_group_size","ft_max_tokens"="max_tokens","ft_kmeans_batch"="kmeans_batch","ft_kmeans_iter"="kmeans_iter",
  aliases_common_prefixed("ft")
)
aliases_resnet <- function() c("resnet_neurons_per_block"="neurons_per_layer","resnet_blocks"="num_hidden_layers",
                                aliases_common_prefixed("resnet"))
aliases_saint <- function() c(
  "saint_d_model"="saint_d_model","saint_heads"="saint_heads","saint_layers"="saint_layers","saint_ff_mult"="saint_ff_mult",
  "saint_dropout"="saint_dropout","saint_token_dropout"="saint_token_dropout","saint_use_cls"="saint_use_cls",
  "saint_use_grouping"="use_grouping","saint_group_trigger"="group_trigger","saint_group_method"="group_method",
  "saint_init_group_size"="init_group_size","saint_max_tokens"="max_tokens","saint_kmeans_batch"="kmeans_batch","saint_kmeans_iter"="kmeans_iter",
  aliases_common_prefixed("saint")
)
aliases_tabnet <- function() c(aliases_common_prefixed("tabnet"))
aliases_node   <- function() c(aliases_common_prefixed("node"))
aliases_deepfm <- function() c(aliases_common_prefixed("deepfm"))
aliases_dcnv2  <- function() c(aliases_common_prefixed("dcn"))
aliases_nam <- function() c(
  "nam_hidden"="nam_hidden","nam_activation"="nam_activation","nam_add_linear"="nam_add_linear","nam_l1"="nam_l1",
  "nam_use_grouping"="use_grouping","nam_group_trigger"="group_trigger","nam_group_method"="group_method",
  "nam_init_group_size"="init_group_size","nam_max_tokens"="max_tokens","nam_kmeans_batch"="kmeans_batch","nam_kmeans_iter"="kmeans_iter",
  aliases_common_prefixed("nam")
)
aliases_moe <- function() c(
  "moe_n_experts"="moe_n_experts","moe_expert_hidden"="moe_expert_hidden","moe_gate_hidden"="moe_gate_hidden",
  "moe_temperature"="moe_temperature","moe_sparse_topk"="moe_sparse_topk","moe_entropy_reg"="moe_entropy_reg",
  "moe_use_grouping"="use_grouping","moe_group_trigger"="group_trigger","moe_group_method"="group_method",
  "moe_init_group_size"="init_group_size","moe_max_tokens"="max_tokens","moe_kmeans_batch"="kmeans_batch","moe_kmeans_iter"="kmeans_iter",
  aliases_common_prefixed("moe")
)
aliases_gp_dkl <- function() c(
  "gp_use_variational"="gp_use_variational","gp_num_inducing"="gp_num_inducing","gp_feature_dim"="gp_feature_dim",
  "gp_kernel"="gp_kernel","gp_ard"="gp_ard","gp_lr_mult"="gp_lr_mult",
  "gp_rff_features"="rff_features","gp_rff_lengthscale"="rff_lengthscale","gp_rff_deep_hidden"="rff_deep_hidden",
  aliases_common_prefixed("gp")
)

# =========================================================
# Safe torch fit call + shim
# =========================================================
torch_fit_filtered <- function(args) {
  sig <- try(names(formals(torch_fit_model)), silent = TRUE)
  if (!inherits(sig, "try-error") && length(sig)) args <- args[intersect(names(args), sig)]
  do.call(torch_fit_model, args)
}

torch_fit_model <- function(X, y, ...) {
  mod  <- get_dl_module()
  #mod <- init_dp_module(prefer_gpu = TRUE)
  dots <- list(...)

  # Duplicate kwargs guard
  dups <- names(dots)[duplicated(names(dots))]
  if (length(dups)) {
    stop(
      "Duplicate kwargs passed to torch_fit_model(): ",
      paste(unique(dups), collapse = ", "),
      call. = FALSE
    )
  }

  # Legacy → current renames
  if (!is.null(dots$deep_learning_model) && is.null(dots$model_type)) {
    mt <- tolower(as.character(dots$deep_learning_model))
    dots$model_type <- switch(mt, "ft"="ft_transformer","dcn"="dcnv2", mt)
    dots$deep_learning_model <- NULL
  }
  if (!is.null(dots$attention_on_final_layer) && is.null(dots$final_attention)) {
    dots$final_attention <- isTRUE(dots$attention_on_final_layer); dots$attention_on_final_layer <- NULL
  }
  if (!is.null(dots$batch_normalization) && is.null(dots$batch_norm)) {
    dots$batch_norm <- isTRUE(dots$batch_normalization); dots$batch_normalization <- NULL
  }
  if (!is.null(dots$l2_regularizer_dp) && is.null(dots$l2_weight_decay)) {
    dots$l2_weight_decay <- as.numeric(dots$l2_regularizer_dp); dots$l2_regularizer_dp <- NULL
  }
  if (!is.null(dots$dropout_rate) && is.null(dots$dropout)) {
    dots$dropout <- as.numeric(dots$dropout_rate); dots$dropout_rate <- NULL
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

  args <- c(list(X = X, y = y), dots)
  do.call(mod$fit_model, args)  # expects (model, history)
}

# =========================================================
# Canonizers + wrappers
# =========================================================
canon_cnn_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_cnn()); C <- normalize_common(dots)
  list(
    neurons_per_layer = coerce_int_vec(dots$neurons_per_layer %||% c(64L,64L,64L)),
    kernel_size       = coerce_int_vec(dots$kernel_size %||% 3L),
    dense_layers_cnn  = coerce_int_vec(dots$dense_layers_cnn %||% c(256L,128L,64L)),
    cnn_use_max_pool  = coerce_logi(dots$cnn_use_max_pool %||% FALSE),
    cnn_pool_kernel   = coerce_int_vec(dots$cnn_pool_kernel %||% 2L),
    cnn_pool_stride   = coerce_int_vec(dots$cnn_pool_stride %||% 2L),
    cnn_pool_padding  = coerce_int_vec(dots$cnn_pool_padding %||% 0L),
    separable         = coerce_logi(dots$separable %||% TRUE),
    dilations         = coerce_int_vec(dots$dilations %||% c(1L,2L,4L)),
    use_se            = coerce_logi(dots$use_se %||% TRUE),
    norm_type         = as.character(dots$norm_type %||% "group"),
    pool_type         = as.character(dots$pool_type %||% "conv"),
    use_global_pool   = coerce_logi(dots$use_global_pool %||% FALSE),
    # common
    learning_rate    = C$learning_rate, epochs = C$epochs, batch_size = C$batch_size,
    l2_weight_decay  = C$l2_weight_decay, dropout = C$dropout, batch_norm = C$batch_norm,
    validation_split = C$validation_split, compile_model = C$compile_model,
    deterministic = C$deterministic, random_seed = C$random_seed, device = C$device,
    use_amp = C$use_amp, max_grad_norm = C$max_grad_norm, optimizer_name = C$optimizer_name,
    heteroscedastic = C$heteroscedastic, auto_class_weights = C$auto_class_weights
  )
}
dl_cnn <- function(X, y, ..., model_type = "cnn") {
  cfg <- canon_cnn_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="cnn"), cfg))
}

canon_mlp_args <- function(dots, with_attention = FALSE) {
  dots <- apply_aliases(dots, aliases_mlp()); C <- normalize_common(dots)
  list(
    neurons_per_layer = coerce_int_vec(dots$neurons_per_layer %||% c(256L,128L)),
    final_attention = coerce_logi(dots$final_attention %||% with_attention),
    attention_across_multiple_layers = coerce_logi(dots$attention_across_multiple_layers %||% FALSE),
    # common
    learning_rate = C$learning_rate, epochs = C$epochs, batch_size = C$batch_size,
    l2_weight_decay = C$l2_weight_decay, dropout = C$dropout, batch_norm = C$batch_norm,
    validation_split = C$validation_split, compile_model = C$compile_model,
    deterministic = C$deterministic, random_seed = C$random_seed, device = C$device,
    use_amp = C$use_amp, max_grad_norm = C$max_grad_norm, optimizer_name = C$optimizer_name,
    heteroscedastic = C$heteroscedastic, auto_class_weights = C$auto_class_weights
  )
}
dl_mlp <- function(X, y, ..., model_type = c("mlp","mlp_with_attention")) {
  model_type <- match.arg(model_type); cfg <- canon_mlp_args(list(...), with_attention = (model_type=="mlp_with_attention"))
  torch_fit_filtered(c(list(X=X, y=y, model_type=model_type), cfg))
}

canon_ft_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_ft()); C <- normalize_common(dots)
  list(
    ft_d_model = coerce_int_vec(dots$ft_d_model %||% 192L),
    ft_heads = coerce_int_vec(dots$ft_heads %||% 8L),
    ft_layers= coerce_int_vec(dots$ft_layers %||% 4L),
    ft_ff_mult = coerce_int_vec(dots$ft_ff_mult %||% 4L),
    ft_dropout = coerce_num(dots$ft_dropout %||% 0.1),
    ft_token_dropout = coerce_num(dots$ft_token_dropout %||% 0.1),
    ft_use_cls = coerce_logi(dots$ft_use_cls %||% TRUE),
    ft_scalar_tokenizer = coerce_logi(dots$ft_scalar_tokenizer %||% TRUE),
    use_grouping = coerce_logi(dots$use_grouping %||% TRUE),
    group_trigger = coerce_int_vec(dots$group_trigger %||% 2048L),
    group_method = as.character(dots$group_method %||% "auto"),
    init_group_size = coerce_int_vec(dots$init_group_size %||% 64L),
    max_tokens = coerce_int_vec(dots$max_tokens %||% 1024L),
    kmeans_batch = coerce_int_vec(dots$kmeans_batch %||% 4096L),
    kmeans_iter  = coerce_int_vec(dots$kmeans_iter %||% 100L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=C$batch_size,
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_ft_transformer <- function(X, y, ..., model_type = "ft_transformer") {
  cfg <- canon_ft_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="ft_transformer"), cfg))
}

canon_resnet_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_resnet()); C <- normalize_common(dots)
  npl <- coerce_int_vec(dots$neurons_per_layer %||% c(64L,64L,64L))
  n_layers <- coerce_int_vec(dots$num_hidden_layers %||% length(npl))
  list(
    neurons_per_layer = npl, num_hidden_layers = n_layers,
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=C$batch_size,
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_resnet <- function(X, y, ..., model_type = "resnet") {
  cfg <- canon_resnet_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="resnet"), cfg))
}

canon_saint_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_saint()); C <- normalize_common(dots)
  list(
    saint_d_model=coerce_int_vec(dots$saint_d_model %||% 128L),
    saint_heads=coerce_int_vec(dots$saint_heads %||% 8L),
    saint_layers=coerce_int_vec(dots$saint_layers %||% 4L),
    saint_ff_mult=coerce_int_vec(dots$saint_ff_mult %||% 4L),
    saint_dropout=coerce_num(dots$saint_dropout %||% 0.1),
    saint_token_dropout=coerce_num(dots$saint_token_dropout %||% 0.1),
    saint_use_cls=coerce_logi(dots$saint_use_cls %||% TRUE),
    use_grouping=coerce_logi(dots$use_grouping %||% TRUE),
    group_trigger=coerce_int_vec(dots$group_trigger %||% 2048L),
    group_method=as.character(dots$group_method %||% "auto"),
    init_group_size=coerce_int_vec(dots$init_group_size %||% 64L),
    max_tokens=coerce_int_vec(dots$max_tokens %||% 1024L),
    kmeans_batch=coerce_int_vec(dots$kmeans_batch %||% 4096L),
    kmeans_iter=coerce_int_vec(dots$kmeans_iter %||% 100L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=C$batch_size,
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_saint <- function(X, y, ..., model_type = "saint") {
  cfg <- canon_saint_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="saint"), cfg))
}

canon_tabnet_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_tabnet()); C <- normalize_common(dots)
  list(
    tabnet_steps=coerce_int_vec(dots$tabnet_steps %||% 5L),
    tabnet_feature_dim=coerce_int_vec(dots$tabnet_feature_dim %||% 64L),
    tabnet_output_dim=coerce_int_vec(dots$tabnet_output_dim %||% 64L),
    tabnet_gamma=coerce_num(dots$tabnet_gamma %||% 1.5),
    tabnet_lambda_sparse=coerce_num(dots$tabnet_lambda_sparse %||% 1e-4),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 512L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_tabnet <- function(X, y, ..., model_type = "tabnet") {
  cfg <- canon_tabnet_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="tabnet"), cfg))
}

canon_node_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_node()); C <- normalize_common(dots)
  list(
    node_trees=coerce_int_vec(dots$node_trees %||% 8L),
    node_depth=coerce_int_vec(dots$node_depth %||% 3L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 1024L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_node <- function(X, y, ..., model_type = "node") {
  cfg <- canon_node_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="node"), cfg))
}

canon_deepfm_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_deepfm()); C <- normalize_common(dots)
  list(
    deepfm_k=coerce_int_vec(dots$deepfm_k %||% 16L),
    deepfm_hidden=coerce_int_vec(dots$deepfm_hidden %||% c(128L,64L)),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 1024L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_deepfm <- function(X, y, ..., model_type = "deepfm") {
  cfg <- canon_deepfm_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="deepfm"), cfg))
}

canon_dcnv2_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_dcnv2()); C <- normalize_common(dots)
  list(
    dcn_layers=coerce_int_vec(dots$dcn_layers %||% 3L),
    dcn_hidden=coerce_int_vec(dots$dcn_hidden %||% c(256L,128L)),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 1024L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_dcnv2 <- function(X, y, ..., model_type = "dcnv2") {
  cfg <- canon_dcnv2_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="dcnv2"), cfg))
}

canon_nam_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_nam()); C <- normalize_common(dots)
  list(
    nam_hidden=coerce_int_vec(dots$nam_hidden %||% c(32L,16L)),
    nam_activation=as.character(dots$nam_activation %||% "relu"),
    nam_add_linear=coerce_logi(dots$nam_add_linear %||% TRUE),
    nam_l1=coerce_num(dots$nam_l1 %||% 1e-4),
    use_grouping=coerce_logi(dots$use_grouping %||% TRUE),
    group_trigger=coerce_int_vec(dots$group_trigger %||% 2048L),
    group_method=as.character(dots$group_method %||% "auto"),
    init_group_size=coerce_int_vec(dots$init_group_size %||% 64L),
    max_tokens=coerce_int_vec(dots$max_tokens %||% 1024L),
    kmeans_batch=coerce_int_vec(dots$kmeans_batch %||% 4096L),
    kmeans_iter=coerce_int_vec(dots$kmeans_iter %||% 100L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 512L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_nam <- function(X, y, ..., model_type = "nam") {
  cfg <- canon_nam_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="nam"), cfg))
}

canon_moe_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_moe()); C <- normalize_common(dots)
  list(
    moe_n_experts=coerce_int_vec(dots$moe_n_experts %||% 4L),
    moe_expert_hidden=coerce_int_vec(dots$moe_expert_hidden %||% c(128L,64L)),
    moe_gate_hidden=coerce_int_vec(dots$moe_gate_hidden %||% 128L),
    moe_temperature=coerce_num(dots$moe_temperature %||% 1.0),
    moe_sparse_topk= if (is.null(dots$moe_sparse_topk)) NULL else coerce_int_vec(dots$moe_sparse_topk),
    moe_entropy_reg=coerce_num(dots$moe_entropy_reg %||% 0.0),
    use_grouping=coerce_logi(dots$use_grouping %||% TRUE),
    group_trigger=coerce_int_vec(dots$group_trigger %||% 2048L),
    group_method=as.character(dots$group_method %||% "auto"),
    init_group_size=coerce_int_vec(dots$init_group_size %||% 64L),
    max_tokens=coerce_int_vec(dots$max_tokens %||% 1024L),
    kmeans_batch=coerce_int_vec(dots$kmeans_batch %||% 4096L),
    kmeans_iter=coerce_int_vec(dots$kmeans_iter %||% 100L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 512L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_moe <- function(X, y, ..., model_type = "moe") {
  cfg <- canon_moe_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="moe"), cfg))
}

canon_gp_dkl_args <- function(dots) {
  dots <- apply_aliases(dots, aliases_gp_dkl()); C <- normalize_common(dots)
  list(
    gp_use_variational=coerce_logi(dots$gp_use_variational %||% TRUE),
    gp_num_inducing=coerce_int_vec(dots$gp_num_inducing %||% 512L),
    gp_feature_dim=coerce_int_vec(dots$gp_feature_dim %||% 64L),
    gp_kernel=canonicalize_gp_kernel(dots$gp_kernel %||% "rbf"),
    gp_ard=coerce_logi(dots$gp_ard %||% TRUE),
    gp_lr_mult=coerce_num(dots$gp_lr_mult %||% 0.5),
    rff_features=coerce_int_vec(dots$rff_features %||% 1024L),
    rff_lengthscale=coerce_num(dots$rff_lengthscale %||% 1.0),
    rff_deep_hidden=coerce_int_vec(dots$rff_deep_hidden %||% 128L),
    # common
    learning_rate=C$learning_rate, epochs=C$epochs, batch_size=coerce_int_vec(dots$batch_size %||% 512L),
    l2_weight_decay=C$l2_weight_decay, dropout=C$dropout, batch_norm=C$batch_norm,
    validation_split=C$validation_split, compile_model=C$compile_model,
    deterministic=C$deterministic, random_seed=C$random_seed, device=C$device,
    use_amp=C$use_amp, max_grad_norm=C$max_grad_norm, optimizer_name=C$optimizer_name,
    heteroscedastic=C$heteroscedastic, auto_class_weights=C$auto_class_weights
  )
}
dl_gp_dkl <- function(X, y, ..., model_type = "gp_dkl") {
  cfg <- canon_gp_dkl_args(list(...))
  torch_fit_filtered(c(list(X=X, y=y, model_type="gp_dkl"), cfg))
}

# =========================================================
# Registry + coordinator
# =========================================================
dl_registry <- new.env(parent = emptyenv())
register_model <- function(key, fn, allowed_args = NULL) {
  dl_registry[[key]] <- list(fn = fn, allowed = allowed_args)
}
get_model <- function(key) dl_registry[[key]]

register_model("cnn", dl_cnn)
register_model("mlp", function(...) dl_mlp(..., model_type = "mlp"))
register_model("mlp_with_attention", function(...) dl_mlp(..., model_type = "mlp_with_attention"))
register_model("ft_transformer", dl_ft_transformer)
register_model("resnet", dl_resnet)
register_model("saint", dl_saint)
register_model("tabnet", dl_tabnet)
register_model("node", dl_node)
register_model("deepfm", dl_deepfm)
register_model("dcnv2", dl_dcnv2)
register_model("nam", dl_nam)
register_model("moe", dl_moe)
register_model("gp_dkl", dl_gp_dkl)

dl_fit <- function(X, y, model = c("cnn","mlp","mlp_with_attention","ft_transformer",
                                   "resnet","saint","tabnet","node","deepfm","dcnv2",
                                   "nam","moe","gp_dkl"), ...) {
  model <- match.arg(model)
  dots  <- list(...)
  X <- assert_matrix(X)
  task <- infer_task(y)
  y_enc <- encode_y(y, task)
  dots$task <- task
  reg <- get_model(model)
  if (is.null(reg)) stop("Model '", model, "' is not registered yet.")
  do.call(reg$fn, c(list(X = X, y = y_enc), dots))
}

# =========================================================
# Boot statistic: train + predict using dl_fit
# =========================================================
train_predict_deeplearning <- function(data, indices, model_type, dl_args, test_geno = NULL) {
  # data: cbind(y_scaled, X); indices: bootstrap sample
  stopifnot(ncol(data) >= 2)
  X_all <- as.matrix(data[, -1, drop = FALSE])
  y_all <- as.numeric(data[, 1])

  X_tr <- X_all[indices, , drop = FALSE]
  y_tr <- y_all[indices]

  # Fit
  fit <- tryCatch({
    do.call(dl_fit, c(list(X = X_tr, y = y_tr, model = tolower(model_type)), dl_args))
  }, error = function(e) stop("dl_fit failed inside bootstrap: ", e$message, call. = FALSE))

  # Predict target matrix
  mod <- get_dl_module()
  #mod <- init_dp_module(prefer_gpu = TRUE)
  X_targ <- if (is.null(test_geno)) X_all else as.matrix(test_geno)
  preds <- mod$predict(fit[[1L]], X_targ,
                       device = if (!is.null(dl_args$device)) as.character(dl_args$device) else NULL)
  as.numeric(preds)
}

# =========================================================
# Main coordinator with boot + CV options
# =========================================================
deep_learning_model <- function(pheno_object = NULL,
                                y = NULL,
                                omics_data = NULL,
                                crossval = FALSE,
                                tst = NULL,
                                geno_omic_object = NULL,
                                geno_omic_test_object = NULL,
                                response = NULL,
                                gen_name = NULL,
                                message = TRUE,
                                scaling = TRUE,
                                centering = FALSE,
                                omic_count = NULL,
                                para_tunning = FALSE,
                                param_grid = NULL,
                                early_stop = TRUE,
                                model_type = "mlp_with_attention",
                                CI_width_thresholds = c(0.33, 0.66),
                                high_reliability_thres = 0.9,
                                low_reliability_thres = 0.5,
                                n_components = 20,
                                threshold = 100,
                                target = "test_set",
                                iqr_multiplier = 1.5,
                                interval_width_high_threshold = NULL,
                                interval_width_low_threshold = NULL,
                                interval_width_moderate_threshold = NULL,
                                n_bootstrap = 100,
                                system_database = FALSE,
                                # common knobs
                                optimizer_name = "adam",
                                use_amp        = TRUE,
                                max_grad_norm  = 1.0,
                                auto_class_weights = FALSE,
                                # CNN
                                cnn_neurons_per_layer = as.integer(c(64,64,64)),
                                cnn_kernel_size = 3L,
                                cnn_dense_layers = as.integer(c(256,128,64)),
                                cnn_use_max_pool = TRUE,
                                cnn_pool_kernel = 2L, cnn_pool_stride = 2L, cnn_pool_padding = 0L,
                                cnn_learning_rate = 1e-3, cnn_separable=TRUE, cnn_dilations=c(1L,2L,4L),
                                cnn_l2_regularizer_dp = 1e-4, cnn_dropout_rate = 0.25,
                                cnn_use_se=TRUE, cnn_norm_type="group", cnn_pool_type="conv", cnn_use_global_pool=FALSE,
                                # ResNet
                                resnet_neurons_per_block = as.integer(c(256,128,64)),
                                resnet_blocks = 3L, resnet_learning_rate = 1e-3,
                                # FT
                                ft_d_model = 192L, ft_heads = 8L, ft_layers = 4L, ft_ff_mult = 4L,
                                ft_dropout = 0.1, ft_token_dropout = 0.1, ft_use_cls = TRUE,
                                # SAINT
                                saint_d_model = 128L, saint_heads = 8L, saint_layers = 4L, saint_ff_mult = 4L,
                                saint_dropout = 0.1, saint_token_dropout = 0.1, saint_use_cls = TRUE,
                                # Grouping
                                use_grouping   = FALSE, group_trigger  = 2048, group_method   = "auto",
                                init_group_size = 64, max_tokens = 1024, kmeans_batch = 4096, kmeans_iter = 100,
                                # TabNet
                                tabnet_steps = 5L, tabnet_feature_dim = 64L, tabnet_output_dim = 64L,
                                tabnet_gamma = 1.5, tabnet_lambda_sparse = 1e-4,
                                # NODE
                                node_trees = 8L, node_depth = 3L,
                                # DeepFM
                                deepfm_k = 16L, deepfm_hidden = as.integer(c(128,64)),
                                # DCNv2
                                dcn_layers = 3L, dcn_hidden = as.integer(c(256,128)),
                                # NAM
                                nam_hidden = as.integer(c(32,16)), nam_activation = "relu", nam_add_linear = TRUE, nam_l1 = 1e-4,
                                # MoE
                                moe_n_experts = 4L, moe_expert_hidden = as.integer(c(128,64)), moe_gate_hidden = 128L,
                                moe_temperature = 1.0, moe_sparse_topk = NA, moe_entropy_reg = 0.0,
                                # GP-DKL / RFF
                                gp_use_variational = TRUE, gp_num_inducing = 256L, gp_feature_dim = 64L,
                                gp_kernel = "rbf", gp_ard = TRUE, gp_lr_mult = 0.5,
                                rff_features = 1024, rff_lengthscale = 1.0, rff_deep_hidden = c(128),
                                # generic training
                                epochs = 30, batch_size = 64 , dropout = 0.2,
                                l2_weight_decay = 1e-4, l2_regularizer_dp = 0.001, dropout_rate = 0.5,
                                batch_norm = TRUE, validation_split = 0.2, compile_model = FALSE,
                                deterministic = TRUE, random_seed = 123, device = NULL,
                                # MLP
                                mlp_neurons_per_layer = as.integer(c(128,64)),
                                mlp_learning_rate = 1e-3,
                                final_attention = TRUE,
                                attention_across_multiple_layers = TRUE,
                                heteroscedastic = FALSE,
                                ...) {

  msg <- "\n==================================================\n"
  prefer_gpu <- !identical(device, "cpu")
  ensure_pydeps(prefer_gpu = prefer_gpu, cuda_version = "auto")
  mod <- get_dl_module()
  #mod <- init_dp_module(prefer_gpu = TRUE)
  # ---------- Prepare X/y + scaling ----------
  if (is.null(geno_omic_object) && is.null(pheno_object) && isFALSE(crossval)) {
    stop(paste(msg, "provide matrix of the predictors and the data.frame of the Y variable."), call. = FALSE)
  }

  if (!is.null(geno_omic_object)) {
    GID <- rownames(geno_omic_object)
    scaler <- caret::preProcess(geno_omic_object, method = c("center","scale"))
    geno_omic_object <- stats::predict(scaler, geno_omic_object)
    cols_with_na <- which(colSums(is.na(geno_omic_object)) > 0)
    if (length(cols_with_na)) {
      geno_omic_object <- geno_omic_object[, -cols_with_na, drop = FALSE]
      if (!is.null(geno_omic_test_object)) {
        geno_omic_test_object <- geno_omic_test_object[, -cols_with_na, drop = FALSE]
      }
    }
  }

  if (!is.null(geno_omic_test_object)) {
    test_label <- rownames(geno_omic_test_object)
    geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
    GID <- rownames(geno_omic_test_object)
    scaler <- caret::preProcess(geno_omic_test_object, method = c("center","scale"))
    geno_omic_test_object <- stats::predict(scaler, geno_omic_test_object)
  }

  if (!is.null(omics_data) && isTRUE(crossval)) {
    scaler <- caret::preProcess(omics_data, method = c("center", "scale"))
    omics_data <- stats::predict(scaler, omics_data)
    cols_with_na <- which(colSums(is.na(omics_data)) > 0)
    if (length(cols_with_na)) omics_data <- omics_data[, -cols_with_na, drop = FALSE]

    # ADD THIS
    omics_data <- as.matrix(omics_data)
    storage.mode(omics_data) <- "double"
  }


  if (isTRUE(crossval)) { para_tunning <- FALSE; param_grid <- NULL }

  # y and scaling
  if (isFALSE(crossval)) {
    y_train <- as.numeric(pheno_object[, response])
  } else {
    y_train <- as.numeric(y)
  }
  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center","scale"))
  y_train_scaled <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[,1]

  # Build matrices
  if (isTRUE(crossval)) {
    if (is.null(tst)) stop("crossval=TRUE but 'tst' indices are NULL", call. = FALSE)
    stopifnot(is.matrix(omics_data))
    stopifnot(is.numeric(tst), all(tst >= 1), all(tst <= nrow(omics_data)))

    X_tr <- as.matrix(omics_data[-tst, , drop = FALSE])
    X_te <- as.matrix(omics_data[ tst, , drop = FALSE])
    y_tr <- y_train_scaled[-tst]
    y_te <- y_train_scaled[ tst]
  } else {
    X_tr <- as.matrix(geno_omic_object)
    X_te <- if (!is.null(geno_omic_test_object)) as.matrix(geno_omic_test_object) else NULL
    y_tr <- y_train_scaled
  }

  # ---------- Common knobs (forwarded to dl_fit) ----------
  args_common <- compact(list(
    optimizer_name = optimizer_name, use_amp = use_amp, max_grad_norm = max_grad_norm,
    auto_class_weights = auto_class_weights, compile_model = compile_model,
    deterministic = deterministic, random_seed = random_seed, device = device,
    batch_norm = batch_norm, validation_split = validation_split
  ))

  # ---------- Model-specific knobs ----------
  args_model <- switch(tolower(model_type),
                       "cnn" = compact(list(
                         cnn_neurons_per_layer = cnn_neurons_per_layer,
                         cnn_kernel_size = cnn_kernel_size, cnn_dense_layers = cnn_dense_layers,
                         cnn_use_max_pool = cnn_use_max_pool, cnn_pool_kernel = cnn_pool_kernel,
                         cnn_pool_stride = cnn_pool_stride, cnn_pool_padding = cnn_pool_padding,
                         cnn_separable = cnn_separable, cnn_dilations = cnn_dilations, cnn_use_se = cnn_use_se,
                         cnn_norm_type = cnn_norm_type, cnn_pool_type = cnn_pool_type, cnn_use_global_pool = cnn_use_global_pool,
                         cnn_learning_rate = cnn_learning_rate, cnn_epochs = epochs, cnn_batch_size = batch_size,
                         cnn_l2_regularizer_dp = cnn_l2_regularizer_dp, cnn_dropout_rate = cnn_dropout_rate
                       )),
                       "mlp" = compact(list(
                         mlp_neurons_per_layer = mlp_neurons_per_layer, mlp_final_attention = FALSE,
                         mlp_attention_across_layers = attention_across_multiple_layers,
                         mlp_learning_rate = mlp_learning_rate, mlp_epochs = epochs, mlp_batch_size = batch_size,
                         mlp_dropout_rate = dropout_rate
                       )),
                       "mlp_with_attention" = compact(list(
                         mlp_neurons_per_layer = mlp_neurons_per_layer, mlp_final_attention = TRUE,
                         mlp_attention_across_layers = attention_across_multiple_layers,
                         mlp_learning_rate = mlp_learning_rate, mlp_epochs = epochs, mlp_batch_size = batch_size,
                         mlp_dropout_rate = dropout_rate
                       )),
                       "ft_transformer" = compact(list(
                         ft_d_model = ft_d_model, ft_heads = ft_heads, ft_layers = ft_layers, ft_ff_mult = ft_ff_mult,
                         ft_dropout = ft_dropout, ft_token_dropout = ft_token_dropout, ft_use_cls = ft_use_cls,
                         ft_use_grouping = use_grouping, ft_group_trigger = group_trigger, ft_group_method = group_method,
                         ft_init_group_size = init_group_size, ft_max_tokens = max_tokens, ft_kmeans_batch = kmeans_batch, ft_kmeans_iter = kmeans_iter,
                         ft_epochs = epochs, ft_batch_size = batch_size
                       )),
                       "resnet" = compact(list(
                         resnet_neurons_per_block = resnet_neurons_per_block, resnet_blocks = resnet_blocks,
                         resnet_learning_rate = resnet_learning_rate, resnet_epochs = epochs,
                         resnet_batch_size = batch_size, resnet_dropout_rate = dropout_rate
                       )),
                       "saint" = compact(list(
                         saint_d_model = saint_d_model, saint_heads = saint_heads, saint_layers = saint_layers, saint_ff_mult = saint_ff_mult,
                         saint_dropout = saint_dropout, saint_token_dropout = saint_token_dropout, saint_use_cls = saint_use_cls,
                         saint_use_grouping = use_grouping, saint_group_trigger = group_trigger, saint_group_method = group_method,
                         saint_init_group_size = init_group_size, saint_max_tokens = max_tokens, saint_kmeans_batch = kmeans_batch, saint_kmeans_iter = kmeans_iter,
                         saint_epochs = epochs, saint_batch_size = batch_size
                       )),
                       "tabnet" = compact(list(
                         tabnet_steps = tabnet_steps, tabnet_feature_dim = tabnet_feature_dim, tabnet_output_dim = tabnet_output_dim,
                         tabnet_gamma = tabnet_gamma, tabnet_lambda_sparse = tabnet_lambda_sparse,
                         tabnet_epochs = epochs, tabnet_batch_size = batch_size
                       )),
                       "node" = compact(list(
                         node_trees = node_trees, node_depth = node_depth, node_epochs = epochs, node_batch_size = batch_size
                       )),
                       "deepfm" = compact(list(
                         deepfm_k = deepfm_k, deepfm_hidden = deepfm_hidden, deepfm_epochs = epochs, deepfm_batch_size = batch_size
                       )),
                       "dcnv2" = compact(list(
                         dcn_layers = dcn_layers, dcn_hidden = dcn_hidden, dcn_epochs = epochs, dcn_batch_size = batch_size
                       )),
                       "nam" = compact(list(
                         nam_hidden = nam_hidden, nam_activation = nam_activation, nam_add_linear = nam_add_linear, nam_l1 = nam_l1,
                         nam_use_grouping = use_grouping, nam_group_trigger = group_trigger, nam_group_method = group_method,
                         nam_init_group_size = init_group_size, nam_max_tokens = max_tokens,
                         nam_kmeans_batch = kmeans_batch, nam_kmeans_iter = kmeans_iter,
                         nam_epochs = epochs, nam_batch_size = batch_size
                       )),
                       "moe" = compact(list(
                         moe_n_experts = moe_n_experts, moe_expert_hidden = moe_expert_hidden, moe_gate_hidden = moe_gate_hidden,
                         moe_temperature = moe_temperature, moe_sparse_topk = moe_sparse_topk, moe_entropy_reg = moe_entropy_reg,
                         moe_use_grouping = use_grouping, moe_group_trigger = group_trigger, moe_group_method = group_method,
                         moe_init_group_size = init_group_size, moe_max_tokens = max_tokens,
                         moe_kmeans_batch = kmeans_batch, moe_kmeans_iter = kmeans_iter,
                         moe_epochs = epochs, moe_batch_size = batch_size
                       )),
                       "gp_dkl" = compact(list(
                         gp_use_variational = gp_use_variational, gp_num_inducing = gp_num_inducing, gp_feature_dim = gp_feature_dim,
                         gp_kernel = gp_kernel, gp_ard = gp_ard, gp_lr_mult = gp_lr_mult,
                         gp_rff_features = rff_features, gp_rff_lengthscale = rff_lengthscale, gp_rff_deep_hidden = rff_deep_hidden,
                         gp_epochs = epochs, gp_batch_size = batch_size
                       )),
                       stop("Unsupported model_type: ", model_type, call. = FALSE)
  )
  dl_args <- c(args_model, args_common)

  # ---------- Cross-validation short path ----------
  if (isTRUE(crossval)) {
    fit_obj <- do.call(dl_fit, c(list(X = X_tr, y = y_tr, model = tolower(model_type)), dl_args))
    pred <- as.numeric(mod$predict(fit_obj[[1L]], X_te,
                                   device = if (!is.null(device)) as.character(device) else NULL))
    pred <- revert_scaling(pred, y_scaler)
    return(pred)
  }

  # ---------- Bootstrap path ----------
  data_label_geno <- cbind(y_train_scaled, X_tr)
  # helper to unscale predictions
  revert_scaling_ml <- function(x, scaler) as.numeric(x * scaler$std + scaler$mean)

  if (!is.null(X_te)) {
    boot_results <- boot::boot(
      data = data_label_geno,
      statistic = train_predict_deeplearning,
      R = n_bootstrap,
      model_type = model_type,
      dl_args = dl_args,
      test_geno = X_te
    )
    # revert all bootstrap predictions
    boot_results$t <- apply(boot_results$t, 2, function(col_vec) revert_scaling_ml(col_vec, y_scaler))
    pred_variances <- apply(boot_results$t, 2, var)
    pred_SE <- apply(boot_results$t, 2, sd)
    AI_pred <- apply(boot_results$t, 2, mean)
    genetic_var <- var(AI_pred)
    AI_pred_reverted <- AI_pred

    result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
      boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
    )
    result_rel <- reliability_thresholds(
      prediction_error_var = pred_variances, genetic_var = genetic_var,
      high_reliability_thres = high_reliability_thres, low_reliability_thres = low_reliability_thres
    )

    train_test_label <- if (!is.null(geno_omic_test_object)) {
      ifelse(rownames(geno_omic_test_object) %in% test_label, "Test", "Train")
    } else rep("Train", nrow(geno_omic_object))

    AI_preds <- data.frame(
      name = GID, Predicted_value = AI_pred_reverted, Train_Test_Label = train_test_label,
      Standard_error = pred_SE, PEV = pred_variances,
      lower_bound = result_rel_MPIW$lower_bound, upper_bound = result_rel_MPIW$upper_bound,
      Uncertainty = result_rel_MPIW$Uncertainty, Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
      Reliability = result_rel$reliability, Reliability_remarks = result_rel$remarks,
      Reliability_percentage = result_rel$reliability_percentage,
      stringsAsFactors = FALSE
    )
    names(AI_preds)[1] <- c(gen_name)

  } else {
    boot_results <- boot::boot(
      data = data_label_geno,
      statistic = train_predict_deeplearning,
      R = n_bootstrap,
      model_type = model_type,
      dl_args = dl_args,
      test_geno = NULL
    )
    pred_variances <- apply(boot_results$t, 2, var)
    pred_SE <- apply(boot_results$t, 2, sd)
    AI_pred <- apply(boot_results$t, 2, mean)
    genetic_var <- var(AI_pred)
    AI_pred_reverted <- revert_scaling_ml(AI_pred, y_scaler)

    result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
      boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
    )
    result_rel <- reliability_thresholds(
      prediction_error_var = pred_variances, genetic_var = genetic_var,
      high_reliability_thres = high_reliability_thres, low_reliability_thres = low_reliability_thres
    )

    train_test_label <- rep("Train", nrow(geno_omic_object))
    AI_preds <- data.frame(
      name = GID, Predicted_value = AI_pred_reverted, Standard_error = pred_SE,
      Train_Test_Label = train_test_label, PEV = pred_variances,
      lower_bound = result_rel_MPIW$lower_bound, upper_bound = result_rel_MPIW$upper_bound,
      Uncertainty = result_rel_MPIW$Uncertainty, Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
      Reliability = result_rel$reliability, Reliability_remarks = result_rel$remarks,
      Reliability_percentage = result_rel$reliability_percentage,
      stringsAsFactors = FALSE
    )
    names(AI_preds)[1] <- c(gen_name)
  }

  # Parameter summary (simple)
  model_para <- data.frame(
    stat = c("model_type","epochs","batch_size"),
    summary = c(model_type, epochs, batch_size),
    stringsAsFactors = FALSE
  )

  list(
    model_parameters = model_para,
    predicted_values = AI_preds,
    diagnostic_plots = diagnostic_plot_true_prediction(
      boot_results = boot_results, GID_names = GID, CI_width_thresholds = CI_width_thresholds,
      predictions = AI_pred_reverted, standard_errors = pred_SE, prediction_error_var = pred_variances,
      genetic_var = genetic_var, confidence_level = 0.95, model_for_CI_cal = "ML",
      high_reliability_thres = high_reliability_thres, low_reliability_thres = low_reliability_thres,
      system_database = system_database
    )
  )
}




#
#
# # =========================================================
# # Helpers
# # =========================================================
# `%||%` <- function(a, b) if (!is.null(a)) a else b
#
# coerce_int_vec <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
#   if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
#   as.integer(x)
# }
# coerce_num <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
#   if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
#   as.numeric(x)
# }
# coerce_logi <- function(x, allow_null = TRUE, name = deparse(substitute(x))) {
#   if (is.null(x)) return(if (allow_null) NULL else stop(name, " cannot be NULL"))
#   isTRUE(x)
# }
#
# assert_matrix <- function(X, name = "X") {
#   if (!is.matrix(X)) {
#     X2 <- try(as.matrix(X), silent = TRUE)
#     if (inherits(X2, "try-error") || !is.matrix(X2)) {
#       stop(name, " must be coercible to a numeric matrix")
#     }
#     X <- X2
#   }
#   storage.mode(X) <- "double"
#   X
# }
#
# infer_task <- function(y) {
#   yv <- if (is.factor(y)) as.integer(y) else y
#   if (is.numeric(yv) && length(unique(na.omit(yv))) > 10) "regression" else "classification"
# }
#
# # 0-based labels for classification
# encode_y <- function(y, task) {
#   if (task == "classification") {
#     if (is.character(y)) y <- factor(y)
#     if (is.factor(y)) {
#       as.integer(y) - 1L
#     } else {
#       lv <- sort(unique(y))
#       match(y, lv) - 1L
#     }
#   } else {
#     as.numeric(y)
#   }
# }
#
# # Synonym mapper for common knobs -> canonical
# normalize_common <- function(args) {
#   # Accept common synonyms
#   if (!is.null(args$l2_regularizer_dp) && is.null(args$l2_weight_decay)) {
#     args$l2_weight_decay <- args$l2_regularizer_dp
#   }
#   if (!is.null(args$dropout_rate) && is.null(args$dropout)) {
#     args$dropout <- args$dropout_rate
#   }
#   if (!is.null(args$batch_normalization) && is.null(args$batch_norm)) {
#     args$batch_norm <- args$batch_normalization
#   }
#   if (!is.null(args$lr) && is.null(args$learning_rate)) {
#     args$learning_rate <- args$lr
#   }
#
#   list(
#     learning_rate    = coerce_num(args$learning_rate %||% 1e-3),
#     epochs           = coerce_int_vec(args$epochs %||% 32L),
#     batch_size       = coerce_int_vec(args$batch_size %||% 64L),
#     l2_weight_decay  = coerce_num(args$l2_weight_decay %||% 1e-3),
#     dropout          = coerce_num(args$dropout %||% 0.5),
#     batch_norm       = coerce_logi(args$batch_norm %||% TRUE),
#     validation_split = coerce_num(args$validation_split %||% 0.2),
#     compile_model    = coerce_logi(args$compile_model %||% TRUE),
#     deterministic    = coerce_logi(args$deterministic %||% FALSE),
#     random_seed      = if (is.null(args$random_seed)) NULL else coerce_int_vec(args$random_seed),
#     device           = if (is.null(args$device)) NULL else as.character(args$device),
#     use_amp          = coerce_logi(args$use_amp %||% TRUE),
#     max_grad_norm    = if (is.null(args$max_grad_norm)) NULL else coerce_num(args$max_grad_norm),
#     optimizer_name   = as.character(args$optimizer_name %||% "adam"),
#     heteroscedastic  = coerce_logi(args$heteroscedastic %||% FALSE),
#     auto_class_weights = coerce_logi(args$auto_class_weights %||% FALSE)
#   )
# }
#
# # ---------- Alias helper ----------
# apply_aliases <- function(dots, aliases) {
#   for (src in names(aliases)) {
#     dst <- aliases[[src]]
#     if (!is.null(dots[[src]]) && is.null(dots[[dst]])) {
#       dots[[dst]] <- dots[[src]]
#     }
#     dots[[src]] <- NULL
#   }
#   dots
# }
#
# # =========================================================
# # Alias maps (prefixed to canonical) for ALL MODELS
# # =========================================================
#
# # Common canonical keys used across models
# .COMMON_CANON <- c(
#   "learning_rate","epochs","batch_size","l2_weight_decay","dropout",
#   "batch_norm","validation_split","compile_model","deterministic",
#   "random_seed","device","use_amp","max_grad_norm","optimizer_name",
#   "heteroscedastic","auto_class_weights","task"
# )
#
# # ---------- common prefixed to canonical this is very essential  ----------
# aliases_common_prefixed <- function(prefix) {
#   stats::setNames(
#     c("learning_rate","epochs","batch_size","l2_weight_decay","l2_weight_decay",
#       "dropout","dropout","batch_norm","validation_split","compile_model",
#       "deterministic","random_seed","device","use_amp","max_grad_norm",
#       "optimizer_name","heteroscedastic","auto_class_weights"),
#     c(paste0(prefix,"_learning_rate"),
#       paste0(prefix,"_epochs"),
#       paste0(prefix,"_batch_size"),
#       paste0(prefix,"_l2_weight_decay"),
#       paste0(prefix,"_l2_regularizer_dp"),
#       paste0(prefix,"_dropout"),
#       paste0(prefix,"_dropout_rate"),
#       paste0(prefix,"_batch_norm"),
#       paste0(prefix,"_validation_split"),
#       paste0(prefix,"_compile_model"),
#       paste0(prefix,"_deterministic"),
#       paste0(prefix,"_random_seed"),
#       paste0(prefix,"_device"),
#       paste0(prefix,"_use_amp"),
#       paste0(prefix,"_max_grad_norm"),
#       paste0(prefix,"_optimizer_name"),
#       paste0(prefix,"_heteroscedastic"),
#       paste0(prefix,"_auto_class_weights"))
#   )
# }
#
# # ---------- per-model alias maps ----------
# aliases_cnn <- function() c(
#   "cnn_neurons_per_layer" = "neurons_per_layer",
#   "cnn_kernel_size"       = "kernel_size",
#   "cnn_dense_layers"      = "dense_layers_cnn",
#   "cnn_use_max_pool"      = "cnn_use_max_pool",
#   "use_max_pool"          = "cnn_use_max_pool",
#   "cnn_pool_kernel"       = "cnn_pool_kernel",
#   "cnn_pool_stride"       = "cnn_pool_stride",
#   "cnn_pool_padding"      = "cnn_pool_padding",
#   "cnn_separable"         = "separable",
#   "cnn_dilations"         = "dilations",
#   "cnn_use_se"            = "use_se",
#   "cnn_norm_type"         = "norm_type",
#   "cnn_pool_type"         = "pool_type",
#   "cnn_use_global_pool"   = "use_global_pool",
#   aliases_common_prefixed("cnn")
# )
#
#
# aliases_mlp <- function() c(
#   "mlp_neurons_per_layer"       = "neurons_per_layer",
#   "mlp_final_attention"         = "final_attention",
#   "mlp_attention_across_layers" = "attention_across_multiple_layers",
#   aliases_common_prefixed("mlp")
# )
#
# aliases_ft <- function() c(
#   "ft_d_model"          = "ft_d_model",
#   "ft_heads"            = "ft_heads",
#   "ft_layers"           = "ft_layers",
#   "ft_ff_mult"          = "ft_ff_mult",
#   "ft_dropout"          = "ft_dropout",
#   "ft_token_dropout"    = "ft_token_dropout",
#   "ft_use_cls"          = "ft_use_cls",
#   "ft_scalar_tokenizer" = "ft_scalar_tokenizer",
#   # grouping (prefixed to canonical)
#   "ft_use_grouping"     = "use_grouping",
#   "ft_group_trigger"    = "group_trigger",
#   "ft_group_method"     = "group_method",
#   "ft_init_group_size"  = "init_group_size",
#   "ft_max_tokens"       = "max_tokens",
#   "ft_kmeans_batch"     = "kmeans_batch",
#   "ft_kmeans_iter"      = "kmeans_iter",
#   aliases_common_prefixed("ft")
# )
#
# aliases_resnet <- function() c(
#   "resnet_neurons_per_block" = "neurons_per_layer",
#   "resnet_blocks"            = "num_hidden_layers",
#   aliases_common_prefixed("resnet")
# )
#
# aliases_saint <- function() c(
#   "saint_d_model"       = "saint_d_model",
#   "saint_heads"         = "saint_heads",
#   "saint_layers"        = "saint_layers",
#   "saint_ff_mult"       = "saint_ff_mult",
#   "saint_dropout"       = "saint_dropout",
#   "saint_token_dropout" = "saint_token_dropout",
#   "saint_use_cls"       = "saint_use_cls",
#   "saint_use_grouping"  = "use_grouping",
#   "saint_group_trigger" = "group_trigger",
#   "saint_group_method"  = "group_method",
#   "saint_init_group_size" = "init_group_size",
#   "saint_max_tokens"    = "max_tokens",
#   "saint_kmeans_batch"  = "kmeans_batch",
#   "saint_kmeans_iter"   = "kmeans_iter",
#   aliases_common_prefixed("saint")
# )
#
# aliases_tabnet <- function() c(
#   # tabnet and others are already canonical, only common prefixed knobs were added:
#   aliases_common_prefixed("tabnet")
# )
#
# aliases_node <- function() c(
#   aliases_common_prefixed("node")
# )
#
# aliases_deepfm <- function() c(
#   aliases_common_prefixed("deepfm")
# )
#
# aliases_dcnv2 <- function() c(
#   # "dcn_" is chosen prefix for DCNv2 commons
#   aliases_common_prefixed("dcn")
# )
#
# aliases_nam <- function() c(
#   "nam_hidden"        = "nam_hidden",
#   "nam_activation"    = "nam_activation",
#   "nam_add_linear"    = "nam_add_linear",
#   "nam_l1"            = "nam_l1",
#   "nam_use_grouping"  = "use_grouping",
#   "nam_group_trigger" = "group_trigger",
#   "nam_group_method"  = "group_method",
#   "nam_init_group_size" = "init_group_size",
#   "nam_max_tokens"    = "max_tokens",
#   "nam_kmeans_batch"  = "kmeans_batch",
#   "nam_kmeans_iter"   = "kmeans_iter",
#   aliases_common_prefixed("nam")
# )
#
# aliases_moe <- function() c(
#   "moe_n_experts"     = "moe_n_experts",
#   "moe_expert_hidden" = "moe_expert_hidden",
#   "moe_gate_hidden"   = "moe_gate_hidden",
#   "moe_temperature"   = "moe_temperature",
#   "moe_sparse_topk"   = "moe_sparse_topk",
#   "moe_entropy_reg"   = "moe_entropy_reg",
#   "moe_use_grouping"  = "use_grouping",
#   "moe_group_trigger" = "group_trigger",
#   "moe_group_method"  = "group_method",
#   "moe_init_group_size" = "init_group_size",
#   "moe_max_tokens"    = "max_tokens",
#   "moe_kmeans_batch"  = "kmeans_batch",
#   "moe_kmeans_iter"   = "kmeans_iter",
#   aliases_common_prefixed("moe")
# )
#
# aliases_gp_dkl <- function() c(
#   "gp_use_variational" = "gp_use_variational",
#   "gp_num_inducing"    = "gp_num_inducing",
#   "gp_feature_dim"     = "gp_feature_dim",
#   "gp_kernel"          = "gp_kernel",
#   "gp_ard"             = "gp_ard",
#   "gp_lr_mult"         = "gp_lr_mult",
#   "gp_rff_features"    = "rff_features",
#   "gp_rff_lengthscale" = "rff_lengthscale",
#   "gp_rff_deep_hidden" = "rff_deep_hidden",
#   aliases_common_prefixed("gp")
# )
#
# # ---- GP-DKL kernel choices & aliases ----
# GP_KERNEL_CHOICES <- c(
#   "rbf", "matern12", "matern32", "matern52",
#   "rq", "linear", "polynomial",
#   "cosine", "periodic", "spectral_mixture"
# )
#
# # optional aliases (case-insensitive)
# GP_KERNEL_ALIASES <- c(
#   "se" = "rbf", "sqexp" = "rbf", "exp_quad" = "rbf",
#   "cos" = "cosine",
#   "matern" = "matern32"  # default if user just says "matern"
# )
#
# canonicalize_gp_kernel <- function(k) {
#   if (is.null(k) || length(k) == 0) return("rbf")
#   k <- tolower(trimws(as.character(k)))
#
#   # apply alias only if present
#   if (exists("GP_KERNEL_ALIASES", inherits = TRUE) &&
#       k %in% names(GP_KERNEL_ALIASES)) {
#     k <- unname(GP_KERNEL_ALIASES[k])
#   }
#
#   # validate against allowed set
#   if (!(k %in% GP_KERNEL_CHOICES)) {
#     stop(sprintf(
#       "Unknown gp_kernel='%s'. Allowed: %s",
#       k, paste(GP_KERNEL_CHOICES, collapse = ", ")
#     ), call. = FALSE)
#   }
#   k
# }
#
#
# # =========================================================
# # Registry
# # =========================================================
# dl_registry <- new.env(parent = emptyenv())
# register_model <- function(key, fn, allowed_args) {
#   dl_registry[[key]] <- list(fn = fn, allowed = allowed_args)
# }
# get_model <- function(key) dl_registry[[key]]
#
# torch_fit_filtered <- function(args) {
#   sig <- try(names(formals(torch_fit_model)), silent = TRUE)
#   if (!inherits(sig, "try-error") && length(sig)) {
#     args <- args[intersect(names(args), sig)]
#   }
#   do.call(torch_fit_model, args)
# }
#
# # =========================================================
# # Model-specific wrappers
# # =========================================================
#
# # ---------- CNN ----------
# dl_cnn <- function(X, y, ..., model_type = "cnn") {
#   dots <- apply_aliases(list(...), aliases_cnn())
#
#   conv_filters    <- coerce_int_vec(dots$neurons_per_layer %||% c(64L,64L,64L))
#   kernel_size     <- coerce_int_vec(dots$kernel_size %||% 3L)
#   dense_layers    <- coerce_int_vec(dots$dense_layers_cnn %||% c(256L,128L,64L))
#   cnn_use_max     <- coerce_logi(dots$cnn_use_max_pool %||% FALSE)   # canonical only
#   pool_k          <- coerce_int_vec(dots$cnn_pool_kernel %||% 2L)
#   pool_s          <- coerce_int_vec(dots$cnn_pool_stride %||% 2L)
#   pool_p          <- coerce_int_vec(dots$cnn_pool_padding %||% 0L)
#   separable       <- coerce_logi(dots$separable %||% TRUE)
#   dilations       <- coerce_int_vec(dots$dilations %||% c(1L,2L,4L))
#   use_se          <- coerce_logi(dots$use_se %||% TRUE)
#   norm_type       <- as.character(dots$norm_type %||% "group")
#   pool_type       <- as.character(dots$pool_type %||% "conv")
#   use_global_pool <- coerce_logi(dots$use_global_pool %||% FALSE)
#
#   C <- normalize_common(dots)
#
#   # backend call WITHOUT `use_max_pool = ...`
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "cnn",
#     neurons_per_layer = conv_filters,
#     kernel_size      = kernel_size,
#     dense_layers_cnn = dense_layers,
#     cnn_use_max_pool = cnn_use_max,
#     cnn_pool_kernel  = pool_k,
#     cnn_pool_stride  = pool_s,
#     cnn_pool_padding = pool_p,
#     separable        = separable,
#     dilations        = dilations,
#     use_se           = use_se,
#     norm_type        = norm_type,
#     pool_type        = pool_type,
#     use_global_pool  = use_global_pool,
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
#
# # function to Call torch_fit_model while auto-dropping unknown args
# torch_fit_filtered <- function(args) {
#   sig <- try(names(formals(torch_fit_model)), silent = TRUE)
#   if (!inherits(sig, "try-error") && length(sig)) {
#     args <- args[intersect(names(args), sig)]
#   }
#   do.call(torch_fit_model, args)
# }
#
# # ---------- MLP (with/without attention) ----------
# dl_mlp <- function(X, y, ..., model_type = c("mlp","mlp_with_attention")) {
#   model_type <- match.arg(model_type)
#   dots <- apply_aliases(list(...), aliases_mlp())
#
#   hidden    <- coerce_int_vec(dots$neurons_per_layer %||% c(256L,128L))
#   final_att <- coerce_logi(dots$final_attention %||% (model_type == "mlp_with_attention"))
#   across    <- coerce_logi(dots$attention_across_multiple_layers %||% FALSE)
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = model_type,
#     neurons_per_layer = hidden,
#     final_attention = final_att,
#     attention_across_multiple_layers = across,
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   )
#   )
# }
#
# # ---------- FT-Transformer ----------
# dl_ft_transformer <- function(X, y, ..., model_type = "ft_transformer") {
#   dots <- apply_aliases(list(...), aliases_ft())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "ft_transformer",
#     ft_d_model       = coerce_int_vec(dots$ft_d_model %||% 192L),
#     ft_heads         = coerce_int_vec(dots$ft_heads %||% 8L),
#     ft_layers        = coerce_int_vec(dots$ft_layers %||% 4L),
#     ft_ff_mult       = coerce_int_vec(dots$ft_ff_mult %||% 4L),
#     ft_dropout       = coerce_num(dots$ft_dropout %||% 0.1),
#     ft_token_dropout = coerce_num(dots$ft_token_dropout %||% 0.1),
#     ft_use_cls       = coerce_logi(dots$ft_use_cls %||% TRUE),
#     ft_scalar_tokenizer = coerce_logi(dots$ft_scalar_tokenizer %||% TRUE),
#
#     use_grouping     = coerce_logi(dots$use_grouping %||% TRUE),
#     group_trigger    = coerce_int_vec(dots$group_trigger %||% 2048L),
#     group_method     = as.character(dots$group_method %||% "auto"),
#     init_group_size  = coerce_int_vec(dots$init_group_size %||% 64L),
#     max_tokens       = coerce_int_vec(dots$max_tokens %||% 1024L),
#     kmeans_batch     = coerce_int_vec(dots$kmeans_batch %||% 4096L),
#     kmeans_iter      = coerce_int_vec(dots$kmeans_iter %||% 100L),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   )
#   )
# }
#
# # ---------- ResNet ----------
# dl_resnet <- function(X, y, ..., model_type = "resnet") {
#   dots <- apply_aliases(list(...), aliases_resnet())
#   C <- normalize_common(dots)
#   blocks <- if (!is.null(dots$num_hidden_layers)) coerce_int_vec(dots$num_hidden_layers) else NULL
#   hidden <- coerce_int_vec(dots$neurons_per_layer %||% c(64L,64L,64L))
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "resnet",
#     num_hidden_layers = blocks,
#     neurons_per_layer = hidden,
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- SAINT ----------
# dl_saint <- function(X, y, ..., model_type = "saint") {
#   dots <- apply_aliases(list(...), aliases_saint())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "saint",
#     saint_d_model       = coerce_int_vec(dots$saint_d_model %||% 128L),
#     saint_heads         = coerce_int_vec(dots$saint_heads %||% 8L),
#     saint_layers        = coerce_int_vec(dots$saint_layers %||% 4L),
#     saint_ff_mult       = coerce_int_vec(dots$saint_ff_mult %||% 4L),
#     saint_dropout       = coerce_num(dots$saint_dropout %||% 0.1),
#     saint_token_dropout = coerce_num(dots$saint_token_dropout %||% 0.1),
#     saint_use_cls       = coerce_logi(dots$saint_use_cls %||% TRUE),
#
#     use_grouping     = coerce_logi(dots$use_grouping %||% TRUE),
#     group_trigger    = coerce_int_vec(dots$group_trigger %||% 2048L),
#     group_method     = as.character(dots$group_method %||% "auto"),
#     init_group_size  = coerce_int_vec(dots$init_group_size %||% 64L),
#     max_tokens       = coerce_int_vec(dots$max_tokens %||% 1024L),
#     kmeans_batch     = coerce_int_vec(dots$kmeans_batch %||% 4096L),
#     kmeans_iter      = coerce_int_vec(dots$kmeans_iter %||% 100L),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- TabNet ----------
# dl_tabnet <- function(X, y, ..., model_type = "tabnet") {
#   dots <- apply_aliases(list(...), aliases_tabnet())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "tabnet",
#     tabnet_steps        = coerce_int_vec(dots$tabnet_steps %||% 5L),
#     tabnet_feature_dim  = coerce_int_vec(dots$tabnet_feature_dim %||% 64L),
#     tabnet_output_dim   = coerce_int_vec(dots$tabnet_output_dim %||% 64L),
#     tabnet_gamma        = coerce_num(dots$tabnet_gamma %||% 1.5),
#     tabnet_lambda_sparse= coerce_num(dots$tabnet_lambda_sparse %||% 1e-4),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- NODE ----------
# dl_node <- function(X, y, ..., model_type = "node") {
#   dots <- apply_aliases(list(...), aliases_node())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "node",
#     node_trees = coerce_int_vec(dots$node_trees %||% 8L),
#     node_depth = coerce_int_vec(dots$node_depth %||% 3L),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- DeepFM ----------
# dl_deepfm <- function(X, y, ..., model_type = "deepfm") {
#   dots <- apply_aliases(list(...), aliases_deepfm())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "deepfm",
#     deepfm_k      = coerce_int_vec(dots$deepfm_k %||% 16L),
#     deepfm_hidden = coerce_int_vec(dots$deepfm_hidden %||% c(128L,64L)),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- DCNv2 ----------
# dl_dcnv2 <- function(X, y, ..., model_type = "dcnv2") {
#   dots <- apply_aliases(list(...), aliases_dcnv2())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "dcnv2",
#     dcn_layers = coerce_int_vec(dots$dcn_layers %||% 3L),
#     dcn_hidden = coerce_int_vec(dots$dcn_hidden %||% c(256L,128L)),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- NAM ----------
# dl_nam <- function(X, y, ..., model_type = "nam") {
#   dots <- apply_aliases(list(...), aliases_nam())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "nam",
#     nam_hidden     = coerce_int_vec(dots$nam_hidden %||% c(32L,16L)),
#     nam_activation = as.character(dots$nam_activation %||% "relu"),
#     nam_add_linear = coerce_logi(dots$nam_add_linear %||% TRUE),
#     nam_l1         = coerce_num(dots$nam_l1 %||% 1e-4),
#
#     use_grouping     = coerce_logi(dots$use_grouping %||% TRUE),
#     group_trigger    = coerce_int_vec(dots$group_trigger %||% 2048L),
#     group_method     = as.character(dots$group_method %||% "auto"),
#     init_group_size  = coerce_int_vec(dots$init_group_size %||% 64L),
#     max_tokens       = coerce_int_vec(dots$max_tokens %||% 1024L),
#     kmeans_batch     = coerce_int_vec(dots$kmeans_batch %||% 4096L),
#     kmeans_iter      = coerce_int_vec(dots$kmeans_iter %||% 100L),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- MoE ----------
# dl_moe <- function(X, y, ..., model_type = "moe") {
#   dots <- apply_aliases(list(...), aliases_moe())
#   C <- normalize_common(dots)
#
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "moe",
#     moe_n_experts     = coerce_int_vec(dots$moe_n_experts %||% 4L),
#     moe_expert_hidden = coerce_int_vec(dots$moe_expert_hidden %||% c(128L,64L)),
#     moe_gate_hidden   = coerce_int_vec(dots$moe_gate_hidden %||% 128L),
#     moe_temperature   = coerce_num(dots$moe_temperature %||% 1.0),
#     moe_sparse_topk   = if (is.null(dots$moe_sparse_topk)) NULL else coerce_int_vec(dots$moe_sparse_topk),
#     moe_entropy_reg   = coerce_num(dots$moe_entropy_reg %||% 0.0),
#
#     use_grouping     = coerce_logi(dots$use_grouping %||% TRUE),
#     group_trigger    = coerce_int_vec(dots$group_trigger %||% 2048L),
#     group_method     = as.character(dots$group_method %||% "auto"),
#     init_group_size  = coerce_int_vec(dots$init_group_size %||% 64L),
#     max_tokens       = coerce_int_vec(dots$max_tokens %||% 1024L),
#     kmeans_batch     = coerce_int_vec(dots$kmeans_batch %||% 4096L),
#     kmeans_iter      = coerce_int_vec(dots$kmeans_iter %||% 100L),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # ---------- GP-DKL (+ RFF fallback) ----------
# dl_gp_dkl <- function(X, y, ..., model_type = "gp_dkl") {
#   dots <- apply_aliases(list(...), aliases_gp_dkl())
#   C <- normalize_common(dots)
#   kern <- canonicalize_gp_kernel(dots$gp_kernel %||% "rbf")
#   torch_fit_filtered(list(
#     X = X, y = y, model_type = "gp_dkl",
#     gp_use_variational = coerce_logi(dots$gp_use_variational %||% TRUE),
#     gp_num_inducing    = coerce_int_vec(dots$gp_num_inducing %||% 512L),
#     gp_feature_dim     = coerce_int_vec(dots$gp_feature_dim %||% 64L),
#     #gp_kernel          = as.character(dots$gp_kernel %||% "rbf"),
#     gp_kernel          = kern,
#     gp_ard             = coerce_logi(dots$gp_ard %||% TRUE),
#     gp_lr_mult         = coerce_num(dots$gp_lr_mult %||% 0.5),
#
#     rff_features    = coerce_int_vec(dots$rff_features %||% 1024L),
#     rff_lengthscale = coerce_num(dots$rff_lengthscale %||% 1.0),
#     rff_deep_hidden = coerce_int_vec(dots$rff_deep_hidden %||% c(128L)),
#
#     # common
#     learning_rate    = C$learning_rate,
#     epochs           = C$epochs,
#     batch_size       = C$batch_size,
#     l2_weight_decay  = C$l2_weight_decay,
#     dropout          = C$dropout,
#     batch_norm       = C$batch_norm,
#     validation_split = C$validation_split,
#     compile_model    = C$compile_model,
#     deterministic    = C$deterministic,
#     random_seed      = C$random_seed,
#     device           = C$device,
#     use_amp          = C$use_amp,
#     max_grad_norm    = C$max_grad_norm,
#     optimizer_name   = C$optimizer_name,
#     heteroscedastic  = C$heteroscedastic,
#     auto_class_weights = C$auto_class_weights
#   ))
# }
#
# # =========================================================
# # Register models + canonical allowed args per model
# # =========================================================
# register_model("cnn", dl_cnn,
#                 allowed_args = c("neurons_per_layer","kernel_size","dense_layers_cnn",
#                                  "cnn_use_max_pool","cnn_pool_kernel","cnn_pool_stride","cnn_pool_padding",
#                                  "separable","dilations","use_se","norm_type","pool_type","use_global_pool")
# )
#
# register_model("mlp", function(...) dl_mlp(..., model_type = "mlp"),
#                 allowed_args = c("neurons_per_layer","final_attention","attention_across_multiple_layers")
# )
# register_model("mlp_with_attention", function(...) dl_mlp(..., model_type = "mlp_with_attention"),
#                 allowed_args = c("neurons_per_layer","final_attention","attention_across_multiple_layers")
# )
#
# register_model("ft_transformer", dl_ft_transformer,
#                 allowed_args = c("ft_d_model","ft_heads","ft_layers","ft_ff_mult","ft_dropout","ft_token_dropout",
#                                  "ft_use_cls","ft_scalar_tokenizer","use_grouping","group_trigger","group_method",
#                                  "init_group_size","max_tokens","kmeans_batch","kmeans_iter")
# )
#
# register_model("resnet", dl_resnet,
#                 allowed_args = c("neurons_per_layer","num_hidden_layers")
# )
#
# register_model("saint", dl_saint,
#                 allowed_args = c("saint_d_model","saint_heads","saint_layers","saint_ff_mult","saint_dropout","saint_token_dropout",
#                                  "saint_use_cls","use_grouping","group_trigger","group_method","init_group_size","max_tokens",
#                                  "kmeans_batch","kmeans_iter")
# )
#
# register_model("tabnet", dl_tabnet,
#                 allowed_args = c("tabnet_steps","tabnet_feature_dim","tabnet_output_dim","tabnet_gamma","tabnet_lambda_sparse")
# )
#
# register_model("node", dl_node,
#                 allowed_args = c("node_trees","node_depth")
# )
#
# register_model("deepfm", dl_deepfm,
#                 allowed_args = c("deepfm_k","deepfm_hidden")
# )
#
# register_model("dcnv2", dl_dcnv2,
#                 allowed_args = c("dcn_layers","dcn_hidden")
# )
#
# register_model("nam", dl_nam,
#                 allowed_args = c("nam_hidden","nam_activation","nam_add_linear","nam_l1",
#                                  "use_grouping","group_trigger","group_method","init_group_size","max_tokens",
#                                  "kmeans_batch","kmeans_iter")
# )
#
# register_model("moe", dl_moe,
#                 allowed_args = c("moe_n_experts","moe_expert_hidden","moe_gate_hidden","moe_temperature",
#                                  "moe_sparse_topk","moe_entropy_reg","use_grouping","group_trigger","group_method",
#                                  "init_group_size","max_tokens","kmeans_batch","kmeans_iter")
# )
#
# register_model("gp_dkl", dl_gp_dkl,
#                 allowed_args = c("gp_use_variational","gp_num_inducing","gp_feature_dim","gp_kernel","gp_ard","gp_lr_mult",
#                                  "rff_features","rff_lengthscale","rff_deep_hidden")
# )
#
# # =========================================================
# # Coordinator with centralized allowed-args guard (alias-aware)
# # =========================================================
# dl_fit <- function(
#     X, y,
#     model = c("cnn","mlp","mlp_with_attention","ft_transformer",
#               "resnet","saint","tabnet","node","deepfm","dcnv2",
#               "nam","moe","gp_dkl"),
#     ...
# ) {
#   model <- match.arg(model)
#   dots  <- list(...)
#
#   X <- assert_matrix(X)
#   task <- infer_task(y)
#   y_enc <- encode_y(y, task)
#   dots$task <- task
#
#   reg <- get_model(model)
#   if (is.null(reg)) stop("Model '", model, "' is not registered yet.")
#
#   # Build alias sets for the selected model
#   allowed_prefixed <- switch(
#     model,
#     cnn   = c(names(aliases_cnn())),
#     mlp   = c(names(aliases_mlp())),
#     mlp_with_attention = c(names(aliases_mlp())),
#     ft_transformer = c(names(aliases_ft())),
#     resnet= c(names(aliases_resnet())),
#     saint = c(names(aliases_saint())),
#     tabnet= c(names(aliases_tabnet())),
#     node  = c(names(aliases_node())),
#     deepfm= c(names(aliases_deepfm())),
#     dcnv2 = c(names(aliases_dcnv2())),
#     nam   = c(names(aliases_nam())),
#     moe   = c(names(aliases_moe())),
#     gp_dkl= c(names(aliases_gp_dkl())),
#     character()
#   )
#
#   allowed_canon <- unique(c(reg$allowed %||% character(), .COMMON_CANON))
#   allowed_all   <- unique(c(allowed_canon, allowed_prefixed))
#
#   extra <- setdiff(names(dots), allowed_all)
#   if (length(extra)) {
#     warning(sprintf(
#       "Ignoring %d unsupported argument(s) for model '%s': %s",
#       length(extra), model, paste(extra, collapse = ", ")
#     ), call. = FALSE)
#   }
#
#   # Dispatch
#   do.call(reg$fn, c(list(X = X, y = y_enc), dots))
# }
