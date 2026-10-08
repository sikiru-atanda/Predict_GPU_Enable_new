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
    if (is.character(y) || is.logical(y)) y <- factor(y)
    if (is.factor(y)) return(as.integer(y) - 1L)
    lv <- sort(unique(y))
    if (length(lv) && identical(as.numeric(lv), seq(0, length(lv) - 1))) {
      return(as.integer(y))
    }
    if (length(lv) && identical(as.numeric(lv), seq_len(length(lv)))) {
      return(as.integer(y) - 1L)
    }
    match(y, lv) - 1L
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
    # identity entries (e.g. "ft_layers"="ft_layers") are already canonical;
    # removing src there would drop the user's value
    if (identical(src, dst)) next
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
  if (!inherits(sig, "try-error") && length(sig) && !"..." %in% sig) {
    args <- args[intersect(names(args), sig)]
  }
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

  # Legacy -> current renames
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
  # CNN pool aliases (old -> new)
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

gp_dl_registered_model_types <- function() {
  sort(ls(envir = dl_registry, all.names = TRUE))
}

gp_dl_supports_local_target_uncertainty <- function(model_type) {
  model_type <- tolower(trimws(as.character(model_type %||% "")[1L]))
  nzchar(model_type) && exists(
    model_type,
    envir = dl_registry,
    inherits = FALSE
  )
}

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
  fam <- gp_resolve_response_family(dl_args$response_family %||% "gaussian", y = y_tr)
  if (identical(fam, "gaussian")) {
    preds <- mod$predict(fit[[1L]], X_targ,
                         device = if (!is.null(dl_args$device)) as.character(dl_args$device) else NULL)
    return(as.numeric(preds))
  }
  preds <- torch_predict(
    fit[[1L]], X_targ,
    device = if (!is.null(dl_args$device)) as.character(dl_args$device) else NULL,
    type = "prob"
  )
  if (identical(fam, "binary")) {
    return(as.numeric(preds))
  }
  as.numeric(preds)
}

gp_dl_y_cache_signature <- function(y) {
  if (is.factor(y) || is.character(y) || is.logical(y)) {
    y_chr <- as.character(y)
    non_na <- y_chr[!is.na(y_chr)]
    return(paste(
      length(y_chr),
      sum(is.na(y_chr)),
      paste(sort(unique(non_na)), collapse = "\r"),
      sum(nchar(non_na, type = "bytes")),
      sep = "::"
    ))
  }
  y_num <- as.numeric(y)
  paste(
    length(y_num),
    sum(is.na(y_num)),
    format(sum(y_num, na.rm = TRUE), digits = 12),
    format(sum(y_num * y_num, na.rm = TRUE), digits = 12),
    sep = "::"
  )
}

gp_dl_cv_cache_key <- function(omics_data, y_raw, tst, scaling = TRUE, centering = FALSE,
                               response_family = "auto", target_data = NULL) {
  x <- as.matrix(omics_data)
  target <- if (is.null(target_data)) {
    matrix(numeric(), nrow = 0L, ncol = ncol(x))
  } else {
    as.matrix(target_data)
  }
  paste(
    nrow(x),
    ncol(x),
    format(sum(x, na.rm = TRUE), digits = 12),
    format(sum(x * x, na.rm = TRUE), digits = 12),
    gp_dl_y_cache_signature(y_raw),
    length(tst),
    sum(as.integer(tst)),
    isTRUE(scaling),
    isTRUE(centering),
    as.character(response_family[[1]] %||% "auto"),
    nrow(target),
    ncol(target),
    format(sum(target, na.rm = TRUE), digits = 12),
    format(sum(target * target, na.rm = TRUE), digits = 12),
    sep = "::"
  )
}

gp_dl_get_cached_cv_split <- function(cache_key, build_fun) {
  key <- paste0("cvsplit::", cache_key)
  if (exists(key, envir = dl_env, inherits = FALSE)) {
    return(get(key, envir = dl_env, inherits = FALSE))
  }
  val <- build_fun()
  assign(key, val, envir = dl_env)
  val
}

gp_dl_get_cached_cv_payload <- function(omics_data,
                                        y_raw,
                                        tst,
                                        fam,
                                        scaling = TRUE,
                                        centering = FALSE,
                                        target_data = NULL) {
  cache_key <- gp_dl_cv_cache_key(
    omics_data = omics_data,
    y_raw = y_raw,
    tst = tst,
    scaling = scaling,
    centering = centering,
    response_family = fam,
    target_data = target_data
  )
  gp_dl_get_cached_cv_split(cache_key, function() {
    x_raw <- as.matrix(omics_data)
    x_prep <- gp_ml_preprocess_predictors(
      x_raw[-tst, , drop = FALSE],
      scaling = scaling,
      centering = centering
    )
    x_all <- as.matrix(gp_ml_apply_preprocessor(x_raw, x_prep))
    storage.mode(x_all) <- "double"
    x_target <- if (is.null(target_data)) {
      x_all
    } else {
      as.matrix(gp_ml_apply_preprocessor(target_data, x_prep))
    }
    storage.mode(x_target) <- "double"

    if (identical(fam, "gaussian")) {
      y_train <- as.numeric(y_raw)
      y_scaler <- gp_fast_y_scale_fit(y_train)
      y_processed <- y_scaler$scaled
      class_levels <- NULL
    } else {
      class_levels <- if (exists("gp_py_ml_class_levels", mode = "function")) {
        gp_py_ml_class_levels(y_raw, fam)
      } else {
        gp_response_class_levels(y_raw, fam)
      }
      if (identical(fam, "binary") && length(class_levels) != 2L) {
        stop("Binary deep-learning support requires exactly two observed classes.", call. = FALSE)
      }
      y_processed <- as.integer(factor(as.character(y_raw), levels = class_levels)) - 1L
      y_scaler <- NULL
    }

    list(
      omics_data = x_all,
      X_tr = unname(as.matrix(x_all[-tst, , drop = FALSE])),
      X_te = unname(as.matrix(x_all[tst, , drop = FALSE])),
      X_target = unname(x_target),
      y_tr = y_processed[-tst],
      y_te = y_processed[tst],
      y_scaler = y_scaler,
      y_train_processed = y_processed,
      class_levels = class_levels
    )
  })
}

# =========================================================
# Main coordinator with boot + CV options
# =========================================================
deep_learning_model <- function(pheno_object = NULL,
                                y = NULL,
                                omics_data = NULL,
                                omics = NULL,
                                crossval = FALSE,
                                response_family = "auto",
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
                                 n_bootstrap = 30,
                                 dl_internal_calibration = TRUE,
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
                                deterministic = TRUE, random_seed = 123,
                                dl_n_seeds = NULL, dl_seeds = NULL,
                                dl_seed_aggregation = "mean", device = NULL,
                                # MLP
                                mlp_neurons_per_layer = as.integer(c(128,64)),
                                mlp_learning_rate = 1e-3,
                                final_attention = TRUE,
                                attention_across_multiple_layers = TRUE,
                                heteroscedastic = FALSE,
                                ...) {
  rng_seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng_seed_before <- if (rng_seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    gp_clear_invalid_random_seed(rng_seed_before, rng_seed_existed),
    add = TRUE
  )

  msg <- ""
  if (is.null(omics_data) && !is.null(omics)) {
    omics_data <- omics
  }
  # ---------- Prepare X/y + scaling ----------
  if (is.null(geno_omic_object) && is.null(pheno_object) && isFALSE(crossval)) {
    stop(paste(msg, "provide matrix of the predictors and the data.frame of the Y variable."), call. = FALSE)
  }

  y_train_raw <- if (isFALSE(crossval)) pheno_object[, response] else y
  fam <- gp_resolve_response_family(response_family, y = y_train_raw)
  if (identical(fam, "ordinal")) {
    stop("Deep-learning public classification support is currently implemented for gaussian, binary, and multiclass traits; ordinal is planned next.", call. = FALSE)
  }
  dl_seed_aggregation <- gp_dl_seed_aggregation(dl_seed_aggregation)
  dl_seed_manifest <- gp_dl_seed_manifest(
    dl_n_seeds = dl_n_seeds,
    dl_seeds = dl_seeds,
    random_seed = random_seed
  )
  dl_primary_training_seed <- dl_seed_manifest$training_seed[[1L]]
  dl_fit_plan <- NULL
  cv_payload <- NULL

  if (!is.null(geno_omic_object)) {
    GID <- rownames(geno_omic_object)
    predictor_prep <- gp_ml_preprocess_predictors(geno_omic_object, scaling = scaling, centering = centering)
    geno_omic_object <- predictor_prep$data
  } else {
    predictor_prep <- NULL
  }

  if (!is.null(geno_omic_test_object)) {
    test_label <- rownames(geno_omic_test_object)
    if (!is.null(predictor_prep)) {
      geno_omic_test_object <- gp_ml_apply_preprocessor(geno_omic_test_object, predictor_prep)
    }
    geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
    GID <- rownames(geno_omic_test_object)
  }

  if (!is.null(omics_data) && isTRUE(crossval)) {
    if (is.null(tst)) stop("crossval=TRUE but 'tst' indices are NULL", call. = FALSE)
    omics_data <- as.matrix(omics_data)
    stopifnot(is.numeric(tst), all(tst >= 1), all(tst <= nrow(omics_data)))
    cv_payload <- gp_dl_get_cached_cv_payload(
      omics_data = omics_data,
      y_raw = y_train_raw,
      tst = tst,
      fam = fam,
      scaling = scaling,
      centering = centering
    )
    omics_data <- cv_payload$omics_data
  }


  if (isTRUE(crossval)) { para_tunning <- FALSE; param_grid <- NULL }
  tuning_summary <- NULL

  # y and scaling
  if (!is.null(cv_payload)) {
    y_scaler <- cv_payload$y_scaler
    y_train_processed <- cv_payload$y_train_processed
    class_levels <- cv_payload$class_levels
  } else if (identical(fam, "gaussian")) {
    y_train <- as.numeric(y_train_raw)
    y_scaler <- gp_fast_y_scale_fit(y_train)
    y_train_processed <- y_scaler$scaled
    class_levels <- NULL
  } else {
    class_levels <- if (exists("gp_py_ml_class_levels", mode = "function")) {
      gp_py_ml_class_levels(y_train_raw, fam)
    } else {
      gp_response_class_levels(y_train_raw, fam)
    }
    if (identical(fam, "binary") && length(class_levels) != 2L) {
      stop("Binary deep-learning support requires exactly two observed classes.", call. = FALSE)
    }
    y_train <- as.character(y_train_raw)
    y_train_processed <- as.integer(factor(y_train, levels = class_levels)) - 1L
    y_scaler <- NULL
  }

  # Build matrices
  if (isTRUE(crossval)) {
    if (is.null(tst)) stop("crossval=TRUE but 'tst' indices are NULL", call. = FALSE)
    stopifnot(is.matrix(omics_data))
    stopifnot(is.numeric(tst), all(tst >= 1), all(tst <= nrow(omics_data)))
    cv_split <- if (!is.null(cv_payload)) {
      cv_payload
    } else {
      cache_key <- gp_dl_cv_cache_key(
        omics_data = omics_data,
        y_raw = y_train_raw,
        tst = tst,
        scaling = scaling,
        centering = centering,
        response_family = fam
      )
      gp_dl_get_cached_cv_split(cache_key, function() {
        list(
          X_tr = as.matrix(omics_data[-tst, , drop = FALSE]),
          X_te = as.matrix(omics_data[tst, , drop = FALSE]),
          y_tr = y_train_processed[-tst],
          y_te = y_train_processed[tst]
        )
      })
    }
    X_tr <- cv_split$X_tr
    X_te <- cv_split$X_te
    y_tr <- cv_split$y_tr
    y_te <- cv_split$y_te
  } else {
    X_tr <- as.matrix(geno_omic_object)
    X_te <- if (!is.null(geno_omic_test_object)) as.matrix(geno_omic_test_object) else NULL
    y_tr <- y_train_processed
  }

  dl_arg_source_env <- environment()
  gp_deep_learning_current_args <- function(overrides = list()) {
    value <- function(name) {
      if (!is.null(overrides) && name %in% names(overrides)) {
        return(overrides[[name]])
      }
      get(name, envir = dl_arg_source_env, inherits = TRUE)
    }

    args_common <- compact(list(
      optimizer_name = value("optimizer_name"),
      use_amp = value("use_amp"),
      max_grad_norm = value("max_grad_norm"),
      auto_class_weights = value("auto_class_weights"),
      compile_model = value("compile_model"),
      deterministic = value("deterministic"),
      random_seed = if (!is.null(overrides) && "random_seed" %in% names(overrides)) {
        value("random_seed")
      } else {
        dl_primary_training_seed
      },
      device = value("device"),
      batch_norm = value("batch_norm"),
      validation_split = value("validation_split")
    ))

    model_key <- tolower(as.character(value("model_type") %||% model_type)[1L])
    args_model <- switch(
      model_key,
      "cnn" = compact(list(
        cnn_neurons_per_layer = value("cnn_neurons_per_layer"),
        cnn_kernel_size = value("cnn_kernel_size"),
        cnn_dense_layers = value("cnn_dense_layers"),
        cnn_use_max_pool = value("cnn_use_max_pool"),
        cnn_pool_kernel = value("cnn_pool_kernel"),
        cnn_pool_stride = value("cnn_pool_stride"),
        cnn_pool_padding = value("cnn_pool_padding"),
        cnn_separable = value("cnn_separable"),
        cnn_dilations = value("cnn_dilations"),
        cnn_use_se = value("cnn_use_se"),
        cnn_norm_type = value("cnn_norm_type"),
        cnn_pool_type = value("cnn_pool_type"),
        cnn_use_global_pool = value("cnn_use_global_pool"),
        cnn_learning_rate = value("cnn_learning_rate"),
        cnn_epochs = value("epochs"),
        cnn_batch_size = value("batch_size"),
        cnn_l2_regularizer_dp = value("cnn_l2_regularizer_dp"),
        cnn_dropout_rate = value("cnn_dropout_rate")
      )),
      "mlp" = compact(list(
        mlp_neurons_per_layer = value("mlp_neurons_per_layer"),
        mlp_final_attention = FALSE,
        mlp_attention_across_layers = value("attention_across_multiple_layers"),
        mlp_learning_rate = value("mlp_learning_rate"),
        mlp_epochs = value("epochs"),
        mlp_batch_size = value("batch_size"),
        mlp_dropout_rate = value("dropout_rate")
      )),
      "mlp_with_attention" = compact(list(
        mlp_neurons_per_layer = value("mlp_neurons_per_layer"),
        mlp_final_attention = TRUE,
        mlp_attention_across_layers = value("attention_across_multiple_layers"),
        mlp_learning_rate = value("mlp_learning_rate"),
        mlp_epochs = value("epochs"),
        mlp_batch_size = value("batch_size"),
        mlp_dropout_rate = value("dropout_rate")
      )),
      "ft_transformer" = compact(list(
        ft_d_model = value("ft_d_model"),
        ft_heads = value("ft_heads"),
        ft_layers = value("ft_layers"),
        ft_ff_mult = value("ft_ff_mult"),
        ft_dropout = value("ft_dropout"),
        ft_token_dropout = value("ft_token_dropout"),
        ft_use_cls = value("ft_use_cls"),
        ft_use_grouping = value("use_grouping"),
        ft_group_trigger = value("group_trigger"),
        ft_group_method = value("group_method"),
        ft_init_group_size = value("init_group_size"),
        ft_max_tokens = value("max_tokens"),
        ft_kmeans_batch = value("kmeans_batch"),
        ft_kmeans_iter = value("kmeans_iter"),
        ft_epochs = value("epochs"),
        ft_batch_size = value("batch_size")
      )),
      "resnet" = compact(list(
        resnet_neurons_per_block = value("resnet_neurons_per_block"),
        resnet_blocks = value("resnet_blocks"),
        resnet_learning_rate = value("resnet_learning_rate"),
        resnet_epochs = value("epochs"),
        resnet_batch_size = value("batch_size"),
        resnet_dropout_rate = value("dropout_rate")
      )),
      "saint" = compact(list(
        saint_d_model = value("saint_d_model"),
        saint_heads = value("saint_heads"),
        saint_layers = value("saint_layers"),
        saint_ff_mult = value("saint_ff_mult"),
        saint_dropout = value("saint_dropout"),
        saint_token_dropout = value("saint_token_dropout"),
        saint_use_cls = value("saint_use_cls"),
        saint_use_grouping = value("use_grouping"),
        saint_group_trigger = value("group_trigger"),
        saint_group_method = value("group_method"),
        saint_init_group_size = value("init_group_size"),
        saint_max_tokens = value("max_tokens"),
        saint_kmeans_batch = value("kmeans_batch"),
        saint_kmeans_iter = value("kmeans_iter"),
        saint_epochs = value("epochs"),
        saint_batch_size = value("batch_size")
      )),
      "tabnet" = compact(list(
        tabnet_steps = value("tabnet_steps"),
        tabnet_feature_dim = value("tabnet_feature_dim"),
        tabnet_output_dim = value("tabnet_output_dim"),
        tabnet_gamma = value("tabnet_gamma"),
        tabnet_lambda_sparse = value("tabnet_lambda_sparse"),
        tabnet_epochs = value("epochs"),
        tabnet_batch_size = value("batch_size")
      )),
      "node" = compact(list(
        node_trees = value("node_trees"),
        node_depth = value("node_depth"),
        node_epochs = value("epochs"),
        node_batch_size = value("batch_size")
      )),
      "deepfm" = compact(list(
        deepfm_k = value("deepfm_k"),
        deepfm_hidden = value("deepfm_hidden"),
        deepfm_epochs = value("epochs"),
        deepfm_batch_size = value("batch_size")
      )),
      "dcnv2" = compact(list(
        dcn_layers = value("dcn_layers"),
        dcn_hidden = value("dcn_hidden"),
        dcn_epochs = value("epochs"),
        dcn_batch_size = value("batch_size")
      )),
      "nam" = compact(list(
        nam_hidden = value("nam_hidden"),
        nam_activation = value("nam_activation"),
        nam_add_linear = value("nam_add_linear"),
        nam_l1 = value("nam_l1"),
        nam_use_grouping = value("use_grouping"),
        nam_group_trigger = value("group_trigger"),
        nam_group_method = value("group_method"),
        nam_init_group_size = value("init_group_size"),
        nam_max_tokens = value("max_tokens"),
        nam_kmeans_batch = value("kmeans_batch"),
        nam_kmeans_iter = value("kmeans_iter"),
        nam_epochs = value("epochs"),
        nam_batch_size = value("batch_size")
      )),
      "moe" = compact(list(
        moe_n_experts = value("moe_n_experts"),
        moe_expert_hidden = value("moe_expert_hidden"),
        moe_gate_hidden = value("moe_gate_hidden"),
        moe_temperature = value("moe_temperature"),
        moe_sparse_topk = value("moe_sparse_topk"),
        moe_entropy_reg = value("moe_entropy_reg"),
        moe_use_grouping = value("use_grouping"),
        moe_group_trigger = value("group_trigger"),
        moe_group_method = value("group_method"),
        moe_init_group_size = value("init_group_size"),
        moe_max_tokens = value("max_tokens"),
        moe_kmeans_batch = value("kmeans_batch"),
        moe_kmeans_iter = value("kmeans_iter"),
        moe_epochs = value("epochs"),
        moe_batch_size = value("batch_size")
      )),
      "gp_dkl" = compact(list(
        gp_use_variational = value("gp_use_variational"),
        gp_num_inducing = value("gp_num_inducing"),
        gp_feature_dim = value("gp_feature_dim"),
        gp_kernel = value("gp_kernel"),
        gp_ard = value("gp_ard"),
        gp_lr_mult = value("gp_lr_mult"),
        gp_rff_features = value("rff_features"),
        gp_rff_lengthscale = value("rff_lengthscale"),
        gp_rff_deep_hidden = value("rff_deep_hidden"),
        gp_epochs = value("epochs"),
        gp_batch_size = value("batch_size")
      )),
      stop("Unsupported model_type: ", model_key, call. = FALSE)
    )

    dl_args <- c(args_model, args_common)
    dl_args$response_family <- fam
    dl_args
  }

  if (isFALSE(crossval) && isTRUE(para_tunning)) {
    if (is.null(param_grid) || !length(param_grid)) {
      stop(paste(msg, "param_grid must be provided when tuning is enabled."), call. = FALSE)
    }
    base_names <- names(formals(deep_learning_model))
    base_args <- mget(
      base_names,
      envir = environment(),
      inherits = TRUE,
      ifnotfound = as.list(rep(list(NULL), length(base_names)))
    )
    base_args$pheno_object <- NULL
    base_args$geno_omic_object <- NULL
    base_args$geno_omic_test_object <- NULL
    base_args$response <- NULL
    base_args$gen_name <- NULL
    base_args$y <- y_train_raw
    base_args$omics_data <- X_tr
    base_args$crossval <- TRUE
    base_args$para_tunning <- FALSE
    base_args$param_grid <- NULL

    tune_res <- gp_grid_tune_cv(
      y = y_train_raw,
      response_family = fam,
      param_grid = param_grid,
      nfolds = 5L,
      random_state = random_seed %||% 123L,
      predict_fun = function(tst, params) {
        args <- utils::modifyList(base_args, params)
        args$tst <- tst
        pred <- do.call(deep_learning_model, args)
        list(pred = pred, prob = attr(pred, "probabilities"))
      },
      batch_predict_fun = function(specs) {
        cv_payloads <- lapply(specs, function(spec) {
          gp_dl_get_cached_cv_payload(
            omics_data = X_tr,
            y_raw = y_train_raw,
            tst = spec$tst,
            fam = fam,
            scaling = scaling,
            centering = centering
          )
        })
        bridge_jobs <- lapply(seq_along(specs), function(i) {
          spec <- specs[[i]]
          payload <- cv_payloads[[i]]
          list(
            id = paste(spec$param_index, spec$fold_index, sep = "_"),
            model_type = spec$params$model_type %||% model_type,
            X_train = payload$X_tr,
            y_train = payload$y_tr,
            X_test = payload$X_te,
            response_family = fam,
            dl_args = gp_deep_learning_current_args(spec$params),
            class_levels = payload$class_levels
          )
        })
        preds <- gp_dl_bridge_fit_predict_seed_batch(
          bridge_jobs,
          seed_manifest = dl_seed_manifest
        )$predictions
        Map(function(pred, payload) {
          if (identical(fam, "gaussian")) {
            pred <- revert_scaling(as.numeric(pred), payload$y_scaler)
          }
          list(pred = pred, prob = attr(pred, "probabilities"))
        }, preds, cv_payloads)
      }
    )
    tuning_summary <- tune_res
    for (nm in names(tune_res$best_params)) {
      assign(nm, tune_res$best_params[[nm]], envir = environment())
    }
    para_tunning <- FALSE
  }

  dl_args <- gp_deep_learning_current_args()

  dl_internal_calibration <- gp_dl_calibration_enabled(dl_internal_calibration)
  if (isFALSE(crossval)) {
    dl_fit_plan <- dl_seed_plan(
      n_bootstrap = n_bootstrap,
      dl_n_seeds = nrow(dl_seed_manifest),
      dl_seeds = dl_seed_manifest$training_seed,
      random_seed = random_seed,
      calibration_folds = if (dl_internal_calibration && identical(fam, "gaussian") &&
          length(y_train_raw) >= 6L) 5L else 0L,
      dl_seed_aggregation = dl_seed_aggregation
    )
    if (isTRUE(message)) {
      total_fits <- dl_fit_plan$fit_counts$model_fits[
        dl_fit_plan$fit_counts$component == "total"
      ]
      base::message(
        "DL true-prediction plan: ", n_bootstrap, " bootstrap samples x ",
        nrow(dl_seed_manifest), " training seeds, ",
        dl_fit_plan$fit_counts$data_resamples_or_folds[
          dl_fit_plan$fit_counts$component == "heldout_risk_calibration"
        ], " internal calibration folds; ",
        format(total_fits, scientific = FALSE, trim = TRUE),
        " planned model fits."
      )
    }
  }

  cv_risk_calibration <- NULL
  if (isFALSE(crossval) && dl_internal_calibration &&
      identical(fam, "gaussian") && length(y_train_raw) >= 6L) {
    seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    restore_seed <- seed_exists && is.integer(old_seed) && length(old_seed) > 1L
    gp_set_seed(random_seed %||% 123L)
    base_names <- names(formals(deep_learning_model))
    base_args <- mget(
      base_names,
      envir = environment(),
      inherits = TRUE,
      ifnotfound = as.list(rep(list(NULL), length(base_names)))
    )
    base_args$pheno_object <- NULL
    base_args$geno_omic_object <- NULL
    base_args$geno_omic_test_object <- NULL
    base_args$response <- NULL
    base_args$gen_name <- NULL
    base_args$y <- y_train_raw
    base_args$omics_data <- X_tr
    base_args$crossval <- TRUE
    base_args$para_tunning <- FALSE
    base_args$param_grid <- NULL

    cv_folds <- gp_tuning_folds(
      y = y_train_raw,
      response_family = fam,
      nfolds = 5L,
      random_state = random_seed %||% 123L
    )
    cv_pred <- rep(NA_real_, length(y_train_raw))
    cv_fold_id <- rep(NA_integer_, length(y_train_raw))
    uncertainty_target_x <- X_te %||% X_tr
    cv_target_predictions <- matrix(
      NA_real_,
      nrow = length(cv_folds),
      ncol = nrow(uncertainty_target_x)
    )
    cv_payloads <- lapply(cv_folds, function(fold_tst) {
      gp_dl_get_cached_cv_payload(
        omics_data = X_tr,
        y_raw = y_train_raw,
        tst = fold_tst,
        fam = fam,
        scaling = scaling,
        centering = centering,
        target_data = uncertainty_target_x
      )
    })
    cv_jobs <- lapply(seq_along(cv_payloads), function(i) {
      payload <- cv_payloads[[i]]
      cv_fold_id[cv_folds[[i]]] <<- i
      list(
        id = paste0("risk_", i),
        model_type = model_type,
        X_train = payload$X_tr,
        y_train = payload$y_tr,
        X_test = rbind(payload$X_te, payload$X_target),
        response_family = fam,
        dl_args = dl_args
      )
    })
    cv_pred_list <- tryCatch(
      gp_dl_bridge_fit_predict_seed_batch(
        cv_jobs,
        seed_manifest = dl_seed_manifest
      )$predictions,
      error = function(e) {
        warning("Batched DL risk calibration failed; falling back to serial folds: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(cv_pred_list)) {
      for (i in seq_along(cv_folds)) {
        values <- revert_scaling(as.numeric(cv_pred_list[[i]]), cv_payloads[[i]]$y_scaler)
        n_holdout <- length(cv_folds[[i]])
        cv_pred[cv_folds[[i]]] <- values[seq_len(n_holdout)]
        cv_target_predictions[i, ] <- values[n_holdout + seq_len(nrow(uncertainty_target_x))]
      }
    } else {
      for (i in seq_along(cv_folds)) {
        payload <- cv_payloads[[i]]
        fold_pred <- gp_dl_bridge_fit_predict_seed_serial(
          model_type = model_type,
          X_train = payload$X_tr,
          y_train = payload$y_tr,
          X_test = rbind(payload$X_te, payload$X_target),
          response_family = fam,
          dl_args = dl_args,
          class_levels = payload$class_levels %||% NULL,
          seed_manifest = dl_seed_manifest
        )$prediction
        values <- revert_scaling(as.numeric(fold_pred), payload$y_scaler)
        n_holdout <- length(cv_folds[[i]])
        cv_pred[cv_folds[[i]]] <- values[seq_len(n_holdout)]
        cv_target_predictions[i, ] <- values[n_holdout + seq_len(nrow(uncertainty_target_x))]
      }
    }
    use_local_target_uncertainty <-
      gp_dl_supports_local_target_uncertainty(model_type)
    cv_risk_calibration <- gp_ml_cv_rank_risk_calibration(
      predicted_cv = cv_pred,
      observed_y = y_train_raw,
      fold_id = cv_fold_id,
      target_fold_predictions = cv_target_predictions,
      heldout_predictors = if (use_local_target_uncertainty) {
        X_tr
      } else NULL,
      target_predictors = if (use_local_target_uncertainty) {
        uncertainty_target_x
      } else NULL,
      heldout_ids = if (use_local_target_uncertainty) {
        rownames(X_tr)
      } else NULL,
      target_ids = if (use_local_target_uncertainty) {
        rownames(uncertainty_target_x)
      } else NULL
    )
    if (restore_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  # ---------- Cross-validation short path ----------
  if (isTRUE(crossval)) {
    if (identical(fam, "gaussian") && nrow(dl_seed_manifest) == 1L) {
      pred <- gp_dl_bridge_fit_predict_gaussian_fast(
        model_type = model_type,
        X_train = X_tr,
        y_train = y_tr,
        X_test = X_te,
        dl_args = dl_args
      )
      pred <- revert_scaling(as.numeric(pred), y_scaler)
    } else {
      pred <- gp_dl_bridge_fit_predict_seed_serial(
        model_type = model_type,
        X_train = X_tr,
        y_train = y_tr,
        X_test = X_te,
        response_family = fam,
        dl_args = dl_args,
        class_levels = class_levels,
        seed_manifest = dl_seed_manifest
      )$prediction
      if (identical(fam, "gaussian")) {
        pred <- revert_scaling(as.numeric(pred), y_scaler)
      }
    }
    return(pred)
  }

  # ---------- Bootstrap path ----------
  data_label_geno <- cbind(y_tr, X_tr)
  revert_scaling_ml <- function(x, scaler) {
    if (identical(fam, "gaussian")) {
      return(as.numeric(x * scaler$std + scaler$mean))
    }
    if (identical(fam, "binary")) {
      return(pmax(0, pmin(1, as.numeric(x))))
    }
    x
  }

  if (!is.null(X_te)) {
    # Published prediction: networks trained on ALL training lines, one per
    # training seed, averaged (seed ensemble). The bootstrap refits below
    # describe uncertainty only.
    full_data_ensemble <- tryCatch(
      gp_dl_full_data_seed_ensemble(
        model_type = model_type, X_train = X_tr, y_train = y_tr, X_pred = X_te,
        response_family = fam, dl_args = dl_args,
        seeds = dl_seed_manifest$training_seed, class_levels = class_levels
      ),
      error = function(e) {
        warning("DL full-data seed ensemble failed; using the bootstrap mean: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    boot_payload <- gp_dl_bridge_bootstrap(
      model_type = model_type,
      X_train = X_tr,
      y_train = y_tr,
      X_pred = X_te,
      n_bootstrap = n_bootstrap,
      response_family = fam,
      dl_args = dl_args,
      seed = random_seed %||% 123L,
      training_seeds = dl_seed_manifest$training_seed,
      seed_aggregation = dl_seed_aggregation
    )
    dl_seed_tables <- gp_dl_seed_output_tables(
      boot_payload = boot_payload,
      ids = GID,
      response_family = fam,
      y_scaler = y_scaler,
      class_levels = class_levels,
      gen_name = gen_name
    )
    dl_seed_tables$manifest <- dl_seed_manifest
    boot_results <- list(t = boot_payload$bootstrap)
    if (identical(fam, "multiclass")) {
      boot_results$t <- t(vapply(seq_len(nrow(boot_results$t)), function(i) {
        p <- matrix(boot_results$t[i, ], nrow = nrow(X_te), ncol = length(class_levels))
        p <- pmax(0, pmin(1, p))
        as.numeric(p)
      }, numeric(nrow(X_te) * length(class_levels))))
    } else {
      boot_results$t <- gp_dl_bootstrap_response_scale(
        boot_results$t,
        response_family = fam,
        y_scaler = y_scaler
      )
    }
    result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
      boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
    )

    train_test_label <- if (!is.null(geno_omic_test_object)) {
      ifelse(rownames(geno_omic_test_object) %in% test_label, "Test", "Train")
    } else rep("Train", nrow(geno_omic_object))
    if (identical(fam, "multiclass")) {
      ml_metrics <- gp_multiclass_bootstrap_metrics(
        boot_matrix = boot_results$t,
        n_obs = if (!is.null(X_te)) nrow(X_te) else nrow(X_tr),
        class_levels = class_levels,
        train_test_label = train_test_label
      )
      pred_variances <- ml_metrics$prediction_error_var
      pred_SE <- ml_metrics$standard_error
      genetic_var <- ml_metrics$genetic_var
      prob_mean <- ml_metrics$mean_prob
      if (is.matrix(full_data_ensemble) && identical(dim(full_data_ensemble), dim(prob_mean))) {
        prob_mean <- full_data_ensemble
      }
      AI_pred_reverted <- class_levels[max.col(prob_mean, ties.method = "first")]
      lower_bound <- ml_metrics$lower_conf
      upper_bound <- ml_metrics$upper_conf
    } else {
      ml_metrics <- gp_ml_bootstrap_target_metrics(
        boot_matrix = boot_results$t,
        train_test_label = train_test_label,
        observed_y = if (identical(fam, "gaussian")) y_train else NULL
      )
      pred_variances <- ml_metrics$prediction_error_var
      pred_SE <- ml_metrics$standard_error
      reference_variance <- ml_metrics$reference_variance
      genetic_var <- reference_variance
      AI_pred_reverted <- ml_metrics$predicted_mean
      if (is.numeric(full_data_ensemble) && !is.matrix(full_data_ensemble) &&
          length(full_data_ensemble) == length(AI_pred_reverted)) {
        AI_pred_reverted <- gp_dl_bootstrap_response_scale(
          matrix(full_data_ensemble, nrow = 1L), response_family = fam, y_scaler = y_scaler
        )[1L, ]
      }
      lower_bound <- result_rel_MPIW$lower_bound
      upper_bound <- result_rel_MPIW$upper_bound
    }

    if (identical(fam, "gaussian")) {
      AI_preds <- gp_ml_gaussian_prediction_output(
        ids = GID,
        predicted_value = AI_pred_reverted,
        train_test_label = train_test_label,
        pred_se = pred_SE,
        pred_variances = pred_variances,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        uncertainty = result_rel_MPIW$Uncertainty,
        uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        reference_variance = reference_variance,
        observed_y = y_train,
        boot_results = boot_results,
        risk_calibration = cv_risk_calibration,
        predictor_matrix = X_te %||% X_tr,
        gen_name = gen_name,
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres,
        confidence_level = 0.95,
        require_heldout_calibration = TRUE
      )
    } else if (identical(fam, "binary")) {
      cls <- gp_classification_prediction_summary(
        prob = AI_pred_reverted,
        class_levels = class_levels,
        high_confidence = high_reliability_thres,
        low_confidence = low_reliability_thres
      )
      AI_preds <- data.frame(
        name = GID, Predicted_value = AI_pred_reverted, Train_Test_Label = train_test_label,
        Standard_error = pred_SE, PEV = pred_variances, Genetic_variance = rep(genetic_var, length(AI_pred_reverted)),
        lower_bound = lower_bound, upper_bound = upper_bound,
        Uncertainty = result_rel_MPIW$Uncertainty, Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        stringsAsFactors = FALSE
      )
      AI_preds <- cbind(AI_preds, cls)
    } else {
      AI_preds <- gp_multiclass_prediction_output(
        ids = GID,
        prob_mean = prob_mean,
        train_test_label = train_test_label,
        pred_se = pred_SE,
        pred_variances = pred_variances,
        genetic_var = genetic_var,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        uncertainty = result_rel_MPIW$Uncertainty,
        uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        class_levels = class_levels,
        gen_name = gen_name,
        high_confidence = high_reliability_thres,
        low_confidence = low_reliability_thres
      )
    }
    names(AI_preds)[1] <- c(gen_name)

  } else {
    boot_payload <- gp_dl_bridge_bootstrap(
      model_type = model_type,
      X_train = X_tr,
      y_train = y_tr,
      X_pred = X_tr,
      n_bootstrap = n_bootstrap,
      response_family = fam,
      dl_args = dl_args,
      seed = random_seed %||% 123L,
      training_seeds = dl_seed_manifest$training_seed,
      seed_aggregation = dl_seed_aggregation
    )
    dl_seed_tables <- gp_dl_seed_output_tables(
      boot_payload = boot_payload,
      ids = GID,
      response_family = fam,
      y_scaler = y_scaler,
      class_levels = class_levels,
      gen_name = gen_name
    )
    dl_seed_tables$manifest <- dl_seed_manifest
    boot_results <- list(t = boot_payload$bootstrap)
    if (!identical(fam, "multiclass")) {
      boot_results$t <- gp_dl_bootstrap_response_scale(
        boot_results$t,
        response_family = fam,
        y_scaler = y_scaler
      )
    }
    if (identical(fam, "multiclass")) {
      ml_metrics <- gp_multiclass_bootstrap_metrics(
        boot_matrix = boot_results$t,
        n_obs = nrow(geno_omic_object),
        class_levels = class_levels,
        train_test_label = rep("Train", nrow(geno_omic_object))
      )
      pred_variances <- ml_metrics$prediction_error_var
      pred_SE <- ml_metrics$standard_error
      genetic_var <- ml_metrics$genetic_var
      prob_mean <- ml_metrics$mean_prob
      AI_pred_reverted <- class_levels[max.col(prob_mean, ties.method = "first")]
    } else {
      ml_metrics <- gp_ml_bootstrap_target_metrics(
        boot_matrix = boot_results$t,
        train_test_label = rep("Train", nrow(geno_omic_object)),
        observed_y = if (identical(fam, "gaussian")) y_train else NULL
      )
      pred_variances <- ml_metrics$prediction_error_var
      pred_SE <- ml_metrics$standard_error
      reference_variance <- ml_metrics$reference_variance
      genetic_var <- reference_variance
      AI_pred_reverted <- ml_metrics$predicted_mean
    }
    result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
      boot_results = boot_results, CI_width_thresholds = CI_width_thresholds
    )

    train_test_label <- rep("Train", nrow(geno_omic_object))
    lower_bound <- if (identical(fam, "multiclass")) ml_metrics$lower_conf else result_rel_MPIW$lower_bound
    upper_bound <- if (identical(fam, "multiclass")) ml_metrics$upper_conf else result_rel_MPIW$upper_bound
    if (identical(fam, "gaussian")) {
      AI_preds <- gp_ml_gaussian_prediction_output(
        ids = GID,
        predicted_value = AI_pred_reverted,
        train_test_label = train_test_label,
        pred_se = pred_SE,
        pred_variances = pred_variances,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        uncertainty = result_rel_MPIW$Uncertainty,
        uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        reference_variance = reference_variance,
        observed_y = y_train,
        boot_results = boot_results,
        risk_calibration = cv_risk_calibration,
        predictor_matrix = X_tr,
        gen_name = gen_name,
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres,
        confidence_level = 0.95,
        require_heldout_calibration = TRUE
      )
    } else if (identical(fam, "binary")) {
      cls <- gp_classification_prediction_summary(
        prob = AI_pred_reverted,
        class_levels = class_levels,
        high_confidence = high_reliability_thres,
        low_confidence = low_reliability_thres
      )
      AI_preds <- data.frame(
        name = GID, Predicted_value = AI_pred_reverted, Standard_error = pred_SE,
        Train_Test_Label = train_test_label, PEV = pred_variances, Genetic_variance = rep(genetic_var, length(AI_pred_reverted)),
        lower_bound = lower_bound, upper_bound = upper_bound,
        Uncertainty = result_rel_MPIW$Uncertainty, Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        stringsAsFactors = FALSE
      )
      AI_preds <- cbind(AI_preds, cls)
    } else {
      AI_preds <- gp_multiclass_prediction_output(
        ids = GID,
        prob_mean = prob_mean,
        train_test_label = train_test_label,
        pred_se = pred_SE,
        pred_variances = pred_variances,
        genetic_var = genetic_var,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        uncertainty = result_rel_MPIW$Uncertainty,
        uncertainty_remarks = result_rel_MPIW$reliability_remarks,
        class_levels = class_levels,
        gen_name = gen_name,
        high_confidence = high_reliability_thres,
        low_confidence = low_reliability_thres
      )
    }
    names(AI_preds)[1] <- c(gen_name)
  }

  if (identical(fam, "gaussian") && !dl_internal_calibration) {
    AI_preds[["Prediction_uncertainty_source"]] <-
      "unavailable_internal_calibration_disabled_by_user"
  }

  # Parameter summary (simple)
  model_para <- data.frame(
    stat = c(
      "model_type", "epochs", "batch_size", "n_bootstrap",
      "dl_n_seeds", "dl_seeds", "dl_seed_aggregation",
      "dl_internal_calibration",
      "dl_planned_model_fits"
    ),
    summary = c(
      model_type,
      epochs,
      batch_size,
      n_bootstrap,
      nrow(dl_seed_manifest),
      paste(dl_seed_manifest$training_seed, collapse = ","),
      dl_seed_aggregation,
      dl_internal_calibration,
      if (is.null(dl_fit_plan)) NA_real_ else
        dl_fit_plan$fit_counts$model_fits[
          dl_fit_plan$fit_counts$component == "total"
        ]
    ),
    stringsAsFactors = FALSE
  )
  if (!is.null(tuning_summary)) {
    tune_rows <- data.frame(
      stat = c("tune_metric", "tune_best_score"),
      summary = c(tuning_summary$metric, as.character(tuning_summary$best_score)),
      stringsAsFactors = FALSE
    )
    if (length(tuning_summary$best_params)) {
      best_param_rows <- data.frame(
        stat = paste0("best_", names(tuning_summary$best_params)),
        summary = vapply(tuning_summary$best_params, function(x) paste(x, collapse = ","), character(1)),
        stringsAsFactors = FALSE
      )
      tune_rows <- rbind(tune_rows, best_param_rows)
    }
    model_para <- rbind(model_para, tune_rows)
  }

  variance_components <- if (identical(fam, "gaussian")) {
    gp_ml_gaussian_variance_components(
      boot_matrix = boot_results$t,
      train_test_label = train_test_label,
      observed_y = y_train,
      predicted_value = AI_pred_reverted,
      reference_variance = reference_variance,
      prediction_error_var = AI_preds[["PEV"]],
      marginal_prediction_mse = attr(
        AI_preds,
        "uncertainty_calibration"
      )$calibration_mean_pev
    )
  } else {
    gp_predictive_model_variance_components(
      AI_preds,
      response_family = fam
    )
  }

  list(
    model_parameters = model_para,
    predicted_values = AI_preds,
    variance_components = variance_components,
    Variance_components = variance_components,
    dl_seed_manifest = dl_seed_tables$manifest,
    dl_computation_plan = dl_fit_plan$fit_counts,
    dl_seed_predictions = dl_seed_tables$predictions,
    dl_seed_variability = dl_seed_tables$variability,
    diagnostic_plots = if (identical(fam, "gaussian") &&
        !is.null(gp_ml_heldout_residual_summary(
          cv_risk_calibration,
          confidence_level = 0.95
        ))) {
      diagnostic_plot_true_prediction(
        boot_results = boot_results, GID_names = GID, CI_width_thresholds = CI_width_thresholds,
        predictions = AI_pred_reverted, standard_errors = pred_SE, prediction_error_var = pred_variances,
        reference_variance = reference_variance, lower_bound = lower_bound, upper_bound = upper_bound,
        observed_y = y_train, train_test_label = train_test_label,
        confidence_level = 0.95, risk_calibration = cv_risk_calibration,
        require_heldout_calibration = TRUE, model_for_CI_cal = "ML",
        response_family = fam,
        high_reliability_thres = high_reliability_thres, low_reliability_thres = low_reliability_thres,
        system_database = system_database
      )
    } else if (identical(fam, "binary")) {
      diagnostic_plot_true_prediction(
        boot_results = boot_results, GID_names = GID, CI_width_thresholds = CI_width_thresholds,
        predictions = AI_pred_reverted, standard_errors = pred_SE, prediction_error_var = pred_variances,
        genetic_var = genetic_var, confidence_level = 0.95, model_for_CI_cal = "ML",
        response_family = fam,
        high_reliability_thres = high_reliability_thres, low_reliability_thres = low_reliability_thres,
        system_database = system_database
      )
    } else {
      NULL
    }
  )
}

gp_met_deep_learning_model <- function(model_type = "mlp",
                                       pheno_object = NULL,
                                       response = NULL,
                                       gen_name = NULL,
                                       heter_groups = NULL,
                                       response_family = "gaussian",
                                       gmatrix = NULL,
                                       omic1_kernel = NULL,
                                       omic2_kernel = NULL,
                                       omic3_kernel = NULL,
                                       kernel_list = NULL,
                                       omics_kernel_label = list(
                                         omic1_kernel = NULL,
                                         omic2_kernel = NULL,
                                         omic3_kernel = NULL
                                       ),
                                       met_kernel_var_explained = 0.95,
                                       met_kernel_min_ev = 1e-8,
                                       met_kernel_max_pcs = NULL,
                                       compile_model = FALSE,
                                       deterministic = TRUE,
                                       random_seed = 123,
                                       device = NULL,
                                       use_amp = FALSE,
                                       batch_norm = FALSE,
                                       validation_split = 0.2,
                                       epochs = 30,
                                       batch_size = 64,
                                       mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                       mlp_learning_rate = 1e-3,
                                       ft_d_model = 64L,
                                       ft_heads = 4L,
                                       ft_layers = 2L,
                                       ft_ff_mult = 2L,
                                       ft_dropout = 0.1,
                                       ft_token_dropout = 0.1,
                                       ft_use_cls = TRUE,
                                       saint_d_model = 64L,
                                       saint_heads = 4L,
                                       saint_layers = 2L,
                                       saint_ff_mult = 2L,
                                       saint_dropout = 0.1,
                                       saint_token_dropout = 0.1,
                                       saint_use_cls = TRUE,
                                       tabnet_steps = 3L,
                                       tabnet_feature_dim = 16L,
                                       tabnet_output_dim = 16L,
                                       tabnet_gamma = 1.3,
                                       tabnet_lambda_sparse = 1e-4,
                                       moe_n_experts = 4L,
                                       moe_expert_hidden = as.integer(c(128L, 64L)),
                                       moe_gate_hidden = 128L,
                                       moe_temperature = 1.0,
                                       moe_sparse_topk = NA,
                                       moe_entropy_reg = 0.0,
                                       dropout = 0.2,
                                       l2_weight_decay = 1e-4,
                                       final_attention = FALSE,
                                       attention_across_multiple_layers = FALSE,
                                       optimizer_name = "adam",
                                       max_grad_norm = NULL,
                                       heteroscedastic = FALSE,
                                        n_bootstrap = 30,
                                        internal_cv_nfolds = 5L,
                                        dl_internal_calibration = TRUE,
                                        CI_width_thresholds = c(0.33, 0.66),
                                       high_reliability_thres = 0.9,
                                       low_reliability_thres = 0.5,
                                       system_database = FALSE) {
  dl_internal_calibration <- gp_dl_calibration_enabled(dl_internal_calibration)
  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  if (!fam %in% c("gaussian", "binary", "multiclass")) {
    stop(model_type, " MET path currently supports gaussian, binary, and multiclass traits only.", call. = FALSE)
  }
  feat_res <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )
  feature_table <- feat_res$feature_table
  if (is.null(feature_table) || !nrow(feature_table)) {
    stop("mlp MET path requires at least one GRM/kernel source.", call. = FALSE)
  }

  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  rownames(x_all) <- long_dat$.met_row_id

  train_idx <- which(!is.na(long_dat[[response]]))
  test_idx <- which(is.na(long_dat[[response]]))
  if (!length(train_idx)) {
    stop("mlp MET path requires non-missing training observations.", call. = FALSE)
  }
  if (identical(fam, "gaussian")) {
    n_bootstrap <- suppressWarnings(as.integer(n_bootstrap[[1L]]))
    if (!is.finite(n_bootstrap) || n_bootstrap < 2L) {
      stop(model_type, " MET Gaussian uncertainty requires n_bootstrap >= 2.",
           call. = FALSE)
    }
  }
  class_levels <- if (exists("gp_py_ml_class_levels", mode = "function")) {
    gp_py_ml_class_levels(long_dat[[response]][train_idx], fam)
  } else if (identical(fam, "gaussian")) {
    NULL
  } else {
    sort(unique(as.character(stats::na.omit(long_dat[[response]][train_idx]))))
  }

  pheno_train <- data.frame(
    MET_ID = long_dat$.met_row_id[train_idx],
    y = long_dat[[response]][train_idx],
    stringsAsFactors = FALSE
  )
  x_prep <- caret::preProcess(x_all[train_idx, , drop = FALSE], method = c("center", "scale"))
  train_x <- stats::predict(x_prep, x_all[train_idx, , drop = FALSE])
  all_x <- stats::predict(x_prep, x_all)
  y_train_raw <- long_dat[[response]][train_idx]
  if (identical(fam, "gaussian")) {
    y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train_raw)), method = c("center", "scale"))
    y_train_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train_raw)))[, 1]
  } else {
    y_scaler <- NULL
    y_train_fit <- gp_dl_bridge_prepare_target(y_train_raw, fam, class_levels = class_levels)$y
  }

  dl_args <- compact(list(
    response_family = fam,
    compile_model = compile_model,
    deterministic = deterministic,
    random_seed = random_seed,
    device = device,
    use_amp = use_amp,
    batch_norm = batch_norm,
    validation_split = validation_split,
    epochs = epochs,
    batch_size = batch_size,
    mlp_neurons_per_layer = mlp_neurons_per_layer,
    mlp_learning_rate = mlp_learning_rate,
    ft_d_model = ft_d_model,
    ft_heads = ft_heads,
    ft_layers = ft_layers,
    ft_ff_mult = ft_ff_mult,
    ft_dropout = ft_dropout,
    ft_token_dropout = ft_token_dropout,
    ft_use_cls = ft_use_cls,
    saint_d_model = saint_d_model,
    saint_heads = saint_heads,
    saint_layers = saint_layers,
    saint_ff_mult = saint_ff_mult,
    saint_dropout = saint_dropout,
    saint_token_dropout = saint_token_dropout,
    saint_use_cls = saint_use_cls,
    tabnet_steps = tabnet_steps,
    tabnet_feature_dim = tabnet_feature_dim,
    tabnet_output_dim = tabnet_output_dim,
    tabnet_gamma = tabnet_gamma,
    tabnet_lambda_sparse = tabnet_lambda_sparse,
    moe_n_experts = moe_n_experts,
    moe_expert_hidden = moe_expert_hidden,
    moe_gate_hidden = moe_gate_hidden,
    moe_temperature = moe_temperature,
    moe_sparse_topk = moe_sparse_topk,
    moe_entropy_reg = moe_entropy_reg,
    dropout = dropout,
    dropout_rate = dropout,
    l2_weight_decay = l2_weight_decay,
    l2_regularizer_dp = l2_weight_decay,
    final_attention = final_attention,
    attention_across_multiple_layers = attention_across_multiple_layers,
    optimizer_name = optimizer_name,
    max_grad_norm = max_grad_norm,
    heteroscedastic = heteroscedastic
  ))

  pred_bridge <- gp_dl_bridge_fit_predict(
    model_type = tolower(model_type),
    X_train = as.matrix(train_x),
    y_train = y_train_fit,
    X_test = as.matrix(all_x),
    response_family = fam,
    dl_args = dl_args,
    class_levels = class_levels
  )
  if (identical(fam, "gaussian")) {
    pred_scaled <- as.numeric(pred_bridge)
    pred_all <- as.numeric(pred_scaled * y_scaler$std + y_scaler$mean)
    bootstrap_payload <- gp_dl_bridge_bootstrap(
      model_type = tolower(model_type),
      X_train = as.matrix(train_x),
      y_train = y_train_fit,
      X_pred = as.matrix(all_x),
      n_bootstrap = n_bootstrap,
      response_family = fam,
      dl_args = dl_args,
      seed = random_seed %||% 123L
    )
    bootstrap_matrix <- as.matrix(bootstrap_payload$bootstrap)
    bootstrap_matrix <- bootstrap_matrix * as.numeric(y_scaler$std) +
      as.numeric(y_scaler$mean)
    heldout_calibration <- if (dl_internal_calibration) gp_met_grouped_crossfit_calibration(
      y = long_dat[[response]],
      observed_idx = train_idx,
      group_id = long_dat[[gen_name]],
      environment = long_dat[[heter_groups]],
      n_targets = nrow(x_all),
      predictor_matrix = x_all,
      row_id = long_dat$.met_row_id,
      fit_predict = function(fit_idx, prediction_idx, fold_index) {
        fold_x_prep <- caret::preProcess(
          x_all[fit_idx, , drop = FALSE],
          method = c("center", "scale")
        )
        fold_train_x <- stats::predict(
          fold_x_prep,
          x_all[fit_idx, , drop = FALSE]
        )
        fold_pred_x <- stats::predict(
          fold_x_prep,
          x_all[prediction_idx, , drop = FALSE]
        )
        fold_y_raw <- as.numeric(long_dat[[response]][fit_idx])
        fold_y_mean <- mean(fold_y_raw)
        fold_y_sd <- stats::sd(fold_y_raw)
        if (!is.finite(fold_y_sd) || fold_y_sd <= 0) fold_y_sd <- 1
        fold_args <- utils::modifyList(
          dl_args,
          list(random_seed = as.integer((random_seed %||% 123L) + fold_index))
        )
        fold_prediction <- gp_dl_bridge_fit_predict(
          model_type = tolower(model_type),
          X_train = as.matrix(fold_train_x),
          y_train = (fold_y_raw - fold_y_mean) / fold_y_sd,
          X_test = as.matrix(fold_pred_x),
          response_family = fam,
          dl_args = fold_args,
          class_levels = NULL
        )
        as.numeric(fold_prediction) * fold_y_sd + fold_y_mean
      },
      nfolds = internal_cv_nfolds,
      seed = random_seed %||% 123L
    ) else NULL
    train_test_label <- ifelse(
      seq_len(nrow(long_dat)) %in% test_idx,
      "Test",
      "Train"
    )
    uncertainty <- gp_met_bootstrap_uncertainty_columns(
      predicted_value = pred_all,
      bootstrap_matrix = bootstrap_matrix,
      observed_y = long_dat[[response]],
      train_test_label = train_test_label,
      heldout_calibration = heldout_calibration,
      n_bootstrap = n_bootstrap,
      confidence_level = 0.95
    )
    if (!dl_internal_calibration) {
      uncertainty[["Prediction_uncertainty_source"]] <-
        "unavailable_internal_calibration_disabled_by_user"
    }
    pred <- data.frame(
      MET_ID = long_dat$.met_row_id,
      Predicted_value = pred_all,
      Train_Test_Label = train_test_label,
      Observed_value = as.numeric(long_dat[[response]]),
      stringsAsFactors = FALSE
    )
    pred <- cbind(pred, uncertainty)
  } else if (identical(fam, "binary")) {
    pred_prob <- as.numeric(pred_bridge)
    pred_prob <- pmax(0, pmin(1, pred_prob))
    cls <- gp_classification_prediction_summary(
      prob = pred_prob,
      class_levels = class_levels,
      high_confidence = high_reliability_thres,
      low_confidence = low_reliability_thres
    )
    prob_df <- data.frame(
      Prob_1 = 1 - pred_prob,
      Prob_2 = pred_prob,
      stringsAsFactors = FALSE
    )
    names(prob_df) <- gp_prob_column_names(class_levels)
    pred <- data.frame(
      MET_ID = long_dat$.met_row_id,
      Predicted_value = pred_prob,
      Train_Test_Label = ifelse(seq_len(nrow(long_dat)) %in% test_idx, "Test", "Train"),
      cls,
      stringsAsFactors = FALSE
    )
    pred <- cbind(pred, prob_df)
  } else {
    prob_mat <- attr(pred_bridge, "probabilities")
    if (is.null(prob_mat)) {
      prob_mat <- matrix(NA_real_, nrow = length(pred_bridge), ncol = length(class_levels))
    }
    prob_mat <- as.matrix(prob_mat)
    colnames(prob_mat) <- class_levels
    cls <- gp_classification_prediction_summary(
      prob = prob_mat,
      class_levels = class_levels,
      high_confidence = high_reliability_thres,
      low_confidence = low_reliability_thres
    )
    prob_df <- as.data.frame(prob_mat, stringsAsFactors = FALSE)
    names(prob_df) <- gp_prob_column_names(class_levels)
    pred <- data.frame(
      MET_ID = long_dat$.met_row_id,
      Predicted_value = cls$Predicted_class,
      Train_Test_Label = ifelse(seq_len(nrow(long_dat)) %in% test_idx, "Test", "Train"),
      cls,
      stringsAsFactors = FALSE
    )
    pred <- cbind(pred, prob_df)
  }
  meta <- long_dat[, c(".met_row_id", gen_name, heter_groups), drop = FALSE]
  names(meta)[1] <- "MET_ID"
  pred <- merge(pred, meta, by = "MET_ID", all.x = TRUE, sort = FALSE)
  pred <- pred[, c(gen_name, heter_groups, setdiff(names(pred), c("MET_ID", gen_name, heter_groups))), drop = FALSE]

  model_para <- data.frame(
    stat = c(
      "model_type", "epochs", "batch_size", "met_feature_source",
      "met_bootstrap", "met_bootstrap_replicates",
      "met_uncertainty_calibration", "dl_internal_calibration",
      "met_internal_calibration_requested_folds"
    ),
    summary = c(
      model_type, as.character(epochs), as.character(batch_size),
      "kernel_eigen", as.character(identical(fam, "gaussian")),
      if (identical(fam, "gaussian")) as.character(n_bootstrap) else "0",
       if (identical(fam, "gaussian")) {
         if (dl_internal_calibration) "CV1_by_GID_grouped_internal_crossfit" else
           "disabled_by_user"
      } else {
        "not_applicable"
      },
      as.character(dl_internal_calibration),
      if (identical(fam, "gaussian") && dl_internal_calibration) {
        as.character(internal_cv_nfolds)
      } else "0"
    ),
    stringsAsFactors = FALSE
  )
  model_para <- rbind(
    model_para,
    data.frame(stat = "met_kernel_var_explained", summary = as.character(met_kernel_var_explained), stringsAsFactors = FALSE),
    data.frame(stat = "met_num_feature_blocks", summary = as.character(nrow(feat_res$feature_summary)), stringsAsFactors = FALSE),
    data.frame(stat = "met_total_kernel_pcs", summary = as.character(ncol(feature_table)), stringsAsFactors = FALSE),
    data.frame(stat = "response_family", summary = fam, stringsAsFactors = FALSE)
  )
  met_kernel_names <- as.character(feat_res$feature_summary$source)
  met_kernel_names[met_kernel_names == "GRM"] <- "gmatrix"
  model_para <- rbind(
    model_para,
    gp_multi_kernel_parameter_rows(
      kernel_names = met_kernel_names,
      strategy = "concatenated_kernel_eigenfeature_blocks"
    )
  )

  list(
    model_parameters = model_para,
    predicted_values = pred,
    diagnostic_plots = NULL,
    met_long_data = long_dat,
    met_feature_summary = feat_res$feature_summary
  )
}

gp_met_dl_cv_predict <- function(model_type = "mlp",
                                 pheno_data,
                                 response,
                                 gen_name,
                                 heter_groups,
                                 response_family = "gaussian",
                                 gmatrix = NULL,
                                 omic1_kernel = NULL,
                                 omic2_kernel = NULL,
                                 omic3_kernel = NULL,
                                 kernel_list = NULL,
                                 omics_kernel_label = list(
                                   omic1_kernel = NULL,
                                   omic2_kernel = NULL,
                                   omic3_kernel = NULL
                                 ),
                                 tst,
                                 model_params = list(),
                                 met_kernel_var_explained = 0.95,
                                 met_kernel_min_ev = 1e-8,
                                 met_kernel_max_pcs = NULL) {
  feat_res <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )
  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feat_res$feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  fam <- gp_resolve_response_family(response_family, y = long_dat[[response]])
  y_all <- if (identical(fam, "gaussian")) as.numeric(long_dat[[response]]) else long_dat[[response]]
  train_idx <- setdiff(which(!is.na(y_all)), tst)
  if (!length(train_idx)) {
    stop("No observed training rows available for MET mlp CV.", call. = FALSE)
  }
  deep_learning_model(
    y = y_all,
    omics_data = as.matrix(x_all),
    tst = tst,
    crossval = TRUE,
    response_family = fam,
    model_type = model_type,
    compile_model = model_params$compile_model %||% FALSE,
    deterministic = model_params$deterministic %||% TRUE,
    random_seed = model_params$random_seed %||% 123L,
    device = model_params$device %||% NULL,
    use_amp = model_params$use_amp %||% FALSE,
    batch_norm = model_params$batch_norm %||% FALSE,
    validation_split = model_params$validation_split %||% 0.2,
    epochs = model_params$epochs %||% 30L,
    batch_size = model_params$batch_size %||% 64L,
    mlp_neurons_per_layer = model_params$mlp_neurons_per_layer %||% as.integer(c(128L, 64L)),
    mlp_learning_rate = model_params$mlp_learning_rate %||% 1e-3,
    ft_d_model = model_params$ft_d_model %||% 64L,
    ft_heads = model_params$ft_heads %||% 4L,
    ft_layers = model_params$ft_layers %||% 2L,
    ft_ff_mult = model_params$ft_ff_mult %||% 2L,
    ft_dropout = model_params$ft_dropout %||% 0.1,
    ft_token_dropout = model_params$ft_token_dropout %||% 0.1,
    ft_use_cls = model_params$ft_use_cls %||% TRUE,
    saint_d_model = model_params$saint_d_model %||% 64L,
    saint_heads = model_params$saint_heads %||% 4L,
    saint_layers = model_params$saint_layers %||% 2L,
    saint_ff_mult = model_params$saint_ff_mult %||% 2L,
    saint_dropout = model_params$saint_dropout %||% 0.1,
    saint_token_dropout = model_params$saint_token_dropout %||% 0.1,
    saint_use_cls = model_params$saint_use_cls %||% TRUE,
    tabnet_steps = model_params$tabnet_steps %||% 3L,
    tabnet_feature_dim = model_params$tabnet_feature_dim %||% 16L,
    tabnet_output_dim = model_params$tabnet_output_dim %||% 16L,
    tabnet_gamma = model_params$tabnet_gamma %||% 1.3,
    tabnet_lambda_sparse = model_params$tabnet_lambda_sparse %||% 1e-4,
    moe_n_experts = model_params$moe_n_experts %||% 4L,
    moe_expert_hidden = model_params$moe_expert_hidden %||% as.integer(c(128L, 64L)),
    moe_gate_hidden = model_params$moe_gate_hidden %||% 128L,
    moe_temperature = model_params$moe_temperature %||% 1.0,
    moe_sparse_topk = model_params$moe_sparse_topk %||% NA,
    moe_entropy_reg = model_params$moe_entropy_reg %||% 0.0,
    dropout = model_params$dropout %||% 0.2,
    dropout_rate = model_params$dropout_rate %||% (model_params$dropout %||% 0.2),
    l2_weight_decay = model_params$l2_weight_decay %||% 1e-4,
    l2_regularizer_dp = model_params$l2_regularizer_dp %||% (model_params$l2_weight_decay %||% 1e-4),
    final_attention = model_params$final_attention %||% FALSE,
    attention_across_multiple_layers = model_params$attention_across_multiple_layers %||% FALSE,
    optimizer_name = model_params$optimizer_name %||% "adam",
    max_grad_norm = model_params$max_grad_norm %||% NULL,
    heteroscedastic = model_params$heteroscedastic %||% FALSE
  )
}

gp_met_dl_cv_predict_batch <- function(model_type = "mlp",
                                       pheno_data,
                                       response,
                                       gen_name,
                                       heter_groups,
                                       response_family = "gaussian",
                                       gmatrix = NULL,
                                       omic1_kernel = NULL,
                                       omic2_kernel = NULL,
                                       omic3_kernel = NULL,
                                       kernel_list = NULL,
                                       omics_kernel_label = list(
                                         omic1_kernel = NULL,
                                         omic2_kernel = NULL,
                                         omic3_kernel = NULL
                                       ),
                                       folds,
                                       model_params = list(),
                                       met_kernel_var_explained = 0.95,
                                       met_kernel_min_ev = 1e-8,
                                       met_kernel_max_pcs = NULL) {
  if (is.null(folds) || !length(folds)) {
    return(list())
  }
  feat_res <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )
  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feat_res$feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  fam <- gp_resolve_response_family(response_family, y = long_dat[[response]])
  y_all <- if (identical(fam, "gaussian")) as.numeric(long_dat[[response]]) else long_dat[[response]]
  observed_idx <- which(!is.na(y_all))
  class_levels <- if (exists("gp_py_ml_class_levels", mode = "function")) {
    gp_py_ml_class_levels(y_all[observed_idx], fam)
  } else {
    NULL
  }
  dl_args <- model_params %||% list()
  dl_args$response_family <- fam
  jobs <- lapply(seq_along(folds), function(i) {
    tst <- as.integer(folds[[i]])
    train_idx <- setdiff(observed_idx, tst)
    if (!length(train_idx)) {
      stop("No observed training rows available for MET ", model_type, " CV.", call. = FALSE)
    }
    list(
      id = i,
      model_type = tolower(as.character(model_type)[1L]),
      X_train = x_all[train_idx, , drop = FALSE],
      y_train = y_all[train_idx],
      X_test = x_all[tst, , drop = FALSE],
      response_family = fam,
      dl_args = dl_args,
      class_levels = class_levels
    )
  })
  gp_dl_bridge_fit_predict_batch(jobs)
}

AI_mlp_MET <- function(...) gp_met_deep_learning_model(model_type = "mlp", ...)
AI_ft_transformer_MET <- function(...) gp_met_deep_learning_model(model_type = "ft_transformer", ...)
AI_saint_MET <- function(...) gp_met_deep_learning_model(model_type = "saint", ...)
AI_tabnet_MET <- function(...) gp_met_deep_learning_model(model_type = "tabnet", ...)
AI_moe_MET <- function(...) gp_met_deep_learning_model(model_type = "moe", ...)

gp_met_mlp_cv_predict <- function(...) gp_met_dl_cv_predict(model_type = "mlp", ...)
gp_met_ft_transformer_cv_predict <- function(...) gp_met_dl_cv_predict(model_type = "ft_transformer", ...)
gp_met_saint_cv_predict <- function(...) gp_met_dl_cv_predict(model_type = "saint", ...)
gp_met_tabnet_cv_predict <- function(...) gp_met_dl_cv_predict(model_type = "tabnet", ...)
gp_met_moe_cv_predict <- function(...) gp_met_dl_cv_predict(model_type = "moe", ...)

gp_multitrait_dl_supported_models <- function() c("mlp", "ft_transformer")

gp_scale_multitrait_matrix <- function(y_mat) {
  centers <- apply(y_mat, 2, function(v) mean(v, na.rm = TRUE))
  scales <- apply(y_mat, 2, function(v) stats::sd(v, na.rm = TRUE))
  centers[!is.finite(centers)] <- 0
  scales[!is.finite(scales) | scales <= 0] <- 1
  y_scaled <- sweep(y_mat, 2, centers, FUN = "-")
  y_scaled <- sweep(y_scaled, 2, scales, FUN = "/")
  list(y_scaled = y_scaled, center = centers, scale = scales)
}

gp_revert_multitrait_matrix <- function(y_mat, center, scale) {
  out <- sweep(y_mat, 2, scale, FUN = "*")
  sweep(out, 2, center, FUN = "+")
}

gp_multitrait_dl_summary_statistics <- function(pred_long, response, model_type) {
  counts <- do.call(rbind, lapply(response, function(tr) {
    dat <- pred_long[pred_long$Trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      observed_count = sum(dat$Train_Test_Label == "Train", na.rm = TRUE),
      predicted_count = sum(dat$Train_Test_Label == "Test", na.rm = TRUE),
      total_count = nrow(dat),
      stringsAsFactors = FALSE
    )
  }))
  data.frame(
    stat = c(
      "mode",
      "response_family",
      "model_type",
      "n_traits",
      "traits",
      paste0("observed_count_", counts$trait),
      paste0("predicted_count_", counts$trait)
    ),
    summary = c(
      "multi_trait_dl",
      "gaussian",
      model_type,
      length(response),
      paste(response, collapse = ", "),
      counts$observed_count,
      counts$predicted_count
    ),
    stringsAsFactors = FALSE
  )
}

gp_multitrait_dl_diagnostic_plot <- function(pred_long) {
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~Trait, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Multi-trait Gaussian DL observed vs predicted",
      #subtitle = "Observed cells are used for masked training; missing cells are final predictions",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_multitrait_local_margin <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  if (!n) {
    return(numeric())
  }
  if (n == 1L) {
    return(Inf)
  }
  ord <- order(values)
  sorted <- values[ord]
  out_sorted <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    left_gap <- if (i > 1L) abs(sorted[i] - sorted[i - 1L]) else Inf
    right_gap <- if (i < n) abs(sorted[i + 1L] - sorted[i]) else Inf
    out_sorted[i] <- min(left_gap, right_gap, na.rm = TRUE)
  }
  out <- rep(NA_real_, n)
  out[ord] <- out_sorted
  out
}

gp_multitrait_risk_remark <- function(risk, low_cut = 1 / 3, high_cut = 2 / 3) {
  ifelse(
    risk >= high_cut,
    "High Risk",
    ifelse(risk <= low_cut, "Low Risk", "Moderate Risk")
  )
}

gp_multitrait_confidence_remark <- function(risk_remark) {
  ifelse(
    risk_remark == "Low Risk",
    "High Confidence",
    ifelse(risk_remark == "High Risk", "Low Confidence", "Moderate Confidence")
  )
}

gp_multitrait_dl_true_prediction_risk <- function(pred_long,
                                                  cv_pred_long,
                                                  gen_name = "GID") {
  if (is.null(pred_long) || !nrow(pred_long) || is.null(cv_pred_long) || !nrow(cv_pred_long)) {
    return(pred_long)
  }

  cv_obs <- cv_pred_long[!is.na(cv_pred_long$Observed_value), , drop = FALSE]
  if (!nrow(cv_obs)) {
    return(pred_long)
  }

  cv_obs$absolute_error <- abs(cv_obs$Predicted_value - cv_obs$Observed_value)
  cv_obs$pred_margin <- ave(cv_obs$Predicted_value, cv_obs$Trait, FUN = gp_multitrait_local_margin)

  rank_summary <- gp_multitrait_dl_rank_risk_summary(cv_obs)
  rank_rows <- rank_summary$combined_risk
  if (is.null(rank_rows) || !nrow(rank_rows)) {
    return(pred_long)
  }

  cv_obs <- merge(
    cv_obs,
    rank_rows[, c("GID", "Trait", "fold", "rep", "model", "rank_shift", "rank_instability_risk", "rank_instability_risk_remarks"), drop = FALSE],
    by = c("GID", "Trait", "fold", "rep", "model"),
    all.x = TRUE,
    sort = FALSE
  )

  train_row_risk <- do.call(rbind, lapply(split(cv_obs, list(cv_obs$GID, cv_obs$Trait), drop = TRUE), function(dat) {
    mean_risk <- mean(dat$rank_instability_risk, na.rm = TRUE)
    data.frame(
      GID = dat$GID[[1]],
      Trait = dat$Trait[[1]],
      train_rank_risk = mean_risk,
      train_prediction_confidence = 1 - mean_risk,
      stringsAsFactors = FALSE
    )
  }))
  rownames(train_row_risk) <- NULL

  trait_profiles <- lapply(split(cv_obs, cv_obs$Trait), function(dat) {
    low_cut <- stats::quantile(dat$rank_instability_risk, probs = 1 / 3, na.rm = TRUE, names = FALSE, type = 7)
    high_cut <- stats::quantile(dat$rank_instability_risk, probs = 2 / 3, na.rm = TRUE, names = FALSE, type = 7)
    margin_lo <- stats::quantile(dat$pred_margin, probs = 1 / 3, na.rm = TRUE, names = FALSE, type = 7)
    margin_hi <- stats::quantile(dat$pred_margin, probs = 2 / 3, na.rm = TRUE, names = FALSE, type = 7)
    list(
      reference = dat[, c("Predicted_value", "rank_instability_risk", "pred_margin"), drop = FALSE],
      low_cut = as.numeric(low_cut),
      high_cut = as.numeric(high_cut),
      margin_lo = as.numeric(margin_lo),
      margin_hi = as.numeric(margin_hi)
    )
  })

  out <- pred_long
  out$Prediction_confidence <- NA_real_
  out$Prediction_confidence_remarks <- NA_character_
  out$Prediction_risk_score <- NA_real_
  out$Prediction_risk_remarks <- NA_character_
  out$Prediction_risk_basis <- NA_character_
  out$Rank_instability_risk <- NA_real_
  out$Rank_instability_risk_remarks <- NA_character_
  out$Rank_instability_risk_basis <- NA_character_

  train_idx <- which(out$Train_Test_Label == "Train")
  if (length(train_idx)) {
    train_merge <- merge(
      out[train_idx, c(gen_name, "Trait"), drop = FALSE],
      train_row_risk,
      by.x = c(gen_name, "Trait"),
      by.y = c("GID", "Trait"),
      all.x = TRUE,
      sort = FALSE
    )
    out$Prediction_risk_score[train_idx] <- train_merge$train_rank_risk
    out$Rank_instability_risk[train_idx] <- train_merge$train_rank_risk
    out$Prediction_confidence[train_idx] <- train_merge$train_prediction_confidence
    out$Prediction_risk_basis[train_idx] <- "Heldout_Training_CV"
    out$Rank_instability_risk_basis[train_idx] <- "Heldout_Training_CV"
  }

  test_idx <- which(out$Train_Test_Label == "Test")
  if (length(test_idx)) {
    for (tr in unique(out$Trait[test_idx])) {
      idx <- test_idx[out$Trait[test_idx] == tr]
      prof <- trait_profiles[[tr]]
      if (is.null(prof) || is.null(prof$reference) || !nrow(prof$reference)) {
        next
      }
      ref_pred <- prof$reference$Predicted_value
      ref_risk <- prof$reference$rank_instability_risk
      pred_vals <- out$Predicted_value[idx]
      pred_margin <- gp_multitrait_local_margin(pred_vals)
      k <- min(7L, length(ref_pred))
      test_risk <- vapply(seq_along(idx), function(i) {
        dist <- abs(ref_pred - pred_vals[i])
        ord <- order(dist, na.last = NA)
        ord <- ord[seq_len(min(k, length(ord)))]
        w <- 1 / (dist[ord] + 1e-6)
        base_risk <- stats::weighted.mean(ref_risk[ord], w = w, na.rm = TRUE)
        if (is.finite(pred_margin[i])) {
          if (pred_margin[i] <= prof$margin_lo) {
            base_risk <- base_risk + 0.08
          } else if (pred_margin[i] >= prof$margin_hi) {
            base_risk <- base_risk - 0.05
          }
        }
        pmax(0, pmin(1, base_risk))
      }, numeric(1))

      out$Prediction_risk_score[idx] <- test_risk
      out$Rank_instability_risk[idx] <- test_risk
      out$Prediction_confidence[idx] <- 1 - test_risk
      out$Prediction_risk_basis[idx] <- "Calibrated_Final_Prediction_CV"
      out$Rank_instability_risk_basis[idx] <- "Calibrated_Final_Prediction_CV"
    }
  }

  for (tr in unique(out$Trait)) {
    idx <- which(out$Trait == tr & is.finite(out$Prediction_risk_score))
    prof <- trait_profiles[[tr]]
    low_cut <- if (!is.null(prof)) prof$low_cut else 1 / 3
    high_cut <- if (!is.null(prof)) prof$high_cut else 2 / 3
    remarks <- gp_multitrait_risk_remark(out$Prediction_risk_score[idx], low_cut = low_cut, high_cut = high_cut)
    out$Prediction_risk_remarks[idx] <- remarks
    out$Rank_instability_risk_remarks[idx] <- remarks
    out$Prediction_confidence_remarks[idx] <- gp_multitrait_confidence_remark(remarks)
  }

  out
}

gp_multitrait_dl_gaussian_model <- function(model_type = c("mlp", "ft_transformer"),
                                            pheno_object,
                                            response,
                                            geno_omic_object,
                                            gen_name,
                                            message = FALSE,
                                            scaling = TRUE,
                                            centering = FALSE,
                                            compile_model = FALSE,
                                            deterministic = TRUE,
                                            random_seed = 123L,
                                            device = NULL,
                                            use_amp = FALSE,
                                            batch_norm = FALSE,
                                            validation_split = 0.2,
                                            epochs = 30L,
                                            batch_size = 64L,
                                            mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                            mlp_learning_rate = 1e-3,
                                            ft_d_model = 64L,
                                            ft_heads = 4L,
                                            ft_layers = 2L,
                                            ft_ff_mult = 2L,
                                            ft_dropout = 0.1,
                                            ft_token_dropout = 0.1,
                                            ft_use_cls = TRUE,
                                            dropout = 0.1,
                                            l2_weight_decay = 1e-4,
                                            final_attention = FALSE,
                                            attention_across_multiple_layers = FALSE,
                                            optimizer_name = "adam",
                                            max_grad_norm = NULL,
                                             internal_cv_nfolds = 5L,
                                             internal_cv_replication = 3L,
                                             dl_internal_calibration = TRUE,
                                             system_database = FALSE) {
  model_type <- match.arg(model_type)
  dl_internal_calibration <- gp_dl_calibration_enabled(dl_internal_calibration)
  if (length(response) < 2L) {
    stop("Multi-trait DL requires at least two response columns.", call. = FALSE)
  }
  if (is.null(geno_omic_object) || !nrow(geno_omic_object)) {
    stop("Multi-trait DL requires a processed geno/omic feature matrix.", call. = FALSE)
  }
  ph <- pheno_object[, c(gen_name, response), drop = FALSE]
  keep_ids <- intersect(as.character(ph[[gen_name]]), rownames(geno_omic_object))
  if (!length(keep_ids)) {
    stop("No overlapping genotype IDs between phenotypes and processed features.", call. = FALSE)
  }
  ph <- ph[match(keep_ids, as.character(ph[[gen_name]])), , drop = FALSE]
  x_all <- as.matrix(geno_omic_object[keep_ids, , drop = FALSE])
  storage.mode(x_all) <- "double"
  y_mat <- as.matrix(ph[, response, drop = FALSE])
  storage.mode(y_mat) <- "double"
  train_rows <- stats::complete.cases(y_mat)
  if (!any(train_rows)) {
    stop("Multi-trait DL requires at least one genotype with all requested traits observed.", call. = FALSE)
  }
  if (isTRUE(scaling) || isTRUE(centering)) {
    pp_methods <- c(if (isTRUE(centering)) "center", if (isTRUE(scaling)) "scale")
    x_prep <- caret::preProcess(x_all[train_rows, , drop = FALSE], method = pp_methods)
    x_all <- as.matrix(stats::predict(x_prep, x_all))
  }

  y_train_mat <- y_mat[train_rows, , drop = FALSE]
  x_train <- x_all[train_rows, , drop = FALSE]
  ph_train <- ph[train_rows, , drop = FALSE]
  y_scale_res <- gp_scale_multitrait_matrix(y_train_mat)
  y_fit <- y_scale_res$y_scaled

  dl_args <- compact(list(
    compile_model = compile_model,
    deterministic = deterministic,
    random_seed = random_seed,
    device = device,
    use_amp = use_amp,
    batch_norm = batch_norm,
    validation_split = validation_split,
    epochs = epochs,
    batch_size = batch_size,
    mlp_neurons_per_layer = mlp_neurons_per_layer,
    mlp_learning_rate = mlp_learning_rate,
    ft_d_model = ft_d_model,
    ft_heads = ft_heads,
    ft_layers = ft_layers,
    ft_ff_mult = ft_ff_mult,
    ft_dropout = ft_dropout,
    ft_token_dropout = ft_token_dropout,
    ft_use_cls = ft_use_cls,
    dropout = dropout,
    dropout_rate = dropout,
    l2_weight_decay = l2_weight_decay,
    l2_regularizer_dp = l2_weight_decay,
    final_attention = final_attention,
    attention_across_multiple_layers = attention_across_multiple_layers,
    optimizer_name = optimizer_name,
    max_grad_norm = max_grad_norm,
    multitask = TRUE
  ))

  pred_scaled <- gp_dl_bridge_fit_predict(
    model_type = model_type,
    X_train = x_train,
    y_train = y_fit,
    X_test = x_all,
    response_family = "multitask_regression",
    dl_args = dl_args
  )
  pred_scaled <- as.matrix(pred_scaled)
  colnames(pred_scaled) <- response
  pred_mat <- gp_revert_multitrait_matrix(pred_scaled, y_scale_res$center, y_scale_res$scale)
  colnames(pred_mat) <- response

  ids <- as.character(ph[[gen_name]])
  pred_long <- do.call(rbind, lapply(seq_along(response), function(j) {
    tr <- response[[j]]
    obs <- y_mat[, j]
    data.frame(
      GID = ids,
      Trait = tr,
      Observed_value = obs,
      Predicted_value = pred_mat[, j],
      Train_Test_Label = ifelse(is.na(obs), "Test", "Train"),
      stringsAsFactors = FALSE
    )
  }))
  names(pred_long)[1] <- gen_name

  n_cv_folds <- max(2L, min(as.integer(internal_cv_nfolds %||% 5L), nrow(y_train_mat)))
  cv_calibration_error <- NULL
  cv_calibration <- if (dl_internal_calibration) tryCatch(
    gp_multitrait_dl_gaussian_cv(
      model_type = model_type,
      # Supply every prediction target to the resampling layer. The CV helper
      # removes rows with all traits missing when it forms training folds, but
      # retains them in its all-target prediction view so true Test GIDs can
      # receive target-specific repeated-prediction SE/PEV.
      pheno_object = ph,
      response = response,
      geno_omic_object = geno_omic_object[as.character(ph[[gen_name]]), , drop = FALSE],
      gen_name = gen_name,
      cross_validation_meth = "K-Folds",
      nfolds = n_cv_folds,
      replication = max(1L, as.integer(internal_cv_replication %||% 3L)),
      eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
      scaling = scaling,
      centering = centering,
      compile_model = compile_model,
      deterministic = deterministic,
      random_seed = random_seed,
      device = device,
      use_amp = use_amp,
      batch_norm = batch_norm,
      validation_split = validation_split,
      epochs = epochs,
      batch_size = batch_size,
      mlp_neurons_per_layer = mlp_neurons_per_layer,
      mlp_learning_rate = mlp_learning_rate,
      ft_d_model = ft_d_model,
      ft_heads = ft_heads,
      ft_layers = ft_layers,
      ft_ff_mult = ft_ff_mult,
      ft_dropout = ft_dropout,
      ft_token_dropout = ft_token_dropout,
      ft_use_cls = ft_use_cls,
      dropout = dropout,
      l2_weight_decay = l2_weight_decay,
      final_attention = final_attention,
      attention_across_multiple_layers = attention_across_multiple_layers,
      optimizer_name = optimizer_name,
      max_grad_norm = max_grad_norm
    ),
    error = function(e) {
      cv_calibration_error <<- conditionMessage(e)
      warning(
        "Multi-trait DL uncertainty calibration failed: ",
        cv_calibration_error,
        call. = FALSE
      )
      NULL
    }
  ) else NULL
  cv_pred_long <- NULL
  target_pred_long <- NULL
  if (!is.null(cv_calibration$cv_results) && length(cv_calibration$cv_results)) {
    cv_pred_long <- do.call(rbind, lapply(cv_calibration$cv_results, function(x) x$ypred_cv_Reps_all))
    target_pred_long <- do.call(rbind, lapply(
      cv_calibration$cv_results,
      function(x) x$target_prediction_Reps_all
    ))
    pred_long <- gp_multitrait_dl_true_prediction_risk(
      pred_long = pred_long,
      cv_pred_long = cv_pred_long,
      gen_name = gen_name
    )
    pred_long <- gp_ml_apply_cv_uncertainty(
      pred = pred_long,
      cv_pred = cv_pred_long,
      id_col = gen_name,
      group_col = "Trait",
      target_pred = target_pred_long
    )
  } else {
    n_pred <- nrow(pred_long)
    pred_long[["Standard_error"]] <- rep(NA_real_, n_pred)
    pred_long[["PEV"]] <- rep(NA_real_, n_pred)
    pred_long[["Prediction_error_variance"]] <- rep(NA_real_, n_pred)
    pred_long[["Prediction_uncertainty_source"]] <- rep(
      if (dl_internal_calibration) "unavailable_internal_cv_calibration_failed" else
        "unavailable_internal_calibration_disabled_by_user",
      n_pred
    )
    pred_long[["PEV_basis"]] <- rep(
      paste(
        gp_ml_unavailable_target_uncertainty_estimand(),
        if (!dl_internal_calibration) {
          "; internal calibration disabled by user"
        } else if (!is.null(cv_calibration_error) && nzchar(cv_calibration_error)) {
          paste0("; internal calibration error: ", cv_calibration_error)
        } else {
          "; internal calibration returned no fold predictions"
        }
      ),
      n_pred
    )
    pred_long[["Reliability"]] <- rep(NA_real_, n_pred)
    pred_long[["Reliability_variance_input"]] <- rep(NA_real_, n_pred)
    pred_long[["Reliability_basis"]] <- rep(
      paste(
        "Reliability unavailable because target-specific predictive PEV",
        "was not produced; not genetic reliability"
      ),
      n_pred
    )
  }

  pred_wide <- data.frame(
    stats::setNames(list(ids), gen_name),
    as.data.frame(pred_mat, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )
  obs_counts <- data.frame(
    trait = response,
    observed_count = colSums(!is.na(y_mat)),
    predicted_count = colSums(is.na(y_mat)),
    stringsAsFactors = FALSE
  )
  model_para <- data.frame(
    stat = c(
      "mode", "response_family", "model_type", "n_traits", "traits",
      "internal_cv_calibration_status", "internal_cv_calibration_error",
      "dl_internal_calibration", "internal_cv_requested_model_fits"
    ),
    summary = c(
      "multi_trait_dl", "gaussian", model_type, length(response),
      paste(response, collapse = ", "),
      if (!dl_internal_calibration) "disabled_by_user" else
        if (!is.null(cv_pred_long)) "available" else "unavailable",
      cv_calibration_error %||% NA_character_,
      as.character(dl_internal_calibration),
      if (dl_internal_calibration) {
        as.character(n_cv_folds * max(1L, as.integer(internal_cv_replication %||% 3L)))
      } else "0"
    ),
    stringsAsFactors = FALSE
  )
  variance_components <- gp_predictive_model_variance_components(
    pred_long,
    response_family = "gaussian"
  )
  predictive_matrices <- gp_ml_predictive_matrix_bundle(
    pred = pred_long,
    id_col = gen_name,
    dimension_col = "Trait",
    error_pred = cv_pred_long,
    random_seed = random_seed
  )

  c(list(
    model_parameters = model_para,
    predicted_values = pred_long,
    multitrait_prediction_wide = pred_wide,
    multitrait_trait_counts = obs_counts,
    trained_model = NULL,
    diagnostic_plots = gp_multitrait_dl_diagnostic_plot(pred_long),
    variance_components = variance_components,
    Variance_components = variance_components,
    prediction_uncertainty_source = if (!dl_internal_calibration) {
      "unavailable_internal_calibration_disabled_by_user"
    } else if (!is.null(cv_pred_long)) "heldout_internal_cv" else "unavailable",
    predictive_covariance_basis = if (!dl_internal_calibration) {
      "descriptive predictions; heldout prediction errors unavailable; not genetic or residual covariance"
    } else {
      "descriptive predictions and heldout prediction errors; not genetic or residual covariance"
    }
  ), predictive_matrices)
}

gp_multitrait_dl_cv_plot <- function(pred_long) {
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = factor(fold))
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~Trait, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Multi-trait Gaussian DL CV observed vs predicted",
      subtitle = "Held-out genotype predictions by fold",
      x = "Observed value",
      y = "Predicted value",
      color = "Fold"
    )
}

gp_multitrait_dl_uncertainty_summary <- function(pred_long) {
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    empty_combined <- data.frame(
      GID = character(),
      Trait = character(),
      Observed_value = numeric(),
      Predicted_value = numeric(),
      fold = integer(),
      rep = integer(),
      model = character(),
      residual = numeric(),
      absolute_error = numeric(),
      squared_error = numeric(),
      stringsAsFactors = FALSE
    )
    empty_agg <- data.frame(
      trait = character(),
      model = character(),
      mean_absolute_error = numeric(),
      root_mean_squared_error = numeric(),
      residual_sd = numeric(),
      observed_sd = numeric(),
      relative_rmse = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(combined_uncertainty = empty_combined, aggregated_uncertainty = empty_agg))
  }
  obs$residual <- obs$Predicted_value - obs$Observed_value
  obs$absolute_error <- abs(obs$residual)
  obs$squared_error <- obs$residual^2
  names(obs)[names(obs) == "Trait"] <- "trait"
  names(obs)[names(obs) == "Observed_value"] <- "y"
  names(obs)[names(obs) == "Predicted_value"] <- "yhat"
  agg <- do.call(rbind, lapply(split(obs, list(obs$trait, obs$model), drop = TRUE), function(dat) {
    obs_sd <- stats::sd(dat$y, na.rm = TRUE)
    rmse <- sqrt(mean(dat$squared_error, na.rm = TRUE))
    data.frame(
      trait = dat$trait[[1]],
      model = dat$model[[1]],
      mean_absolute_error = mean(dat$absolute_error, na.rm = TRUE),
      root_mean_squared_error = rmse,
      residual_sd = stats::sd(dat$residual, na.rm = TRUE),
      observed_sd = obs_sd,
      relative_rmse = if (is.finite(obs_sd) && obs_sd > 0) rmse / obs_sd else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL
  list(combined_uncertainty = obs, aggregated_uncertainty = agg)
}

gp_multitrait_dl_rank_risk_summary <- function(pred_long) {
  pred_long <- as.data.frame(pred_long, stringsAsFactors = FALSE)
  if (!"GID" %in% names(pred_long)) {
    non_id_cols <- c(
      "Trait", "model", "Observed_value", "Predicted_value",
      "Train_Test_Label", "fold", "rep", "cv_role",
      "Prediction_SE", "Prediction_SD", "Prediction_Lower", "Prediction_Upper"
    )
    id_candidates <- setdiff(names(pred_long), non_id_cols)
    if (length(id_candidates)) {
      names(pred_long)[names(pred_long) == id_candidates[[1L]]] <- "GID"
    } else {
      pred_long$GID <- as.character(seq_len(nrow(pred_long)))
    }
  }
  if (!"model" %in% names(pred_long)) {
    pred_long$model <- "unknown"
  }
  if (!"fold" %in% names(pred_long)) {
    pred_long$fold <- 1L
  }
  if (!"rep" %in% names(pred_long)) {
    pred_long$rep <- 1L
  }
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    empty_combined <- data.frame(
      GID = character(),
      Trait = character(),
      fold = integer(),
      rep = integer(),
      model = character(),
      rank_shift = numeric(),
      rank_instability_risk = numeric(),
      rank_instability_risk_remarks = character(),
      stringsAsFactors = FALSE
    )
    empty_agg <- data.frame(
      trait = character(),
      model = character(),
      mean_rank_shift = numeric(),
      mean_rank_instability_risk = numeric(),
      low_risk_percentage = numeric(),
      moderate_risk_percentage = numeric(),
      high_risk_percentage = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(combined_risk = empty_combined, aggregated_risk = empty_agg))
  }

  risk_rows <- do.call(rbind, lapply(split(obs, list(obs[["Trait"]], obs[["fold"]], obs[["rep"]]), drop = TRUE), function(dat) {
    if (!nrow(dat)) return(NULL)
    obs_rank <- rank(-dat[["Observed_value"]], ties.method = "average")
    pred_rank <- rank(-dat[["Predicted_value"]], ties.method = "average")
    rank_shift <- abs(obs_rank - pred_rank)
    if (length(rank_shift) == 1L) {
      risk <- 0.5
      remarks <- "Moderate Risk"
    } else {
      q <- stats::quantile(rank_shift, probs = c(1/3, 2/3), na.rm = TRUE, type = 7)
      denom <- max(rank_shift, na.rm = TRUE)
      risk <- if (is.finite(denom) && denom > 0) rank_shift / denom else rep(0, length(rank_shift))
      remarks <- ifelse(
        rank_shift <= q[[1]], "Low Risk",
        ifelse(rank_shift >= q[[2]], "High Risk", "Moderate Risk")
      )
    }
    data.frame(
      GID = as.character(dat[["GID"]]),
      Trait = as.character(dat[["Trait"]]),
      fold = dat[["fold"]],
      rep = dat[["rep"]],
      model = as.character(dat[["model"]]),
      rank_shift = rank_shift,
      rank_instability_risk = risk,
      rank_instability_risk_remarks = remarks,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(risk_rows) || !nrow(risk_rows)) {
    return(list(
      combined_risk = data.frame(),
      aggregated_risk = data.frame()
    ))
  }
  agg <- do.call(rbind, lapply(split(risk_rows, list(risk_rows[["Trait"]], risk_rows[["model"]]), drop = TRUE), function(dat) {
    data.frame(
      trait = dat[["Trait"]][[1]],
      model = dat[["model"]][[1]],
      mean_rank_shift = mean(dat[["rank_shift"]], na.rm = TRUE),
      mean_rank_instability_risk = mean(dat[["rank_instability_risk"]], na.rm = TRUE),
      low_risk_percentage = mean(dat[["rank_instability_risk_remarks"]] == "Low Risk", na.rm = TRUE) * 100,
      moderate_risk_percentage = mean(dat[["rank_instability_risk_remarks"]] == "Moderate Risk", na.rm = TRUE) * 100,
      high_risk_percentage = mean(dat[["rank_instability_risk_remarks"]] == "High Risk", na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL
  list(combined_risk = risk_rows, aggregated_risk = agg)
}

gp_multitrait_dl_risk_plot <- function(risk_summary) {
  if (is.null(risk_summary) || !is.list(risk_summary) || is.null(risk_summary$aggregated_risk) || !nrow(risk_summary$aggregated_risk)) {
    return(NULL)
  }
  ggplot2::ggplot(
    risk_summary$aggregated_risk,
    ggplot2::aes(x = model, y = mean_rank_instability_risk, fill = model)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~trait, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Multi-trait Gaussian DL mean rank instability risk",
      subtitle = "Lower values indicate lower held-out rank displacement within trait",
      x = "Model",
      y = "Mean rank instability risk"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

gp_multitrait_dl_metric_table <- function(pred_long,
                                          eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                          model_type = "mlp",
                                          rep = 1L) {
  out <- list()
  idx <- 1L
  for (tr in unique(pred_long$Trait)) {
    dat_tr <- pred_long[pred_long$Trait == tr & !is.na(pred_long$Observed_value), , drop = FALSE]
    if (!nrow(dat_tr)) next
    for (fd in sort(unique(dat_tr$fold))) {
      dat_fd <- dat_tr[dat_tr$fold == fd, , drop = FALSE]
      if (!nrow(dat_fd)) next
      for (metric in eval_metrics) {
        val <- tryCatch(
          evaluation_metrics(
            y_observed = dat_fd$Observed_value,
            y_predicted = dat_fd$Predicted_value,
            eval_metrics = metric,
            response_family = "gaussian"
          ),
          error = function(e) NA_real_
        )
        out[[idx]] <- data.frame(
          trait = tr,
          model = model_type,
          rep = rep,
          fold = fd,
          metric = metric,
          value = as.numeric(val),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (!length(out)) {
    return(data.frame(
      trait = character(),
      model = character(),
      rep = integer(),
      fold = integer(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

gp_multitrait_dl_cv_process <- function(cv_results) {
  pred_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, function(x) x$ypred_cv_Reps_all))
  if (length(pred_parts)) {
    pred_all <- do.call(rbind, pred_parts)
    counts <- stats::aggregate(
      Predicted_value ~ Trait + Train_Test_Label,
      data = pred_all,
      FUN = length
    )
    names(counts) <- c("Trait", "Train_Test_Label", "count")
  } else {
    pred_all <- data.frame(
      Trait = character(),
      Train_Test_Label = character(),
      Predicted_value = numeric(),
      stringsAsFactors = FALSE
    )
    counts <- data.frame(
      Trait = character(),
      Train_Test_Label = character(),
      count = integer(),
      stringsAsFactors = FALSE
    )
  }
  metric_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, function(x) x$eval_metrics_reps))
  if (length(metric_parts)) {
    metrics_all <- do.call(rbind, metric_parts)
    if (!is.null(metrics_all) && nrow(metrics_all) > 0) {
      metric_summary <- tryCatch(
        stats::aggregate(
          value ~ trait + model + metric,
          data = metrics_all,
          FUN = function(v) mean(v, na.rm = TRUE)
        ),
        error = function(e) {
          data.frame(
            trait = character(),
            model = character(),
            metric = character(),
            value = numeric(),
            stringsAsFactors = FALSE
          )
        }
      )
    } else {
      metrics_all <- data.frame(
        trait = character(),
        model = character(),
        rep = integer(),
        fold = integer(),
        metric = character(),
        value = numeric(),
        stringsAsFactors = FALSE
      )
      metric_summary <- data.frame(
        trait = character(),
        model = character(),
        metric = character(),
        value = numeric(),
        stringsAsFactors = FALSE
      )
    }
  } else {
    metrics_all <- data.frame(
      trait = character(),
      model = character(),
      rep = integer(),
      fold = integer(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
    metric_summary <- data.frame(
      trait = character(),
      model = character(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  uncertainty_summary <- gp_multitrait_dl_uncertainty_summary(pred_all)
  rank_risk_summary <- gp_multitrait_dl_rank_risk_summary(pred_all)
  list(
    multitrait_metric_summary = metric_summary,
    multitrait_prediction_counts = counts,
    multitrait_cv_plots = gp_multitrait_dl_cv_plot(pred_all),
    multitrait_uncertainty_summaries = uncertainty_summary,
    multitrait_rank_risk_summaries = rank_risk_summary,
    multitrait_rank_risk_plot = gp_multitrait_dl_risk_plot(rank_risk_summary)
  )
}

gp_multitrait_dl_retry_fit <- function(fit_call,
                                       max_attempts = 2L,
                                       context = "multi-trait DL resampling fit") {
  max_attempts <- max(1L, as.integer(max_attempts)[1L])
  errors <- character()
  for (attempt in seq_len(max_attempts)) {
    value <- tryCatch(
      fit_call(),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      }
    )
    if (!is.null(value)) {
      return(value)
    }
    if (attempt < max_attempts) {
      warning(
        context,
        " failed on attempt ", attempt,
        "; retrying once in a fresh backend process.",
        call. = FALSE
      )
    }
  }
  stop(
    context,
    " failed after ", max_attempts, " attempts: ",
    paste(unique(errors), collapse = " | "),
    call. = FALSE
  )
}

gp_multitrait_dl_gaussian_cv <- function(model_type = c("mlp", "ft_transformer"),
                                         pheno_object,
                                         response,
                                         geno_omic_object,
                                         gen_name,
                                         cross_validation_meth = "K-Folds",
                                         nfolds = 5L,
                                         replication = 1L,
                                         eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                         scaling = TRUE,
                                         centering = FALSE,
                                         compile_model = FALSE,
                                         deterministic = TRUE,
                                         random_seed = 123L,
                                         device = NULL,
                                         use_amp = FALSE,
                                         batch_norm = FALSE,
                                         validation_split = 0.2,
                                         epochs = 30L,
                                         batch_size = 64L,
                                         mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                         mlp_learning_rate = 1e-3,
                                         ft_d_model = 64L,
                                         ft_heads = 4L,
                                         ft_layers = 2L,
                                         ft_ff_mult = 2L,
                                         ft_dropout = 0.1,
                                         ft_token_dropout = 0.1,
                                         ft_use_cls = TRUE,
                                         dropout = 0.1,
                                         l2_weight_decay = 1e-4,
                                         final_attention = FALSE,
                                         attention_across_multiple_layers = FALSE,
                                         optimizer_name = "adam",
                                         max_grad_norm = NULL,
                                         feature_score_metadata = NULL,
                                         feature_k = NULL,
                                         feature_scoring_cv = "fixed",
                                         feature_scoring_model = "Ridge_Regression",
                                         feature_source_map = NULL,
                                         feature_ridge_lambda = 1,
                                         feature_bayes_nIter = 1500L,
                                         feature_bayes_burnIn = 500L,
                                         feature_bayes_thin = 5L,
                                         ntree = 500L,
                                         mtry = NULL,
                                         nodesize = NULL,
                                         rf_n_jobs = 1L) {
  model_type <- match.arg(model_type)
  cv_token <- normalize_cv_token(cross_validation_meth)
  if (!cv_token %in% c("k_folds", "repeated_k_folds")) {
    stop("Multi-trait DL CV currently supports K-Folds only.", call. = FALSE)
  }
  ph <- pheno_object[, c(gen_name, response), drop = FALSE]
  keep_ids <- intersect(as.character(ph[[gen_name]]), rownames(geno_omic_object))
  if (!length(keep_ids)) {
    stop("No overlapping genotype IDs between phenotypes and processed features.", call. = FALSE)
  }
  target_ph <- ph[match(keep_ids, as.character(ph[[gen_name]])), , drop = FALSE]
  target_x_all0 <- as.matrix(geno_omic_object[keep_ids, , drop = FALSE])
  storage.mode(target_x_all0) <- "double"
  target_y_mat0 <- as.matrix(target_ph[, response, drop = FALSE])
  storage.mode(target_y_mat0) <- "double"
  target_gids <- as.character(target_ph[[gen_name]])
  ph <- target_ph
  x_all0 <- target_x_all0
  y_mat0 <- target_y_mat0
  keep_rows <- rowSums(!is.na(y_mat0)) > 0
  ph <- ph[keep_rows, , drop = FALSE]
  x_all0 <- x_all0[keep_rows, , drop = FALSE]
  y_mat0 <- y_mat0[keep_rows, , drop = FALSE]
  gids <- as.character(ph[[gen_name]])
  feature_active <- !is.null(feature_score_metadata) && !is.null(feature_k)
  feature_scoring_cv <- gp_feature_normalize_cv_policy(feature_scoring_cv, allow_both = FALSE)
  cv_results <- vector("list", length = replication)

  for (rep_idx in seq_len(replication)) {
    feature_records <- list()
    folds <- gp_tuning_folds(
      y = seq_len(nrow(y_mat0)),
      response_family = "gaussian",
      nfolds = nfolds,
      random_state = as.integer((random_seed %||% 123L) + rep_idx - 1L)
    )
    target_blocks <- vector("list", length = length(folds))
    pred_blocks <- vector("list", length = length(folds))
    for (fold_idx in seq_along(folds)) {
      test_idx <- sort(unique(as.integer(folds[[fold_idx]])))
      train_idx <- setdiff(seq_len(nrow(y_mat0)), test_idx)
      x_fold <- x_all0
      x_target_fold <- target_x_all0
      if (isTRUE(feature_active)) {
        feature_view <- gp_feature_multitrait_view(
          predictor_data = x_all0,
          pheno_data = ph,
          response = response,
          gen_name = gen_name,
          k = feature_k,
          selection_mode = feature_scoring_cv,
          feature_score_metadata = feature_score_metadata,
          test_rows = test_idx,
          scoring_model = feature_scoring_model,
          seed = as.integer((random_seed %||% 123L) + rep_idx * 10000L + fold_idx),
          source_block = feature_source_map,
          response_family = "gaussian",
          replication = rep_idx,
          fold = fold_idx,
          ridge_lambda = feature_ridge_lambda,
          bayes_nIter = feature_bayes_nIter,
          bayes_burnIn = feature_bayes_burnIn,
          bayes_thin = feature_bayes_thin,
          ntree = ntree,
          mtry = mtry,
          nodesize = nodesize,
          rf_n_jobs = rf_n_jobs
        )
        x_fold <- feature_view$predictor_data
        selected_columns <- intersect(colnames(x_fold), colnames(target_x_all0))
        if (length(selected_columns) != ncol(x_fold)) {
          stop("Could not apply fold-selected DL features to all prediction targets.", call. = FALSE)
        }
        x_target_fold <- target_x_all0[, selected_columns, drop = FALSE]
        feature_records[[length(feature_records) + 1L]] <- feature_view$metadata
      }
      x_train <- x_fold[train_idx, , drop = FALSE]
      x_target <- x_target_fold
      if (isTRUE(scaling) || isTRUE(centering)) {
        x_prep <- gp_fast_center_scale_fit(
          x_train,
          center = isTRUE(centering),
          scale = isTRUE(scaling)
        )
        x_train <- as.matrix(x_prep$data)
        x_target <- as.matrix(gp_fast_center_scale_apply(x_target, x_prep))
      }
      y_train_mat <- y_mat0[train_idx, , drop = FALSE]
      y_scale_res <- gp_scale_multitrait_matrix(y_train_mat)
      y_train_fit <- y_scale_res$y_scaled

      dl_args <- compact(list(
        compile_model = compile_model,
        deterministic = deterministic,
        random_seed = as.integer((random_seed %||% 123L) + rep_idx - 1L),
        device = device,
        use_amp = use_amp,
        batch_norm = batch_norm,
        validation_split = validation_split,
        epochs = epochs,
        batch_size = batch_size,
        mlp_neurons_per_layer = mlp_neurons_per_layer,
        mlp_learning_rate = mlp_learning_rate,
        ft_d_model = ft_d_model,
        ft_heads = ft_heads,
        ft_layers = ft_layers,
        ft_ff_mult = ft_ff_mult,
        ft_dropout = ft_dropout,
        ft_token_dropout = ft_token_dropout,
        ft_use_cls = ft_use_cls,
        dropout = dropout,
        dropout_rate = dropout,
        l2_weight_decay = l2_weight_decay,
        l2_regularizer_dp = l2_weight_decay,
        final_attention = final_attention,
        attention_across_multiple_layers = attention_across_multiple_layers,
        optimizer_name = optimizer_name,
        max_grad_norm = max_grad_norm,
        multitask = TRUE
      ))

      pred_scaled <- gp_multitrait_dl_retry_fit(
        fit_call = function() gp_dl_bridge_fit_predict(
          model_type = model_type,
          X_train = x_train,
          y_train = y_train_fit,
          X_test = x_target,
          response_family = "multitask_regression",
          dl_args = dl_args
        ),
        max_attempts = 2L,
        context = paste(
          "Multi-trait", model_type,
          "resampling rep", rep_idx,
          "fold", fold_idx
        )
      )
      pred_scaled <- as.matrix(pred_scaled)
      colnames(pred_scaled) <- response
      target_pred_mat <- gp_revert_multitrait_matrix(pred_scaled, y_scale_res$center, y_scale_res$scale)
      colnames(target_pred_mat) <- response
      evaluation_positions <- match(gids[test_idx], target_gids)

      pred_blocks[[fold_idx]] <- do.call(rbind, lapply(seq_along(response), function(j) {
        tr <- response[[j]]
        data.frame(
          GID = gids[test_idx],
          Trait = tr,
          model = model_type,
          Observed_value = y_mat0[test_idx, j],
          Predicted_value = target_pred_mat[evaluation_positions, j],
          Train_Test_Label = "Test",
          fold = fold_idx,
          rep = rep_idx,
          cv_role = "test",
          stringsAsFactors = FALSE
        )
      }))
      target_blocks[[fold_idx]] <- do.call(rbind, lapply(seq_along(response), function(j) {
        tr <- response[[j]]
        data.frame(
          GID = target_gids,
          Trait = tr,
          model = model_type,
          Observed_value = target_y_mat0[, j],
          Predicted_value = target_pred_mat[, j],
          Train_Test_Label = ifelse(is.na(target_y_mat0[, j]), "Test", "Train"),
          fold = fold_idx,
          rep = rep_idx,
          cv_role = "target_resampling",
          stringsAsFactors = FALSE
        )
      }))
      if (isTRUE(feature_active)) {
        pred_blocks[[fold_idx]][["feature_k"]] <- as.integer(feature_k)
        pred_blocks[[fold_idx]][["feature_scoring_model"]] <- feature_scoring_model
        pred_blocks[[fold_idx]][["feature_scoring_cv"]] <- feature_scoring_cv
      }
    }

    pred_long <- do.call(rbind, pred_blocks)
    names(pred_long)[1] <- gen_name
    target_pred_long <- do.call(rbind, target_blocks)
    names(target_pred_long)[1] <- gen_name
    metric_table <- gp_multitrait_dl_metric_table(
      pred_long = pred_long,
      eval_metrics = eval_metrics,
      model_type = model_type,
      rep = rep_idx
    )
    if (isTRUE(feature_active) && nrow(metric_table)) {
      metric_table[["feature_k"]] <- as.integer(feature_k)
      metric_table[["feature_scoring_model"]] <- feature_scoring_model
      metric_table[["feature_scoring_cv"]] <- feature_scoring_cv
    }
    cv_results[[rep_idx]] <- list(
      trait = paste(response, collapse = ","),
      model = model_type,
      rep = rep_idx,
      eval_metrics_reps = metric_table,
      ypred_cv_Reps_all = pred_long,
      target_prediction_Reps_all = target_pred_long,
      yprob_cv_Reps_all = NULL,
      feature_selection_metadata = gp_feature_bind_metadata(feature_records)
    )
  }

  processed <- gp_multitrait_dl_cv_process(cv_results)
  processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(cv_results, `[[`, "feature_selection_metadata")
  )
  list(cv_results = cv_results, cv_results_processed = processed)
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
