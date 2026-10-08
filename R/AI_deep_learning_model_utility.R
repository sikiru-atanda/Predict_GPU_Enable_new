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
#   # --- Legacy -> current arg renames ---
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

  # --- Legacy -> current arg renames ---
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

  # OPTIONAL: auto-dedupe (keep last) if you'd rather not hard-stop:
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
# - For heteroscedastic reg: returns mean/variance from the model's second channel.
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

gp_fast_center_scale_fit <- function(x, center = TRUE, scale = TRUE) {
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  x_center <- if (isTRUE(center)) colMeans(x_mat, na.rm = TRUE) else rep(0, ncol(x_mat))
  x_scale <- if (isTRUE(scale)) apply(x_mat, 2, stats::sd, na.rm = TRUE) else rep(1, ncol(x_mat))
  x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  x_scaled <- sweep(sweep(x_mat, 2, x_center, "-"), 2, x_scale, "/")
  list(data = x_scaled, center = x_center, scale = x_scale)
}

gp_fast_center_scale_apply <- function(x, fit) {
  x_mat <- as.matrix(x)
  storage.mode(x_mat) <- "double"
  sweep(sweep(x_mat, 2, fit$center, "-"), 2, fit$scale, "/")
}

gp_fast_y_scale_fit <- function(y) {
  y_num <- as.numeric(y)
  mu <- mean(y_num, na.rm = TRUE)
  sdv <- stats::sd(y_num, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) {
    sdv <- 1
  }
  list(
    scaled = (y_num - mu) / sdv,
    mean = mu,
    std = sdv
  )
}

gp_dl_bridge_fit_params <- function(model_type, dl_args = list()) {
  model_key <- tolower(as.character(model_type %||% "mlp"))
  args <- dl_args %||% list()
  args$response_family <- NULL
  out <- switch(
    model_key,
    "cnn" = canon_cnn_args(args),
    "mlp" = canon_mlp_args(args, with_attention = FALSE),
    "mlp_with_attention" = canon_mlp_args(args, with_attention = TRUE),
    "ft_transformer" = canon_ft_args(args),
    "resnet" = canon_resnet_args(args),
    "saint" = canon_saint_args(args),
    "tabnet" = canon_tabnet_args(args),
    "node" = canon_node_args(args),
    "deepfm" = canon_deepfm_args(args),
    "dcnv2" = canon_dcnv2_args(args),
    "nam" = canon_nam_args(args),
    "moe" = canon_moe_args(args),
    "gp_dkl" = canon_gp_dkl_args(args),
    args
  )
  # NULL or a single NA means "not set": leave it to the Python default
  Filter(function(v) !is.null(v) && !(length(v) == 1L && is.atomic(v) && is.na(v)), out)
}

gp_dl_payload_matrix <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  if (!is.null(dimnames(x))) {
    dimnames(x) <- NULL
  }
  x
}

gp_dl_bridge_script_path <- function() {
  pyfile <- system.file("python/dl_bridge.py", package = "PredictProR")
  if (!nzchar(pyfile)) {
    candidates <- c("inst/python/dl_bridge.py", "python/dl_bridge.py", "../inst/python/dl_bridge.py")
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) {
      pyfile <- normalizePath(hit[1], winslash = "/", mustWork = TRUE)
    }
  }
  if (!nzchar(pyfile)) {
    stop("dl_bridge.py not found under inst/python/.", call. = FALSE)
  }
  pyfile
}

gp_detect_dl_python <- function() {
  py <- Sys.getenv("PREDICTPRO_DL_PYTHON", unset = "")
  if (nzchar(py) && file.exists(py)) {
    return(normalizePath(py, winslash = "/", mustWork = TRUE))
  }
  py <- Sys.getenv("PREDICTPRO_PYTHON_DL", unset = "")
  if (nzchar(py) && file.exists(py)) {
    return(normalizePath(py, winslash = "/", mustWork = TRUE))
  }
  preferred <- gp_preferred_python(purpose = "dl")
  if (!is.null(preferred) && nzchar(preferred)) {
    return(normalizePath(preferred, winslash = "/", mustWork = TRUE))
  }
  NULL
}

gp_dl_subprocess_env <- function() {
  paste0(
    "CUBLAS_WORKSPACE_CONFIG=",
    gp_configure_torch_runtime_env()
  )
}

gp_dl_run_cli <- function(command,
                          args,
                          python_bin = NULL,
                          max_attempts = 1L,
                          run_fn = system2) {
  python_bin <- python_bin %||% gp_detect_dl_python()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop(
      "No configured Python runtime was found for Python DL models.\n",
      "Set PREDICTPRO_DL_PYTHON or configure a preferred DL Python runtime.",
      call. = FALSE
    )
  }
  max_attempts <- suppressWarnings(as.integer(max_attempts)[1L])
  if (!is.finite(max_attempts) || max_attempts < 1L || max_attempts > 2L) {
    stop("max_attempts must be one or two.", call. = FALSE)
  }
  if (!is.function(run_fn)) {
    stop("run_fn must be a function.", call. = FALSE)
  }
  cmd_args <- c(gp_dl_bridge_script_path(), command, args)
  gp_dl_subprocess_env()
  for (attempt in seq_len(max_attempts)) {
    out <- suppressWarnings(run_fn(
      python_bin,
      args = gp_quote_system_args(cmd_args),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- suppressWarnings(as.integer(attr(out, "status") %||% 0L)[1L])
    if (identical(status, 0L)) return(invisible(out))
    if (attempt < max_attempts) {
      # Training bridges write only to a call-private temporary output
      # directory. Remove known partial products before repeating the exact
      # deterministic command; retain params.json, which is an input.
      out_flag <- match("--out-dir", args)
      if (!is.na(out_flag) && out_flag < length(args)) {
        out_dir <- args[[out_flag + 1L]]
        partial_products <- file.path(
          out_dir,
          c("predictions.csv", "bootstrap.csv", "standard_error.csv",
            "probabilities.csv", "meta.json", "embeddings.csv")
        )
        unlink(partial_products[file.exists(partial_products)], force = TRUE)
      }
      warning(
        "DL Python bridge exited with status ", status,
        "; retrying the same deterministic command once.",
        call. = FALSE
      )
      next
    }
    stop(
      paste0(
        "DL Python bridge failed with status ", status, ".\n",
        paste(out, collapse = "\n")
      ),
      call. = FALSE
    )
  }
}

gp_dl_write_matrix_csv <- function(x, path) {
  gp_bridge_csv_write(as.data.frame(unname(gp_dl_payload_matrix(x)), check.names = FALSE), path, quote = FALSE)
  invisible(path)
}

gp_dl_write_vector_csv <- function(x, path) {
  gp_bridge_csv_write(data.frame(value = x, stringsAsFactors = FALSE), path, quote = TRUE)
  invisible(path)
}

gp_dl_write_target_csv <- function(x, path) {
  x_mat <- if (is.matrix(x) || is.data.frame(x)) as.matrix(x) else matrix(x, ncol = 1L)
  gp_bridge_csv_write(as.data.frame(unname(x_mat), check.names = FALSE), path, quote = TRUE)
  invisible(path)
}

gp_dl_write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = TRUE), path, useBytes = TRUE)
  invisible(path)
}

gp_dl_write_seed_json <- function(x, path) {
  seeds <- gp_dl_validate_seed_vector(x, name = "training_seeds")
  gp_dl_write_json(unname(as.list(seeds)), path)
}

gp_dl_bridge_read_cli_result <- function(out_dir) {
  pred_path <- file.path(out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("DL Python bridge did not write predictions.csv.", call. = FALSE)
  }
  pred_df <- utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE)
  prob_path <- file.path(out_dir, "probabilities.csv")
  prob <- if (file.exists(prob_path)) {
    as.matrix(utils::read.csv(prob_path, stringsAsFactors = FALSE, check.names = FALSE))
  } else {
    NULL
  }
  meta_path <- file.path(out_dir, "meta.json")
  meta <- if (file.exists(meta_path)) {
    tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = TRUE), error = function(e) list())
  } else {
    list()
  }
  list(
    predictions = if (ncol(pred_df) == 1L) pred_df[[1L]] else as.matrix(pred_df),
    probabilities = prob,
    classes = meta$classes %||% NULL,
    meta = meta
  )
}

gp_dl_bridge_response_family <- function(response_family = NULL, y = NULL) {
  fam <- tolower(trimws(as.character(response_family %||% "auto")[1L]))
  if (fam %in% c("multitask", "multi_task", "multioutput", "multi_output", "multitask_regression")) {
    return("multitask_regression")
  }
  gp_resolve_response_family(response_family, y = y)
}

gp_dl_bridge_class_levels <- function(y, fam, class_levels = NULL) {
  if (!identical(fam, "binary") && !identical(fam, "multiclass")) {
    return(NULL)
  }
  if (!is.null(class_levels) && length(class_levels)) {
    return(as.character(class_levels))
  }
  if (exists("gp_py_ml_class_levels", mode = "function")) {
    return(gp_py_ml_class_levels(y, fam))
  }
  if (is.factor(y) || is.ordered(y)) {
    return(as.character(levels(y)))
  }
  sort(unique(as.character(stats::na.omit(y))))
}

gp_dl_bridge_prepare_target <- function(y, fam, class_levels = NULL) {
  if (!identical(fam, "binary") && !identical(fam, "multiclass")) {
    return(list(y = y, class_levels = NULL))
  }
  levs <- gp_dl_bridge_class_levels(y, fam, class_levels)
  if (!length(levs)) {
    stop("Classification deep-learning target has no observed class levels.", call. = FALSE)
  }
  y_chr <- as.character(y)
  y_int <- suppressWarnings(as.integer(y_chr))
  non_na <- !is.na(y)
  encoded_levels <- seq_along(levs) - 1L
  already_encoded <- any(non_na) &&
    all(!is.na(y_int[non_na])) &&
    all(y_int[non_na] %in% encoded_levels) &&
    !any(y_chr[non_na] %in% levs)
  if (already_encoded) {
    return(list(y = as.integer(y_int), class_levels = levs))
  }
  y_encoded <- as.integer(factor(y_chr, levels = levs)) - 1L
  bad <- non_na & is.na(y_encoded)
  if (any(bad)) {
    stop("Classification deep-learning target contains labels outside class_levels.", call. = FALSE)
  }
  list(y = y_encoded, class_levels = levs)
}

gp_dl_bridge_fit_predict_gaussian_fast <- function(model_type,
                                                   X_train,
                                                   y_train,
                                                   X_test,
                                                   dl_args = list()) {
  res <- gp_dl_bridge_fit_predict_raw(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    response_family = "gaussian",
    dl_args = dl_args
  )
  as.numeric(res$predictions)
}

gp_dl_bridge_fit_predict_raw <- function(model_type,
                                         X_train,
                                         y_train,
                                         X_test,
                                         response_family = "gaussian",
                                         dl_args = list(),
                                         class_levels = NULL) {
  fam <- gp_dl_bridge_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("DL fit-predict", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  fit_params <- gp_dl_bridge_fit_params(model_type = model_type, dl_args = dl_args)
  x_train_payload <- unname(gp_dl_payload_matrix(X_train))
  x_test_payload <- unname(gp_dl_payload_matrix(X_test))
  target <- gp_dl_bridge_prepare_target(y_train, fam, class_levels = class_levels)
  y_payload <- target$y
  class_levels <- target$class_levels
  write_payload <- function(input_dir) {
    gp_dl_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
    gp_dl_write_target_csv(y_payload, file.path(input_dir, "y_train.csv"))
    gp_dl_write_matrix_csv(x_test_payload, file.path(input_dir, "x_test.csv"))
    if (!is.null(class_levels)) {
      gp_dl_write_json(as.character(class_levels), file.path(input_dir, "class_levels.json"))
    }
    invisible(input_dir)
  }
  payload_cache <- gp_bridge_payload_cache_prepare(
    prefix = "DL",
    purpose = paste("fit_predict", fam, sep = "::"),
    key_parts = list(
      response_family = fam,
      x_train = x_train_payload,
      y_train = y_payload,
      x_test = x_test_payload,
      class_levels = as.character(class_levels %||% character())
    ),
    cells = gp_bridge_payload_cache_cells(x_train_payload, y_payload, x_test_payload),
    build_fun = write_payload,
    min_cells = 20000L,
    max_entries = 12L
  )
  if (is.null(payload_cache)) {
    bridge_dir <- tempfile("predictpror_dl_bridge_")
    input_dir <- file.path(bridge_dir, "input")
    out_dir <- file.path(bridge_dir, "out")
    dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(bridge_dir, recursive = TRUE, force = TRUE), add = TRUE)
    write_payload(input_dir)
  } else {
    input_dir <- payload_cache$input_dir
    out_dir <- tempfile("predictpror_dl_bridge_out_")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  x_train_csv <- file.path(input_dir, "x_train.csv")
  y_train_csv <- file.path(input_dir, "y_train.csv")
  x_test_csv <- file.path(input_dir, "x_test.csv")
  levels_json <- file.path(input_dir, "class_levels.json")
  params_json <- file.path(out_dir, "params.json")
  gp_dl_write_json(Filter(Negate(is.null), fit_params), params_json)
  args <- c(
    "--model-type", tolower(as.character(model_type)),
    "--task", fam,
    "--x-train-csv", x_train_csv,
    "--y-train-csv", y_train_csv,
    "--x-test-csv", x_test_csv,
    "--params-json", params_json,
    "--out-dir", out_dir
  )
  if (!is.null(class_levels)) {
    args <- c(args, "--class-levels-json", levels_json)
  }
  # One fit-predict runs per CV fold; retry a transient native exit once, as
  # the bootstrap and multi-trait DL paths already do.
  gp_dl_run_cli("fit-predict", args, max_attempts = 2L)
  invisible(gc(verbose = FALSE))

  gp_dl_bridge_read_cli_result(out_dir)
}

gp_dl_bridge_fit_predict_batch_raw <- function(jobs) {
  if (is.null(jobs) || !length(jobs)) {
    return(list())
  }
  batch_dir <- tempfile("predictpror_dl_batch_")
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(batch_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_jobs <- vector("list", length(jobs))
  out_dirs <- character(length(jobs))

  for (i in seq_along(jobs)) {
    job <- jobs[[i]]
    model_type <- tolower(as.character(job$model_type %||% "mlp")[1L])
    fam <- gp_dl_bridge_response_family(job$response_family %||% "gaussian", y = job$y_train)
    filtered <- gp_py_ml_filter_missing_training(
      X_train = job$X_train,
      y_train = job$y_train,
      response_family = fam,
      context = paste("DL batch fit-predict", model_type)
    )
    job$X_train <- filtered$X_train
    job$y_train <- filtered$y_train
    fit_params <- gp_dl_bridge_fit_params(model_type = model_type, dl_args = job$dl_args %||% list())
    x_train_payload <- unname(gp_dl_payload_matrix(job$X_train))
    x_test_payload <- unname(gp_dl_payload_matrix(job$X_test %||% job$X_pred))
    target <- gp_dl_bridge_prepare_target(job$y_train, fam, class_levels = job$class_levels %||% NULL)
    y_payload <- target$y
    class_levels <- target$class_levels
    write_payload <- function(input_dir) {
      gp_dl_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
      gp_dl_write_target_csv(y_payload, file.path(input_dir, "y_train.csv"))
      gp_dl_write_matrix_csv(x_test_payload, file.path(input_dir, "x_test.csv"))
      if (!is.null(class_levels)) {
        gp_dl_write_json(as.character(class_levels), file.path(input_dir, "class_levels.json"))
      }
      invisible(input_dir)
    }
    payload_cache <- gp_bridge_payload_cache_prepare(
      prefix = "DL",
      purpose = paste("fit_predict", fam, sep = "::"),
      key_parts = list(
        response_family = fam,
        x_train = x_train_payload,
        y_train = y_payload,
        x_test = x_test_payload,
        class_levels = as.character(class_levels %||% character())
      ),
      cells = gp_bridge_payload_cache_cells(x_train_payload, y_payload, x_test_payload),
      build_fun = write_payload,
      min_cells = 20000L,
      max_entries = 12L
    )
    if (is.null(payload_cache)) {
      input_dir <- file.path(batch_dir, "input", sprintf("job_%05d", i))
      dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
      write_payload(input_dir)
    } else {
      input_dir <- payload_cache$input_dir
    }
    out_dir <- file.path(batch_dir, "out", sprintf("job_%05d", i))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    params_json <- file.path(out_dir, "params.json")
    gp_dl_write_json(Filter(Negate(is.null), fit_params), params_json)

    out_dirs[[i]] <- out_dir
    manifest_jobs[[i]] <- list(
      id = as.character(job$id %||% i),
      model_type = model_type,
      task = fam,
      x_train_csv = file.path(input_dir, "x_train.csv"),
      y_train_csv = file.path(input_dir, "y_train.csv"),
      x_test_csv = file.path(input_dir, "x_test.csv"),
      params_json = params_json,
      class_levels_json = if (!is.null(class_levels)) file.path(input_dir, "class_levels.json") else NULL,
      out_dir = out_dir
    )
  }

  manifest_json <- file.path(batch_dir, "manifest.json")
  gp_dl_write_json(list(jobs = manifest_jobs), manifest_json)
  gp_dl_run_cli("batch-fit-predict", c("--manifest-json", manifest_json))
  lapply(out_dirs, gp_dl_bridge_read_cli_result)
}

gp_dl_bridge_format_prediction <- function(res,
                                           fam,
                                           class_levels = NULL) {
  if (identical(fam, "gaussian")) {
    return(as.numeric(res$predictions))
  }

  if (identical(fam, "multitask_regression")) {
    pred <- as.matrix(res$predictions)
    storage.mode(pred) <- "double"
    return(pred)
  }

  prob <- if (is.null(res$probabilities)) NULL else {
    prob_mat <- as.matrix(res$probabilities)
    storage.mode(prob_mat) <- "double"
    prob_mat
  }

  if (identical(fam, "binary")) {
    pred <- if (!is.null(prob)) as.numeric(prob[, ncol(prob), drop = TRUE]) else as.numeric(res$predictions)
    pred <- pmax(0, pmin(1, pred))
    if (!is.null(prob)) {
      if (!is.null(class_levels) && length(class_levels) == ncol(prob)) {
        colnames(prob) <- class_levels
      } else if (is.null(colnames(prob))) {
        colnames(prob) <- c("0", "1")[seq_len(ncol(prob))]
      }
      attr(pred, "probabilities") <- prob
    }
    return(pred)
  }

  if (!is.null(prob)) {
    classes <- as.character(class_levels %||% colnames(prob) %||% res$classes %||% paste0("class_", seq_len(ncol(prob))))
    if (length(classes) == ncol(prob)) {
      colnames(prob) <- classes
      pred <- classes[max.col(prob, ties.method = "first")]
      attr(pred, "probabilities") <- prob
      return(pred)
    }
  }

  pred <- as.character(res$predictions)
  if (!is.null(prob)) {
    attr(pred, "probabilities") <- prob
  }
  pred
}

gp_dl_bridge_fit_predict <- function(model_type,
                                     X_train,
                                     y_train,
                                     X_test,
                                     response_family = "gaussian",
                                     dl_args = list(),
                                     class_levels = NULL) {
  fam <- gp_dl_bridge_response_family(response_family, y = y_train)
  class_levels <- gp_dl_bridge_class_levels(y_train, fam, class_levels = class_levels)
  res <- gp_dl_bridge_fit_predict_raw(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    response_family = fam,
    dl_args = dl_args,
    class_levels = class_levels
  )

  gp_dl_bridge_format_prediction(res, fam = fam, class_levels = class_levels)
}

gp_dl_bridge_fit_predict_batch <- function(jobs) {
  if (is.null(jobs) || !length(jobs)) {
    return(list())
  }
  prepared <- lapply(jobs, function(job) {
    fam <- gp_dl_bridge_response_family(job$response_family %||% "gaussian", y = job$y_train)
    job$response_family <- fam
    job$class_levels <- gp_dl_bridge_class_levels(job$y_train, fam, class_levels = job$class_levels %||% NULL)
    job
  })
  raw <- gp_dl_bridge_fit_predict_batch_raw(prepared)
  Map(function(result, job) {
    gp_dl_bridge_format_prediction(
      result,
      fam = job$response_family,
      class_levels = job$class_levels %||% NULL
    )
  }, raw, prepared)
}

gp_dl_bridge_bootstrap <- function(model_type,
                                   X_train,
                                   y_train,
                                   X_pred,
                                   n_bootstrap,
                                   response_family = "gaussian",
                                   dl_args = list(),
                                   seed = 123L,
                                   training_seeds = NULL,
                                   seed_aggregation = "mean") {
  fam <- gp_dl_bridge_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("DL bootstrap", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  fit_params <- gp_dl_bridge_fit_params(model_type = model_type, dl_args = dl_args)
  x_train_payload <- unname(gp_dl_payload_matrix(X_train))
  x_pred_payload <- unname(gp_dl_payload_matrix(X_pred))
  y_payload <- if (identical(fam, "gaussian")) y_train else y_train
  write_payload <- function(input_dir) {
    gp_dl_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
    gp_dl_write_target_csv(y_payload, file.path(input_dir, "y_train.csv"))
    gp_dl_write_matrix_csv(x_pred_payload, file.path(input_dir, "x_pred.csv"))
    invisible(input_dir)
  }
  payload_cache <- gp_bridge_payload_cache_prepare(
    prefix = "DL",
    purpose = paste("bootstrap", fam, sep = "::"),
    key_parts = list(
      response_family = fam,
      x_train = x_train_payload,
      y_train = y_payload,
      x_pred = x_pred_payload
    ),
    cells = gp_bridge_payload_cache_cells(x_train_payload, y_payload, x_pred_payload),
    build_fun = write_payload,
    min_cells = 20000L,
    max_entries = 12L
  )
  if (is.null(payload_cache)) {
    bridge_dir <- tempfile("predictpror_dl_bootstrap_")
    input_dir <- file.path(bridge_dir, "input")
    out_dir <- file.path(bridge_dir, "out")
    dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(bridge_dir, recursive = TRUE, force = TRUE), add = TRUE)
    write_payload(input_dir)
  } else {
    input_dir <- payload_cache$input_dir
    out_dir <- tempfile("predictpror_dl_bootstrap_out_")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  x_train_csv <- file.path(input_dir, "x_train.csv")
  y_train_csv <- file.path(input_dir, "y_train.csv")
  x_pred_csv <- file.path(input_dir, "x_pred.csv")
  params_json <- file.path(out_dir, "params.json")
  training_seeds_json <- file.path(out_dir, "training_seeds.json")
  gp_dl_write_json(Filter(Negate(is.null), fit_params), params_json)
  training_seeds <- gp_dl_validate_seed_vector(
    training_seeds %||% fit_params$random_seed %||% seed %||% 123L,
    name = "training_seeds"
  )
  gp_dl_write_seed_json(training_seeds, training_seeds_json)
  gp_dl_run_cli(
    "bootstrap-fit-predict",
    c(
      "--model-type", tolower(as.character(model_type)),
      "--task", fam,
      "--x-train-csv", x_train_csv,
      "--y-train-csv", y_train_csv,
      "--x-pred-csv", x_pred_csv,
      "--params-json", params_json,
      "--n-bootstrap", as.character(as.integer(n_bootstrap)),
      "--seed", as.character(as.integer(seed %||% 123L)),
      "--training-seeds-json", training_seeds_json,
      "--seed-aggregation", gp_dl_seed_aggregation(seed_aggregation),
      "--out-dir", out_dir
    ),
    max_attempts = 2L
  )
  boot_path <- file.path(out_dir, "bootstrap.csv")
  if (!file.exists(boot_path)) {
    stop("DL Python bridge did not write bootstrap.csv.", call. = FALSE)
  }
  boot <- as.matrix(utils::read.csv(boot_path, stringsAsFactors = FALSE, check.names = FALSE))
  seed_predictions_path <- file.path(out_dir, "seed_predictions.csv")
  seed_predictions <- if (file.exists(seed_predictions_path)) {
    as.matrix(utils::read.csv(
      seed_predictions_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  } else {
    NULL
  }
  seed_variability_path <- file.path(out_dir, "seed_variability.csv")
  seed_variability <- if (file.exists(seed_variability_path)) {
    utils::read.csv(
      seed_variability_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }
  meta_path <- file.path(out_dir, "meta.json")
  meta <- if (file.exists(meta_path)) {
    tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = TRUE), error = function(e) list())
  } else {
    list()
  }
  invisible(gc(verbose = FALSE))

  list(
    bootstrap = as.matrix(boot),
    classes = meta$classes %||% NULL,
    task = fam,
    seed_manifest = data.frame(
      seed_index = seq_along(training_seeds),
      training_seed = as.integer(training_seeds),
      source = rep("resolved_by_R_api", length(training_seeds)),
      stringsAsFactors = FALSE
    ),
    seed_predictions = seed_predictions,
    seed_variability = seed_variability,
    seed_aggregation = meta$seed_aggregation %||% seed_aggregation,
    n_model_fits = meta$n_model_fits %||%
      (as.double(n_bootstrap) * length(training_seeds))
  )
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
#   # ---- y: make sure it's numeric-coded even if character ----
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
#   msg <- ""
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
