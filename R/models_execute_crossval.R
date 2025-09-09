
resolve_model_name <- function(x, name_lookup) {
  if (is.null(x) || length(x) == 0) return(NULL)
  res <- name_lookup[x]
  res <- res[!is.na(res)]
  if (length(res) == 0) return(NULL)
  unname(res)
}

replace_with_canonical <- function(x, canonical_names, friendly_names) {
  if (is.null(x) || length(x) == 0) return(NULL)

  # Normalize to lowercase
  k <- tolower(trimws(x))

  # Create lowercase lookup
  keys <- tolower(c(canonical_names, friendly_names))
  vals <- rep(canonical_names, 2)
  lookup <- setNames(vals, keys)

  # Lookup canonical names
  res <- lookup[k]

  # Preserve unknowns as original
  res[is.na(res)] <- x[is.na(res)]

  unname(res)
}

replace_with_friendly_name <- function(x, friendly_name_lookup) {
  if (is.null(x) || length(x) == 0) return(NULL)

  k <- tolower(trimws(x))
  names(friendly_name_lookup) <- tolower(names(friendly_name_lookup))

  res <- friendly_name_lookup[k]
  res[is.na(res)] <- x[is.na(res)]

  unname(res)
}
# =========================
# Smart parallel helpers
# =========================
`%||%` <- function(a, b) if (is.null(a)) b else a

# Detect the Python binary reticulate is using (or from RETICULATE_PYTHON)
gp_detect_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  conf <- try(reticulate::py_config(), silent = TRUE)
  if (!inherits(conf, "try-error") && !is.null(conf$python) && nzchar(conf$python)) {
    return(conf$python)
  }
  py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(py)) return(py)
  NULL
}

# Initialize that Python once in the current R session
gp_init_python_once <- function(py_bin) {
  if (!requireNamespace("reticulate", quietly = TRUE) || is.null(py_bin) || !nzchar(py_bin)) return(invisible(FALSE))
  if (!reticulate::py_available(initialize = FALSE)) {
    reticulate::use_python(py_bin, required = TRUE)
    invisible(TRUE)
  } else {
    # Already initialized; if different, we can't switch now, just warn.
    conf <- reticulate::py_config()
    same <- tryCatch(
      identical(normalizePath(conf$python, winslash = "/"), normalizePath(py_bin, winslash = "/")),
      error = function(e) FALSE
    )
    if (!same) message(sprintf("[reticulate] Already initialized to '%s'; requested '%s'.", conf$python, py_bin))
    invisible(FALSE)
  }
}

sp_available_cores <- function(num_cores = NULL) {
  if (!is.null(num_cores) && is.finite(num_cores) && num_cores >= 1) return(as.integer(num_cores))
  cores <- tryCatch({
    if (requireNamespace("parallelly", quietly = TRUE)) parallelly::availableCores()
    else parallel::detectCores(logical = TRUE)
  }, error = function(e) parallel::detectCores(logical = TRUE))
  max(1L, as.integer(cores))
}

sp_estimate_bytes <- function(x) {
  safe_sz <- function(obj) {
    sz <- tryCatch(utils::object.size(obj), error = function(e) NA)
    as.numeric(sz)
  }
  if (is.list(x)) sum(vapply(x, safe_sz, numeric(1), USE.NAMES = FALSE), na.rm = TRUE) else safe_sz(x)
}

# Optional user hints to mark models that do their own threading:
# options(gp.parallel.hints = list(xgboost="internal", ranger="internal"))
sp_model_has_internal_parallel <- function(model, params = list()) {
  hints <- getOption("gp.parallel.hints", NULL)
  if (is.list(hints)) {
    key <- tolower(as.character(model))
    if (key %in% names(hints)) {
      v <- tolower(as.character(hints[[key]]))
      if (v %in% c("internal","yes","true","1")) return(TRUE)
      if (v %in% c("external","no","false","0")) return(FALSE)
    }
  }
  thread_keys <- c("nthread","nthreads","n_threads","num_threads","num.threads",
                   "n_jobs","threads","omp_threads","openmp_threads","num_workers")
  vals <- unlist(params, recursive = TRUE, use.names = TRUE)
  if (length(vals)) {
    nk <- names(vals); hit <- which(tolower(nk) %in% thread_keys)
    if (length(hit)) {
      to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))
      tv <- to_num(vals[hit])
      if (any(is.finite(tv) & tv > 1)) return(TRUE)
    }
  }
  FALSE
}

# sp_decide_policy <- function(models,
#                               model_params_list = NULL,
#                               globals = list(),
#                               n_tasks = NULL,
#                               user_mode = c("auto","future","sequential","base_parallel","foreach"),
#                               num_cores = NULL,
#                               globals_max_GB = 4,
#                               verbose = FALSE,
#                               sys_name = Sys.info()[["sysname"]],
#                               prefer_fork = TRUE) {
#
#   user_mode <- match.arg(user_mode)
#   workers   <- sp_available_cores(num_cores)
#
#   gbytes <- sp_estimate_bytes(globals)
#   g_gb   <- gbytes / (1024^3)
#
#   params_for <- function(idx) {
#     if (is.null(model_params_list)) return(list())
#     if (is.list(model_params_list) && length(model_params_list) == length(models)) return(model_params_list[[idx]])
#     if (is.list(model_params_list) && !is.null(model_params_list[[as.character(models[idx])]]))
#       return(model_params_list[[as.character(models[idx])]])
#     list()
#   }
#   internal_flags <- vapply(seq_along(models), function(i) {
#     sp_model_has_internal_parallel(models[[i]], params_for(i))
#   }, logical(1))
#
#   any_internal <- any(internal_flags)
#   few_tasks    <- is.null(n_tasks) || !is.finite(n_tasks) || n_tasks < 2
#
#   backend <- "future"
#   plan    <- if (sys_name == "Windows" || !prefer_fork) "multisession" else "multicore"
#   reason  <- NULL
#
#   if (user_mode == "sequential") {
#     backend <- "sequential"; plan <- "sequential"; reason <- "user_forced_sequential"
#   } else if (user_mode == "base_parallel") {
#     backend <- "base_parallel"; plan <- "base_parallel"; reason <- "user_forced_base_parallel"
#   } else if (user_mode == "foreach") {
#     backend <- "foreach"; plan <- "foreach"; reason <- "user_forced_foreach"
#   } else if (user_mode == "future") {
#     if (g_gb > globals_max_GB * 0.8) { backend <- "sequential"; plan <- "sequential"; reason <- sprintf("globals %.2fGB exceed limit %.2fGB", g_gb, globals_max_GB) }
#     else if (few_tasks || workers < 2) { backend <- "sequential"; plan <- "sequential"; reason <- "too_few_tasks_or_workers" }
#     else if (any_internal) { backend <- "sequential"; plan <- "sequential"; reason <- "internal_parallel_detected" }
#   } else { # auto
#     if (g_gb > globals_max_GB * 0.8) {
#       backend <- "sequential"; plan <- "sequential"; reason <- sprintf("globals %.2fGB exceed limit %.2fGB", g_gb, globals_max_GB)
#     } else if (few_tasks || workers < 2) {
#       backend <- "sequential"; plan <- "sequential"; reason <- "too_few_tasks_or_workers"
#     } else if (any_internal) {
#       backend <- "sequential"; plan <- "sequential"; reason <- "internal_parallel_detected"
#     } else {
#       backend <- "future"; plan <- if (sys_name == "Windows" || !prefer_fork) "multisession" else "multicore"; reason <- "future_is_beneficial"
#     }
#   }
#
#   if (!identical(backend, "future") || is.null(n_tasks) || !is.finite(n_tasks) || workers < 2) {
#     chunk_size <- NULL
#   } else {
#     chunk_size <- ceiling(n_tasks / workers)
#     chunk_size <- max(1L, min(chunk_size, max(1L, floor(n_tasks / (workers * 1.5)))))
#   }
#
#   if (isTRUE(verbose)) {
#     message(sprintf(
#       "[smart-parallel] backend=%s, plan=%s, workers=%d, globals=%.2f GB, tasks=%s, any_internal=%s, reason=%s",
#       backend, plan, workers, g_gb, if (is.null(n_tasks)) "NA" else n_tasks, any_internal, reason %||% "n/a"
#     ))
#   }
#
#   list(
#     backend    = backend,
#     plan       = plan,
#     workers    = sp_available_cores(num_cores),
#     chunk_size = chunk_size,
#     reason     = reason,
#     internal_flags = internal_flags,
#     globals_gb = g_gb
#   )
# }

# Compose multiple initializers for PSOCK workers (optional helper)
compose_initializers <- function(...) {
  inits <- list(...)
  function() {
    for (f in inits) if (is.function(f)) try(f(), silent = TRUE)
  }
}

# Optional: force single-threaded BLAS/OMP inside workers (use only if you decide to parallelize CPU-heavy tasks)
init_single_thread_blas <- function() {
  # Environment variables many libs respect
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1L),  silent = TRUE)
  }
}

detect_hidden_threads <- function() {
  vals <- c(
    suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", NA))),
    suppressWarnings(as.integer(Sys.getenv("OPENBLAS_NUM_THREADS", NA))),
    suppressWarnings(as.integer(Sys.getenv("MKL_NUM_THREADS", NA))),
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      c(tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_),
        tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NA_integer_))
    } else integer(0)
  )
  max_val <- suppressWarnings(max(vals, na.rm = TRUE))
  if (!is.finite(max_val)) NA_integer_ else max_val
}
# Decide smart parallel/sequential policy per task list
sp_decide_policy <- function(models,
                             model_params_list,
                             globals,
                             n_tasks,
                             user_mode = c("auto","future","sequential","base_parallel","foreach"),
                             num_cores = NULL,
                             globals_max_GB = 4,
                             # NEW: RKHS controls
                             rkhs_mem_threshold_GB = getOption("gp.rkhs.mem.threshold.GB", 1.0),
                             always_seq_rkhs      = getOption("gp.always.seq.RKHS", TRUE),
                             verbose = TRUE,
                             sys_name = Sys.info()[["sysname"]],
                             prefer_fork = TRUE) {

  user_mode <- match.arg(user_mode)
  n_tasks   <- as.integer(n_tasks %||% length(models))

  # ---- helpers -------------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a)) b else a
  size_bytes <- function(x) if (is.null(x)) 0 else as.numeric(utils::object.size(x))

  # ---- early exit: user forced sequential ---------------------------------
  if (identical(user_mode, "sequential")) {
    return(list(
      backend = "sequential",
      plan = NULL,
      workers = 1L,
      chunk_size = NULL,
      internal_flags = rep(TRUE, length(models))
    ))
  }

  # ---- heuristics: models that should default to sequential ----------------
  deep_models <- c("cnn","resnet","ft_transformer","saint","tabnet","node",
                   "deepfm","dcnv2","nam","moe","gp_dkl","mlp_with_attention","mlp")

  # 1) Internal-parallel / Python DL -> sequential
  seq_deep <- models %in% deep_models

  # 2) ASReml present + GBLUP -> sequential (threaded lib; avoid oversubscription)
  has_asreml_prep <- !is.null(globals$asreml)
  seq_asreml <- (models == "GBLUP") & has_asreml_prep

  # 3) XGBoost with nthread>1 -> sequential
  seq_xgb <- vapply(seq_along(models), function(i) {
    if (!identical(models[i], "Xgboost")) return(FALSE)
    p  <- model_params_list[[i]]
    nt <- p$xgb_nthread %||% p$nthread %||% p$threads %||% NA_integer_
    isTRUE(!is.na(nt) && nt > 1L)
  }, logical(1))

  hidden_thr <- detect_hidden_threads()

  #seq_RKHS <- vapply(seq_along(models), function(i) identical(models[i], "RKHS"), logical(1))

  # Memory estimate for RKHS kernel: ~8 * n^2 bytes
  n <- tryCatch(NROW(globals$pheno_data), error = function(e) NA_integer_)
  rkhs_kernel_GB <- if (is.finite(n)) (8 * as.numeric(n) * as.numeric(n)) / (1024^3) else NA_real_

  # Policies:
  #  A) Default: always sequential for RKHS (safe)
  #  B) If user disables always_seq_rkhs, then sequential if hidden threads >1 OR kernel too big
  if (isTRUE(always_seq_rkhs)) {
    seq_RKHS <- (models == "RKHS")
  } else {
    seq_rkhs_hidden <- (models == "RKHS") & is.finite(hidden_thr) & hidden_thr > 1L
    seq_rkhs_mem    <- (models == "RKHS") & is.finite(rkhs_kernel_GB) & rkhs_kernel_GB > rkhs_mem_threshold_GB
    seq_RKHS <- seq_rkhs_hidden | seq_rkhs_mem
  }
  internal_flags <- seq_deep | seq_asreml | seq_xgb | seq_RKHS


  # ---- globals payload guard ----------------------------------------------
  tot_GB <- sum(c(
    size_bytes(globals$pheno_data),
    size_bytes(globals$omics_data),
    size_bytes(globals$bayes),
    size_bytes(globals$asreml)
  )) / (1024^3)

  too_big_to_ship <- is.finite(tot_GB) && (tot_GB > globals_max_GB)
  if (isTRUE(too_big_to_ship) && isTRUE(verbose)) {
    message(sprintf("[policy] Globals ~ %s GB > %s GB → force sequential.", safe_num(tot_GB), safe_num(globals_max_GB)))
  }
  if (isTRUE(too_big_to_ship)) {
    internal_flags[] <- TRUE
  }

  # ---- very small parallel workload? just go sequential --------------------
  n_par_tasks <- sum(!internal_flags)
  if (n_par_tasks < 2L) {
    internal_flags[] <- TRUE
    n_par_tasks <- 0L
  }

  # ---- choose backend ------------------------------------------------------
  is_windows <- identical(tolower(sys_name), "windows")
  n_cores <- if (is.null(num_cores)) {
    max(1L, parallel::detectCores(logical = TRUE) - 1L)
  } else {
    max(1L, as.integer(num_cores))
  }
  workers <- if (n_par_tasks > 0) max(1L, min(n_cores, n_par_tasks)) else 1L

  if (identical(user_mode, "base_parallel")) {
    backend <- "base_parallel"; plan <- NULL
  } else if (identical(user_mode, "foreach")) {
    backend <- "foreach"; plan <- NULL
  } else if (identical(user_mode, "future") || identical(user_mode, "auto")) {
    backend <- "future"
    # Use strategy function (not a string)
    if (!is_windows && isTRUE(prefer_fork)) {
      plan <- future::multicore
    } else {
      plan <- future::multisession
    }
  } else {
    backend <- "future"; plan <- future::multisession
  }

  chunk_size <- if (n_par_tasks > 0 && workers > 0) ceiling(n_par_tasks / workers) else NULL

  if (isTRUE(verbose)) {
    plan_name <- if (is.null(plan)) "NA" else if (identical(plan, future::multicore)) "multicore" else "multisession"

    # optional tails for extra fields (only shown when finite)
    tail_hidden <- if (exists("hidden_thr") && is.numeric(hidden_thr) && is.finite(hidden_thr))
      sprintf(" | hidden_thr=%s GB", safe_num(hidden_thr)) else ""

    tail_rkhs <- if (exists("rkhs_kernel_GB") && is.numeric(rkhs_kernel_GB) && is.finite(rkhs_kernel_GB))
      sprintf(" | rkhs=%s GB", safe_num(rkhs_kernel_GB)) else ""

    message(sprintf(
      "[policy] backend=%s plan=%s workers=%d chunk=%s | seq=%d par=%d | globals=%s GB%s%s",
      backend,
      plan_name,
      workers,
      if (is.null(chunk_size)) "NA" else as.character(chunk_size),
      sum(internal_flags), sum(!internal_flags),
      safe_num(tot_GB),
      tail_hidden,
      tail_rkhs
    ))

    if (any(seq_deep))
      message(sprintf("         deep/tabular forced seq: %s", paste(which(seq_deep), collapse = ",")))
    if (any(seq_asreml))
      message(sprintf("         GBLUP(asreml) forced seq: %s", paste(which(seq_asreml), collapse = ",")))
    if (any(seq_xgb))
      message(sprintf("         Xgboost(nthread>1) seq: %s", paste(which(seq_xgb), collapse = ",")))
    if (any(models == "RKHS"))
      message(sprintf("         RKHS policy: always_seq_rkhs=%s, threshold=%s GB",
                      as.character(always_seq_rkhs), safe_num(rkhs_mem_threshold_GB)))
  }

  list(
    backend = backend,
    plan = plan,
    workers = workers,
    chunk_size = chunk_size,
    internal_flags = internal_flags
  )
}

sp_apply <- function(X, FUN, decision, packages = NULL, seed = TRUE, initializer = NULL) {
  backend <- decision$backend
  packages <- unique(c(packages %||% character(), "reticulate"))  # ensure reticulate on workers

  if (identical(backend, "sequential")) {
    if (is.function(initializer)) initializer()
    return(lapply(X, FUN))
  }

  if (identical(backend, "base_parallel")) {
    sys_name <- Sys.info()[["sysname"]]
    if (sys_name != "Windows" && .Platform$OS.type == "unix") {
      # Forked → inherits Python, but calling initializer is harmless.
      return(parallel::mclapply(X, function(x) { if (is.function(initializer)) initializer(); FUN(x) },
                                mc.cores = max(1L, decision$workers)))
    } else {
      cl <- parallel::makePSOCKcluster(max(1L, decision$workers))
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      if (length(packages)) {
        parallel::clusterCall(cl, function(pkgs) lapply(pkgs, requireNamespace, quietly = TRUE), packages)
      }
      if (is.function(initializer)) parallel::clusterCall(cl, initializer)
      return(parallel::parLapply(cl, X, FUN))
    }
  }

  if (identical(backend, "foreach")) {
    if (!requireNamespace("foreach", quietly = TRUE)) stop("foreach package required for backend='foreach'")
    if (!requireNamespace("doParallel", quietly = TRUE)) stop("doParallel package required for backend='foreach'")
    cl <- parallel::makePSOCKcluster(max(1L, decision$workers))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
    if (is.function(initializer)) parallel::clusterCall(cl, initializer)
    `%dopar%` <- foreach::`%dopar%`
    return(
      foreach::foreach(i = X, .packages = packages) %dopar% {
        FUN(i)  # initializer already ran via clusterCall
      }
    )
  }

  # future backend
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    stop("future/future.apply packages required for backend='future'")
  }
  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
  future::plan(decision$plan, workers = max(1L, decision$workers))
  innerFUN <- if (is.function(initializer)) {
    function(x) { initializer(); FUN(x) }
  } else FUN
  future.apply::future_lapply(
    X, innerFUN,
    future.seed = isTRUE(seed),
    future.chunk.size = decision$chunk_size,
    future.packages = packages
  )
}

# sp_apply <- function(X, FUN, decision, packages = NULL, seed = TRUE) {
#   backend <- decision$backend
#
#   if (identical(backend, "sequential")) {
#     return(lapply(X, FUN))
#   }
#
#   if (identical(backend, "base_parallel")) {
#     sys_name <- Sys.info()[["sysname"]]
#     if (sys_name != "Windows" && .Platform$OS.type == "unix") {
#       return(parallel::mclapply(X, FUN, mc.cores = max(1L, decision$workers)))
#     } else {
#       cl <- parallel::makePSOCKcluster(max(1L, decision$workers))
#       on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
#       if (!is.null(packages) && length(packages)) {
#         parallel::clusterCall(cl, function(pkgs) lapply(pkgs, requireNamespace, quietly = TRUE), packages)
#       }
#       return(parallel::parLapply(cl, X, FUN))
#     }
#   }
#
#   if (identical(backend, "foreach")) {
#     if (!requireNamespace("foreach", quietly = TRUE)) stop("foreach package required for backend='foreach'")
#     if (!requireNamespace("doParallel", quietly = TRUE)) stop("doParallel package required for backend='foreach'")
#     cl <- parallel::makePSOCKcluster(max(1L, decision$workers))
#     on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
#     doParallel::registerDoParallel(cl)
#     `%dopar%` <- foreach::`%dopar%`
#     return(foreach::foreach(i = X, .packages = packages %||% character()) %dopar% FUN(i))
#   }
#
#   if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
#     stop("future/future.apply packages required for backend='future'")
#   }
#   old_plan <- future::plan()
#   on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
#   future::plan(decision$plan, workers = max(1L, decision$workers))
#   future.apply::future_lapply(
#     X, FUN,
#     future.seed = isTRUE(seed),
#     future.chunk.size = decision$chunk_size,
#     future.packages = packages %||% character()
#   )
# }

# =========================
# Model name canonicalization (friendly ↔ canonical)
# =========================
map_to_canonical <- function(x, canonical_names, friendly_names) {
  if (is.null(x) || !length(x)) return(character())
  key  <- tolower(c(canonical_names, friendly_names))
  val  <- c(canonical_names, canonical_names)
  LUT  <- stats::setNames(val, key)
  out  <- vapply(x, function(s) {
    k <- tolower(as.character(s))
    if (!is.na(LUT[k])) LUT[k] else as.character(s)
  }, character(1))
  unname(out)
}

# Orientation tables (extend if you use more names)
maximize_metrics <- c("accuracy","kendalls_tau","spearman","pearson","correlation","rsq","r2")
minimize_metrics <- c("mean_squared_error", "bias",
                      "root_mean_squared_error", "relative_squared_error",
                      "mean_absolute_error", "mean_absolute_percent_error")

# Worst-case numeric sentinels (finite)
metric_worst_value <- function(metric, y_true = NULL) {
  metric <- tolower(metric)
  if (metric %in% maximize_metrics) return(-1e12)        # very bad for "maximize"
  if (metric %in% minimize_metrics) {
    # scale worst value ~ variance when available, else big constant
    v <- tryCatch(stats::var(y_true, na.rm = TRUE), error = function(e) NA_real_)
    v <- if (is.finite(v) && v > 0) v else 1
    return(1e12 * v)                                       # very bad for "minimize"
  }
  # Unknown metric: play safe and return big penalty
  1e12
}

# Robust metric that NEVER returns NA/Inf; falls back to worst-case finite value
safe_metric_value <- function(y_true, y_pred, metric) {
  # keep finite pairs
  idx <- is.finite(y_true) & is.finite(y_pred)
  if (!any(idx)) return(metric_worst_value(metric, y_true))
  yt <- y_true[idx]; yp <- y_pred[idx]

  needs_corr <- tolower(metric) %in% c("kendalls_tau","accuracy", "spearman","pearson","correlation","rsq","r2")
  # too few points or zero variance → worst value
  if (length(yt) < if (needs_corr) 2L else 1L) return(metric_worst_value(metric, y_true))
  if (needs_corr && (stats::sd(yt) == 0 || stats::sd(yp) == 0)) return(metric_worst_value(metric, y_true))

  # compute; if anything errors, return worst-case
  out <- tryCatch(
    evaluation_metrics(y_observed = yt, y_predicted = yp, eval_metrics = metric),
    error = function(e) NA_real_
  )
  if (!is.finite(out)) out <- metric_worst_value(metric, y_true)
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
# ---------- Robust hold-out maker with retries ----------
# Create strata for regression: 10 quantile bins (works for "stratified" hold-out)
ho_strata <- function(pd, trait, sampling_method) {
  if (!grepl("stratified", sampling_method %||% "", ignore.case = TRUE)) return(NULL)
  y_raw <- pd[[trait]]
  if (is.factor(y_raw) || (is.character(y_raw) && length(unique(y_raw)) <= 20)) {
    return(as.integer(factor(y_raw)))
  }
  y <- suppressWarnings(as.numeric(y_raw))
  nunq <- length(unique(stats::na.omit(y)))
  if (!is.finite(nunq) || nunq < 10L) {
    return(as.integer(factor(y_raw)))
  }
  qs <- unique(stats::quantile(y, probs = seq(0, 1, length.out = 11), na.rm = TRUE, type = 8))
  cut(y, qs, include.lowest = TRUE, labels = FALSE)
}

ho_valid <- function(y, tst, min_test = 5L) {
  if (!length(tst) || length(tst) < min_test) return(FALSE)
  if (length(tst) >= length(y)) return(FALSE)
  s <- suppressWarnings(stats::sd(y[tst]))
  isTRUE(is.finite(s) && s > 0)
}

# Try up to max_attempts to get a test split with enough size and non-zero variance
make_holdout <- function(y, strata = NULL, test_size = 0.2, seed = NULL,
                         min_test = 5L, max_attempts = 30L) {
  n <- length(y); idx_all <- seq_len(n); test_size <- as.numeric(test_size %||% 0.2)
  for (att in seq_len(max_attempts)) {
    if (!is.null(seed)) set.seed(seed + att - 1L)
    if (!is.null(strata)) {
      # per-stratum sampling
      split_ix <- split(idx_all, strata)
      tst <- unlist(lapply(split_ix, function(ix) {
        if (!length(ix)) return(integer(0))
        k <- max(1L, floor(length(ix) * test_size))
        if (length(ix) <= k) ix else sample(ix, size = k)
      }), use.names = FALSE)
    } else {
      k <- max(1L, floor(n * test_size))
      tst <- sample(idx_all, size = k)
    }
    tst <- sort(unique(tst))
    if (ho_valid(y, tst, min_test = min_test)) {
      return(list(test = tst, train = setdiff(idx_all, tst), attempts = att))
    }
  }
  # Fallback: simple split (may still be degenerate; your metric guard should handle it)
  if (!is.null(seed)) set.seed(seed + max_attempts + 1L)
  k   <- max(1L, floor(n * test_size))
  tst <- sort(sample(idx_all, size = k))
  list(test = tst, train = setdiff(idx_all, tst), attempts = max_attempts)
}

# Wrapper: try your existing hold_out_stratified_and_un(), fall back to robust splitter
get_holdout_indices <- function(pheno_data, gen_name, trait, test_size, seed, sampling_method,
                                min_test = 5L, verbose = FALSE) {
  y <- suppressWarnings(as.numeric(pheno_data[[trait]]))
  # try the package splitter first
  out <- try(hold_out_stratified_and_un(
    pheno_data = pheno_data, gen_name = gen_name, response = trait,
    test_size = test_size, random_state = seed, replication = 1L,
    sampling_method = sampling_method
  ), silent = TRUE)

  if (!inherits(out, "try-error") && length(out) && length(out[[1]]) > 0) {
    tst <- out[[1]]
    if (ho_valid(y, tst, min_test = min_test)) return(out)
    if (isTRUE(verbose)) message(sprintf("[holdout] fallback: invalid test split (n=%d, sd=%.4f)",
                                         length(tst), suppressWarnings(stats::sd(y[tst]))))
  } else if (isTRUE(verbose)) {
    message("[holdout] fallback: package splitter errored; switching to robust sampler.")
  }

  strata <- ho_strata(pheno_data, trait, sampling_method)
  sp <- make_holdout(y, strata = strata, test_size = test_size, seed = seed,
                      min_test = min_test, max_attempts = 30L)
  if (isTRUE(verbose)) message(sprintf("[holdout] robust sampler attempts=%d, test=%d", sp$attempts, length(sp$test)))
  list(sp$test)
}

# Assign predictions to test indices robustly (length or name alignment)
assign_preds <- function(ypred_cv, tst, preds, pheno_data, verbose = FALSE) {
  if (is.null(preds)) return(ypred_cv)
  if (length(preds) == length(tst)) {
    ypred_cv$yhat[tst] <- as.double(preds)
    return(ypred_cv)
  }
  idx <- NULL
  if (!is.null(names(preds)) && !is.null(rownames(pheno_data))) {
    idx <- match(names(preds), rownames(pheno_data))
  }
  if (!is.null(idx) && any(!is.na(idx))) {
    ok <- !is.na(idx)
    ypred_cv$yhat[idx[ok]] <- as.double(preds[ok])
    return(ypred_cv)
  }
  m <- min(length(preds), length(tst))
  if (m > 0) ypred_cv$yhat[tst[seq_len(m)]] <- as.double(preds[seq_len(m)])
  if (isTRUE(verbose)) message(sprintf("[holdout] partial assignment: %d of %d", m, length(tst)))
  ypred_cv
}


# ---------- SAFE METRICS LAYER (paste once at top of models_execute_crossval) ----------

# capture the original evaluation_metrics function (from your package or global env)
orig_eval_metrics <- tryCatch(get("evaluation_metrics", mode = "function", inherits = TRUE),
                              error = function(e) NULL)

metric_worst_value <- function(metric, y_true = NULL) {
  m <- tolower(metric)
  if (m %in% maximize_metrics) return(-1e12)  # worst for "maximize"
  if (m %in% minimize_metrics) {
    v <- tryCatch(stats::var(y_true, na.rm = TRUE), error = function(e) NA_real_)
    v <- if (is.finite(v) && v > 0) v else 1
    return(1e12 * v)                           # worst for "minimize"
  }
  1e12
}

safe_metric_value <- local({
  # bind the original into the closure so we don't recurse when we shadow the name below
  orig <- orig_eval_metrics
  function(y_true, y_pred, metric) {
    # keep finite pairs only
    idx <- is.finite(y_true) & is.finite(y_pred)
    if (!any(idx)) return(metric_worst_value(metric, y_true))
    yt <- y_true[idx]; yp <- y_pred[idx]

    # metrics that truly need correlation-like conditions (NOT accuracy)
    needs_corr <- tolower(metric) %in% c("kendalls_tau","spearman","pearson","correlation","rsq","r2")
    if (length(yt) < if (needs_corr) 2L else 1L) return(metric_worst_value(metric, y_true))
    if (needs_corr && (stats::sd(yt) == 0 || stats::sd(yp) == 0)) return(metric_worst_value(metric, y_true))

    # compute via the original if present; otherwise return a neutral penalty
    out <- tryCatch(
      if (is.function(orig)) orig(y_observed = yt, y_predicted = yp, eval_metrics = metric) else NA_real_,
      error = function(e) NA_real_
    )
    if (!is.finite(out)) out <- metric_worst_value(metric, y_true)
    out
  }
})

# CRITICAL: shadow any direct calls to evaluation_metrics() inside this function scope
evaluation_metrics <- function(y_observed, y_predicted, eval_metrics, ...) {
  safe_metric_value(y_true = y_observed, y_pred = y_predicted, metric = eval_metrics)
}
# ---------- END SAFE METRICS LAYER ----------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Is hold-out viable for this dataset/config?
is_holdout_feasible <- function(pheno_data, response, sampling_method, test_size,
                                min_total_n = getOption("gp.min_total_n_holdout", 40L),
                                min_test = NULL, min_train = NULL,
                                min_per_stratum = 2L, max_bins = 5L) {
  n <- NROW(pheno_data)
  if (n < min_total_n) return(list(ok=FALSE, reason=sprintf("n=%d < min_total_n=%d", n, min_total_n)))
  if (is.null(min_test))  min_test  <- max(5L, ceiling(0.05 * n))
  if (is.null(min_train)) min_train <- max(5L, ceiling(0.10 * n))
  exp_test  <- floor(n * test_size)
  exp_train <- n - exp_test
  if (exp_test  < min_test)  return(list(ok=FALSE, reason=sprintf("expected test=%d < min_test=%d", exp_test,  min_test)))
  if (exp_train < min_train) return(list(ok=FALSE, reason=sprintf("expected train=%d < min_train=%d", exp_train, min_train)))

  y <- pheno_data[[response]]
  if (grepl("stratified", sampling_method %||% "", ignore.case=TRUE)) {
    # Build strata
    if (is.factor(y) || (is.character(y) && length(unique(y)) <= 20L)) {
      groups <- as.integer(factor(y, exclude=NULL))
    } else {
      ynum <- suppressWarnings(as.numeric(y))
      bins <- max(2L, min(max_bins, floor(length(unique(ynum[is.finite(ynum)]))/2)))
      qs   <- unique(quantile(ynum, probs = seq(0,1,length.out=bins+1), na.rm=TRUE))
      if (length(qs) <= 2L) return(list(ok=FALSE, reason="cannot form stratification bins"))
      groups <- findInterval(ynum, qs, all.inside=TRUE)
    }
    tab <- as.integer(table(groups))
    # Need enough per stratum for both train & test
    req_test  <- ceiling(min_per_stratum / test_size)
    req_train <- ceiling(min_per_stratum / (1 - test_size))
    if (any(tab < pmax(req_test, req_train))) {
      return(list(ok=FALSE, reason="some strata too small for requested test_size"))
    }
  }
  list(ok=TRUE, reason="ok")
}

# Choose a safe nfolds for (stratified) K-fold when hold-out is not viable
pick_safe_kfold <- function(pheno_data, response, sampling_method,
                            nfolds_requested = 5L, min_test = NULL,
                            min_per_stratum = 2L, max_bins = 5L) {
  n <- NROW(pheno_data)
  if (is.null(min_test)) min_test <- max(5L, ceiling(0.05 * n))
  nfolds <- min(nfolds_requested, max(2L, floor(n / min_test)))  # each fold’s test size >= min_test

  if (grepl("stratified", sampling_method %||% "", ignore.case=TRUE)) {
    y <- pheno_data[[response]]
    if (is.factor(y) || (is.character(y) && length(unique(y)) <= 20L)) {
      groups <- as.integer(factor(y, exclude=NULL))
    } else {
      ynum <- suppressWarnings(as.numeric(y))
      bins <- max(2L, min(max_bins, floor(length(unique(ynum[is.finite(ynum)]))/2)))
      qs   <- unique(quantile(ynum, probs = seq(0,1,length.out=bins+1), na.rm=TRUE))
      groups <- if (length(qs) > 2L) findInterval(ynum, qs, all.inside=TRUE) else rep(1L, n)
    }
    tab <- as.integer(table(groups))
    # Each stratum must have at least nfolds rows to distribute one per fold
    nfolds <- min(nfolds, max(2L, min(tab)))
  }
  max(2L, nfolds)
}


safe_num <- function(x, digits = 2, na = "NA") {
  if (length(x) == 0L || is.null(x) || !is.numeric(x) || !is.finite(x)) return(na)
  formatC(x, format = "f", digits = digits)
}

#' Title
#'
#' @param model
#' @param y
#' @param omics_data
#' @param tst
#' @param additional_params
#'
#' @return
#' @export
#'
#' @examples
predict_with_model <- function(model = NULL,
                               y = NULL,
                               omics_data = NULL,
                               tst = NULL,
                               additional_params = NULL) {
  # Generalized function to handle predictions for various models

  switch(model,
         "Xgboost" = AI_xgboost_cv(y = y, omics = omics_data, tst = tst, eta = additional_params$eta,
                                   nrounds = additional_params$nrounds, max_depth = additional_params$max_depth,
                                   scaling = additional_params$scaling, omic_count = additional_params$omic_count,
                                   centering = additional_params$centering, xgb_gamma = additional_params$xgb_gamma,
                                   colsample_bytree = additional_params$colsample_bytree,
                                   subsample = additional_params$subsample, min_child_weight = additional_params$min_child_weight,
                                   xgb_alpha = additional_params$xgb_alpha, xgb_lambda = additional_params$xgb_lambda,
                                   xgb_booster = additional_params$xgb_booster,
                                   xgb_rate_drop = additional_params$xgb_rate_drop,xgb_skip_drop = additional_params$xgb_skip_drop,
                                   xgb_objective = additional_params$xgb_objective,xgb_sample_type = additional_params$xgb_sample_type,
                                   xgb_normalize_type = additional_params$xgb_normalize_type,
                                   early_stop_for_iteration_xgb = additional_params$early_stop_for_iteration_xgb),
         #unified deep learning branch
         "cnn" = ,
         "resnet" = ,
         "ft_transformer" = ,
         "saint" = ,
         "tabnet" = ,
         "node" = ,
         "deepfm" = ,
         "dcnv2" = ,
         "nam" = ,
         "moe" = ,
         "gp_dkl" = ,
         "mlp_with_attention" = ,
         "mlp" = deep_learning_model(
           y = y,
           omics = omics_data,
           tst = tst,
           scaling = additional_params$scaling,
           centering = additional_params$centering,
           crossval = additional_params$crossval,
           omic_count = additional_params$omic_count,
           early_stop = additional_params$early_stop,
           deep_learning_model = additional_params$deep_learning_model,

           optimizer_name = additional_params$optimizer_name,
           use_amp = additional_params$use_amp,
           max_grad_norm = additional_params$max_grad_norm,
           auto_class_weights = additional_params$auto_class_weights,

           ##### CNN
           cnn_neurons_per_layer = additional_params$cnn_neurons_per_layer,
           cnn_kernel_size = additional_params$cnn_kernel_size,
           cnn_dense_layers = additional_params$cnn_dense_layers,
           cnn_use_max_pool = additional_params$cnn_use_max_pool,
           cnn_pool_kernel = additional_params$cnn_pool_kernel,
           cnn_pool_stride = additional_params$cnn_pool_stride,
           cnn_pool_padding = additional_params$cnn_pool_padding,
           cnn_learning_rate = additional_params$cnn_learning_rate,
           cnn_separable = additional_params$cnn_separable,
           cnn_dilations = additional_params$cnn_dilations,
           cnn_use_se = additional_params$cnn_use_se,
           cnn_norm_type = additional_params$cnn_norm_type,
           cnn_pool_type = additional_params$cnn_pool_type,
           cnn_use_global_pool = additional_params$cnn_use_global_pool,

           ##### ResNet
           resnet_neurons_per_block = additional_params$resnet_neurons_per_block,
           resnet_blocks = additional_params$resnet_blocks,
           resnet_learning_rate = additional_params$resnet_learning_rate,

           ##### FT Transformer
           ft_d_model = additional_params$ft_d_model,
           ft_heads = additional_params$ft_heads,
           ft_layers = additional_params$ft_layers,
           ft_ff_mult = additional_params$ft_ff_mult,
           ft_dropout = additional_params$ft_dropout,
           ft_token_dropout = additional_params$ft_token_dropout,
           ft_use_cls = additional_params$ft_use_cls,

           ##### SAINT
           saint_d_model = additional_params$saint_d_model,
           saint_heads = additional_params$saint_heads,
           saint_layers = additional_params$saint_layers,
           saint_ff_mult = additional_params$saint_ff_mult,
           saint_dropout = additional_params$saint_dropout,
           saint_token_dropout = additional_params$saint_token_dropout,
           saint_use_cls = additional_params$saint_use_cls,

           ##### Grouping controls
           use_grouping = additional_params$use_grouping,
           group_trigger = additional_params$group_trigger,
           group_method = additional_params$group_method,
           init_group_size = additional_params$init_group_size,
           max_tokens = additional_params$max_tokens,
           kmeans_batch = additional_params$kmeans_batch,
           kmeans_iter = additional_params$kmeans_iter,

           ##### TabNet
           tabnet_steps = additional_params$tabnet_steps,
           tabnet_feature_dim = additional_params$tabnet_feature_dim,
           tabnet_output_dim = additional_params$tabnet_output_dim,
           tabnet_gamma = additional_params$tabnet_gamma,
           tabnet_lambda_sparse = additional_params$tabnet_lambda_sparse,

           ##### NODE
           node_trees = additional_params$node_trees,
           node_depth = additional_params$node_depth,

           ##### DeepFM
           deepfm_k = additional_params$deepfm_k,
           deepfm_hidden = additional_params$deepfm_hidden,

           ##### DCNv2
           dcn_layers = additional_params$dcn_layers,
           dcn_hidden = additional_params$dcn_hidden,

           ##### NAM
           nam_hidden = additional_params$nam_hidden,
           nam_activation = additional_params$nam_activation,
           nam_add_linear = additional_params$nam_add_linear,
           nam_l1 = additional_params$nam_l1,

           ##### Mixture-of-Experts
           moe_n_experts = additional_params$moe_n_experts,
           moe_expert_hidden = additional_params$moe_expert_hidden,
           moe_gate_hidden = additional_params$moe_gate_hidden,
           moe_temperature = additional_params$moe_temperature,
           moe_sparse_topk = additional_params$moe_sparse_topk,
           moe_entropy_reg = additional_params$moe_entropy_reg,

           ##### GP / RFF
           gp_use_variational = additional_params$gp_use_variational,
           gp_num_inducing = additional_params$gp_num_inducing,
           gp_feature_dim = additional_params$gp_feature_dim,
           gp_kernel = additional_params$gp_kernel,
           gp_ard = additional_params$gp_ard,
           gp_lr_mult = additional_params$gp_lr_mult,
           rff_features = additional_params$rff_features,
           rff_lengthscale = additional_params$rff_lengthscale,
           rff_deep_hidden = additional_params$rff_deep_hidden,

           ##### General deep learning params
           model_type = additional_params$model_type,
           epochs = additional_params$epochs,
           batch_size = additional_params$batch_size,
           dropout = additional_params$dropout,
           l2_weight_decay = additional_params$l2_weight_decay,
           l2_regularizer_dp = additional_params$l2_regularizer_dp,
           dropout_rate = additional_params$dropout_rate,
           batch_norm = additional_params$batch_norm,
           validation_split = additional_params$validation_split,
           compile_model = additional_params$compile_model,
           deterministic = additional_params$deterministic,
           random_seed = additional_params$random_seed,
           device = additional_params$device,

           ##### MLP and Attention
           mlp_neurons_per_layer = additional_params$mlp_neurons_per_layer,
           mlp_learning_rate = additional_params$mlp_learning_rate,
           final_attention = additional_params$final_attention,
           attention_across_multiple_layers = additional_params$attention_across_multiple_layers,
           heteroscedastic = additional_params$heteroscedastic
         ),
         "RandomForest" = AI_randomforest_cv(y = y, omics = omics_data, tst = tst,
                                             scaling = additional_params$scaling,
                                             centering = additional_params$centering, ntree = additional_params$ntree,
                                             omic_count = additional_params$omic_count),
         "PartialLeastSquare" = AI_pls_cv(y = y, omics = omics_data, tst = tst,
                                          scaling = additional_params$scaling,
                                          centering = additional_params$centering, ncomp = additional_params$ncomp,
                                          omic_count = additional_params$omic_count),
         "Ridge_Regression" = AI_ridge_regression_cv(y = y, omics = omics_data, tst = tst,
                                                     scaling = additional_params$scaling,
                                                     centering = additional_params$centering,
                                                     omic_count = additional_params$omic_count),
         "Lasso" = AI_lasso_cv(y = y, omics = omics_data, tst = tst, scaling = additional_params$scaling,
                               centering = additional_params$centering,
                               omic_count = additional_params$omic_count),
         "SupportVectorMachine" = AI_svm_cv(y = y, omics = omics_data, tst = tst,
                                            scaling = additional_params$scaling,
                                            centering = additional_params$centering,
                                            C_value = additional_params$C_value,
                                            degree_value = additional_params$degree_value,
                                            scale_value = additional_params$scale_value,
                                            offset_value = additional_params$offset_value,
                                            omic_count = additional_params$omic_count),
         "K-NearestNeighbors" = AI_knn_cv(y = y, omics = omics_data, tst = tst,
                                          scaling = additional_params$scaling,
                                          centering = additional_params$centering, k = additional_params$k,
                                          omic_count = additional_params$omic_count),
         "GBLUP" = asreml_mod_cv(asreml_models_prep_cv = additional_params$asreml_models_prep_cv,
                                 pheno_data = additional_params$pheno_data,
                                 response = additional_params$response,
                                 heter_groups = additional_params$heter_groups,
                                 gen_name = additional_params$gen_name, tst = tst),
         "Bayes" = bayes_mod_cv(y = y, ETA = additional_params$ETA, weights = additional_params$weights,
                                bayes_para = additional_params$bayes_para, tst = tst,
                                bayes_model = additional_params$bayes_model, bayes_trait = additional_params$bayes_trait)
  )
}

full_dp_models <- setNames(
  c("Conv1DNet", "TabTransformer", "TabAttention", "TabNet", "LightTreeNet",
    "FactorNet", "CrossNet", "NeuralAdditive", "MixtureOfExperts", "GPNet",
    "DenseAttentionNet", "DenseNeuralNet", "ResNet"),
  c("cnn", "ft_transformer", "saint", "tabnet", "node",
    "deepfm", "dcnv2", "nam", "moe", "gp_dkl", "mlp_with_attention", "mlp", "resnet")

)

# robust CV token normalizer  # <<< CHANGED
normalize_cv_token <- function(x) {
  if (is.null(x) || !nzchar(x)) return(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", x)
  x <- sub("_$", "", x)
  x
}

# canonicalize model names (friendly -> canonical)  # <<< CHANGED
map_to_canonical <- function(x, canonical, friendly) {
  if (is.null(x) || !length(x)) return(character())
  key <- tolower(c(canonical, friendly))
  val <- c(canonical, canonical)
  lut <- stats::setNames(val, key)
  unname(vapply(x, function(s) lut[[tolower(as.character(s))]] %||% as.character(s), character(1)))
}

# ======================================================================
models_execute_crossval <- function(pheno_data = NULL,
                                    test_set = NULL,
                                    response = NULL,
                                    gen_name = NULL,
                                    test_size = NULL,
                                    random_state = NULL,
                                    replication = NULL,
                                    weights = NULL,
                                    selected_raw,
                                    gam_method = NULL,
                                    max_features = 50,
                                    k_value = 5,
                                    var_explained = 0.9,
                                    engine = NULL,
                                    model_prep_all_bayes_cv = NULL,
                                    asreml_models_prep_cv = NULL,
                                    ml_dat_res = NULL,
                                    heter_groups = NULL,
                                    verbose = FALSE,
                                    num_cores = NULL,
                                    nfolds = 5,
                                    cross_validation_meth = NULL,
                                    sampling_method = NULL,
                                    eval_metrics = NULL,
                                    bayes_model = NULL,
                                    GS_model_cv = NULL,
                                    scaling = FALSE,
                                    centering = TRUE,
                                    eta = 0.1,
                                    nrounds = 100,
                                    max_depth = 6,
                                    xgb_gamma = 4,
                                    subsample = 0.5,
                                    colsample_bytree = 1,
                                    xgb_alpha = 0.001,
                                    xgb_lambda = 1,
                                    min_child_weight = 1,
                                    early_stop_for_iteration_xgb = TRUE,
                                    xgb_booster = "dart",
                                    xgb_rate_drop = 0.1,
                                    xgb_skip_drop = 0.5,
                                    xgb_objective = "reg:squarederror",
                                    xgb_sample_type = "uniform",
                                    xgb_normalize_type = "tree",
                                    ncomp = 3,
                                    ntree = 500,
                                    k = 5,
                                    svm_kernel = "Gaussian",
                                    sigma_value  = 0.1,
                                    C_value  = 1,
                                    degree_value = 3,
                                    scale_value  = 1,
                                    offset_value = 1,
                                    num_hidden_layers = 1,
                                    neurons_per_layer = 64,
                                    learning_rate_dp = 0.001,
                                    early_stop = TRUE,
                                    crossval = TRUE,
                                    optimizer_name = "adam",
                                    use_amp        = TRUE,
                                    max_grad_norm  = 1.0,
                                    auto_class_weights = FALSE,
                                    ##### cnn
                                    cnn_neurons_per_layer = as.integer(c(64, 64, 64)),
                                    cnn_kernel_size = 3L,
                                    cnn_dense_layers = as.integer(c(256, 128, 64)),
                                    cnn_use_max_pool = FALSE,
                                    cnn_pool_kernel = 2L,
                                    cnn_pool_stride = 2L,
                                    cnn_pool_padding = 0L,
                                    cnn_learning_rate = 1e-3,
                                    cnn_separable=TRUE,
                                    cnn_dilations=c(1,2,4),
                                    cnn_use_se=TRUE,
                                    cnn_norm_type="group",
                                    cnn_pool_type="conv",
                                    cnn_use_global_pool=FALSE,
                                    ##### resnet
                                    resnet_neurons_per_block = as.integer(c(256, 128, 64)),
                                    resnet_blocks = 3,
                                    resnet_learning_rate = 1e-3,
                                    #### ft_transformer
                                    ft_d_model = 192L,
                                    ft_heads = 8L,
                                    ft_layers = 3L,
                                    ft_ff_mult = 4L,
                                    ft_dropout = 0.1,
                                    ft_token_dropout = 0.0,
                                    ft_use_cls = TRUE,
                                    #### saint
                                    saint_d_model = 128L,
                                    saint_heads = 8L,
                                    saint_layers = 3L,
                                    saint_ff_mult = 4L,
                                    saint_dropout = 0.1,
                                    saint_token_dropout = 0.0,
                                    saint_use_cls       = TRUE,
                                    ###### Grouping controls
                                    use_grouping   = FALSE,
                                    group_trigger  = 2048,
                                    group_method   = "auto",
                                    init_group_size = 64,
                                    max_tokens      = 1024,
                                    kmeans_batch    = 4096,
                                    kmeans_iter     = 100,
                                    #### tabnet
                                    tabnet_steps = 5L,
                                    tabnet_feature_dim = 64L,
                                    tabnet_output_dim = 64L,
                                    tabnet_gamma = 1.5,
                                    tabnet_lambda_sparse = 1e-4,
                                    #### node
                                    node_trees = 8L,
                                    node_depth = 3L,
                                    #### deepfm
                                    deepfm_k = 16L,
                                    deepfm_hidden = as.integer(c(128, 64)),
                                    #### dcnv2
                                    dcn_layers = 3L,
                                    dcn_hidden = as.integer(c(256, 128, 64)),
                                    #### nam
                                    nam_hidden = as.integer(c(32, 16)),
                                    nam_activation = "relu",
                                    nam_add_linear = TRUE,
                                    nam_l1 = 1e-4,
                                    ### moe
                                    moe_n_experts = 4L,
                                    moe_expert_hidden = as.integer(c(128, 64)),
                                    moe_gate_hidden = 128L,
                                    moe_temperature = 1.0,
                                    moe_sparse_topk = NA,
                                    moe_entropy_reg = 0.0,
                                    ### gp_dkl/RFF knobs
                                    gp_use_variational = TRUE,
                                    gp_num_inducing = 256L,
                                    gp_feature_dim = 64L,
                                    gp_kernel = "rbf",
                                    gp_ard = TRUE,
                                    gp_lr_mult = 0.5,
                                    rff_features       = 1024,
                                    rff_lengthscale    = 1.0,
                                    rff_deep_hidden    = c(128),
                                    ##### General dp
                                    model_type = "resnet",
                                    epochs = 10,
                                    batch_size = 64 ,
                                    dropout = 0.2,
                                    l2_weight_decay = 1e-4,
                                    l2_regularizer_dp = 0.001,
                                    dropout_rate = 0.5,
                                    batch_norm = TRUE,
                                    validation_split = 0.2,
                                    compile_model = FALSE,
                                    deterministic = TRUE,
                                    random_seed = 123,
                                    device = NULL,
                                    #### mlp and attention
                                    mlp_neurons_per_layer = as.integer(c(128, 64)),
                                    mlp_learning_rate = 1e-3,
                                    final_attention = TRUE,
                                    attention_across_multiple_layers = TRUE,
                                    heteroscedastic = TRUE,
                                    docker_nd_usage = FALSE,
                                    globals_max_GB = 4,
                                    # Parallel controls
                                    parallel_mode = c("auto","future","sequential","base_parallel","foreach"),
                                    parallel_backend_prefer_fork = TRUE,
                                    # NEW: user-controlled sequential models
                                    sequential_models = NULL,
                                    ...) {

  on.exit(try(future::plan("sequential"), silent = TRUE), add = TRUE)
  parallel_mode <- match.arg(parallel_mode)

  # --- helpers --------------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # capture Python binding for reticulate workers (if you use Python)
  gp_py_bin <- if (exists("gp_detect_python", mode = "function")) gp_detect_python() else NULL
  init_py <- function() if (exists("gp_init_python_once", mode = "function")) gp_init_python_once(gp_py_bin)

  msg <- "\n==================================================\n"
  if (!is.null(cross_validation_meth) & length(cross_validation_meth) > 1) {
    stop(paste(msg, 'use only one cross_validation method at a time.'), call. = FALSE)
  }

  # sampling defaults
  patterns <- c("stratified", "Repeated")
  if (is.null(sampling_method) | is.null(replication)) {
    matche_strings <- sapply(patterns, function(pattern) {
      length(grep(pattern, cross_validation_meth, ignore.case = TRUE)) > 0
    }, simplify = FALSE)
    names(matche_strings) <- patterns
    present_patterns <- names(matche_strings)[unlist(matche_strings)]
    if ("stratified" %in% present_patterns) sampling_method <- "stratified"
  }

  # data prep
  if (!is.null(ml_dat_res)) {
    omics_data <-  ml_dat_res[["merged_data"]][["merge_data"]]
    if (!"merged_data_test" %in% names(ml_dat_res) && !is.null(test_set)) {
      omics_data <- omics_data[!rownames(omics_data) %in% test_set, ]
    }
    omic_count <- if ("omic_count" %in% names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL
  }

  # parameters bag (unchanged)
  additional_params <- list(
    ETA = NULL, weights = weights, bayes_para = NULL, bayes_model = NULL, bayes_trait = NULL,
    scaling = scaling, centering = centering,
    eta = eta, nrounds = nrounds, max_depth = max_depth, xgb_gamma = xgb_gamma,
    early_stop_for_iteration_xgb = early_stop_for_iteration_xgb, colsample_bytree = colsample_bytree,
    subsample = subsample, ntree = ntree, xgb_alpha = xgb_alpha, xgb_lambda = xgb_lambda,
    min_child_weight = min_child_weight, xgb_booster = xgb_booster, xgb_rate_drop = xgb_rate_drop,
    xgb_skip_drop = xgb_skip_drop, xgb_objective = xgb_objective, xgb_sample_type = xgb_sample_type,
    xgb_normalize_type = xgb_normalize_type, ncomp = ncomp, C_value = C_value,
    degree_value = degree_value, scale_value = scale_value, offset_value = offset_value,
    k = k, omic_count = if (exists("omic_count")) omic_count else NULL,
    asreml_models_prep_cv = asreml_models_prep_cv, gen_name = gen_name, pheno_data = pheno_data,
    max_features = max_features, heter_groups = heter_groups, crossval = crossval, early_stop = early_stop,
    model_type = model_type, optimizer_name = optimizer_name, use_amp = use_amp, max_grad_norm = max_grad_norm,
    auto_class_weights = auto_class_weights,
    cnn_neurons_per_layer = cnn_neurons_per_layer, cnn_kernel_size = cnn_kernel_size,
    cnn_dense_layers = cnn_dense_layers, cnn_use_max_pool = cnn_use_max_pool, cnn_pool_kernel = cnn_pool_kernel,
    cnn_pool_stride = cnn_pool_stride, cnn_pool_padding = cnn_pool_padding, cnn_learning_rate = cnn_learning_rate,
    cnn_separable = cnn_separable, cnn_dilations = cnn_dilations, cnn_use_se = cnn_use_se, cnn_norm_type = cnn_norm_type,
    cnn_pool_type = cnn_pool_type, cnn_use_global_pool = cnn_use_global_pool,
    resnet_neurons_per_block = resnet_neurons_per_block, resnet_blocks = resnet_blocks,
    resnet_learning_rate = resnet_learning_rate,
    ft_d_model = ft_d_model, ft_heads = ft_heads, ft_layers = ft_layers, ft_ff_mult = ft_ff_mult,
    ft_dropout = ft_dropout, ft_token_dropout = ft_token_dropout, ft_use_cls = ft_use_cls,
    saint_d_model = saint_d_model, saint_heads = saint_heads, saint_layers = saint_layers, saint_ff_mult = saint_ff_mult,
    saint_dropout = saint_dropout, saint_token_dropout = saint_token_dropout, saint_use_cls = saint_use_cls,
    use_grouping = use_grouping, group_trigger = group_trigger, group_method = group_method,
    init_group_size = init_group_size, max_tokens = max_tokens, kmeans_batch = kmeans_batch, kmeans_iter = kmeans_iter,
    tabnet_steps = tabnet_steps, tabnet_feature_dim = tabnet_feature_dim, tabnet_output_dim = tabnet_output_dim,
    tabnet_gamma = tabnet_gamma, tabnet_lambda_sparse = tabnet_lambda_sparse,
    node_trees = node_trees, node_depth = node_depth,
    deepfm_k = deepfm_k, deepfm_hidden = deepfm_hidden,
    dcn_layers = dcn_layers, dcn_hidden = dcn_hidden,
    nam_hidden = nam_hidden, nam_activation = nam_activation, nam_add_linear = nam_add_linear, nam_l1 = nam_l1,
    moe_n_experts = moe_n_experts, moe_expert_hidden = moe_expert_hidden, moe_gate_hidden = moe_gate_hidden,
    moe_temperature = moe_temperature, moe_sparse_topk = moe_sparse_topk, moe_entropy_reg = moe_entropy_reg,
    gp_use_variational = gp_use_variational, gp_num_inducing = gp_num_inducing, gp_feature_dim = gp_feature_dim,
    gp_kernel = gp_kernel, gp_ard = gp_ard, gp_lr_mult = gp_lr_mult,
    rff_features = rff_features, rff_lengthscale = rff_lengthscale, rff_deep_hidden = rff_deep_hidden,
    epochs = epochs, batch_size = batch_size, dropout = dropout, l2_weight_decay = l2_weight_decay,
    l2_regularizer_dp = l2_regularizer_dp, dropout_rate = dropout_rate, batch_norm = batch_norm,
    validation_split = validation_split, compile_model = compile_model, deterministic = deterministic,
    random_seed = random_seed, device = device,
    mlp_neurons_per_layer = mlp_neurons_per_layer, mlp_learning_rate = mlp_learning_rate,
    final_attention = final_attention, attention_across_multiple_layers = attention_across_multiple_layers,
    heteroscedastic = heteroscedastic
  )

  # model families
  AI_valid_models <- c("Xgboost","RandomForest","PartialLeastSquare",
                       "SupportVectorMachine","K-NearestNeighbors","Lasso","Ridge_Regression")
  bayes_valid_models <- c("BRR","BayesA","BayesB","BayesC","BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR","RKHS")

  canonical_names <- c("cnn","ft_transformer","saint","tabnet","node","deepfm","dcnv2",
                       "nam","moe","gp_dkl","mlp_with_attention","mlp","resnet")
  friendly_names  <- c("Conv1DNet","TabTransformer","TabAttention","TabNet","LightTreeNet","FactorNet",
                       "CrossNet","NeuralAdditive","MixtureOfExperts","GPNet","DenseAttentionNet","DenseNeuralNet","ResNet")

  # always canonicalize inputs  # <<< CHANGED
  name_lookup <- setNames(
    rep(canonical_names, 2),
    c(canonical_names, friendly_names)
  )

  friendly_name_lookup <- setNames(
    rep(friendly_names, 2),
    c(friendly_names, canonical_names)
  )
  #GS_model_cv_can <- map_to_canonical(GS_model_cv, canonical_names, friendly_names)

  dp_models <-  resolve_model_name(GS_model_cv, name_lookup = name_lookup)
  if(!is.null(dp_models) & length(dp_models)>0){
    AI_valid_models <- c(AI_valid_models, dp_models)

    GS_model_cv_can <-  replace_with_canonical(GS_model_cv, canonical_names = canonical_names,
                                           friendly_names = friendly_names)

    GS_model_cv <- replace_with_friendly_name(GS_model_cv, friendly_name_lookup)
    # model_type <- replace_with_canonical(model_type, canonical_names = canonical_names,
    #                                      friendly_names = friendly_names)
  }else{
    GS_model_cv_can <- GS_model_cv
  }
  if(!is.null(sequential_models)){
  sequential_models <- map_to_canonical(sequential_models, canonical_names, friendly_names)
  }
  AI_valid_models <- unique(c(AI_valid_models, canonical_names))


  ###################
  # CV method tokens (transparent; no conversion)  # <<< CHANGED
  cv_token <- normalize_cv_token(cross_validation_meth)
  holds_tokens <- c("hold_out","stratified_hold_out","repeated_hold_out","repeated_stratified_hold_out")
  kfold_tokens <- c("k_folds","stratified_k_folds","repeated_k_folds","repeated_stratified_k_folds")
  multi_tokens <- c("cv0", "cv1","cv2","repeated_cv0", "repeated_cv1","repeated_cv2")

  is_holdout <- cv_token %in% holds_tokens
  is_kfold   <- cv_token %in% kfold_tokens
  is_multi   <- cv_token %in% multi_tokens
  if (!(is_holdout || is_kfold || is_multi)) {
    stop(paste0("Unsupported cross-validation method: ", cross_validation_meth,
                "\nAllowed: Hold-Out, K-Folds families, or CV0/CV1/CV2 variants."), call. = FALSE)
  }


  sys_name <- Sys.info()[["sysname"]]
  if (docker_nd_usage) sys_name <- "Windows"

  # tasks (keep friendly label for display)  # <<< CHANGED
  tasks_full <- expand.grid(response = response,
                            replication = seq_len(replication),
                            modell = GS_model_cv_can,
                            modell_label = GS_model_cv,
                            stringsAsFactors = FALSE)

  idx_seq_user <- which(tasks_full$modell %in% sequential_models)
  idx_rest     <- setdiff(seq_len(nrow(tasks_full)), idx_seq_user)

  globals_for_size <- list(
    pheno_data = pheno_data,
    omics_data = if (exists("omics_data")) omics_data else NULL,
    bayes = model_prep_all_bayes_cv,
    asreml = asreml_models_prep_cv
  )


  # per-task
  run_one_task <- function(task_row) {
    if (is.function(init_py)) init_py()
    trait       <- as.character(task_row$response)
    rep_i       <- as.integer(task_row$replication)
    model_can   <- as.character(task_row$modell)        # canonical for routing  # <<< CHANGED
    model_label <- as.character(task_row$modell_label)  # friendly for display  # <<< CHANGED
#####

    repp <- rep_i
    base_seed <- if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1)
      as.integer(random_state) else 123L
    new_seed <- (base_seed + rep_i * 10000L) %% .Machine$integer.max


    # splits
    if (is_holdout) {
      test_set_val <- hold_out_stratified_and_un(
        pheno_data = pheno_data, gen_name = gen_name, response = trait,
        test_size = test_size, random_state = new_seed, replication = repp,
        sampling_method = sampling_method
      )

    } else if (is_kfold) {
      test_set_val <- kfolds_stratified_un(
        pheno_data = pheno_data, gen_name = gen_name, response = trait,
        nfolds = nfolds, random_state = new_seed, replication = repp,
        sampling_method = sampling_method
      )
    } else { # CV1/CV2
      #CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
      CVi <- as.integer(sub("^.*?([0-2])$", "\\1", tolower(cross_validation_meth)))
      test_set_val <- CV0_CV1_CV2_for_multi_environment(
        pheno_data = pheno_data, gen_name = gen_name,
        heter_groups = heter_groups,
        CV = CVi, nfolds = nfolds,
        random_state = new_seed, replication = repp, message = verbose
      )
    }

    len_y <- nrow(pheno_data)
    y <- as.double(pheno_data[[trait]])
    handle_error <- FALSE

    ypred_cv <- data.frame(y = y, yhat = rep(NA_real_, len_y), cv_role = "train", check.names = FALSE)
    results_eval_metrics_reps <- NULL

    if (is_multi) {
      if (is.null(heter_groups)) stop(paste(msg, "For CV1/CV2 you must supply `heter_groups`."), call. = FALSE)
      ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
    }

    # helper to assign predictions safely per-fold  # <<< CHANGED
    assign_preds <- function(tst, preds) {
      if (is.null(preds)) return()
      # try length match first
      if (length(preds) == length(tst)) {
        ypred_cv$yhat[tst] <<- as.double(preds)
        return()
      }
      # try name/rownames alignment
      idx <- NULL
      if (!is.null(names(preds))) {
        nm <- names(preds)
        if (!is.null(rownames(pheno_data))) {
          idx <- match(nm, rownames(pheno_data))
        }
      }
      if (!is.null(idx) && any(!is.na(idx))) {
        ypred_cv$yhat[idx[!is.na(idx)]] <<- as.double(preds[!is.na(idx)])
        return()
      }
      # fallback: partial assign
      m <- min(length(preds), length(tst))
      if (m > 0) ypred_cv$yhat[tst[seq_len(m)]] <<- as.double(preds[seq_len(m)])
      if (isTRUE(verbose)) message(sprintf("[cv] warning: preds length %d != tst length %d; assigned %d", length(preds), length(tst), m))
    }

    # ---- Hold-Out
    if (is_holdout) {
      tst <- test_set_val[[repp]]
      yNA <- y; yNA[tst] <- NA
      ypred_cv$cv_role[tst] <- "test"

      if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
        tryCatch({
          if (model_can == "GBLUP_BRR") {
            mGB <- "BRR"
            additional_params$bayes_model <- mGB
            additional_params$bayes_trait <- trait
            additional_params$ETA        <- model_prep_all_bayes_cv[[mGB]][["bayes_ETA"]][["ETA"]]
            additional_params$bayes_para  <- model_prep_all_bayes_cv[[mGB]][["bayes_para"]]
          } else {
            additional_params$bayes_model <- model_can
            additional_params$bayes_trait <- trait
            additional_params$ETA        <- model_prep_all_bayes_cv[[model_can]][["bayes_ETA"]][["ETA"]]
            additional_params$bayes_para  <- model_prep_all_bayes_cv[[model_can]][["bayes_para"]]
          }
          assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = additional_params))
        }, error = function(e) { message(paste("Error in Bayes", model_can, "for", trait, ":", e$message)); handle_error <<- TRUE })
      }

      if (!is.null(engine) && model_can == "GBLUP" && engine == "asreml") {
        tryCatch({
          additional_params$response <- trait
          assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = additional_params))
        }, error = function(e) { message(paste("Error in GBLUP(asreml) for", trait, ":", e$message)); handle_error <<- TRUE })
      }

      if (model_can %in% AI_valid_models) {
        tryCatch({
          assign_preds(tst, predict_with_model(
            model = model_can, y = yNA,
            omics_data = if (exists("omics_data")) omics_data else NULL,
            tst = tst, additional_params = additional_params))
        }, error = function(e) { message(paste("Error in AI model", model_can, "for", trait, ":", e$message)); handle_error <<- TRUE })
      }

      #test_idx <- which(ypred_cv$cv_role == "test" & !is.na(ypred_cv$yhat))
      test_idx <- which(ypred_cv$cv_role == "test")
      y_t <- ypred_cv$y[test_idx]
      y_p <- ypred_cv$yhat[test_idx]
      ok <- is.finite(y_t) & is.finite(y_p)
      n_ok <- sum(ok)
      results_eval_metrics_reps <- data.frame(Rep = repp, check.names = FALSE)
      if (n_ok < 2 || isTRUE(stats::sd(y_p[ok]) == 0)) {
        # degenerate predictions (too few or constant ŷ) → deterministic, finite penalties
        for (m in eval_metrics) {
          results_eval_metrics_reps[[m]] <- metric_worst_value(m, y_true = y_t[ok])
        }
      } else {
        # normal path (still guarded against errors inside evaluation_metrics)
        for (m in eval_metrics) {
          results_eval_metrics_reps[[m]] <- safe_metric_value(
            y_true = y_t, y_pred = y_p, metric = m
          )
        }
      }
      # for (m in eval_metrics) {
      #
      #   results_eval_metrics_reps[[m]] <- tryCatch(
      #     evaluation_metrics(y_observed = ypred_cv$y[test_idx],
      #                        y_predicted = ypred_cv$yhat[test_idx],
      #                        eval_metrics = m),
      #
      #     error = function(e) NA_real_
      #   )
      # }
    }

    # ---- K-Folds
    if (is_kfold) {
      group <- test_set_val[[repp]]  # vector length n
      for (j in seq_len(nfolds)) {
        tst <- which(group == j)
        if (!length(tst)) next
        yNA <- y; yNA[tst] <- NA
        ypred_cv$cv_role[tst] <- "test"

        if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            if (model_can == "GBLUP_BRR") {
              mGB <- "BRR"
              additional_params$bayes_model <- mGB
              additional_params$bayes_trait <- trait
              additional_params$ETA        <- model_prep_all_bayes_cv[[mGB]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para  <- model_prep_all_bayes_cv[[mGB]][["bayes_para"]]
            } else {
              additional_params$bayes_model <- model_can
              additional_params$bayes_trait <- trait
              additional_params$ETA        <- model_prep_all_bayes_cv[[model_can]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para  <- model_prep_all_bayes_cv[[model_can]][["bayes_para"]]
            }
            assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = additional_params))
          }, error = function(e) { message(paste("Error in Bayes", model_can, "for", trait, "fold", j, ":", e$message)); handle_error <<- TRUE })
        }

        if (!is.null(engine) && model_can == "GBLUP" && engine == "asreml") {
          tryCatch({
            additional_params$response <- trait
            assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = additional_params))
          }, error = function(e) { message(paste("Error in GBLUP(asreml) for", trait, "fold", j, ":", e$message)); handle_error <<- TRUE })
        }

        if (model_can %in% AI_valid_models) {
          tryCatch({
            assign_preds(tst, predict_with_model(
              model = model_can, y = yNA,
              omics_data = if (exists("omics_data")) omics_data else NULL,
              tst = tst, additional_params = additional_params))
          }, error = function(e) { message(paste("Error in AI model", model_can, "for", trait, "fold", j, ":", e$message)); handle_error <<- TRUE })
        }

        if (isTRUE(verbose)) {
          n_now <- sum(!is.na(ypred_cv$yhat[tst]))
          if (n_now != length(tst)) message(sprintf("[cv] fold %d/%d: predicted %d of %d", j, nfolds, n_now, length(tst)))
        }
      }

      oof <- which(ypred_cv$cv_role == "test" & !is.na(ypred_cv$yhat))
      results_eval_metrics_reps <- data.frame(Rep = repp, check.names = FALSE)
      for (m in eval_metrics) {
        results_eval_metrics_reps[[m]] <- tryCatch(
          # evaluation_metrics(y_observed = ypred_cv$y[oof],
          #                    y_predicted = ypred_cv$yhat[oof],
          #                    eval_metrics = m),
          safe_metric_value(
            y_true = ypred_cv$y[oof],
            y_pred = ypred_cv$yhat[oof],
            metric = m
          ),
          error = function(e) NA_real_
        )
      }
      if (isTRUE(verbose) && length(oof) < sum(ypred_cv$cv_role == "test"))
        message(sprintf("[cv] OOF coverage: %d/%d rows got predictions.",
                        length(oof), sum(ypred_cv$cv_role == "test")))
    }

    ## --- MULTI-ENVIRONMENT CV (CV0/CV1/CV2) -----------------------------------
    ## --- MULTI-ENVIRONMENT CV (CV0/CV1/CV2) -----------------------------------
    if (!is_holdout && !is_kfold) {
      if (is.null(heter_groups)) stop("For CV0/CV1/CV2 you must supply `heter_groups`.", call. = FALSE)
      if (!heter_groups %in% names(ypred_cv)) {
        ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
      }

      group <- test_set_val[[repp]]  # integer vector of length nrow(pheno_data)
      if (length(group) != nrow(pheno_data)) {
        stop(sprintf("Split length mismatch: got %d, expected %d.",
                     length(group), nrow(pheno_data)))
      }

      # unique folds; in CV2, 0 marks train-only singletons
      uf <- sort(unique(group))
      uf <- uf[uf > 0L]   # <- IMPORTANT: skip 0

      if (!length(uf) && isTRUE(verbose)) {
        message("[cv-me] No test folds (>0) produced by the ME scheme; all rows remain train.")
      }

      for (j in uf) {
        tst <- which(group == j)
        if (!length(tst)) next

        # mask current fold as test
        yNA <- y
        yNA[tst] <- NA_real_
        ypred_cv$cv_role[tst] <- "test"

        # ---- route to your models (same as in your k-fold block) ----
        if (model_can %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            if (model_can == "GBLUP_BRR") {
              mGB <- "BRR"
              additional_params$bayes_model <- mGB
              additional_params$bayes_trait <- trait
              additional_params$ETA        <- model_prep_all_bayes_cv[[mGB]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para  <- model_prep_all_bayes_cv[[mGB]][["bayes_para"]]
            } else {
              additional_params$bayes_model <- model_can
              additional_params$bayes_trait <- trait
              additional_params$ETA        <- model_prep_all_bayes_cv[[model_can]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para  <- model_prep_all_bayes_cv[[model_can]][["bayes_para"]]
            }
            assign_preds(tst, predict_with_model("Bayes", y = yNA, tst = tst, additional_params = additional_params))
          }, error = function(e) {
            message(sprintf("Error in Bayes %s for %s (fold %d): %s", model_can, trait, j, e$message))
            handle_error <<- TRUE
          })
        }

        if (!is.null(engine) && model_can == "GBLUP" && engine == "asreml") {
          tryCatch({
            additional_params$response <- trait
            assign_preds(tst, predict_with_model("GBLUP", tst = tst, additional_params = additional_params))
          }, error = function(e) {
            message(sprintf("Error in GBLUP(asreml) for %s (fold %d): %s", trait, j, e$message))
            handle_error <<- TRUE
          })
        }

        if (model_can %in% AI_valid_models) {
          tryCatch({
            assign_preds(
              tst,
              predict_with_model(
                model = model_can, y = yNA,
                omics_data = if (exists("omics_data")) omics_data else NULL,
                tst = tst, additional_params = additional_params
              )
            )
          }, error = function(e) {
            message(sprintf("Error in AI model %s for %s (fold %d): %s", model_can, trait, j, e$message))
            handle_error <<- TRUE
          })
        }

        if (isTRUE(verbose)) {
          n_now <- sum(!is.na(ypred_cv$yhat[tst]))
          if (n_now != length(tst))
            message(sprintf("[cv-me] fold %d: predicted %d of %d", j, n_now, length(tst)))
        }
      }


      # Keep predictions ONLY on test rows (mask any accidental train preds)
      ypred_cv$yhat[ypred_cv$cv_role != "test"] <- NA_real_

      # Sanity: ensure we actually have test rows
      if (!any(ypred_cv$cv_role == "test")) {
        stop("[cv] No 'test' rows were marked. Check that group_vec has valid fold indices.", call. = FALSE)
      }

      # Metrics per environment + overall
      ENV <- sort(unique(ypred_cv[[heter_groups]]))
      by_env <- do.call(
        rbind,
        lapply(ENV, function(env_i) {
          idx <- which(ypred_cv$cv_role == "test" & ypred_cv[[heter_groups]] == env_i)
          row <- list(Rep = repp, Env = env_i)
          for (m in eval_metrics) {
            ok <- is.finite(ypred_cv$y[idx]) & is.finite(ypred_cv$yhat[idx])
            row[[m]] <- if (any(ok)) safe_metric_value(
              y_true = ypred_cv$y[idx][ok], y_pred = ypred_cv$yhat[idx][ok], metric = m
            ) else NA_real_
          }
          as.data.frame(row, check.names = FALSE)
        })
      )
      oof <- which(ypred_cv$cv_role == "test")
      ok  <- is.finite(ypred_cv$y[oof]) & is.finite(ypred_cv$yhat[oof])
      overall <- data.frame(Rep = repp, Env = "ALL", check.names = FALSE)
      for (m in eval_metrics) {
        overall[[m]] <- if (any(ok)) safe_metric_value(ypred_cv$y[oof][ok], ypred_cv$yhat[oof][ok], m) else NA_real_
      }
      results_eval_metrics_reps <- rbind(by_env, overall)

      # # cv_info
      # nfolds_eff <- length(folds)
      # cv_info_list <- list(
      #   method      = cross_validation_meth,
      #   token       = cv_token,
      #   nfolds      = nfolds_eff,
      #   replication = repp,
      #   n_test      = sum(ypred_cv$cv_role == "test"),
      #   n_train     = sum(ypred_cv$cv_role == "train")
      # )

      # Make sure your final return uses cv_info_list
      #cv_info <- cv_info_list
    }
    ## --- END MULTI-ENVIRONMENT CV ---------------------------------------------


    # cv_group <- if (!is_holdout) test_set_val[[repp]] else rep(NA_integer_, nrow(pheno_data))
    # nfolds_eff <- if (!is_holdout) suppressWarnings(max(cv_group, na.rm = TRUE)) else NA_integer_
    # if (!is.finite(nfolds_eff)) nfolds_eff <- NA_integer_

    list(trait = trait,
         rep = rep_i,
         model = model_label,  # show friendly name  # <<< CHANGED
         eval_metrics_reps = results_eval_metrics_reps,
         ypred_cv_Reps_all = ypred_cv,
         cv_info = list(method = cross_validation_meth,
                        token = cv_token,
                        #nfolds = nfolds_eff,
                        nfolds = if (is_kfold) nfolds else NA_integer_,
                        replication = repp,
                        n_test = sum(ypred_cv$cv_role == "test"),
                        n_train = sum(ypred_cv$cv_role == "train")))
  }

  results <- list()

  # 1) user-forced sequential
  if (length(idx_seq_user) > 0) {
    res_seq_user <- lapply(idx_seq_user, function(i) run_one_task(tasks_full[i, ]))
    results <- c(results, res_seq_user)
  }

  # 2) smart policy for the rest
  if (length(idx_rest) > 0) {
    tasks_rest <- tasks_full[idx_rest, , drop = FALSE]

    decision_rest <- sp_decide_policy(
      models = tasks_rest$modell,
      model_params_list = replicate(nrow(tasks_rest), additional_params, simplify = FALSE),
      globals = globals_for_size,
      n_tasks = nrow(tasks_rest),
      user_mode = parallel_mode,
      num_cores = num_cores,
      globals_max_GB = globals_max_GB,
      verbose = isTRUE(verbose),
      sys_name = sys_name,
      prefer_fork = isTRUE(parallel_backend_prefer_fork)
    )

    int_flags <- decision_rest$internal_flags
    idx_seq2  <- which(int_flags)
    idx_par   <- which(!int_flags)

    if (length(idx_seq2) > 0) {
      res_seq2 <- lapply(idx_seq2, function(k) run_one_task(tasks_rest[k, ]))
      results <- c(results, res_seq2)
    }

    # # When you call sp_apply(...), compose both initializers only if I want to keep worker single thread
    # init_all <- compose_initializers(init_py, init_single_thread_blas)
    # res_par <- sp_apply(..., initializer = init_all)

    if (length(idx_par) > 0) {
      res_par <- sp_apply(
        X = idx_par,
        FUN = function(k) run_one_task(tasks_rest[k, ]),
        decision = decision_rest,
        packages = c("dplyr"),
        seed = TRUE,
        initializer = init_py
      )
      results <- c(results, res_par)
    }
  }

  # keep results if any predictions exist (transparent for Hold-Out/K-fold)
  filtered_results <- Filter(function(res) {
    ok <- !is.null(res$eval_metrics_reps) && !is.null(res$ypred_cv_Reps_all)
    if (!ok) return(FALSE)
    df <- res$ypred_cv_Reps_all
    "yhat" %in% names(df) && any(!is.na(df$yhat))
  }, results)

  filtered_results <- Filter(function(res) {
    if (is.null(res$eval_metrics_reps)) return(FALSE)
    mcols <- setdiff(names(res$eval_metrics_reps), "Rep")
    all(is.finite(unlist(res$eval_metrics_reps[mcols], use.names = FALSE)))
  }, results)

  try(future::plan("sequential"), silent = TRUE)
  return(filtered_results)
}




#'
#' resolve_model_name <- function(x, name_lookup) {
#'   if (is.null(x) || length(x) == 0) return(NULL)
#'   res <- name_lookup[x]
#'   res <- res[!is.na(res)]
#'   if (length(res) == 0) return(NULL)
#'   unname(res)
#' }
#'
#' replace_with_canonical <- function(x, canonical_names, friendly_names) {
#'   if (is.null(x) || length(x) == 0) return(NULL)
#'
#'   # Normalize to lowercase
#'   k <- tolower(trimws(x))
#'
#'   # Create lowercase lookup
#'   keys <- tolower(c(canonical_names, friendly_names))
#'   vals <- rep(canonical_names, 2)
#'   lookup <- setNames(vals, keys)
#'
#'   # Lookup canonical names
#'   res <- lookup[k]
#'
#'   # Preserve unknowns as original
#'   res[is.na(res)] <- x[is.na(res)]
#'
#'   unname(res)
#' }
#' #' @export
#' log_thread_env_vars <- function() {
#'   # now log what’s actually set
#'   vars <- c(
#'     "OPENBLAS_NUM_THREADS",
#'     "OMP_NUM_THREADS",
#'     "MKL_NUM_THREADS",
#'     "NUMEXPR_NUM_THREADS",
#'     "TF_INTRA_OP_PARALLELISM_THREADS",
#'     "TF_INTER_OP_PARALLELISM_THREADS",
#'     "TF_ENABLE_ONEDNN_OPTS"
#'   )
#'   for (v in vars) {
#'     message(sprintf("→ %s = %s", v, Sys.getenv(v, unset = "<unset>")))
#'   }
#'
#'   invisible(NULL)
#' }
#'
#'
#' set_parallel_plan <- function(n_trait,
#'                               n_model = 1,
#'                               replication = 1,
#'                               num_cores = NULL,
#'                               globals_max_GB = 4,
#'                               docker_override = FALSE,
#'                               mode = "cross_validation",
#'                               sys_name) {
#'   task_type <- if (mode == "cross_validation") "Cross Validation" else "True Prediction"
#'   message(sprintf(
#'     "\nSetting up parallel plan for %s\n",
#'     task_type
#'   ))
#'   sys_name  <- if (docker_override) "Windows" else Sys.info()[["sysname"]]
#'   plan_type <- if (sys_name == "Windows") "multisession" else "multicore"
#'   message(sprintf("→ OS: %s (docker_override=%s) → using future plan '%s'",
#'                   sys_name, docker_override, plan_type))
#'
#'   phys <- parallel::detectCores(logical = FALSE)
#'   avail <- if (is.null(num_cores)) floor(phys * 0.5) else floor(num_cores)
#'   message(sprintf("→ Core budget: avail = %d (num_cores=%s)", avail,
#'                   if (is.null(num_cores)) "auto" else num_cores))
#'
#'   max_needed <- n_trait * n_model * replication
#'   workers <- max(1L, min(avail, phys, n_trait * n_model * replication))
#'   message(sprintf("→ Workload: traits × models × reps = %d; spawning %d worker(s)",
#'                   max_needed, workers))
#'
#'   options(future.globals.maxSize = globals_max_GB * 1024^3L,
#'           future.rng.onMisuse    = "ignore")
#'   message(sprintf("→ future.globals.maxSize = %.1f GB", globals_max_GB))
#'
#'   if (workers > 1L) {
#'     future::plan(plan_type, workers = workers, gc = TRUE)
#'     message(sprintf("→ future plan set to '%s' with %d workers", plan_type, workers))
#'   } else {
#'     future::plan("sequential")
#'     message("→ Using sequential plan (workers = 1)")
#'   }
#'
#'   # now log all your thread‑control env vars in one shot
#'   #log_thread_env_vars()
#'
#'   invisible(workers)
#' }
#'
#' #' Compute per-future thread budget and apply it globally
#' #'
#' #' @return integer number of threads to use per model
#' set_per_worker_threads <- function() {
#'   # how many futures will be running in parallel?
#'   # workers <- future::plan()$workers
#'   # if (is.null(workers) || workers < 1L) workers <- 1L
#'   #
#'   # cores   <- parallel::detectCores(logical = TRUE)
#'   workers <- future::nbrOfWorkers()
#'   # how many hardware threads on the box?
#'   cores   <- parallel::detectCores(logical = TRUE)
#'   # give each future an equal share
#'   intra   <- max(1L, floor(cores / workers))
#'
#'   # enforce in OpenMP / MKL / BLAS
#'   ## These also good whene using Eigen, TBB, etc.
#'   Sys.setenv(OMP_NUM_THREADS = intra)
#'   Sys.setenv(MKL_NUM_THREADS = intra)
#'
#'
#'   invisible(intra)
#'
#' }
#'
#' #' Title
#' #'
#' #' @param model
#' #' @param y
#' #' @param omics_data
#' #' @param tst
#' #' @param additional_params
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' predict_with_model <- function(model = NULL,
#'                                y = NULL,
#'                                omics_data = NULL,
#'                                tst = NULL,
#'                                additional_params = NULL) {
#'   # Generalized function to handle predictions for various models
#'
#'   switch(model,
#'          "Xgboost" = AI_xgboost_cv(y = y, omics = omics_data, tst = tst, eta = additional_params$eta,
#'                                    nrounds = additional_params$nrounds, max_depth = additional_params$max_depth,
#'                                    scaling = additional_params$scaling, omic_count = additional_params$omic_count,
#'                                    centering = additional_params$centering, xgb_gamma = additional_params$xgb_gamma,
#'                                    colsample_bytree = additional_params$colsample_bytree,
#'                                    subsample = additional_params$subsample, min_child_weight = additional_params$min_child_weight,
#'                                    xgb_alpha = additional_params$xgb_alpha, xgb_lambda = additional_params$xgb_lambda,
#'                                    xgb_booster = additional_params$xgb_booster,
#'                                    xgb_rate_drop = additional_params$xgb_rate_drop,xgb_skip_drop = additional_params$xgb_skip_drop,
#'                                    xgb_objective = additional_params$xgb_objective,xgb_sample_type = additional_params$xgb_sample_type,
#'                                    xgb_normalize_type = additional_params$xgb_normalize_type,
#'                                    early_stop_for_iteration_xgb = additional_params$early_stop_for_iteration_xgb),
#'          #unified deep learning branch
#'          "cnn" = ,
#'          "resnet" = ,
#'          "ft_transformer" = ,
#'          "saint" = ,
#'          "tabnet" = ,
#'          "node" = ,
#'          "deepfm" = ,
#'          "dcnv2" = ,
#'          "nam" = ,
#'          "moe" = ,
#'          "gp_dkl" = ,
#'          "mlp_with_attention" = ,
#'          "mlp" = deep_learning_model(
#'            y = y,
#'            omics = omics_data,
#'            tst = tst,
#'            scaling = additional_params$scaling,
#'            centering = additional_params$centering,
#'            crossval = additional_params$crossval,
#'            omic_count = additional_params$omic_count,
#'            early_stop = additional_params$early_stop,
#'            deep_learning_model = additional_params$deep_learning_model,
#'
#'            optimizer_name = additional_params$optimizer_name,
#'            use_amp = additional_params$use_amp,
#'            max_grad_norm = additional_params$max_grad_norm,
#'            auto_class_weights = additional_params$auto_class_weights,
#'
#'            ##### CNN
#'            cnn_neurons_per_layer = additional_params$cnn_neurons_per_layer,
#'            cnn_kernel_size = additional_params$cnn_kernel_size,
#'            cnn_dense_layers = additional_params$cnn_dense_layers,
#'            cnn_use_max_pool = additional_params$cnn_use_max_pool,
#'            cnn_pool_kernel = additional_params$cnn_pool_kernel,
#'            cnn_pool_stride = additional_params$cnn_pool_stride,
#'            cnn_pool_padding = additional_params$cnn_pool_padding,
#'            cnn_learning_rate = additional_params$cnn_learning_rate,
#'            cnn_separable = additional_params$cnn_separable,
#'            cnn_dilations = additional_params$cnn_dilations,
#'            cnn_use_se = additional_params$cnn_use_se,
#'            cnn_norm_type = additional_params$cnn_norm_type,
#'            cnn_pool_type = additional_params$cnn_pool_type,
#'            cnn_use_global_pool = additional_params$cnn_use_global_pool,
#'
#'            ##### ResNet
#'            resnet_neurons_per_block = additional_params$resnet_neurons_per_block,
#'            resnet_blocks = additional_params$resnet_blocks,
#'            resnet_learning_rate = additional_params$resnet_learning_rate,
#'
#'            ##### FT Transformer
#'            ft_d_model = additional_params$ft_d_model,
#'            ft_heads = additional_params$ft_heads,
#'            ft_layers = additional_params$ft_layers,
#'            ft_ff_mult = additional_params$ft_ff_mult,
#'            ft_dropout = additional_params$ft_dropout,
#'            ft_token_dropout = additional_params$ft_token_dropout,
#'            ft_use_cls = additional_params$ft_use_cls,
#'
#'            ##### SAINT
#'            saint_d_model = additional_params$saint_d_model,
#'            saint_heads = additional_params$saint_heads,
#'            saint_layers = additional_params$saint_layers,
#'            saint_ff_mult = additional_params$saint_ff_mult,
#'            saint_dropout = additional_params$saint_dropout,
#'            saint_token_dropout = additional_params$saint_token_dropout,
#'            saint_use_cls = additional_params$saint_use_cls,
#'
#'            ##### Grouping controls
#'            use_grouping = additional_params$use_grouping,
#'            group_trigger = additional_params$group_trigger,
#'            group_method = additional_params$group_method,
#'            init_group_size = additional_params$init_group_size,
#'            max_tokens = additional_params$max_tokens,
#'            kmeans_batch = additional_params$kmeans_batch,
#'            kmeans_iter = additional_params$kmeans_iter,
#'
#'            ##### TabNet
#'            tabnet_steps = additional_params$tabnet_steps,
#'            tabnet_feature_dim = additional_params$tabnet_feature_dim,
#'            tabnet_output_dim = additional_params$tabnet_output_dim,
#'            tabnet_gamma = additional_params$tabnet_gamma,
#'            tabnet_lambda_sparse = additional_params$tabnet_lambda_sparse,
#'
#'            ##### NODE
#'            node_trees = additional_params$node_trees,
#'            node_depth = additional_params$node_depth,
#'
#'            ##### DeepFM
#'            deepfm_k = additional_params$deepfm_k,
#'            deepfm_hidden = additional_params$deepfm_hidden,
#'
#'            ##### DCNv2
#'            dcn_layers = additional_params$dcn_layers,
#'            dcn_hidden = additional_params$dcn_hidden,
#'
#'            ##### NAM
#'            nam_hidden = additional_params$nam_hidden,
#'            nam_activation = additional_params$nam_activation,
#'            nam_add_linear = additional_params$nam_add_linear,
#'            nam_l1 = additional_params$nam_l1,
#'
#'            ##### Mixture-of-Experts
#'            moe_n_experts = additional_params$moe_n_experts,
#'            moe_expert_hidden = additional_params$moe_expert_hidden,
#'            moe_gate_hidden = additional_params$moe_gate_hidden,
#'            moe_temperature = additional_params$moe_temperature,
#'            moe_sparse_topk = additional_params$moe_sparse_topk,
#'            moe_entropy_reg = additional_params$moe_entropy_reg,
#'
#'            ##### GP / RFF
#'            gp_use_variational = additional_params$gp_use_variational,
#'            gp_num_inducing = additional_params$gp_num_inducing,
#'            gp_feature_dim = additional_params$gp_feature_dim,
#'            gp_kernel = additional_params$gp_kernel,
#'            gp_ard = additional_params$gp_ard,
#'            gp_lr_mult = additional_params$gp_lr_mult,
#'            rff_features = additional_params$rff_features,
#'            rff_lengthscale = additional_params$rff_lengthscale,
#'            rff_deep_hidden = additional_params$rff_deep_hidden,
#'
#'            ##### General deep learning params
#'            model_type = additional_params$model_type,
#'            epochs = additional_params$epochs,
#'            batch_size = additional_params$batch_size,
#'            dropout = additional_params$dropout,
#'            l2_weight_decay = additional_params$l2_weight_decay,
#'            l2_regularizer_dp = additional_params$l2_regularizer_dp,
#'            dropout_rate = additional_params$dropout_rate,
#'            batch_norm = additional_params$batch_norm,
#'            validation_split = additional_params$validation_split,
#'            compile_model = additional_params$compile_model,
#'            deterministic = additional_params$deterministic,
#'            random_seed = additional_params$random_seed,
#'            device = additional_params$device,
#'
#'            ##### MLP and Attention
#'            mlp_neurons_per_layer = additional_params$mlp_neurons_per_layer,
#'            mlp_learning_rate = additional_params$mlp_learning_rate,
#'            final_attention = additional_params$final_attention,
#'            attention_across_multiple_layers = additional_params$attention_across_multiple_layers,
#'            heteroscedastic = additional_params$heteroscedastic
#'          ),
#'          "RandomForest" = AI_randomforest_cv(y = y, omics = omics_data, tst = tst,
#'                                              scaling = additional_params$scaling,
#'                                              centering = additional_params$centering, ntree = additional_params$ntree,
#'                                              omic_count = additional_params$omic_count),
#'          "PartialLeastSquare" = AI_pls_cv(y = y, omics = omics_data, tst = tst,
#'                                           scaling = additional_params$scaling,
#'                                           centering = additional_params$centering, ncomp = additional_params$ncomp,
#'                                           omic_count = additional_params$omic_count),
#'          "Ridge_Regression" = AI_ridge_regression_cv(y = y, omics = omics_data, tst = tst,
#'                                                      scaling = additional_params$scaling,
#'                                                      centering = additional_params$centering,
#'                                                      omic_count = additional_params$omic_count),
#'          "Lasso" = AI_lasso_cv(y = y, omics = omics_data, tst = tst, scaling = additional_params$scaling,
#'                                centering = additional_params$centering,
#'                                omic_count = additional_params$omic_count),
#'          "SupportVectorMachine" = AI_svm_cv(y = y, omics = omics_data, tst = tst,
#'                                             scaling = additional_params$scaling,
#'                                             centering = additional_params$centering,
#'                                             C_value = additional_params$C_value,
#'                                             degree_value = additional_params$degree_value,
#'                                             scale_value = additional_params$scale_value,
#'                                             offset_value = additional_params$offset_value,
#'                                             omic_count = additional_params$omic_count),
#'          "K-NearestNeighbors" = AI_knn_cv(y = y, omics = omics_data, tst = tst,
#'                                           scaling = additional_params$scaling,
#'                                           centering = additional_params$centering, k = additional_params$k,
#'                                           omic_count = additional_params$omic_count),
#'          "GBLUP" = asreml_mod_cv(asreml_models_prep_cv = additional_params$asreml_models_prep_cv,
#'                                  pheno_data = additional_params$pheno_data,
#'                                  response = additional_params$response,
#'                                  heter_groups = additional_params$heter_groups,
#'                                  gen_name = additional_params$gen_name, tst = tst),
#'          "Bayes" = bayes_mod_cv(y = y, ETA = additional_params$ETA, weights = additional_params$weights,
#'                                 bayes_para = additional_params$bayes_para, tst = tst,
#'                                 bayes_model = additional_params$bayes_model, bayes_trait = additional_params$bayes_trait)
#'   )
#' }
#'
#' ###########################
#'
#' #' Title
#' #'
#' #' @param pheno_data
#' #' @param test_set
#' #' @param response
#' #' @param gen_name
#' #' @param test_size
#' #' @param random_state
#' #' @param replication
#' #' @param weights
#' #' @param model_prep_all_bayes_cv
#' #' @param asreml_models_prep_cv
#' #' @param ml_dat_res
#' #' @param heter_groups
#' #' @param verbose
#' #' @param num_cores
#' #' @param nfolds
#' #' @param cross_validation_meth
#' #' @param sampling_method
#' #' @param eval_metrics
#' #' @param GS_model_cv
#' #' @param scaling
#' #' @param centering
#' #' @param eta
#' #' @param nrounds
#' #' @param max_depth
#' #' @param gamma
#' #' @param subsample
#' #' @param colsample_bytree
#' #' @param ncomp
#' #' @param ntree
#' #' @param k
#' #' @param c
#' #' @param ...
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' models_execute_crossval <- function(pheno_data = NULL,
#'                                     test_set = NULL,
#'                                     response = NULL,
#'                                     gen_name = NULL,
#'                                     test_size = NULL,
#'                                     random_state = NULL,
#'                                     replication = NULL,
#'                                     weights = NULL,
#'                                     selected_raw, # nd_mods
#'                                     gam_method = NULL, # nd_mods
#'                                     max_features = 50,
#'                                     k_value = 5,
#'                                     var_explained = 0.9,
#'                                     engine = NULL,
#'                                     model_prep_all_bayes_cv = NULL,
#'                                     asreml_models_prep_cv = NULL,
#'                                     ml_dat_res = NULL,
#'                                     heter_groups = NULL,
#'                                     verbose = FALSE,
#'                                     num_cores = NULL,
#'                                     nfolds = 5,
#'                                     cross_validation_meth = NULL, ## this handle th ETA for bayes model
#'                                     sampling_method = NULL,
#'                                     eval_metrics = NULL,
#'                                     bayes_model = NULL,
#'                                     GS_model_cv = NULL,
#'                                     scaling = FALSE,
#'                                     centering = TRUE,
#'                                     eta = 0.1, ## xgboost
#'                                     nrounds = 100, ## xgboost
#'                                     max_depth = 6, ## xgboost
#'                                     xgb_gamma = 4, ## xgboost
#'                                     subsample = 0.5, ## xgboost
#'                                     colsample_bytree = 1, ## xgboost
#'                                     xgb_alpha = 0.001, ## xgboost linear
#'                                     xgb_lambda = 1, ## xgboost linear
#'                                     min_child_weight = 1, ## xgboost
#'                                     early_stop_for_iteration_xgb = TRUE,
#'                                     xgb_booster = "dart", #"gbtree",
#'                                     xgb_rate_drop = 0.1,
#'                                     xgb_skip_drop = 0.5,
#'                                     xgb_objective = "reg:squarederror",
#'                                     xgb_sample_type = "uniform",
#'                                     xgb_normalize_type = "tree",
#'                                     ncomp = 3, #### pls
#'                                     ntree = 500, ### random forest
#'                                     k = 5, ## for knn
#'                                     svm_kernel = "Gaussian", # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
#'                                     sigma_value  = 0.1,       # Default sigma value for RBF kernel
#'                                     C_value  = 1,             # Default cost parameter
#'                                     degree_value = 3,        # Default degree for polynomial kernel
#'                                     scale_value  = 1,         # Default scale for polynomial kernel
#'                                     offset_value = 1,
#'                                     num_hidden_layers = 1,
#'                                     neurons_per_layer = 64,
#'                                     learning_rate_dp = 0.001,
#'
#'                                     early_stop = TRUE,
#'
#'                                     crossval = TRUE,
#'                                     optimizer_name = "adam",
#'                                     use_amp        = TRUE,
#'                                     max_grad_norm  = 1.0,
#'                                     auto_class_weights = FALSE,
#'                                     ##### cnn
#'                                     cnn_neurons_per_layer = as.integer(c(64, 64, 64)),
#'                                     cnn_kernel_size = 3L,
#'                                     cnn_dense_layers = as.integer(c(256, 128, 64)),
#'                                     cnn_use_max_pool = FALSE,
#'                                     cnn_pool_kernel = 2L,
#'                                     cnn_pool_stride = 2L,
#'                                     cnn_pool_padding = 0L,
#'                                     cnn_learning_rate = 1e-3,
#'                                     cnn_separable=TRUE,
#'                                     cnn_dilations=c(1,2,4),
#'                                     cnn_use_se=TRUE,
#'                                     cnn_norm_type="group",
#'                                     cnn_pool_type="conv",
#'                                     cnn_use_global_pool=FALSE,
#'                                     ##### resnet
#'                                     resnet_neurons_per_block = as.integer(c(256, 128, 64)),
#'                                     resnet_blocks = 3,
#'                                     resnet_learning_rate = 1e-3,
#'                                     #### ft_transformer
#'                                     ft_d_model = 192L,
#'                                     ft_heads = 8L,
#'                                     ft_layers = 3L,
#'                                     ft_ff_mult = 4L,
#'                                     ft_dropout = 0.1,
#'                                     ft_token_dropout = 0.0,
#'                                     ft_use_cls = TRUE,
#'                                     #### saint
#'                                     saint_d_model = 128L,
#'                                     saint_heads = 8L,
#'                                     saint_layers = 3L,
#'                                     saint_ff_mult = 4L,
#'                                     saint_dropout = 0.1,
#'                                     saint_token_dropout = 0.0,
#'                                     saint_use_cls       = TRUE,
#'                                     ###### Grouping controls (FT/SAINT and NAM/MoE)
#'                                     use_grouping   = FALSE,
#'                                     group_trigger  = 2048,
#'                                     group_method   = "auto",
#'                                     init_group_size = 64,
#'                                     max_tokens      = 1024,
#'                                     kmeans_batch    = 4096,
#'                                     kmeans_iter     = 100,
#'                                     #### tabnet
#'                                     tabnet_steps = 5L,
#'                                     tabnet_feature_dim = 64L,
#'                                     tabnet_output_dim = 64L,
#'                                     tabnet_gamma = 1.5,
#'                                     tabnet_lambda_sparse = 1e-4,
#'                                     #### node
#'                                     node_trees = 8L,
#'                                     node_depth = 3L,
#'                                     #### deepfm
#'                                     deepfm_k = 16L,
#'                                     deepfm_hidden = as.integer(c(128, 64)),
#'                                     #### dcnv2
#'                                     dcn_layers = 3L,
#'                                     dcn_hidden = as.integer(c(256, 128, 64)),
#'                                     #### nam
#'                                     nam_hidden = as.integer(c(32, 16)),
#'                                     nam_activation = "relu",
#'                                     nam_add_linear = TRUE,
#'                                     nam_l1 = 1e-4,
#'                                     ### moe
#'                                     moe_n_experts = 4L,
#'                                     moe_expert_hidden = as.integer(c(128, 64)),
#'                                     moe_gate_hidden = 128L,
#'                                     moe_temperature = 1.0,
#'                                     moe_sparse_topk = NA,
#'                                     moe_entropy_reg = 0.0,
#'                                     ### gp_dkl/RFF knobs
#'                                     gp_use_variational = TRUE,
#'                                     gp_num_inducing = 256L,
#'                                     gp_feature_dim = 64L,
#'                                     gp_kernel = "rbf",
#'                                     gp_ard = TRUE,
#'                                     gp_lr_mult = 0.5,
#'                                     rff_features       = 1024,
#'                                     rff_lengthscale    = 1.0,
#'                                     rff_deep_hidden    = c(128),
#'                                     ##### General dp
#'                                     model_type = "resnet",
#'                                     epochs = 10,
#'                                     batch_size = 64 ,
#'                                     dropout = 0.2,
#'                                     l2_weight_decay = 1e-4,
#'                                     l2_regularizer_dp = 0.001,
#'                                     dropout_rate = 0.5,
#'                                     batch_norm = TRUE,
#'                                     validation_split = 0.2,
#'                                     compile_model = FALSE,
#'                                     deterministic = TRUE,
#'                                     random_seed = 123,
#'                                     device = NULL,
#'                                     #### mlp and attention
#'                                     mlp_neurons_per_layer = as.integer(c(128, 64)),
#'                                     mlp_learning_rate = 1e-3,
#'                                     final_attention = TRUE,
#'                                     attention_across_multiple_layers = TRUE,
#'                                     heteroscedastic = TRUE,
#'                                     docker_nd_usage = FALSE,
#'                                     globals_max_GB = 4,
#'                                     ...){
#'
#'  #browser()
#'   on.exit(future::plan("sequential"), add = TRUE)
#'   msg <- "\n==================================================\n"
#'
#'   if(!is.null(cross_validation_meth) & length(cross_validation_meth)>1){
#'     stop(paste(msg, 'use only one cross_validation method at a time.'), call. = FALSE)
#'   }
#'   ##### if user choose repeated and stratified CV strategy but
#'   ## forget to choose sampling stratgy or replication is not defined.
#'   patterns <- c("stratified", "Repeated")
#'
#'   # Use sapply to apply grep to each pattern and return a named logical vector indicating presence.
#'   if(is.null(sampling_method) | is.null(replication)){
#'     matche_strings <- sapply(patterns, function(pattern) {
#'       length(grep(pattern, cross_validation_meth, ignore.case = TRUE)) > 0
#'     }, simplify = FALSE)
#'
#'     # Name the list elements with the patterns.
#'     names(matche_strings) <- patterns
#'
#'     # Filter to get only the patterns that were found.
#'     present_patterns <- names(matche_strings)[unlist(matche_strings)]
#'
#'     if("stratified"%in%present_patterns) sampling_method <- "stratified"
#'
#'
#'   }
#'
#'   # Prepare additional parameters for model prediction
#'   ETA <-  NULL
#'   bayes_model <- NULL
#'   bayes_trait <- NULL
#'   bayes_para <- NULL
#'
#'   results_seq_dl <- NULL
#'
#'   if(!is.null(ml_dat_res)){
#'   omics_data <-  ml_dat_res[["merged_data"]][["merge_data"]]
#'
#'   if(!"merged_data_test"%in%names(ml_dat_res)){
#'   if(!is.null(test_set)) {
#'
#'     omics_data <-  omics_data[!rownames(omics_data) %in% test_set, ]
#'
#'   }
#'
#'   }
#'
#'
#'   if("omic_count"%in%names(ml_dat_res)){
#'     omic_count <- ml_dat_res[["omic_count"]]
#'
#'   }else{
#'     omic_count <- NULL
#'   }
#'
#'   #print(omics_data)
#'   }
#'
#'
#'   additional_params <- list(ETA = ETA, weights = weights, bayes_para = bayes_para,
#'                             bayes_model = bayes_model, bayes_trait = bayes_trait,
#'                             scaling = scaling,centering = centering,
#'                             eta = eta, nrounds = nrounds,
#'                             max_depth = max_depth, xgb_gamma = xgb_gamma,
#'                             early_stop_for_iteration_xgb = early_stop_for_iteration_xgb,
#'                             colsample_bytree = colsample_bytree, subsample = subsample, ntree = ntree,
#'                             xgb_alpha = xgb_alpha, xgb_lambda = xgb_lambda,
#'                             min_child_weight = min_child_weight,
#'                             xgb_booster = xgb_booster, xgb_rate_drop = xgb_rate_drop,
#'                             xgb_skip_drop = xgb_skip_drop,
#'                             xgb_objective = xgb_objective,
#'                             xgb_sample_type = xgb_sample_type,
#'                             xgb_normalize_type = xgb_normalize_type,
#'                             ncomp = ncomp, C_value = C_value,degree_value = degree_value,
#'                             scale_value = scale_value, offset_value = offset_value,
#'                             k = k, omic_count = omic_count,
#'                             asreml_models_prep_cv = asreml_models_prep_cv, ## asreml
#'                             gen_name = gen_name,
#'                             pheno_data = pheno_data,
#'
#'                             max_features = max_features,
#'                             heter_groups = heter_groups,
#'
#'                             crossval = crossval,
#'                             early_stop = early_stop,
#'                             model_type = model_type,
#'                             #####
#'                             optimizer_name = optimizer_name,
#'                             use_amp        = use_amp,
#'                             max_grad_norm  = max_grad_norm,
#'                             auto_class_weights = auto_class_weights,
#'                             cnn_neurons_per_layer = cnn_neurons_per_layer,
#'                             cnn_kernel_size = cnn_kernel_size,
#'                             cnn_dense_layers = cnn_dense_layers,
#'                             cnn_use_max_pool = cnn_use_max_pool,
#'                             cnn_pool_kernel = cnn_pool_kernel,
#'                             cnn_pool_stride = cnn_pool_stride,
#'                             cnn_pool_padding = cnn_pool_padding,
#'                             cnn_learning_rate = cnn_learning_rate,
#'                             cnn_separable = cnn_separable,
#'                             cnn_dilations = cnn_dilations,
#'                             cnn_use_se = cnn_use_se,
#'                             cnn_norm_type = cnn_norm_type,
#'                             cnn_pool_type = cnn_pool_type,
#'                             cnn_use_global_pool = cnn_use_global_pool,
#'                             ##### resnet
#'                             resnet_neurons_per_block = resnet_neurons_per_block,
#'                             resnet_blocks = resnet_blocks,
#'                             resnet_learning_rate = resnet_learning_rate,
#'                             #### ft_transformer
#'                             ft_d_model = ft_d_model,
#'                             ft_heads = ft_heads,
#'                             ft_layers = ft_layers,
#'                             ft_ff_mult = ft_ff_mult,
#'                             ft_dropout = ft_dropout,
#'                             ft_token_dropout = ft_token_dropout,
#'                             ft_use_cls = ft_use_cls,
#'                             #### saint
#'                             saint_d_model = saint_d_model,
#'                             saint_heads = saint_heads,
#'                             saint_layers = saint_layers,
#'                             saint_ff_mult = saint_ff_mult,
#'                             saint_dropout = saint_dropout,
#'                             saint_token_dropout = saint_token_dropout,
#'                             saint_use_cls = saint_use_cls,
#'                             ###### Grouping controls (FT/SAINT and NAM/MoE)
#'                             use_grouping = use_grouping,
#'                             group_trigger = group_trigger,
#'                             group_method = group_method,
#'                             init_group_size = init_group_size,
#'                             max_tokens = max_tokens,
#'                             kmeans_batch = kmeans_batch,
#'                             kmeans_iter = kmeans_iter,
#'                             #### tabnet
#'                             tabnet_steps = tabnet_steps,
#'                             tabnet_feature_dim = tabnet_feature_dim,
#'                             tabnet_output_dim = tabnet_output_dim,
#'                             tabnet_gamma = tabnet_gamma,
#'                             tabnet_lambda_sparse = tabnet_lambda_sparse,
#'                             #### node
#'                             node_trees = node_trees,
#'                             node_depth = node_depth,
#'                             #### deepfm
#'                             deepfm_k = deepfm_k,
#'                             deepfm_hidden = deepfm_hidden,
#'                             #### dcnv2
#'                             dcn_layers = dcn_layers,
#'                             dcn_hidden = dcn_hidden,
#'                             #### nam
#'                             nam_hidden = nam_hidden,
#'                             nam_activation = nam_activation,
#'                             nam_add_linear = nam_add_linear,
#'                             nam_l1 = nam_l1,
#'                             ### moe
#'                             moe_n_experts = moe_n_experts,
#'                             moe_expert_hidden = moe_expert_hidden,
#'                             moe_gate_hidden = moe_gate_hidden,
#'                             moe_temperature = moe_temperature,
#'                             moe_sparse_topk = moe_sparse_topk,
#'                             moe_entropy_reg = moe_entropy_reg,
#'                             ### gp_dkl/RFF knobs
#'                             gp_use_variational = gp_use_variational,
#'                             gp_num_inducing = gp_num_inducing,
#'                             gp_feature_dim = gp_feature_dim,
#'                             gp_kernel = gp_kernel,
#'                             gp_ard = gp_ard,
#'                             gp_lr_mult = gp_lr_mult,
#'                             rff_features = rff_features,
#'                             rff_lengthscale = rff_lengthscale,
#'                             rff_deep_hidden = rff_deep_hidden,
#'                             ##### General dp
#'                             model_type = model_type,
#'                             epochs = epochs,
#'                             batch_size = batch_size,
#'                             dropout = dropout,
#'                             l2_weight_decay = l2_weight_decay,
#'                             l2_regularizer_dp = l2_regularizer_dp,
#'                             dropout_rate = dropout_rate,
#'                             batch_norm = batch_norm,
#'                             validation_split = validation_split,
#'                             compile_model = compile_model,
#'                             deterministic = deterministic,
#'                             random_seed = random_seed,
#'                             device = device,
#'                             #### mlp and attention
#'                             mlp_neurons_per_layer = mlp_neurons_per_layer,
#'                             mlp_learning_rate = mlp_learning_rate,
#'                             final_attention = final_attention,
#'                             attention_across_multiple_layers = attention_across_multiple_layers,
#'                             heteroscedastic = heteroscedastic
#'
#'                             )
#'
#'
#'
#'   AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare",
#'                        "SupportVectorMachine", "K-NearestNeighbors", "Lasso",
#'                        "Ridge_Regression")
#'
#'
#'   bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
#'   bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
#'
#'
#'   # Canonical names your
#'   canonical_names <- c(
#'     "cnn", "ft_transformer", "saint", "tabnet", "node",
#'     "deepfm", "dcnv2", "nam", "moe", "gp_dkl",
#'     "mlp_with_attention", "mlp", "resnet"
#'   )
#'
#'   # Friendly names for users
#'   friendly_names <- c(
#'     "Conv1DNet", "TabTransformer", "TabAttention", "TabNet", "LightTreeNet",
#'     "FactorNet", "CrossNet", "NeuralAdditive", "MixtureOfExperts", "GPNet",
#'     "DenseAttentionNet", "DenseNeuralNet", "ResNet"
#'   )
#'
#'   # Build lookup table, both canonical and friendly point to canonical
#'   name_lookup <- setNames(
#'     rep(canonical_names, 2),
#'     c(canonical_names, friendly_names)
#'   )
#'
#'
#'   dp_models <-  resolve_model_name(GS_model_cv, name_lookup = name_lookup)
#'   if(!is.null(dp_models) & length(dp_models)>0){
#'     AI_valid_models <- c(AI_valid_models, dp_models)
#'
#'     GS_model_cv <-  replace_with_canonical(GS_model_cv, canonical_names = canonical_names,
#'                                            friendly_names = friendly_names)
#'     model_type <- replace_with_canonical(model_type, canonical_names = canonical_names,
#'                                          friendly_names = friendly_names)
#'   }
#'
#'   # create mapping
#'
#'   asreml_model <- "GBLUP"
#'
#'   n_trait <- length(response)
#'   n_model <- length(GS_model_cv)
#'
#'
#'
#'   holds_out_methods_avail <- c("Hold_Out",
#'                                "Stratified_Hold_Out",
#'                                "Repeated_Hold_Out",
#'                                "Repeated_Stratified_Hold_Out"
#'                                #"Leave_one_Out"
#'                                )
#'
#'   Kfolds_methods_avail <- c("K-Folds",
#'                             "Stratified_K-Folds",
#'                             "Repeated_K-Folds",
#'                             "Repeated_Stratified_K-Folds")
#'
#'
#'   # Create a named vector to map hold-out methods to K-folds methods
#'   method_mapping <- setNames(Kfolds_methods_avail, holds_out_methods_avail)
#'
#'   # Function to convert hold-out method to K-folds method
#'   convert_to_kfolds <- function(method) {
#'     if (method %in% names(method_mapping)) {
#'       return(method_mapping[method])
#'     } else {
#'       stop(paste(msg, "Provided method is not available in hold-out methods."), call. = FALSE)
#'     }
#'   }
#'
#'   if(cross_validation_meth %in% holds_out_methods_avail){
#'
#'     cross_validation_meth <- as.character(convert_to_kfolds(cross_validation_meth))
#'
#'     if(is.null(nfolds)) nfolds <- 5
#'
#'   }
#'
#'
#'   CVs_multi_envs_methods_avail <- c("CV1",
#'                                     "CV2",
#'                                     "Repeated_CV1",
#'                                     "Repeated_CV2")
#'
#'
#'   #ypred_cv_Reps_all <- list()
#'   ###############################################################
#'   # Main logic
#'   sys_name <- Sys.info()["sysname"]
#'   if(docker_nd_usage) sys_name <- "Windows"
#'   workers <- set_parallel_plan(n_trait = n_trait, n_model = n_model,
#'                     replication = replication,num_cores = num_cores,
#'                     sys_name = sys_name,
#'                     globals_max_GB = globals_max_GB,
#'                     docker_override = docker_nd_usage)
#'
#'
#'   # Create a list of all combinations of response variables and replications
#'   tasks_full <- expand.grid(response = response,
#'                        replication = seq_len(replication),
#'                        modell = GS_model_cv,
#'                        stringsAsFactors = FALSE)
#'
#'   tasks <- subset(tasks_full, ! modell %in% dp_models)
#'
#'   chunk_size <- if (workers > 1) ceiling(nrow(tasks) / workers) else NULL
#'
#'   results <- future.apply::future_lapply(seq_len(nrow(tasks)),
#'                                          future.packages   = c("dplyr"),
#'                                          future.seed       = TRUE,
#'                                          future.chunk.size = chunk_size,
#'                                          function(i) {
#'     task_row <- tasks[i, ]
#'
#'     trait <- as.character(task_row$response)
#'     rep <- as.integer(task_row$replication)
#'     model <- as.character(task_row$modell)
#'
#'
#'     # if(any(model%in%dp_models)){
#'     #  additional_params$deep_learning_model <- model
#'     #  #model_use <- model
#'     #  model <- "deep_learning_model"
#'     # }
#'
#'     repp <- 1
#'
#'     if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1) {
#'       base_seed <- as.integer(random_state)
#'       new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
#'     } else {
#'       base_seed <- 123L
#'       new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
#'     }
#'
#'     if (cross_validation_meth %in% c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail)) {
#'       if (cross_validation_meth %in% holds_out_methods_avail) {
#'         test_set_val <- hold_out_stratified_and_un(pheno_data = pheno_data,
#'                                                    gen_name = gen_name,
#'                                                    response = trait,
#'                                                    test_size = test_size,
#'                                                    random_state = new_seed,
#'                                                    replication = repp,
#'                                                    sampling_method = sampling_method)
#'       } else if (cross_validation_meth %in% Kfolds_methods_avail) {
#'         test_set_val <- kfolds_stratified_un(pheno_data = pheno_data,
#'                                              gen_name = gen_name,
#'                                              response = trait,
#'                                              #test_size = test_size,
#'                                              nfolds = nfolds,
#'                                              random_state = new_seed,
#'                                              replication = repp,
#'                                              sampling_method = sampling_method)
#'       } else if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'         CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
#'         test_set_val <- CV0_CV1_CV2_for_multi_environment(pheno_data = pheno_data,
#'                                                       gen_name = gen_name,
#'                                                       response = trait,
#'                                                       #test_size = test_size,
#'                                                       CV = CV,
#'                                                       nfolds = nfolds,
#'                                                       heter_groups = heter_groups,
#'                                                       random_state = new_seed,
#'                                                       replication = repp
#'                                                       #sampling_method = sampling_method
#'                                                       )
#'       }
#'     } else {
#'       stop(paste(msg, "Unsupported cross-validation method specified. Choose from: ",
#'                  paste(c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail), collapse = ", ")), call. = FALSE)
#'     }
#'
#'     len_y <- nrow(pheno_data)
#'     y <- as.double(pheno_data[[trait]])
#'
#'     if (cross_validation_meth %in% c(Kfolds_methods_avail, holds_out_methods_avail)) {
#'       ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
#'       colnames(ypred_cv) <- c("y", "yhat")
#'       ypred_cv <- as.data.frame(ypred_cv)
#'       ypred_cv[, "y"] <- as.double(pheno_data[[trait]])
#'
#'       results_eval_metrics_reps <- matrix(NA, nrow = repp, ncol = length(eval_metrics)+1)
#'       results_eval_metrics_reps[, 1] <- 1
#'       rownames(results_eval_metrics_reps) <- paste("REP", 1, sep = "_")
#'       colnames(results_eval_metrics_reps) <- c("Rep", eval_metrics)
#'       results_eval_metrics_reps <- as.data.frame(results_eval_metrics_reps)
#'     }
#'
#'     if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'       if (is.null(heter_groups)) {
#'         stop(message(paste(msg, 'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)
#'       }
#'       ENV <- as.character(unique(pheno_data[[heter_groups]]))
#'       ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
#'       colnames(ypred_cv) <- c("y", "yhat", heter_groups)
#'       ypred_cv <- as.data.frame(ypred_cv)
#'       ypred_cv[["y"]] <- as.double(pheno_data[[trait]])
#'       ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
#'
#'       results_eval_metrics_reps <- matrix(NA, nrow = length(ENV), ncol = length(eval_metrics)+2)
#'       results_eval_metrics_reps[, 1] <- rep(1, length(ENV))
#'       results_eval_metrics_reps[, 2] <- ENV
#'       colnames(results_eval_metrics_reps) <- c("Rep", heter_groups, eval_metrics)
#'       results_eval_metrics_reps_use <- results_eval_metrics_reps
#'       results_eval_metrics_reps <- data.frame()
#'     }
#'
#'     group <- test_set_val[[repp]]
#'
#'     handle_error <- FALSE
#'
#'     if (!cross_validation_meth %in% holds_out_methods_avail) {
#'       for (j in 1:nfolds) {
#'         yNA <- y
#'         for (g in 1:len_y) {
#'           if (group[g] == j) { yNA[g] <- NA }
#'         }
#'         tst <- which(is.na(yNA))
#'
#'         if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
#'           tryCatch({
#'             if (model == "GBLUP_BRR") {
#'               model_GBLUP <- "BRR"
#'               additional_params$bayes_model <-  model_GBLUP #model
#'               additional_params$bayes_trait <- trait
#'               additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
#'               additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
#'               #rm(model_GBLUP)
#'             } else {
#'               additional_params$bayes_model <- model
#'               additional_params$bayes_trait <- trait
#'               additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
#'               additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
#'             }
#'             ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
#'                                                         tst = tst, additional_params = additional_params)
#'           }, error = function(e) {
#'             message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'
#'         if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
#'           tryCatch({
#'             additional_params$response <- trait
#'             preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
#'             if (!is.null(preds)) {
#'               ypred_cv[tst, "yhat"] <- preds
#'             } else {
#'               stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
#'             }
#'           }, error = function(e) {
#'             message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'
#'         if (model %in% c(AI_valid_models)) {
#'           tryCatch({
#'             ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = omics_data,
#'                                                         tst = tst, additional_params = additional_params)
#'
#'           }, error = function(e) {
#'             message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'         ##
#'       }
#'     } else {
#'       if (cross_validation_meth %in% holds_out_methods_avail) {
#'         tst <- test_set_val[[repp]]
#'         yNA <- y
#'         yNA[tst] <- NA
#'
#'         if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
#'           tryCatch({
#'             if (model == "GBLUP_BRR") {
#'               model_GBLUP <- "BRR"
#'               additional_params$bayes_model <- model_GBLUP#model
#'               additional_params$bayes_trait <- trait
#'               additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
#'               additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
#'               #rm(model_GBLUP)
#'             } else {
#'               additional_params$bayes_model <- model
#'               additional_params$bayes_trait <- trait
#'               additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
#'               additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
#'             }
#'             ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
#'                                                         tst = tst, additional_params = additional_params)
#'           }, error = function(e) {
#'             message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'
#'         if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
#'           tryCatch({
#'             additional_params$response <- trait
#'             preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
#'             if (!is.null(preds)) {
#'               ypred_cv[tst, "yhat"] <- preds
#'             } else {
#'               stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
#'             }
#'           }, error = function(e) {
#'             message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'
#'         if (model %in% c(AI_valid_models)) {
#'           tryCatch({
#'             ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = omics_data,
#'                                                         tst = tst, additional_params = additional_params)
#'
#'             }, error = function(e) {
#'             message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
#'             handle_error <<- TRUE
#'           })
#'         }
#'         ##
#'       }
#'     }
#'
#'     if (isTRUE(handle_error)) {
#'       #return(NULL)
#'       ypred_cv <- NULL
#'       results_eval_metrics_reps <- NULL
#'     }
#'
#'     if (cross_validation_meth %in% CVs_multi_envs_methods_avail & isFALSE(handle_error)) {
#'       ypred_cv <- as.data.frame(ypred_cv)
#'       ypred_cv[, 'y'] <- as.double(ypred_cv[, 'y'])
#'       ypred_cv[, 'yhat'] <- as.double(ypred_cv[, 'yhat'])
#'
#'       results_eval_metrics_reps_use <- as.data.frame(results_eval_metrics_reps_use)
#'
#'       for (eva in 1:length(eval_metrics)) {
#'         sik <- ypred_cv |>
#'           dplyr::group_by(!!dplyr::sym(heter_groups)) |>
#'           dplyr::summarise(
#'             eval_metric = tryCatch(
#'               evaluation_metrics(yhat, y, eval_metrics = eval_metrics[eva]),
#'               error = function(e) NULL
#'             )
#'           )
#'
#'         if (!is.null(sik)) {
#'           for (ii in 1:nrow(sik)) {
#'             if (!is.null(sik$eval_metric[ii])) {
#'               results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups] == as.character(sik[[heter_groups]])[ii], c("Rep", eval_metrics[eva])] <- c(repp, sik$eval_metric[ii])
#'             }
#'           }
#'         }
#'       }
#'
#'       results_eval_metrics_reps <- rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)
#'     } else {
#'       if (!cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'         for (eva in 1:length(eval_metrics)) {
#'           tryCatch({
#'             results_eval_metrics_reps[repp, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
#'                                                                                      y_predicted = ypred_cv[tst, "yhat"],
#'                                                                                      eval_metrics = eval_metrics[eva])
#'           }, error = function(e) {
#'             results_eval_metrics_reps[repp, eval_metrics[eva]] <- NA
#'           })
#'         }
#'       }
#'     }
#'
#'     #if(model == "deep_learning_model") model <- as.character(task_row$modell)
#'
#'     list(trait = trait,
#'          rep = rep,
#'          model = model,
#'          eval_metrics_reps = results_eval_metrics_reps,
#'          ypred_cv_Reps_all = ypred_cv)
#'   })
#'
#'   future::plan("sequential")
#'
#'   dp_present <- intersect(GS_model_cv, dp_models)
#'
#'   if (length(dp_present) > 0) {
#'
#'     # build *just* the DP grid
#'     tasks_dl <- expand.grid(
#'       response    = response,
#'       replication = seq_len(replication),
#'       modell      = dp_present,
#'       stringsAsFactors = FALSE
#'     )
#'     tasks_by_trait <- split(tasks_dl, tasks_dl$response)
#'
#'   # 3) sequential loop by trait
#'   results_seq_dl <- lapply(tasks_by_trait, function(task_df) {
#'     lapply(seq_len(nrow(task_df)), function(i) {
#'
#'       # ─── extract and rename ─────────────────────────────────
#'       row    <- task_df[i, ]
#'       trait  <- row$response
#'       rep_i  <- as.integer(row$replication)
#'       model  <- as.character(row$modell)
#'
#'       # ─── log current state ───────────────────────────────────
#'       # message(sprintf(
#'       #   "Processing trait='%s', replication=%d, model='%s'",
#'       #   trait, rep_i, model
#'       # ))
#'
#'       # ─── force DL branch ───────────────────────────────────
#'       #additional_params$deep_learning_model <- model
#'       #additional_params$model_type <- model
#'       #model_flag <- "deep_learning_model"       # branch flag
#'       repp       <- 1
#'
#'       # ─── reproducible seed ─────────────────────────────────
#'       base_seed <- if (
#'         !is.null(random_state) &&
#'         is.numeric(random_state) &&
#'         length(random_state) == 1
#'       ) as.integer(random_state) else 123L
#'       new_seed <- (base_seed + rep_i * 10000L) %% .Machine$integer.max
#'
#'       # ─── train/test split ──────────────────────────────────
#'       if (cross_validation_meth %in% holds_out_methods_avail) {
#'         test_set_val <- hold_out_stratified_and_un(
#'           pheno_data      = pheno_data,
#'           gen_name        = gen_name,
#'           response        = trait,
#'           test_size       = test_size,
#'           random_state    = new_seed,
#'           replication     = repp,
#'           sampling_method = sampling_method
#'         )
#'
#'       } else if (cross_validation_meth %in% Kfolds_methods_avail) {
#'         test_set_val <- kfolds_stratified_un(
#'           pheno_data      = pheno_data,
#'           gen_name        = gen_name,
#'           response        = trait,
#'           #test_size       = test_size,
#'           nfolds          = nfolds,
#'           random_state    = new_seed,
#'           replication     = repp
#'           #sampling_method = sampling_method
#'         )
#'
#'       } else if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'         CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
#'         test_set_val <- CV0_CV1_CV2_for_multi_environment(
#'           pheno_data      = pheno_data,
#'           gen_name        = gen_name,
#'           response        = trait,
#'           #test_size       = test_size,
#'           CV              = CV,
#'           nfolds          = nfolds,
#'           heter_groups    = heter_groups,
#'           random_state    = new_seed,
#'           replication     = repp
#'           #sampling_method = sampling_method
#'         )
#'
#'       } else {
#'         stop(
#'           paste(
#'             msg,
#'             "Unsupported cross-validation method. Choose from:",
#'             paste(c(holds_out_methods_avail,
#'                     Kfolds_methods_avail,
#'                     CVs_multi_envs_methods_avail),
#'                   collapse = ", ")
#'           ),
#'           call. = FALSE
#'         )
#'       }
#'
#'       # ─── init y + containers ───────────────────────────────
#'       len_y <- nrow(pheno_data)
#'       y     <- as.double(pheno_data[[trait]])
#'
#'       if (cross_validation_meth %in% c(Kfolds_methods_avail, holds_out_methods_avail)) {
#'         ypred_cv <- data.frame(
#'           y    = y,
#'           yhat = rep(NA_real_, len_y)
#'         )
#'         results_eval_metrics_reps <- data.frame(
#'           Rep = 1,
#'           matrix(
#'             NA_real_,
#'             nrow = 1,
#'             ncol = length(eval_metrics),
#'             dimnames = list(NULL, eval_metrics)
#'           )
#'         )
#'       }
#'
#'       if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'         if (is.null(heter_groups)) {
#'           stop(paste(msg, "For CV1/CV2 you must supply `heter_groups`."), call. = FALSE)
#'         }
#'         ENV <- unique(as.character(pheno_data[[heter_groups]]))
#'         ypred_cv <- data.frame(
#'           y    = y,
#'           yhat = rep(NA_real_, len_y),
#'           env  = as.character(pheno_data[[heter_groups]])
#'         )
#'         results_eval_metrics_reps_use <- data.frame(
#'           Rep = rep(1, length(ENV)),
#'           env = ENV,
#'           matrix(NA_real_,
#'                  nrow = length(ENV),
#'                  ncol = length(eval_metrics),
#'                  dimnames = list(NULL, eval_metrics))
#'         )
#'         results_eval_metrics_reps <- data.frame()  # to be filled later
#'       }
#'
#'       # ─── grab the single fold/hold‑out group ────────────────
#'       group <- test_set_val[[1]]    # always [[1]] because repp == 1
#'
#'       handle_error <- FALSE
#'
#'       if (!cross_validation_meth %in% holds_out_methods_avail) {
#'         for (j in 1:nfolds) {
#'           yNA <- y
#'           for (g in 1:len_y) {
#'             if (group[g] == j) { yNA[g] <- NA }
#'           }
#'           tst <- which(is.na(yNA))
#'
#'           if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
#'             tryCatch({
#'               if (model == "GBLUP_BRR") {
#'                 model_GBLUP <- "BRR"
#'                 additional_params$bayes_model <-  model_GBLUP #model
#'                 additional_params$bayes_trait <- trait
#'                 additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
#'                 additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
#'                 #rm(model_GBLUP)
#'               } else {
#'                 additional_params$bayes_model <- model
#'                 additional_params$bayes_trait <- trait
#'                 additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
#'                 additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
#'               }
#'               ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
#'                                                           tst = tst, additional_params = additional_params)
#'             }, error = function(e) {
#'               message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'
#'           if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
#'             tryCatch({
#'               additional_params$response <- trait
#'               preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
#'               if (!is.null(preds)) {
#'                 ypred_cv[tst, "yhat"] <- preds
#'               } else {
#'                 stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
#'               }
#'             }, error = function(e) {
#'               message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'
#'           # if (model %in% c(AI_valid_models)) {
#'           if (model %in% dp_models) {
#'             tryCatch({
#'               ypred_cv[tst, "yhat"] <- predict_with_model(
#'                 model = model,
#'                 y = yNA,
#'                 omics_data = omics_data,
#'                 tst = tst,
#'                 additional_params = additional_params
#'               )
#'             }, error = function(e) {
#'               message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'           ##
#'         }
#'       } else {
#'         if (cross_validation_meth %in% holds_out_methods_avail) {
#'           tst <- test_set_val[[1]]
#'           yNA <- y
#'           yNA[tst] <- NA
#'
#'           if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
#'             tryCatch({
#'               if (model == "GBLUP_BRR") {
#'                 model_GBLUP <- "BRR"
#'                 additional_params$bayes_model <- model_GBLUP#model
#'                 additional_params$bayes_trait <- trait
#'                 additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
#'                 additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
#'                 #rm(model_GBLUP)
#'               } else {
#'                 additional_params$bayes_model <- model
#'                 additional_params$bayes_trait <- trait
#'                 additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
#'                 additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
#'               }
#'               ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
#'                                                           tst = tst, additional_params = additional_params)
#'             }, error = function(e) {
#'               message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'
#'           if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
#'             tryCatch({
#'               additional_params$response <- trait
#'               preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
#'               if (!is.null(preds)) {
#'                 ypred_cv[tst, "yhat"] <- preds
#'               } else {
#'                 stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
#'               }
#'             }, error = function(e) {
#'               message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'
#'           # if (model %in% c(AI_valid_models)) {
#'           if (model %in% dp_models) {
#'             tryCatch({
#'               ypred_cv[tst, "yhat"] <- predict_with_model(
#'                 model = model,
#'                 y = yNA,
#'                 omics_data = omics_data,
#'                 tst = tst,
#'                 additional_params = additional_params)
#'
#'             }, error = function(e) {
#'               message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
#'               handle_error <<- TRUE
#'             })
#'           }
#'           ##
#'         }
#'       }
#'
#'       if (isTRUE(handle_error)) {
#'         #return(NULL)
#'         ypred_cv <- NULL
#'         results_eval_metrics_reps <- NULL
#'       }
#'
#'       if (cross_validation_meth %in% CVs_multi_envs_methods_avail & isFALSE(handle_error)) {
#'         ypred_cv <- as.data.frame(ypred_cv)
#'         ypred_cv[, 'y'] <- as.double(ypred_cv[, 'y'])
#'         ypred_cv[, 'yhat'] <- as.double(ypred_cv[, 'yhat'])
#'
#'         results_eval_metrics_reps_use <- as.data.frame(results_eval_metrics_reps_use)
#'
#'         for (eva in 1:length(eval_metrics)) {
#'           sik <- ypred_cv |>
#'             dplyr::group_by(!!dplyr::sym(heter_groups)) |>
#'             dplyr::summarise(
#'               eval_metric = tryCatch(
#'                 evaluation_metrics(yhat, y, eval_metrics = eval_metrics[eva]),
#'                 error = function(e) NULL
#'               )
#'             )
#'
#'           if (!is.null(sik)) {
#'             for (ii in 1:nrow(sik)) {
#'               if (!is.null(sik$eval_metric[ii])) {
#'                 results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups] == as.character(sik[[heter_groups]])[ii], c("Rep", eval_metrics[eva])] <- c(repp, sik$eval_metric[ii])
#'               }
#'             }
#'           }
#'         }
#'
#'         results_eval_metrics_reps <- rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)
#'       } else {
#'         if (!cross_validation_meth %in% CVs_multi_envs_methods_avail) {
#'           for (eva in 1:length(eval_metrics)) {
#'             tryCatch({
#'               results_eval_metrics_reps[repp, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
#'                                                                                        y_predicted = ypred_cv[tst, "yhat"],
#'                                                                                        eval_metrics = eval_metrics[eva])
#'             }, error = function(e) {
#'               results_eval_metrics_reps[repp, eval_metrics[eva]] <- NA
#'             })
#'           }
#'         }
#'       }
#'
#'       # if(model == "deep_learning_model") model <- as.character(task_row$modell)
#'
#'       list(
#'         trait = trait,
#'         rep = rep_i,
#'         model = model,
#'         eval_metrics_reps = results_eval_metrics_reps,
#'         ypred_cv_Reps_all = ypred_cv
#'       )
#'     })
#'   })
#'
#'   # flatten nested list
#'   results_seq_dl_flat <- unlist(results_seq_dl, recursive = FALSE)
#'
#'   } else {
#'     # no DP models requested → skip that pass
#'     results_seq_dl_flat <- list()
#'   }
#'
#'   # combine into one flat list
#'   final_results <- c(results, results_seq_dl_flat)
#'
#'   # now filter those that have no NAs
#'   filtered_results <- Filter(
#'     function(res) {
#'       !is.null(res$eval_metrics_reps) &&
#'         !is.null(res$ypred_cv_Reps_all) &&
#'         all(!is.na(unlist(res$eval_metrics_reps))) &&
#'         all(!is.na(unlist(res$ypred_cv_Reps_all)))
#'     },
#'     final_results
#'   )
#'
#'   future::plan("sequential")
#'   return(filtered_results)
#' }
