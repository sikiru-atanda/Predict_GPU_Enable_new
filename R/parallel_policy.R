
# Initialize BLAS settings
init_single_thread_blas <- function() {
  blas_threads <- as.integer(Sys.getenv("GP_BLAS_THREADS", "1"))
  Sys.setenv(
    OMP_NUM_THREADS = blas_threads,
    OPENBLAS_NUM_THREADS = blas_threads,
    MKL_NUM_THREADS = blas_threads,
    BLIS_NUM_THREADS = blas_threads,
    VECLIB_MAXIMUM_THREADS = blas_threads
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch({
      RhpcBLASctl::blas_set_num_threads(blas_threads)
      RhpcBLASctl::omp_set_num_threads(blas_threads)
    }, error = function(e) logger::log_warn("Failed to set BLAS threads: {conditionMessage(e)}"))
  } else if (Sys.info()[["sysname"]] != "Windows") {
    logger::log_debug("RhpcBLASctl not available, BLAS threading may be uncontrolled")
  }
}

# Detect hidden threads
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

gp_parallel_globals_size_bytes <- function(globals) {
  if (!is.list(globals) || !length(globals)) {
    return(0)
  }
  sum(vapply(globals, function(x) {
    if (is.null(x)) 0 else as.numeric(utils::object.size(x))
  }, numeric(1)), na.rm = TRUE)
}

gp_parallel_has_nonnull_global <- function(globals, pattern) {
  if (!is.list(globals) || !length(globals) || is.null(names(globals))) {
    return(FALSE)
  }
  idx <- grepl(pattern, names(globals), ignore.case = TRUE)
  any(idx & !vapply(globals, is.null, logical(1)))
}

gp_parallel_phenotype_n <- function(globals) {
  if (!is.list(globals) || !length(globals)) {
    return(NA_integer_)
  }
  candidates <- globals[grepl("pheno", names(globals), ignore.case = TRUE)]
  for (candidate in candidates) {
    if (is.data.frame(candidate) || is.matrix(candidate)) {
      return(NROW(candidate))
    }
    if (is.list(candidate)) {
      nested <- candidate[vapply(candidate, function(x) is.data.frame(x) || is.matrix(x), logical(1))]
      if (length(nested)) {
        return(NROW(nested[[1L]]))
      }
    }
  }
  NA_integer_
}

gp_parallel_memory_budget_gb <- function() {
  configured <- suppressWarnings(as.numeric(Sys.getenv("GP_PAR_MEMORY_BUDGET_GB", "")))
  if (is.finite(configured) && configured > 0) {
    return(configured)
  }
  fraction <- suppressWarnings(as.numeric(Sys.getenv("GP_PAR_MEMORY_FRACTION", "0.7")))
  if (!is.finite(fraction) || fraction <= 0 || fraction > 1) {
    fraction <- 0.7
  }
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(Inf)
  }
  available <- tryCatch(as.numeric(ps::ps_system_memory()[["avail"]]), error = function(e) NA_real_)
  if (!is.finite(available) || available <= 0) {
    return(Inf)
  }
  (available / 1024^3) * fraction
}

gp_parallel_memory_worker_cap <- function(globals_gb, memory_budget_gb, worker_overhead_gb = NULL) {
  globals_gb <- suppressWarnings(as.numeric(globals_gb)[1L])
  memory_budget_gb <- suppressWarnings(as.numeric(memory_budget_gb)[1L])
  if (!is.finite(globals_gb) || globals_gb < 0 || !is.finite(memory_budget_gb)) {
    return(.Machine$integer.max)
  }
  payload_multiplier <- suppressWarnings(as.numeric(Sys.getenv("GP_PAR_PAYLOAD_MULTIPLIER", "4")))
  if (!is.finite(payload_multiplier) || payload_multiplier < 1) {
    payload_multiplier <- 4
  }
  worker_overhead_gb <- suppressWarnings(as.numeric(
    worker_overhead_gb %||% Sys.getenv("GP_PAR_WORKER_OVERHEAD_GB", "0.75")
  )[1L])
  if (!is.finite(worker_overhead_gb) || worker_overhead_gb < 0) {
    worker_overhead_gb <- 0.75
  }
  per_worker_gb <- globals_gb * payload_multiplier + worker_overhead_gb
  if (per_worker_gb <= 0) {
    return(.Machine$integer.max)
  }
  max(1L, as.integer(floor(memory_budget_gb / per_worker_gb)))
}

# Models fitted in R (BGLR / ASReml) need no Python child process per worker.
gp_parallel_r_only_models <- function() {
  c("BayesA", "BayesB", "BayesC", "BL", "BRR", "GBLUP_BRR", "RKHS", "GBLUP")
}

# Measured memory of one idle parallel worker (R + PredictProR loaded) plus,
# when Python-backed models run, one Python child with numpy/sklearn loaded.
# Measured once per session; NULL when it cannot be measured.
gp_parallel_measure_worker_gb <- function(python_models = TRUE) {
  key <- paste0("worker_memory_measured::python=", isTRUE(python_models))
  if (exists(key, envir = PredictProR_runtime_cache, inherits = FALSE)) {
    return(get(key, envir = PredictProR_runtime_cache, inherits = FALSE))
  }
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(NULL)
  }
  # A working worker has PredictProR's imports loaded (BGLR, Matrix,
  # data.table, dplyr, ...), so load the same set before measuring.
  imports <- tryCatch({
    raw <- utils::packageDescription("PredictProR")[["Imports"]] %||% ""
    pkgs <- trimws(sub("\\(.*$", "", strsplit(raw, ",")[[1L]]))
    pkgs[nzchar(pkgs)]
  }, error = function(e) character())
  r_gb <- tryCatch({
    cl <- parallel::makePSOCKcluster(1L)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mem <- parallel::clusterCall(cl, function(libs, pkgs) {
      .libPaths(unique(c(libs, .libPaths())))
      suppressMessages(suppressWarnings(
        for (p in unique(c("PredictProR", pkgs))) requireNamespace(p, quietly = TRUE)
      ))
      invisible(gc())
      if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
      info <- ps::ps_memory_info()
      # private bytes on Windows (committed memory), RSS elsewhere
      as.numeric(info[["mem_private"]] %||% info[["rss"]])
    }, unique(c(gp_predictpror_loaded_lib(), .libPaths())), imports)[[1L]]
    as.numeric(mem) / 1024^3
  }, error = function(e) NA_real_)
  if (!is.finite(r_gb) || r_gb <= 0) {
    return(NULL)
  }
  py_gb <- 0
  if (isTRUE(python_models)) {
    py_bin <- tryCatch(gp_detect_ml_python(), error = function(e) NULL)
    py_gb <- 0.5
    if (!is.null(py_bin) && file.exists(py_bin)) {
      out <- tryCatch(suppressWarnings(system2(
        py_bin,
        c("-c", shQuote(paste(
          "import numpy, sklearn",
          "try:\n import psutil, os; print(psutil.Process(os.getpid()).memory_info().rss)\nexcept Exception:\n print('NA')",
          sep = "\n"
        ))),
        stdout = TRUE, stderr = FALSE
      )), error = function(e) character())
      val <- suppressWarnings(as.numeric(utils::tail(out, 1L)))
      if (length(val) && is.finite(val) && val > 0) py_gb <- val / 1024^3
    }
  }
  # 25% headroom for transient allocations during a fit.
  measured <- (r_gb + py_gb) * 1.25
  assign(key, measured, envir = PredictProR_runtime_cache)
  measured
}

# Per-worker memory allowance, in order of precedence:
#   1. worker_memory_gb supplied by the user (model_execute argument);
#   2. GP_PAR_WORKER_OVERHEAD_GB environment variable;
#   3. measured on this machine (once per session) when `measure` is TRUE;
#   4. 0.75 GiB fallback.
gp_parallel_resolve_worker_memory_gb <- function(models = character(),
                                                 worker_memory_gb = NULL,
                                                 measure = TRUE) {
  user <- suppressWarnings(as.numeric(worker_memory_gb)[1L])
  if (length(user) && is.finite(user) && user > 0) {
    return(list(gb = user, source = "user_worker_memory_gb"))
  }
  env <- suppressWarnings(as.numeric(Sys.getenv("GP_PAR_WORKER_OVERHEAD_GB", "")))
  if (is.finite(env) && env >= 0) {
    return(list(gb = env, source = "env_GP_PAR_WORKER_OVERHEAD_GB"))
  }
  if (isTRUE(measure)) {
    python_models <- any(!as.character(models) %in% gp_parallel_r_only_models())
    measured <- tryCatch(gp_parallel_measure_worker_gb(python_models), error = function(e) NULL)
    if (!is.null(measured) && is.finite(measured) && measured > 0) {
      return(list(gb = measured, source = if (python_models) "measured_r_worker_plus_python" else "measured_r_worker"))
    }
  }
  list(gb = 0.75, source = "default_0.75")
}

gp_parallel_worker_memory_gb <- function(globals_gb,
                                         payload_multiplier = NULL,
                                         worker_overhead_gb = NULL) {
  globals_gb <- suppressWarnings(as.numeric(globals_gb)[1L])
  if (!is.finite(globals_gb) || globals_gb < 0) globals_gb <- 0
  payload_multiplier <- payload_multiplier %||%
    suppressWarnings(as.numeric(Sys.getenv("GP_PAR_PAYLOAD_MULTIPLIER", "4")))
  worker_overhead_gb <- worker_overhead_gb %||%
    suppressWarnings(as.numeric(Sys.getenv("GP_PAR_WORKER_OVERHEAD_GB", "0.75")))
  if (!is.finite(payload_multiplier) || payload_multiplier < 1) payload_multiplier <- 4
  if (!is.finite(worker_overhead_gb) || worker_overhead_gb < 0) worker_overhead_gb <- 0.75
  globals_gb * payload_multiplier + worker_overhead_gb
}

#' Estimate parallel worker memory for PredictProR workloads
#'
#' Reports the serialized size of each supplied object and recommends a worker
#' cap using the same cross-platform memory policy as `model_execute()`.
#' The estimate includes a configurable payload expansion factor and fixed R
#' process/package overhead for every worker.
#'
#' @param objects Named list of phenotype, marker, kernel, omic, cache, or
#'   prepared-model objects used by worker tasks.
#' @param n_tasks Number of independent tasks.
#' @param workers Requested workers. Defaults to detected logical cores minus
#'   one, bounded by `n_tasks`.
#' @param memory_budget_gb Optional worker-memory budget in GiB. The default is
#'   derived from `GP_PAR_MEMORY_BUDGET_GB` or available system memory.
#' @param payload_multiplier Multiplier for transient copies and model-specific
#'   transformations. Defaults to `GP_PAR_PAYLOAD_MULTIPLIER` (4).
#' @param worker_overhead_gb Fixed memory allowance for each R worker. Defaults
#'   to `GP_PAR_WORKER_OVERHEAD_GB` (0.75 GiB: an idle R worker plus a Python child process).
#'
#' @return A list containing an object inventory, payload estimates, requested
#'   and recommended worker counts, and projected worker memory.
#' @export
parallel_memory_plan <- function(objects,
                                 n_tasks = 1L,
                                 workers = NULL,
                                 memory_budget_gb = NULL,
                                 payload_multiplier = NULL,
                                 worker_overhead_gb = NULL) {
  if (!is.list(objects)) {
    stop("objects must be a named list.", call. = FALSE)
  }
  if (is.null(names(objects))) {
    names(objects) <- paste0("object_", seq_along(objects))
  }
  n_tasks <- suppressWarnings(as.integer(n_tasks)[1L])
  if (!is.finite(n_tasks) || n_tasks < 1L) {
    stop("n_tasks must be a positive integer.", call. = FALSE)
  }
  requested <- if (is.null(workers)) {
    gp_detect_cores(logical = TRUE, reserve = 1L)
  } else {
    suppressWarnings(as.integer(workers)[1L])
  }
  if (!is.finite(requested) || requested < 1L) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  requested <- min(requested, n_tasks)
  budget <- memory_budget_gb %||% gp_parallel_memory_budget_gb()
  budget <- suppressWarnings(as.numeric(budget)[1L])
  sizes <- vapply(objects, function(x) {
    if (is.null(x)) 0 else as.numeric(utils::object.size(x))
  }, numeric(1))
  inventory <- data.frame(
    name = names(objects),
    class = vapply(objects, function(x) paste(class(x), collapse = "/"), character(1)),
    size_mb = sizes / 1024^2,
    stringsAsFactors = FALSE
  )
  inventory <- inventory[order(inventory$size_mb, decreasing = TRUE), , drop = FALSE]
  rownames(inventory) <- NULL
  payload_gb <- sum(sizes) / 1024^3
  per_worker_gb <- gp_parallel_worker_memory_gb(
    globals_gb = payload_gb,
    payload_multiplier = payload_multiplier,
    worker_overhead_gb = worker_overhead_gb
  )
  cap <- if (is.finite(budget) && per_worker_gb > 0) {
    max(1L, as.integer(floor(budget / per_worker_gb)))
  } else {
    .Machine$integer.max
  }
  recommended <- max(1L, min(requested, cap))
  list(
    inventory = inventory,
    payload_gb = payload_gb,
    per_worker_gb = per_worker_gb,
    memory_budget_gb = budget,
    workers_requested = requested,
    memory_worker_cap = cap,
    workers_recommended = recommended,
    projected_worker_memory_gb = recommended * per_worker_gb
  )
}

gp_parallel_deep_models <- function() {
  c("cnn", "resnet", "ft_transformer", "saint", "tabnet", "node",
    "deepfm", "dcnv2", "nam", "moe", "gp_dkl", "mlp_with_attention", "mlp")
}

gp_parallel_gp_models <- function() {
  if (exists("gp_lowrank_supported_models", mode = "function")) {
    return(gp_lowrank_supported_models())
  }
  c("KRR", "GP", "GP_FA", "LowRankGP")
}

gp_parallel_param_threads <- function(model, params = list()) {
  model <- as.character(model %||% "")
  p <- params %||% list()
  keys <- unique(c(
    "n_jobs", "nthread", "threads", "thread_count",
    "xgb_nthread", "lightgbm_nthread", "catboost_thread_count",
    "rf_n_jobs"
  ))
  vals <- suppressWarnings(vapply(keys, function(k) {
    if (is.null(p[[k]])) {
      NA_real_
    } else {
      as.numeric(p[[k]])
    }
  }, numeric(1)))
  val <- suppressWarnings(max(vals, na.rm = TRUE))
  # Tree models with an unset or single-thread setting are raised to this
  # process's share of the cores at fit time (gp_py_ml_normalize_model_params),
  # so the policy must count them as multi-threaded: run in a worker pool they
  # gained nothing (mice 5-fold CV: 600 s parallel vs 595 s sequential) and
  # held ~3x the memory.
  if (model %in% c("RandomForest", "Xgboost", "CatBoost", "LightGBM") &&
      (!is.finite(val) || val <= 1) &&
      !tolower(Sys.getenv("PREDICTPRO_ML_AUTO_THREADS", "true")) %in% c("false", "0", "no")) {
    return(as.integer(tryCatch(gp_ml_thread_share(), error = function(e) 1L)))
  }
  if (!is.finite(val)) {
    return(NA_integer_)
  }
  as.integer(val)
}

gp_parallel_extract_model_policy_params <- function(model, params = list()) {
  model <- as.character(model %||% "")
  p <- params %||% list()
  keep <- if (model %in% gp_parallel_deep_models()) {
    c("device")
  } else if (model %in% gp_parallel_gp_models()) {
    c("device", "gp_device", "gp_backend", "gp_output_level", "gp_varcomp_mode",
      "gp_return_se", "gp_full_vc", "gp_iters", "gp_lr", "gp_exact_fast_cv")
  } else {
    switch(
      model,
      "Xgboost" = c("xgb_nthread", "nthread", "threads", "device", "xgb_booster", "booster"),
      "CatBoost" = c("catboost_thread_count", "thread_count", "n_jobs", "device"),
      "LightGBM" = c("lightgbm_nthread", "nthread", "n_jobs", "device"),
      "RandomForest" = c("rf_n_jobs", "n_jobs", "threads"),
      "SupportVectorMachine" = c("n_jobs", "threads"),
      "Ridge_Regression" = c("n_jobs", "threads"),
      "Lasso" = c("n_jobs", "threads"),
      "PartialLeastSquare" = c("n_jobs", "threads"),
      character()
    )
  }
  out <- p[intersect(names(p), keep)]
  if (!length(out)) {
    return(list())
  }
  out
}

gp_parallel_device_pref <- function(params = list()) {
  p <- params %||% list()
  dev <- p$gp_device %||% p$device %||% p$torch_device %||% p$xgb_device %||% NULL
  if (is.null(dev)) {
    return(NULL)
  }
  tolower(as.character(dev[[1]]))
}

gp_parallel_xgb_gpu_selected <- function(model, params = list(), num_gpus = 0L) {
  if (!identical(as.character(model %||% ""), "Xgboost")) {
    return(FALSE)
  }
  if (as.integer(num_gpus %||% 0L) <= 0L) {
    return(FALSE)
  }
  p <- params %||% list()
  booster <- tolower(as.character(p$xgb_booster %||% p$booster %||% "gbtree"))
  if (identical(booster, "gblinear")) {
    return(FALSE)
  }
  device_pref <- gp_parallel_device_pref(p)
  !identical(device_pref, "cpu")
}

gp_parallel_model_capability_table <- function(models,
                                               model_params_list,
                                               num_gpus = 0L,
                                               hidden_thr = NA_integer_,
                                               has_asreml_prep = FALSE,
                                               rkhs_kernel_GB = NA_real_,
                                               rkhs_mem_threshold_GB = 1,
                                               always_seq_rkhs = TRUE,
                                               gpu_usage = NA_real_,
                                               gpu_busy_threshold = gp_gpu_busy_threshold()) {
  deep_models <- gp_parallel_deep_models()
  gp_models <- gp_parallel_gp_models()
  gpu_busy <- gp_gpu_is_busy(gpu_usage, gpu_busy_threshold)
  out <- lapply(seq_along(models), function(i) {
    model <- as.character(models[[i]])
    params <- model_params_list[[i]] %||% list()
    threads <- gp_parallel_param_threads(model, params)
    thread_sensitive <- isTRUE(!is.na(threads) && threads > 1L)
    device_pref <- gp_parallel_device_pref(params)
    gpu_capable <- model %in% deep_models || model %in% gp_models || identical(model, "Xgboost")
    gpu_selected <- FALSE
    gpu_exclusive <- FALSE
    device_route <- NA_character_
    reasons <- character()

    if (model %in% deep_models) {
      gpu_selected <- as.integer(num_gpus %||% 0L) > 0L && !identical(device_pref, "cpu")
      gpu_exclusive <- TRUE
      if (gpu_selected) {
        reasons <- c(reasons, "deep_gpu_selected")
      } else {
        reasons <- c(reasons, "deep_cpu_or_no_gpu")
      }
    }

    if (model %in% gp_models) {
      if (identical(device_pref, "cpu")) {
        device_route <- "cpu"
        reasons <- c(reasons, "gp_user_cpu")
      } else if (as.integer(num_gpus %||% 0L) <= 0L) {
        device_route <- "cpu"
        reasons <- c(reasons, "gp_cpu_no_gpu")
      } else if (isTRUE(gpu_busy)) {
        device_route <- "cpu"
        reasons <- c(reasons, "gp_cpu_gpu_busy")
      } else {
        gpu_selected <- TRUE
        gpu_exclusive <- TRUE
        device_route <- "auto"
        reasons <- c(reasons, "gp_gpu_selected")
      }
    }

    if (identical(model, "Xgboost")) {
      if (gp_parallel_xgb_gpu_selected(model, params, num_gpus = num_gpus)) {
        gpu_selected <- TRUE
        gpu_exclusive <- TRUE
        reasons <- c(reasons, "xgboost_gpu_selected")
      }
      if (thread_sensitive) {
        reasons <- c(reasons, "xgboost_internal_threads")
      }
    }
    if (identical(model, "CatBoost") && thread_sensitive) {
      reasons <- c(reasons, "catboost_internal_threads")
    }
    if (identical(model, "LightGBM") && thread_sensitive) {
      reasons <- c(reasons, "lightgbm_internal_threads")
    }
    if (identical(model, "RandomForest") && thread_sensitive) {
      reasons <- c(reasons, "randomforest_internal_jobs")
    }

    rkhs_hidden <- identical(model, "RKHS") && is.finite(hidden_thr) && hidden_thr > 1L
    rkhs_mem <- identical(model, "RKHS") && is.finite(rkhs_kernel_GB) && rkhs_kernel_GB > rkhs_mem_threshold_GB
    rkhs_forced_seq <- if (isTRUE(always_seq_rkhs)) {
      identical(model, "RKHS")
    } else {
      rkhs_hidden || rkhs_mem
    }
    if (rkhs_hidden) reasons <- c(reasons, "rkhs_hidden_threads")
    if (rkhs_mem) reasons <- c(reasons, "rkhs_memory_heavy")
    if (identical(model, "RKHS") && isTRUE(always_seq_rkhs)) reasons <- c(reasons, "rkhs_forced_seq")

    asreml_forced_seq <- identical(model, "GBLUP") && isTRUE(has_asreml_prep)
    if (asreml_forced_seq) reasons <- c(reasons, "asreml_prepared")

    single_gpu_exclusive <- isTRUE(gpu_exclusive) && as.integer(num_gpus %||% 0L) <= 1L
    if (single_gpu_exclusive) reasons <- c(reasons, "single_gpu_exclusive")

    forced_seq <- single_gpu_exclusive || asreml_forced_seq || thread_sensitive || rkhs_forced_seq

    data.frame(
      model = model,
      internal_threads = if (is.na(threads)) NA_integer_ else as.integer(threads),
      thread_sensitive = isTRUE(thread_sensitive),
      gpu_capable = isTRUE(gpu_capable),
      gpu_selected = isTRUE(gpu_selected),
      gpu_exclusive = isTRUE(gpu_exclusive),
      device_route = device_route,
      forced_sequential = isTRUE(forced_seq),
      reasons = paste(unique(reasons), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

gp_parallel_score_weights <- function() {
  list(
    tasks = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_TASKS", "0.9"))),
    workers = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WORKERS", "0.5"))),
    globals = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_GLOBALS", "0.35"))),
    windows_future = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WINDOWS_FUTURE", "0.9"))),
    windows_future_startup = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WINDOWS_FUTURE_STARTUP", "2.0"))),
    windows_mirai_penalty = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WINDOWS_MIRAI_PENALTY", "0.45"))),
    mirai_worker_pool = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_MIRAI_WORKER_POOL", "0.8"))),
    gpu_mirai = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_GPU_MIRAI", "3.0"))),
    gpu_future_penalty = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_GPU_FUTURE_PENALTY", "1.75"))),
    hidden_threads = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_HIDDEN_THREADS", "0.4"))),
    mori = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_MORI", "0.6"))),
    tiny_seq = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_TINY_SEQ", "2.2"))),
    threaded_seq = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_THREADED_SEQ", "1.8"))),
    overhead = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_OVERHEAD", "0.45"))),
    worker_startup = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WORKER_STARTUP", "3.0"))),
    fork_worker_startup = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_FORK_WORKER_STARTUP", "1.2"))),
    windows_worker_startup = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_WINDOWS_WORKER_STARTUP", "3.0"))),
    base_parallel_unix = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_BASE_PARALLEL_UNIX", "0.9"))),
    foreach_penalty = suppressWarnings(as.numeric(Sys.getenv("GP_PAR_SCORE_FOREACH_PENALTY", "0.35")))
  )
}

gp_parallel_backend_available <- function(backend, plan = NULL) {
  backend <- as.character(backend %||% "")
  if (identical(backend, "sequential")) {
    return(TRUE)
  }
  if (identical(backend, "future")) {
    return(requireNamespace("future", quietly = TRUE) &&
             requireNamespace("future.apply", quietly = TRUE))
  }
  if (identical(backend, "mirai")) {
    return(requireNamespace("mirai", quietly = TRUE) &&
             requireNamespace("purrr", quietly = TRUE))
  }
  if (identical(backend, "base_parallel")) {
    return(TRUE)
  }
  if (identical(backend, "foreach")) {
    return(requireNamespace("foreach", quietly = TRUE) &&
             requireNamespace("doParallel", quietly = TRUE))
  }
  FALSE
}

gp_parallel_small_cpu_queue <- function(signals) {
  !isTRUE(signals$any_gpu_selected) &&
    is.finite(signals$n_par_tasks) &&
    signals$n_par_tasks <= max(4L, 2L * as.integer(signals$workers %||% 1L))
}

gp_parallel_worker_startup_penalty <- function(backend, plan, signals, weights) {
  if (!isTRUE(gp_parallel_small_cpu_queue(signals))) {
    return(0)
  }
  if (isTRUE(signals$is_windows)) {
    return(weights$windows_worker_startup)
  }
  future_multicore <- requireNamespace("future", quietly = TRUE) && identical(plan, future::multicore)
  fork_style <- isTRUE(future_multicore) || identical(backend, "base_parallel")
  if (isTRUE(fork_style)) {
    return(weights$fork_worker_startup)
  }
  weights$worker_startup
}

gp_parallel_backend_plan <- function(backend,
                                     sys_name = Sys.info()[["sysname"]],
                                     prefer_fork = TRUE) {
  if (!identical(backend, "future")) {
    return(NULL)
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    return(NULL)
  }
  is_windows <- identical(tolower(sys_name), "windows")
  if (!is_windows && isTRUE(prefer_fork)) {
    return(future::multicore)
  }
  future::multisession
}

gp_parallel_mori_signal <- function(backend,
                                    plan = NULL,
                                    globals_gb = 0,
                                    workers = 1L,
                                    n_tasks = 1L) {
  if (!isTRUE(gp_mori_enabled()) || !isTRUE(gp_mori_available())) {
    return(0)
  }
  if (!isTRUE(gp_mori_backend_eligible(backend = backend, plan = plan))) {
    return(0)
  }
  globals_mb <- max(0, as.numeric(globals_gb) * 1024)
  min_mb <- gp_mori_min_mb()
  if (!is.finite(globals_mb) || !is.finite(min_mb) || globals_mb < min_mb) {
    return(0)
  }
  fanout <- max(1L, min(as.integer(workers %||% 1L), as.integer(n_tasks %||% 1L)))
  log1p(globals_mb / max(min_mb, 1e-6)) * log1p(fanout)
}

# Rough seconds for ONE CV fold fit, by model family. Coefficients were
# calibrated on the wheat benchmark (n_train ~ 225, 8,068 markers; i7-12700H):
# BayesB 13.5 s, GBLUP_BRR 0.5 s, RandomForest(300 trees) ~ 48 s, Python
# CLI ML ~ 3 s start-up. Only the ratio of work to start-up cost matters.
# Size of what serialising `fun` ships to a worker: every object in its
# enclosing frames up to the global or a namespace environment. Serialisation
# does not share vectors, so object.size() of each object (shared or not) is a
# close proxy (mice CV: 1.81 GB vs 1.79 GB serialised).
gp_parallel_closure_payload_gb <- function(fun) {
  if (!is.function(fun)) return(NA_real_)
  env <- environment(fun)
  total <- 0
  seen <- 0L
  while (!is.null(env) && !identical(env, globalenv()) && !identical(env, emptyenv()) &&
         !isNamespace(env) && !identical(env, baseenv()) && seen < 50L) {
    for (nm in ls(env, all.names = TRUE)) {
      obj <- get(nm, envir = env, inherits = FALSE)
      if (!is.function(obj)) total <- total + as.numeric(utils::object.size(obj))
    }
    env <- parent.env(env)
    seen <- seen + 1L
  }
  total / 1024^3
}

gp_parallel_estimate_fold_seconds <- function(model, n_train, n_markers, nIter = 6000, ntree = 500,
                                              iterations = 300) {
  model <- as.character(model)
  n <- max(as.numeric(n_train), 1)
  p <- max(as.numeric(n_markers), 1)
  nIter <- max(as.numeric(nIter %||% 6000), 1)
  py_start <- 3
  if (model %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")) return(2.5e-9 * nIter * n * p + 0.2)
  if (model %in% c("GBLUP_BRR", "RKHS")) return(3.3e-9 * nIter * n * n + 0.3)
  if (identical(model, "GBLUP")) return(1e-8 * n^3 + 1)
  if (identical(model, "RandomForest")) return(py_start + 8.8e-8 * max(as.numeric(ntree %||% 500), 1) * n * p)
  if (model %in% c("Xgboost", "CatBoost", "LightGBM")) return(py_start + 2e-9 * max(as.numeric(iterations %||% 300), 1) * n * p)
  if (model %in% c("Ridge_Regression", "Lasso", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors")) {
    return(py_start + 1e-7 * n * p)
  }
  if (model %in% gp_parallel_gp_models()) return(10)
  if (model %in% gp_parallel_deep_models()) return(30)
  py_start + 1e-7 * n * p
}

# Expected wall-clock for running the parallel-eligible tasks on `workers`
# workers versus sequentially; worker start-up on Windows (PSOCK, package
# load) is much more expensive than a Unix fork.
gp_parallel_cost_estimate <- function(models, task_cost, workers, globals_gb, is_windows) {
  if (!length(models)) {
    return(list(est_seq = 0, est_par = 0, task_seconds = numeric()))
  }
  folds_per_task <- max(1, as.numeric(task_cost$folds_per_task %||% 1))
  secs <- vapply(models, function(m) {
    folds_per_task * gp_parallel_estimate_fold_seconds(
      m, task_cost$n_train, task_cost$n_markers,
      nIter = task_cost$nIter, ntree = task_cost$ntree, iterations = task_cost$iterations
    )
  }, numeric(1))
  workers <- max(1L, as.integer(workers))
  startup <- if (isTRUE(is_windows)) 6 + 1.5 * workers else 1 + 0.2 * workers
  startup <- startup + 8 * max(as.numeric(globals_gb), 0)
  est_seq <- sum(secs)
  est_par <- max(est_seq / workers, max(secs)) + startup
  list(est_seq = est_seq, est_par = est_par, task_seconds = secs)
}

gp_parallel_candidate_backends <- function(user_mode = "auto") {
  user_mode <- match.arg(user_mode, c("auto", "future", "sequential", "base_parallel", "foreach", "mirai"))
  if (!identical(user_mode, "auto")) {
    return(user_mode)
  }
  c("sequential", "future", "mirai", "base_parallel", "foreach")
}

gp_parallel_score_backend <- function(backend, signals, weights) {
  backend <- as.character(backend %||% "")
  eligible <- isTRUE(signals$available[[backend]])
  plan <- signals$plans[[backend]]
  reasons <- character()
  score <- -Inf

  if (!eligible) {
    return(list(
      backend = backend,
      eligible = FALSE,
      plan = plan,
      score = score,
      reasons = "backend_unavailable"
    ))
  }

  tasks_signal <- log1p(max(signals$n_par_tasks, 0))
  workers_signal <- log1p(max(signals$workers, 1L))
  globals_signal <- log1p(max(signals$tot_GB, 0) + 1e-8)
  hidden_signal <- if (is.finite(signals$hidden_thr) && signals$hidden_thr > 1L) log1p(signals$hidden_thr - 1L) else 0
  mori_signal <- gp_parallel_mori_signal(
    backend = backend,
    plan = plan,
    globals_gb = signals$tot_GB,
    workers = signals$workers,
    n_tasks = signals$n_par_tasks
  )
  tasks_per_worker <- signals$n_par_tasks / max(signals$workers, 1L)
  small_cpu_queue <- gp_parallel_small_cpu_queue(signals)

  if (identical(backend, "sequential")) {
    score <- 0
    if (signals$n_par_tasks <= 2L || tasks_per_worker <= 1) {
      score <- score + weights$tiny_seq * (2 - min(tasks_per_worker, 2))
      reasons <- c(reasons, "tiny_workload")
    }
    if (signals$thread_sensitive_parallel) {
      score <- score + weights$threaded_seq
      reasons <- c(reasons, "thread_sensitive_parallel")
    }
    if (small_cpu_queue) {
      score <- score + weights$tiny_seq
      reasons <- c(reasons, "small_cpu_queue")
      if (signals$is_windows) {
        reasons <- c(reasons, "small_windows_cpu_queue")
      }
    }
    if (signals$any_gpu_selected && signals$num_gpus <= 1L) {
      score <- score + weights$gpu_mirai
      reasons <- c(reasons, "single_gpu_model")
    }
    score <- score - weights$overhead * tasks_signal
    return(list(
      backend = backend,
      eligible = TRUE,
      plan = plan,
      score = score,
      reasons = unique(reasons)
    ))
  }

  score <- 0.25 + weights$tasks * tasks_signal + weights$workers * workers_signal
  score <- score + weights$globals * globals_signal
  score <- score + weights$mori * mori_signal
  if (mori_signal > 0) {
    reasons <- c(reasons, "mori_shared_globals")
  }

  if (identical(backend, "future")) {
    reasons <- c(reasons, "general_psock_backend")
    if (signals$is_windows) {
      score <- score + weights$windows_future - weights$windows_future_startup
      reasons <- c(reasons, "windows_psock", "windows_future_startup_penalty")
    }
    if (signals$any_gpu_selected && signals$num_gpus > 1L) {
      score <- score - weights$gpu_future_penalty
      reasons <- c(reasons, "gpu_penalty")
    }
    if (signals$thread_sensitive_parallel) {
      score <- score - (weights$hidden_threads * hidden_signal * 0.5)
      reasons <- c(reasons, "thread_sensitive_penalty")
    }
  } else if (identical(backend, "mirai")) {
    reasons <- c(reasons, "low_overhead_worker_pool")
    if (!small_cpu_queue && !signals$any_gpu_selected) {
      score <- score + weights$mirai_worker_pool
      reasons <- c(reasons, "reusable_worker_pool")
    }
    if (signals$is_windows) {
      score <- score - weights$windows_mirai_penalty
      reasons <- c(reasons, "windows_socket_penalty")
    }
    if (signals$any_gpu_selected && signals$num_gpus > 1L) {
      score <- score + weights$gpu_mirai
      reasons <- c(reasons, "multi_gpu_model_bonus")
    }
    if (signals$n_par_tasks <= 2L) {
      score <- score - weights$overhead
      reasons <- c(reasons, "small_queue_penalty")
    }
  } else if (identical(backend, "base_parallel")) {
    reasons <- c(reasons, "base_parallel_fallback")
    if (!signals$is_windows) {
      score <- score + weights$base_parallel_unix
      reasons <- c(reasons, "unix_fork_bonus")
    } else {
      score <- score - weights$overhead
      reasons <- c(reasons, "windows_psock_penalty")
    }
  } else if (identical(backend, "foreach")) {
    reasons <- c(reasons, "foreach_fallback")
    score <- score - weights$foreach_penalty
  }

  startup_penalty <- gp_parallel_worker_startup_penalty(backend, plan, signals, weights)
  if (startup_penalty > 0) {
    score <- score - startup_penalty
    reasons <- c(reasons, "small_cpu_queue_penalty")
    if (signals$is_windows) {
      reasons <- c(reasons, "small_windows_cpu_queue_penalty")
    }
  }

  if (signals$thread_sensitive_parallel) {
    score <- score - weights$hidden_threads * hidden_signal
  }

  list(
    backend = backend,
    eligible = TRUE,
    plan = plan,
    score = score,
    reasons = unique(reasons)
  )
}

gp_parallel_choose_backend <- function(user_mode = "auto",
                                       signals,
                                       prefer_fork = TRUE,
                                       exclude = character()) {
  candidates <- setdiff(gp_parallel_candidate_backends(user_mode = user_mode), exclude)
  plans <- lapply(candidates, function(backend) {
    gp_parallel_backend_plan(
      backend = backend,
      sys_name = signals$sys_name,
      prefer_fork = prefer_fork
    )
  })
  names(plans) <- candidates
  available <- vapply(candidates, function(backend) {
    gp_parallel_backend_available(backend, plan = plans[[backend]])
  }, logical(1))

  signals$plans <- plans
  signals$available <- available
  weights <- gp_parallel_score_weights()
  scored <- lapply(candidates, function(backend) gp_parallel_score_backend(backend, signals, weights))

  table <- data.frame(
    backend = vapply(scored, `[[`, character(1), "backend"),
    eligible = vapply(scored, `[[`, logical(1), "eligible"),
    score = vapply(scored, `[[`, numeric(1), "score"),
    plan = vapply(scored, function(x) {
      if (is.null(x$plan)) {
        "NA"
      } else if (requireNamespace("future", quietly = TRUE) && identical(x$plan, future::multicore)) {
        "multicore"
      } else if (requireNamespace("future", quietly = TRUE) && identical(x$plan, future::multisession)) {
        "multisession"
      } else {
        "custom"
      }
    }, character(1)),
    reasons = vapply(scored, function(x) paste(x$reasons %||% character(), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )

  eligible_rows <- which(table$eligible)
  best_idx <- if (length(eligible_rows)) {
    eligible_rows[[which.max(table$score[eligible_rows])]]
  } else {
    1L
  }

  list(
    backend = table$backend[[best_idx]],
    plan = plans[[table$backend[[best_idx]]]],
    score_table = table,
    decision_reason = table$reasons[[best_idx]],
    backend_score = table$score[[best_idx]]
  )
}

# Decide parallel/sequential policy
sp_decide_policy <- function(models,
                             model_params_list,
                             globals,
                             n_tasks,
                             user_mode = c("auto", "future", "sequential", "base_parallel", "foreach", "mirai"),
                             num_cores = NULL,
                             globals_max_GB = as.numeric(Sys.getenv("GP_GLOBALS_MAX_GB", "4")),
                             rkhs_mem_threshold_GB = getOption("gp.rkhs.mem.threshold.GB", 1.0),
                             always_seq_rkhs = getOption("gp.always.seq.RKHS", TRUE),
                             verbose = TRUE,
                             sys_name = Sys.info()[["sysname"]],
                             prefer_fork = TRUE,
                             worker_memory_gb = NULL,
                             memory_budget_gb = NULL,
                             task_cost = NULL,
                             payload_gb = NULL) {
  user_mode <- match.arg(user_mode)
  n_tasks <- as.integer(n_tasks %||% length(models))

  `%||%` <- function(a, b) if (is.null(a)) b else a
  safe_num <- function(x) if (is.finite(x)) sprintf("%.2f", x) else "NA"

  if (identical(user_mode, "sequential")) {
    return(list(
      backend = "sequential",
      plan = NULL,
      prefer_fork = isTRUE(prefer_fork),
      workers = 1L,
      chunk_size = NULL,
      internal_flags = rep(TRUE, length(models)),
      num_gpus = 0L,
      gpu_ids = character(0),
      backend_score = Inf,
      decision_reason = "user_forced_sequential",
      backend_candidates = data.frame(
        backend = "sequential",
        eligible = TRUE,
        score = Inf,
        plan = "NA",
        reasons = "user_forced_sequential",
        stringsAsFactors = FALSE
      )
    ))
  }

  deep_models <- gp_parallel_deep_models()
  num_gpus <- detect_num_gpus()
  gpu_usage <- tryCatch(get_gpu_usage(), error = function(e) NA_real_)
  gpu_busy_threshold <- gp_gpu_busy_threshold()
  has_asreml_prep <- gp_parallel_has_nonnull_global(globals, "asreml")

  hidden_thr <- detect_hidden_threads()

  n <- gp_parallel_phenotype_n(globals)
  rkhs_kernel_GB <- if (is.finite(n)) (8 * as.numeric(n) * as.numeric(n)) / (1024^3) else NA_real_
  capability_table <- gp_parallel_model_capability_table(
    models = models,
    model_params_list = model_params_list,
    num_gpus = num_gpus,
    hidden_thr = hidden_thr,
    has_asreml_prep = has_asreml_prep,
    rkhs_kernel_GB = rkhs_kernel_GB,
    rkhs_mem_threshold_GB = rkhs_mem_threshold_GB,
    always_seq_rkhs = always_seq_rkhs,
    gpu_usage = gpu_usage,
    gpu_busy_threshold = gpu_busy_threshold
  )
  internal_flags <- capability_table$forced_sequential
  seq_deep <- models %in% deep_models & internal_flags
  seq_asreml <- models == "GBLUP" & grepl("asreml_prepared", capability_table$reasons, fixed = TRUE)
  seq_xgb <- models == "Xgboost" & capability_table$thread_sensitive
  seq_catboost <- models == "CatBoost" & capability_table$thread_sensitive
  seq_lightgbm <- models == "LightGBM" & capability_table$thread_sensitive
  seq_randomforest <- models == "RandomForest" & capability_table$thread_sensitive
  seq_RKHS <- models == "RKHS" & internal_flags

  tot_GB <- gp_parallel_globals_size_bytes(globals) / (1024^3)
  # What each worker really receives (the task closure's frames), when known.
  payload_gb <- suppressWarnings(as.numeric(payload_gb %||% NA_real_)[1L])
  if (is.finite(payload_gb) && (!is.finite(tot_GB) || payload_gb > tot_GB)) tot_GB <- payload_gb

  if (is.finite(tot_GB) && tot_GB > globals_max_GB && isTRUE(verbose)) {
    logger::log_warn("Globals ~ {safe_num(tot_GB)} GB > {safe_num(globals_max_GB)} GB. Consider using file paths or database keys")
  }
  if (is.finite(tot_GB) && tot_GB > globals_max_GB) {
    internal_flags[] <- TRUE
  }

  n_par_tasks <- sum(!internal_flags)
  if (n_par_tasks < 2L) {
    internal_flags[] <- TRUE
    n_par_tasks <- 0L
  }

  is_windows <- identical(tolower(sys_name), "windows")
  n_cores <- if (is.null(num_cores)) {
    gp_detect_cores(logical = TRUE, reserve = 1L)
  } else {
    max(1L, as.integer(num_cores))
  }
  workers <- if (n_par_tasks > 0) max(1L, min(n_cores, n_par_tasks)) else 1L
  workers_requested <- workers
  user_budget <- suppressWarnings(as.numeric(memory_budget_gb)[1L])
  memory_budget_GB <- if (length(user_budget) && is.finite(user_budget) && user_budget > 0) {
    user_budget
  } else {
    gp_parallel_memory_budget_gb()
  }
  # Measure the real worker footprint only when more than one worker could
  # run, so sequential runs never pay for starting a probe worker.
  worker_memory <- gp_parallel_resolve_worker_memory_gb(
    models = models[!internal_flags],
    worker_memory_gb = worker_memory_gb,
    measure = workers > 1L
  )
  # Per-worker memory = overhead + GP_PAR_PAYLOAD_MULTIPLIER (default 4) x the
  # payload each worker receives. Calibrated on system-level memory, mice CV
  # (1.0 GB payload): BGLR workers used ~4.4 GB (3.8x + 0.6 GB), mixed
  # BGLR/Python CV ~3.5x including the Python children -- so Python models are
  # no longer doubled on top. The former 1.25x let 4 BGLR workers take 16.7 GB
  # and leave 1 GB of available memory.
  data_GB <- tot_GB
  memory_worker_cap <- gp_parallel_memory_worker_cap(data_GB, memory_budget_GB, worker_memory$gb)
  per_worker_memory_GB <- gp_parallel_worker_memory_gb(data_GB, worker_overhead_gb = worker_memory$gb)
  memory_limited <- workers > memory_worker_cap
  if (isTRUE(memory_limited)) {
    workers <- max(1L, min(workers, memory_worker_cap))
  }
  memory_forced_sequential <- isTRUE(memory_limited) && n_par_tasks > 0L && workers < 2L
  if (isTRUE(memory_forced_sequential)) {
    internal_flags[] <- TRUE
    n_par_tasks <- 0L
    workers <- 1L
  }

  gpu_ids <- strsplit(Sys.getenv("GP_CUDA_VISIBLE_DEVICES", paste0("0:", max(0, num_gpus - 1))), ":")[[1]]
  if (num_gpus > 1L && any(capability_table$gpu_selected, na.rm = TRUE)) {
    workers <- max(1L, min(workers, length(gpu_ids)))
  }

  thread_sensitive_parallel <- any(capability_table$thread_sensitive | (models == "RKHS" & internal_flags), na.rm = TRUE)
  signals <- list(
    sys_name = sys_name,
    is_windows = is_windows,
    n_tasks = n_tasks,
    n_par_tasks = n_par_tasks,
    workers = workers,
    tot_GB = tot_GB,
    hidden_thr = hidden_thr,
    num_gpus = num_gpus,
    gpu_usage = gpu_usage,
    gpu_busy_threshold = gpu_busy_threshold,
    gpu_busy = gp_gpu_is_busy(gpu_usage, gpu_busy_threshold),
    any_deep = any(models %in% deep_models),
    any_gpu_selected = any(capability_table$gpu_selected, na.rm = TRUE),
    thread_sensitive_parallel = thread_sensitive_parallel
  )

  no_parallel_reason <- if (isTRUE(memory_forced_sequential)) {
    "memory_limited_workers"
  } else {
    "no_parallel_eligible_tasks"
  }
  if (n_par_tasks <= 0L) {
    forced_reason_strings <- capability_table$reasons[internal_flags]
    forced_reason_strings <- forced_reason_strings[nzchar(forced_reason_strings %||% "")]
    forced_reason_tokens <- unique(unlist(strsplit(forced_reason_strings, ";", fixed = TRUE)))
    forced_reason_tokens <- forced_reason_tokens[nzchar(forced_reason_tokens)]
    deep_forced_single_gpu <- num_gpus == 1L &&
      any(models %in% deep_models) &&
      any(capability_table$gpu_selected, na.rm = TRUE) &&
      all(internal_flags)
    if (isTRUE(memory_forced_sequential)) {
      no_parallel_reason <- "memory_limited_workers"
    } else if (isTRUE(deep_forced_single_gpu)) {
      no_parallel_reason <- "single_gpu_exclusive_deep_models"
    } else {
      specific_reasons <- c(
        "rkhs_hidden_threads",
        "xgboost_internal_threads",
        "catboost_internal_threads",
        "lightgbm_internal_threads",
        "randomforest_internal_jobs"
      )
      matched_specific <- intersect(specific_reasons, forced_reason_tokens)
      if (length(matched_specific) == 1L) {
        no_parallel_reason <- matched_specific[[1L]]
      } else if (length(matched_specific) > 1L) {
        no_parallel_reason <- "multiple_sequential_constraints"
      }
    }
  }

  if (n_par_tasks <= 0L) {
    backend_choice <- list(
      backend = "sequential",
      plan = NULL,
      backend_score = Inf,
      decision_reason = no_parallel_reason,
      score_table = data.frame(
        backend = "sequential",
        eligible = TRUE,
        score = Inf,
        plan = "NA",
        reasons = no_parallel_reason,
        stringsAsFactors = FALSE
      )
    )
  } else if (!identical(user_mode, "auto")) {
    backend_choice <- list(
      backend = user_mode,
      plan = gp_parallel_backend_plan(user_mode, sys_name = sys_name, prefer_fork = prefer_fork),
      backend_score = Inf,
      decision_reason = paste0("user_forced_", user_mode),
      score_table = data.frame(
        backend = user_mode,
        eligible = TRUE,
        score = Inf,
        plan = if (identical(user_mode, "future")) {
          if (!is_windows && isTRUE(prefer_fork)) "multicore" else "multisession"
        } else "NA",
        reasons = paste0("user_forced_", user_mode),
        stringsAsFactors = FALSE
      )
    )
  } else {
    cost <- NULL
    use_cost_model <- !is.null(task_cost) &&
      !tolower(Sys.getenv("GP_PAR_COST_MODEL", "true")) %in% c("false", "0", "no")
    if (isTRUE(use_cost_model)) {
      cost <- gp_parallel_cost_estimate(
        models = models[!internal_flags], task_cost = task_cost, workers = workers,
        globals_gb = tot_GB, is_windows = is_windows
      )
    }
    if (!is.null(cost) && cost$est_par < 0.75 * cost$est_seq) {
      # Expected work clearly outweighs worker start-up. Prefer base_parallel:
      # it ships the shared data once per worker (mirai ships it once per
      # task; future duplicated closure data and hit its globals limit) and
      # forks on Unix. Users can still force future/mirai/foreach.
      backend_choice <- if (gp_parallel_backend_available(
        "base_parallel", plan = gp_parallel_backend_plan("base_parallel", sys_name = sys_name, prefer_fork = prefer_fork)
      )) {
        list(
          backend = "base_parallel",
          plan = gp_parallel_backend_plan("base_parallel", sys_name = sys_name, prefer_fork = prefer_fork),
          backend_score = Inf,
          decision_reason = "base_parallel_preferred",
          score_table = data.frame(backend = "base_parallel", eligible = TRUE, score = Inf, plan = "NA",
                                   reasons = "base_parallel_preferred", stringsAsFactors = FALSE)
        )
      } else {
        gp_parallel_choose_backend(
          user_mode = user_mode, signals = signals, prefer_fork = prefer_fork,
          exclude = "sequential"
        )
      }
      backend_choice$decision_reason <- paste0(
        sprintf("cost_model_parallel(est_seq=%.0fs,est_par=%.0fs)", cost$est_seq, cost$est_par),
        ";", backend_choice$decision_reason %||% ""
      )
    } else if (!is.null(cost)) {
      backend_choice <- list(
        backend = "sequential",
        plan = NULL,
        backend_score = Inf,
        decision_reason = sprintf("cost_model_sequential(est_seq=%.0fs,est_par=%.0fs)", cost$est_seq, cost$est_par),
        score_table = data.frame(backend = "sequential", eligible = TRUE, score = Inf, plan = "NA",
                                 reasons = "cost_model_sequential", stringsAsFactors = FALSE)
      )
    } else {
      backend_choice <- gp_parallel_choose_backend(
        user_mode = user_mode,
        signals = signals,
        prefer_fork = prefer_fork
      )
    }
  }

  backend <- backend_choice$backend
  plan <- backend_choice$plan

  if (isTRUE(memory_limited) && !isTRUE(memory_forced_sequential)) {
    backend_choice$decision_reason <- paste(
      backend_choice$decision_reason %||% "",
      "memory_capped_workers",
      sep = if (nzchar(backend_choice$decision_reason %||% "")) ";" else ""
    )
  }

  chunk_size <- if (n_par_tasks > 0 && workers > 0) ceiling(n_par_tasks / workers) else NULL

  if (isTRUE(verbose)) {
    plan_name <- if (is.null(plan)) "NA" else if (identical(plan, future::multicore)) "multicore" else "multisession"
    tail_hidden <- if (is.finite(hidden_thr)) sprintf(" | hidden_thr=%s", safe_num(hidden_thr)) else ""
    tail_rkhs <- if (is.finite(rkhs_kernel_GB)) sprintf(" | rkhs=%s GB", safe_num(rkhs_kernel_GB)) else ""
    tail_gpus <- if (num_gpus > 0) sprintf(" | num_gpus=%d", num_gpus) else ""
    tail_gpu_usage <- if (is.finite(gpu_usage)) sprintf(" | gpu_usage=%.1f%%", gpu_usage) else ""
    tail_memory <- if (is.finite(memory_budget_GB)) {
      sprintf(" | memory_budget=%s GB worker_cap=%d", safe_num(memory_budget_GB), memory_worker_cap)
    } else ""
    logger::log_info(
      "Policy: backend={backend} plan={plan_name} workers={workers} chunk={if (is.null(chunk_size)) 'NA' else chunk_size} | seq={sum(internal_flags)} par={sum(!internal_flags)} | globals={safe_num(tot_GB)} GB{tail_hidden}{tail_rkhs}{tail_gpus}{tail_gpu_usage}{tail_memory} | reason={backend_choice$decision_reason %||% 'NA'}"
    )
    if (any(seq_deep))
      logger::log_info("Deep/GPU-exclusive forced seq: {paste(which(seq_deep), collapse = ',')}")
    if (any(seq_asreml))
      logger::log_info("GBLUP(asreml) forced seq: {paste(which(seq_asreml), collapse = ',')}")
    if (any(seq_xgb))
      logger::log_info("Xgboost(nthread>1) seq: {paste(which(seq_xgb), collapse = ',')}")
    if (any(seq_catboost))
      logger::log_info("CatBoost(thread_count>1) seq: {paste(which(seq_catboost), collapse = ',')}")
    if (any(seq_lightgbm))
      logger::log_info("LightGBM(nthread>1) seq: {paste(which(seq_lightgbm), collapse = ',')}")
    if (any(seq_randomforest))
      logger::log_info("RandomForest(n_jobs>1) seq: {paste(which(seq_randomforest), collapse = ',')}")
    if (any(models == "RKHS"))
      logger::log_info("RKHS policy: always_seq_rkhs={always_seq_rkhs}, threshold={safe_num(rkhs_mem_threshold_GB)} GB")
    if (is.data.frame(capability_table) && nrow(capability_table)) {
      logger::log_debug("Model capability summary: {paste(sprintf('%s[%s]', capability_table$model, capability_table$reasons), collapse = ' | ')}")
    }
  }

  list(
    backend = backend,
    plan = plan,
    prefer_fork = isTRUE(prefer_fork),
    workers = workers,
    workers_requested = workers_requested,
    globals_gb = tot_GB,
    per_worker_memory_gb = per_worker_memory_GB,
    # estimated seconds of each parallel task (same order as the parallel
    # tasks) so the scheduler can start the heaviest first
    task_seconds = if (!is.null(task_cost)) {
      gp_parallel_cost_estimate(models[!internal_flags], task_cost, workers, tot_GB, is_windows)$task_seconds
    } else NULL,
    worker_memory_gb = worker_memory$gb,
    worker_memory_source = worker_memory$source,
    projected_worker_memory_gb = workers * per_worker_memory_GB,
    memory_budget_gb = memory_budget_GB,
    memory_worker_cap = memory_worker_cap,
    n_tasks = n_tasks,
    chunk_size = chunk_size,
    internal_flags = internal_flags,
    num_gpus = num_gpus,
    gpu_ids = gpu_ids,
    gpu_usage = gpu_usage,
    gpu_busy_threshold = gpu_busy_threshold,
    gpu_busy = gp_gpu_is_busy(gpu_usage, gpu_busy_threshold),
    mori_enabled = gp_mori_enabled() && gp_mori_available(),
    mori_backend_eligible = gp_mori_backend_eligible(backend = backend, plan = plan),
    backend_score = backend_choice$backend_score %||% NA_real_,
    decision_reason = backend_choice$decision_reason %||% NA_character_,
    backend_candidates = backend_choice$score_table %||% data.frame(),
    model_capabilities = capability_table,
    gpu_parallel_enabled = isTRUE(any(capability_table$gpu_selected, na.rm = TRUE) && num_gpus > 1L)
  )
}
