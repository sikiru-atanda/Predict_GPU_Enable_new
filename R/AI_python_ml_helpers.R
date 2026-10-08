gp_ml_bridge_script_path <- function() {
  candidates <- c("inst/python/ml_bridge.py", "python/ml_bridge.py", "../inst/python/ml_bridge.py")
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) {
    return(normalizePath(hit[1], winslash = "/", mustWork = TRUE))
  }
  pyfile <- system.file("python/ml_bridge.py", package = "PredictProR")
  if (nzchar(pyfile)) {
    pyfile <- normalizePath(pyfile, winslash = "/", mustWork = TRUE)
  }
  if (!nzchar(pyfile)) {
    stop("ml_bridge.py not found under inst/python/.", call. = FALSE)
  }
  pyfile
}

PredictProR_ml_bridge_env <- new.env(parent = emptyenv())

gp_ml_bridge_module_key <- function(python_bin, bridge_path) {
  paste(
    normalizePath(python_bin, winslash = "/", mustWork = FALSE),
    normalizePath(bridge_path, winslash = "/", mustWork = FALSE),
    sep = "::"
  )
}

gp_get_ml_bridge_module <- function(python_bin = NULL) {
  python_bin <- python_bin %||% gp_detect_ml_python()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop(
      "No configured Python runtime was found for Python ML models.\n",
      "Set PREDICTPRO_PYTHON to a Python with numpy, pandas and scikit-learn ",
      "(plus xgboost/catboost/lightgbm for those models).",
      call. = FALSE
    )
  }
  bridge_path <- gp_ml_bridge_script_path()
  gp_init_python_once(python_bin)
  key <- gp_ml_bridge_module_key(python_bin, bridge_path)
  if (exists(key, envir = PredictProR_ml_bridge_env, inherits = FALSE)) {
    return(get(key, envir = PredictProR_ml_bridge_env, inherits = FALSE))
  }
  importlib <- reticulate::import("importlib.util", delay_load = FALSE, convert = FALSE)
  spec <- importlib$spec_from_file_location("predictpror_ml_bridge", bridge_path)
  mod <- importlib$module_from_spec(spec)
  spec$loader$exec_module(mod)
  assign(key, mod, envir = PredictProR_ml_bridge_env)
  mod
}

gp_detect_ml_python <- function() {
  # Called for every ML fit; cache per session (keyed by the configuring env
  # vars) so a PREDICTPRO_PYTHON probe is not repeated for every CV fold.
  cache_key <- paste(
    "ml_python", Sys.getenv("PREDICTPRO_ML_PYTHON", unset = ""),
    Sys.getenv("PREDICTPRO_PYTHON_ML", unset = ""), Sys.getenv("PREDICTPRO_PYTHON", unset = ""),
    sep = "::"
  )
  if (exists(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE)) {
    cached <- get(cache_key, envir = PredictProR_runtime_cache, inherits = FALSE)
    if (file.exists(cached)) return(cached)
  }
  found <- gp_detect_ml_python_uncached()
  if (!is.null(found)) {
    assign(cache_key, found, envir = PredictProR_runtime_cache)
  }
  found
}

gp_detect_ml_python_uncached <- function() {
  explicit <- Sys.getenv("PREDICTPRO_ML_PYTHON", unset = "")
  if (!nzchar(explicit)) {
    explicit <- Sys.getenv("PREDICTPRO_PYTHON_ML", unset = "")
  }
  if (nzchar(explicit) && file.exists(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }

  generic <- Sys.getenv("PREDICTPRO_PYTHON", unset = "")
  if (nzchar(generic) && file.exists(generic)) {
    generic <- normalizePath(generic, winslash = "/", mustWork = TRUE)
    generic_probe <- tryCatch(
      gp_python_probe(generic, modules = c("numpy", "sklearn")),
      error = function(e) NULL
    )
    if (isTRUE(generic_probe[["sklearn"]])) {
      return(generic)
    }
  }

  # The best-scoring candidate can still lack numpy/sklearn (e.g. a bare system
  # python3 on Linux); auto-select it only if the bridge can actually run.
  preferred <- gp_preferred_python(purpose = "ml")
  if (!is.null(preferred) && nzchar(preferred)) {
    preferred_probe <- tryCatch(
      gp_python_probe(preferred, modules = c("numpy", "sklearn")),
      error = function(e) NULL
    )
    if (isTRUE(preferred_probe[["numpy"]]) && isTRUE(preferred_probe[["sklearn"]])) {
      return(normalizePath(preferred, winslash = "/", mustWork = TRUE))
    }
  }

  if (nzchar(generic) && file.exists(generic)) {
    return(generic)
  }
  NULL
}

gp_ml_run_cli <- function(command, args, python_bin = NULL) {
  python_bin <- python_bin %||% gp_detect_ml_python()
  if (is.null(python_bin) || !nzchar(python_bin)) {
    stop(
      "No configured Python runtime was found for Python ML models.\n",
      "Set PREDICTPRO_PYTHON to a Python with numpy, pandas and scikit-learn ",
      "(plus xgboost/catboost/lightgbm for those models).",
      call. = FALSE
    )
  }
  cmd_args <- c(gp_ml_bridge_script_path(), command, args)
  out <- system2(
    python_bin,
    args = gp_quote_system_args(cmd_args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) {
    stop(
      paste("ML Python bridge failed.", paste(out, collapse = "\n"), sep = "\n"),
      call. = FALSE
    )
  }
  invisible(out)
}

gp_ml_write_matrix_csv <- function(x, path) {
  gp_bridge_csv_write(as.data.frame(unname(as.matrix(x)), check.names = FALSE), path, quote = FALSE)
  invisible(path)
}

gp_ml_write_vector_csv <- function(x, path) {
  gp_bridge_csv_write(data.frame(value = x, stringsAsFactors = FALSE), path, quote = TRUE)
  invisible(path)
}

gp_ml_write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = TRUE), path, useBytes = TRUE)
  invisible(path)
}

gp_ml_read_json <- function(path, default = NULL) {
  if (is.null(path) || !file.exists(path)) {
    return(default)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

gp_ml_cli_backend <- function() {
  "cli"
}

gp_ml_temp_run_dir <- function(prefix = "predictpror_ml_bridge_") {
  path <- tempfile(pattern = prefix)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

gp_ml_read_prediction_payload <- function(out_dir) {
  pred_path <- file.path(out_dir, "predictions.csv")
  if (!file.exists(pred_path)) {
    stop("ML Python bridge did not write predictions.csv.", call. = FALSE)
  }
  pred_df <- utils::read.csv(pred_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!ncol(pred_df)) {
    stop("ML Python bridge wrote an empty predictions.csv.", call. = FALSE)
  }

  prob_path <- file.path(out_dir, "probabilities.csv")
  prob <- NULL
  if (file.exists(prob_path)) {
    prob_df <- utils::read.csv(prob_path, check.names = FALSE, stringsAsFactors = FALSE)
    prob <- as.matrix(prob_df)
    storage.mode(prob) <- "double"
  }

  meta <- gp_ml_read_json(file.path(out_dir, "meta.json"), default = list())
  classes <- meta$classes %||% NULL
  if (!is.null(classes) && !length(classes)) {
    classes <- NULL
  }

  list(
    predictions = unname(as.vector(pred_df[[1L]])),
    probabilities = prob,
    classes = classes,
    selected_alpha = meta$selected_alpha %||% NULL
  )
}

gp_py_ml_training_label_n <- function(y_train) {
  if (is.matrix(y_train) || is.data.frame(y_train)) {
    return(nrow(y_train))
  }
  length(y_train)
}

gp_py_ml_is_multitask_response_family <- function(response_family) {
  fam <- tolower(trimws(as.character(response_family %||% "")[1L]))
  fam %in% c("multitask", "multi_task", "multioutput", "multi_output", "multitask_regression")
}

gp_py_ml_nonmissing_training_rows <- function(y_train, response_family = "gaussian") {
  if (gp_py_ml_is_multitask_response_family(response_family)) {
    if (is.matrix(y_train) || is.data.frame(y_train)) {
      y_num <- as.data.frame(lapply(as.data.frame(y_train, stringsAsFactors = FALSE), function(col) suppressWarnings(as.numeric(col))))
      finite <- is.finite(as.matrix(y_num))
      return(rowSums(finite) > 0L)
    }
    return(is.finite(suppressWarnings(as.numeric(y_train))))
  }

  fam <- gp_resolve_response_family(response_family, y = y_train)
  if (is.matrix(y_train) || is.data.frame(y_train)) {
    y_df <- as.data.frame(y_train, stringsAsFactors = FALSE)
    if (identical(fam, "gaussian")) {
      y_num <- as.data.frame(lapply(y_df, function(col) suppressWarnings(as.numeric(col))))
      return(stats::complete.cases(y_num) & apply(as.matrix(y_num), 1L, function(row) all(is.finite(row))))
    }
    return(stats::complete.cases(y_df))
  }
  if (identical(fam, "gaussian")) {
    return(is.finite(suppressWarnings(as.numeric(y_train))))
  }
  !is.na(y_train)
}

gp_py_ml_filter_missing_training <- function(X_train,
                                             y_train,
                                             response_family = "gaussian",
                                             context = "ML Python bridge") {
  n_x <- nrow(X_train)
  n_y <- gp_py_ml_training_label_n(y_train)
  if (!identical(as.integer(n_x), as.integer(n_y))) {
    stop(
      sprintf(
        "%s received %d training predictor rows but %d training labels.",
        context,
        as.integer(n_x),
        as.integer(n_y)
      ),
      call. = FALSE
    )
  }

  keep <- gp_py_ml_nonmissing_training_rows(y_train, response_family = response_family)
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) {
    stop(
      sprintf("%s has no non-missing training labels after trait-specific missingness filtering.", context),
      call. = FALSE
    )
  }
  if (all(keep)) {
    return(list(X_train = X_train, y_train = y_train, dropped = 0L))
  }

  filtered_y <- if (is.matrix(y_train) || is.data.frame(y_train)) {
    y_train[keep, , drop = FALSE]
  } else {
    y_train[keep]
  }
  list(
    X_train = X_train[keep, , drop = FALSE],
    y_train = filtered_y,
    dropped = sum(!keep)
  )
}

gp_py_ml_normalize_model_params <- function(model_type, model_params = list()) {
  model_key <- tolower(as.character(model_type %||% ""))
  params <- Filter(Negate(is.null), model_params %||% list())
  if (identical(model_key, "randomforest")) {
    params$n_jobs <- as.integer(params$n_jobs %||% params$rf_n_jobs %||% 1L)
    params$rf_n_jobs <- NULL
  }
  # Thread-capable tree models default to 1 thread, which left the other
  # cores idle whenever the fit ran alone. A thread count of 1 (the default)
  # is raised to this process's share of the machine; values > 1 are kept.
  # PREDICTPRO_ML_AUTO_THREADS=false keeps the supplied values exactly.
  thread_key <- switch(
    model_key,
    randomforest = "n_jobs",
    xgboost = "xgb_nthread",
    catboost = "catboost_thread_count",
    lightgbm = "lightgbm_nthread",
    NULL
  )
  if (!is.null(thread_key) &&
      !tolower(Sys.getenv("PREDICTPRO_ML_AUTO_THREADS", "true")) %in% c("false", "0", "no")) {
    current <- suppressWarnings(as.integer(params[[thread_key]] %||% 1L))
    if (!is.finite(current) || current <= 1L) {
      params[[thread_key]] <- gp_ml_thread_share()
    }
  }
  params
}

# Threads one model fit may use: all logical cores but one in the main
# process, a worker's share (cores / workers) inside a parallel worker.
gp_ml_thread_share <- function() {
  cores <- tryCatch(gp_detect_cores(logical = TRUE, reserve = 1L), error = function(e) 1L)
  if (identical(Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER", ""), "1")) {
    n_workers <- suppressWarnings(as.integer(Sys.getenv("PREDICTPRO_PARALLEL_WORKERS", "")))
    if (!is.finite(n_workers) || n_workers < 1L) return(1L)
    return(max(1L, as.integer(floor(cores / n_workers))))
  }
  max(1L, as.integer(cores))
}

gp_ml_fit_predict_cli <- function(model_type,
                                  X_train,
                                  y_train,
                                  X_test,
                                  response_family = "gaussian",
                                  model_params = list(),
                                  class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML fit-predict", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  params <- gp_py_ml_normalize_model_params(model_type, model_params)
  if (identical(model_type, "xgboost") && is.null(params$device)) {
    params$device <- gp_py_ml_xgb_device(
      params,
      n_rows = nrow(X_train) + nrow(X_test),
      n_cols = ncol(X_train)
    )
  }

  x_train_payload <- unname(as.matrix(X_train))
  x_test_payload <- unname(as.matrix(X_test))
  y_payload <- gp_py_ml_prepare_target(y_train, fam)
  write_payload <- function(input_dir) {
    gp_ml_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
    gp_ml_write_vector_csv(y_payload, file.path(input_dir, "y_train.csv"))
    gp_ml_write_matrix_csv(x_test_payload, file.path(input_dir, "x_test.csv"))
    if (!is.null(class_levels)) {
      gp_ml_write_json(as.character(class_levels), file.path(input_dir, "class_levels.json"))
    }
    invisible(input_dir)
  }
  payload_cache <- gp_bridge_payload_cache_prepare(
    prefix = "ML",
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
    max_entries = 16L
  )
  if (is.null(payload_cache)) {
    run_dir <- gp_ml_temp_run_dir()
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
    input_dir <- run_dir
    out_dir <- file.path(run_dir, "out")
    write_payload(input_dir)
  } else {
    input_dir <- payload_cache$input_dir
    out_dir <- gp_ml_temp_run_dir(prefix = "predictpror_ml_bridge_out_")
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  x_train_csv <- file.path(input_dir, "x_train.csv")
  y_train_csv <- file.path(input_dir, "y_train.csv")
  x_test_csv <- file.path(input_dir, "x_test.csv")
  params_json <- file.path(out_dir, "params.json")
  levels_json <- file.path(input_dir, "class_levels.json")
  gp_ml_write_json(params, params_json)

  args <- c(
    "--model", model_type,
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

  gp_ml_run_cli("fit-predict", args)
  gp_ml_read_prediction_payload(out_dir)
}

gp_ml_fit_predict_batch_cli <- function(jobs) {
  if (is.null(jobs) || !length(jobs)) {
    return(list())
  }
  batch_dir <- gp_ml_temp_run_dir(prefix = "predictpror_ml_batch_")
  on.exit(unlink(batch_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_jobs <- vector("list", length(jobs))
  out_dirs <- character(length(jobs))

  for (i in seq_along(jobs)) {
    job <- jobs[[i]]
    if (is.function(job$build)) {
      # Materialise this job's matrices only now; they are released below once
      # its payload is on disk, so one job's copies exist at a time.
      built <- job$build()
      job$X_train <- built$X_train
      job$X_test <- built$X_test
      job$build <- NULL
      rm(built)
    }
    model_type <- tolower(as.character(job$model_type %||% job$model %||% "")[1L])
    if (!nzchar(model_type)) {
      stop("ML batch job ", i, " is missing `model_type`.", call. = FALSE)
    }
    fam <- gp_resolve_response_family(job$response_family %||% "gaussian", y = job$y_train)
    params <- gp_py_ml_normalize_model_params(model_type, job$model_params %||% job$params %||% list())
    X_train <- job$X_train
    X_test <- job$X_test %||% job$X_pred
    y_train <- job$y_train
    filtered <- gp_py_ml_filter_missing_training(
      X_train = X_train,
      y_train = y_train,
      response_family = fam,
      context = paste("ML batch fit-predict", model_type)
    )
    X_train <- filtered$X_train
    y_train <- filtered$y_train
    class_levels <- job$class_levels %||% gp_py_ml_class_levels(y_train, fam)
    if (identical(model_type, "xgboost") && is.null(params$device)) {
      params$device <- gp_py_ml_xgb_device(
        params,
        n_rows = nrow(X_train) + nrow(X_test),
        n_cols = ncol(X_train)
      )
    }

    x_train_payload <- unname(as.matrix(X_train))
    x_test_payload <- unname(as.matrix(X_test))
    y_payload <- gp_py_ml_prepare_target(y_train, fam)
    write_payload <- function(input_dir) {
      gp_ml_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
      gp_ml_write_vector_csv(y_payload, file.path(input_dir, "y_train.csv"))
      gp_ml_write_matrix_csv(x_test_payload, file.path(input_dir, "x_test.csv"))
      if (!is.null(class_levels)) {
        gp_ml_write_json(as.character(class_levels), file.path(input_dir, "class_levels.json"))
      }
      invisible(input_dir)
    }
    payload_cache <- gp_bridge_payload_cache_prepare(
      prefix = "ML",
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
      max_entries = 16L
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
    gp_ml_write_json(params, params_json)

    out_dirs[[i]] <- out_dir
    manifest_jobs[[i]] <- Filter(Negate(is.null), list(
      id = as.character(job$id %||% i),
      model = model_type,
      task = fam,
      x_train_csv = file.path(input_dir, "x_train.csv"),
      y_train_csv = file.path(input_dir, "y_train.csv"),
      x_test_csv = file.path(input_dir, "x_test.csv"),
      params_json = params_json,
      class_levels_json = if (!is.null(class_levels)) file.path(input_dir, "class_levels.json") else NULL,
      out_dir = out_dir
    ))
    rm(job, X_train, X_test, filtered, x_train_payload, x_test_payload, write_payload)
  }

  manifest_json <- file.path(batch_dir, "manifest.json")
  gp_ml_write_json(list(jobs = manifest_jobs), manifest_json)
  gp_ml_run_cli("batch-fit-predict", c("--manifest-json", manifest_json))
  lapply(out_dirs, gp_ml_read_prediction_payload)
}

gp_ml_bootstrap_cli <- function(model_type,
                                X_train,
                                y_train,
                                X_pred,
                                response_family = "gaussian",
                                model_params = list(),
                                class_levels = NULL,
                                n_bootstrap = 30L,
                                seed = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML bootstrap", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  params <- gp_py_ml_normalize_model_params(model_type, model_params)
  if (identical(model_type, "xgboost") && is.null(params$device)) {
    params$device <- gp_py_ml_xgb_device(
      params,
      n_rows = nrow(X_train) + nrow(X_pred),
      n_cols = ncol(X_train)
    )
  }

  x_train_payload <- unname(as.matrix(X_train))
  x_pred_payload <- unname(as.matrix(X_pred))
  y_payload <- gp_py_ml_prepare_target(y_train, fam)
  bootstrap_seed <- as.integer(seed %||% gp_ml_random_state())
  write_payload <- function(input_dir) {
    gp_ml_write_matrix_csv(x_train_payload, file.path(input_dir, "x_train.csv"))
    gp_ml_write_vector_csv(y_payload, file.path(input_dir, "y_train.csv"))
    gp_ml_write_matrix_csv(x_pred_payload, file.path(input_dir, "x_pred.csv"))
    if (!is.null(class_levels)) {
      gp_ml_write_json(as.character(class_levels), file.path(input_dir, "class_levels.json"))
    }
    invisible(input_dir)
  }
  payload_cache <- gp_bridge_payload_cache_prepare(
    prefix = "ML",
    purpose = paste("bootstrap", fam, sep = "::"),
    key_parts = list(
      response_family = fam,
      x_train = x_train_payload,
      y_train = y_payload,
      x_pred = x_pred_payload,
      class_levels = as.character(class_levels %||% character())
    ),
    cells = gp_bridge_payload_cache_cells(x_train_payload, y_payload, x_pred_payload),
    build_fun = write_payload,
    min_cells = 20000L,
    max_entries = 16L
  )
  if (is.null(payload_cache)) {
    run_dir <- gp_ml_temp_run_dir()
    on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
    input_dir <- run_dir
    out_dir <- file.path(run_dir, "out")
    write_payload(input_dir)
  } else {
    input_dir <- payload_cache$input_dir
    out_dir <- gp_ml_temp_run_dir(prefix = "predictpror_ml_bootstrap_out_")
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  x_train_csv <- file.path(input_dir, "x_train.csv")
  y_train_csv <- file.path(input_dir, "y_train.csv")
  x_pred_csv <- file.path(input_dir, "x_pred.csv")
  params_json <- file.path(out_dir, "params.json")
  levels_json <- file.path(input_dir, "class_levels.json")
  gp_ml_write_json(params, params_json)

  args <- c(
    "--model", model_type,
    "--task", fam,
    "--x-train-csv", x_train_csv,
    "--y-train-csv", y_train_csv,
    "--x-pred-csv", x_pred_csv,
    "--params-json", params_json,
    "--n-bootstrap", as.character(as.integer(n_bootstrap %||% 30L)),
    "--n-jobs", as.character(gp_ml_bootstrap_jobs(
      n_bootstrap %||% 30L,
      data_gb = 8 * (length(x_train_payload) + length(x_pred_payload)) / 1024^3
    )),
    "--seed", as.character(bootstrap_seed),
    "--out-dir", out_dir
  )
  if (!is.null(class_levels)) {
    args <- c(args, "--class-levels-json", levels_json)
  }

  gp_ml_run_cli("bootstrap-fit-predict", args)
  boot_path <- file.path(out_dir, "bootstrap.csv")
  if (!file.exists(boot_path)) {
    stop("ML Python bridge did not write bootstrap.csv.", call. = FALSE)
  }
  boot <- utils::read.csv(boot_path, check.names = FALSE, stringsAsFactors = FALSE)
  as.matrix(boot)
}

ensure_ml_pydeps <- function(prefer_gpu = TRUE) {
  gp_ml_run_cli("setup-deps", if (isTRUE(prefer_gpu)) c("--prefer-gpu") else character())
  invisible(TRUE)
}

gp_ml_feature_importance_cli <- function(model_type,
                                         X_train,
                                         y_train,
                                         response_family = "gaussian",
                                         model_params = list(),
                                         class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML feature importance", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  params <- gp_py_ml_normalize_model_params(model_type, model_params)
  run_dir <- gp_ml_temp_run_dir(prefix = "predictpror_ml_importance_")
  on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
  out_dir <- file.path(run_dir, "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  x_train_csv <- file.path(run_dir, "x_train.csv")
  y_train_csv <- file.path(run_dir, "y_train.csv")
  params_json <- file.path(out_dir, "params.json")
  levels_json <- file.path(run_dir, "class_levels.json")
  gp_ml_write_matrix_csv(unname(as.matrix(X_train)), x_train_csv)
  gp_ml_write_vector_csv(gp_py_ml_prepare_target(y_train, fam), y_train_csv)
  gp_ml_write_json(params, params_json)
  if (!is.null(class_levels)) {
    gp_ml_write_json(as.character(class_levels), levels_json)
  }

  args <- c(
    "--model", model_type,
    "--task", fam,
    "--x-train-csv", x_train_csv,
    "--y-train-csv", y_train_csv,
    "--params-json", params_json,
    "--out-dir", out_dir
  )
  if (!is.null(class_levels)) {
    args <- c(args, "--class-levels-json", levels_json)
  }

  gp_ml_run_cli("feature-importance", args)
  imp_path <- file.path(out_dir, "feature_importance.csv")
  if (!file.exists(imp_path)) {
    stop("ML Python bridge did not write feature_importance.csv.", call. = FALSE)
  }
  imp <- utils::read.csv(imp_path, check.names = FALSE, stringsAsFactors = FALSE)
  as.numeric(imp[[1L]])
}

gp_py_ml_feature_importance <- function(model_type,
                                        X_train,
                                        y_train,
                                        response_family = "gaussian",
                                        model_params = list(),
                                        class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML feature importance", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  class_levels <- class_levels %||% gp_py_ml_class_levels(y_train, fam)
  params <- gp_py_ml_normalize_model_params(model_type, model_params)

  if (identical(tolower(model_type), "lasso")) {   # fitted in R by glmnet
    return(gp_lasso_feature_importance(X_train, y_train, response_family = fam, model_params = params))
  }
  if (identical(gp_ml_cli_backend(), "reticulate")) {
    mod <- gp_get_ml_bridge_module()
    py_res <- mod$feature_importance_payload(
      model = model_type,
      task = fam,
      X_train = unname(as.matrix(X_train)),
      y_train = gp_py_ml_prepare_target(y_train, fam),
      params = params,
      class_levels = class_levels
    )
    return(as.numeric(reticulate::py_to_r(py_res$get("importance"))))
  }

  gp_ml_feature_importance_cli(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    model_params = params,
    class_levels = class_levels
  )
}

gp_py_ml_class_levels <- function(y, response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family, y = y)
  if (identical(fam, "gaussian")) {
    return(NULL)
  }
  if (is.factor(y) || is.ordered(y)) {
    return(as.character(levels(y)))
  }
  sort(unique(as.character(stats::na.omit(y))))
}

gp_py_ml_prepare_target <- function(y, response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family, y = y)
  if (identical(fam, "gaussian")) {
    return(as.numeric(y))
  }
  as.character(y)
}

gp_py_ml_xgb_device <- function(model_params = NULL,
                                n_rows = NULL,
                                n_cols = NULL) {
  params <- model_params %||% list()
  booster <- tolower(as.character(params$xgb_booster %||% params$booster %||% "gbtree")[[1]])
  if (identical(booster, "gblinear")) {
    return(NULL)
  }

  explicit_device <- params$device %||%
    params$xgb_device %||%
    Sys.getenv("PREDICTPRO_XGB_DEVICE", unset = "")
  if (!is.null(explicit_device) && nzchar(as.character(explicit_device[[1]]))) {
    dev <- tolower(as.character(explicit_device[[1]]))
    if (dev %in% c("cuda", "gpu")) {
      return("cuda")
    }
    if (identical(dev, "cpu")) {
      return(NULL)
    }
  }

  auto_gpu <- tolower(Sys.getenv("PREDICTPRO_XGB_AUTO_GPU", unset = "true"))
  if (auto_gpu %in% c("0", "false", "no", "off")) {
    return(NULL)
  }

  # Measured (8 CPU threads vs CUDA, 100 trees, depth 6): CPU is ~3x faster
  # at 231 x 8,068 (1.9M cells); the GPU only wins around 2,000 x 20,000
  # (40M cells). Per-process CUDA start-up adds a few seconds per fit.
  min_cells <- suppressWarnings(as.numeric(Sys.getenv("PREDICTPRO_XGB_GPU_MIN_CELLS", "2e7")))
  n_cells <- suppressWarnings(as.numeric(n_rows %||% NA_real_) * as.numeric(n_cols %||% NA_real_))
  if (is.finite(min_cells) && min_cells > 0 && is.finite(n_cells) && n_cells < min_cells) {
    return(NULL)
  }

  n_gpu <- tryCatch(detect_physical_num_gpus(), error = function(e) NA_integer_)
  if (!is.finite(n_gpu) || as.integer(n_gpu) <= 0L) {
    return(NULL)
  }
  gpu_usage <- tryCatch(get_gpu_usage(), error = function(e) NA_real_)
  if (isTRUE(gp_gpu_is_busy(gpu_usage))) {
    return(NULL)
  }
  "cuda"
}

gp_py_ml_fit_predict_raw_reticulate <- function(model_type,
                                                X_train,
                                                y_train,
                                                X_test,
                                                response_family = "gaussian",
                                                model_params = list(),
                                                class_levels = NULL,
                                                prefer_gpu = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML reticulate fit-predict", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  class_levels <- class_levels %||% gp_py_ml_class_levels(y_train, fam)
  params <- gp_py_ml_normalize_model_params(model_type, model_params)
  if (identical(model_type, "xgboost") && is.null(params$device)) {
    params$device <- gp_py_ml_xgb_device(
      params,
      n_rows = nrow(X_train) + nrow(X_test),
      n_cols = ncol(X_train)
    )
  }
  mod <- gp_get_ml_bridge_module()
  py_res <- mod$fit_predict_payload(
    model = model_type,
    task = fam,
    X_train = unname(as.matrix(X_train)),
    y_train = gp_py_ml_prepare_target(y_train, fam),
    X_test = unname(as.matrix(X_test)),
    params = params,
    class_levels = class_levels
  )
  pred_df <- reticulate::py_to_r(py_res$get("predictions"))
  prob_obj <- py_res$get("probabilities")
  prob_df <- if (reticulate::py_is_null_xptr(prob_obj) || is.null(prob_obj)) NULL else reticulate::py_to_r(prob_obj)
  classes_obj <- py_res$get("classes")
  classes <- if (reticulate::py_is_null_xptr(classes_obj) || is.null(classes_obj)) NULL else reticulate::py_to_r(classes_obj)
  rm(py_res)
  invisible(gc(verbose = FALSE))

  list(
    predictions = unname(as.vector(pred_df)),
    probabilities = if (is.null(prob_df)) NULL else as.matrix(prob_df),
    classes = classes %||% NULL
  )
}

gp_py_ml_fit_predict_raw <- function(model_type,
                                     X_train,
                                     y_train,
                                     X_test,
                                     response_family = "gaussian",
                                     model_params = list(),
                                     class_levels = NULL,
                                     prefer_gpu = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  model_params <- gp_py_ml_normalize_model_params(model_type, model_params)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML fit-predict", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  class_levels <- class_levels %||% gp_py_ml_class_levels(y_train, fam)
  if (identical(tolower(model_type), "lasso")) {   # fitted in R by glmnet
    return(gp_lasso_fit_predict(X_train, y_train, X_test, response_family = fam, model_params = model_params))
  }
  if (identical(gp_ml_cli_backend(), "reticulate")) {
    return(gp_py_ml_fit_predict_raw_reticulate(
      model_type = model_type,
      X_train = X_train,
      y_train = y_train,
      X_test = X_test,
      response_family = fam,
      model_params = model_params,
      class_levels = class_levels,
      prefer_gpu = prefer_gpu
    ))
  }
  gp_ml_fit_predict_cli(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    response_family = fam,
    model_params = model_params,
    class_levels = class_levels
  )
}

gp_py_ml_format_prediction <- function(result,
                                       response_family = "gaussian",
                                       class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family)
  preds <- result$predictions

  if (identical(fam, "gaussian")) {
    out <- as.numeric(preds)
    if (!is.null(result$selected_alpha)) attr(out, "selected_alpha") <- as.numeric(result$selected_alpha)
    return(out)
  }

  prob <- NULL
  if (!is.null(result$probabilities)) {
    prob_raw <- result$probabilities
    if (is.list(prob_raw) && !is.data.frame(prob_raw)) {
      prob_raw <- do.call(rbind, lapply(prob_raw, function(x) unlist(x, use.names = FALSE)))
    }
    prob <- as.matrix(prob_raw)
    storage.mode(prob) <- "double"
    prob_levels <- as.character(result$classes %||% class_levels %||% colnames(prob))
    if (length(prob_levels) != ncol(prob)) {
      if (!is.null(class_levels) && length(class_levels) == ncol(prob)) {
        prob_levels <- as.character(class_levels)
      } else if (ncol(prob) == 1L && !is.null(class_levels) && length(class_levels) >= 1L) {
        prob_levels <- as.character(class_levels[[min(2L, length(class_levels))]])
      } else {
        prob_levels <- paste0("class", seq_len(ncol(prob)))
      }
    }
    colnames(prob) <- prob_levels
    if (!is.null(class_levels) && length(class_levels) > 0L && !identical(fam, "binary")) {
      missing_levels <- setdiff(as.character(class_levels), colnames(prob))
      if (length(missing_levels)) {
        missing_mat <- matrix(0, nrow = nrow(prob), ncol = length(missing_levels))
        colnames(missing_mat) <- missing_levels
        prob <- cbind(prob, missing_mat)
      }
      prob <- prob[, as.character(class_levels), drop = FALSE]
    }
  }

  if (identical(fam, "binary")) {
    pred <- gp_ml_binary_prob_column(prob, class_levels = class_levels %||% colnames(prob))
    if (!is.null(prob)) {
      attr(pred, "probabilities") <- prob
    }
    return(pred)
  }

  pred <- as.character(preds)
  if (!is.null(prob)) {
    attr(pred, "probabilities") <- prob
  }
  pred
}

gp_py_ml_fit_predict <- function(model_type,
                                 X_train,
                                 y_train,
                                 X_test,
                                 response_family = "gaussian",
                                 model_params = list(),
                                 class_levels = NULL,
                                 prefer_gpu = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  class_levels <- class_levels %||% gp_py_ml_class_levels(y_train, fam)
  res <- gp_py_ml_fit_predict_raw(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    response_family = fam,
    model_params = model_params,
    class_levels = class_levels,
    prefer_gpu = prefer_gpu
  )
  gp_py_ml_format_prediction(
    result = res,
    response_family = fam,
    class_levels = class_levels
  )
}

gp_py_ml_fit_predict_batch <- function(jobs) {
  if (is.null(jobs) || !length(jobs)) {
    return(list())
  }
  prepared <- lapply(jobs, function(job) {
    fam <- gp_resolve_response_family(job$response_family %||% "gaussian", y = job$y_train)
    # Lazy jobs (`build()` returns X_train/X_test) are materialised and
    # filtered one at a time by gp_ml_fit_predict_batch_cli().
    if (!is.function(job$build)) {
      filtered <- gp_py_ml_filter_missing_training(
        X_train = job$X_train,
        y_train = job$y_train,
        response_family = fam,
        context = paste("ML batch fit-predict", job$model_type %||% job$model %||% "")
      )
      job$X_train <- filtered$X_train
      job$y_train <- filtered$y_train
    }
    job$response_family <- fam
    job$class_levels <- job$class_levels %||% gp_py_ml_class_levels(job$y_train, fam)
    job$model_params <- gp_py_ml_normalize_model_params(job$model_type %||% job$model, job$model_params %||% job$params %||% list())
    job
  })
  # Lasso jobs are fitted in R by glmnet (one at a time, built lazily); the
  # rest go to the Python bridge in one call.
  is_lasso <- vapply(prepared, function(job) {
    identical(tolower(as.character(job$model_type %||% job$model %||% "")[1L]), "lasso")
  }, logical(1L))
  raw <- vector("list", length(prepared))
  for (i in which(is_lasso)) {
    job <- prepared[[i]]
    if (is.function(job$build)) {
      built <- job$build()
      job$X_train <- built$X_train
      job$X_test <- built$X_test
      rm(built)
    }
    filtered <- gp_py_ml_filter_missing_training(job$X_train, job$y_train, response_family = job$response_family,
                                                 context = "ML batch fit-predict lasso")
    raw[[i]] <- gp_lasso_fit_predict(filtered$X_train, filtered$y_train, job$X_test,
                                     response_family = job$response_family, model_params = job$model_params)
    rm(job, filtered)
  }
  if (any(!is_lasso)) {
    raw[!is_lasso] <- gp_ml_fit_predict_batch_cli(prepared[!is_lasso])
  }
  Map(function(result, job) {
    gp_py_ml_format_prediction(
      result = result,
      response_family = job$response_family,
      class_levels = job$class_levels
    )
  }, raw, prepared)
}

gp_py_ml_bootstrap_predict <- function(data_label_geno,
                                       indices,
                                       test_geno = NULL,
                                       model_type,
                                       response_family = "gaussian",
                                       model_params = list(),
                                       class_levels = NULL) {
  fam <- gp_resolve_response_family(response_family, y = data_label_geno[, 1])
  train_x <- data_label_geno[indices, -1, drop = FALSE]
  train_y <- data_label_geno[indices, 1]
  pred_x <- test_geno %||% data_label_geno[, -1, drop = FALSE]

  res <- gp_py_ml_fit_predict_raw(
    model_type = model_type,
    X_train = train_x,
    y_train = train_y,
    X_test = pred_x,
    response_family = fam,
    model_params = model_params,
    prefer_gpu = identical(model_type, "xgboost")
  )

  if (identical(fam, "multiclass")) {
    prob <- as.matrix(res$probabilities)
    prob_levels <- as.character(res$classes %||% class_levels %||% colnames(prob))
    if (length(prob_levels) != ncol(prob)) {
      if (!is.null(class_levels) && length(class_levels) == ncol(prob)) {
        prob_levels <- as.character(class_levels)
      } else {
        prob_levels <- paste0("class", seq_len(ncol(prob)))
      }
    }
    colnames(prob) <- prob_levels
    target_levels <- class_levels %||% prob_levels
    missing_levels <- setdiff(as.character(target_levels), colnames(prob))
    if (length(missing_levels)) {
      missing_mat <- matrix(0, nrow = nrow(prob), ncol = length(missing_levels))
      colnames(missing_mat) <- missing_levels
      prob <- cbind(prob, missing_mat)
    }
    prob <- prob[, target_levels, drop = FALSE]
    return(as.numeric(t(prob)))
  }

  formatted <- gp_py_ml_format_prediction(
    result = res,
    response_family = fam,
    class_levels = class_levels
  )
  as.numeric(formatted)
}

gp_py_ml_bootstrap_matrix_reticulate <- function(model_type,
                                                 X_train,
                                                 y_train,
                                                 X_pred,
                                                 response_family = "gaussian",
                                                 model_params = list(),
                                                 class_levels = NULL,
                                                 n_bootstrap = 30L,
                                                 seed = NULL,
                                                 prefer_gpu = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML reticulate bootstrap", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  params <- gp_py_ml_normalize_model_params(model_type, model_params)
  if (identical(model_type, "xgboost") && is.null(params$device)) {
    params$device <- gp_py_ml_xgb_device(
      params,
      n_rows = nrow(X_train) + nrow(X_pred),
      n_cols = ncol(X_train)
    )
  }
  mod <- gp_get_ml_bridge_module()
  py_res <- mod$bootstrap_fit_predict_payload(
    model = model_type,
    task = fam,
    X_train = unname(as.matrix(X_train)),
    y_train = gp_py_ml_prepare_target(y_train, fam),
    X_pred = unname(as.matrix(X_pred)),
    n_bootstrap = as.integer(n_bootstrap %||% 30L),
    seed = as.integer(seed %||% gp_ml_random_state()),
    params = params,
    class_levels = class_levels
  )
  boot <- reticulate::py_to_r(py_res$get("bootstrap"))
  rm(py_res)
  invisible(gc(verbose = FALSE))
  as.matrix(boot)
}

gp_py_ml_bootstrap_matrix <- function(model_type,
                                      X_train,
                                      y_train,
                                      X_pred,
                                      response_family = "gaussian",
                                      model_params = list(),
                                      class_levels = NULL,
                                      n_bootstrap = 30L,
                                     seed = NULL,
                                     prefer_gpu = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_train)
  model_params <- gp_py_ml_normalize_model_params(model_type, model_params)
  filtered <- gp_py_ml_filter_missing_training(
    X_train = X_train,
    y_train = y_train,
    response_family = fam,
    context = paste("ML bootstrap", model_type)
  )
  X_train <- filtered$X_train
  y_train <- filtered$y_train
  if (identical(tolower(model_type), "lasso")) {   # fitted in R by glmnet
    return(gp_lasso_bootstrap_matrix(X_train, y_train, X_pred, response_family = fam,
                                     model_params = model_params, n_bootstrap = n_bootstrap, seed = seed))
  }
  if (identical(gp_ml_cli_backend(), "reticulate")) {
    return(gp_py_ml_bootstrap_matrix_reticulate(
      model_type = model_type,
      X_train = X_train,
      y_train = y_train,
      X_pred = X_pred,
      response_family = fam,
      model_params = model_params,
      class_levels = class_levels,
      n_bootstrap = n_bootstrap,
      seed = seed,
      prefer_gpu = prefer_gpu
    ))
  }
  gp_ml_bootstrap_cli(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_pred = X_pred,
    response_family = fam,
    model_params = model_params,
    class_levels = class_levels,
    n_bootstrap = n_bootstrap,
    seed = seed
  )
}
