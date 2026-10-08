# Policy decisions below must not depend on this machine's free RAM.
withr::local_envvar(c(GP_PAR_MEMORY_BUDGET_GB = "16"))   # reset when this file finishes

test_that("python runtime pickup sets env var without eagerly initializing reticulate", {
  skip_if_not_installed("reticulate")

  old_reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
  old_predictpro_python <- Sys.getenv("PREDICTPRO_PYTHON", unset = NA_character_)
  on.exit({
    if (is.na(old_reticulate_python)) Sys.unsetenv("RETICULATE_PYTHON") else Sys.setenv(RETICULATE_PYTHON = old_reticulate_python)
    if (is.na(old_predictpro_python)) Sys.unsetenv("PREDICTPRO_PYTHON") else Sys.setenv(PREDICTPRO_PYTHON = old_predictpro_python)
  }, add = TRUE)

  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.unsetenv("PREDICTPRO_PYTHON")

  preferred <- PredictProR:::gp_preferred_python()
  skip_if(is.null(preferred), "No preferred Python candidate found on this machine")

  local_mocked_bindings(
    gp_conda_python_candidates = function(env_names) character(),
    .package = "PredictProR"
  )

  was_initialized <- reticulate::py_available(initialize = FALSE)
  configured <- PredictProR:::gp_configure_python_runtime(force = TRUE)
  expect_identical(configured, preferred)
  expect_identical(Sys.getenv("RETICULATE_PYTHON"), preferred)
  expect_identical(reticulate::py_available(initialize = FALSE), was_initialized)
})

test_that("DL subprocess environment always supplies deterministic cuBLAS configuration", {
  old_value <- Sys.getenv("CUBLAS_WORKSPACE_CONFIG", unset = NA_character_)
  on.exit({
    if (is.na(old_value)) {
      Sys.unsetenv("CUBLAS_WORKSPACE_CONFIG")
    } else {
      Sys.setenv(CUBLAS_WORKSPACE_CONFIG = old_value)
    }
  }, add = TRUE)

  Sys.unsetenv("CUBLAS_WORKSPACE_CONFIG")
  expect_identical(
    PredictProR:::gp_dl_subprocess_env(),
    "CUBLAS_WORKSPACE_CONFIG=:4096:8"
  )

  Sys.setenv(CUBLAS_WORKSPACE_CONFIG = ":16:8")
  expect_identical(
    PredictProR:::gp_dl_subprocess_env(),
    "CUBLAS_WORKSPACE_CONFIG=:16:8"
  )
})

test_that("preferred python selection is workload-aware for DL", {
  rm(
    list = ls(envir = PredictProR:::PredictProR_runtime_cache),
    envir = PredictProR:::PredictProR_runtime_cache
  )

  # Candidates must exist on disk; use stand-in files so this runs on any OS.
  root <- withr::local_tempdir()
  general_py <- file.path(root, "Python", "python.exe")
  dl_py <- file.path(root, ".virtualenvs", "r-torchgwas", "Scripts", "python.exe")
  for (p in c(general_py, dl_py)) {
    dir.create(dirname(p), recursive = TRUE)
    file.create(p)
  }
  general_py <- normalizePath(general_py, winslash = "/")
  dl_py <- normalizePath(dl_py, winslash = "/")

  local_mocked_bindings(
    gp_python_candidates = function() c(general_py, dl_py),
    gp_python_module_score = function(py_bin, purpose = c("general", "ml", "dl"), modules = c("numpy", "sklearn", "xgboost", "lightgbm", "catboost", "torch")) {
      purpose <- match.arg(purpose)
      if (grepl("r-torchgwas", py_bin, fixed = TRUE)) {
        return(if (identical(purpose, "dl")) 30 else 5)
      }
      return(if (identical(purpose, "dl")) 10 else 12)
    },
    .package = "PredictProR"
  )

  expect_identical(PredictProR:::gp_preferred_python(purpose = "general"), general_py)
  expect_identical(PredictProR:::gp_preferred_python(purpose = "dl"), dl_py)
  rm(
    list = ls(envir = PredictProR:::PredictProR_runtime_cache),
    envir = PredictProR:::PredictProR_runtime_cache
  )
})

test_that("ML python detection falls back from generic DL runtime without sklearn", {
  ml_py <- tempfile("ml-python-", fileext = ".exe")
  dl_py <- tempfile("dl-python-", fileext = ".exe")
  file.create(ml_py)
  file.create(dl_py)

  withr::local_envvar(c(
    PREDICTPRO_PYTHON = normalizePath(dl_py, winslash = "/", mustWork = TRUE),
    PREDICTPRO_ML_PYTHON = NA
  ))

  local_mocked_bindings(
    gp_preferred_python = function(purpose = c("general", "ml", "dl", "gp")) {
      purpose <- match.arg(purpose)
      if (identical(purpose, "ml")) {
        return(normalizePath(ml_py, winslash = "/", mustWork = TRUE))
      }
      NULL
    },
    gp_python_probe = function(py_bin, modules = c("numpy", "sklearn")) {
      list(numpy = TRUE, sklearn = identical(
        normalizePath(py_bin, winslash = "/", mustWork = TRUE),
        normalizePath(ml_py, winslash = "/", mustWork = TRUE)
      ))
    },
    .package = "PredictProR"
  )

  actual <- PredictProR:::gp_detect_ml_python()
  expected <- normalizePath(ml_py, winslash = "/", mustWork = TRUE)
  if (!identical(actual, expected)) {
    stop(
      "Expected ML Python '", expected, "' but detected '", actual, "'.",
      call. = FALSE
    )
  }
  expect_true(TRUE)
})

test_that("ML python detection rejects an auto-found interpreter without numpy/sklearn", {
  rm(
    list = ls(envir = PredictProR:::PredictProR_runtime_cache),
    envir = PredictProR:::PredictProR_runtime_cache
  )
  bare_py <- tempfile("bare-python-", fileext = ".exe")
  file.create(bare_py)
  withr::local_envvar(c(PREDICTPRO_PYTHON = NA, PREDICTPRO_ML_PYTHON = NA, PREDICTPRO_PYTHON_ML = NA))
  local_mocked_bindings(
    gp_preferred_python = function(purpose = c("general", "ml", "dl", "gp")) {
      normalizePath(bare_py, winslash = "/", mustWork = TRUE)
    },
    gp_python_probe = function(py_bin, modules = c("numpy", "sklearn")) {
      list(numpy = FALSE, sklearn = FALSE)
    },
    .package = "PredictProR"
  )
  expect_null(PredictProR:::gp_detect_ml_python_uncached())
  expect_error(
    PredictProR:::gp_ml_run_cli("fit-predict", character()),
    "numpy, pandas and scikit-learn"
  )
})

test_that("XGBoost CLI device policy avoids reticulate initialization", {
  skip_if_not_installed("reticulate")

  old_device <- Sys.getenv("PREDICTPRO_XGB_DEVICE", unset = NA_character_)
  old_auto <- Sys.getenv("PREDICTPRO_XGB_AUTO_GPU", unset = NA_character_)
  on.exit({
    if (is.na(old_device)) Sys.unsetenv("PREDICTPRO_XGB_DEVICE") else Sys.setenv(PREDICTPRO_XGB_DEVICE = old_device)
    if (is.na(old_auto)) Sys.unsetenv("PREDICTPRO_XGB_AUTO_GPU") else Sys.setenv(PREDICTPRO_XGB_AUTO_GPU = old_auto)
  }, add = TRUE)
  Sys.unsetenv("PREDICTPRO_XGB_DEVICE")
  Sys.setenv(PREDICTPRO_XGB_AUTO_GPU = "true")

  local_mocked_bindings(
    detect_physical_num_gpus = function() 1L,
    get_gpu_usage = function() 5,
    .package = "PredictProR"
  )

  was_initialized <- reticulate::py_available(initialize = FALSE)
  expect_null(PredictProR:::gp_py_ml_xgb_device(list(xgb_booster = "gbtree"), n_rows = 30L, n_cols = 8L))
  # CPU is faster below ~20M cells; GPU only for large matrices
  expect_null(PredictProR:::gp_py_ml_xgb_device(list(xgb_booster = "gbtree"), n_rows = 5000L, n_cols = 500L))
  expect_identical(PredictProR:::gp_py_ml_xgb_device(list(xgb_booster = "gbtree"), n_rows = 5000L, n_cols = 5000L), "cuda")
  expect_identical(reticulate::py_available(initialize = FALSE), was_initialized)
})

test_that("deep-model policy forces sequential execution on a single GPU", {
  local_mocked_bindings(
    detect_num_gpus = function() 1L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("mlp", "RandomForest"),
    model_params_list = list(list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:10),
      omics_data = matrix(1, nrow = 2, ncol = 2),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )

  expect_identical(res$backend, "sequential")
  expect_identical(res$internal_flags, c(TRUE, TRUE))
  expect_identical(res$num_gpus, 1L)
  expect_identical(res$decision_reason, "single_gpu_exclusive_deep_models")
})

test_that("deep-model policy allows parallel backend when multiple GPUs are present", {
  local_mocked_bindings(
    detect_num_gpus = function() 2L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("mlp", "tabnet"),
    model_params_list = list(list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:10),
      omics_data = matrix(1, nrow = 2, ncol = 2),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 4L
  )

  expect_identical(res$backend, "mirai")
  expect_identical(res$internal_flags, c(FALSE, FALSE))
  expect_identical(res$num_gpus, 2L)
  expect_true(res$workers >= 1L)
})

test_that("policy forces Xgboost with internal threads and RKHS into sequential execution", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 4L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("Xgboost", "RKHS", "RandomForest"),
    model_params_list = list(
      list(xgb_nthread = 4L),
      list(),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 3L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 3L
  )

  expect_identical(res$internal_flags, c(TRUE, TRUE, TRUE))
  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "multiple_sequential_constraints")
})

test_that("policy forces CatBoost and LightGBM with internal threads into sequential execution", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("CatBoost", "LightGBM", "RandomForest"),
    model_params_list = list(
      list(catboost_thread_count = 4L),
      list(lightgbm_nthread = 3L),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 3L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 3L
  )

  expect_identical(res$internal_flags, c(TRUE, TRUE, TRUE))
  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "multiple_sequential_constraints")
})

test_that("cross-validation orchestration passes overload-sensitive models into policy decision", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  omics <- matrix(
    1:16,
    nrow = 4,
    dimnames = list(pheno$GID, paste0("m", 1:4))
  )

  seen_decision_args <- list()

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, sequential_models = NULL, decision_args = NULL, ...) {
      out <- lapply(seq_along(model_ids), run_task)
      seen_decision_args[[length(seen_decision_args) + 1L]] <<- list(
        model_ids = model_ids,
        sequential_models = sequential_models,
        decision_args = decision_args
      )
      attr(out, "gp_execution_policy") <- list(
        backend = "sequential",
        workers = 1L,
        plan = NULL,
        chunk_size = NULL,
        num_gpus = 0L,
        model_ids = model_ids,
        user_sequential_models = sequential_models,
        internal_flags = rep(TRUE, length(model_ids)),
        internal_sequential_count = length(model_ids),
        parallel_count = 0L
      )
      out
    },
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, omics_data = NULL, tst = NULL, additional_params = NULL) {
      rep(mean(y[-tst], na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = "Yield",
      gen_name = "GID",
      test_size = 0.25,
      random_state = 1L,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = "mean_squared_error",
      GS_model_cv = c("Xgboost", "RKHS", "mlp"),
      ml_dat_res = list(
        merged_data = list(merge_data = omics),
        omic_count = 1L
      ),
      model_prep_all_bayes_cv = list(dummy = TRUE),
      xgb_nthread = 4L,
      parallel_mode = "auto"
    ),
    verbose = FALSE
  )

  expect_true(length(res) >= 3L)
  expect_length(seen_decision_args, 1L)
  expect_equal(sort(seen_decision_args[[1L]]$model_ids), sort(c("Xgboost", "RKHS", "mlp")))
  expect_true(is.null(seen_decision_args[[1L]]$sequential_models) || !("mlp" %in% seen_decision_args[[1L]]$sequential_models))
  expect_identical(seen_decision_args[[1L]]$decision_args$model_params_list[[1L]]$xgb_nthread, 4L)
  expect_null(seen_decision_args[[1L]]$decision_args$model_params_list[[2L]]$xgb_nthread)
  expect_null(seen_decision_args[[1L]]$decision_args$model_params_list[[3L]]$xgb_nthread)
  expect_identical(seen_decision_args[[1L]]$decision_args$user_mode, "auto")
  expect_identical(seen_decision_args[[1L]]$decision_args$num_cores, NULL)
  expect_identical(attr(res, "gp_execution_policy")$backend, "sequential")
})

test_that("best-model orchestration skips Python initializer for direct CLI models", {
  skip_if_not_installed("reticulate")

  old_reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = NA_character_)
  old_predictpro_python <- Sys.getenv("PREDICTPRO_PYTHON", unset = NA_character_)
  on.exit({
    if (is.na(old_reticulate_python)) Sys.unsetenv("RETICULATE_PYTHON") else Sys.setenv(RETICULATE_PYTHON = old_reticulate_python)
    if (is.na(old_predictpro_python)) Sys.unsetenv("PREDICTPRO_PYTHON") else Sys.setenv(PREDICTPRO_PYTHON = old_predictpro_python)
  }, add = TRUE)

  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.unsetenv("PREDICTPRO_PYTHON")

  preferred <- PredictProR:::gp_preferred_python()
  skip_if(is.null(preferred), "No preferred Python candidate found on this machine")

  seen <- NULL

  local_mocked_bindings(
    gp_execute_partitioned_tasks = function(model_ids, run_task, sequential_models = NULL, decision_args = NULL, initializer = NULL, ...) {
      seen <<- list(
        model_ids = model_ids,
        sequential_models = sequential_models,
        decision_args = decision_args,
        initializer = initializer
      )
      list(list(ok = TRUE))
    },
    .package = "PredictProR"
  )

  was_initialized <- reticulate::py_available(initialize = FALSE)

  out <- PredictProR:::gp_execute_best_model_tasks(
    best_models = data.frame(model = "mlp", trait = "Yield", stringsAsFactors = FALSE),
    run_one_best = function(i) list(id = i),
    sequential_models = "mlp",
    canonical_names = c("mlp", "tabnet"),
    friendly_names = c("DenseNeuralNet", "TabNet"),
    parallel_mode = "auto",
    num_cores = 2L,
    globals_max_GB = 1,
    verbose = FALSE,
    parallel_backend_prefer_fork = FALSE,
    pheno_clean = list(pheno_clean_data = data.frame(GID = c("g1", "g2"), Yield = c(1, 2))),
    ml_dat_res = list(),
    gmatrix_kernel_model_ready_list = list(),
    geno_omic_model_ready_list = list(),
    response = "Yield",
    init_py = function() PredictProR:::gp_init_python_once(PredictProR:::gp_detect_python())
  )

  expect_length(out, 1L)
  expect_identical(seen$model_ids, "mlp")
  expect_identical(seen$sequential_models, "mlp")
  expect_null(seen$initializer)
  expect_identical(Sys.getenv("RETICULATE_PYTHON", unset = ""), "")
  expect_identical(reticulate::py_available(initialize = FALSE), was_initialized)
})

test_that("partitioned task initializer runs for user-sequential tasks without pinning BLAS in the main process", {
  # BLAS/OpenMP pinning is for parallel workers only: pinning the main
  # process leaked OMP_NUM_THREADS = 1 into the user's session.
  withr::local_envvar(OMP_NUM_THREADS = NA)
  seen <- character()

  local_mocked_bindings(
    init_single_thread_blas = function() {
      Sys.setenv(OMP_NUM_THREADS = "1")
      seen <<- c(seen, "blas")
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_execute_partitioned_tasks(
    model_ids = c("mlp"),
    run_task = function(i) {
      seen <<- c(seen, paste0("task", i, ":", Sys.getenv("OMP_NUM_THREADS", unset = "")))
      list(i = i)
    },
    sequential_models = c("mlp"),
    initializer = function() {
      seen <<- c(seen, "user")
    }
  )

  expect_length(out, 1L)
  expect_identical(seen, c("user", "task1:"))
  expect_identical(Sys.getenv("OMP_NUM_THREADS", unset = ""), "")
})

test_that("user-sequential execution policy keeps detected GPU count", {
  local_mocked_bindings(
    detect_num_gpus = function() 1L,
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_execute_partitioned_tasks(
    model_ids = c("mlp"),
    run_task = function(i) list(i = i),
    sequential_models = c("mlp")
  )

  expect_identical(attr(out, "gp_execution_policy")$backend, "sequential")
  expect_identical(attr(out, "gp_execution_policy")$num_gpus, 1L)
})

test_that("partitioned execution restores original task order", {
  local_mocked_bindings(
    sp_decide_policy = function(...) {
      list(
        backend = "future",
        workers = 2L,
        plan = "multisession",
        chunk_size = 1L,
        backend_score = 1,
        decision_reason = "test_mixed_partitions",
        num_gpus = 0L,
        internal_flags = c(FALSE, TRUE, FALSE, TRUE)
      )
    },
    gp_task_device_routes = function(model_ids, decision, task_indices) {
      rep(NA_character_, length(model_ids))
    },
    sp_apply = function(X, FUN, ...) lapply(X, FUN),
    .package = "PredictProR"
  )

  model_ids <- c("RandomForest", "RKHS", "Xgboost", "GBLUP")
  out <- PredictProR:::gp_execute_partitioned_tasks(
    model_ids = model_ids,
    run_task = function(i) list(task_index = i, model = model_ids[[i]]),
    decision_args = list()
  )

  expect_identical(
    vapply(out, `[[`, integer(1L), "task_index"),
    seq_along(model_ids)
  )
  expect_identical(
    vapply(out, `[[`, character(1L), "model"),
    model_ids
  )
  expect_identical(
    attr(out, "gp_execution_policy")[["internal_flags"]],
    c(FALSE, TRUE, FALSE, TRUE)
  )
})

test_that("GPU detection is quiet when torch is unavailable", {
  skip_if_not_installed("reticulate")

  local_mocked_bindings(
    gp_detect_python = function() "python",
    gp_init_python_once = function(py_bin) invisible(TRUE),
    detect_physical_num_gpus = function() NA_integer_,
    .package = "PredictProR"
  )
  local_mocked_bindings(
    py_module_available = function(module) FALSE,
    import = function(...) stop("ModuleNotFoundError: No module named 'torch'"),
    .package = "reticulate"
  )

  expect_silent({
    n <- PredictProR:::detect_num_gpus()
    expect_identical(n, 0L)
  })
})

test_that("python runtime status reports configured-but-uninitialized runtime without probing torch", {
  skip_if_not_installed("reticulate")

  withr::local_envvar(c(
    RETICULATE_PYTHON = NA_character_,
    PREDICTPRO_PYTHON = NA_character_
  ))

  local_mocked_bindings(
    gp_preferred_python = function(...) "C:/python.exe",
    gp_detect_python = function() "C:/python.exe",
    .package = "PredictProR"
  )
  local_mocked_bindings(
    py_available = function(initialize = FALSE) FALSE,
    py_module_available = function(module) stop("should not probe modules when initialize = FALSE"),
    import = function(...) stop("should not import modules when initialize = FALSE"),
    .package = "reticulate"
  )

  st <- PredictProR:::gp_python_runtime_status(initialize = FALSE)
  expect_identical(st$preferred_python, "C:/python.exe")
  expect_identical(st$active_python, "C:/python.exe")
  expect_false(st$reticulate_initialized)
  expect_true(st$python_configured)
  expect_false(st$torch_available)
  expect_identical(st$num_gpus, 0L)
})

test_that("python runtime status honors dedicated DL Python on servers", {
  skip_if_not_installed("reticulate")

  withr::local_envvar(c(
    PREDICTPRO_DL_PYTHON = "C:/server/predictdl/python.exe",
    RETICULATE_PYTHON = NA_character_,
    PREDICTPRO_PYTHON = NA_character_
  ))

  local_mocked_bindings(
    gp_preferred_python = function(...) NULL,
    gp_detect_python = function() stop("generic Python detection should not override DL Python"),
    .package = "PredictProR"
  )
  local_mocked_bindings(
    py_available = function(initialize = FALSE) FALSE,
    py_module_available = function(module) stop("should not probe modules when initialize = FALSE"),
    import = function(...) stop("should not import modules when initialize = FALSE"),
    .package = "reticulate"
  )

  st <- PredictProR:::gp_python_runtime_status(initialize = FALSE, purpose = "dl")
  expect_identical(st$active_python, "C:/server/predictdl/python.exe")
  expect_true(st$python_configured)
  expect_false(st$reticulate_initialized)
})

test_that("runtime metadata captures thread policy fields", {
  meta <- PredictProR:::gp_runtime_metadata(
      execution_policy = list(
        backend = "sequential",
        workers = 1L,
        plan = NULL,
        chunk_size = NULL,
        num_gpus = 0L,
        decision_reason = "unit_test",
        backend_score = 1.23,
        user_sequential_models = "mlp",
        internal_sequential_count = 1L,
        parallel_count = 0L,
        model_ids = "mlp"
    ),
    python_path = "python",
    preferred_python = "python",
    context = "cross_validation"
  )

    expect_true(all(c(
      "policy_decision_reason",
      "policy_backend_score",
      "policy_hidden_threads",
      "policy_blas_threads_env",
      "reticulate_available",
    "python_configured",
    "python_exists",
    "python_initializable",
    "torch_available",
    "cuda_available",
    "runtime_num_gpus",
    "runtime_gpu_name"
  ) %in% meta$key))
})

test_that("model_execute does not call setup_predictdl_env when deep runtime is already configured", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  omics <- matrix(
    1:16,
    nrow = 4,
    dimnames = list(pheno$GID, paste0("m", 1:4))
  )

  local_mocked_bindings(
    gp_detect_python = function() "python",
    setup_predictdl_env = function(...) stop("setup_predictdl_env should not be called"),
    phenotype_to_model = function(...) {
      ph <- pheno
      attr(ph, "cleared") <- "for_model_fit"
      list(
        pheno_clean_data = ph,
        test_set = NULL
      )
    },
    gp_prepare_model_input_objects = function(ctx) {
      ph <- pheno
      attr(ph, "cleared") <- "for_model_fit"
      list(
        low_call_rate_inds_removed = NULL,
        pheno_clean = list(pheno_clean_data = ph, test_set = NULL),
        geno_res = list(),
        omic1_res = list(),
        omic2_res = list(),
        omic3_res = list(),
        geno_omic_model_ready_list = list(),
        gmatrix_kernel_model_ready_list = list(),
        test_set = NULL,
        ml_dat_res = list(
          pheno_clean_data = ph,
          merged_data = list(merge_data = omics),
          omic_count = 1L
        )
      )
    },
    gp_prepare_cross_validation_artifacts = function(ctx) {
      list(
        pheno_data = pheno,
        model_prep_all_bayes_cv = NULL,
        asreml_models_prep_cv = NULL
      )
    },
    gp_run_cross_validation_pipeline = function(ctx) {
      list(
        cv_results = list(),
        cv_results_processed = list(
          best_models_list = list(mean_squared_error = data.frame(trait = "Yield", model = "DenseNeuralNet", mean_squared_error = 1)),
          plot_reps_list = list(mean_squared_error = list(ggplot_boxplot_reps = NULL)),
          plot_mean_list = list(mean_squared_error = list(ggplot_lineplot_mean = NULL)),
          run_metadata = data.frame(key = "context", value = "cross_validation", stringsAsFactors = FALSE)
        ),
        cv_results_predicted_vs_observed = NULL,
        run_metadata = data.frame(key = "context", value = "cross_validation", stringsAsFactors = FALSE),
        best_models = data.frame(trait = "Yield", model = "DenseNeuralNet", mean_squared_error = 1),
        best_models_ggplot_rep = NULL,
        best_models_ggplot_mean = NULL
      )
    },
    results_handling = function(..., run_metadata = NULL) {
      list(run_metadata = run_metadata)
    },
    gp_execute_best_model_tasks = function(...) list(),
    gp_finalize_model_execute_results = function(..., run_metadata = NULL) {
      list(run_metadata = run_metadata)
    },
    .package = "PredictProR"
  )

  expect_silent(
    PredictProR::model_execute(
      pheno_data = pheno,
      omic1_data = omics,
      response = "Yield",
      gen_name = "GID",
      GS_model_cv = "mlp",
      cross_validation = TRUE,
      cv_evaluation_only = TRUE,
      system_database = TRUE,
      eval_metrics = "mean_squared_error",
      metric_for_ranking = "mean_squared_error",
      verbose = FALSE
    )
  )
})
