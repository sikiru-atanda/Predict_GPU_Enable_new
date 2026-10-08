# Policy tests check backend logic, not this machine's free RAM: with little
# free memory the gate correctly picks sequential (memory_limited_workers).
withr::local_envvar(c(GP_PAR_MEMORY_BUDGET_GB = "16"))   # reset when this file finishes
# Many tests below use tree models as generic parallel-eligible tasks. By
# default those are auto-threaded (and therefore run sequentially with all
# cores); the documented opt-out keeps them single-threaded here. The default
# behaviour is tested explicitly at the end of this file.
withr::local_envvar(c(PREDICTPRO_ML_AUTO_THREADS = "false"))   # reset when this file finishes

test_that("parallel runtime environment propagates project library paths", {
  skip_if_not_installed("withr")

  temp_lib <- tempfile("predictpror-lib path ")
  dir.create(temp_lib, recursive = TRUE)
  old_libs <- .libPaths()
  old_env <- Sys.getenv(c("R_LIBS", "PREDICTPRO_R_LIBPATHS"), unset = NA_character_)
  withr::defer(.libPaths(old_libs))
  withr::defer({
    if (is.na(old_env[["R_LIBS"]])) {
      Sys.unsetenv("R_LIBS")
    } else {
      Sys.setenv(R_LIBS = old_env[["R_LIBS"]])
    }
    if (is.na(old_env[["PREDICTPRO_R_LIBPATHS"]])) {
      Sys.unsetenv("PREDICTPRO_R_LIBPATHS")
    } else {
      Sys.setenv(PREDICTPRO_R_LIBPATHS = old_env[["PREDICTPRO_R_LIBPATHS"]])
    }
  })

  .libPaths(c(temp_lib, old_libs))
  env <- PredictProR:::gp_parallel_runtime_env()
  env <- env[setdiff(names(env), c("PREDICTPRO_GP_PYTHON", "PREDICTPRO_PYTHON", "RETICULATE_PYTHON"))]
  normalized_lib <- normalizePath(temp_lib, winslash = "/", mustWork = TRUE)

  expect_true(grepl(normalized_lib, env[["R_LIBS"]], fixed = TRUE))
  expect_true(grepl(normalized_lib, env[["PREDICTPRO_R_LIBPATHS"]], fixed = TRUE))

  .libPaths(old_libs)
  PredictProR:::gp_parallel_init_runtime(env)
  expect_true(normalized_lib %in% normalizePath(.libPaths(), winslash = "/", mustWork = FALSE))
})

test_that("auto policy chooses sequential when too few tasks remain parallel-eligible", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("RandomForest"),
    model_params_list = list(list()),
    globals = list(
      pheno_data = data.frame(y = 1:10),
      omics_data = matrix(1, nrow = 2, ncol = 2),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 1L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )

  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "no_parallel_eligible_tasks")
  expect_true(is.data.frame(res$backend_candidates))
  expect_true("sequential" %in% res$backend_candidates$backend)
})

test_that("auto policy avoids PSOCK startup cost for tiny CPU queues on Windows", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("RandomForest", "Ridge_Regression", "Lasso"),
    model_params_list = list(list(), list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(rnorm(400), nrow = 40),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 3L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 3L,
    sys_name = "Windows",
    prefer_fork = FALSE
  )

  expect_identical(res$backend, "sequential")
  expect_true(grepl("small_windows_cpu_queue", res$decision_reason, fixed = TRUE))
  expect_true(all(c("future", "mirai", "sequential") %in% res$backend_candidates$backend))
})

test_that("auto policy avoids worker startup cost for tiny CPU queues on Unix", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("RandomForest", "Ridge_Regression", "Lasso"),
    model_params_list = list(list(), list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(rnorm(400), nrow = 40),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 3L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 3L,
    sys_name = "Linux",
    prefer_fork = TRUE
  )

  if (!identical(res$backend, "sequential")) {
    stop(
      paste("Expected sequential backend, got", res$backend, "with reason", res$decision_reason),
      call. = FALSE
    )
  }
  expect_true(grepl("small_cpu_queue", res$decision_reason, fixed = TRUE))
})

test_that("auto policy prefers mirai for larger CPU queues on Windows", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  models <- rep(c("RandomForest", "Ridge_Regression", "Lasso", "PartialLeastSquare"), 2L)
  res <- PredictProR:::sp_decide_policy(
    models = models,
    model_params_list = rep(list(list()), length(models)),
    globals = list(
      pheno_data = data.frame(y = seq_len(200L)),
      omics_data = matrix(rnorm(4000), nrow = 200L),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = length(models),
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L,
    sys_name = "Windows",
    prefer_fork = FALSE
  )

  if (!identical(res$backend, "mirai")) {
    stop(
      paste("Expected mirai backend, got", res$backend, "with reason", res$decision_reason),
      call. = FALSE
    )
  }
  expect_true(grepl("low_overhead_worker_pool", res$decision_reason, fixed = TRUE))
  expect_gt(
    res$backend_candidates$score[res$backend_candidates$backend == "mirai"],
    res$backend_candidates$score[res$backend_candidates$backend == "future"]
  )
})

test_that("mori signal biases large shared-global CPU queues toward mirai", {
  old_use <- Sys.getenv("GP_USE_MORI", unset = NA_character_)
  old_min <- Sys.getenv("GP_MORI_MIN_MB", unset = NA_character_)
  on.exit({
    if (is.na(old_use)) Sys.unsetenv("GP_USE_MORI") else Sys.setenv(GP_USE_MORI = old_use)
    if (is.na(old_min)) Sys.unsetenv("GP_MORI_MIN_MB") else Sys.setenv(GP_MORI_MIN_MB = old_min)
  }, add = TRUE)
  Sys.setenv(GP_USE_MORI = "true", GP_MORI_MIN_MB = "0")

  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    gp_mori_available = function() TRUE,
    .package = "PredictProR"
  )

  models <- rep(c("RandomForest", "Ridge_Regression", "Lasso", "PartialLeastSquare"), 3L)
  res <- PredictProR:::sp_decide_policy(
    models = models,
    model_params_list = rep(list(list()), length(models)),
    globals = list(
      pheno_data = data.frame(y = seq_len(400L)),
      omics_data = matrix(rnorm(40000), nrow = 400L),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = length(models),
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 3L,
    sys_name = "Windows",
    prefer_fork = FALSE
  )

  if (!identical(res$backend, "mirai")) {
    stop(
      paste("Expected mori-biased mirai backend, got", res$backend, "with reason", res$decision_reason),
      call. = FALSE
    )
  }
  expect_true(grepl("mori_shared_globals", res$decision_reason, fixed = TRUE))
})

test_that("auto policy prefers mirai for multi-GPU deep workloads", {
  local_mocked_bindings(
    detect_num_gpus = function() 2L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("mlp", "ft_transformer"),
    model_params_list = list(list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:20),
      omics_data = matrix(rnorm(200), nrow = 20),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 4L,
    sys_name = "Windows"
  )

  expect_identical(res$backend, "mirai")
  expect_true(grepl("multi_gpu_model_bonus", res$decision_reason, fixed = TRUE))
  expect_gt(
    res$backend_candidates$score[res$backend_candidates$backend == "mirai"],
    res$backend_candidates$score[res$backend_candidates$backend == "future"]
  )
  expect_true(isTRUE(res$gpu_parallel_enabled))
})

test_that("forced backend overrides scoring while preserving explanation", {
  res <- PredictProR:::sp_decide_policy(
    models = c("RandomForest", "Ridge_Regression"),
    model_params_list = list(list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:20),
      omics_data = matrix(rnorm(200), nrow = 20),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "future",
    verbose = FALSE,
    num_cores = 2L,
    sys_name = "Windows",
    prefer_fork = FALSE
  )

  expect_identical(res$backend, "future")
  expect_identical(res$decision_reason, "user_forced_future")
  expect_true(is.infinite(res$backend_score))
  expect_identical(nrow(res$backend_candidates), 1L)
})

test_that("capability table captures internal-threaded models and RKHS BLAS sensitivity", {
  caps <- PredictProR:::gp_parallel_model_capability_table(
    models = c("Xgboost", "RKHS", "RandomForest"),
    model_params_list = list(
      list(xgb_nthread = 4L),
      list(),
      list()
    ),
    num_gpus = 0L,
    hidden_thr = 4L,
    has_asreml_prep = FALSE,
    rkhs_kernel_GB = 2,
    rkhs_mem_threshold_GB = 1,
    always_seq_rkhs = FALSE
  )

  expect_true(caps$forced_sequential[[1]])
  expect_true(caps$forced_sequential[[2]])
  expect_false(caps$forced_sequential[[3]])
  expect_true(grepl("xgboost_internal_threads", caps$reasons[[1]], fixed = TRUE))
  expect_true(grepl("rkhs_hidden_threads", caps$reasons[[2]], fixed = TRUE))
})

test_that("policy reports RKHS-specific sequential reason when RKHS is the only blocker", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 4L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("RKHS", "RandomForest"),
    model_params_list = list(list(), list()),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )

  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "rkhs_hidden_threads")
})

test_that("policy reports Xgboost-specific sequential reason when Xgboost is the only blocker", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("Xgboost", "RandomForest"),
    model_params_list = list(
      list(xgb_nthread = 4L),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )

  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "xgboost_internal_threads")
})

test_that("policy reports CatBoost and LightGBM specific sequential reasons", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res_cat <- PredictProR:::sp_decide_policy(
    models = c("CatBoost", "RandomForest"),
    model_params_list = list(
      list(catboost_thread_count = 4L),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )
  expect_identical(res_cat$decision_reason, "catboost_internal_threads")

  res_lgbm <- PredictProR:::sp_decide_policy(
    models = c("LightGBM", "RandomForest"),
    model_params_list = list(
      list(lightgbm_nthread = 4L),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(1, nrow = 4, ncol = 4),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L
  )
  expect_identical(res_lgbm$decision_reason, "lightgbm_internal_threads")
})

test_that("policy reports multiple constraints when more than one sequential cause remains", {
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

  expect_identical(res$backend, "sequential")
  expect_identical(res$decision_reason, "multiple_sequential_constraints")
})

test_that("single-GPU xgboost is treated as GPU-exclusive and kept sequential", {
  local_mocked_bindings(
    detect_num_gpus = function() 1L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )

  res <- PredictProR:::sp_decide_policy(
    models = c("Xgboost", "RandomForest"),
    model_params_list = list(
      list(xgb_booster = "dart"),
      list()
    ),
    globals = list(
      pheno_data = data.frame(y = 1:40),
      omics_data = matrix(rnorm(400), nrow = 40),
      bayes = NULL,
      asreml = NULL
    ),
    n_tasks = 2L,
    user_mode = "auto",
    verbose = FALSE,
    num_cores = 2L,
    sys_name = "Windows"
  )

  expect_true(res$model_capabilities$forced_sequential[[1]])
  expect_true(grepl("single_gpu_exclusive", res$model_capabilities$reasons[[1]], fixed = TRUE))
  expect_false(isTRUE(res$gpu_parallel_enabled))
})

test_that("policy parameter extraction keeps only orchestration-relevant fields", {
  out <- PredictProR:::gp_parallel_extract_model_policy_params(
    "Xgboost",
    list(
      xgb_nthread = 4L,
      xgb_booster = "dart",
      device = "cuda",
      learning_rate = 0.1,
      random_seed = 123L
    )
  )

  expect_identical(sort(names(out)), sort(c("xgb_nthread", "xgb_booster", "device")))
  expect_false("learning_rate" %in% names(out))
  expect_false("random_seed" %in% names(out))
})

test_that("mirai production policy exposes bounded queue settings", {
  old <- Sys.getenv("GP_MIRAI_QUEUE_MEMORY_MB", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("GP_MIRAI_QUEUE_MEMORY_MB") else Sys.setenv(GP_MIRAI_QUEUE_MEMORY_MB = old)
  }, add = TRUE)

  Sys.setenv(GP_MIRAI_QUEUE_MEMORY_MB = "256")
  expect_identical(PredictProR:::gp_mirai_queue_memory_mb(workers = 4L, n_tasks = 12L), 256)
  expect_identical(PredictProR:::gp_mirai_queue_memory_bytes(workers = 4L, n_tasks = 12L), 256 * 1024^2)
})

test_that("mirai backend uses bounded non-blocking submission before collecting results", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("purrr")

  submitted <- 0L
  collected_after_submissions <- integer()
  daemon_memory <- NULL
  fake_submit <- function(expr, .args = list(), .compute = NULL) {
    submitted <<- submitted + 1L
    out <- new.env(parent = emptyenv())
    makeActiveBinding("value", function() {
      collected_after_submissions <<- c(collected_after_submissions, submitted)
      paste0("result_", length(collected_after_submissions))
    }, out)
    out
  }

  local_mocked_bindings(
    start_metrics_endpoint = function(...) NULL,
    stop_metrics_endpoint = function(...) invisible(NULL),
    stop_daemons = function(...) invisible(NULL),
    get_gpu_usage = function(...) NA_real_,
    gp_mirai_start_daemons = function(n, .compute = NULL, .expr = NULL, sync = FALSE, memory = NULL) {
      daemon_memory <<- memory
      invisible(NULL)
    },
    gp_mirai_submit = fake_submit,
    gp_mirai_collect = function(x) x$value,
    gp_mirai_is_error_value = function(x) FALSE,
    gp_mirai_queue_depth = function() 0L,
    gp_mirai_queue_memory_bytes = function(...) 64 * 1024^2,
    .package = "PredictProR"
  )

  out <- PredictProR:::sp_apply(
    X = 1:3,
    FUN = function(i) i,
    decision = list(
      backend = "mirai",
      workers = 3L,
      num_gpus = 0L,
      gpu_parallel_enabled = FALSE,
      gpu_ids = character()
    ),
    packages = character(),
    seed = FALSE
  )

  expect_length(out, 3L)
  expect_identical(collected_after_submissions, rep(3L, 3L))
  expect_identical(daemon_memory, 64 * 1024^2)
})

test_that("base_parallel PSOCK workers initialize runtime environment", {
  skip_if_not(identical(.Platform$OS.type, "windows"), "PSOCK branch is Windows-specific here")

  out <- PredictProR:::sp_apply(
    X = 1:2,
    FUN = function(i) i + 1L,
    decision = list(
      backend = "base_parallel",
      workers = 1L
    ),
    packages = character()
  )

  expect_identical(unlist(out), 2:3)
})

test_that("parallel policy sizes every supplied global and detects current ASReml names", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    gp_parallel_memory_budget_gb = function() Inf,
    .package = "PredictProR"
  )

  custom_matrix <- matrix(1, nrow = 50L, ncol = 20L)
  asreml_models_prep_cv <- list(fitted = TRUE)
  globals <- list(
    pheno_clean = list(pheno_clean_data = data.frame(y = seq_len(50L))),
    custom_matrix = custom_matrix,
    asreml_models_prep_cv = asreml_models_prep_cv
  )
  expected_gb <- PredictProR:::gp_parallel_globals_size_bytes(globals) / 1024^3

  out <- PredictProR:::sp_decide_policy(
    models = c("GBLUP", "RandomForest"),
    model_params_list = list(list(), list()),
    globals = globals,
    n_tasks = 2L,
    user_mode = "auto",
    num_cores = 2L,
    sys_name = "Windows",
    prefer_fork = FALSE,
    verbose = FALSE
  )

  expect_equal(out$globals_gb, expected_gb, tolerance = 1e-12)
  expect_true(out$model_capabilities$forced_sequential[out$model_capabilities$model == "GBLUP"])
  expect_match(
    out$model_capabilities$reasons[out$model_capabilities$model == "GBLUP"],
    "asreml_prepared",
    fixed = TRUE
  )
})

test_that("worker fanout respects the available-memory budget on every OS", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    gp_parallel_memory_budget_gb = function() 0.005,
    .package = "PredictProR"
  )

  globals <- list(arbitrary_large_input = matrix(0, nrow = 1000L, ncol = 1000L))
  for (sys_name in c("Windows", "Linux", "Darwin")) {
    out <- PredictProR:::sp_decide_policy(
      models = c("RandomForest", "Lasso"),
      model_params_list = list(list(), list()),
      globals = globals,
      n_tasks = 2L,
      user_mode = "auto",
      num_cores = 2L,
      globals_max_GB = 1,
      sys_name = sys_name,
      prefer_fork = FALSE,
      verbose = FALSE
    )

    expect_identical(out$backend, "sequential", info = sys_name)
    expect_identical(out$decision_reason, "memory_limited_workers", info = sys_name)
    expect_identical(out$workers, 1L, info = sys_name)
    expect_true(out$memory_worker_cap < out$workers_requested, info = sys_name)
    expect_false(out$prefer_fork, info = sys_name)
  }
})

test_that("sp_apply gives every backend deterministic task-specific RNG streams", {
  seed <- 20260811L
  sequential <- PredictProR:::sp_apply(
    X = 1:3,
    FUN = function(i) stats::runif(4L),
    decision = list(backend = "sequential", workers = 1L),
    packages = character(),
    seed = seed
  )
  base_parallel <- PredictProR:::sp_apply(
    X = 1:3,
    FUN = function(i) stats::runif(4L),
    decision = list(backend = "base_parallel", workers = 2L),
    packages = character(),
    seed = seed
  )

  expect_equal(base_parallel, sequential, tolerance = 0)
  expect_false(identical(sequential[[1L]], sequential[[2L]]))
})

test_that("future global limits follow the per-worker memory budget", {
  derived <- PredictProR:::gp_future_globals_limit_bytes(
    decision = list(memory_budget_gb = 2, workers = 2L),
    current_option = NULL
  )
  expect_equal(derived$bytes, 1024^3, tolerance = 0)
  expect_identical(derived$source, "predictpror_per_worker_memory_budget")
  expect_true(derived$temporary)

  explicit <- PredictProR:::gp_future_globals_limit_bytes(
    decision = list(memory_budget_gb = 2, workers = 2L),
    current_option = 64 * 1024^2
  )
  expect_equal(explicit$bytes, 64 * 1024^2, tolerance = 0)
  expect_identical(explicit$source, "user_option")
  expect_false(explicit$temporary)
})

test_that("PSOCK ports are process-distinct unless the user sets one", {
  skip_if_not_installed("withr")
  first <- PredictProR:::gp_psock_port_candidates(pid = 1000L)
  second <- PredictProR:::gp_psock_port_candidates(pid = 1001L)

  expect_length(first, 64L)
  expect_length(unique(first), 64L)
  expect_true(all(first >= 20000L & first <= 59999L))
  expect_false(identical(first[[1L]], second[[1L]]))

  withr::local_envvar(c(R_PARALLEL_PORT = "24567"))
  expect_identical(PredictProR:::gp_psock_cluster_port(), 24567L)
})

test_that("base and foreach PSOCK paths use the collision-resistant helper", {
  source_text <- paste(deparse(body(PredictProR:::sp_apply)), collapse = "\n")
  expect_match(source_text, "gp_make_psock_cluster", fixed = TRUE)
  expect_false(grepl("parallel::makePSOCKcluster", source_text, fixed = TRUE))
})

test_that("sp_apply restores an absent future globals option", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  old_option <- getOption("future.globals.maxSize", NULL)
  on.exit(options(future.globals.maxSize = old_option), add = TRUE)
  options(future.globals.maxSize = NULL)

  observed <- PredictProR:::sp_apply(
    X = 1:2,
    FUN = function(i) i + 1L,
    decision = list(
      backend = "future", workers = 1L, plan = future::multisession,
      chunk_size = 1L, memory_budget_gb = 0.75
    ),
    packages = character(),
    seed = 20260902L
  )

  expect_identical(unlist(observed, use.names = FALSE), 2:3)
  expect_null(getOption("future.globals.maxSize", NULL))
})

test_that("gp_set_seed repairs an invalid global RNG binding", {
  set.seed(991L)
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)

  assign(".Random.seed", NULL, envir = .GlobalEnv)
  expect_silent(PredictProR:::gp_set_seed(123L))
  expect_true(is.integer(.Random.seed))
  observed <- stats::runif(3L)

  PredictProR:::gp_set_seed(123L)
  expect_identical(stats::runif(3L), observed)
})

test_that("PSOCK workers apply the configured BLAS thread cap without mutating the parent", {
  skip_if_not_installed("withr")
  skip_if_not(identical(.Platform$OS.type, "windows"), "PSOCK branch is Windows-specific here")
  withr::local_envvar(c(GP_BLAS_THREADS = "1", OMP_NUM_THREADS = "7"))

  worker_value <- PredictProR:::sp_apply(
    X = 1L,
    FUN = function(i) Sys.getenv("OMP_NUM_THREADS"),
    decision = list(backend = "base_parallel", workers = 1L),
    packages = character(),
    seed = FALSE
  )

  expect_identical(worker_value[[1L]], "1")
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "7")
})

test_that("auto-threaded tree models run sequentially with all cores; BGLR models stay parallel", {
  withr::local_envvar(c(PREDICTPRO_ML_AUTO_THREADS = "true", GP_PAR_MEMORY_BUDGET_GB = "64"))
  globals <- list(x = matrix(0, 10, 10))
  for (m in c("RandomForest", "Xgboost", "CatBoost", "LightGBM")) {
    caps <- PredictProR:::sp_decide_policy(
      models = rep(m, 8L), model_params_list = replicate(8L, list(), simplify = FALSE),
      globals = globals, n_tasks = 8L, user_mode = "auto", num_cores = 8L, verbose = FALSE
    )
    expect_identical(caps$backend, "sequential", info = m)
    expect_true(all(caps$capability_table$thread_sensitive), info = m)
  }
  bglr <- PredictProR:::sp_decide_policy(
    models = rep("BayesB", 8L), model_params_list = replicate(8L, list(), simplify = FALSE),
    globals = globals, n_tasks = 8L, user_mode = "auto", num_cores = 8L, verbose = FALSE
  )
  expect_false(any(bglr$capability_table$thread_sensitive))
})

test_that("the policy sizes what a task closure actually ships", {
  big <- matrix(0, 2000, 2000)        # ~30 MB captured by the closure below
  task <- local({ data <- big; function(i) sum(data[i, ]) })
  expect_gt(PredictProR:::gp_parallel_closure_payload_gb(task), 0.025)
  res <- PredictProR:::sp_decide_policy(
    models = rep("BayesB", 4L), model_params_list = replicate(4L, list(), simplify = FALSE),
    globals = list(), n_tasks = 4L, user_mode = "auto", num_cores = 4L, verbose = FALSE,
    payload_gb = 3
  )
  expect_gte(res$per_worker_memory_gb, 3 * 4)   # payload x GP_PAR_PAYLOAD_MULTIPLIER (4)
})
