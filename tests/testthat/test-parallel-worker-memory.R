test_that("worker memory resolves user value, then env var, then measurement, then fallback", {
  withr::local_envvar(GP_PAR_WORKER_OVERHEAD_GB = NA)
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb("BayesB", worker_memory_gb = 2.5, measure = FALSE),
    list(gb = 2.5, source = "user_worker_memory_gb")
  )
  withr::local_envvar(GP_PAR_WORKER_OVERHEAD_GB = "1.2")
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb("BayesB", measure = FALSE),
    list(gb = 1.2, source = "env_GP_PAR_WORKER_OVERHEAD_GB")
  )
  # the user argument still wins over the env var
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb("BayesB", worker_memory_gb = 3, measure = FALSE)$gb,
    3
  )
  withr::local_envvar(GP_PAR_WORKER_OVERHEAD_GB = NA)
  local_mocked_bindings(
    gp_parallel_measure_worker_gb = function(python_models = TRUE) if (python_models) 1.4 else 0.6,
    .package = "PredictProR"
  )
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb(c("BayesB", "RandomForest"), measure = TRUE),
    list(gb = 1.4, source = "measured_r_worker_plus_python")
  )
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb(c("BayesB", "GBLUP_BRR"), measure = TRUE),
    list(gb = 0.6, source = "measured_r_worker")
  )
  expect_identical(
    PredictProR:::gp_parallel_resolve_worker_memory_gb("BayesB", measure = FALSE),
    list(gb = 0.75, source = "default_0.75")
  )
})

test_that("worker memory and memory budget cap the number of workers", {
  withr::local_envvar(GP_PAR_PAYLOAD_MULTIPLIER = "1")
  expect_identical(PredictProR:::gp_parallel_memory_worker_cap(0.5, 8, worker_overhead_gb = 1.5), 4L)
  expect_identical(PredictProR:::gp_parallel_memory_worker_cap(0.5, 8, worker_overhead_gb = 3.5), 2L)
  # tiny shared data: the worker allowance alone still caps the workers
  expect_identical(PredictProR:::gp_parallel_memory_worker_cap(0, 4, worker_overhead_gb = 1), 4L)
})

test_that("the CV policy honours user worker_memory_gb and memory_budget_gb", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    .package = "PredictProR"
  )
  globals <- list(params = list(kernel = diag(10)))
  policy <- PredictProR:::sp_decide_policy(
    models = rep("BayesB", 8),
    model_params_list = rep(list(list()), 8),
    globals = globals,
    n_tasks = 8L,
    user_mode = "base_parallel",
    num_cores = 8L,
    verbose = FALSE,
    sys_name = "Windows",
    worker_memory_gb = 2,
    memory_budget_gb = 6
  )
  expect_identical(policy$worker_memory_source, "user_worker_memory_gb")
  expect_equal(policy$worker_memory_gb, 2)
  expect_equal(policy$memory_budget_gb, 6)
  expect_identical(as.integer(policy$workers), 2L)
})

test_that("model_execute rejects invalid memory settings", {
  expect_error(PredictProR::model_execute(worker_memory_gb = -1), "worker_memory_gb must be NULL or one positive number")
  expect_error(PredictProR::model_execute(memory_budget_gb = "8"), "memory_budget_gb must be NULL or one positive number")
})

test_that("ML bootstrap process count respects override, worker flag, cores and replicates", {
  withr::local_envvar(PREDICTPRO_BOOTSTRAP_JOBS = "3", PREDICTPRO_IN_PARALLEL_WORKER = NA)
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(100), 3L)
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(2), 2L)
  withr::local_envvar(PREDICTPRO_BOOTSTRAP_JOBS = NA, PREDICTPRO_IN_PARALLEL_WORKER = "1",
                      PREDICTPRO_PARALLEL_WORKERS = NA)
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(100), 1L)
  local_mocked_bindings(
    gp_detect_cores = function(...) 7L,
    gp_parallel_memory_budget_gb = function() 2,
    .package = "PredictProR"
  )
  # a worker that shares the machine with 2 others gets floor(7 / 3) = 2 jobs
  withr::local_envvar(PREDICTPRO_PARALLEL_WORKERS = "3")
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(100), 1L)  # memory: (2 / 3) / 0.5 = 1
  local_mocked_bindings(gp_parallel_memory_budget_gb = function() 12, .package = "PredictProR")
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(100), 2L)
  withr::local_envvar(PREDICTPRO_IN_PARALLEL_WORKER = NA)
  local_mocked_bindings(gp_parallel_memory_budget_gb = function() 2, .package = "PredictProR")
  # memory: 2 GB / (0.5 + 0) = 4 jobs; cores 7; replicates 100
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(100), 4L)
  expect_identical(PredictProR:::gp_ml_bootstrap_jobs(3), 3L)
})