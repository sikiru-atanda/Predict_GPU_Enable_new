test_that("parallel memory plan inventories the complete serialized payload", {
  pheno <- data.frame(GID = paste0("g", 1:100), y = rnorm(100))
  K <- diag(100)
  dimnames(K) <- list(pheno$GID, pheno$GID)
  params <- list(pheno_data = pheno, gmatrix_model_ready = K)

  plan <- parallel_memory_plan(
    objects = list(params = params, pheno_data = pheno),
    n_tasks = 8L,
    workers = 4L,
    memory_budget_gb = 1,
    payload_multiplier = 1.25,
    worker_overhead_gb = 0.1
  )

  expect_identical(plan$inventory$name[[1L]], "params")
  expect_gt(plan$payload_gb, as.numeric(object.size(K)) / 1024^3)
  expect_equal(plan$per_worker_gb, plan$payload_gb * 1.25 + 0.1, tolerance = 1e-12)
  expect_lte(plan$workers_recommended, plan$workers_requested)
  expect_lte(plan$projected_worker_memory_gb, plan$memory_budget_gb + 1e-12)
})

test_that("parallel memory cap includes worker process overhead", {
  withr::local_envvar(c(
    GP_PAR_PAYLOAD_MULTIPLIER = "1",
    GP_PAR_WORKER_OVERHEAD_GB = "0.25"
  ))

  expect_identical(
    PredictProR:::gp_parallel_memory_worker_cap(0.25, 1),
    2L
  )
})

test_that("cross-validation policy exposes projected worker memory", {
  local_mocked_bindings(
    detect_num_gpus = function() 0L,
    detect_hidden_threads = function() 1L,
    gp_parallel_memory_budget_gb = function() 2,
    .package = "PredictProR"
  )
  withr::local_envvar(c(
    GP_PAR_PAYLOAD_MULTIPLIER = "1.25",
    GP_PAR_WORKER_OVERHEAD_GB = "0.1"
  ))
  globals <- list(params = list(kernel = diag(100)))
  result <- PredictProR:::sp_decide_policy(
    models = rep("RandomForest", 8),
    model_params_list = rep(list(list()), 8),
    globals = globals,
    n_tasks = 8L,
    user_mode = "base_parallel",
    num_cores = 4L,
    verbose = FALSE,
    sys_name = "Linux",
    prefer_fork = TRUE
  )

  expect_true(is.finite(result$per_worker_memory_gb))
  expect_equal(
    result$projected_worker_memory_gb,
    result$workers * result$per_worker_memory_gb,
    tolerance = 1e-12
  )
  expect_lte(result$projected_worker_memory_gb, result$memory_budget_gb + 1e-12)
})
