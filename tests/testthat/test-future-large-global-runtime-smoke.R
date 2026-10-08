test_that("future executes a shared payload within the PredictProR memory budget", {
  skip_if_not(
    identical(Sys.getenv("PREDICTPROR_RUN_PARALLEL_MEMORY_SMOKE", unset = "0"), "1"),
    "Set PREDICTPROR_RUN_PARALLEL_MEMORY_SMOKE=1 for the 128 MiB future smoke."
  )
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("withr")
  withr::local_envvar(c(GP_PAR_MEMORY_BUDGET_GB = "2"))
  old_option <- getOption("future.globals.maxSize", NULL)
  on.exit(options(future.globals.maxSize = old_option), add = TRUE)
  options(future.globals.maxSize = NULL)

  count <- as.integer(128 * 1024^2 / 8)
  payload <- numeric(count)
  for (start in seq.int(1L, count, by = 1000000L)) {
    end <- min(count, start + 1000000L - 1L)
    index <- seq.int(start, end)
    payload[index] <- sin(index * 0.000013) + cos(index * 0.000007)
  }
  task_fun <- local({
    shared <- payload
    function(i) {
      index <- seq.int(as.integer(i), length(shared), by = 257L)
      c(task = as.integer(i), value = sum(shared[index]),
        random = stats::runif(1L))
    }
  })

  sequential <- PredictProR:::sp_apply(
    X = 1:4, FUN = task_fun,
    decision = list(backend = "sequential", workers = 1L),
    packages = character(), seed = 20260902L
  )
  decision <- PredictProR:::sp_decide_policy(
    models = rep("Lasso", 4L),
    model_params_list = rep(list(list()), 4L),
    globals = list(payload = payload), n_tasks = 4L,
    user_mode = "future", num_cores = 2L, verbose = FALSE
  )
  expect_identical(decision$backend, "future")
  expect_identical(decision$workers, 2L)
  expect_equal(decision$memory_budget_gb, 2, tolerance = 0)

  observed <- PredictProR:::sp_apply(
    X = 1:4, FUN = task_fun, decision = decision,
    packages = character(), seed = 20260902L
  )

  expect_equal(observed, sequential, tolerance = 0)
  expect_null(getOption("future.globals.maxSize", NULL))
})

