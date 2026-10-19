# base_parallel's PSOCK path sent the task data to every worker twice: the
# package-loading helper was a closure created inside sp_apply(), so it
# carried sp_apply()'s frame (and the task function's data) with it.

test_that("PSOCK workers receive the task function's data once", {
  big <- stats::runif(2e6)  # 16 MB the task closes over
  task <- function(i) i + length(big)
  sent <- numeric()
  setup <- PredictProR:::gp_parallel_psock_worker_setup
  local_mocked_bindings(
    gp_make_psock_cluster = function(workers) structure(list(), class = "fake_cluster"),
    .package = "PredictProR"
  )
  local_mocked_bindings(
    clusterCall = function(cl, fun, ...) {
      sent <<- c(sent, length(serialize(list(fun, list(...)), NULL)))
      # run only the task set-up here; the worker runtime set-up would change this session
      if (identical(fun, setup)) list(fun(...)) else list(NULL)
    },
    parLapplyLB = function(cl, X, fun, ...) lapply(X, fun),
    stopCluster = function(cl) invisible(NULL),
    .package = "parallel"
  )
  withr::defer(if (exists(".gp_sp_task_fun", envir = globalenv())) rm(".gp_sp_task_fun", envir = globalenv()))

  out <- PredictProR:::sp_apply(
    as.list(1:3), task,
    list(backend = "base_parallel", workers = 2L, prefer_fork = FALSE),
    packages = "stats"
  )
  expect_equal(unlist(out), 1:3 + 2e6)
  expect_gt(sum(sent), 16e6)        # the data did go to the worker
  expect_lt(sum(sent), 1.3 * 16e6)  # once
})

# sp_apply() wrapped each seeded task in a closure made in its own frame. An
# argument not yet evaluated there (fail_fast, on the PSOCK path) is a promise
# whose environment is the caller's frame, and that frame went to every worker
# with the task: a 3-trait CV's 61 MB task travelled as 1,026 MB per worker.
test_that("a seeded task goes to workers without its caller's frame", {
  task <- local(function(i) i + 1, envir = new.env(parent = globalenv()))
  sent <- numeric()
  setup <- PredictProR:::gp_parallel_psock_worker_setup
  local_mocked_bindings(
    gp_make_psock_cluster = function(workers) structure(list(), class = "fake_cluster"),
    .package = "PredictProR"
  )
  local_mocked_bindings(
    clusterCall = function(cl, fun, ...) {
      sent <<- c(sent, length(serialize(list(fun, list(...)), NULL)))
      if (identical(fun, setup)) list(fun(...)) else list(NULL)
    },
    parLapplyLB = function(cl, X, fun, ...) lapply(X, fun),
    stopCluster = function(cl) invisible(NULL),
    .package = "parallel"
  )
  withr::defer(if (exists(".gp_sp_task_fun", envir = globalenv())) rm(".gp_sp_task_fun", envir = globalenv()))
  caller <- function() {
    unused <- stats::runif(2e6)  # 16 MB the task never reads
    flag <- FALSE
    PredictProR:::sp_apply(as.list(1:3), task, list(backend = "base_parallel", workers = 2L, prefer_fork = FALSE),
                           packages = "stats", seed = TRUE, fail_fast = flag)
  }
  out <- caller()
  expect_equal(unlist(out), 2:4)
  # 17.2 MB with the caller's frame; ~1 MB without it under load_all(), which
  # keeps each package function's source file
  expect_lt(sum(sent), 4e6)
})

# The per-task closure built in gp_execute_partitioned_tasks() carried that
# function's frame: the policy's `globals` list and its mori-prepared copy,
# the run's data twice more (heavy CV: 314 MB shipped per worker, ~100 MB of
# it the data, measured 2026-10-10).
test_that("the task function the partitioned scheduler hands to workers carries the data once", {
  big <- stats::runif(2e6)  # 16 MB the task closes over
  # the task's frame holds the data; its parent is the global env, so
  # serialising the task does not walk into this test's frame as well
  run_task <- local(function(i) i + length(data), envir = list2env(list(data = big), parent = globalenv()))
  shipped <- NA_real_
  local_mocked_bindings(
    sp_apply = function(X, FUN, decision, ...) {
      shipped <<- length(serialize(FUN, NULL))
      lapply(X, FUN)
    },
    detect_num_gpus = function() 0L,
    .package = "PredictProR"
  )
  out <- PredictProR:::gp_execute_partitioned_tasks(
    model_ids = c("BayesB", "BayesB"),
    run_task = run_task,
    decision_args = list(
      model_params_list = list(list(), list()), globals = list(data = big),
      user_mode = "base_parallel", num_cores = 2L, verbose = FALSE, worker_memory_gb = 0.5
    ),
    seed = FALSE
  )
  expect_equal(unlist(out), 1:2 + 2e6)
  expect_gt(shipped, 16e6)
  expect_lt(shipped, 1.3 * 16e6)
})
