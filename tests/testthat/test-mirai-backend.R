# sp_apply(backend = "mirai") used to start daemons with sync = TRUE, which in
# mirai >= 2.x evaluates every task in the calling process, and passed the
# daemon set-up as `.expr`, which daemons() ignores: tasks ran sequentially in
# the session, and no daemon got the runtime env or its own GPU.

mirai_decision <- function(workers = 2L, gpu_ids = character()) {
  list(
    backend = "mirai",
    workers = workers,
    num_gpus = length(gpu_ids),
    gpu_parallel_enabled = length(gpu_ids) > 1L,
    gpu_ids = gpu_ids
  )
}

skip_if_no_mirai <- function() {
  skip_on_cran()
  skip_if_not_installed("mirai", "2.7.0")
  skip_if_not_installed("purrr")
}

test_that("mirai tasks run in set-up daemon processes, not in the session", {
  skip_if_no_mirai()
  withr::local_envvar(PREDICTPRO_IN_PARALLEL_WORKER = NA)
  out <- PredictProR:::sp_apply(
    as.list(1:4),
    function(i) {
      Sys.sleep(0.3)
      list(i = i, pid = Sys.getpid(), worker = Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER"))
    },
    mirai_decision(workers = 2L),
    packages = character(),
    seed = FALSE
  )
  pids <- vapply(out, function(r) as.numeric(r$pid), numeric(1))
  expect_identical(vapply(out, function(r) r$i, integer(1)), 1:4)
  expect_false(any(pids == Sys.getpid()))
  expect_length(unique(pids), 2L)
  # set by gp_parallel_init_worker_runtime() on every daemon
  expect_true(all(vapply(out, function(r) r$worker, character(1)) == "1"))
  expect_false(identical(Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER"), "1"))
})

test_that("mirai results equal the sequential results for the same seed", {
  skip_if_no_mirai()
  task <- function(i) c(i, stats::runif(2))
  seq_out <- PredictProR:::sp_apply(as.list(1:5), task, list(backend = "sequential", workers = 1L),
                                    packages = character(), seed = 42L)
  mirai_out <- PredictProR:::sp_apply(as.list(1:5), task, mirai_decision(workers = 2L),
                                      packages = character(), seed = 42L)
  expect_identical(mirai_out, seq_out)
})

test_that("the task function is sent to each daemon once, not with every task", {
  skip_if_no_mirai()
  big <- stats::runif(2e6)  # 16 MB the task closes over
  sent <- list()
  real_submit <- PredictProR:::gp_mirai_submit
  local_mocked_bindings(
    gp_mirai_submit = function(task_x, .compute = NULL) {
      sent[[length(sent) + 1L]] <<- task_x
      real_submit(task_x, .compute = .compute)
    },
    .package = "PredictProR"
  )
  out <- PredictProR:::sp_apply(as.list(1:4), function(i) i + length(big), mirai_decision(workers = 2L),
                                packages = character(), seed = TRUE)
  expect_equal(unlist(out), 1:4 + 2e6)
  expect_length(sent, 4L)
  expect_lt(max(vapply(sent, function(x) length(serialize(x, NULL)), numeric(1))), 1e4)
})

test_that("each GPU daemon sees only its own GPU, and every profile is stopped", {
  skip_if_no_mirai()
  withr::local_envvar(GP_CUDA_VISIBLE_DEVICES = NA, CUDA_VISIBLE_DEVICES = NA)
  out <- PredictProR:::sp_apply(
    as.list(1:4),
    function(i) list(cuda = Sys.getenv("CUDA_VISIBLE_DEVICES"), pid = Sys.getpid()),
    mirai_decision(workers = 2L, gpu_ids = c("0", "1")),
    packages = character(),
    seed = FALSE
  )
  expect_identical(vapply(out, function(r) r$cuda, character(1)), c("0", "1", "0", "1"))
  expect_false(any(vapply(out, function(r) as.numeric(r$pid), numeric(1)) == Sys.getpid()))
  expect_identical(Sys.getenv("CUDA_VISIBLE_DEVICES"), "")
  for (profile in c("default", "gpu1", "gpu2")) {
    expect_false(mirai::daemons_set(.compute = profile), info = profile)
  }
})

test_that("a daemon that fails its set-up stops the run", {
  skip_if_no_mirai()
  real_env <- PredictProR:::gp_parallel_runtime_env
  local_mocked_bindings(
    gp_parallel_runtime_env = function() {
      env <- real_env()
      env[["PREDICTPRO_EXPECTED_VERSION"]] <- "0.0.1"
      env
    },
    .package = "PredictProR"
  )
  expect_error(
    PredictProR:::sp_apply(as.list(1:2), function(i) i, mirai_decision(workers = 2L),
                           packages = character(), seed = FALSE),
    "this session runs 0.0.1"
  )
  expect_false(mirai::daemons_set())
})

test_that("a task error is logged and returns NULL, or stops the run with fail_fast", {
  skip_if_no_mirai()
  task <- function(i) if (i == 2L) stop("task two failed") else i
  out <- PredictProR:::sp_apply(as.list(1:3), task, mirai_decision(workers = 2L),
                                packages = character(), seed = FALSE)
  expect_identical(out, list(1L, NULL, 3L))
  expect_error(
    PredictProR:::sp_apply(as.list(1:3), task, mirai_decision(workers = 2L),
                           packages = character(), seed = FALSE, fail_fast = TRUE),
    "fail_fast"
  )
  expect_false(mirai::daemons_set())
})

test_that("mirai submits the most expensive tasks first and returns results in the original order", {
  skip_if_no_mirai()
  submitted <- integer()
  local_mocked_bindings(
    start_metrics_endpoint = function(...) NULL,
    stop_metrics_endpoint = function(...) invisible(NULL),
    gp_mirai_stop_daemons = function(...) invisible(NULL),
    gp_mirai_start_daemons = function(...) invisible(NULL),
    gp_mirai_init_daemons = function(...) invisible(NULL),
    gp_mirai_submit = function(task_x, .compute = NULL) {
      submitted <<- c(submitted, task_x)
      list(value = task_x * 10L, profile = .compute)
    },
    gp_mirai_collect = function(x) x,
    gp_mirai_is_error_value = function(x) FALSE,
    gp_mirai_queue_depth = function() 0L,
    .package = "PredictProR"
  )
  decision <- mirai_decision(workers = 2L, gpu_ids = c("0", "1"))
  decision$task_seconds <- c(1, 30, 5, 60)
  out <- PredictProR:::sp_apply(as.list(1:4), identity, decision, packages = character(), seed = FALSE)
  expect_identical(submitted, c(4L, 2L, 3L, 1L))
  expect_identical(vapply(out, `[[`, integer(1), "value"), c(10L, 20L, 30L, 40L))
  # the two heaviest go to different GPUs
  expect_identical(vapply(out, `[[`, character(1), "profile"), c("gpu2", "gpu2", "gpu1", "gpu1"))
})

test_that("the dispatcher queue cap is given to mirai in MB", {
  captured <- NULL
  local_mocked_bindings(
    gp_mirai_queue_memory_mb = function(...) 64,
    gp_mirai_start_daemons = function(n, .compute = NULL, memory = NULL) {
      captured <<- memory
      stop("stop after capture")
    },
    gp_mirai_stop_daemons = function(...) invisible(NULL),
    .package = "PredictProR"
  )
  expect_error(
    PredictProR:::sp_apply(as.list(1:2), identity, mirai_decision(workers = 2L),
                           packages = character(), seed = FALSE),
    "stop after capture"
  )
  expect_identical(captured, 64)
})

test_that("GPU ids cover every GPU", {
  withr::local_envvar(GP_CUDA_VISIBLE_DEVICES = NA)
  expect_identical(PredictProR:::gp_parallel_gpu_ids(4L), c("0", "1", "2", "3"))
  expect_identical(PredictProR:::gp_parallel_gpu_ids(1L), "0")
  expect_identical(PredictProR:::gp_parallel_gpu_ids(0L), "0")
  withr::local_envvar(GP_CUDA_VISIBLE_DEVICES = "2:3")
  expect_identical(PredictProR:::gp_parallel_gpu_ids(4L), c("2", "3"))
})
