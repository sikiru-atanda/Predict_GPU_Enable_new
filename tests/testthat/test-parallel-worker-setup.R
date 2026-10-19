# Worker set-up cost and fork safety.
# - base_parallel on Unix forked once per task (mclapply(mc.preschedule =
#   FALSE)), so every task re-ran the worker set-up and the initializer;
# - the set-up started embedded Python whenever PREDICTPRO_PYTHON was set
#   (about 4 s per task on Linux);
# - a fork after mirai daemons, in a session that had forked before, crashed
#   the forked workers.
# The fork tests run in a fresh R process: earlier tests in this session may
# have run mirai daemons, after which sp_apply() no longer forks.

run_fresh <- function(fun) {
  callr::r(fun, libpath = .libPaths(), env = c(callr::rcmd_safe_env(), PREDICTPRO_IN_PARALLEL_WORKER = ""))
}

test_that("base_parallel forks each worker once and runs the initializer once per worker", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  res <- run_fresh(function() {
    suppressPackageStartupMessages(library(PredictProR))
    logger::log_threshold(logger::WARN)
    marks <- tempfile(); dir.create(marks)
    out <- PredictProR:::sp_apply(
      as.list(1:12),
      function(i) list(i = i, pid = Sys.getpid(), worker = Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER")),
      list(backend = "base_parallel", workers = 2L, prefer_fork = TRUE),
      packages = character(), seed = FALSE,
      initializer = function() file.create(file.path(marks, Sys.getpid()))
    )
    list(out = out, marks = list.files(marks), main_pid = Sys.getpid(),
         main_worker_flag = Sys.getenv("PREDICTPRO_IN_PARALLEL_WORKER"),
         fork_env_left = length(ls(PredictProR:::gp_parallel_fork_env)))
  })
  pids <- unique(vapply(res$out, function(r) as.numeric(r$pid), numeric(1)))
  expect_identical(vapply(res$out, function(r) r$i, integer(1)), 1:12)
  expect_lte(length(pids), 2L)
  expect_false(any(pids == res$main_pid))
  expect_true(all(vapply(res$out, function(r) r$worker, character(1)) == "1"))
  expect_setequal(res$marks, as.character(pids))
  expect_false(identical(res$main_worker_flag, "1"))
  expect_identical(res$fork_env_left, 0L)
})

test_that("base_parallel fork results equal the sequential ones, heaviest tasks first", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  same <- run_fresh(function() {
    suppressPackageStartupMessages(library(PredictProR))
    logger::log_threshold(logger::WARN)
    task <- function(i) c(i, stats::runif(2))
    seq_out <- PredictProR:::sp_apply(as.list(1:6), task, list(backend = "sequential", workers = 1L),
                                      packages = character(), seed = 7L)
    decision <- list(backend = "base_parallel", workers = 3L, prefer_fork = TRUE,
                     task_seconds = c(1, 9, 2, 8, 3, 7))
    fork_out <- PredictProR:::sp_apply(as.list(1:6), task, decision, packages = character(), seed = 7L)
    identical(fork_out, seq_out)
  })
  expect_true(same)
})

test_that("after mirai daemons, base_parallel uses PSOCK workers instead of crashing forked ones", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  skip_if_not_installed("mirai", "2.7.0")
  skip_if_not_installed("purrr")
  res <- run_fresh(function() {
    suppressPackageStartupMessages(library(PredictProR))
    logger::log_threshold(logger::WARN)
    f <- function(i) { s <- 0; for (k in 1:200) s <- s + sum(sqrt(seq_len(2e4) + k)); i }
    d <- function(b) list(backend = b, workers = 4L, prefer_fork = TRUE, num_gpus = 0L,
                          gpu_parallel_enabled = FALSE, gpu_ids = character())
    run <- function(b) tryCatch(unlist(PredictProR:::sp_apply(as.list(1:8), f, d(b), packages = character())),
                                error = function(e) conditionMessage(e))
    before <- PredictProR:::gp_fork_safe_after_mirai()
    a <- run("base_parallel"); b <- run("mirai"); c <- run("base_parallel")
    list(before = before, after = PredictProR:::gp_fork_safe_after_mirai(), a = a, b = b, c = c,
         prefer_fork_after = PredictProR:::gp_parallel_safe_fork_preference(TRUE, FALSE))
  })
  expect_true(res$before)
  expect_false(res$after)
  expect_false(res$prefer_fork_after)
  for (step in c("a", "b", "c")) expect_identical(res[[step]], 1:8, info = step)
})

# Every worker loaded reticulate (0.58 s on Windows, 0.43 s of it Matrix),
# while BGLR tasks never use it and code that embeds Python loads it itself.
test_that("PSOCK workers load reticulate only when a task uses it", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  decision <- list(backend = "base_parallel", workers = 2L, prefer_fork = FALSE)
  loaded <- function(i) c(reticulate = isNamespaceLoaded("reticulate"), logger = isNamespaceLoaded("logger"))
  out <- PredictProR:::sp_apply(as.list(1:2), loaded, decision, packages = character(), seed = FALSE)
  expect_false(any(vapply(out, `[[`, logical(1), "reticulate")))
  expect_true(all(vapply(out, `[[`, logical(1), "logger")))
  used <- PredictProR:::sp_apply(as.list(1:2), function(i) is.function(reticulate::py_available),
                                 decision, packages = character(), seed = FALSE)
  expect_true(all(unlist(used)))
})

test_that("parallel workers do not start embedded Python because PREDICTPRO_PYTHON is set", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  skip_if(reticulate::py_available(initialize = FALSE), "Python already embedded in this session")
  py <- Sys.which(c("python3", "python"))
  py <- py[nzchar(py)]
  skip_if(!length(py), "no Python on PATH")
  withr::local_envvar(PREDICTPRO_PYTHON = unname(py[[1L]]))
  backends <- c(if (.Platform$OS.type == "unix") "base_parallel",
                if (requireNamespace("mirai", quietly = TRUE)) "mirai")
  for (backend in backends) {
    decision <- list(backend = backend, workers = 2L, prefer_fork = TRUE,
                     num_gpus = 0L, gpu_parallel_enabled = FALSE, gpu_ids = character())
    out <- PredictProR:::sp_apply(as.list(1:4), function(i) reticulate::py_available(initialize = FALSE),
                                  decision, packages = character(), seed = FALSE)
    expect_false(any(unlist(out)), info = backend)
  }
})
