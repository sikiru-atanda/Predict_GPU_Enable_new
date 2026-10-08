# Parallel workers must load the same PredictProR as the session. They used to
# resolve it from the default library, silently running an older copy (e.g.
# 0.20.152 workers under a 0.20.173 session, without the thread-share logic,
# so parallel CV ran ~5x slower than sequential).
test_that("worker runtime env puts the session's PredictProR library first", {
  lib <- withr::local_tempdir()
  local_mocked_bindings(gp_predictpror_loaded_lib = function() normalizePath(lib, winslash = "/"), .package = "PredictProR")
  env <- PredictProR:::gp_parallel_runtime_env()
  first <- strsplit(env[["R_LIBS"]], .Platform$path.sep, fixed = TRUE)[[1L]][1L]
  expect_identical(first, normalizePath(lib, winslash = "/"))
})

test_that("worker runtime env always carries the session's version, and the build stamp of an install", {
  env <- PredictProR:::gp_parallel_runtime_env()
  expect_identical(env[["PREDICTPRO_EXPECTED_VERSION"]], as.character(getNamespaceVersion("PredictProR")))
  built <- PredictProR:::gp_predictpror_identity()$built
  if (nzchar(built)) expect_identical(env[["PREDICTPRO_EXPECTED_BUILD"]], built)
})

test_that("a worker refuses to run a different PredictProR version", {
  current <- as.character(getNamespaceVersion("PredictProR"))
  expect_identical(PredictProR:::gp_parallel_check_worker_version(current, ""), unname(current))
  expect_identical(PredictProR:::gp_parallel_check_worker_version("", ""), unname(current))
  expect_error(PredictProR:::gp_parallel_check_worker_version("0.0.1", ""), "this session runs 0.0.1")
})

# Same version number is not enough: a reinstall with changed code keeps the
# number but gets a new build stamp.
test_that("a worker refuses a different install of the same version", {
  current <- as.character(getNamespaceVersion("PredictProR"))
  local_mocked_bindings(gp_predictpror_identity = function() {
    list(version = current, built = "R 4.3.3; ; 2026-01-01 00:00:00 UTC; windows", path = tempdir())
  }, .package = "PredictProR")
  expect_error(
    PredictProR:::gp_parallel_check_worker_version(current, "R 4.3.3; ; 2026-10-04 12:00:00 UTC; windows"),
    "different install of PredictProR"
  )
  expect_identical(
    PredictProR:::gp_parallel_check_worker_version(current, "R 4.3.3; ; 2026-01-01 00:00:00 UTC; windows"),
    current
  )
})

# Helpers that set their own future plan (feature selection, KNN imputation,
# VCF QC) run through gp_future_lapply_session: workers get the session's
# library and check the version before every task.
test_that("session future_lapply runs tasks and stops on a worker version mismatch", {
  skip_if_not_installed("future.apply")
  old_plan <- future::plan()
  withr::defer(future::plan(old_plan))
  future::plan(future::multisession, workers = 2L)

  out <- PredictProR:::gp_future_lapply_session(1:3, function(x) x * 2)
  expect_identical(unlist(out), c(2, 4, 6))

  local_mocked_bindings(gp_parallel_runtime_env = function() c(PREDICTPRO_EXPECTED_VERSION = "0.0.1"), .package = "PredictProR")
  expect_error(suppressWarnings(PredictProR:::gp_future_lapply_session(1:2, function(x) x)), "this session runs 0.0.1")
})

test_that("the session wrapper does not ship the caller's data", {
  big <- runif(1e6)  # 8 MB in the caller's frame
  task <- function(x) x
  environment(task) <- globalenv()
  wrapped <- PredictProR:::gp_session_checked_fun(task, c(A = "1"), Sys.getpid())
  expect_lt(length(serialize(wrapped, NULL)), 2e6)
  expect_identical(wrapped(3), 3)  # in-process: runs without touching the env
})
