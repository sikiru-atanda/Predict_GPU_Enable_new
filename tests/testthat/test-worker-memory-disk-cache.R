# The worker-memory check starts an R worker and loads PredictProR's imports
# (2.3-2.6 s on a Windows test machine, more with a Python child) in the
# first parallel run of every session. The measurement is kept on disk until
# PredictProR or one of its imports is reinstalled, R changes or the ML
# interpreter changes, for at most 30 days.

# A private cache file, a counted fake measurement, and a file that stands
# for an installed package's DESCRIPTION (reinstalling rewrites it).
local_memory_cache <- function(env = parent.frame()) {
  dir <- normalizePath(withr::local_tempdir(.local_envir = env), winslash = "/")
  mc <- new.env()
  mc$file <- file.path(dir, "python-discovery.rds")
  mc$desc <- file.path(dir, "DESCRIPTION")
  writeLines("Package: fake", mc$desc)
  mc$py <- file.path(dir, "python.exe")
  writeLines("fake interpreter", mc$py)
  mc$runs <- 0L
  mc$value <- 0.5
  local_mocked_bindings(
    gp_python_discovery_cache_file = function() mc$file,
    gp_parallel_worker_memory_paths = function(py_bin = NULL) c(mc$desc, py_bin),
    gp_parallel_measure_worker_gb_run = function(python_models, py_bin = NULL) {
      mc$runs <- mc$runs + 1L
      mc$value
    },
    gp_detect_ml_python = function() mc$py,
    .package = "PredictProR",
    .env = env
  )
  withr::local_envvar(PREDICTPRO_PROBE_CACHE = "true", .local_envir = env)
  cache <- PredictProR:::PredictProR_runtime_cache
  saved <- mget(ls(cache), envir = cache)
  withr::defer({
    rm(list = ls(cache), envir = cache)
    list2env(saved, envir = cache)
  }, envir = env)
  memory_new_session()
  mc
}

# What a new R session starts with: no measurement in memory.
memory_new_session <- function() {
  cache <- PredictProR:::PredictProR_runtime_cache
  rm(list = grep("^worker_memory_measured::", ls(cache), value = TRUE), envir = cache)
}

measure <- function(python = FALSE) PredictProR:::gp_parallel_measure_worker_gb(python)

test_that("the worker memory measured in one session is not measured again in the next", {
  mc <- local_memory_cache()
  expect_identical(measure(), 0.5)
  expect_identical(mc$runs, 1L)
  expect_true(file.exists(mc$file))
  memory_new_session()
  expect_identical(measure(), 0.5)
  expect_identical(mc$runs, 1L)
})

test_that("reinstalling a package makes the next session measure again", {
  mc <- local_memory_cache()
  measure()
  Sys.setFileTime(mc$desc, Sys.time() + 5)
  mc$value <- 0.7
  memory_new_session()
  expect_identical(measure(), 0.7)
  expect_identical(mc$runs, 2L)
})

test_that("R-only and Python measurements are kept apart", {
  mc <- local_memory_cache()
  mc$value <- 0.5
  measure(FALSE)
  mc$value <- 0.9
  measure(TRUE)
  expect_identical(mc$runs, 2L)
  memory_new_session()
  expect_identical(measure(FALSE), 0.5)
  expect_identical(measure(TRUE), 0.9)
  expect_identical(mc$runs, 2L)
})

test_that("another ML interpreter is measured again", {
  mc <- local_memory_cache()
  measure(TRUE)
  other <- file.path(dirname(mc$py), "python3.exe")
  writeLines("another interpreter", other)
  mc$py <- other
  memory_new_session()
  measure(TRUE)
  expect_identical(mc$runs, 2L)
})

test_that("a measurement made under another R or PredictProR is not used", {
  mc <- local_memory_cache()
  measure()
  entries <- readRDS(mc$file)
  key <- grep("^worker_memory::", names(entries), value = TRUE)
  expect_length(key, 1L)
  entries[[key]]$extra$r <- "R version 0.0.0"
  saveRDS(entries, mc$file)
  memory_new_session()
  measure()
  expect_identical(mc$runs, 2L)
})

test_that("a measurement older than 30 days is measured again", {
  mc <- local_memory_cache()
  measure()
  entries <- readRDS(mc$file)
  key <- grep("^worker_memory::", names(entries), value = TRUE)
  entries[[key]]$saved <- Sys.time() - 31 * 86400
  saveRDS(entries, mc$file)
  memory_new_session()
  measure()
  expect_identical(mc$runs, 2L)
})

test_that("a failed measurement is not kept on disk and is not retried in the same session", {
  mc <- local_memory_cache()
  mc$value <- NULL
  expect_null(measure())
  expect_null(measure())
  expect_identical(mc$runs, 1L)
  memory_new_session()
  mc$value <- 0.6
  expect_identical(measure(), 0.6)
  expect_identical(mc$runs, 2L)
})

test_that("with PREDICTPRO_PROBE_CACHE=false every session measures", {
  mc <- local_memory_cache()
  withr::local_envvar(PREDICTPRO_PROBE_CACHE = "false")
  measure()
  memory_new_session()
  measure()
  expect_identical(mc$runs, 2L)
  expect_false(file.exists(mc$file))
})

test_that("the paths a measurement depends on are PredictProR's and its imports' installs", {
  paths <- PredictProR:::gp_parallel_worker_memory_paths()
  # the loaded PredictProR's: an install, or the source tree under load_all()
  own <- file.path(normalizePath(find.package("PredictProR"), winslash = "/"), "DESCRIPTION")
  expect_identical(paths[[1L]], own)
  expect_true(any(grepl("/BGLR/DESCRIPTION$", paths)))
  expect_true(all(file.exists(paths)))
  py <- tempfile(fileext = ".exe")
  expect_identical(utils::tail(PredictProR:::gp_parallel_worker_memory_paths(py), 1L),
                   PredictProR:::gp_python_path(py, mustWork = FALSE))
})
