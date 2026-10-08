test_that("policy calibration benchmark writes report tables", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  script <- file.path(root, "tools", "benchmark_parallel_policy_calibration.R")
  if (!file.exists(script)) {
    testthat::skip("policy calibration benchmark helper is not included in source-tarball checks")
  }

  env <- new.env(parent = globalenv())
  source(script, local = env)

  out_dir <- tempfile("policy-calibration-")
  res <- env$run_policy_calibration_benchmark(
    out_dir = out_dir,
    workloads = c("tiny_cpu"),
    backends = c("sequential"),
    n_tasks = 2L,
    repeats = 1L,
    workers = 1L,
    sleep_seconds = 0,
    compute_size = 50L,
    root = root
  )

  expect_true(file.exists(res$timings_file))
  expect_true(file.exists(res$summary_file))
  expect_true(file.exists(res$decisions_file))

  timings <- utils::read.csv(res$timings_file, stringsAsFactors = FALSE)
  summary <- utils::read.csv(res$summary_file, stringsAsFactors = FALSE)
  decisions <- utils::read.csv(res$decisions_file, stringsAsFactors = FALSE)

  expect_true(all(c("workload", "backend", "elapsed_seconds", "replicate") %in% names(timings)))
  expect_true(all(c("workload", "fastest_backend", "auto_backend", "auto_matches_fastest") %in% names(summary)))
  expect_true(all(c("workload", "auto_backend", "decision_reason") %in% names(decisions)))
  expect_identical(summary$workload, "tiny_cpu")
  expect_identical(summary$fastest_backend, "sequential")
})
