# 0.20.220: the ridge penalty grid no longer stops at 1e3, settings reach Python at full precision,
# and the ML bridge script is the package's own copy. Evidence: NEWS 0.20.220.

test_that("settings written for the Python bridges keep full precision", {
  x <- list(alpha = 1.2345e-4, small = 8.7096e-6, lr = 0.123456789, grid = c(1e-6, 2.5e-5, 7.25e-3))
  for (writer in list(gp_ml_write_json, gp_dl_write_json, gp_bridge_write_json)) {
    path <- withr::local_tempfile(fileext = ".json")
    writer(x, path)
    back <- jsonlite::fromJSON(path)
    expect_equal(back$alpha, x$alpha, tolerance = 0)
    expect_equal(back$small, x$small, tolerance = 0)
    expect_equal(back$lr, x$lr, tolerance = 0)
    expect_equal(back$grid, x$grid, tolerance = 0)
  }
})

test_that("the ML bridge script is the package's copy, not one in the working directory", {
  decoy_root <- withr::local_tempdir()
  for (sub in c("inst/python", "python")) {
    dir.create(file.path(decoy_root, sub), recursive = TRUE)
    writeLines("raise SystemExit('decoy')", file.path(decoy_root, sub, "ml_bridge.py"))
  }
  withr::local_dir(decoy_root)
  path <- gp_ml_bridge_script_path()
  expect_identical(path, normalizePath(system.file("python/ml_bridge.py", package = "PredictProR"),
                                       winslash = "/", mustWork = TRUE))
  expect_false(startsWith(path, normalizePath(decoy_root, winslash = "/")))
})

test_that("ridge's default penalty search reaches past the old 1e3 cap", {
  skip_if_no_ml_python()
  # Pure-noise trait on 2,000 standardized predictors: the best ridge penalty is as large as the
  # grid allows. The old grid stopped at 1e3, the value it chose on every simulated marker panel
  # of 2,000 or more columns.
  withr::local_seed(20261016)
  n <- 120; p <- 2000
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  pred <- gp_py_ml_fit_predict("ridge", X_train = X[1:100, ], y_train = y[1:100],
                               X_test = X[101:120, ], response_family = "gaussian")
  alpha <- attr(pred, "selected_alpha", exact = TRUE)
  expect_true(is.numeric(alpha) && length(alpha) == 1L)
  expect_gt(alpha, 1e3)
  expect_length(pred, 20L)
  expect_true(all(is.finite(pred)))
})

test_that("a user-supplied ridge grid with small values reaches Python intact", {
  skip_if_no_ml_python()
  withr::local_seed(7)
  X <- matrix(rnorm(60 * 5), 60, 5)
  y <- drop(X %*% c(1, -1, 0.5, 0, 0)) + rnorm(60, sd = 0.1)
  # 10^seq(-6, 0, 25) holds values such as 1.78e-5 that the old JSON writer wrote as 0, which
  # scikit-learn rejects ("alphas[...] == 0.0, must be > 0.0")
  pred <- gp_py_ml_fit_predict("ridge", X_train = X[1:50, ], y_train = y[1:50], X_test = X[51:60, ],
                               response_family = "gaussian",
                               model_params = list(ridge_alphas = 10^seq(-6, 0, length.out = 25)))
  expect_length(pred, 10L)
  expect_true(all(is.finite(pred)))
})
