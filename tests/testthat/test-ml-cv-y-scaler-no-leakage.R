# TDD (audit finding): the classical-ML CV functions fit the response (y)
# centring/scaling on the FULL response vector before the train/test split, so
# held-out test responses leaked into the scaler used for the training fit and
# the reverse transform. gp_ml_cv_train_y_scaler must fit on training rows only.

.resolve_fn_yscaler <- function(name) {
  if (exists(name, mode = "function")) return(get(name))
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(tryCatch(getFromNamespace(name, "PredictProR"), error = function(e) NULL))
  }
  NULL
}

test_that("CV response scaler is fit on training rows only (no test-response leakage)", {
  fn <- .resolve_fn_yscaler("gp_ml_cv_train_y_scaler")
  if (is.null(fn)) skip("gp_ml_cv_train_y_scaler not available")
  set.seed(1)
  # Training rows 1:30 ~ N(0,1); held-out test rows 31:40 ~ N(1000,1). A leaky
  # full-vector fit would centre near ~250; the train-only fit centres near 0.
  y <- c(rnorm(30, 0, 1), rnorm(10, 1000, 1))
  tst <- 31:40
  sc <- fn(y, tst)

  expect_lt(abs(as.numeric(sc$mean[[1]])), 1)   # train mean ~0, NOT ~250 (no leakage)
  expect_lt(as.numeric(sc$std[[1]]), 3)         # train sd ~1, NOT ~480

  # The scaler must round-trip training values (apply then reverse = identity).
  scaled <- stats::predict(sc, as.data.frame(as.matrix(y)))[, 1]
  reverted <- scaled * as.numeric(sc$std[[1]]) + as.numeric(sc$mean[[1]])
  expect_equal(reverted[1:5], y[1:5], tolerance = 1e-8)
})

test_that("CV response scaler tolerates an empty test set", {
  fn <- .resolve_fn_yscaler("gp_ml_cv_train_y_scaler")
  if (is.null(fn)) skip("gp_ml_cv_train_y_scaler not available")
  set.seed(2)
  y <- rnorm(20, 5, 2)
  sc <- fn(y, integer(0))
  expect_equal(as.numeric(sc$mean[[1]]), mean(y), tolerance = 1e-8)
})
