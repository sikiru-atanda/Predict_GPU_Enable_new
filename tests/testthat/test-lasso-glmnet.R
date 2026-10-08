# Lasso is fitted in R by glmnet (sklearn's LassoCV took 82 min per penalty
# search at mice scale; cv.glmnet 15 s with the same penalty and accuracy).
make_lasso_data <- function(n = 120L, p = 60L, seed = 3L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n)
  y <- as.numeric(X[, 1:4] %*% c(1.5, -1, 0.8, 0.5) + rnorm(n, sd = 0.5))
  list(X = X, y = y, X_test = matrix(rnorm(20L * p), 20L))
}

test_that("fixed-penalty Lasso equals a direct glmnet fit", {
  skip_if_not_installed("glmnet")
  d <- make_lasso_data()
  res <- PredictProR:::gp_lasso_fit_predict(d$X, d$y, d$X_test, model_params = list(lasso_alpha = 0.05))
  ref <- glmnet::glmnet(d$X, d$y, alpha = 1, lambda = 0.05, standardize = FALSE)
  expect_equal(res$predictions, as.numeric(predict(ref, newx = d$X_test)), tolerance = 1e-10)
  expect_identical(res$selected_alpha, 0.05)
})

test_that("CV Lasso reports its chosen penalty and is reproducible for a fixed seed", {
  skip_if_not_installed("glmnet")
  d <- make_lasso_data()
  a <- PredictProR:::gp_lasso_fit_predict(d$X, d$y, d$X_test, model_params = list(random_state = 7))
  b <- PredictProR:::gp_lasso_fit_predict(d$X, d$y, d$X_test, model_params = list(random_state = 7))
  expect_true(is.finite(a$selected_alpha) && a$selected_alpha > 0)
  expect_identical(a, b)
  # the true signal is recovered: held-in correlation well above noise
  expect_gt(cor(PredictProR:::gp_lasso_fit_predict(d$X, d$y, d$X, model_params = list(random_state = 7))$predictions, d$y), 0.8)
})

test_that("Lasso bootstrap refits at the given penalty and is reproducible", {
  skip_if_not_installed("glmnet")
  d <- make_lasso_data()
  b1 <- PredictProR:::gp_lasso_bootstrap_matrix(d$X, d$y, d$X_test, model_params = list(lasso_alpha = 0.05),
                                                n_bootstrap = 5L, seed = 11)
  b2 <- PredictProR:::gp_lasso_bootstrap_matrix(d$X, d$y, d$X_test, model_params = list(lasso_alpha = 0.05),
                                                n_bootstrap = 5L, seed = 11)
  expect_identical(dim(b1), c(5L, 20L))
  expect_identical(b1, b2)
  expect_true(all(apply(b1, 2, stats::sd) > 0))   # resamples differ
})

test_that("Lasso rejects classification and routes through glmnet, not Python", {
  skip_if_not_installed("glmnet")
  d <- make_lasso_data()
  expect_error(PredictProR:::gp_lasso_fit_predict(d$X, factor(d$y > 0), d$X_test, response_family = "binary"),
               "only supports gaussian")
  local_mocked_bindings(
    gp_ml_fit_predict_cli = function(...) stop("Python bridge must not be used for Lasso"),
    gp_ml_fit_predict_batch_cli = function(...) stop("Python bridge must not be used for Lasso"),
    .package = "PredictProR"
  )
  single <- PredictProR:::gp_py_ml_fit_predict("lasso", d$X, d$y, d$X_test, model_params = list(lasso_alpha = 0.05))
  expect_equal(as.numeric(single), PredictProR:::gp_lasso_fit_predict(d$X, d$y, d$X_test, model_params = list(lasso_alpha = 0.05))$predictions)
  expect_equal(attr(single, "selected_alpha"), 0.05)
  batch <- PredictProR:::gp_py_ml_fit_predict_batch(list(list(
    model_type = "lasso", y_train = d$y, response_family = "gaussian", model_params = list(lasso_alpha = 0.05),
    build = function() list(X_train = d$X, X_test = d$X_test)
  )))
  expect_equal(as.numeric(batch[[1L]]), as.numeric(single))
})
