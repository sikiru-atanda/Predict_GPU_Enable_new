# Lasso runs in R with glmnet rather than through the Python bridge.
#
# scikit-learn's LassoCV needed 82 min for ONE penalty search on 1,451 mice x
# 10,346 markers; cv.glmnet (strong-rule screening, Fortran) took 15 s and
# chose the same penalty (lambda 0.040 vs alpha 0.041) with the same held-out
# accuracy (0.325 vs 0.320). glmnet's lambda and scikit-learn's alpha share
# one scale: minimise (1/2n)||y - b0 - Xb||^2 + lambda * ||b||_1.
#
# Inputs arrive already scaled by the ML preprocessor, so glmnet is called with
# standardize = FALSE and an intercept, matching the former sklearn Lasso.
# `lasso_alpha` in model_params fixes the penalty (tuning grids, bootstrap
# refits); without it the penalty is chosen by 5-fold cv.glmnet (lambda.min).

gp_lasso_seed <- function(model_params) {
  seed <- suppressWarnings(as.numeric(model_params$random_state %||% gp_ml_random_state())[1L])
  if (!is.finite(seed)) seed <- 42
  seed
}

gp_lasso_check_gaussian <- function(response_family) {
  if (!identical(gp_resolve_response_family(response_family), "gaussian")) {
    stop("Lasso only supports gaussian regression.", call. = FALSE)
  }
  invisible(TRUE)
}

gp_lasso_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- 0
  unname(x)
}

# Fit at a fixed penalty or choose it by CV; returns the fit and its lambda.
gp_lasso_fit <- function(X_train, y_train, model_params = list()) {
  x <- gp_lasso_matrix(X_train)
  y <- as.numeric(y_train)
  alpha <- suppressWarnings(as.numeric(model_params$lasso_alpha %||% NA_real_)[1L])
  if (is.finite(alpha)) {
    fit <- glmnet::glmnet(x, y, family = "gaussian", alpha = 1, lambda = alpha, standardize = FALSE)
    return(list(fit = fit, lambda = alpha))
  }
  nfolds <- max(3L, min(5L, length(y)))
  foldid <- gp_with_pinned_seed(gp_lasso_seed(model_params), sample(rep(seq_len(nfolds), length.out = length(y))))
  cv <- glmnet::cv.glmnet(x, y, family = "gaussian", alpha = 1, foldid = foldid, standardize = FALSE)
  list(fit = cv$glmnet.fit, lambda = cv$lambda.min)
}

# Same return shape as the Python bridge's fit-predict payload.
gp_lasso_fit_predict <- function(X_train, y_train, X_test, response_family = "gaussian", model_params = list()) {
  gp_lasso_check_gaussian(response_family)
  res <- gp_lasso_fit(X_train, y_train, model_params)
  pred <- as.numeric(stats::predict(res$fit, newx = gp_lasso_matrix(X_test), s = res$lambda, exact = FALSE))
  list(predictions = pred, probabilities = NULL, classes = NULL, selected_alpha = res$lambda)
}

# n_bootstrap x n_pred matrix of refits on resampled training rows, at the
# penalty in model_params$lasso_alpha (chosen once by CV when absent).
gp_lasso_bootstrap_matrix <- function(X_train, y_train, X_pred, response_family = "gaussian",
                                      model_params = list(), n_bootstrap = 30L, seed = NULL) {
  gp_lasso_check_gaussian(response_family)
  if (!is.finite(suppressWarnings(as.numeric(model_params$lasso_alpha %||% NA_real_)[1L]))) {
    model_params$lasso_alpha <- gp_lasso_fit(X_train, y_train, model_params)$lambda
  }
  x <- gp_lasso_matrix(X_train)
  y <- as.numeric(y_train)
  xp <- gp_lasso_matrix(X_pred)
  n <- nrow(x)
  idx <- gp_with_pinned_seed(seed %||% gp_lasso_seed(model_params),
                             lapply(seq_len(as.integer(n_bootstrap)), function(i) sample.int(n, n, replace = TRUE)))
  out <- matrix(NA_real_, nrow = length(idx), ncol = nrow(xp))
  for (b in seq_along(idx)) {
    fit <- glmnet::glmnet(x[idx[[b]], , drop = FALSE], y[idx[[b]]], family = "gaussian", alpha = 1,
                          lambda = model_params$lasso_alpha, standardize = FALSE)
    out[b, ] <- as.numeric(stats::predict(fit, newx = xp))
  }
  out
}

# |coefficient| at the chosen penalty, one value per predictor column.
gp_lasso_feature_importance <- function(X_train, y_train, response_family = "gaussian", model_params = list()) {
  gp_lasso_check_gaussian(response_family)
  res <- gp_lasso_fit(X_train, y_train, model_params)
  abs(as.numeric(stats::coef(res$fit, s = res$lambda, exact = FALSE))[-1L])
}
