test_that("grid tuner selects lower-loss parameter setting for gaussian tasks", {
  y <- c(1, 2, 3, 4, 5, 6)

  res <- PredictProR:::gp_grid_tune_cv(
    y = y,
    response_family = "gaussian",
    param_grid = list(alpha = c(1, 2)),
    nfolds = 3L,
    random_state = 1L,
    predict_fun = function(tst, params) {
      pred <- if (params$alpha == 1) y[tst] else rep(mean(y[-tst]), length(tst))
      list(pred = pred, prob = NULL)
    }
  )

  expect_identical(res$best_params$alpha, 1)
  expect_identical(res$metric, "root_mean_squared_error")
})

test_that("python tabular model tuning forwards best parameter set", {
  skip_if_no_ml_python()  # the full-data fit runs through the real bridge
  pheno <-data.frame(GID = c("g1", "g2", "g3", "g4"), Yield = c(1, 2, 3, 4))
  geno <- matrix(1:16, nrow = 4, dimnames = list(pheno$GID, paste0("m", 1:4)))
  seen <- list()

  local_mocked_bindings(
    gp_grid_tune_cv = function(...) {
      list(best_params = list(ridge_alpha = 0.5), metric = "root_mean_squared_error", best_score = 0.1)
    },
    gp_py_ml_bootstrap_matrix = function(model_type, X_train, y_train, X_pred, response_family,
                                         model_params, class_levels, n_bootstrap, seed, prefer_gpu) {
      seen$model_params <<- model_params
      matrix(0, nrow = as.integer(n_bootstrap), ncol = nrow(X_pred))
    },
    diagnostic_plot_true_prediction = function(...) NULL,
    .package = "PredictProR"
  )

  out <- PredictProR:::AI_RidgeRegression_Lasso(
    pheno_object = pheno,
    geno_omic_object = geno,
    response = "Yield",
    gen_name = "GID",
    GS_model = "Ridge_Regression",
    para_tunning = TRUE,
    n_bootstrap = 2L,
    scaling = FALSE,
    centering = FALSE
  )

  expect_true("tune_metric" %in% out$model_parameters$stat)
  expect_identical(seen$model_params$ridge_alpha, 0.5)
})

test_that("auto-tuned Lasso bootstrap reuses the full fit's penalty; calibration folds search their own", {
  set.seed(9)
  geno <- matrix(rnorm(80), nrow = 10, dimnames = list(paste0("g", 1:10), paste0("m", 1:8)))
  pheno <- data.frame(GID = rownames(geno), Yield = rnorm(10))
  seen <- list(batch = list(), boot = NULL)
  local_mocked_bindings(
    gp_py_ml_fit_predict = function(model_type, X_train, y_train, X_test, model_params = list(), ...) {
      out <- rep(0, nrow(X_test))
      attr(out, "selected_alpha") <- 0.042
      out
    },
    gp_py_ml_fit_predict_batch = function(jobs) {
      seen$batch <<- lapply(jobs, function(job) job$model_params)
      lapply(jobs, function(job) rep(0, nrow(job$build()$X_test)))
    },
    gp_py_ml_bootstrap_matrix = function(model_type, X_train, y_train, X_pred, response_family,
                                         model_params, class_levels, n_bootstrap, seed, prefer_gpu) {
      seen$boot <<- model_params
      matrix(0, nrow = as.integer(n_bootstrap), ncol = nrow(X_pred))
    },
    diagnostic_plot_true_prediction = function(...) NULL,
    .package = "PredictProR"
  )
  PredictProR:::AI_RidgeRegression_Lasso(
    pheno_object = pheno, geno_omic_object = geno, response = "Yield", gen_name = "GID",
    GS_model = "Lasso", n_bootstrap = 2L, scaling = FALSE, centering = FALSE
  )
  expect_identical(seen$boot$lasso_alpha, 0.042)
  expect_true(length(seen$batch) > 0L)
  expect_true(all(vapply(seen$batch, function(p) is.null(p$lasso_alpha), logical(1))))
})

test_that("ML bridge payload carries a CV-selected penalty through to the prediction", {
  d <- withr::local_tempdir()
  writeLines(c("prediction", "1.5", "2.5"), file.path(d, "predictions.csv"))
  writeLines('{"classes": [], "task": "gaussian", "selected_alpha": 0.0123}', file.path(d, "meta.json"))
  res <- PredictProR:::gp_ml_read_prediction_payload(d)
  expect_equal(res$selected_alpha, 0.0123)
  pred <- PredictProR:::gp_py_ml_format_prediction(res, response_family = "gaussian")
  expect_equal(as.numeric(pred), c(1.5, 2.5))
  expect_equal(attr(pred, "selected_alpha"), 0.0123)
})

test_that("blend selection avoids full metric recomputation for every rmse alpha", {
  calls <- 0L
  metric_stub <- function(prediction, observed, env = NULL) {
    calls <<- calls + 1L
    prediction <- as.numeric(prediction)
    observed <- as.numeric(observed)
    list(
      rmse = sqrt(mean((prediction - observed)^2)),
      mae = mean(abs(prediction - observed)),
      cor = NA_real_,
      centered_rmse = NA_real_,
      centered_mae = NA_real_,
      centered_cor = NA_real_,
      within_env_spearman = NA_real_
    )
  }

  local_mocked_bindings(
    gp_tuning_score_vectors = metric_stub,
    .package = "PredictProR"
  )

  alpha_grid <- seq(0, 1, by = 0.1)
  gp_prediction <- c(1.0, 2.1, 2.9, 4.2)
  baseline_prediction <- c(2.0, 2.0, 2.0, 2.0)
  observed <- c(1, 2, 3, 4)

  out <- PredictProR:::gp_tuning_select_blend(
    gp_prediction = gp_prediction,
    baseline_prediction = baseline_prediction,
    observed = observed,
    objective = "rmse",
    alpha_grid = alpha_grid
  )

  brute_force <- vapply(alpha_grid, function(alpha) {
    pred <- alpha * gp_prediction + (1 - alpha) * baseline_prediction
    sqrt(mean((pred - observed)^2))
  }, numeric(1L))

  expect_lte(calls, 4L)
  expect_equal(out$blend_rmse, min(brute_force), tolerance = 1e-12)
  expect_equal(out$blend_alpha, alpha_grid[which.min(brute_force)])
})

test_that("candidate row indexing preserves batch candidate slices", {
  rows <- data.frame(
    candidate_id = c(2L, 1L, 2L, 3L),
    value = c("two-a", "one", "two-b", "three"),
    stringsAsFactors = FALSE
  )

  by_id <- PredictProR:::gp_tuning_rows_by_candidate_id(rows)

  expect_identical(by_id[["1"]]$value, "one")
  expect_identical(by_id[["2"]]$value, c("two-a", "two-b"))
  expect_identical(by_id[["3"]]$value, "three")
})
