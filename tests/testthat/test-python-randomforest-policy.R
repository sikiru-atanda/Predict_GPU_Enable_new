test_that("RandomForest policy treats Python RF internal jobs as sequential work", {
  params <- gp_parallel_extract_model_policy_params(
    "RandomForest",
    list(n_jobs = 4L, ntree = 100L)
  )
  cap <- gp_parallel_model_capability_table(
    models = "RandomForest",
    model_params_list = list(params),
    num_gpus = 0L
  )

  expect_equal(params$n_jobs, 4L)
  expect_true(cap$thread_sensitive)
  expect_true(cap$forced_sequential)
  expect_match(cap$reasons, "randomforest_internal_jobs")
})

test_that("RandomForest fit-predict routes through Python bridge", {
  captured <- NULL

  local_mocked_bindings(
    gp_ml_cli_backend = function() "cli",
    gp_ml_fit_predict_cli = function(model_type,
                                     X_train,
                                     y_train,
                                     X_test,
                                     response_family = "gaussian",
                                     model_params = list(),
                                     class_levels = NULL) {
      captured <<- list(
        model_type = model_type,
        X_train = X_train,
        y_train = y_train,
        X_test = X_test,
        response_family = response_family,
        model_params = model_params
      )
      list(predictions = rep(0, nrow(X_test)), probabilities = NULL, classes = NULL)
    },
    .package = "PredictProR"
  )

  x <- matrix(seq_len(12), nrow = 4)
  pred <- PredictProR:::gp_py_ml_fit_predict(
    model_type = "randomforest",
    X_train = x,
    y_train = c(1, NA, 3, 4),
    X_test = matrix(seq_len(6), nrow = 2),
    response_family = "gaussian",
    model_params = list(ntree = 20, n_jobs = 1L)
  )

  expect_length(pred, 2)
  expect_identical(captured$model_type, "randomforest")
  expect_false(anyNA(captured$y_train))
  expect_identical(as.numeric(captured$y_train), c(1, 3, 4))
  expect_equal(captured$X_train, x[c(1, 3, 4), , drop = FALSE])
  # n_jobs = 1 (the default) becomes this process's share of the cores
  expect_equal(captured$model_params$n_jobs, PredictProR:::gp_ml_thread_share())
  withr::local_envvar(PREDICTPRO_ML_AUTO_THREADS = "false")
  PredictProR:::gp_py_ml_fit_predict(
    model_type = "randomforest", X_train = x, y_train = c(1, NA, 3, 4),
    X_test = matrix(seq_len(6), nrow = 2), response_family = "gaussian",
    model_params = list(ntree = 20, n_jobs = 1L)
  )
  expect_equal(captured$model_params$n_jobs, 1L)
})

test_that("RandomForest batch jobs stay on Python batch bridge", {
  captured <- NULL

  local_mocked_bindings(
    gp_ml_fit_predict_batch_cli = function(jobs) {
      captured <<- jobs
      lapply(jobs, function(job) {
        list(predictions = rep(0, nrow(job$X_test)), probabilities = NULL, classes = NULL)
      })
    },
    .package = "PredictProR"
  )

  x <- matrix(seq_len(15), nrow = 5)
  out <- PredictProR:::gp_py_ml_fit_predict_batch(list(list(
    model_type = "randomforest",
    X_train = x,
    y_train = c(1, NA, 3, NA, 5),
    X_test = matrix(seq_len(6), nrow = 2),
    response_family = "gaussian",
    model_params = list(ntree = 20, n_jobs = 1L)
  )))

  expect_length(out, 1)
  expect_identical(captured[[1]]$model_type, "randomforest")
  expect_false(anyNA(captured[[1]]$y_train))
  expect_identical(as.numeric(captured[[1]]$y_train), c(1, 3, 5))
  expect_equal(captured[[1]]$X_train, x[c(1, 3, 5), , drop = FALSE])
})
