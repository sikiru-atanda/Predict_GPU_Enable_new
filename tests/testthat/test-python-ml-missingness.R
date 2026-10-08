test_that("ML fit-predict drops missing training labels before Python bridge", {
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
        X_train = X_train,
        y_train = y_train,
        X_test = X_test,
        response_family = response_family
      )
      list(predictions = rep(0, nrow(X_test)), probabilities = NULL, classes = NULL)
    },
    .package = "PredictProR"
  )

  X_train <- matrix(seq_len(12), nrow = 4)
  y_train <- c(1, NA, 3, 4)
  X_test <- matrix(seq_len(6), nrow = 2)

  out <- PredictProR:::gp_py_ml_fit_predict(
    model_type = "svm",
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    response_family = "gaussian"
  )

  expect_length(out, nrow(X_test))
  expect_false(anyNA(captured$y_train))
  expect_identical(as.numeric(captured$y_train), c(1, 3, 4))
  expect_equal(captured$X_train, X_train[c(1, 3, 4), , drop = FALSE])
})

test_that("ML batch fit-predict drops missing labels per job", {
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

  X_train <- matrix(seq_len(15), nrow = 5)
  jobs <- list(list(
    model_type = "svm",
    X_train = X_train,
    y_train = c(1, NA, 3, NA, 5),
    X_test = matrix(seq_len(6), nrow = 2),
    response_family = "gaussian"
  ))

  out <- PredictProR:::gp_py_ml_fit_predict_batch(jobs)

  expect_length(out, 1)
  expect_false(anyNA(captured[[1]]$y_train))
  expect_identical(as.numeric(captured[[1]]$y_train), c(1, 3, 5))
  expect_equal(captured[[1]]$X_train, X_train[c(1, 3, 5), , drop = FALSE])
})

test_that("lazy ML batch jobs write the same payload as eager jobs, one job at a time", {
  withr::local_envvar(c(PREDICTPRO_ML_PAYLOAD_CACHE = "false"))
  written <- list()
  local_mocked_bindings(
    gp_ml_run_cli = function(command, args, python_bin = NULL) {
      manifest <- jsonlite::fromJSON(args[[2]], simplifyVector = FALSE)
      for (job in manifest$jobs) {
        written[[length(written) + 1L]] <<- lapply(
          job[c("x_train_csv", "y_train_csv", "x_test_csv")], readLines
        )
        writeLines(c("prediction", rep("0", length(readLines(job$x_test_csv)) - 1L)),
                   file.path(job$out_dir, "predictions.csv"))
      }
      invisible(character())
    },
    .package = "PredictProR"
  )
  set.seed(4)
  x <- matrix(rnorm(40), 10, dimnames = list(paste0("g", 1:10), paste0("m", 1:4)))
  y <- c(1, NA, 3, 4, 5, 6, NA, 8, 9, 10)
  tst <- c(2, 5)
  eager <- list(model_type = "ridge", X_train = x[-tst, , drop = FALSE], y_train = y[-tst],
                X_test = x[tst, , drop = FALSE], response_family = "gaussian")
  lazy <- list(model_type = "ridge", y_train = y[-tst], response_family = "gaussian",
               build = function() list(X_train = x[-tst, , drop = FALSE], X_test = x[tst, , drop = FALSE]))

  out <- PredictProR:::gp_py_ml_fit_predict_batch(list(eager, lazy))
  expect_length(out, 2L)
  expect_identical(written[[1L]], written[[2L]])
  # missing labels were filtered from the lazily built training rows as well
  expect_identical(length(written[[2L]]$y_train_csv) - 1L, sum(!is.na(y[-tst])))
})

test_that("ML bootstrap drops missing labels before sampling", {
  captured <- NULL

  local_mocked_bindings(
    gp_ml_cli_backend = function() "cli",
    gp_ml_bootstrap_cli = function(model_type,
                                   X_train,
                                   y_train,
                                   X_pred,
                                   response_family = "gaussian",
                                   model_params = list(),
                                   class_levels = NULL,
                                   n_bootstrap = 100L,
                                   seed = NULL) {
      captured <<- list(X_train = X_train, y_train = y_train, X_pred = X_pred)
      matrix(0, nrow = as.integer(n_bootstrap), ncol = nrow(X_pred))
    },
    .package = "PredictProR"
  )

  X_train <- matrix(seq_len(12), nrow = 4)
  y_train <- c(1, NA, 3, 4)
  X_pred <- matrix(seq_len(8), nrow = 2)

  out <- PredictProR:::gp_py_ml_bootstrap_matrix(
    model_type = "svm",
    X_train = X_train,
    y_train = y_train,
    X_pred = X_pred,
    response_family = "gaussian",
    n_bootstrap = 3L
  )

  expect_equal(dim(out), c(3L, nrow(X_pred)))
  expect_false(anyNA(captured$y_train))
  expect_identical(as.numeric(captured$y_train), c(1, 3, 4))
  expect_equal(captured$X_train, X_train[c(1, 3, 4), , drop = FALSE])
})

test_that("DL fit-predict drops missing training labels before Python bridge", {
  captured <- NULL

  local_mocked_bindings(
    gp_dl_run_cli = function(command, args, python_bin = NULL,
                             max_attempts = 1L, run_fn = system2) {
      arg_value <- function(flag) args[[match(flag, args) + 1L]]
      x_train <- utils::read.csv(arg_value("--x-train-csv"), check.names = FALSE)
      y_train <- utils::read.csv(arg_value("--y-train-csv"), check.names = FALSE)
      out_dir <- arg_value("--out-dir")
      captured <<- list(X_train = as.matrix(x_train), y_train = y_train[[1L]])
      utils::write.csv(data.frame(prediction = rep(0, 2)), file.path(out_dir, "predictions.csv"), row.names = FALSE)
      invisible(character())
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::gp_dl_bridge_fit_predict_raw(
    model_type = "mlp",
    X_train = matrix(seq_len(12), nrow = 4),
    y_train = c(1, NA, 3, 4),
    X_test = matrix(seq_len(6), nrow = 2),
    response_family = "gaussian",
    dl_args = list()
  )

  expect_length(res$predictions, 2L)
  expect_false(anyNA(captured$y_train))
  expect_identical(as.numeric(captured$y_train), c(1, 3, 4))
  expect_equal(nrow(captured$X_train), 3L)
})

test_that("multitask DL missingness keeps partially observed trait rows", {
  y_train <- matrix(
    c(
      1, NA,
      NA, 2,
      NaN, NA,
      3, 4
    ),
    nrow = 4L,
    byrow = TRUE
  )

  keep <- PredictProR:::gp_py_ml_nonmissing_training_rows(
    y_train,
    response_family = "multitask_regression"
  )

  expect_identical(keep, c(TRUE, TRUE, FALSE, TRUE))
})
