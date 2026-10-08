test_that("gp_grid_tune_cv can score all fold/candidate jobs through a batch callback", {
  y <- seq_len(18)
  seen <- 0L
  res <- gp_grid_tune_cv(
    y = y,
    response_family = "gaussian",
    param_grid = list(bias = c(0, 3)),
    nfolds = 3L,
    random_state = 99L,
    predict_fun = function(tst, params) {
      stop("serial predict_fun should not be used when batch callback succeeds")
    },
    batch_predict_fun = function(specs) {
      seen <<- length(specs)
      lapply(specs, function(spec) {
        list(pred = y[spec$tst] + as.numeric(spec$params$bias), prob = NULL)
      })
    }
  )

  expect_equal(seen, 6L)
  expect_equal(res$best_params$bias, 0)
  expect_equal(nrow(res$all_results), 2L)
})

test_that("gp_grid_tune_cv falls back to serial scoring if the batch callback fails", {
  y <- seq_len(12)
  serial_calls <- 0L
  res <- NULL
  expect_warning(
    res <- gp_grid_tune_cv(
      y = y,
      response_family = "gaussian",
      param_grid = list(bias = c(0, 2)),
      nfolds = 3L,
      random_state = 101L,
      predict_fun = function(tst, params) {
        serial_calls <<- serial_calls + 1L
        list(pred = y[tst] + as.numeric(params$bias), prob = NULL)
      },
      batch_predict_fun = function(specs) {
        stop("simulated batch failure")
      }
    ),
    regexp = "falling back to serial tuning"
  )

  expect_equal(serial_calls, 6L)
  expect_equal(res$best_params$bias, 0)
})
