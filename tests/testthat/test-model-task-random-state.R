test_that("model task seeds honor the public random_state", {
  expect_identical(PredictProR:::gp_model_task_seed(123L, 1L), 10123L)
  expect_identical(PredictProR:::gp_model_task_seed(456L, 1L), 10456L)
  expect_identical(PredictProR:::gp_model_task_seed(123L, 2L), 20123L)
  expect_false(identical(
    PredictProR:::gp_model_task_seed(123L, 1L),
    PredictProR:::gp_model_task_seed(124L, 1L)
  ))
})

test_that("model task seed validation is bounded and deterministic", {
  expect_identical(
    PredictProR:::gp_model_task_seed(.Machine$integer.max - 10, 1L),
    9990L
  )
  expect_error(PredictProR:::gp_model_task_seed(123L, 0L), "positive integer")
  expect_identical(PredictProR:::gp_model_task_seed(NA_real_, 1L), 10123L)
})
