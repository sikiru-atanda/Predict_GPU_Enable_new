test_that("ML predictor preprocessor preserves columns and only formats finite matrices", {
  x <- matrix(
    c(
      1, NA, 3,
      4, 5, NA,
      7, 8, NA
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:3), c("m1", "m2", "m3"))
  )

  prep <- PredictProR:::gp_ml_preprocess_predictors(x, scaling = FALSE, centering = FALSE)

  expect_equal(dim(prep$data), dim(x))
  expect_identical(colnames(prep$data), colnames(x))
  expect_identical(rownames(prep$data), rownames(x))
  expect_length(prep$removed_cols, 0L)
  expect_false(anyNA(prep$data))
  expect_identical(as.numeric(prep$data[, "m3"]), c(3, 3, 3))

  test_x <- x[c(2, 1), c("m3", "m1", "m2"), drop = FALSE]
  applied <- PredictProR:::gp_ml_apply_preprocessor(test_x, prep)

  expect_equal(dim(applied), c(2L, 3L))
  expect_identical(colnames(applied), colnames(x))
  expect_false(anyNA(applied))
})

test_that("DL CV payload keeps predictor columns with residual missing values", {
  x <- matrix(
    c(
      1, NA, 3,
      4, 5, NA,
      7, 8, 9,
      10, NA, 12
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("g", 1:4), c("m1", "m2", "m3"))
  )
  y <- c(1, 2, 3, 4)

  payload <- PredictProR:::gp_dl_get_cached_cv_payload(
    omics_data = x,
    y_raw = y,
    tst = c(2L, 4L),
    fam = "gaussian",
    scaling = FALSE,
    centering = FALSE
  )

  expect_equal(dim(payload$omics_data), dim(x))
  expect_equal(ncol(payload$X_tr), ncol(x))
  expect_equal(ncol(payload$X_te), ncol(x))
  expect_false(anyNA(payload$omics_data))
  expect_false(anyNA(payload$X_tr))
  expect_false(anyNA(payload$X_te))
})

test_that("DL CV payload fits preprocessing on fold training rows and transforms all uncertainty targets", {
  x <- cbind(m1 = c(1, 100, 3, 200), m2 = c(2, 50, 4, 80))
  target <- rbind(x, c(m1 = 500, m2 = 120))
  payload <- PredictProR:::gp_dl_get_cached_cv_payload(
    omics_data = x,
    y_raw = c(1, 2, 3, 4),
    tst = c(2L, 4L),
    fam = "gaussian",
    scaling = TRUE,
    centering = TRUE,
    target_data = target
  )

  expect_equal(unname(colMeans(payload$X_tr)), c(0, 0), tolerance = 1e-12)
  expect_equal(dim(payload$X_target), dim(target))
  expect_false(anyNA(payload$X_target))
  expect_gt(payload$X_te[1, 1], 50)
})
