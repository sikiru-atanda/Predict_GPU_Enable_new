test_that("BGLR ordinal fits keep the declared class order, not alphabetical order", {
  skip_if_not_installed("BGLR")
  # BGLR::BGLR() runs as.vector(y) and then factor(y, ordered = TRUE), so a
  # response passed as labels is ordered ALPHABETICALLY whatever its factor
  # levels say. "low/medium/high" became high < low < medium: the ordinal
  # thresholds were fitted on a scrambled scale and predictions collapsed.
  set.seed(11)
  n <- 240L
  x <- stats::rnorm(n)
  liability <- 2 * x + stats::rnorm(n)
  lv <- c("low", "medium", "high")
  y <- ordered(lv[findInterval(liability, stats::quantile(liability, c(1/3, 2/3))) + 1L], levels = lv)
  y_fit <- y
  y_fit[1:40] <- NA

  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  fit <- PredictProR:::gp_bayes_bglr_fit(
    fam = "ordinal",
    y = PredictProR:::gp_bayes_prepare_response(y_fit, "ordinal"),
    response_type = "ordinal",
    ETA = list(list(X = cbind(x = x), model = "FIXED")),
    nIter = 1500L, burnIn = 500L, thin = 5L, verbose = FALSE,
    saveAt = "ordinal_order_"
  )

  expect_identical(as.character(fit$levels), lv)
  expect_identical(colnames(fit$probs), lv)
  expect_identical(fit$y[!is.na(y_fit)], as.character(y_fit[!is.na(y_fit)]))
  # probability-weighted class index must rise with the liability driver
  expected_class <- as.vector(fit$probs %*% seq_along(lv))
  expect_gt(stats::cor(expected_class[1:40], x[1:40]), 0.8)
})

test_that("ordinal codes sort correctly for ten or more classes", {
  lv <- sprintf("c%d", 1:12)
  y <- ordered(lv[c(1, 2, 10, 12, 9)], levels = lv)
  codes <- PredictProR:::gp_bayes_prepare_response(y, "ordinal")
  expect_identical(attr(codes, "class_levels"), lv)
  # zero-padded codes sort as strings in the declared order ("10" > "02")
  expect_identical(order(as.vector(codes)), order(as.integer(y)))
})
