test_that("quadratic weighted kappa uses the declared ordinal class order", {
  # gp_quadratic_weighted_kappa() sorted labels alphabetically, so for
  # low < medium < high it weighted distances on high, low, medium: a low/high
  # swap cost the same as a one-class miss. QWK is the default ranking metric
  # for ordinal cross-validation.
  lv <- c("low", "medium", "high")
  y <- ordered(rep(lv, each = 20), levels = lv)
  one_off <- ordered(rep(c("medium", "high", "medium"), each = 20), levels = lv)
  swapped <- ordered(rep(c("high", "medium", "low"), each = 20), levels = lv)
  qwk <- function(pred) as.numeric(PredictProR::evaluation_metrics(
    y_observed = y, y_predicted = pred, eval_metrics = "quadratic_weighted_kappa",
    response_family = "ordinal"
  ))
  expect_equal(qwk(y), 1)
  expect_gt(qwk(one_off), qwk(swapped))
  # swapping the two extreme classes is worse than chance agreement
  expect_lt(qwk(swapped), 0)
})
