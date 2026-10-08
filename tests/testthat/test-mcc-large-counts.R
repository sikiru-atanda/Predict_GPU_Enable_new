test_that("MCC stays finite for large confusion-matrix counts", {
  # Confusion-matrix counts come from table() and are integers; the MCC
  # denominator multiplies four of them and overflowed R's integer range
  # (~2.1e9) once each margin passed roughly 215, returning NA. Pooled MET
  # cross-validation hit this and then failed with "no rows to aggregate".
  lv <- c("no", "yes")
  y <- factor(rep(c("no", "yes", "no", "yes"), c(300, 100, 100, 300)), levels = lv)
  y_hat <- factor(rep(c("no", "no", "yes", "yes"), c(300, 100, 100, 300)), levels = lv)
  # tp = tn = 300, fp = fn = 100  ->  MCC = (300^2 - 100^2) / 400^2 = 0.5
  mcc <- PredictProR::evaluation_metrics(
    y_observed = y, y_predicted = y_hat, eval_metrics = "mcc",
    response_family = "binary", positive_class = "yes"
  )
  expect_equal(as.numeric(mcc), 0.5)
})
