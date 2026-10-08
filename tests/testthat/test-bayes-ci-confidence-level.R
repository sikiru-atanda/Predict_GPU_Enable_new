# TDD (audit finding): the Bayesian marker posterior-quantile credible interval
# hardcoded probs = 0.05/0.95 (a 90% interval), ignoring `confidence_level`,
# while kernel/multi-trait paths used yHat +/- z(confidence_level)*SE (95%). The
# same lower_bound/upper_bound columns therefore carried different nominal
# coverage by model type. The Bayes branch must honor confidence_level, taking
# posterior quantiles at (1-cl)/2 and 1-(1-cl)/2.

.resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name))
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(tryCatch(getFromNamespace(name, "PredictProR"), error = function(e) NULL))
  }
  NULL
}

test_that("Bayesian posterior-quantile CI honors confidence_level", {
  reliability_thresholds_MPIW_from_CI <- .resolve_fn("reliability_thresholds_MPIW_from_CI")
  if (is.null(reliability_thresholds_MPIW_from_CI)) {
    skip("reliability_thresholds_MPIW_from_CI not available")
  }
  set.seed(1)
  # 6 individuals x 4000 posterior draws ~ N(mean_i, 1); interval WIDTH depends
  # only on the per-row spread (sd = 1), not the mean.
  means <- c(-2, -1, 0, 1, 2, 3)
  draws <- t(sapply(means, function(m) stats::rnorm(4000, mean = m, sd = 1)))
  mod <- list(model = list(mu = 0))

  width <- function(cl) {
    r <- reliability_thresholds_MPIW_from_CI(
      Predicted_value_for_CI = draws, mod = mod,
      model_for_CI_cal = "Bayes", confidence_level = cl,
      CI_width_thresholds = c(0.33, 0.66)
    )
    mean(r$upper_bound - r$lower_bound)
  }
  w95 <- width(0.95)
  w80 <- width(0.80)

  # A 95% credible interval must be strictly wider than an 80% one.
  expect_gt(w95, w80)
  # And each matches the expected normal-quantile spread of the draws.
  expect_equal(w95, 2 * stats::qnorm(0.975), tolerance = 0.12)  # ~3.92
  expect_equal(w80, 2 * stats::qnorm(0.90), tolerance = 0.12)   # ~2.56
})
