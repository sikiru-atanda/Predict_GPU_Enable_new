# TDD (audit follow-up): the ML bootstrap branch of reliability_thresholds_MPIW_from_CI
# hardcoded probs = 0.05/0.95 (a 90% bootstrap interval), ignoring confidence_level
# -- the same defect just fixed for the Bayesian branch. It must honor
# confidence_level so ML/DL bootstrap intervals match the rest of the package.

.resolve_fn_ml <- function(name) {
  if (exists(name, mode = "function")) return(get(name))
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(tryCatch(getFromNamespace(name, "PredictProR"), error = function(e) NULL))
  }
  NULL
}

test_that("ML bootstrap CI honors confidence_level", {
  fn <- .resolve_fn_ml("reliability_thresholds_MPIW_from_CI")
  if (is.null(fn)) skip("reliability_thresholds_MPIW_from_CI not available")
  set.seed(2)
  # boot_results$t: rows = bootstrap replicates, cols = observations.
  # 6 observations, 4000 reps ~ N(mean_j, 1); interval WIDTH depends on the
  # per-column spread (sd = 1), not the mean.
  means <- c(-2, -1, 0, 1, 2, 3)
  t_mat <- sapply(means, function(m) stats::rnorm(4000, mean = m, sd = 1))  # 4000 x 6
  boot <- list(t = t_mat)

  width <- function(cl) {
    r <- fn(boot_results = boot, model_for_CI_cal = "ML",
            confidence_level = cl, CI_width_thresholds = c(0.33, 0.66))
    mean(r$upper_bound - r$lower_bound)
  }
  w95 <- width(0.95)
  w80 <- width(0.80)

  expect_gt(w95, w80)
  expect_equal(w95, 2 * stats::qnorm(0.975), tolerance = 0.12)  # ~3.92
  expect_equal(w80, 2 * stats::qnorm(0.90), tolerance = 0.12)   # ~2.56
})

test_that("ML bootstrap summaries tolerate non-finite replicate predictions", {
  fn <- .resolve_fn_ml("reliability_thresholds_MPIW_from_CI")
  if (is.null(fn)) skip("reliability_thresholds_MPIW_from_CI not available")
  boot <- list(t = cbind(
    c(1, 1.2, NaN, 0.8),
    c(2, NA_real_, 2.2, 1.8),
    rep(NaN, 4L)
  ))

  out <- fn(
    boot_results = boot,
    model_for_CI_cal = "ML",
    confidence_level = 0.95,
    CI_width_thresholds = c(0.33, 0.66)
  )

  expect_true(all(is.finite(out$lower_bound[1:2])))
  expect_true(all(is.finite(out$upper_bound[1:2])))
  expect_true(all(is.finite(out$standard_errors[1:2])))
  expect_true(is.na(out$lower_bound[3]))
  expect_true(is.na(out$standard_errors[3]))
})
