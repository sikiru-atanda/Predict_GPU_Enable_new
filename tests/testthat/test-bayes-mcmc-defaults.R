# TDD (audit finding): the default MCMC settings (nIter=200, burnIn=50, thin=1)
# produce unconverged posteriors. The package's own warnings flag anything below
# nIter=16000 / burnIn=1600 as sub-optimal, so the defaults should at least meet
# those thresholds (no self-warning) and yield converged runs out of the box.

.resolve_bayes_param_check <- function() {
  if (exists("bayes_parameter_check", mode = "function")) return(get("bayes_parameter_check"))
  if (requireNamespace("PredictProR", quietly = TRUE)) {
    return(tryCatch(getFromNamespace("bayes_parameter_check", "PredictProR"), error = function(e) NULL))
  }
  NULL
}

test_that("default Bayesian MCMC settings meet the package's optimal thresholds", {
  bayes_parameter_check <- .resolve_bayes_param_check()
  if (is.null(bayes_parameter_check)) skip("bayes_parameter_check not available")
  para <- bayes_parameter_check(message = FALSE)
  expect_gte(para$nIter, 16000L)
  expect_gte(para$burnIn, 1600L)
  expect_gte(para$thin, 1L)
  expect_lt(para$burnIn, para$nIter)
})

test_that("explicit Bayesian MCMC settings are still honored and validated", {
  bayes_parameter_check <- .resolve_bayes_param_check()
  if (is.null(bayes_parameter_check)) skip("bayes_parameter_check not available")
  para <- bayes_parameter_check(nIter = 500L, burnIn = 100L, thin = 2L, message = FALSE)
  expect_equal(para$nIter, 500L)
  expect_equal(para$burnIn, 100L)
  expect_equal(para$thin, 2L)
  # burnIn >= nIter must still error.
  expect_error(bayes_parameter_check(nIter = 100L, burnIn = 100L, message = FALSE))
})
