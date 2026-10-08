# Guard: direct GP models (KRR, GP, LowRankGP) can combine several genomic or
# omics kernels but retain a single-environment covariance structure. On a
# multi-environment panel they can collapse the genetic signal into the residual
# (per-env h2 ~0.08-0.16 vs ~0.5-0.7 for GP_FA / ASReml corgh). Warn and steer
# users to GP_FA / ASReml corgh.
resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  ns <- tryCatch(asNamespace("PredictProR"), error = function(e) NULL)
  if (!is.null(ns) && exists(name, envir = ns, mode = "function")) return(get(name, envir = ns, mode = "function"))
  NULL
}

test_that("gp_warn_single_kernel_gp_met warns for direct GP + MET, silent otherwise", {
  fn <- resolve_fn("gp_warn_single_kernel_gp_met")
  skip_if(is.null(fn), "gp_warn_single_kernel_gp_met not available")

  # MET + direct GP -> warn (and name GP_FA in the message)
  expect_warning(fn("GP", NULL, FALSE, is_met = TRUE), "GP_FA")
  expect_warning(fn("KRR", NULL, FALSE, is_met = TRUE), "multiple genomic/omics kernels")
  expect_warning(fn("LowRankGP", NULL, FALSE, is_met = TRUE), "multi-environment")

  # MET + GP_FA (the correct MET model) -> silent
  expect_silent(fn("GP_FA", NULL, FALSE, is_met = TRUE))
  # MET + non-GP model -> silent
  expect_silent(fn("GBLUP", NULL, FALSE, is_met = TRUE))
  # single-environment -> silent regardless of model
  expect_silent(fn("GP", NULL, FALSE, is_met = FALSE))
  expect_silent(fn("KRR", NULL, FALSE, is_met = FALSE))
})

test_that("gp_warn_single_kernel_gp_met uses GS_model_cv when cross_validation = TRUE", {
  fn <- resolve_fn("gp_warn_single_kernel_gp_met")
  skip_if(is.null(fn), "gp_warn_single_kernel_gp_met not available")
  # CV path: the CV model name is what matters
  expect_warning(fn(GS_model = "GP_FA", GS_model_cv = "GP", cross_validation = TRUE, is_met = TRUE),
                 "single-environment covariance structure")
  expect_silent(fn(GS_model = "GP", GS_model_cv = "GP_FA", cross_validation = TRUE, is_met = TRUE))
})
