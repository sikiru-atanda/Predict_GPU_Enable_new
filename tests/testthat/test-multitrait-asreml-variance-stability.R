# 0.20.223: a multi-trait ASReml fit that converges but whose variance components still move by
# more than 1% per continuation update used to stop with an error, so the user got no
# predictions (grouped MT-MET at 4 traits x 3 environments in the trait-count capacity study).
# It now returns the predictions and withholds only the variance output, with a warning.

asreml_ready <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) return(FALSE)
  tryCatch(grepl("No error", asreml::asreml.license.status(quiet = TRUE)$rv_text),
           error = function(e) FALSE)
}

stability_panel <- function(traits = c("T1", "T2"), envs = NULL, n_fam = 8, fam_size = 5, p = 120) {
  withr::local_seed(17)
  n <- n_fam * fam_size
  fam <- rep(seq_len(n_fam), each = fam_size)
  P <- matrix(stats::rbinom(2 * n_fam * p, 1, 0.3), 2 * n_fam, p)
  M <- t(vapply(seq_len(n), function(i) {
    f <- fam[i]
    P[2 * f - 1, ] + ifelse(stats::runif(p) < 0.5, P[2 * f, ], P[2 * f - 1, ])
  }, numeric(p)))
  Z <- scale(M); Z[is.na(Z)] <- 0
  K <- tcrossprod(Z) / p
  ids <- sprintf("L%03d", seq_len(n)); dimnames(K) <- list(ids, ids)
  g <- drop(Z[, 1:20] %*% stats::rnorm(20))
  one <- function(shift) {
    y <- sapply(seq_along(traits), function(j) scale(g)[, 1] * 0.8 + stats::rnorm(n, sd = 0.7) + shift + j)
    d <- data.frame(GID = ids, y, stringsAsFactors = FALSE); names(d)[-1] <- traits; d
  }
  ph <- if (is.null(envs)) one(0) else do.call(rbind, lapply(seq_along(envs), function(e) cbind(Env = envs[e], one(e))))
  ph[ph$GID %in% tail(ids, 5), traits] <- NA
  list(K = K, ph = ph, traits = traits, ids = ids)
}

# the real stability helper, but reporting a converged fit whose components still moved by 2%
mock_unstable <- function() {
  real <- gp_multitrait_asreml_fit_until_stable
  function(model, update_fn, max_updates = 4L, max_percent_change = 1) {
    out <- real(model, update_fn, max_updates = 1L, max_percent_change = Inf)
    out$trace$maximum_percent_change[nrow(out$trace)] <- 2
    out$trace$converged[nrow(out$trace)] <- TRUE
    out$stable <- FALSE
    out
  }
}

test_that("grouped MT-MET ASReml returns predictions when its variance is not stable", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- stability_panel(envs = c("E1", "E2", "E3"))
  local_mocked_bindings(gp_multitrait_asreml_fit_until_stable = mock_unstable())
  expect_warning(
    res <- gp_multitrait_asreml_mtmet_gaussian_model(
      d$ph, d$traits, "GID", heter_groups = "Env", gmatrix = d$K, var_cov_str = "fa1"
    ),
    "converged but did not reach a variance-stable endpoint.*Predictions are returned"
  )
  expect_gt(nrow(res$predicted_values), 0L)
  expect_true(all(is.finite(res$predicted_values$Predicted_value)))
  expect_identical(res$variance_component_status, "not_promoted_nonconverged")
  expect_equal(nrow(res$variance_components), 0L)
  mp <- res$model_parameters
  expect_identical(mp$summary[mp$stat == "variance_endpoint"], "converged_variance_unstable_not_promoted")
})

test_that("multi-kernel single-environment ASReml withholds variance and reliability, not predictions", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- stability_panel(traits = c("T1", "T2", "T3"))
  K2 <- 0.6 * d$K + 0.4 * diag(nrow(d$K)); dimnames(K2) <- dimnames(d$K)
  local_mocked_bindings(gp_multitrait_asreml_fit_until_stable = mock_unstable())
  expect_warning(
    res <- gp_multitrait_asreml_gaussian_model(
      d$ph, d$traits, "GID", gmatrix = d$K, kernel_list = list(second = K2), trait_covariance = "us"
    ),
    "converged but did not reach a variance-stable endpoint"
  )
  pl <- res$predicted_values
  expect_true(all(is.finite(pl$Predicted_value)))
  expect_true(all(is.na(pl$Reliability)))
  expect_match(pl$Reliability_basis[1], "^unavailable: variance components did not reach a stable endpoint")
  expect_identical(res$variance_component_status, "not_promoted_nonconverged")
  expect_equal(nrow(res$variance_components), 0L)
  expect_null(res$Genetic_covariance_traits)
})

test_that("a stable fit is unchanged: variance promoted", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- stability_panel(envs = c("E1", "E2", "E3"))
  res <- suppressWarnings(gp_multitrait_asreml_mtmet_gaussian_model(
    d$ph, d$traits, "GID", heter_groups = "Env", gmatrix = d$K, var_cov_str = "fa1",
    require_variance_stability = FALSE
  ))
  expect_identical(res$variance_component_status, "promoted")
  expect_gt(nrow(res$variance_components), 0L)
})

test_that("model_execute reports the withheld variance as not promoted for grouped MT-MET", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- stability_panel(envs = c("E1", "E2", "E3"))
  local_mocked_bindings(gp_multitrait_asreml_fit_until_stable = mock_unstable())
  withr::local_dir(withr::local_tempdir())
  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = d$ph, gmatrix = d$K, response = d$traits, gen_name = "GID", GS_model = "GBLUP",
    engine = "asreml", response_family = "gaussian", multi_trait_asreml = TRUE,
    heter_groups = "Env", heter_resid = TRUE, var_cov_str = "fa1", message = FALSE
  ))
  preds <- out$model_results$predicted_values
  expect_gt(nrow(preds), 0L)
  expect_true(all(is.finite(preds$Predicted_value)))
  vc <- out$model_results$variance_components
  expect_true(any(grepl("not_promoted", vc$Component)))
  expect_true(all(is.na(vc$Components)))
})
