# 0.20.221: ASReml reliability is genetic, 1 - genetic PEV / (sigma2_g * K_ii), the quantity the
# GP engine reports. It was 1 - total-mean SE^2 / var(y), about twice the genetic value (0.77 vs
# 0.41 on a 150-line panel). PEV and Standard_error keep the prediction's total SE^2; the genetic
# PEV is Reliability_variance_input. The single-environment variance table also carries ASReml's
# standard errors (they were NA).

asreml_ready <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) return(FALSE)
  tryCatch(grepl("No error", asreml::asreml.license.status(quiet = TRUE)$rv_text),
           error = function(e) FALSE)
}

test_that("without genetic inputs, ASReml reliability is NA, never 1 - SE^2 / var(y)", {
  ph <- data.frame(GID = c("a", "b", "c"), y = c(1, 4, NA))
  pred <- data.frame(GID = c("a", "b", "c"), Predicted_value = c(1, 3, 2),
                     Standard_error = c(0.5, 0.6, 0.9))
  out <- asreml_add_uncertainty_fields(pred, ph, "y")
  expect_equal(out$PEV, pred$Standard_error^2)
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(is.na(out$Reliability_reference_variance)))
  expect_match(out$Reliability_basis[1], "^unavailable")
  pred$Reliability_variance_input <- c(0.2, 0.3, 0.8)
  pred$Reliability_reference_variance <- 1
  out <- asreml_add_uncertainty_fields(pred, ph, "y")
  expect_equal(out$Reliability, c(0.8, 0.7, 0.2))
})

test_that("environment genetic variances become a kernels x environments matrix", {
  envs <- c("E1", "E2")
  m <- asreml_met_varG_matrix(list(varG_per_omics = matrix(1:4, 2)), 2L, envs)
  expect_equal(dim(m), c(2L, 2L)); expect_equal(colnames(m), envs)
  m <- asreml_met_varG_matrix(list(varG_per_omics = list(c(1, 2), c(3, 4))), 2L, envs)
  expect_equal(m[2, ], c(E1 = 3, E2 = 4))
  m <- asreml_met_varG_matrix(list(Total_genetic_var = c(5, 6)), 1L, envs)
  expect_equal(m[1, ], c(E1 = 5, E2 = 6))
  expect_null(asreml_met_varG_matrix(list(Total_genetic_var = 5), 1L, envs))
})

test_that("compute_inverse_and_sparse keeps the fitted kernel's diagonal", {
  K <- matrix(c(1.2, 0.3, 0.3, 0.9), 2, dimnames = list(c("a", "b"), c("a", "b")))
  inv <- compute_inverse_and_sparse(K, epsilon = 1e-6, inverse = TRUE)
  expect_equal(attr(inv, "kernel_diag"), c(a = 1.2 + 1e-6, b = 0.9 + 1e-6))
})

test_that("ASReml genetic reliability equals a Henderson mixed-model-equations reference", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  withr::local_seed(11)
  n <- 80; p <- 300
  M <- matrix(stats::rbinom(n * p, 2, 0.3), n, p)
  Z <- scale(M); K <- tcrossprod(Z) / p
  dimnames(K) <- list(sprintf("L%02d", 1:n), sprintf("L%02d", 1:n))
  g <- drop(Z[, 1:20] %*% stats::rnorm(20)); g <- g / stats::sd(g)
  ph <- data.frame(GID = factor(rownames(K), levels = rownames(K)), y = 10 + g + stats::rnorm(n))
  ph$y[71:80] <- NA
  Kinv <- compute_inverse_and_sparse(kernel = K, epsilon = 1e-6, inverse = TRUE)
  assign("rel_test_Kinv", Kinv, envir = .GlobalEnv)
  withr::defer(rm("rel_test_Kinv", envir = .GlobalEnv))
  mod <- asreml::asreml(fixed = y ~ 1, random = ~ vm(GID, rel_test_Kinv), data = ph, trace = FALSE)
  # predict() works at the last iteration's parameters: update until the components settle,
  # so the reference (built from the reported components) describes the same fit
  for (i in 1:6) mod <- suppressWarnings(asreml::update.asreml(mod))
  vc <- summary(mod)$varcomp
  s2g <- vc[grep("^vm", rownames(vc)), "component"]
  s2e <- vc[grep("!R$", rownames(vc)), "component"]
  pred <- asreml_run_with_model_data_context(mod, ph, asreml_predict_pvals(mod, classify = "GID"))
  names(pred)[names(pred) == "std.error"] <- "Standard_error"
  out <- asreml_apply_genetic_reliability(
    pred_df = pred, mod = mod, pheno_data = ph, gen_name = "GID", classify = "GID",
    varG_per_omics = matrix(s2g), names_in_inv_list = "rel_test_Kinv"
  )
  # Henderson: C22 of [X'X/s2e, X'Z/s2e; Z'X/s2e, Z'Z/s2e + K^-1/s2g]^-1 is Var(ghat - g)
  obs <- !is.na(ph$y)
  X <- matrix(1, sum(obs), 1); Zm <- diag(n)[obs, , drop = FALSE]
  Kf <- K + diag(1e-6, n)
  C <- rbind(cbind(crossprod(X) / s2e, crossprod(X, Zm) / s2e),
             cbind(crossprod(Zm, X) / s2e, crossprod(Zm) / s2e + solve(Kf) / s2g))
  pev_ref <- diag(solve(C))[-1]
  rel_ref <- 1 - pev_ref / (s2g * diag(Kf))
  idx <- match(rownames(K), as.character(out$GID))
  expect_equal(out$Reliability_variance_input[idx], unname(pev_ref), tolerance = 1e-4)
  expect_equal(out$Reliability_reference_variance[idx], unname(s2g * diag(Kf)), tolerance = 1e-8)
  expect_equal(out$Reliability[idx], unname(rel_ref), tolerance = 1e-4)
  expect_true(all(out$Reliability_variance_input < out$Standard_error^2))
  expect_match(out$Reliability_basis[1], "^genetic_prior_variance")
})

asreml_e2e_panel <- function(n_env = 1L) {
  withr::local_seed(3, .local_envir = parent.frame())
  n <- 90; p <- 300
  M <- matrix(stats::rbinom(n * p, 2, 0.3), n, p,
              dimnames = list(sprintf("L%02d", 1:n), sprintf("m%03d", 1:p)))
  Z <- scale(M); G <- tcrossprod(Z) / p; diag(G) <- diag(G) + 1e-3
  g <- drop(Z[, 1:20] %*% stats::rnorm(20)); g <- g / stats::sd(g)
  envs <- paste0("E", seq_len(n_env))
  ph <- do.call(rbind, lapply(seq_along(envs), function(e) {
    data.frame(GID = rownames(M), Env = envs[e], Yield = 5 * e + g + stats::rnorm(n))
  }))
  ph$Yield[ph$GID %in% rownames(M)[81:90]] <- NA
  list(G = G, ph = ph)
}

test_that("single-environment ASReml output: genetic reliability and variance-component SEs", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- asreml_e2e_panel(1L)
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(model_execute(
    pheno_data = d$ph[, c("GID", "Yield")], gmatrix = d$G, gen_name = "GID", response = "Yield",
    response_family = "gaussian", fixed = ~1, random = ~GID, GS_model = "GBLUP", engine = "asreml",
    cross_validation = FALSE, parallel_mode = "sequential", system_database = FALSE,
    message = FALSE, verbose = FALSE
  ))
  mr <- res$model_results_by_trait$Yield
  if (!"predicted_values" %in% names(mr)) mr <- mr[[1L]]
  vc <- mr$variance_components
  se <- stats::setNames(vc$Standard_error, rownames(vc))
  expect_true(all(is.finite(se[c("genetic_variance", "residual_variance", "heritability")])))
  expect_true(all(se[c("genetic_variance", "residual_variance", "heritability")] > 0))
  pv <- mr$predicted_values
  expect_equal(pv$PEV, pv$Standard_error^2, tolerance = 1e-8)
  expect_true(all(is.finite(pv$Reliability)))
  expect_equal(pv$Reliability,
               pmax(0, pmin(1, 1 - pv$Reliability_variance_input / pv$Reliability_reference_variance)),
               tolerance = 1e-8)
  expect_true(all(pv$Reliability_variance_input < pv$PEV))
  expect_true(all(grepl("^genetic_prior_variance", pv$Reliability_basis)))
  # test lines have no phenotype, so their genetic reliability is lower than the training lines'
  expect_lt(mean(pv$Reliability[pv$Train_Test_Label == "Test"]),
            mean(pv$Reliability[pv$Train_Test_Label == "Train"]))
  # the contract holds the identity: a reliability computed from anything else is rejected
  expect_true(all(validate_prediction_output(mr, gen_name = "GID", strict = FALSE)$valid))
  bad <- mr
  bad$predicted_values$Reliability[1] <- 1 - bad$predicted_values$PEV[1] /
    stats::var(d$ph$Yield, na.rm = TRUE)
  report <- validate_prediction_output(bad, gen_name = "GID", strict = FALSE)
  expect_false(all(report$valid))
  expect_match(paste(report$issues, collapse = " "), "genetic Reliability must equal")
})

test_that("multi-environment ASReml output: genotype-by-environment genetic reliability", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- asreml_e2e_panel(2L)
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(model_execute(
    pheno_data = d$ph, gmatrix = d$G, gen_name = "GID", response = "Yield",
    response_family = "gaussian", fixed = ~Env, random = ~GID + GID:Env, heter_groups = "Env",
    GS_model = "GBLUP", engine = "asreml", cross_validation = FALSE,
    parallel_mode = "sequential", system_database = FALSE, message = FALSE, verbose = FALSE
  ))
  mr <- res$model_results_by_trait$Yield
  if (!"predicted_values" %in% names(mr)) mr <- mr[[1L]]
  ge <- if ("Env" %in% names(mr$predicted_values)) mr$predicted_values else mr$across_environment_predicted_values
  expect_true("Env" %in% names(ge))
  ok <- is.finite(ge$Reliability)
  expect_true(any(ok))
  expect_true(all(grepl("^genetic_prior_variance", ge$Reliability_basis[ok])))
  expect_equal(ge$Reliability[ok],
               pmax(0, pmin(1, 1 - ge$Reliability_variance_input[ok] / ge$Reliability_reference_variance[ok])),
               tolerance = 1e-8)
  expect_true(all(ge$Reliability_variance_input[ok] < ge$PEV[ok]))
})

test_that("the genetic-PEV prediction retries once with more prediction workspace", {
  calls <- list()
  fake_predict <- function(mod, classify, only, sed, pworkspace = NULL, ...) {
    calls[[length(calls) + 1L]] <<- pworkspace
    if (length(calls) == 1L) stop("Insufficient workspace available when calculating predictions.")
    list(pvals = data.frame(GID = c("a", "b"), predicted.value = c(0.1, -0.1), std.error = c(0.5, 0.6)))
  }
  withr::local_options(PredictProR.asreml_predict_impl = fake_predict)
  mod <- list(call = quote(asreml(data = ph)))
  pv <- asreml_genetic_pev_table(mod, data.frame(GID = c("a", "b")), "GID", "vm(GID, K)",
                                 pworkspace = 1e6)
  expect_equal(pv$std.error, c(0.5, 0.6))
  expect_equal(calls[[1L]], 1e6)
  expect_identical(calls[[2L]], "1gb")
  # any other error is not retried
  calls <- list()
  withr::local_options(PredictProR.asreml_predict_impl = function(...) stop("singular"))
  expect_error(asreml_genetic_pev_table(mod, data.frame(GID = "a"), "GID", "vm(GID, K)"), "singular")
})

test_that("the total-mean ASReml prediction retries once with more prediction workspace", {
  calls <- list()
  fake_predict <- function(mod, classify, sed, pworkspace = NULL, ...) {
    calls[[length(calls) + 1L]] <<- pworkspace
    if (length(calls) == 1L) stop("Insufficient workspace available when calculating predictions.")
    list(pvals = data.frame(GID = c("a", "b"), predicted.value = c(1, 2), std.error = c(0.5, 0.6)))
  }
  withr::local_options(PredictProR.asreml_predict_impl = fake_predict)
  pv <- asreml_predict_pvals(list(), classify = "GID", pworkspace = 1e6)
  expect_equal(pv$std.error, c(0.5, 0.6))
  expect_equal(calls[[1L]], 1e6)
  expect_identical(calls[[2L]], "1gb")
  withr::local_options(PredictProR.asreml_predict_impl = function(...) stop("singular"))
  expect_error(asreml_predict_pvals(list(), classify = "GID"), "singular")
})
