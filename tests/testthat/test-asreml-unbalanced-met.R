# Unbalanced MET (lines missing from some environments) used to fail in ASReml
# output processing: kernel coefficients paired EBVs with phenotype rows by
# position ("requires one finite EBV per kernel row").
test_that("ASReml MET runs on unbalanced data and labels coefficients by genotype x environment", {
  skip_if_not_installed("asreml")
  set.seed(31)
  n <- 40L; p <- 120L
  G <- matrix(rbinom(n * p, 2, 0.35), n, dimnames = list(sprintf("L%02d", 1:n), sprintf("m%03d", 1:p)))
  u <- as.numeric(scale(G[, 1:10] %*% rnorm(10)))
  envs <- c("E1", "E2", "E3")
  ph <- do.call(rbind, lapply(seq_along(envs), function(k) {
    keep <- setdiff(seq_len(n), sample(n, 8L))          # each environment misses 8 different lines
    data.frame(GID = rownames(G)[keep], Env = envs[k], Yield = 5 + k + u[keep] + rnorm(length(keep), sd = 0.5))
  }))
  ph$Yield[ph$GID %in% c("L01", "L02", "L03")] <- NA
  withr::local_dir(withr::local_tempdir())
  out <- tryCatch(suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~ Env,
    random = ~ GID + GID:Env, heter_groups = "Env", heter_resid = TRUE, var_cov_str = "corgh",
    GS_model = "GBLUP", engine = "asreml", gmatrix_method = "VanRaden", qc_filtering = FALSE,
    impute = FALSE, ld_prunning_qc = FALSE, parallel_mode = "sequential", system_database = TRUE, message = FALSE
  )), error = function(e) e)
  if (inherits(out, "error") && grepl("License|licen", conditionMessage(out))) skip("ASReml license unavailable")
  expect_false(inherits(out, "error"), info = if (inherits(out, "error")) conditionMessage(out) else "")
  pv <- out$model_results$predicted_values
  expect_true(is.data.frame(pv) && nrow(pv) > 0L)
  expect_true(all(is.finite(pv$Predicted_value)))
  find_coef <- function(x, d = 0) {
    if (is.data.frame(x) && all(c("x_variables", "coeff") %in% names(x))) return(x)
    if (is.list(x) && d < 6) for (e in x) { r <- find_coef(e, d + 1); if (!is.null(r)) return(r) }
    NULL
  }
  cf <- find_coef(out)
  if (!is.null(cf)) {
    expect_true("Env" %in% names(cf))
    expect_false(anyDuplicated(paste(cf$x_variables, cf$Env)) > 0)
  }
})

test_that("invalid-correlation failures are recognised", {
  expect_true(PredictProR:::asreml_is_invalid_correlation_error(
    "Reconstructed correlation matrix is invalid: matrix is not positive semidefinite within tolerance"))
  expect_true(PredictProR:::asreml_is_invalid_correlation_error("A valid correlation matrix with a unit diagonal is required."))
  expect_false(PredictProR:::asreml_is_invalid_correlation_error("Expected 3 correlation parameters for 3 environments, found 2."))
})

# corgh on 60 lines x 3 unbalanced environments converges to genetic
# correlations that are not positive semidefinite. The predictions used to be
# lost ("ASReml output processing failed"); they are kept, with variance
# components reported as unavailable and a warning.
test_that("ASReml MET keeps predictions when the genetic correlation estimate is invalid", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  # Same data as sim_met(n = 60, p = 300) in inst/examples/scenarios/00_setup.R.
  n <- 60L; p <- 300L; envs <- c("E1", "E2", "E3")
  ids <- sprintf("L%03d", seq_len(n))
  set.seed(11)
  freq <- stats::runif(p, 0.1, 0.9)
  G <- matrix(2L * stats::rbinom(n * p, 1, rep(freq, each = n)), n, p, dimnames = list(ids, sprintf("m%04d", seq_len(p))))
  G <- G[, apply(G, 2, stats::var) > 0, drop = FALSE]
  gv <- function(s) { set.seed(s); q <- sample(ncol(G), 30); as.numeric(scale(G[, q, drop = FALSE] %*% stats::rnorm(30))) }
  g_common <- gv(12)
  g_env <- sapply(seq_along(envs), function(k) sqrt(0.6) * g_common + sqrt(0.4) * gv(11 + 100 + k))
  set.seed(13)
  ph <- expand.grid(GID = ids, Env = envs, stringsAsFactors = FALSE)
  env_effect <- stats::setNames(seq(0, by = 2, length.out = length(envs)), envs)
  ph$Yield <- 10 + env_effect[ph$Env] + g_env[cbind(match(ph$GID, ids), match(ph$Env, envs))] +
    stats::rnorm(nrow(ph), sd = 0.7)
  ph <- ph[sort(sample(nrow(ph), round(0.85 * nrow(ph)))), ]
  ph$Weight <- 1
  test_ids <- sample(unique(ph$GID), round(0.2 * n))
  ph$Yield[ph$GID %in% test_ids] <- NA
  rownames(ph) <- NULL
  p_freq <- colMeans(G) / 2
  Z <- sweep(G, 2, 2 * p_freq)
  K <- tcrossprod(Z) / (2 * sum(p_freq * (1 - p_freq))) + diag(1e-4, n)
  withr::local_dir(withr::local_tempdir())
  warns <- character()
  out <- withCallingHandlers(
    tryCatch(suppressMessages(PredictProR::model_execute(
      pheno_data = ph, gmatrix = K, gen_name = "GID", response = "Yield", heter_groups = "Env",
      fixed = ~ Env, random = ~ GID + GID:Env, heter_resid = TRUE, var_cov_str = "corgh",
      GS_model = "GBLUP", engine = "asreml", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
      parallel_mode = "sequential", system_database = TRUE, message = FALSE
    )), error = function(e) e),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(out, "error") && grepl("License|licen", conditionMessage(out))) skip("ASReml license unavailable")
  expect_false(inherits(out, "error"), info = if (inherits(out, "error")) conditionMessage(out) else "")
  pv <- out$model_results$predicted_values
  expect_true(is.data.frame(pv) && nrow(pv) > 0L && all(is.finite(pv$Predicted_value)))
  # When this data set does give a valid correlation estimate (ASReml version
  # differences), there is nothing more to check.
  if (any(grepl("positive semidefinite", warns))) {
    expect_true(all(is.na(pv$Reliability)))
    status <- PredictProR:::gp_vc_find_scalar(out, "variance_component_status")
    expect_match(status, "unavailable")
  }
})

# Large MET skips the kernel coefficients (more genotype x environment rows
# than PREDICTPRO_ASREML_COEF_MAX_ROWS). The skipped slot was dropped from the
# coefficient list, the renaming then failed on an empty list, and the whole
# fitted model was lost (G2F, 23,110 rows, after a 22-minute fit).
test_that("ASReml MET keeps its predictions when kernel coefficients are skipped", {
  skip_if_not_installed("asreml")
  set.seed(32)
  n <- 30L; p <- 100L
  G <- matrix(rbinom(n * p, 2, 0.35), n, dimnames = list(sprintf("L%02d", 1:n), sprintf("m%03d", 1:p)))
  u <- as.numeric(scale(G[, 1:10] %*% rnorm(10)))
  ph <- expand.grid(GID = rownames(G), Env = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
  ph$Yield <- 5 + as.integer(factor(ph$Env)) + u[match(ph$GID, rownames(G))] + rnorm(nrow(ph), sd = 0.5)
  ph$Yield[ph$GID %in% c("L01", "L02")] <- NA
  withr::local_envvar(PREDICTPRO_ASREML_COEF_MAX_ROWS = "1")
  withr::local_dir(withr::local_tempdir())
  out <- tryCatch(suppressMessages(suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~ Env,
    random = ~ GID + GID:Env, heter_groups = "Env", GS_model = "GBLUP", engine = "asreml",
    gmatrix_method = "VanRaden", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE
  ))), error = function(e) e)
  if (inherits(out, "error") && grepl("License|licen", conditionMessage(out))) skip("ASReml license unavailable")
  expect_false(inherits(out, "error"), info = if (inherits(out, "error")) conditionMessage(out) else "")
  pv <- out$model_results$predicted_values
  expect_true(is.data.frame(pv) && nrow(pv) == nrow(ph))
})
