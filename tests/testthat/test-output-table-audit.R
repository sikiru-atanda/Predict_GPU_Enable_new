# gp_audit_output_tables() is what the scenario suite runs on every result.
# Each bug class found by the 2026-10-19 audit must make it fail; a correct
# result must pass every check.

audit_fixture <- function() {
  set.seed(3)
  envs <- c("E1", "E2", "E3")
  ids <- sprintf("L%02d", 1:20)
  pv <- expand.grid(GID = ids, Env = envs, stringsAsFactors = FALSE)
  g <- rnorm(20)
  pv$Observed_value <- 10 + c(E1 = 0, E2 = 2, E3 = 4)[pv$Env] + g[match(pv$GID, ids)] + rnorm(60, sd = 0.5)
  pv$Predicted_value <- pv$Observed_value + rnorm(60, sd = 0.3)
  pv$Train_Test_Label <- "Train"
  pv$Standard_error <- 0.4
  pv$PEV <- 0.16
  pv$Prediction_interval_method <- "Gaussian approximation from model-reported standard error"
  pv$Prediction_interval_nominal_coverage <- 0.95
  pv$lower_bound <- pv$Predicted_value - qnorm(0.975) * 0.4
  pv$upper_bound <- pv$Predicted_value + qnorm(0.975) * 0.4
  pv$Reliability_variance_input <- 0.1
  pv$Reliability_reference_variance <- 1
  pv$Reliability <- 0.9
  pv$Reliability_basis <- "genetic_prior_variance: 1 - genetic PEV / (sum of sigma2_g * K_ii)"
  pv$Reliability_remarks <- "Reliable"
  std <- PredictProR::prediction_output_standard(response_family = "gaussian", task = "multi_environment",
                                                 gen_name = "GID", heter_groups = "Env")
  for (col in setdiff(std$prediction_columns, names(pv))) {
    pv[[col]] <- if (grepl("remarks|basis|method|source", col, ignore.case = TRUE)) NA_character_ else NA_real_
  }
  pv <- pv[, std$prediction_columns]
  vc <- data.frame(
    Component = c(paste0("genetic_variance_", envs), paste0("residual_variance_", envs), paste0("heritability_", envs)),
    Components = c(rep(1, 3), rep(0.25, 3), rep(0.8, 3)),
    Standard_error = 0.1, stringsAsFactors = FALSE
  )
  rownames(vc) <- vc$Component
  G <- matrix(0.6, 3, 3, dimnames = list(envs, envs)); diag(G) <- 1
  long <- function(M, name, comp) data.frame(Matrix = name, Component = comp, Row = rep(rownames(M), 3),
                                             Column = rep(colnames(M), each = 3), Value = as.vector(M),
                                             Estimation_method = "model_reported", stringsAsFactors = FALSE)
  list(model_results = structure(list(
    model_parameters = data.frame(stat = "gp_model", summary = "GBLUP"),
    predicted_values = pv, diagnostic_plots = list(), variance_components = vc,
    covariance_components = long(G, "Genetic_covariance_environments", "environment_covariance"),
    correlation_components = long(cov2cor(G), "Genetic_correlation_environments", "environment_correlation")
  ), class = c("predictpror_model_result", "list")))
}

failed_checks <- function(x) {
  a <- PredictProR:::gp_audit_output_tables(x)
  unique(a$check[a$ok %in% FALSE])
}

test_that("a correct result passes every table check", {
  a <- PredictProR:::gp_audit_output_tables(audit_fixture())
  expect_gt(nrow(a), 20)
  expect_identical(a$check[a$ok %in% FALSE], character())
})

test_that("a doubled intercept and a missing environment effect fail the level checks", {
  x <- audit_fixture()
  x$model_results$predicted_values$Predicted_value <- x$model_results$predicted_values$Predicted_value + 10
  expect_true("prediction level matches observed (training rows)" %in% failed_checks(x))
  y <- audit_fixture()
  pv <- y$model_results$predicted_values
  y$model_results$predicted_values$Predicted_value <- pv$Predicted_value - c(E1 = 0, E2 = 2, E3 = 4)[pv$Env]
  expect_true("prediction level matches observed in each environment" %in% failed_checks(y))
})

test_that("collapsed labels, broken correlations and odd environment names fail", {
  x <- audit_fixture()
  x$model_results$covariance_components$Row <- "not_identifiable"
  expect_true("real row/column labels" %in% failed_checks(x))
  y <- audit_fixture()
  y$model_results$correlation_components$Value[1] <- 0.7
  expect_true(any(c("correlation diagonal = 1", "symmetric") %in% failed_checks(y)))
  z <- audit_fixture()
  z$model_results$variance_components$Component <- sub("_E1$", "_common", z$model_results$variance_components$Component)
  expect_true("per-environment rows name prediction environments" %in% failed_checks(z))
})

test_that("identity, range and boundary violations fail", {
  x <- audit_fixture()
  x$model_results$predicted_values$PEV[1] <- 1
  expect_true("PEV = SE^2" %in% failed_checks(x))
  y <- audit_fixture()
  y$model_results$variance_components["heritability_E2", "Components"] <- 0.5
  expect_true("heritability = genetic / (genetic + residual)" %in% failed_checks(y))
  z <- audit_fixture()
  z$model_results$variance_components$Boundary <- ""
  z$model_results$variance_components["residual_variance_E2", c("Components", "Boundary")] <- list(NA, "not_estimable")
  z$model_results$variance_components["genetic_variance_E2", "Components"] <- NA
  z$model_results$variance_components["heritability_E2", "Components"] <- NA
  expect_true("no reliability where the residual is at zero" %in% failed_checks(z))
  w <- audit_fixture()
  w$model_results$predicted_values$Standard_error <- 20
  w$model_results$predicted_values$PEV <- 400
  w$model_results$predicted_values$lower_bound <- w$model_results$predicted_values$Predicted_value - qnorm(0.975) * 20
  w$model_results$predicted_values$upper_bound <- w$model_results$predicted_values$Predicted_value + qnorm(0.975) * 20
  expect_true("model-based SE below twice the phenotypic SD" %in% failed_checks(w))
})

test_that("a covariance matrix without its correlation matrix fails", {
  x <- audit_fixture()
  x$model_results$correlation_components <- NULL
  expect_true("every covariance matrix has a correlation matrix" %in% failed_checks(x))
})
