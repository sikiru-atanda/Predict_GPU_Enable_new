test_that("multi-trait model reliability uses the exported genetic variance reference", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), times = 2L),
    Trait = rep(c("Yield", "Height"), each = 2L),
    Predicted_value = c(1.1, 1.3, 2.1, 2.3),
    Train_Test_Label = c("Train", "Test", "Train", "Test"),
    Observed_value = c(1, NA, 2, NA),
    Standard_error = sqrt(c(0.10, 0.20, 0.15, 0.30)),
    PEV = c(0.10, 0.20, 0.15, 0.30),
    Reliability = rep(0.99, 4L),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Trait = c("Yield", "Height"),
    Component = "genetic_variance",
    Components = c(0.50, 0.75),
    Standard_error = c(0.05, 0.06),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      variance_components = vc,
      # Complex Bayesian/ASReml result bundles can contain GP-like metadata;
      # the already-labelled public component table remains authoritative.
      gp_info = list(backend = "mock_gp_like_bundle")
    ),
    gen_name = "GID"
  )$predicted_values

  expected_reference <- c(0.50, 0.50, 0.75, 0.75)
  expected_reliability <- pmax(0, pmin(1, 1 - pred$PEV / expected_reference))
  expect_equal(out$Reliability_reference_variance, expected_reference)
  expect_equal(out$Reliability_variance_input, pred$PEV)
  expect_equal(out$Reliability, expected_reliability)
  expect_true(all(grepl("model-specific genetic variance", out$Reliability_basis, fixed = TRUE)))
})

test_that("specialized multi-trait routes reconcile tables without dropping fitted models", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), times = 2L),
    Trait = rep(c("T1", "T2"), each = 2L),
    Predicted_value = c(1.1, 1.2, 2.1, 2.2),
    Train_Test_Label = rep(c("Train", "Test"), 2L),
    Observed_value = c(1, NA, 2, NA),
    Standard_error = sqrt(c(0.1, 0.2, 0.15, 0.25)),
    PEV = c(0.1, 0.2, 0.15, 0.25),
    Reliability = NA_real_,
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Trait = c("T1", "T2"),
    Component = "genetic_variance",
    Components = c(0.5, 0.8),
    Standard_error = NA_real_,
    stringsAsFactors = FALSE
  )
  fitted_marker <- structure(list(id = "keep_me"), class = "mock_fitted_model")

  out <- PredictProR:::gp_reconcile_multitrait_model_output(
    list(
      predicted_values = pred,
      variance_components = vc,
      Asreml_model = fitted_marker
    ),
    gen_name = "GID"
  )

  expect_identical(out$Asreml_model, fitted_marker)
  expect_equal(
    out$predicted_values$Reliability_reference_variance,
    c(0.5, 0.5, 0.8, 0.8)
  )
  expect_equal(out$predicted_values$Reliability_variance_input, pred$PEV)
  expect_true(all(nzchar(out$predicted_values$Reliability_basis)))
})
