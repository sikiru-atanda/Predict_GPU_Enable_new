test_that("Bayesian MT-MET aggregate reliability retains its arithmetic provenance", {
  pred <- data.frame(
    GID = rep(c("G1", "G2"), each = 2L),
    Env = rep(c("E1", "E2"), times = 2L),
    Trait = "TraitA",
    Predicted_value = c(1.0, 1.4, 0.2, 0.6),
    Train_Test_Label = c("Train", "Train", "Test", "Test"),
    Observed_value = c(0.9, 1.5, NA_real_, NA_real_),
    Standard_error = c(0.2, 0.3, 0.4, 0.5),
    PEV = c(0.04, 0.09, 0.16, 0.25),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::bayes_multitrait_across_environment_prediction(
    pred_long = pred,
    trait_genetic_variance = c(TraitA = 2)
  )

  expect_equal(nrow(out), 2L)
  expect_true(all(c(
    "Reliability_variance_input", "Reliability_reference_variance",
    "Reliability_basis", "PEV_basis", "Prediction_uncertainty_source"
  ) %in% names(out)))
  expect_equal(out$Reliability_variance_input, out$PEV, tolerance = 1e-12)
  expect_equal(out$Reliability_reference_variance, rep(2, 2L), tolerance = 1e-12)
  expect_equal(out$Reliability, pmax(0, pmin(1, 1 - out$PEV / 2)), tolerance = 1e-12)
  expect_true(all(grepl("model-specific genetic variance", out$Reliability_basis, fixed = TRUE)))
  expect_true(all(grepl("independence approximation", out$PEV_basis, fixed = TRUE)))
  expect_identical(
    out$Prediction_uncertainty_source,
    rep("posterior_marginal_independence_approximation", 2L)
  )
})
