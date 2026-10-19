test_that("multi-trait output preserves genetic PEV separately from total mean uncertainty", {
  pheno <- data.frame(GID = c("a", "b", "c"), T1 = c(1, 2, NA), T2 = c(2, 4, NA))
  pred <- data.frame(gid = rep(pheno$GID, each = 2), trait = rep(c("T1", "T2"), 3),
                     Prediction = seq_len(6), SE_latent = rep(sqrt(1.5), 6),
                     Reliability_variance_input = rep(1, 6), PEV = rep(1.5, 6), Prediction_Var_latent = rep(1.5, 6))
  out <- PredictProR:::gp_model_execute_multitrait_predicted_values(
    predictions = pred, pheno_data = pheno, response = c("T1", "T2"), gen_name = "GID")
  expect_equal(out$PEV, rep(1.5, 6))
  expect_equal(out$Reliability_variance_input, rep(1, 6))
  expect_equal(out$Standard_error^2, rep(1.5, 6))
})

test_that("MT reliability uses the row genetic prior and retains its basis", {
  pheno <- data.frame(GID = c("a", "b", "c"), T1 = c(100, 200, NA), T2 = c(2, 4, NA))
  pred <- data.frame(gid = rep(pheno$GID, each = 2), trait = rep(c("T1", "T2"), 3),
                     Prediction = seq_len(6), SE_latent = 2, PEV = 4, Reliability_variance_input = 1,
                     Reliability_reference_variance = c(2, 4, 2, 4, 2, 4),
                     Reliability_basis = "genetic_prior_variance")
  out <- PredictProR:::gp_model_execute_multitrait_predicted_values(
    predictions = pred, pheno_data = pheno, response = c("T1", "T2"), gen_name = "GID")
  expect_equal(out$Reliability, rep(c(.5, .75), 3))
  expect_equal(out$Reliability_reference_variance, pred$Reliability_reference_variance)
  expect_equal(out$Reliability_basis, pred$Reliability_basis)
  # The final formatter must preserve K_ii-specific priors even if the
  # component table contains different, trait-level scalar variances.
  vc <- data.frame(Trait = c("T1", "T2"), Component = "genetic_variance", Components = c(9, 10))
  final <- PredictProR:::gp_standardize_multitrait_prediction_contract(out, vc, "GID")
  expect_equal(final$Reliability_reference_variance, pred$Reliability_reference_variance)
  expect_equal(final$Reliability_variance_input, pred$Reliability_variance_input)
  expect_equal(final$Reliability, out$Reliability)
  expect_equal(final$PEV, final$Standard_error^2)
  pred$Reliability_reference_variance <- NULL
  out <- PredictProR:::gp_model_execute_multitrait_predicted_values(
    predictions = pred, pheno_data = pheno, response = c("T1", "T2"), gen_name = "GID")
  expect_true(all(is.na(out$Reliability)))
})

test_that("single-trait adapter prefers explicit genetic PEV over total latent variance", {
  pheno <- data.frame(GID = c("a", "b", "c"), Yield = c(1, 2, NA))
  pred <- data.frame(Genotype = pheno$GID, Prediction = c(1.1, 2.1, 3.1),
                     Prediction_SE_latent = rep(sqrt(1.5), 3),
                     Reliability_variance_input = rep(1, 3), PEV = rep(1.5, 3), Prediction_Var_latent = rep(1.5, 3),
                     Reliability_reference_variance = rep(2, 3))
  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = pred, pheno_data = pheno, response = "Yield", gen_name = "GID", test_set = "c",
    gp_result = list(per_env = data.frame(Env = "ENV1", total_genetic_variance = 2)))
  expect_equal(out$PEV, rep(1.5, 3))
  expect_equal(out$Reliability_variance_input, rep(1, 3))
  expect_equal(out$Standard_error^2, rep(1.5, 3))
  expect_equal(out$Reliability, rep(.5, 3))
  expect_true(all(grepl("^genetic_prior_variance", out$Reliability_basis)))
})

test_that("KRR without genetic PEV cannot report total-variance reliability", {
  pheno <- data.frame(GID = c("a", "b"), Yield = c(1, NA))
  pred <- data.frame(Name = pheno$GID, Prediction = c(1, 2),
    Prediction_SE_latent = c(.5, .8), Prediction_Var_latent = c(.25, .64))
  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    pred, pheno, "Yield", "GID", test_set = "b")
  expect_equal(out$PEV, out$Standard_error^2)
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(is.na(out$Reliability_variance_input)))
  expect_true(all(nzchar(out$Reliability_remarks)))
})

test_that("parameter SE tables and conditioning survive public R formatting", {
  tab <- data.frame(component = c("genetic!T1:T1", "residual!T1:T1"),
    environment = "overall", trait1 = "T1", trait2 = "T1",
    estimate = c(2, 3), std.error = c(.4, .5), inference_status = "available")
  cor_tab <- data.frame(component = "genetic!cor(T1,T2)", environment = "overall",
    estimate = .5, std.error = .1, inference_status = "available")
  summary <- list(inference_status = "available", inference_available = TRUE,
    covariance_parameter_conditioning = "kernel_weights_and_standardization")
  native <- list(result = list(varcomp = tab, correlation_parameters = cor_tab,
    inference_summary = summary), predictions = data.frame())
  native <- PredictProR:::gp_public_normalize_multitrait_uncertainty(native)
  expect_equal(native$covariance_parameters$std.error, c(.4, .5))
  p <- data.frame(GID = letters[1:3], Predicted_value = 1:3,
    Train_Test_Label = c("Train", "Train", "Test"), Standard_error = .5,
    PEV = .25, Reliability_variance_input = .2, Reliability_reference_variance = 2,
    Reliability = .9, Reliability_basis = "genetic_prior_variance")
  public <- PredictProR:::gp_standardize_public_model_result(list(
    predicted_values = p, gp_result = native$result,
    variance_components = data.frame(Component = c("genetic_variance", "residual_variance"),
      Components = c(2,3), Standard_error = c(.4,.5))), model = "GP")
  expect_equal(public$covariance_parameters, tab)
  expect_equal(public$correlation_parameters, cor_tab)
  expect_equal(public$inference_summary, summary)
  # Single-trait/FA summaries must not be shadowed by the model-parameter
  # data frame's character column, which is also named "summary".
  st_native <- native$result
  st_native$inference_summary <- NULL
  st_native$summary <- c(summary, list(fa_derived = list(
    inference_basis = "delta_method_conditional_on_active_set_and_kernel_weights")))
  st <- PredictProR:::gp_standardize_public_model_result(list(
    model_parameters = data.frame(stat = "gp_model", summary = "GP_FA"),
    predicted_values = p, gp_result = st_native,
    variance_components = public$variance_components), model = "GP_FA")
  expect_equal(st$inference_summary$inference_status, "available")
  expect_equal(st$inference_summary$derived_se_method,
    "delta_method_conditional_on_active_set_and_kernel_weights")
})
