test_that("Gaussian single-environment prediction contract is compact and ordered", {
  pred <- data.frame(
    NAME = c("g1", "g2"),
    Predicted_value = c(1.2, 2.4),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.1, NA),
    Standard_error = c(0.1, 0.2),
    Genetic_variance = c(3, 3),
    Reliability_percentage = c(100, 100),
    Prediction_stability = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred, gen_name = "NAME")

  # Output column carries the user-supplied gen_name (here "NAME"), not a
  # hardcoded "GID" -- the user's variable name is respected end-to-end.
  expect_identical(names(out), PredictProR:::gp_gaussian_prediction_columns(gen_name = "NAME"))
  expect_identical(out$NAME, c("g1", "g2"))
  expect_equal(out$PEV, c(0.01, 0.04))
  expect_false("Genetic_variance" %in% names(out))
  expect_false("Reliability_percentage" %in% names(out))
  expect_equal(out$Prediction_stability, c(0.9, 0.8))
})

test_that("specialized result indexing preserves completed export status", {
  handled <- list(
    model_results = list(
      predicted_values = data.frame(
        GID = c("g1", "g2"),
        Predicted_value = c(1.2, 1.8),
        Train_Test_Label = c("Train", "Test"),
        Observed_value = c(1.1, NA_real_),
        Standard_error = c(0.1, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    summary_statistic = list(hybrid_component_summary = data.frame(component = "female")),
    test_diagonistic_plots = structure(list(), class = "gtable"),
    export_status = "Successful"
  )

  out <- PredictProR:::gp_index_specialized_model_execute_result(
    handled_result = handled,
    model_label = "RandomForest",
    trait_name = "multi_trait",
    gen_name = "GID"
  )

  expect_identical(out$export_status, "Successful")
  expect_identical(out$summary_statistic, handled$summary_statistic)
  expect_identical(out$test_diagonistic_plots, handled$test_diagonistic_plots)
})

test_that("Gaussian true-prediction contract derives bounds and reliability remarks", {
  pred <- data.frame(
    NAME = c("g1", "g2", "g3"),
    Predicted_value = c(1.2, 2.4, 3.1),
    Train_Test_Label = c("Train", "Test", "Test"),
    Observed_value = c(1.1, NA, NA),
    Standard_error = c(0.1, 0.2, 0.3),
    Reliability = c(0.95, 0.65, 0.25),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(predicted_values = pred),
    gen_name = "NAME"
  )
  public <- out$predicted_values

  # User-supplied gen_name ("NAME") is preserved in the output column order.
  expect_identical(names(public), PredictProR:::gp_gaussian_prediction_columns(gen_name = "NAME"))
  expect_true(all(is.finite(public$Standard_error)))
  expect_true(all(is.finite(public$PEV)))
  expect_true(all(is.finite(public$lower_bound)))
  expect_true(all(is.finite(public$upper_bound)))
  expect_true(all(nzchar(public$PEV_basis)))
  expect_true(all(public$Prediction_uncertainty_source == "model_reported_prediction_uncertainty"))
  expect_true(all(public$Prediction_interval_method ==
                    "Gaussian approximation from model-reported standard error"))
  expect_true(all(public$Prediction_interval_nominal_coverage == 0.95))
  expect_true(all(is.na(public$Prediction_interval_calibration_n)))
  expect_identical(public$Reliability_remarks, c("Reliable", "Acceptable", "Unreliable"))
})

test_that("public Gaussian results restore known training observations by model keys", {
  pheno <- data.frame(
    Line = c("g1", "g1", "g2", "g3"),
    Site = c("A", "A", "B", "B"),
    Yield = c(1, 3, 5, NA_real_),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    Line = c("g1", "g2", "g3"),
    Site = c("A", "B", "B"),
    Predicted_value = c(1.8, 4.7, 6.1),
    Train_Test_Label = c("Train", "Train", "Test"),
    Observed_value = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(predicted_values = pred),
    gen_name = "Line",
    heter_groups = "Site",
    pheno_data = pheno,
    response = "Yield"
  )

  expect_equal(out$predicted_values$Observed_value, c(2, 5, NA_real_))
})

test_that("observed-value restoration never discloses supplied values on Test rows", {
  pheno <- data.frame(
    GID = c("g1", "g2"),
    Yield = c(1.25, 9.5),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 8.8),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = NA_real_,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(predicted_values = pred),
    gen_name = "GID",
    pheno_data = pheno,
    response = "Yield"
  )

  expect_equal(out$predicted_values$Observed_value, c(1.25, NA_real_))
})

test_that("Gaussian formatter preserves backend-specific interval provenance", {
  pred <- data.frame(
    GID = "g1",
    Predicted_value = 1.2,
    Standard_error = 0.1,
    PEV = 0.01,
    PEV_basis = "posterior target variance",
    Prediction_uncertainty_source = "posterior draws",
    Prediction_interval_method = "posterior quantiles",
    Prediction_interval_nominal_coverage = 0.9,
    lower_bound = 0.9,
    upper_bound = 1.5,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred, confidence_level = 0.95)

  expect_identical(out$PEV_basis, "posterior target variance")
  expect_identical(out$Prediction_uncertainty_source, "posterior draws")
  expect_identical(out$Prediction_interval_method, "posterior quantiles")
  expect_identical(out$Prediction_interval_nominal_coverage, 0.9)
  expect_equal(out$lower_bound, 0.9)
  expect_equal(out$upper_bound, 1.5)
})

test_that("Gaussian formatter aliases supplied stability when reliability is blank", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.2, 2.4),
    Train_Test_Label = c("Train", "Test"),
    Standard_error = c(0.2, 0.3),
    Reliability = c(NA_real_, NA_real_),
    Reliability_remarks = c(NA_character_, NA_character_),
    Prediction_stability = c(0.82, 0.41),
    Prediction_stability_remarks = c("Stable", "Unstable"),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred)

  expect_identical(out$Reliability, c(0.82, 0.41))
  expect_identical(out$Prediction_stability, c(0.82, 0.41))
  expect_identical(out$Prediction_stability_remarks, c("Stable", "Unstable"))
  expect_true(all(grepl("not genetic reliability", out$Reliability_basis, fixed = TRUE)))
  plots <- PredictProR:::gp_gaussian_diagnostic_plots(out)
  expect_match(plots$predicted_vs_reliability$labels$x, "non-genetic", fixed = TRUE)
})

test_that("Gaussian MET contract adds Env and drops constant Trait", {
  pred <- data.frame(
    GID = c("g1", "g1"),
    Env = c("E1", "E2"),
    Trait = c("YLD", "YLD"),
    Predicted_value = c(1.2, 1.5),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.1, NA),
    Standard_error = c(0.1, 0.2),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred)

  expect_identical(names(out), PredictProR:::gp_gaussian_prediction_columns(include_env = TRUE))
})

test_that("Gaussian formatter recomputes flat reliability within environments when PEV varies", {
  pred <- data.frame(
    GID = paste0("g", seq_len(6)),
    Env = rep(c("E1", "E2"), each = 3),
    Predicted_value = seq_len(6),
    Train_Test_Label = "Test",
    Standard_error = sqrt(c(0.1, 0.2, 0.35, 0.08, 0.1, 0.12)),
    PEV = c(0.1, 0.2, 0.35, 0.08, 0.1, 0.12),
    Reliability_reference_variance = c(rep(1.2, 3), rep(0.8, 3)),
    Reliability = c(rep(0, 3), rep(0.75, 3)),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred)

  expect_true("Env" %in% names(out))
  expect_gt(length(unique(round(out$Reliability[out$Env == "E1"], 10L))), 1L)
  expect_gt(length(unique(round(out$Reliability[out$Env == "E2"], 10L))), 1L)
  expect_true(all(out$Reliability > 0 & out$Reliability <= 1))
})

test_that("Gaussian formatter repairs flat reliability after reference variance columns are dropped", {
  pred <- data.frame(
    GID = paste0("g", seq_len(6)),
    Env = "E4",
    Predicted_value = c(5.1, 5.3, 5.84, 5.75, 5.46, 6.19),
    Train_Test_Label = c("Train", "Train", rep("Test", 4)),
    Standard_error = sqrt(c(0.12, 0.18, 0.210, 0.268, 0.2684, 0.0952)),
    PEV = c(0.12, 0.18, 0.210, 0.268, 0.2684, 0.0952),
    Reliability = c(0.41, 0.55, rep(0, 4)),
    Reliability_remarks = c("Unreliable", "Acceptable", rep("Unreliable", 4)),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred)
  tst <- out[out$Train_Test_Label == "Test", , drop = FALSE]

  expect_gt(length(unique(round(tst$Reliability, 10L))), 1L)
  expect_true(all(out$Reliability > 0 & out$Reliability <= 1))
  expect_true(all(nzchar(out$Reliability_remarks)))
})

test_that("GP MET prediction contract preserves Env and across-environment totals", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), each = 2),
    Env = rep(c("E1", "E2"), times = 2),
    Predicted_value = c(1.0, 1.4, 2.0, 2.2),
    Train_Test_Label = c("Train", "Test", "Train", "Train"),
    Observed_value = c(1.1, NA, 2.1, 2.3),
    Standard_error = c(0.2, 0.3, 0.2, 0.25),
    PEV = c(0.04, 0.09, 0.04, 0.0625),
    lower_bound = c(0.6, 0.8, 1.6, 1.7),
    upper_bound = c(1.4, 2.0, 2.4, 2.7),
    Uncertainty = c(0.8, 1.2, 0.8, 1.0),
    Uncertainty_remarks = "model_based",
    Reliability = c(0.9, 0.8, 0.95, 0.92),
    Reliability_remarks = "Reliable",
    stringsAsFactors = FALSE
  )

  public <- PredictProR:::gp_public_prediction_table(pred)
  expect_true("Env" %in% names(public))
  expect_identical(public$Env, pred$Env)

  total <- PredictProR:::gp_bridge_total_prediction_output(pred)
  expect_s3_class(total, "data.frame")
  expect_false("Env" %in% names(total))
  expect_equal(nrow(total), 2L)
  expect_true(all(c("GID", "Predicted_value", "Standard_error", "PEV", "Reliability") %in% names(total)))
  expect_true("Test" %in% total$Train_Test_Label)
})

test_that("GP model-execute MET adapter returns environment rows and across totals", {
  pheno <- data.frame(
    GID = rep(c("g1", "g2"), each = 2),
    Env = rep(c("E1", "E2"), times = 2),
    Yield = c(1.0, NA, 2.0, 2.2),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    GID = pheno$GID,
    Env = pheno$Env,
    Predicted_value = c(1.1, 1.4, 2.1, 2.3),
    Prediction_SE = c(0.2, 0.3, 0.2, 0.25),
    Prediction_variance = c(0.04, 0.09, 0.04, 0.0625),
    stringsAsFactors = FALSE
  )
  K <- diag(2)
  rownames(K) <- colnames(K) <- c("g1", "g2")

  local_mocked_bindings(
    gp_single_trait_model = function(...) {
      list(
        predictions = pred,
        result = list(varcomp = data.frame(component = "genetic", estimate = 1)),
        info = list(backend = "mock")
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_backend_gaussian_model(
    model_name = "KRR",
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    gmatrix = K,
    gp_output_level = "predict_with_se",
    gp_prediction_output = "all"
  )

  public <- PredictProR:::gp_standardize_public_model_result(out, gen_name = "GID")
  expect_true("Env" %in% names(public$predicted_values))
  expect_identical(public$predicted_values$Env, pheno$Env)
  expect_s3_class(public$across_environment_predicted_values, "data.frame")
  expect_false("Env" %in% names(public$across_environment_predicted_values))
  expect_equal(nrow(public$across_environment_predicted_values), length(unique(pheno$GID)))
})

test_that("Gaussian MET true-prediction contract derives remarks for environment and across outputs", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), each = 2),
    Env = rep(c("E1", "E2"), times = 2),
    Predicted_value = c(1.0, 1.4, 2.0, 2.2),
    Train_Test_Label = c("Train", "Test", "Train", "Test"),
    Observed_value = c(1.1, NA, 2.1, NA),
    Standard_error = c(0.2, 0.3, 0.2, 0.25),
    Reliability = c(0.9, 0.8, 0.45, 0.55),
    stringsAsFactors = FALSE
  )
  total <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.2, 2.1),
    Train_Test_Label = c("Test", "Test"),
    Observed_value = c(NA, NA),
    PEV = c(0.04, 0.09),
    Reliability = c(0.75, 0.35),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      across_environment_predicted_values = total
    ),
    gen_name = "GID"
  )

  expect_identical(
    names(out$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns(include_env = TRUE)
  )
  expect_identical(
    names(out$across_environment_predicted_values),
    PredictProR:::gp_gaussian_prediction_columns()
  )
  expect_true(all(is.finite(out$predicted_values$lower_bound)))
  expect_true(all(is.finite(out$predicted_values$upper_bound)))
  expect_true(all(nzchar(out$predicted_values$Reliability_remarks)))
  expect_true(all(is.finite(out$across_environment_predicted_values$Standard_error)))
  expect_true(all(is.finite(out$across_environment_predicted_values$lower_bound)))
  expect_true(all(is.finite(out$across_environment_predicted_values$upper_bound)))
  expect_true(all(nzchar(out$across_environment_predicted_values$Reliability_remarks)))
})

test_that("GP model-execute MT-MET adapter returns the common Gaussian columns", {
  pheno <- data.frame(
    GID = rep(c("g1", "g2"), each = 2),
    Env = rep(c("E1", "E2"), times = 2),
    T1 = c(1.0, NA, 2.0, 2.2),
    T2 = c(2.0, 2.1, NA, 2.8),
    stringsAsFactors = FALSE
  )
  pred <- expand.grid(
    gid = c("g1", "g2"),
    env = c("E1", "E2"),
    trait = c("T1", "T2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred$Prediction <- seq_len(nrow(pred)) / 10
  pred$Prediction_SE_observed <- 0.2
  pred$Prediction_Var_observed <- 0.04
  pred$Prediction_SE_latent <- 0.1
  pred$Prediction_Var_latent <- 0.01

  out <- PredictProR:::gp_model_execute_multitrait_predicted_values(
    predictions = pred,
    pheno_data = pheno,
    response = c("T1", "T2"),
    gen_name = "GID",
    heter_groups = "Env"
  )

  expect_identical(
    names(PredictProR:::gp_format_gaussian_prediction_table(out)),
    PredictProR:::gp_gaussian_prediction_columns(include_env = TRUE, include_trait = TRUE)
  )
  expect_true(all(c("lower_bound", "upper_bound", "Uncertainty", "Reliability") %in% names(out)))
  expect_equal(out$Standard_error, rep(0.1, nrow(out)))
  expect_equal(out$PEV, rep(0.01, nrow(out)))

  total <- PredictProR:::gp_multitrait_across_environment_prediction(out, gen_name = "GID")
  public_total <- PredictProR:::gp_format_total_prediction_table(total, gen_name = "GID", include_trait = TRUE)
  expect_identical(
    names(public_total),
    PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
  )
  expect_false("Env" %in% names(public_total))
  expect_equal(nrow(public_total), 4L)
  expect_true(all(is.finite(public_total$Standard_error)))
  expect_true(all(is.finite(public_total$PEV)))
  expect_true(all(is.finite(public_total$Reliability)))
})

test_that("BGLR joint multi-trait Bayesian helper uses documented Cov and common columns", {
  skip_if_not_installed("BGLR")
  set.seed(27)
  ids <- paste0("g", seq_len(6))
  pheno <- data.frame(
    GID = ids,
    T1 = c(1.1, 1.8, NA, 2.7, 3.0, 3.4),
    T2 = c(2.2, NA, 2.8, 3.3, 3.8, 4.1),
    stringsAsFactors = FALSE
  )
  K <- diag(length(ids))
  rownames(K) <- colnames(K) <- ids
  K2 <- diag(seq(0.7, 1.3, length.out = length(ids)))
  rownames(K2) <- colnames(K2) <- ids

  out <- PredictProR:::bayes_multitrait_joint_fit(
    pheno_data = pheno,
    response = c("T1", "T2"),
    gen_name = "GID",
    heter_groups = NULL,
    kernels = list(gmatrix = K, transcript_kernel = K2),
    GS_model = "RKHS",
    bayes_para = list(nIter = 45L, burnIn = 15L, thin = 2L),
    verbose = FALSE
  )

  pred <- out$bayes_result$Predicted_value
  expect_identical(
    names(PredictProR:::gp_format_gaussian_prediction_table(pred)),
    PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
  )
  expect_true(all(c("T1", "T2") %in% unique(pred$Trait)))
  expect_identical(out$bayes_result$model_parameters$eta_cov_argument[[1]], "Cov")
  expect_identical(
    sort(names(out$bayes_result$Genetic_covariance_by_kernel)),
    c("gmatrix", "transcript_kernel")
  )
  expect_equal(
    Reduce(`+`, out$bayes_result$Genetic_covariance_by_kernel),
    out$bayes_result$Genetic_covariance_traits,
    tolerance = 1e-10
  )
  expect_true(all(out$bayes_result$kernel_variance_components$Independent_kernel_estimate))
  expect_true(all(is.finite(out$bayes_result$kernel_variance_components$Standard_error)))
})

test_that("Joint multi-trait ML/DL keeps Reliability unavailable with PEV and retains stability", {
  pred <- data.frame(
    GID = rep(paste0("g", 1:6), times = 2),
    Trait = rep(c("T1", "T2"), each = 6),
    Predicted_value = c(1.0, 1.4, 1.9, 2.5, 3.0, 3.8, 4.1, 4.4, 4.9, 5.5, 6.0, 6.7),
    Train_Test_Label = rep(c("Train", "Train", "Train", "Train", "Test", "Test"), times = 2),
    Observed_value = c(1.1, 1.1, 2.2, 2.4, NA, NA, 4.0, 4.8, 4.5, 5.7, NA, NA),
    Prediction_confidence = rep(0.6, 12),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(predicted_values = pred),
    gen_name = "GID"
  )

  expect_identical(
    names(out$predicted_values),
    PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
  )
  expect_true(all(is.na(out$predicted_values$Standard_error)))
  expect_true(all(is.na(out$predicted_values$PEV)))
  expect_true(all(is.na(out$predicted_values$lower_bound)))
  expect_true(all(is.na(out$predicted_values$upper_bound)))
  expect_true(all(is.na(out$predicted_values$Reliability)))
  expect_true(all(is.na(out$predicted_values$Reliability_variance_input)))
  expect_true(all(is.finite(out$predicted_values$Prediction_stability)))
  expect_true(all(grepl("not genetic reliability", out$predicted_values$Reliability_basis, fixed = TRUE)))
})

test_that("Joint multi-trait ML risk calibration accepts user ID columns from data.table CV output", {
  skip_if_not_installed("data.table")
  cv_pred <- data.table::data.table(
    ID = rep(paste0("g", 1:18), times = 6),
    Trait = rep(rep(c("Yield", "Protein", "Oil"), each = 18), times = 2),
    model = "RandomForest",
    Observed_value = rep(seq_len(18), times = 6),
    Predicted_value = rep(rev(seq_len(18)), times = 6) + rep(c(0, 0.1, 0.2), each = 18, times = 2),
    Train_Test_Label = "Test",
    fold = rep(rep(1:2, each = 54), each = 1),
    rep = 1L,
    cv_role = "test"
  )

  out <- PredictProR:::gp_multitrait_dl_rank_risk_summary(cv_pred)

  expect_s3_class(out$combined_risk, "data.frame")
  expect_equal(nrow(out$combined_risk), nrow(cv_pred))
  expect_true(all(c("GID", "Trait", "model", "rank_instability_risk") %in% names(out$combined_risk)))
  expect_true(all(is.finite(out$combined_risk$rank_instability_risk)))
  expect_true(all(c("Yield", "Protein", "Oil") %in% out$aggregated_risk$trait))
})

test_that("Joint multi-trait ASReml-style outputs derive reliability from REML trait variances", {
  pred <- data.frame(
    GID = rep(paste0("g", 1:4), times = 2),
    Trait = rep(c("T1", "T2"), each = 4),
    Predicted_value = c(1.0, 1.6, 2.1, 2.7, 3.0, 3.8, 4.2, 5.1),
    Train_Test_Label = rep(c("Train", "Train", "Test", "Test"), times = 2),
    Observed_value = c(1.2, 1.4, NA, NA, 3.1, 3.6, NA, NA),
    Standard_error = c(0.10, 0.14, 0.20, 0.25, 0.16, 0.18, 0.24, 0.30),
    PEV = c(0.01, 0.0196, 0.04, 0.0625, 0.0256, 0.0324, 0.0576, 0.09),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Trait = c("T1", "T1", "T2", "T2"),
    Component = rep(c("genetic_variance", "residual_variance"), times = 2),
    Components = c(0.8, 0.2, 1.4, 0.3),
    Standard_error = c(0.12, 0.05, 0.18, 0.06),
    Estimation_method = "ASReml_REML",
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(predicted_values = pred, variance_components = vc),
    gen_name = "GID"
  )

  expect_equal(out$predicted_values$PEV, pred$PEV)
  expect_equal(out$predicted_values$Standard_error, pred$Standard_error)
  expect_true(all(is.finite(out$predicted_values$Reliability)))
  expect_true(all(out$predicted_values$Reliability >= 0 & out$predicted_values$Reliability <= 1))
  expect_false("Estimation_method" %in% names(out$variance_components))
})

test_that("Variance, covariance, and correlation outputs use common public schemas", {
  pred <- data.frame(
    GID = c("g1", "g2", "g3"),
    Predicted_value = c(1.1, 1.8, 2.5),
    Train_Test_Label = c("Train", "Train", "Test"),
    Observed_value = c(1.0, 2.0, NA),
    Standard_error = c(0.10, 0.20, 0.25),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(1.2, 0.4, 0.75),
    Standard_error = c(0.2, 0.1, 0.05),
    Estimation_method = "REML",
    stringsAsFactors = FALSE
  )
  cov_mat <- matrix(c(1.2, 0.3, 0.3, 0.8), nrow = 2)
  rownames(cov_mat) <- colnames(cov_mat) <- c("T1", "T2")
  cor_mat <- stats::cov2cor(cov_mat)

  out <- PredictProR:::gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      variance_components = vc,
      Genetic_covariance_traits = cov_mat,
      Genetic_correlation_traits = cor_mat
    ),
    gen_name = "GID"
  )

  expect_identical(names(out$variance_components), PredictProR:::gp_variance_component_columns())
  expect_false("Trait" %in% names(out$variance_components))
  expect_false("Env" %in% names(out$variance_components))
  expect_identical(names(out$covariance_components), PredictProR:::gp_covariance_component_columns())
  expect_identical(names(out$correlation_components), PredictProR:::gp_covariance_component_columns())
  expect_equal(nrow(out$covariance_components), 4L)
  expect_equal(nrow(out$correlation_components), 4L)
  expect_true(all(out$covariance_components$Component == "genetic_covariance"))
  expect_true(all(out$correlation_components$Component == "genetic_correlation"))
  expect_true(is.matrix(out$Genetic_covariance_traits))
  expect_true(is.matrix(out$Genetic_correlation_traits))
})

test_that("Single-level covariance matrices are omitted from public outputs", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 1.8),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.0, NA_real_),
    Standard_error = c(0.10, 0.20),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(1.2, 0.8, 0.6),
    Standard_error = c(0.2, 0.1, 0.05),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_public_model_result(
    list(
      predicted_values = pred,
      variance_components = vc,
      env_covariance_fa = matrix(1.2, nrow = 1L),
      env_correlation_fa = matrix(1, nrow = 1L)
    ),
    gen_name = "GID"
  )

  expect_false("Genetic_covariance_environments" %in% names(out))
  expect_false("Genetic_correlation_environments" %in% names(out))
  expect_false("covariance_components" %in% names(out))
  expect_false("correlation_components" %in% names(out))
})

test_that("Generic ASReml covariance and correlation lists are normalized", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c(1.1, 1.8),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.0, NA),
    Standard_error = c(0.10, 0.20),
    stringsAsFactors = FALSE
  )
  vc <- data.frame(
    Components = c(1.1, 0.5, 0.69),
    Standard_error = NA_real_,
    row.names = c("genetic_variance", "residual_variance", "heritability")
  )
  cov_mat <- matrix(c(1.1, 0.2, 0.2, 0.7), nrow = 2)
  rownames(cov_mat) <- colnames(cov_mat) <- c("E1", "E2")
  cor_mat <- stats::cov2cor(cov_mat)

  out <- PredictProR:::gp_standardize_public_model_result(
    list(
      Predicted_value = pred,
      Variance_components = vc,
      Covariance = list(site = cov_mat),
      Correlation = list(site = cor_mat)
    ),
    gen_name = "GID"
  )

  expect_identical(names(out$variance_components), PredictProR:::gp_variance_component_columns())
  expect_false("Env" %in% names(out$variance_components))
  expect_identical(names(out$covariance_components), PredictProR:::gp_covariance_component_columns())
  expect_identical(names(out$correlation_components), PredictProR:::gp_covariance_component_columns())
  expect_false("Estimation_method" %in% names(out$variance_components))
  expect_equal(unique(out$covariance_components$Matrix), "site")
  expect_equal(unique(out$correlation_components$Matrix), "site")
})

test_that("Classification prediction contract uses class probabilities and no Gaussian uncertainty columns", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c("Resistant", "Susceptible"),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c("Resistant", NA),
    Prediction_confidence = c(0.8, 0.7),
    Classification_uncertainty = c(0.2, 0.3),
    Reliability = c(0.8, 0.7),
    Prob_Resistant = c(0.8, 0.3),
    Prob_Susceptible = c(0.2, 0.7),
    Standard_error = c(0.1, 0.2),
    PEV = c(0.01, 0.04),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_classification_prediction_table(pred)

  expect_true(all(c("Probability_Resistant", "Probability_Susceptible") %in% names(out)))
  expect_false("Standard_error" %in% names(out))
  expect_false("PEV" %in% names(out))
  expect_identical(out$Predicted_class, c("Resistant", "Susceptible"))
  expect_identical(out$Observed_class, c("Resistant", NA_character_))
})

test_that("MET classification prediction contract preserves Env", {
  pred <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Predicted_value = c("yes", "no", "yes", "yes"),
    Train_Test_Label = c("Train", "Test", "Train", "Train"),
    Observed_value = c("yes", NA, "no", "yes"),
    Prediction_confidence = c(0.8, 0.7, 0.6, 0.9),
    Classification_uncertainty = c(0.2, 0.3, 0.4, 0.1),
    Reliability = c(0.8, 0.7, 0.6, 0.9),
    Prob_no = c(0.2, 0.7, 0.4, 0.1),
    Prob_yes = c(0.8, 0.3, 0.6, 0.9),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(pred)

  expect_identical(
    names(out),
    c(
      PredictProR:::gp_public_prediction_columns("binary", include_env = TRUE),
      "Probability_no", "Probability_yes"
    )
  )
  expect_identical(out$Env, pred$Env)
  expect_true(all(c("Probability_no", "Probability_yes") %in% names(out)))
  expect_equal(out$Prediction_confidence, apply(as.matrix(out[, c("Probability_no", "Probability_yes")]), 1L, max))
  expect_equal(out$Classification_uncertainty, 1 - out$Prediction_confidence)
  expect_equal(out$Reliability, out$Prediction_confidence)
  expect_false(any(c("Standard_error", "PEV") %in% names(out)))
})

test_that("MET classification summaries include across-environment totals", {
  pheno <- data.frame(
    GID = c("g1", "g1", "g2", "g2"),
    Env = c("E1", "E2", "E1", "E2"),
    Class = factor(c("yes", NA, "no", "yes"), levels = c("no", "yes")),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    GID = pheno$GID,
    Env = pheno$Env,
    Predicted_value = c(0.8, 0.3, 0.6, 0.9),
    Train_Test_Label = c("Train", "Test", "Train", "Train"),
    Predicted_class = c("yes", "no", "yes", "yes"),
    Prediction_confidence = c(0.8, 0.7, 0.6, 0.9),
    Classification_uncertainty = c(0.2, 0.3, 0.4, 0.1),
    Prob_no = c(0.2, 0.7, 0.4, 0.1),
    Prob_yes = c(0.8, 0.3, 0.6, 0.9),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_met_summary_statistics(
    predicted_values = pred,
    pheno_data = pheno,
    response = "Class",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "RandomForest",
    response_family = "binary"
  )

  expect_s3_class(out$MET_prediction_table, "data.frame")
  expect_true("Env" %in% names(out$MET_prediction_table))
  expect_s3_class(out$Total_Predicted_value, "data.frame")
  expect_false("Env" %in% names(out$Total_Predicted_value))
  expect_equal(nrow(out$Total_Predicted_value), length(unique(pheno$GID)))
  expect_true(all(c(
    "GID", "Predicted_class", "Train_Test_Label", "Observed_class",
    "Prediction_confidence", "Classification_uncertainty",
    "Reliability", "Reliability_remarks", "Probability_no", "Probability_yes"
  ) %in% names(out$Total_Predicted_value)))
  expect_true("Test" %in% out$Total_Predicted_value$Train_Test_Label)
})

test_that("classification labels survive generic prediction standardization", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Predicted_value = c("Resistant", "Susceptible"),
    Train_Test_Label = c("Train", "Test"),
    Prob_Resistant = c(0.8, 0.3),
    Prob_Susceptible = c(0.2, 0.7),
    stringsAsFactors = FALSE
  )

  standardized <- PredictProR:::gp_standardize_prediction_table(pred)
  out <- PredictProR:::gp_public_prediction_table(standardized)

  expect_identical(out$Predicted_class, c("Resistant", "Susceptible"))
})

test_that("Variance component contract preserves Bayesian style rows", {
  vc <- data.frame(
    Components = c(4, 1, 0.8),
    Standard_error = c(0.4, 0.1, 0.05),
    row.names = c("genetic_variance", "residual_variance", "heritability")
  )

  out <- PredictProR:::gp_format_variance_components(list(Variance_components = vc))

  expect_identical(names(out), PredictProR:::gp_variance_component_columns())
  expect_identical(rownames(out), c("genetic_variance", "residual_variance", "heritability"))
  expect_identical(out$Component, c("genetic_variance", "residual_variance", "heritability"))
})

test_that("Variance component contract carries explicit component labels", {
  vc <- data.frame(
    Components = c(4, 1, 0.8),
    Standard_error = c(0.4, 0.1, 0.05),
    row.names = c("genetic_variance", "residual_variance", "heritability")
  )

  out <- PredictProR:::gp_format_variance_components(list(Variance_components = vc))

  expect_true(all(c("Component", "Components", "Standard_error") %in% names(out)))
  expect_identical(out$Component, c("genetic_variance", "residual_variance", "heritability"))
})

test_that("GP single-environment variance components prefer REML SE and drop duplicate reporting rows", {
  raw <- list(
    predicted_values = data.frame(
      GID = c("g1", "g2", "g3"),
      Env = c("ENV1", "ENV1", "ENV1"),
      Predicted_value = c(1.1, 1.2, 1.3),
      stringsAsFactors = FALSE
    ),
    gp_result = list(
      var_components = data.frame(
        term = c("G", "G", "G:Env", "G:Env", "Env"),
        kernel = c("gmatrix", "total", "gmatrix", "total", "total"),
        component_type = c("genetic_main", "genetic_main", "interaction", "interaction", "environment_main"),
        estimate = c(4.0, 4.0, 0, 0, 0),
        stringsAsFactors = FALSE
      ),
      varcomp = data.frame(
        component = c("w_g[0]", "sigma2_noise"),
        estimate = c(4.0, 1.0),
        std.error = c(0.4, 0.2),
        bound = c("P", "P"),
        stringsAsFactors = FALSE
      )
    ),
    gp_info = list(gp_varcomp_mode = "reml")
  )

  out <- PredictProR:::gp_public_variance_components(raw)

  expect_identical(names(out), PredictProR:::gp_variance_component_columns())
  expect_identical(rownames(out), c("genetic_variance", "residual_variance", "heritability"))
  expect_equal(out["genetic_variance", "Standard_error"], 0.4)
  expect_equal(out["residual_variance", "Standard_error"], 0.2)
  expect_false(any(grepl("G:Env|Env|G\\.1", rownames(out))))
  expect_true(is.finite(out["heritability", "Standard_error"]))
  expect_false("Estimation_method" %in% names(out))
})

test_that("GP single-environment REML boundary residual is not reported as zero heritability input", {
  raw <- list(
    predicted_values = data.frame(
      GID = c("g1", "g2", "g3"),
      Env = c("ENV1", "ENV1", "ENV1"),
      Predicted_value = c(1.1, 1.2, 1.3),
      stringsAsFactors = FALSE
    ),
    gp_result = list(
      var_components = data.frame(
        term = c("G", "G", "Env"),
        kernel = c("gmatrix", "total", "total"),
        component_type = c("genetic_main", "genetic_main", "environment_main"),
        estimate = c(4.0, 4.0, 0),
        stringsAsFactors = FALSE
      ),
      varcomp = data.frame(
        component = c("vm(GID,G)!var", "Env!var", "Env_ENV1!R"),
        estimate = c(4.0, 0.0, 0.0),
        std.error = c(0.4, NA_real_, NA_real_),
        bound = c("P", "B", "B"),
        stringsAsFactors = FALSE
      )
    ),
    gp_info = list(gp_varcomp_mode = "reml")
  )

  out <- PredictProR:::gp_public_variance_components(raw)

  expect_identical(rownames(out), c("genetic_variance", "residual_variance", "heritability"))
  # Residual (Env_ENV1!R, bound="B") at zero: genetic and residual variance are
  # not separable, so genetic, residual and h^2 are all not estimable.
  expect_true(all(is.na(out[c("genetic_variance", "residual_variance", "heritability"), "Components"])))
  expect_identical(unname(out[, "Boundary"]), rep("not_estimable", 3L))
  expect_false("Estimation_method" %in% names(out))
})

test_that("GP REML variance components do not mix ASReml varcomp with report-bundle summaries", {
  raw <- list(
    predicted_values = data.frame(
      GID = c("g1", "g2", "g3"),
      Predicted_value = c(1.1, 1.2, 1.3),
      stringsAsFactors = FALSE
    ),
    gp_result = list(
      var_components = data.frame(
        term = "G",
        kernel = "total",
        component_type = "genetic_main",
        estimate = 4.0,
        stringsAsFactors = FALSE
      ),
      varcomp = data.frame(
        component = "sigma2_noise",
        estimate = 1.0,
        std.error = 0.2,
        bound = "P",
        stringsAsFactors = FALSE
      )
    ),
    gp_info = list(gp_varcomp_mode = "reml")
  )

  out <- PredictProR:::gp_public_variance_components(raw)

  expect_identical(out$Component, "reml_variance_components_incomplete")
  expect_true(is.na(out$Components[[1]]))
  expect_false("Estimation_method" %in% names(out))
})

test_that("GP MET REML variance components use backend varcomp instead of prediction fallback", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), times = 2),
    Env = rep(c("E1", "E2"), each = 2),
    Predicted_value = c(1, 2, 10, 20),
    Train_Test_Label = "Train",
    Observed_value = c(1.1, 1.9, 9.8, 20.2),
    stringsAsFactors = FALSE
  )
  raw <- list(
    predicted_values = pred,
    gp_result = list(
      varcomp = data.frame(
        component = c("vm(GID,G)!var", "Env_E1!R", "Env_E2!R"),
        estimate = c(2.5, 0.4, 0.8),
        std.error = c(0.3, 0.1, 0.2),
        bound = "P",
        stringsAsFactors = FALSE
      )
    ),
    gp_info = list(gp_varcomp_mode = "reml")
  )

  out <- PredictProR:::gp_public_variance_components(raw)

  expect_true(all(c("vm(GID,G)!var", "Env_E1!R", "Env_E2!R") %in% out$Component))
  expect_false("Estimation_method" %in% names(out))
  expect_false(any(out$Component == "genetic_variance"))
})

test_that("GP MET REML boundary residual rows are not reported as zero variance", {
  pred <- data.frame(
    GID = rep(c("g1", "g2"), times = 2),
    Env = rep(c("E1", "E2"), each = 2),
    Predicted_value = c(1, 2, 10, 20),
    Train_Test_Label = "Train",
    Observed_value = c(1.1, 1.9, 9.8, 20.2),
    stringsAsFactors = FALSE
  )
  raw <- list(
    predicted_values = pred,
    gp_result = list(
      varcomp = data.frame(
        component = c("vm(GID,G)!var", "Env_E1!R", "Env_E2!R"),
        estimate = c(2.5, 0.0, 0.8),
        std.error = c(0.3, NA_real_, 0.2),
        bound = c("P", "B", "P"),
        stringsAsFactors = FALSE
      )
    ),
    gp_info = list(gp_varcomp_mode = "reml")
  )

  out <- PredictProR:::gp_public_variance_components(raw)

  e1 <- out[out$Component == "Env_E1!R", , drop = FALSE]
  e2 <- out[out$Component == "Env_E2!R", , drop = FALSE]
  expect_true(is.na(e1$Components[[1]]))
  expect_equal(e2$Components[[1]], 0.8)
  expect_false("Estimation_method" %in% names(out))
})

test_that("GP bridge defaults request REML variance components", {
  params <- PredictProR:::gp_bridge_model_params("KRR")

  expect_identical(params$varcomp_mode, "reml")
})

test_that("GP variance components are not synthesized from prediction tables", {
  pred <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Predicted_value = c(1.1, 1.8, 3.1, 3.9),
    Train_Test_Label = c("Train", "Train", "Test", "Test"),
    Observed_value = c(1.0, 2.0, NA, NA),
    Standard_error = c(0.1, 0.2, 0.3, 0.4),
    PEV = c(0.01, 0.04, 0.09, 0.16),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_variance_components(list(
    predicted_values = pred,
    gp_result = list(predictions = pred),
    gp_info = list(backend = "mock_gp")
  ))

  expect_true(all(c("Component", "Components", "Standard_error") %in% names(out)))
  expect_identical(out$Component, "variance_components_not_estimated")
  expect_true(is.na(out$Components[[1]]))
  expect_false(any(out$Component %in% c("genetic_variance", "residual_variance", "heritability")))
})

test_that("ML Gaussian bootstrap output separates predictive and genetic estimands", {
  boot <- rbind(
    c(1.0, 2.0, 3.0, 4.0),
    c(1.1, 2.1, 2.9, 4.1),
    c(0.9, 1.9, 3.1, 3.9)
  )

  out <- PredictProR:::gp_ml_gaussian_variance_components(
    boot_matrix = boot,
    train_test_label = c("Train", "Train", "Test", "Test"),
    observed_y = c(1.2, 1.8),
    predicted_value = colMeans(boot),
    prediction_error_var = apply(boot, 2, var)
  )

  expect_identical(names(out), PredictProR:::gp_variance_component_columns())
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable",
    "phenotypic_reference_variance",
    "prediction_variance",
    "estimated_prediction_mse"
  ) %in% out$Component))
  expect_true(all(is.na(out$Components[grepl("not_identifiable$", out$Component)])))
  expect_true(all(is.finite(out$Components[!grepl("not_identifiable$", out$Component)])))
  expect_true(any(is.finite(out$Standard_error)))
  expect_false("Estimation_method" %in% names(out))
})

test_that("Public ML predictive fallback does not synthesize mixed-model components", {
  pred <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Trait = c("YLD", "YLD", "HT", "HT"),
    Predicted_value = c(1.1, 1.8, 3.1, 3.9),
    Train_Test_Label = c("Train", "Test", "Train", "Test"),
    Observed_value = c(1.0, NA, 3.0, NA),
    Standard_error = c(0.1, 0.2, 0.1, 0.2),
    PEV = c(0.01, 0.04, 0.01, 0.04),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_variance_components(
    list(predicted_values = pred),
    model = "RandomForest",
    response_family = "gaussian"
  )

  expect_true(all(c("Trait", "Component", "Components", "Standard_error") %in% names(out)))
  expect_false("Env" %in% names(out))
  expect_true(setequal(unique(out$Trait), c("YLD", "HT")))
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable",
    "prediction_variance",
    "estimated_prediction_mse"
  ) %in% out$Component))
})

test_that("Single-environment trait variance components omit Env", {
  pred <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Trait = c("YLD", "YLD", "HT", "HT"),
    Predicted_value = c(1.1, 1.8, 3.1, 3.9),
    Train_Test_Label = c("Train", "Test", "Train", "Test"),
    Observed_value = c(1.0, NA, 3.0, NA),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_prediction_table_variance_components(pred)

  expect_identical(
    names(out),
    PredictProR:::gp_variance_component_columns(include_trait = TRUE)
  )
  expect_true("Trait" %in% names(out))
  expect_false("Env" %in% names(out))
  expect_true(setequal(unique(out$Trait), c("YLD", "HT")))
})

test_that("Public MET predictive fallback reports environment-specific prediction variance", {
  pred <- data.frame(
    GID = rep(paste0("g", 1:4), times = 2),
    Env = rep(c("E1", "E2"), each = 4),
    Predicted_value = c(1.0, 1.4, 1.8, 2.2, 3.0, 4.0, 5.0, 6.0),
    Train_Test_Label = "Train",
    Observed_value = c(1.1, 1.3, 1.7, 2.3, 2.8, 4.1, 4.9, 6.2),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_prediction_table_variance_components(pred)

  expect_true(all(c("Env", "Component", "Components", "Standard_error") %in% names(out)))
  expect_true(setequal(unique(out$Env), c("E1", "E2")))
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "prediction_variance"
  ) %in% out$Component))
  e1_prediction <- out$Components[out$Env == "E1" & out$Component == "prediction_variance"]
  e2_prediction <- out$Components[out$Env == "E2" & out$Component == "prediction_variance"]
  expect_equal(e1_prediction, stats::var(pred$Predicted_value[pred$Env == "E1"]))
  expect_equal(e2_prediction, stats::var(pred$Predicted_value[pred$Env == "E2"]))
})

test_that("Hybrid prediction contract uses common hybrid columns only", {
  pred <- data.frame(
    HybridID = c("F1_M1", "F2_M2"),
    Female = c("F1", "F2"),
    Male = c("M1", "M2"),
    Env = c("E1", "E1"),
    Predicted_value = c(5.1, 5.4),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(5.0, NA),
    Standard_error = c(0.2, 0.3),
    Prediction_error_variance = c(0.04, 0.09),
    Reliability_reference_variance = c(1, 1),
    Female_GCA = c(0.1, 0.2),
    SCA_effect = c(0.05, -0.02),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_hybrid_prediction_table(pred, hybrid_id = "HybridID")

  expect_identical(
    names(out),
    PredictProR:::gp_hybrid_prediction_columns()
  )
  expect_false("Female_GCA" %in% names(out))
  expect_false("SCA_effect" %in% names(out))
  expect_true(all(is.finite(out$PEV)))
  expect_true(all(is.finite(out$lower_bound)))
  expect_true(all(is.finite(out$upper_bound)))
  expect_true(all(is.finite(out$Reliability)))
  expect_true(all(nzchar(out$Reliability_remarks)))
})

test_that("Hybrid standardizer keeps Reliability unavailable with PEV and retains stability", {
  pred <- data.frame(
    HybridID = paste0("H", 1:6),
    Female = c("F1", "F1", "F2", "F2", "F3", "F3"),
    Male = c("M1", "M2", "M1", "M2", "M1", "M2"),
    Predicted_value = c(4.9, 5.2, 5.8, 6.1, 5.5, 6.4),
    Train_Test_Label = c("Train", "Train", "Train", "Train", "Test", "Test"),
    Observed_value = c(5.0, 5.3, 5.5, 6.2, NA, NA),
    Female_additive_contribution = c(-0.2, -0.2, 0.4, 0.4, 0.1, 0.1),
    Male_additive_contribution = c(-0.1, 0.2, -0.1, 0.2, -0.1, 0.2),
    Hybrid_interaction_contribution = c(0.0, 0.1, 0.2, 0.3, -0.1, 0.4),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_standardize_hybrid_model_result(
    list(predicted_values = pred),
    model = "RandomForest"
  )

  expect_identical(
    names(out$predicted_values),
    PredictProR:::gp_hybrid_prediction_columns()
  )
  expect_true(all(is.na(out$predicted_values$Standard_error)))
  expect_true(all(is.na(out$predicted_values$PEV)))
  expect_true(all(is.na(out$predicted_values$lower_bound)))
  expect_true(all(is.na(out$predicted_values$upper_bound)))
  expect_true(all(is.na(out$predicted_values$Reliability)))
  expect_true(all(is.finite(out$predicted_values$Prediction_stability)))
  expect_true(all(grepl("not genetic reliability", out$predicted_values$Reliability_basis, fixed = TRUE)))
  expect_true(all(is.na(out$predicted_values$Reliability_remarks)))
  expect_false("Estimation_method" %in% names(out$variance_components))
})

test_that("Hybrid variance component contract uses a hybrid-specific schema", {
  vc <- data.frame(
    component = c("female_gca_variance", "male_gca_variance", "sca_variance"),
    estimate = c(1.2, 0.8, 0.4),
    std.error = c(0.1, 0.08, 0.05),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_hybrid_variance_components(vc)

  expect_identical(
    names(out),
    PredictProR:::gp_hybrid_variance_component_columns()
  )
  expect_false("Env" %in% names(out))
  expect_identical(out$Hybrid_component, c("female_gca", "male_gca", "sca"))
  expect_identical(out$Component, vc$component)
  expect_equal(out$Components, vc$estimate)
  expect_false("Estimation_method" %in% names(out))
})

test_that("Hybrid variance component fallback labels empirical decomposition components", {
  pred <- data.frame(
    HybridID = paste0("H", 1:5),
    Female = paste0("F", c(1, 1, 2, 2, 3)),
    Male = paste0("M", c(1, 2, 1, 2, 1)),
    Predicted_value = c(4.9, 5.2, 5.8, 6.1, 5.5),
    Train_Test_Label = c("Train", "Train", "Train", "Train", "Test"),
    Observed_value = c(5.0, 5.3, 5.5, 6.2, NA),
    Female_additive_contribution = c(-0.2, -0.2, 0.4, 0.4, 0.1),
    Male_additive_contribution = c(-0.1, 0.2, -0.1, 0.2, -0.1),
    Hybrid_interaction_contribution = c(0.0, 0.1, 0.2, 0.3, -0.1),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_hybrid_variance_components(
    list(predicted_values = pred),
    model = "RandomForest"
  )

  expect_identical(
    names(out),
    PredictProR:::gp_hybrid_variance_component_columns()
  )
  expect_false("Env" %in% names(out))
  expect_true(all(out$Hybrid_component %in% c("female_gca", "male_gca", "sca", "hybrid_other")))
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable"
  ) %in% out$Component))
  expect_true(any(is.finite(out$Components)))
  expect_false("Estimation_method" %in% names(out))
})

test_that("Hybrid variance component fallback aligns multi-environment rows", {
  pred <- data.frame(
    HybridID = rep(paste0("H", 1:4), times = 2L),
    Female = rep(c("F1", "F1", "F2", "F2"), times = 2L),
    Male = rep(c("M1", "M2", "M1", "M2"), times = 2L),
    Env = rep(c("E1", "E2"), each = 4L),
    Predicted_value = c(4.9, 5.2, 5.8, 6.1, 5.4, 5.7, 6.3, 6.6),
    Train_Test_Label = "Train",
    Observed_value = c(5.0, 5.3, 5.5, 6.2, 5.5, 5.8, 6.1, 6.7),
    Female_additive_contribution = c(-0.2, -0.2, 0.4, 0.4, -0.1, -0.1, 0.5, 0.5),
    Male_additive_contribution = c(-0.1, 0.2, -0.1, 0.2, 0.0, 0.3, 0.0, 0.3),
    Hybrid_interaction_contribution = c(0.0, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.4),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_hybrid_variance_components(
    list(predicted_values = pred),
    model = "RandomForest"
  )

  expect_identical(
    names(out),
    PredictProR:::gp_hybrid_variance_component_columns(include_env = TRUE)
  )
  expect_true(all(c("E1", "E2") %in% out$Env))
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable"
  ) %in% out$Component))
})

test_that("ASReml-style numeric component column is treated as variance estimate", {
  vc <- data.frame(
    component = c(1.2, 0.4),
    std.error = c(0.2, 0.08),
    bound = c("P", "P"),
    stringsAsFactors = FALSE
  )
  rownames(vc) <- c("vm(Female, female_inv)", "units!R")

  out <- PredictProR:::gp_format_variance_components(vc, default_estimation_method = "ASReml_REML")

  expect_identical(out$Component, rownames(vc))
  expect_equal(out$Components, vc$component)
  expect_false("Estimation_method" %in% names(out))
})

test_that("Diagnostic plot contract has the same slots across models", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  out <- PredictProR:::gp_public_diagnostic_plots(list(diagnostic_plots = p))

  expect_identical(
    names(out),
    c("predicted_vs_observed", "predicted_vs_reliability", "prediction_interval", "reliability_summary")
  )
  expect_true(inherits(out$predicted_vs_observed, "ggplot"))
})

test_that("Diagnostic plot contract collapses true-prediction GP plot lists to the standard plot slot", {
  p1 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  p2 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  out <- PredictProR:::gp_public_diagnostic_plots(list(
    diagnostic_plots = list(
      predicted_vs_reliability = p1,
      prediction_inter_vs_reliability = p2
    )
  ))

  expect_identical(
    names(out),
    c("predicted_vs_observed", "predicted_vs_reliability", "prediction_interval", "reliability_summary")
  )
  expect_true(inherits(out$predicted_vs_observed, "gtable"))
  expect_true(all(vapply(out[-1], is.null, logical(1L))))
})

test_that("Single-environment auto formatting drops constant Env and Trait columns", {
  pred <- data.frame(
    GID = c("g1", "g2"),
    Env = c("E1", "E1"),
    Trait = c("YLD", "YLD"),
    Predicted_value = c(1.2, 2.4),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(1.1, NA),
    Standard_error = c(0.1, 0.2),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_format_gaussian_prediction_table(pred)

  expect_identical(names(out), PredictProR:::gp_gaussian_prediction_columns())
  expect_false("Env" %in% names(out))
  expect_false("Trait" %in% names(out))
})
