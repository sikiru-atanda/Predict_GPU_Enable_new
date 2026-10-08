test_that("multi-kernel eigenfeature builder cbinds aligned blocks", {
  ids <- paste0("g", 1:4)
  K1 <- diag(4)
  dimnames(K1) <- list(ids, ids)
  K2 <- matrix(c(
    1, 0.2, 0.1, 0,
    0.2, 1, 0.3, 0.1,
    0.1, 0.3, 1, 0.2,
    0, 0.1, 0.2, 1
  ), 4, 4, byrow = TRUE)
  dimnames(K2) <- list(ids, ids)

  res <- gp_prepare_met_kernel_features(
    gmatrix = K1,
    omic1_kernel = K2,
    omics_kernel_label = list(omic1_kernel = "COP"),
    var_explained = 0.95
  )

  expect_s3_class(res$feature_table, "data.frame")
  expect_equal(rownames(res$feature_table), ids)
  expect_true(any(grepl("^GRM_PC", names(res$feature_table))))
  expect_true(any(grepl("^COP_PC", names(res$feature_table))))
  expect_equal(sort(res$feature_summary$source), c("GRM", "omic1"))
})

test_that("public MET results retain the kernel feature summary", {
  feature_summary <- data.frame(
    source = c("GRM", "Gaussian", "Matern32"),
    n_pcs = c(4L, 3L, 2L),
    stringsAsFactors = FALSE
  )
  out <- gp_standardize_public_model_result(
    list(
      model_parameters = data.frame(
        stat = "model_type", summary = "test", stringsAsFactors = FALSE
      ),
      predicted_values = data.frame(
        GID = c("g1", "g2"),
        Predicted_value = c(0.1, 0.2),
        Train_Test_Label = c("Train", "Test"),
        Observed_value = c(0.15, NA_real_),
        Standard_error = c(0.2, 0.3),
        PEV = c(0.04, 0.09),
        stringsAsFactors = FALSE
      ),
      met_feature_summary = feature_summary
    ),
    gen_name = "GID"
  )

  expect_identical(out$met_feature_summary, feature_summary)
})

test_that("MET kernel eigenfeature builder uses user kernel_list blocks", {
  ids <- paste0("g", 1:4)
  K1 <- diag(4)
  dimnames(K1) <- list(ids, ids)
  K2 <- matrix(c(
    1, 0.4, 0.2, 0.1,
    0.4, 1, 0.3, 0.2,
    0.2, 0.3, 1, 0.5,
    0.1, 0.2, 0.5, 1
  ), 4, 4, byrow = TRUE)
  dimnames(K2) <- list(ids, ids)

  res <- gp_prepare_met_kernel_features(
    kernel_list = list(VanRaden = K1, Yang = K2),
    var_explained = 0.95
  )

  expect_s3_class(res$feature_table, "data.frame")
  expect_equal(rownames(res$feature_table), ids)
  expect_equal(sort(res$feature_summary$source), c("VanRaden", "Yang"))
  expect_true(any(grepl("^VanRaden_PC", names(res$feature_table))))
  expect_true(any(grepl("^Yang_PC", names(res$feature_table))))
})

test_that("validate_met_input_standard accepts kernel_list as the relationship source", {
  ids <- paste0("g", 1:4)
  ph <- expand.grid(
    GID = ids,
    Env = c("E1", "E2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ph$Yield <- seq_len(nrow(ph))
  K <- diag(length(ids))
  dimnames(K) <- list(ids, ids)

  expect_error(
    validate_met_input_standard(
      pheno_data = ph,
      response = "Yield",
      gen_name = "GID",
      heter_groups = "Env",
      met_ml_dl = TRUE,
      GS_model = "RandomForest",
      kernel_list = list(user_kernel = K)
    ),
    NA
  )
})

test_that("MET long data builder expands genotype features over environments", {
  ids <- paste0("g", 1:3)
  feature_table <- data.frame(
    GRM_PC1 = c(0.1, 0.2, 0.3),
    COP_PC1 = c(1, 2, 3),
    row.names = ids,
    check.names = FALSE
  )
  ph <- data.frame(
    GID = c("g1", "g1", "g2", "g3"),
    Env = c("E1", "E2", "E1", "E2"),
    Yield = c(5.1, NA, 4.8, 5.5),
    stringsAsFactors = FALSE
  )

  long_dat <- gp_build_met_long_data(
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    kernel_feature_table = feature_table
  )

  expect_equal(nrow(long_dat), 4L)
  expect_true(all(c("GRM_PC1", "COP_PC1", ".met_row_id") %in% names(long_dat)))
  expect_equal(long_dat$GRM_PC1[long_dat$GID == "g1" & long_dat$Env == "E1"], 0.1)
  expect_equal(long_dat$GRM_PC1[long_dat$GID == "g1" & long_dat$Env == "E2"], 0.1)
})

test_that("MET CV assignments support CV0 CV1 and CV2", {
  ph <- data.frame(
    GID = c("g1", "g1", "g2", "g2", "g3", "g3"),
    Env = c("E1", "E2", "E1", "E2", "E1", "E2"),
    Yield = c(1, 2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )

  cv0 <- gp_met_cv_assignments(ph, "GID", "Env", "CV0", nfolds = 2, random_state = 1, replication = 1)
  cv1 <- gp_met_cv_assignments(ph, "GID", "Env", "CV1", nfolds = 2, random_state = 1, replication = 1)
  cv2 <- gp_met_cv_assignments(ph, "GID", "Env", "CV2", nfolds = 2, random_state = 1, replication = 1)

  expect_length(cv0, 1L)
  expect_length(cv1, 1L)
  expect_length(cv2, 1L)
  expect_equal(length(cv0[[1]]), nrow(ph))
  expect_equal(length(cv1[[1]]), nrow(ph))
  expect_equal(length(cv2[[1]]), nrow(ph))
})

test_that("Gaussian MET ML/DL summaries expose plot-ready uncertainty and across-environment prediction", {
  ph <- data.frame(
    GID = rep(paste0("g", 1:5), each = 3),
    Env = rep(c("E1", "E2", "E3"), times = 5),
    Yield = c(
      4.9, 5.2, NA,
      5.8, NA, 6.1,
      NA, 4.7, 5.0,
      6.2, 6.4, NA,
      5.3, NA, 5.5
    ),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    GID = ph$GID,
    Env = ph$Env,
    Predicted_value = c(
      5.0, 5.1, 5.4,
      5.6, 5.9, 6.0,
      4.8, 4.9, 5.1,
      6.0, 6.3, 6.5,
      5.4, 5.6, 5.7
    ),
    Train_Test_Label = ifelse(is.na(ph$Yield), "Test", "Train"),
    stringsAsFactors = FALSE
  )
  gaussian_cols <- c(
    "Standard_error",
    "PEV",
    "lower_bound",
    "upper_bound",
    "Uncertainty",
    "Uncertainty_remarks",
    "Reliability",
    "Reliability_remarks"
  )

  summary_res <- gp_met_summary_statistics(
    predicted_values = pred,
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "RandomForest",
    response_family = "gaussian"
  )

  expect_true(all(gaussian_cols %in% names(summary_res$MET_prediction_table)))
  expect_true(all(is.na(summary_res$MET_prediction_table$Standard_error)))
  expect_true(all(is.na(summary_res$MET_prediction_table$PEV)))
  expect_true(all(is.na(summary_res$MET_prediction_table$Reliability)))
  expect_true(all(is.finite(summary_res$MET_prediction_table$Prediction_stability)))
  expect_true(all(
    summary_res$MET_prediction_table$Prediction_uncertainty_source ==
      "unavailable_no_heldout_or_model_based_uncertainty"
  ))

  total_pred <- summary_res$across_environment_predicted_values
  expect_s3_class(total_pred, "data.frame")
  expect_equal(nrow(total_pred), length(unique(ph$GID)))
  expect_true(all(gaussian_cols %in% names(total_pred)))
  expect_true(all(is.na(total_pred$Standard_error)))
  expect_true(all(is.na(total_pred$PEV)))
  expect_true(all(is.na(total_pred$Reliability)))
  expect_true(any(summary_res$summary_statistics$summary == "Per-environment and across-environment prediction"))
  expect_true(is.matrix(summary_res$Prediction_covariance_environments))
  expect_true(is.matrix(summary_res$Prediction_error_covariance_environments))
  expect_true(is.matrix(summary_res$Prediction_covariance_environments_SE))
  expect_false("Genetic_covariance_environments" %in% names(summary_res))

  plots <- gp_met_true_prediction_plot(
    predicted_values = summary_res$MET_prediction_table,
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    model_label = "RandomForest",
    response_family = "gaussian"
  )
  expect_type(plots, "list")
  expect_true(all(c(
    "predicted_vs_observed",
    "predicted_vs_reliability",
    "prediction_interval",
    "reliability_summary"
  ) %in% names(plots)))
  expect_s3_class(plots$predicted_vs_reliability, "ggplot")
  expect_true(is.null(plots$prediction_interval))
})

test_that("Gaussian MET ML/DL in-sample predictions do not synthesize individual uncertainty", {
  gids <- paste0("g", 1:6)
  envs <- c("E1", "E2", "E3")
  ph <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  gid_eff <- setNames(c(-0.7, -0.25, 0.05, 0.35, 0.7, 1.15), gids)
  env_eff <- setNames(c(-0.35, 0.1, 0.45), envs)
  ph$Yield <- 5 + gid_eff[ph$GID] + env_eff[ph$Env]
  pred <- data.frame(
    GID = ph$GID,
    Env = ph$Env,
    Predicted_value = 5 + 0.8 * gid_eff[ph$GID] + 0.65 * env_eff[ph$Env] +
      rep(c(-0.42, 0.08, 0.27, -0.12, 0.38, -0.05), each = length(envs)),
    Train_Test_Label = "Train",
    stringsAsFactors = FALSE
  )

  summary_res <- gp_met_summary_statistics(
    predicted_values = pred,
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "RandomForest",
    response_family = "gaussian"
  )

  by_env <- split(summary_res$MET_prediction_table, summary_res$MET_prediction_table$Env)
  expect_true(all(vapply(by_env, function(df) all(is.na(df$Standard_error)), logical(1))))
  expect_true(all(vapply(by_env, function(df) all(is.na(df$Reliability)), logical(1))))
  expect_true(all(vapply(by_env, function(df) all(is.finite(df$Prediction_stability)), logical(1))))

  total_pred <- summary_res$across_environment_predicted_values
  expect_true(all(is.na(total_pred$Standard_error)))
  expect_true(all(is.na(total_pred$Reliability)))
})

test_that("Gaussian MET ML/DL keeps reliability unavailable with PEV while retaining stability", {
  gids <- paste0("g", 1:6)
  envs <- c("E1", "E2")
  pred_obs <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pred_obs$Predicted_value <- rep(c(8.0, 6.8, 5.9, 4.4, 3.5, 2.1), times = length(envs)) +
    rep(c(0, 0.3), each = length(gids))
  pred_obs$Observed_value <- rep(c(1.00, 1.05, 0.95, 1.10, 0.90, 1.02), times = length(envs))
  pred_obs$Train_Test_Label <- "Train"

  out <- gp_met_add_gaussian_uncertainty(
    pred_obs = pred_obs,
    gen_name = "GID",
    heter_groups = "Env",
    ml_dl_estimand = TRUE
  )

  by_env <- split(out, out$Env)
  expect_true(all(vapply(by_env, function(df) all(is.na(df$PEV)), logical(1))))
  expect_true(all(vapply(by_env, function(df) all(is.na(df$Reliability)), logical(1))))
  expect_true(all(vapply(by_env, function(df) all(is.na(df$Reliability_variance_input)), logical(1))))
  expect_true(all(vapply(by_env, function(df) all(is.finite(df$Prediction_stability)), logical(1))))

  total_pred <- gp_met_across_environment_prediction(
    pred_obs = out,
    gen_name = "GID",
    heter_groups = "Env"
  )
  expect_true(all(is.na(total_pred$Reliability)))
})

test_that("Gaussian MET ML/DL uses fitted prediction variance on the marker-adjustment scale", {
  ph <- data.frame(
    GID = rep(paste0("g", 1:5), times = 2),
    Env = rep(c("E1", "E2"), each = 5),
    Yield = c(1.0, 1.4, 1.8, 2.2, 2.6, 4.0, 5.0, 6.0, 7.0, 8.0),
    stringsAsFactors = FALSE
  )
  pred <- data.frame(
    GID = ph$GID,
    Env = ph$Env,
    Predicted_value = c(1.1, 1.5, 1.7, 2.3, 2.7, 3.8, 5.2, 5.8, 7.3, 7.9),
    Train_Test_Label = "Train",
    stringsAsFactors = FALSE
  )

  summary_res <- gp_met_summary_statistics(
    predicted_values = pred,
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model = "DenseNeuralNet",
    response_family = "gaussian"
  )

  env_pred <- summary_res$MET_prediction_table
  expect_true("Prediction_stability_reference_variance" %in% names(env_pred))
  expect_false("Genetic_variance" %in% names(env_pred))
  by_env <- split(env_pred, env_pred$Env)
  expected_ref <- tapply(pred$Predicted_value, pred$Env, stats::var)
  expect_true(all(vapply(names(by_env), function(env) {
    all.equal(unique(round(by_env[[env]]$Prediction_stability_reference_variance, 10)), round(expected_ref[[env]], 10)) == TRUE
  }, logical(1))))
  expect_true(all(is.na(env_pred$Reliability_reference_variance)))
  expect_true(all(is.na(env_pred$Reliability)))
  expect_true(all(is.finite(env_pred$Prediction_stability)))

  across <- summary_res$across_environment_predicted_values
  expect_true("Prediction_stability_reference_variance" %in% names(across))
  expect_false("Genetic_variance" %in% names(across))
  expect_equal(
    unique(round(across$Prediction_stability_reference_variance, 10)),
    round(mean(expected_ref), 10)
  )
  expect_true(all(is.na(across$Reliability)))
})

test_that("Gaussian MET ML/DL computes non-genetic reliability from heldout PEV", {
  pred_obs <- data.frame(
    GID = rep(paste0("g", 1:4), times = 2),
    Env = rep(c("E3", "E4"), each = 4),
    Predicted_value = c(1.0, 1.4, 1.8, 2.2, 3.0, 3.4, 3.8, 4.2),
    Observed_value = c(1.1, 1.2, 2.0, 2.1, NA, NA, NA, NA),
    Train_Test_Label = rep(c("Train", "Test"), each = 4),
    Standard_error = c(0.2, 0.3, 0.4, 0.5, 1.0, 1.2, 1.4, 1.6),
    PEV = c(0.04, 0.09, 0.16, 0.25, 1.0, 1.44, 1.96, 2.56),
    Reliability = rep(0, 8),
    PEV_basis = "heldout cross-fitted predictive MSE; not mixed-model PEV",
    stringsAsFactors = FALSE
  )

  out <- gp_met_add_gaussian_uncertainty(
    pred_obs = pred_obs,
    gen_name = "GID",
    heter_groups = "Env",
    ml_dl_estimand = TRUE
  )

  e4 <- out[out$Env == "E4", , drop = FALSE]
  expect_gt(length(unique(round(e4$PEV, 10))), 1L)
  expect_equal(
    e4$Reliability,
    e4$Reliability_reference_variance /
      (e4$Reliability_reference_variance + e4$PEV)
  )
  expect_equal(e4$Reliability_variance_input, e4$PEV)
  expect_gt(length(unique(round(e4$Prediction_stability, 10))), 1L)
  expect_true(all(e4$Prediction_stability > 0 & e4$Prediction_stability < 1))
})

test_that("Gaussian MET ML/DL does not derive reliability from resampling variance alone", {
  pred_obs <- data.frame(
    GID = paste0("g", 1:6),
    Env = "E1",
    Predicted_value = c(1, 2, 3, 10, 11, 12),
    Observed_value = c(1.1, 1.9, 3.2, NA, NA, NA),
    Train_Test_Label = c("Train", "Train", "Train", "Test", "Test", "Test"),
    Repeated_prediction_variance = rep(0.01, 6),
    stringsAsFactors = FALSE
  )

  out <- gp_met_add_gaussian_uncertainty(
    pred_obs = pred_obs,
    gen_name = "GID",
    heter_groups = "Env",
    ml_dl_estimand = TRUE
  )

  expect_identical(out$Train_Test_Label, pred_obs$Train_Test_Label)
  expect_true(all(is.na(out$PEV)))
  expect_true(all(is.na(out$Reliability_reference_variance)))
  expect_true(all(is.na(out$Reliability)))
  expect_true(all(is.finite(out$Prediction_stability)))
  expect_true(all(grepl(
    "Reliability unavailable",
    out$Reliability_basis,
    fixed = TRUE
  )))
})

test_that("MET CatBoost CV helper returns numeric predictions for held-out rows", {
  py <- skip_if_no_ml_python()
  Sys.unsetenv("RETICULATE_PYTHON")
  Sys.setenv(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli")

  gids <- paste0("g", 1:6)
  envs <- c("E1", "E2", "E3")
  ph <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  set.seed(1)
  M <- matrix(rnorm(length(gids) * 6), nrow = length(gids), dimnames = list(gids, paste0("m", 1:6)))
  K <- tcrossprod(scale(M, center = TRUE, scale = FALSE))
  diag(K) <- diag(K) + 1e-3
  dimnames(K) <- list(gids, gids)
  g_eff <- setNames(rnorm(length(gids), 0, 0.2), gids)
  e_eff <- setNames(c(-0.1, 0.2, 0.4), envs)
  ph$Yield <- 5 + g_eff[ph$GID] + e_eff[ph$Env]
  tst <- which(ph$GID %in% c("g2", "g5") & ph$Env %in% c("E2", "E3"))

  pred <- gp_met_catboost_cv_predict(
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    gmatrix = K,
    tst = tst,
    model_params = list(
      catboost_iterations = 20,
      catboost_depth = 4,
      catboost_learning_rate = 0.1,
      catboost_l2_leaf_reg = 3,
      catboost_thread_count = 1L
    )
  )

  expect_type(pred, "double")
  expect_equal(length(pred), length(tst))
})

test_that("MET uncertainty calibration holds complete GIDs out across environments", {
  gids <- rep(paste0("g", 1:6), each = 3)
  envs <- rep(c("E1", "E2", "E3"), times = 6)
  y <- seq_len(length(gids)) / 3
  fitted_group_sets <- list()

  calibration <- gp_met_grouped_crossfit_calibration(
    y = y,
    observed_idx = seq_along(y),
    group_id = gids,
    environment = envs,
    n_targets = length(y),
    predictor_matrix = cbind(
      gid_score = rep(seq_len(6), each = 3),
      env_score = rep(seq_len(3), times = 6)
    ),
    row_id = paste(gids, envs, sep = "__"),
    fit_predict = function(fit_idx, prediction_idx, fold_index) {
      fitted_group_sets[[fold_index]] <<- unique(gids[fit_idx])
      mean(y[fit_idx]) + 0.07 * prediction_idx + 0.02 * fold_index
    },
    nfolds = 3L,
    seed = 917L
  )

  expect_false(is.null(calibration))
  expect_equal(calibration$heldout_n, length(y))
  expect_equal(dim(calibration$target_fold_predictions), c(3L, length(y)))
  expect_equal(length(fitted_group_sets), 3L)
  expect_true(all(vapply(fitted_group_sets, function(fit_groups) {
    length(setdiff(unique(gids), fit_groups)) == 2L
  }, logical(1L))))
})

test_that("MET bootstrap uncertainty is finite, target-specific, and predictive", {
  gids <- rep(paste0("g", 1:6), each = 3)
  envs <- rep(c("E1", "E2", "E3"), times = 6)
  y <- seq_len(length(gids)) / 3
  predicted <- y + sin(seq_along(y)) / 8
  calibration <- gp_met_grouped_crossfit_calibration(
    y = y,
    observed_idx = seq_along(y),
    group_id = gids,
    environment = envs,
    n_targets = length(y),
    predictor_matrix = cbind(
      gid_score = rep(seq_len(6), each = 3),
      env_score = rep(seq_len(3), times = 6)
    ),
    row_id = paste(gids, envs, sep = "__"),
    fit_predict = function(fit_idx, prediction_idx, fold_index) {
      mean(y[fit_idx]) + 0.07 * prediction_idx +
        0.015 * fold_index * prediction_idx
    },
    nfolds = 3L,
    seed = 917L
  )
  offsets <- seq(-1, 1, length.out = 10)
  bootstrap_matrix <- vapply(offsets, function(offset) {
    predicted + offset * seq(0.03, 0.3, length.out = length(predicted))
  }, numeric(length(predicted)))
  bootstrap_matrix <- t(bootstrap_matrix)

  uncertainty <- gp_met_bootstrap_uncertainty_columns(
    predicted_value = predicted,
    bootstrap_matrix = bootstrap_matrix,
    observed_y = y,
    train_test_label = rep("Train", length(y)),
    heldout_calibration = calibration,
    n_bootstrap = nrow(bootstrap_matrix)
  )

  expect_true(all(is.finite(uncertainty$Standard_error)))
  expect_true(all(is.finite(uncertainty$PEV)))
  expect_equal(uncertainty$PEV, uncertainty$Standard_error^2)
  expect_gt(length(unique(round(uncertainty$Standard_error, 10))), 1L)
  expect_true(all(grepl("not mixed-model PEV", uncertainty$PEV_basis, fixed = TRUE)))
  expect_true(all(uncertainty$Prediction_uncertainty_source ==
    "local_cross_fitted_residual_risk_plus_fold_instability"))
  expect_false(any(grepl("unavailable", uncertainty$Prediction_uncertainty_source, fixed = TRUE)))
})
