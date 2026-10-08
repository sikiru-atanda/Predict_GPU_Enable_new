direct_marker_target_summary <- function(mod, ETA, model_name, bayes_para) {
  eta_types <- vapply(ETA$ETA, function(x) x$model, character(1))
  random_idx <- which(eta_types != "FIXED")

  total_draws <- NULL
  prior_genetic_variance <- rep(0, length(mod$model$yHat))
  tst <- which(is.na(mod$model$y))
  train_idx <- which(!is.na(mod$model$y))
  draw_indices <- gp_bayes_saved_draw_indices(
    bayes_para[["nIter"]],
    bayes_para[["burnIn"]],
    bayes_para[["thin"]]
  )
  for (eta_idx in random_idx) {
    b_file <- mod$output_files_names[grepl(
      paste0("ETA_", eta_idx, "_b\\.bin$"),
      basename(mod$output_files_names)
    )]
    beta <- as.matrix(BGLR::readBinMat(b_file[1]))
    if (max(draw_indices) <= nrow(beta)) {
      beta <- beta[draw_indices, , drop = FALSE]
    }
    X_eta <- as.matrix(ETA$ETA[[eta_idx]]$X)
    eta_draws <- as.matrix(X_eta %*% t(beta))
    total_draws <- if (is.null(total_draws)) eta_draws else total_draws + eta_draws
    var_u_res <- process_var_u_new(
      GS_model = model_name,
      geno_data = X_eta,
      y = mod$model$y,
      B = beta,
      tst = tst
    )
    component_mean <- mean(var_u_res$var_u_omics, na.rm = TRUE)
    row_energy <- rowSums(X_eta ^ 2)
    reference_energy <- mean(row_energy[train_idx], na.rm = TRUE)
    k_diag <- row_energy / reference_energy
    prior_genetic_variance <- prior_genetic_variance + (component_mean * k_diag)
  }

  mu_file <- mod$output_files_names[grepl("mu\\.dat$", basename(mod$output_files_names))]
  mu_draws <- as.double(utils::read.table(mu_file, header = FALSE)[[1]][draw_indices])
  pred_draws <- as.matrix(sweep(total_draws, 2, mu_draws, FUN = "+"))
  pred_mean <- as.double(mod$model$yHat)
  standard_error <- apply(pred_draws, 1, stats::sd, na.rm = TRUE)
  prediction_error_var <- apply(pred_draws, 1, stats::var, na.rm = TRUE)
  fallback_genetic_variance <- mean(
    apply(total_draws[train_idx, , drop = FALSE], 2, stats::var, na.rm = TRUE),
    na.rm = TRUE
  )
  target_genetic_variance <- prior_genetic_variance
  bad_idx <- !is.finite(target_genetic_variance) | target_genetic_variance <= 0
  target_genetic_variance[bad_idx] <- fallback_genetic_variance
  reliability_res <- reliability_thresholds(
    prediction_error_var = prediction_error_var,
    genetic_var = target_genetic_variance
  )

  list(
    predicted_mean = pred_mean,
    standard_error = standard_error,
    prediction_error_var = prediction_error_var,
    target_genetic_variance = target_genetic_variance,
    reliability = reliability_res$reliability
  )
}

expect_bayes_environment_covariance_output <- function(out, labels) {
  expect_true(is.matrix(out$Genetic_covariance_environments))
  expect_true(is.matrix(out$Genetic_correlation_environments))
  expect_true(is.matrix(out$Genetic_covariance_environments_SE))
  expect_true(is.matrix(out$Genetic_correlation_environments_SE))
  expect_identical(rownames(out$Genetic_covariance_environments), labels)
  expect_equal(
    out$Genetic_covariance_environments,
    t(out$Genetic_covariance_environments),
    tolerance = 1e-10
  )
  expect_true(min(eigen(
    out$Genetic_covariance_environments,
    symmetric = TRUE,
    only.values = TRUE
  )$values) > -1e-8)
  expect_equal(unname(diag(out$Genetic_correlation_environments)), rep(1, length(labels)))
  expect_true(all(is.finite(out$Genetic_covariance_environments_SE)))
  expect_true(all(is.finite(out$Genetic_correlation_environments_SE)))
}

test_that("single-environment BayesA shared-fit metrics match direct BGLR summaries", {
  skip_if_not_installed("BGLR")

  set.seed(101)
  pheno <- data.frame(
    GID = paste0("g", 1:8),
    Yield = rnorm(8),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(2, 6)] <- NA_real_

  geno <- matrix(rnorm(8 * 4), nrow = 8, ncol = 4)
  rownames(geno) <- pheno$GID
  colnames(geno) <- paste0("m", 1:4)

  prep <- bayes_finalize_A_B_C_BL_BRR(
    random = ~ GID,
    GS_model = "BayesA",
    response = "Yield",
    pheno_data = pheno,
    geno_data = geno,
    gen_name = "GID",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "BayesA"
  )

  direct <- direct_marker_target_summary(
    mod = mod,
    ETA = prep$bayes_ETA,
    model_name = "BayesA",
    bayes_para = prep$bayes_para
  )
  out <- mod_output_bayes(
    mod = mod,
    ETA = prep$bayes_ETA,
    geno_data = geno,
    gen_name = "GID",
    bayes_para = prep$bayes_para,
    GS_model = "BayesA",
    confidence_level = 0.90
  )

  pred <- out$Predicted_value
  expect_equal(as.double(pred$Predicted_value), direct$predicted_mean, tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Standard_error)), unname(direct$standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$PEV)), unname(direct$prediction_error_var), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Reliability)), unname(direct$reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Genetic_variance)), unname(direct$target_genetic_variance), tolerance = 1e-8)
  expect_true(all(pred$PEV_basis == "BGLR posterior target-draw variance"))
  expect_true(all(pred$Prediction_uncertainty_source == "BGLR posterior target draws"))
  expect_true(all(pred$Prediction_interval_method == "BGLR posterior target-draw quantiles"))
  expect_true(all(pred$Prediction_interval_nominal_coverage == 0.90))
  expect_true(all(pred$Prediction_interval_calibration_n == ncol(direct$pred_draws)))
  expect_gt(mean(pred$Reliability[pred$Train_Test_Label == "Train"], na.rm = TRUE), 0)
})

test_that("single-environment BayesA plus omic source matches direct BGLR summaries", {
  skip_if_not_installed("BGLR")

  set.seed(111)
  pheno <- data.frame(
    GID = paste0("g", 1:8),
    Yield = rnorm(8),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(3, 7)] <- NA_real_

  geno <- matrix(rnorm(8 * 4), nrow = 8, ncol = 4)
  omic1 <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3)
  rownames(geno) <- rownames(omic1) <- pheno$GID
  colnames(geno) <- paste0("m", 1:4)
  colnames(omic1) <- paste0("o", 1:3)

  prep <- bayes_finalize_A_B_C_BL_BRR(
    random = ~ GID,
    GS_model = "BayesA",
    response = "Yield",
    pheno_data = pheno,
    geno_data = geno,
    omic1_data = omic1,
    gen_name = "GID",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "BayesA"
  )

  direct <- direct_marker_target_summary(
    mod = mod,
    ETA = prep$bayes_ETA,
    model_name = "BayesA",
    bayes_para = prep$bayes_para
  )
  out <- mod_output_bayes(
    mod = mod,
    ETA = prep$bayes_ETA,
    geno_data = geno,
    omic1_data = omic1,
    gen_name = "GID",
    bayes_para = prep$bayes_para,
    GS_model = "BayesA"
  )

  pred <- out$Predicted_value
  expect_equal(as.double(pred$Predicted_value), direct$predicted_mean, tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Standard_error)), unname(direct$standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$PEV)), unname(direct$prediction_error_var), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Reliability)), unname(direct$reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Genetic_variance)), unname(direct$target_genetic_variance), tolerance = 1e-8)
  expect_gt(mean(pred$Reliability[pred$Train_Test_Label == "Train"], na.rm = TRUE), 0)
  expect_identical(sort(names(out$Coefficients)), c("Geno_coefficient", "Omics_coefficient"))
  expect_identical(sort(names(out$Estimated_breeding_value)), c("Geno_estimated_breeding_value", "Omics_estimated_breeding_value"))
})

test_that("single-environment Bayes Alphabet models retain every multi-omics block", {
  skip_if_not_installed("BGLR")

  set.seed(112)
  pheno <- data.frame(
    GID = paste0("g", 1:10),
    Yield = rnorm(10),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(3, 8)] <- NA_real_
  geno <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  omic1 <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  omic2 <- matrix(rnorm(10 * 2), nrow = 10, ncol = 2)
  rownames(geno) <- rownames(omic1) <- rownames(omic2) <- pheno$GID
  colnames(geno) <- paste0("m", 1:5)
  colnames(omic1) <- paste0("transcript", 1:3)
  colnames(omic2) <- paste0("metabolite", 1:2)

  for (model in c("BRR", "BayesA", "BayesB", "BayesC", "BL")) {
    set.seed(700 + match(model, c("BRR", "BayesA", "BayesB", "BayesC", "BL")))
    fit <- bayes_finalize_A_B_C_BL_BRR(
      random = ~ GID,
      GS_model = model,
      response = "Yield",
      pheno_data = pheno,
      geno_data = geno,
      omic1_data = omic1,
      omic2_data = omic2,
      gen_name = "GID",
      nIter = 120,
      burnIn = 20,
      thin = 2,
      cross_validation = FALSE,
      system_database = TRUE
    )
    out <- fit$bayes_result
    expect_length(out$M_matrix_model_ready, 3L)
    expect_length(out$Coefficients, 3L)
    expect_length(out$Estimated_breeding_value, 3L)
    expect_true(all(is.finite(out$Predicted_value$Predicted_value)), info = model)
    expect_true(all(c(
      "genetic_variance_geno_data",
      "genetic_variance_omic1_data",
      "genetic_variance_omic2_data",
      "total_genetic_variance",
      "residual_variance",
      "heritability"
    ) %in% rownames(out$Variance_components)), info = model)
    parameter <- stats::setNames(out$model_parameters$summary, out$model_parameters$stat)
    expect_identical(parameter[["bayesian_feature_block_count"]], "3", info = model)
    expect_identical(
      parameter[["bayesian_total_genetic_variance_estimand"]],
      "posterior_variance_of_summed_genetic_effects",
      info = model
    )
  }
})

test_that("single-environment RKHS shared-fit metrics match direct target summaries", {
  skip_if_not_installed("BGLR")

  set.seed(102)
  pheno <- data.frame(
    GID = paste0("g", 1:6),
    Yield = rnorm(6),
    stringsAsFactors = FALSE
  )
  pheno$Yield[c(1, 5)] <- NA_real_

  gmatrix <- diag(6)
  rownames(gmatrix) <- colnames(gmatrix) <- pheno$GID

  prep <- bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "RKHS"
  )

  direct <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "RKHS",
    pheno_data = prep$bayes_ETA$pheno_data,
    gen_name = "GID",
    output_files_names = mod$output_files_names,
    draw_indices = gp_bayes_saved_draw_indices(
      prep$bayes_para[["nIter"]], prep$bayes_para[["burnIn"]], prep$bayes_para[["thin"]]
    )
  )

  out <- mod_output_bayes_gbluBRR_RKHS(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "RKHS",
    gmatrix = gmatrix,
    gen_name = "GID",
    pheno_data = prep$bayes_ETA$pheno_data,
    bayes_para = prep$bayes_para,
    confidence_level = 0.90
  )
  pred <- out$Predicted_value
  expect_equal(as.double(pred$Predicted_value), direct$observation$predicted_mean, tolerance = 1e-8)
  expect_equal(as.double(pred$Standard_error), direct$observation$Standard_error, tolerance = 1e-8)
  expect_equal(as.double(pred$PEV), direct$observation$PEV, tolerance = 1e-8)
  expect_equal(as.double(pred$Reliability), direct$observation$Reliability, tolerance = 1e-8)
  expect_equal(as.double(pred$Genetic_variance), direct$observation$Genetic_variance, tolerance = 1e-8)
  expect_true(all(pred$PEV_basis == "BGLR RKHS analytic posterior prediction variance"))
  expect_true(all(pred$Prediction_uncertainty_source ==
                    "BGLR RKHS analytic posterior prediction covariance"))
  expect_true(all(pred$Prediction_interval_method ==
                    "Gaussian approximation from BGLR RKHS posterior covariance"))
  expect_true(all(pred$Prediction_interval_nominal_coverage == 0.90))
  expect_true(all(is.na(pred$Prediction_interval_calibration_n)))
})

test_that("multi-environment BRR shared-fit totals match direct target summaries and labels", {
  skip_if_not_installed("BGLR")

  set.seed(103)
  pheno <- expand.grid(
    GID = paste0("g", 1:5),
    Env = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- rnorm(nrow(pheno))
  test_gid <- c("g2", "g4")
  pheno$Yield[pheno$GID %in% test_gid] <- NA_real_

  gmatrix <- diag(5)
  rownames(gmatrix) <- colnames(gmatrix) <- unique(pheno$GID)

  prep <- bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "BRR",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "BRR"
  )

  direct <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "BRR",
    pheno_data = prep$bayes_ETA$pheno_data,
    gen_name = "GID",
    heter_groups = "Env",
    output_files_names = mod$output_files_names,
    draw_indices = gp_bayes_saved_draw_indices(
      prep$bayes_para[["nIter"]], prep$bayes_para[["burnIn"]], prep$bayes_para[["thin"]]
    )
  )

  out <- mod_output_bayes_gbluBRR_RKHS(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "BRR",
    gmatrix = gmatrix,
    gen_name = "GID",
    pheno_data = prep$bayes_ETA$pheno_data,
    heter_groups = "Env",
    bayes_para = prep$bayes_para,
    confidence_level = 0.90
  )
  expect_bayes_environment_covariance_output(out, c("E1", "E2"))

  pred_total <- out$Total_Predicted_value
  expect_equal(unname(as.double(pred_total$Predicted_value)), unname(direct$across$predicted_mean), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Standard_error)), unname(direct$across$Standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$PEV)), unname(direct$across$PEV), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Reliability)), unname(direct$across$Reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Genetic_variance)), unname(direct$across$Genetic_variance), tolerance = 1e-8)
  expect_identical(
    as.character(pred_total$Train_Test_Label),
    ifelse(as.character(pred_total$GID) %in% test_gid, "Test", "Train")
  )
})

test_that("multi-environment BRR row-level observation metrics match direct target summaries", {
  skip_if_not_installed("BGLR")

  set.seed(104)
  pheno <- expand.grid(
    GID = paste0("g", 1:5),
    Env = c("E1", "E2", "E3"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- rnorm(nrow(pheno))
  test_gid <- c("g2", "g5")
  pheno$Yield[pheno$GID %in% test_gid] <- NA_real_

  gmatrix <- diag(5)
  rownames(gmatrix) <- colnames(gmatrix) <- unique(pheno$GID)

  prep <- bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "BRR",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "BRR"
  )

  direct <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "BRR",
    pheno_data = prep$bayes_ETA$pheno_data,
    gen_name = "GID",
    heter_groups = "Env",
    output_files_names = mod$output_files_names,
    draw_indices = gp_bayes_saved_draw_indices(
      prep$bayes_para[["nIter"]], prep$bayes_para[["burnIn"]], prep$bayes_para[["thin"]]
    )
  )

  out <- mod_output_bayes_gbluBRR_RKHS(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "BRR",
    gmatrix = gmatrix,
    gen_name = "GID",
    pheno_data = prep$bayes_ETA$pheno_data,
    heter_groups = "Env",
    bayes_para = prep$bayes_para,
    confidence_level = 0.90
  )

  pred <- out$Predicted_value
  expect_equal(unname(as.double(pred$Predicted_value)), unname(direct$observation$predicted_mean), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Standard_error)), unname(direct$observation$Standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$PEV)), unname(direct$observation$PEV), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Reliability)), unname(direct$observation$Reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Genetic_variance)), unname(direct$observation$Genetic_variance), tolerance = 1e-8)
  expect_true(all(pred$PEV_basis == "BGLR posterior target-draw variance"))
  expect_true(all(pred$Prediction_uncertainty_source == "BGLR posterior target draws"))
  expect_true(all(pred$Prediction_interval_method == "BGLR posterior target-draw quantiles"))
  expect_true(all(pred$Prediction_interval_nominal_coverage == 0.90))
  expect_true(all(pred$Prediction_interval_calibration_n ==
                    length(gp_bayes_saved_draw_indices(
                      prep$bayes_para[["nIter"]], prep$bayes_para[["burnIn"]],
                      prep$bayes_para[["thin"]]
                    ))))
  finite_interval <- is.finite(pred$lower_bound) &
    is.finite(pred$Predicted_value) & is.finite(pred$upper_bound)
  expect_true(all(pred$lower_bound[finite_interval] <= pred$Predicted_value[finite_interval]))
  expect_true(all(pred$upper_bound[finite_interval] >= pred$Predicted_value[finite_interval]))
  expect_identical(as.character(pred$Train_Test_Label), ifelse(as.character(pred$GID) %in% test_gid, "Test", "Train"))
})

test_that("multi-environment RKHS shared-fit row-level and across metrics match direct target summaries", {
  skip_if_not_installed("BGLR")

  set.seed(105)
  pheno <- expand.grid(
    GID = paste0("g", 1:5),
    Env = c("E1", "E2", "E3"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- rnorm(nrow(pheno))
  test_gid <- c("g1", "g4")
  pheno$Yield[pheno$GID %in% test_gid] <- NA_real_

  gmatrix <- diag(5)
  rownames(gmatrix) <- colnames(gmatrix) <- unique(pheno$GID)

  prep <- bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  mod <- bayes_mod_execute(
    pheno_data = prep$bayes_ETA$pheno_data,
    response = "Yield",
    ETA = prep$bayes_ETA$ETA,
    bayes_para = prep$bayes_para,
    GS_model = "RKHS"
  )

  direct <- bayes_compute_prediction_targets(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "RKHS",
    pheno_data = prep$bayes_ETA$pheno_data,
    gen_name = "GID",
    heter_groups = "Env",
    output_files_names = mod$output_files_names,
    draw_indices = gp_bayes_saved_draw_indices(
      prep$bayes_para[["nIter"]], prep$bayes_para[["burnIn"]], prep$bayes_para[["thin"]]
    )
  )

  out <- mod_output_bayes_gbluBRR_RKHS(
    mod = mod,
    ETA = prep$bayes_ETA,
    GS_model = "RKHS",
    gmatrix = gmatrix,
    gen_name = "GID",
    pheno_data = prep$bayes_ETA$pheno_data,
    heter_groups = "Env",
    bayes_para = prep$bayes_para
  )
  expect_bayes_environment_covariance_output(out, c("E1", "E2", "E3"))

  pred <- out$Predicted_value
  expect_equal(unname(as.double(pred$Predicted_value)), unname(direct$observation$predicted_mean), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Standard_error)), unname(direct$observation$Standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$PEV)), unname(direct$observation$PEV), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Reliability)), unname(direct$observation$Reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred$Genetic_variance)), unname(direct$observation$Genetic_variance), tolerance = 1e-8)
  expect_identical(as.character(pred$Train_Test_Label), ifelse(as.character(pred$GID) %in% test_gid, "Test", "Train"))
  expect_gt(mean(pred$Reliability, na.rm = TRUE), 0.3)
  expect_gt(mean(pred$Reliability[pred$Train_Test_Label == "Train"], na.rm = TRUE), 0.3)
  expect_gt(mean(pred$Reliability[pred$Train_Test_Label == "Test"], na.rm = TRUE), 0.3)

  pred_total <- out$Total_Predicted_value
  expect_equal(unname(as.double(pred_total$Predicted_value)), unname(direct$across$predicted_mean), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Standard_error)), unname(direct$across$Standard_error), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$PEV)), unname(direct$across$PEV), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Reliability)), unname(direct$across$Reliability), tolerance = 1e-8)
  expect_equal(unname(as.double(pred_total$Genetic_variance)), unname(direct$across$Genetic_variance), tolerance = 1e-8)
  expect_gt(mean(pred_total$Reliability, na.rm = TRUE), 0.4)
  expect_gt(mean(pred_total$Reliability[pred_total$Train_Test_Label == "Train"], na.rm = TRUE), 0.4)
  expect_gt(mean(pred_total$Reliability[pred_total$Train_Test_Label == "Test"], na.rm = TRUE), 0.4)
  expect_identical(
    as.character(pred_total$Train_Test_Label),
    ifelse(as.character(pred_total$GID) %in% test_gid, "Test", "Train")
  )
})

test_that("multi-environment RKHS defaults to BRR main effect and RKHS interaction", {
  set.seed(106)
  pheno <- expand.grid(
    GID = paste0("g", 1:4),
    Env = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  pheno$Yield <- rnorm(nrow(pheno))

  gmatrix <- diag(4)
  rownames(gmatrix) <- colnames(gmatrix) <- unique(pheno$GID)

  prep <- bayes_finalize_RKHS_GBLUPBRR(
    random = ~ GID + GID:Env,
    GS_model = "RKHS",
    response = "Yield",
    pheno_data = pheno,
    gmatrix = gmatrix,
    gen_name = "GID",
    heter_groups = "Env",
    nIter = 120,
    burnIn = 20,
    thin = 2,
    cross_validation = TRUE
  )

  expect_identical(
    vapply(prep$bayes_ETA$ETA, function(x) x$model, character(1)),
    c("BRR", "RKHS")
  )
})

test_that("Bayesian summary statistics include target-specific metric definitions", {
  mod <- list(model = list(y = c(1, 2, NA), yHat = c(1.1, 2.1, 2.9), varE = 0.5, ETA = list()))
  model_result <- list(
    Predicted_value = data.frame(
      GID = c("g1", "g2", "g3"),
      Predicted_value = c(1.1, 2.1, 2.9),
      PEV = c(0.1, 0.2, 0.3),
      Reliability = c(0.9, 0.8, 0.7),
      Genetic_variance = c(1, 1, 1)
    ),
    Total_Predicted_value = data.frame(
      GID = c("g1", "g2", "g3"),
      Predicted_value = c(1, 2, 3),
      PEV = c(0.1, 0.1, 0.1),
      Reliability = c(0.9, 0.9, 0.9),
      Genetic_variance = c(1, 1, 1)
    )
  )

  res <- summary_statistics_bayes(
    mod = mod,
    GS_model = "RKHS",
    model_result = model_result,
    gen_name = "GID",
    heter_groups = "Env"
  )

  stats_df <- res$summary_statistics
  expect_true(all(c(
    "Prediction_Target",
    "Genetic_Variance_Definition",
    "PEV_Definition",
    "Reliability_Definition",
    "Across_Environment_Target"
  ) %in% stats_df$stat))
})
