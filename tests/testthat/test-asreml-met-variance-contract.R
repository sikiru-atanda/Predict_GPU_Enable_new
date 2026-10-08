make_asreml_met_corh_summary <- function(rho = 0.25) {
  envs <- c("E1", "E2", "E3")
  rows <- c(
    "corh(Env):vm(GID, gmatrix_inv)!Env!cor",
    paste0("corh(Env):vm(GID, gmatrix_inv)!Env_", envs),
    paste0("Env_", envs, "!R")
  )
  component <- c(rho, 1.2, 0.8, 2.0, 0.4, 0.5, 0.6)
  structure(
    list(
      varcomp = data.frame(
        component = component,
        std.error = c(0.03, 0.12, 0.08, 0.20, 0.04, 0.05, 0.06),
        z.ratio = NA_real_,
        bound = "P",
        `%ch` = 0.1,
        row.names = rows,
        check.names = FALSE
      ),
      mf = data.frame(
        Env = factor(rep(envs, each = 2L), levels = envs)
      )
    ),
    class = "summary.asreml"
  )
}

make_asreml_met_weighted_summary <- function(rho = 0.25) {
  envs <- c("E1", "E2", "E3")
  rows <- c(
    "corh(Env):vm(GID, gmatrix_inv)!Env!cor",
    paste0("corh(Env):vm(GID, gmatrix_inv)!Env_", envs),
    paste0("Env_", envs, "!R")
  )
  structure(
    list(
      varcomp = data.frame(
        component = c(rho, 1.2, 0.8, 2.0, 1, 1, 1),
        std.error = c(0.03, 0.12, 0.08, 0.20, NA, NA, NA),
        z.ratio = NA_real_,
        bound = c(rep("P", 4L), rep("F", 3L)),
        `%ch` = 0.1,
        row.names = rows,
        check.names = FALSE
      ),
      mf = data.frame(
        GID = paste0("g", seq_len(9L)),
        Env = factor(rep(envs, each = 3L), levels = envs),
        Yield = rep(c(1, 2, NA_real_), 3L),
        .PredictProR_stage2_precision_weight = c(
          1, 2, 1,
          2, 4, 1,
          0.5, 1, 1
        ),
        check.names = FALSE
      )
    ),
    class = "summary.asreml"
  )
}

test_that("ASReml MET reliability requires one positive genetic variance per environment", {
  env <- c("E1", "E2", "E3")
  expect_true(PredictProR:::asreml_met_reliability_variance_valid(
    c(E1 = 1.2, E2 = 0.8, E3 = 0.5), env
  ))
  expect_false(PredictProR:::asreml_met_reliability_variance_valid(NULL, env))
  expect_false(PredictProR:::asreml_met_reliability_variance_valid(c(1, 0.8), env))
  expect_false(PredictProR:::asreml_met_reliability_variance_valid(c(1, NA, 0.5), env))
  expect_false(PredictProR:::asreml_met_reliability_variance_valid(c(1, 0, 0.5), env))

  ebv <- data.frame(
    Env = env,
    Prediction_error_variance = c(0.2, 0.3, 0.1),
    stringsAsFactors = FALSE
  )
  applied <- PredictProR:::asreml_apply_met_reliability(
    ebv_df = ebv,
    heter_groups = "Env",
    heter_grp = env,
    variance_by_env = c(1.2, 0.8, 0.5)
  )
  expect_equal(applied$Reliability,
               1 - c(0.2 / 1.2, 0.3 / 0.8, 0.1 / 0.5))
  expect_equal(applied$Reliability_variance_input,
               ebv$Prediction_error_variance)
  expect_equal(applied$Reliability_reference_variance, c(1.2, 0.8, 0.5))
  expect_true(all(grepl(
    "model-specific genetic variance",
    applied$Reliability_basis,
    fixed = TRUE
  )))
})

test_that("ASReml public MET reliability replaces the phenotypic fallback", {
  pred <- data.frame(
    GID = rep(c("G1", "G2"), each = 3L),
    Env = rep(c("E1", "E2", "E3"), 2L),
    Predicted_value = seq_len(6L) / 10,
    Train_Test_Label = rep(c("Train", "Test"), each = 3L),
    Observed_value = c(1:3, NA, NA, NA),
    Standard_error = sqrt(c(0.2, 0.3, 0.1, 0.4, 0.2, 0.3)),
    PEV = c(0.2, 0.3, 0.1, 0.4, 0.2, 0.3),
    Reliability = rep(0.75, 6L),
    Reliability_variance_input = rep(NA_real_, 6L),
    Reliability_reference_variance = rep(2, 6L),
    Reliability_basis = rep("observed_response_variance", 6L),
    stringsAsFactors = FALSE
  )
  applied <- PredictProR:::asreml_apply_met_reliability(
    ebv_df = pred,
    heter_groups = "Env",
    heter_grp = c("E1", "E2", "E3"),
    variance_by_env = c(1.2, 0.8, 0.5)
  )
  expected_reference <- rep(c(1.2, 0.8, 0.5), 2L)
  expected <- pmax(0, pmin(1, 1 - pred$PEV / expected_reference))
  expect_equal(applied$Reliability_variance_input, pred$PEV)
  expect_equal(applied$Reliability_reference_variance, expected_reference)
  expect_equal(applied$Reliability, expected)
  expect_true(all(grepl(
    "model-specific genetic variance",
    applied$Reliability_basis,
    fixed = TRUE
  )))

  public <- PredictProR:::gp_format_gaussian_prediction_table(
    applied, gen_name = "GID", heter_groups = "Env"
  )
  expect_equal(public$Reliability_variance_input, pred$PEV)
  expect_equal(public$Reliability_reference_variance, expected_reference)
  expect_equal(public$Reliability, expected)
  expect_true(all(grepl(
    "model-specific genetic variance",
    public$Reliability_basis,
    fixed = TRUE
  )))
})

test_that("ASReml MET exposes fitted genetic and residual matrices", {
  fit <- make_asreml_met_corh_summary()
  extracted <- asreml_herit_varCov_new(
    model = fit,
    heter_groups = "Env",
    var_cov_str = "corh",
    heter_resid = TRUE,
    names_in_inv_list = "gmatrix_inv",
    inter_gen_pos = 2L,
    gen_pos = 1L
  )

  expected_vg <- c(E1 = 1.2, E2 = 0.8, E3 = 2.0)
  expect_equal(diag(extracted$Genetic_covariance_environments), expected_vg)
  expect_equal(unname(diag(extracted$Genetic_correlation_environments)), rep(1, 3L))
  expect_equal(
    extracted$Genetic_covariance_environments[1L, 2L],
    0.25 * sqrt(1.2 * 0.8),
    tolerance = 1e-12
  )
  expect_equal(unname(diag(extracted$Residual_covariance_environments)), c(0.4, 0.5, 0.6))
  expect_equal(unname(extracted$Residual_correlation_environments), diag(3L))
  expect_equal(unname(extracted$varG_per_omics_SE[1L, ]), c(0.12, 0.08, 0.20))
  expect_equal(extracted$Residual_Var_SE, c(E1 = 0.04, E2 = 0.05, E3 = 0.06))
  expect_true(all(is.finite(extracted$Heritability)))

  formatted <- asreml_variance_components(extracted)
  expect_true(is.list(formatted))
  vc <- formatted$Variance_components
  expect_equal(vc["genetic_variance_E1", "Components"], 1.2)
  expect_equal(vc["genetic_variance_E1", "Standard_error"], 0.12)
  expect_equal(vc["residual_variance_E1", "Components"], 0.4)
  expect_equal(vc["residual_variance_E1", "Standard_error"], 0.04)
  expect_true(all(vc$Estimation_method == "ASReml_REML_corh"))
})

test_that("weighted ASReml MET reports response-scale residual variance by environment", {
  extracted <- asreml_herit_varCov_new(
    model = make_asreml_met_weighted_summary(),
    heter_groups = "Env",
    var_cov_str = "corh",
    heter_resid = TRUE,
    names_in_inv_list = "gmatrix_inv",
    inter_gen_pos = 2L,
    gen_pos = 1L,
    response = "Yield",
    gen_name = "GID"
  )

  expected_ve <- c(E1 = 0.75, E2 = 0.375, E3 = 1.5)
  expect_equal(extracted$Residual_scale_parameter, c(E1 = 1, E2 = 1, E3 = 1))
  expect_equal(extracted$Residual_Var, expected_ve)
  expect_true(all(is.na(extracted$Residual_Var_SE)))
  expect_equal(unname(diag(extracted$Residual_covariance_environments)), unname(expected_ve))
  expect_equal(
    as.numeric(extracted$Heritability),
    c(1.2 / (1.2 + 0.75), 0.8 / (0.8 + 0.375), 2 / (2 + 1.5))
  )

  detail <- extracted$Residual_variance_by_observation
  expect_equal(detail$Observation_index, c(1L, 2L, 4L, 5L, 7L, 8L))
  expect_equal(detail$Residual_variance,
               c(1, 0.5, 0.5, 0.25, 2, 1))
  expect_equal(detail$Residual_variance,
               detail$Residual_scale_parameter / detail$Precision_weight)
  expect_false(any(detail$GID %in% c("g3", "g6", "g9")))

  env_summary <- extracted$Residual_variance_by_environment
  expect_equal(env_summary$Env, names(expected_ve))
  expect_equal(env_summary$Residual_variance, unname(expected_ve))
  expect_equal(env_summary$N_observations, rep(2L, 3L))

  formatted <- asreml_variance_components(extracted)$Variance_components
  expect_equal(formatted["residual_variance_E1", "Components"], expected_ve[["E1"]])
  expect_true(is.na(formatted["residual_variance_E1", "Standard_error"]))
  expect_identical(
    rownames(formatted)[1:3],
    c(
      "mean_genetic_variance_across_environments",
      "mean_residual_variance_across_environments",
      "heritability_from_mean_variances"
    )
  )
})

test_that("ASReml MET canonical matrices survive the public result contract", {
  extracted <- asreml_herit_varCov_new(
    model = make_asreml_met_weighted_summary(),
    heter_groups = "Env",
    var_cov_str = "corh",
    heter_resid = TRUE,
    names_in_inv_list = "gmatrix_inv",
    inter_gen_pos = 2L,
    gen_pos = 1L,
    response = "Yield",
    gen_name = "GID"
  )
  formatted <- asreml_variance_components(extracted)
  pred <- data.frame(
    GID = c("g1", "g2"),
    Env = c("E1", "E2"),
    Predicted_value = c(0.2, 0.4),
    Train_Test_Label = c("Train", "Test"),
    Observed_value = c(0.1, NA_real_),
    Standard_error = c(0.2, 0.3)
  )

  public <- PredictProR:::gp_standardize_public_model_result(
    list(
      Predicted_value = pred,
      Variance_components = formatted$Variance_components,
      Covariance = formatted$Covariance,
      Correlation = formatted$Correlation,
      Genetic_covariance_environments = extracted$Genetic_covariance_environments,
      Genetic_correlation_environments = extracted$Genetic_correlation_environments,
      Residual_covariance_environments = extracted$Residual_covariance_environments,
      Residual_correlation_environments = extracted$Residual_correlation_environments,
      Residual_variance_by_environment = extracted$Residual_variance_by_environment,
      Residual_variance_by_observation = extracted$Residual_variance_by_observation
    ),
    gen_name = "GID",
    heter_groups = "Env"
  )

  expect_equal(
    public$Genetic_covariance_environments,
    extracted$Genetic_covariance_environments
  )
  expect_equal(
    public$Genetic_correlation_environments,
    extracted$Genetic_correlation_environments
  )
  expect_true(nrow(public$covariance_components) > 0L)
  expect_true(nrow(public$correlation_components) > 0L)
  expect_true(any(public$variance_components$Component ==
                    "mean_genetic_variance_across_environments"))
  expect_equal(public$Residual_variance_by_environment,
               extracted$Residual_variance_by_environment)
  expect_equal(public$Residual_variance_by_observation,
               extracted$Residual_variance_by_observation)
  expect_false(any(grepl("not_identifiable", public$variance_components$Component)))
})

test_that("invalid ASReml MET correlation estimates are rejected", {
  envs <- c("E1", "E2", "E3")
  rows <- c(
    paste0(
      "corgh(Env):vm(GID, gmatrix_inv)!Env!",
      c("E2:!Env!E1", "E3:!Env!E1", "E3:!Env!E2"),
      ".cor"
    ),
    paste0("corgh(Env):vm(GID, gmatrix_inv)!Env_", envs),
    paste0("Env_", envs, "!R")
  )
  bad <- structure(
    list(
      varcomp = data.frame(
        component = c(0.9, 0.9, -0.9, 1, 1, 1, 0.5, 0.5, 0.5),
        std.error = 0.1,
        bound = "P",
        row.names = rows
      ),
      mf = data.frame(Env = factor(rep(envs, each = 2L), levels = envs))
    ),
    class = "summary.asreml"
  )

  expect_error(
    asreml_herit_varCov_new(
      model = bad,
      heter_groups = "Env",
      var_cov_str = "corgh",
      heter_resid = TRUE,
      names_in_inv_list = "gmatrix_inv",
      inter_gen_pos = 2L,
      gen_pos = 1L
    ),
    "Reconstructed correlation matrix is invalid"
  )
})

test_that("ASReml output is withheld for untrustworthy fits", {
  stable <- list(
    converge = TRUE,
    ifault = 0L,
    vparameters.pc = c(vg = 0.4, rho = 8),
    vparameters.con = c("P", "B")
  )
  expect_no_error(PredictProR:::gp_asreml_assert_trustworthy_fit(stable))

  not_converged <- stable
  not_converged$converge <- FALSE
  expect_error(
    PredictProR:::gp_asreml_assert_trustworthy_fit(not_converged),
    "did not converge"
  )

  changing <- stable
  changing$vparameters.pc <- c(vg = 1.2, rho = 0)
  expect_error(
    PredictProR:::gp_asreml_assert_trustworthy_fit(changing),
    "changing by more than 1%"
  )

  fault <- stable
  fault$ifault <- 2L
  expect_error(
    PredictProR:::gp_asreml_assert_trustworthy_fit(fault),
    "ifault = 2"
  )
})

test_that("ASReml refinement continues when converge is true but components are unstable", {
  mod <- list(
    converge = TRUE,
    ifault = 0L,
    vparameters.pc = c(vg = 2.4, rho = 0),
    vparameters.con = c("P", "B")
  )
  update_count <- 0L
  update_fn <- function(x) {
    update_count <<- update_count + 1L
    x$vparameters.pc[["vg"]] <- x$vparameters.pc[["vg"]] / 2
    x
  }

  refined <- PredictProR:::gp_asreml_refine_fit(
    mod,
    max_updates = 8L,
    update_fn = update_fn
  )

  expect_equal(update_count, 2L)
  expect_lte(abs(refined$vparameters.pc[["vg"]]), 1)
  expect_no_error(PredictProR:::gp_asreml_assert_trustworthy_fit(refined))
})

test_that("ASReml refinement is bounded and preserves fail-loud output checks", {
  mod <- list(
    converge = TRUE,
    ifault = 0L,
    vparameters.pc = c(vg = 2),
    vparameters.con = c("P")
  )
  update_count <- 0L
  unchanged_update <- function(x) {
    update_count <<- update_count + 1L
    x
  }

  refined <- PredictProR:::gp_asreml_refine_fit(
    mod,
    max_updates = 3L,
    update_fn = unchanged_update
  )

  expect_equal(update_count, 3L)
  expect_error(
    PredictProR:::gp_asreml_assert_trustworthy_fit(refined),
    "changing by more than 1%"
  )
})
