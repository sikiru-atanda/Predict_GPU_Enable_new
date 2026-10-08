test_that("GP-ASReml benchmark scenarios cover the requested design matrix", {
  scenarios <- PredictProR:::gp_asreml_benchmark_scenarios()

  expect_identical(
    scenarios$scenario,
    c(
      "single_trait_single_environment",
      "single_trait_met_single_kernel",
      "single_trait_met_multi_kernel",
      "multi_trait_single_environment",
      "multi_trait_multi_environment"
    )
  )
  expect_identical(scenarios$n_kernels, c(1L, 1L, 2L, 1L, 1L))
  expect_identical(scenarios$n_traits, c(1L, 1L, 1L, 2L, 2L))
})

test_that("benchmark predictions normalize GP and ASReml aliases identically", {
  gp <- data.frame(
    Name = c("g1", "g2"),
    Env = c("E1", "E1"),
    Prediction = c(1.2, 2.3),
    Prediction_SE_observed = c(0.2, 0.3)
  )
  asr <- data.frame(
    GID = c("g1", "g2"),
    Env = c("E1", "E1"),
    predicted.value = c(1.1, 2.2),
    std.error = c(0.25, 0.35),
    check.names = FALSE
  )

  gp_out <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    gp,
    scenario = "single_trait_met_single_kernel",
    dataset = "synthetic_1",
    engine = "PredictProR_GP",
    model_family = "GP_FA",
    observed = c(1, 2),
    default_trait = "Yield"
  )
  asr_out <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    asr,
    scenario = "single_trait_met_single_kernel",
    dataset = "synthetic_1",
    engine = "ASReml_R",
    model_family = "FA1_GBLUP",
    observed = c(1, 2),
    default_trait = "Yield"
  )

  expect_identical(names(gp_out), names(asr_out))
  expect_identical(names(gp_out), PredictProR:::gp_asreml_benchmark_prediction_columns())
  expect_equal(gp_out$PEV, gp_out$Standard_error^2)
  expect_equal(asr_out$PEV, asr_out$Standard_error^2)
  expect_true(PredictProR:::gp_validate_asreml_benchmark_predictions(gp_out)$valid)
  expect_true(PredictProR:::gp_validate_asreml_benchmark_predictions(asr_out)$valid)
})

test_that("benchmark validation rejects duplicated keys and non-finite predictions", {
  out <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    data.frame(GID = c("g1", "g2"), Prediction = c(1, 2)),
    scenario = "single_trait_single_environment",
    dataset = "synthetic_1",
    engine = "PredictProR_GP",
    model_family = "GP",
    observed = c(1, 2),
    default_env = "E1",
    default_trait = "Yield"
  )
  bad <- rbind(out, out[1L, , drop = FALSE])
  bad$Predicted_value[[2L]] <- Inf
  validation <- PredictProR:::gp_validate_asreml_benchmark_predictions(bad)

  expect_false(validation$valid)
  expect_match(paste(validation$errors, collapse = " "), "non-finite")
  expect_match(paste(validation$errors, collapse = " "), "duplicated")
})

test_that("benchmark normalization uses observed-scale variance consistently", {
  out <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    data.frame(
      GID = paste0("g", seq_len(3L)),
      Prediction = 1:3,
      Prediction_SE_observed = rep(2, 3L),
      Prediction_Var_observed = rep(4, 3L),
      Prediction_SE_latent = rep(1, 3L),
      PEV = rep(1, 3L)
    ),
    scenario = "multi_trait_single_environment",
    dataset = "synthetic_1",
    engine = "PredictProR_GP",
    model_family = "MT_GP",
    observed = 1:3,
    default_env = "E1",
    default_trait = "Trait1"
  )

  expect_equal(out$Standard_error, rep(2, 3L))
  expect_equal(out$PEV, rep(4, 3L))
  expect_true(PredictProR:::gp_validate_asreml_benchmark_predictions(out)$valid)

  out$PEV[[1L]] <- 1
  validation <- PredictProR:::gp_validate_asreml_benchmark_predictions(out)
  expect_false(validation$valid)
  expect_match(paste(validation$errors, collapse = " "), "Standard_error\\^2")
})

test_that("covariance diagnostics check PSD and truth error", {
  truth <- matrix(c(1, 0.4, 0.4, 0.8), 2L)
  good <- PredictProR:::gp_asreml_benchmark_covariance_diagnostics(truth, truth)
  bad <- PredictProR:::gp_asreml_benchmark_covariance_diagnostics(
    matrix(c(1, 2, 2, 1), 2L), truth
  )

  expect_true(good$positive_semidefinite)
  expect_equal(good$relative_frobenius_truth_error, 0)
  expect_false(bad$positive_semidefinite)
})

test_that("benchmark accuracy is stratified by applicable scenario dimensions", {
  base <- expand.grid(
    GID = paste0("g", seq_len(4L)),
    Env = c("E1", "E2"),
    Trait = c("T1", "T2"),
    stringsAsFactors = FALSE
  )
  base$Observed_value <- seq_len(nrow(base))

  make_predictions <- function(engine, offset) {
    PredictProR:::gp_asreml_benchmark_standardize_predictions(
      data.frame(
        GID = base$GID,
        Env = base$Env,
        Trait = base$Trait,
        Predicted_value = base$Observed_value + offset,
        stringsAsFactors = FALSE
      ),
      scenario = "multi_trait_multi_environment",
      dataset = "synthetic_1",
      engine = engine,
      model_family = if (engine == "PredictProR_GP") "MT_MET_GP" else "CORGH_GBLUP",
      observed = base$Observed_value
    )
  }
  predictions <- rbind(
    make_predictions("PredictProR_GP", 1),
    make_predictions("ASReml_R", 2)
  )
  tables <- PredictProR:::gp_asreml_benchmark_metric_tables(predictions)

  expect_identical(
    names(tables),
    c("overall", "by_environment", "by_trait", "by_environment_trait", "stratified")
  )
  expect_identical(names(tables$overall), c(
    "scenario", "dataset", "engine", "n", "RMSE", "MAE", "Pearson", "Spearman"
  ))
  expect_equal(nrow(tables$by_environment), 0L)
  expect_equal(nrow(tables$by_trait), 4L)
  expect_equal(nrow(tables$by_environment_trait), 8L)
  expect_true(all(tables$by_trait$n == 8L))
  expect_true(all(tables$by_environment_trait$n == 4L))
  expect_equal(
    tables$overall$RMSE[tables$overall$engine == "PredictProR_GP"],
    1
  )
  expect_equal(
    tables$by_environment_trait$RMSE[
      tables$by_environment_trait$engine == "ASReml_R"
    ],
    rep(2, 4L)
  )
  expect_setequal(
    unique(tables$stratified$metric_scope),
    c("overall", "trait", "environment_trait")
  )
})

test_that("single-trait MET accuracy is reported separately by environment", {
  base <- expand.grid(
    GID = paste0("g", seq_len(4L)),
    Env = c("E1", "E2"),
    stringsAsFactors = FALSE
  )
  predictions <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    data.frame(
      GID = base$GID,
      Env = base$Env,
      Prediction = seq_len(nrow(base)) + 0.5
    ),
    scenario = "single_trait_met_single_kernel",
    dataset = "synthetic_1",
    engine = "PredictProR_GP",
    model_family = "GP_FA",
    observed = seq_len(nrow(base)),
    default_trait = "Yield"
  )
  tables <- PredictProR:::gp_asreml_benchmark_metric_tables(predictions)

  expect_equal(nrow(tables$by_environment), 2L)
  expect_equal(tables$by_environment$n, c(4L, 4L))
  expect_equal(tables$by_environment$RMSE, c(0.5, 0.5))
  expect_false("Trait" %in% names(tables$by_environment))
})

test_that("stratified metrics omit inapplicable scenario slices", {
  predictions <- PredictProR:::gp_asreml_benchmark_standardize_predictions(
    data.frame(GID = paste0("g", seq_len(4L)), Prediction = 2:5),
    scenario = "single_trait_single_environment",
    dataset = "synthetic_1",
    engine = "PredictProR_GP",
    model_family = "GP",
    observed = 1:4,
    default_env = "E1",
    default_trait = "Yield"
  )
  tables <- PredictProR:::gp_asreml_benchmark_metric_tables(predictions)

  expect_equal(nrow(tables$overall), 1L)
  expect_equal(nrow(tables$by_environment), 0L)
  expect_equal(nrow(tables$by_trait), 0L)
  expect_equal(nrow(tables$by_environment_trait), 0L)
  expect_identical(names(tables$by_environment), c(
    "scenario", "dataset", "engine", "Env", "n", "RMSE", "MAE", "Pearson", "Spearman"
  ))
})
