test_that("ASReml cross-validation uses trait-specific test indices under unbalanced traits", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(10, NA, 30, 40),
    Height = c(NA, 20, 30, 40),
    stringsAsFactors = FALSE
  )

  seen_tst <- list()

  local_mocked_bindings(
    hold_out_stratified_and_un = function(pheno_data, ...) list(1L),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) lapply(seq_along(model_ids), run_task),
    gp_filter_crossval_results = function(results) results,
    predict_with_model = function(model = NULL, y = NULL, tst = NULL, additional_params = NULL, ...) {
      seen_tst[[additional_params$response]] <<- tst
      rep(0.5, length(tst))
    },
    .package = "PredictProR"
  )

  res <- PredictProR:::models_execute_crossval(
    pheno_data = pheno,
    params = list(
      response = c("Yield", "Height"),
      gen_name = "GID",
      test_size = 0.5,
      random_state = 1,
      replication = 1L,
      cross_validation_meth = "hold_out",
      eval_metrics = "mean_squared_error",
      GS_model_cv = "GBLUP",
      engine = "asreml"
    ),
    verbose = FALSE
  )

  expect_length(res, 2)
  expect_identical(seen_tst[["Yield"]], 1L)
  expect_identical(seen_tst[["Height"]], 2L)
})

test_that("stacked multi-trait ASReml data are unit-major and trait-minor", {
  pheno <- data.frame(
    GID = c("g2", "g1", "g3"),
    Yield = c(2, NA, 3),
    Height = c(20, 10, NA),
    stringsAsFactors = FALSE
  )

  long <- PredictProR:::gp_multitrait_asreml_long_data(
    pheno_object = pheno,
    response = c("Yield", "Height"),
    gen_name = "GID"
  )

  expect_identical(long$GID, rep(pheno$GID, each = 2L))
  expect_identical(long$Trait, rep(c("Yield", "Height"), times = nrow(pheno)))
  expect_identical(long$Observed_value, c(2, 20, NA, 10, 3, NA))
})

test_that("ASReml best-model routing preserves trait-specific true-prediction labels", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    Height = c(NA, 2, 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    pheno_clean <- PredictProR::phenotype_to_model(
      pheno_data = pheno,
      response = c("Yield", "Height"),
      gen_name = "GID"
    ),
    "trait-specific inferred testing sets"
  )

  seen_test_ids <- list()
  seen_summary_test_ids <- list()

  local_mocked_bindings(
    asreml_utilis_new = function(...) list(model = structure(list(), class = "mock_asreml")),
    asreml_varcomp_table = function(...) {
      data.frame(
        Components = 1,
        bound = "P",
        row.names = "genetic_variance",
        stringsAsFactors = FALSE
      )
    },
    asreml_mod_output_new = function(pheno_data = NULL, response = NULL, gen_name = NULL, ...) {
      lbl <- ifelse(is.na(pheno_data[[response]]), "Test", "Train")
      seen_test_ids[[response]] <<- as.character(pheno_data[[gen_name]][lbl == "Test"])
      list(
        Predicted_value = data.frame(
          GID = pheno_data[[gen_name]],
          Predicted_value = seq_len(nrow(pheno_data)),
          Standard_error = rep(0.1, nrow(pheno_data)),
          Prediction_error_variance = rep(0.01, nrow(pheno_data)),
          Train_Test_Label = lbl,
          stringsAsFactors = FALSE
        ),
        Variance_components = data.frame(
          Components = 1,
          row.names = "genetic_variance",
          stringsAsFactors = FALSE
        ),
        Asreml_model = list(model = "mock")
      )
    },
    summary_statistics_asreml = function(pheno_data = NULL, response = NULL, ...) {
      seen_summary_test_ids[[response]] <<- as.character(pheno_data$GID[is.na(pheno_data[[response]])])
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )

  base_ctx <- list(
    GS_model = "GBLUP",
    engine = "asreml",
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    AI_valid_models = character(),
    canonical_names = character(),
    pheno_clean = pheno_clean,
    ml_dat_res = list(),
    gen_name = "GID",
    heter_groups = NULL,
    fixed = NULL,
    random = ~ GID,
    cova = NULL,
    weights = NULL,
    gmatrix_model_ready = diag(3),
    omic1_kernel_model_ready = NULL,
    omic2_kernel_model_ready = NULL,
    omic3_kernel_model_ready = NULL,
    heter_resid = FALSE,
    var_cov_str = NULL,
    pworkspace = 1e06,
    workspace = 1e08,
    maxit = 10,
    inverse = TRUE,
    epsilon = 1e-6,
    eval_metrics = "mean_squared_error",
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.7,
    low_reliability_thres = 0.4,
    system_database = FALSE,
    friendly_name_lookup = c(GBLUP = "GBLUP"),
    msg = "",
    cross_validation = FALSE
  )
  rownames(base_ctx$gmatrix_model_ready) <- colnames(base_ctx$gmatrix_model_ready) <- c("g1", "g2", "g3")

  yield_res <- PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Yield")))
  height_res <- PredictProR:::gp_run_best_model_task(utils::modifyList(base_ctx, list(response = "Height")))

  expect_identical(unname(seen_test_ids[["Yield"]]), "g2")
  expect_identical(unname(seen_test_ids[["Height"]]), "g1")
  expect_identical(unname(seen_summary_test_ids[["Yield"]]), "g2")
  expect_identical(unname(seen_summary_test_ids[["Height"]]), "g1")
  expect_identical(
    as.character(yield_res$res_model_output$Predicted_value$Train_Test_Label),
    c("Train", "Test", "Train")
  )
  expect_identical(
    as.character(height_res$res_model_output$Predicted_value$Train_Test_Label),
    c("Test", "Train", "Train")
  )
})
