test_that("GP CV routes each trait response and keeps public GP labels", {
  ph <- data.frame(
    GID = paste0("g", seq_len(6)),
    Trait1 = c(NA, 2, 3, 4, 5, 6),
    Trait2 = c(1, 2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  K <- diag(6)
  rownames(K) <- colnames(K) <- ph$GID
  captured_response <- character()
  captured_model <- character()

  local_mocked_bindings(
    stop_daemons = function(...) NULL,
    gp_detect_python = function(...) NULL,
    gp_init_python_once = function(...) NULL,
    gp_backend_prepare_kernel_cache = function(...) NULL,
    gp_backend_cv_batch_enabled = function() TRUE,
    gp_backend_cv_predict_batch = function(model_name, y, folds, additional_params = NULL) {
      captured_model <<- c(captured_model, model_name)
      captured_response <<- c(captured_response, additional_params$response)
      testthat::expect_length(additional_params$response, 1L)
      lapply(folds, function(idx) rep(mean(y, na.rm = TRUE), length(idx)))
    },
    gp_parallel_extract_model_policy_params = function(...) list(),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    hold_out_stratified_and_un = function(pheno_data, ...) {
      list(seq_len(min(2L, nrow(pheno_data))))
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::models_execute_crossval(
    pheno_data = ph,
    params = list(
      response = c("Trait1", "Trait2"),
      gen_name = "GID",
      test_size = 0.25,
      random_state = 1,
      replication = 1,
      GS_model_cv = "Kernel-GBLUP",
      cross_validation_meth = "Stratified_Hold_Out",
      nfolds = 2,
      eval_metrics = "mean_squared_error",
      response_family = "gaussian",
      gmatrix_model_ready = K,
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )

  expect_setequal(captured_response, c("Trait1", "Trait2"))
  expect_true(all(captured_model == "KRR"))
  expect_equal(length(out), 2L)
  expect_true(all(vapply(out, function(x) identical(x$model, "Kernel-GBLUP"), logical(1))))
  expect_true(all(vapply(out, function(x) identical(x$model_canonical, "KRR"), logical(1))))
})

test_that("GP-only CV evaluation keeps user-facing plots by default", {
  ctx <- list(
    cv_generate_plots = NULL,
    cv_evaluation_only = TRUE,
    GS_model_cv = PredictProR:::gp_lowrank_supported_models(),
    gp_valid_models = PredictProR:::gp_lowrank_supported_models()
  )

  expect_true(PredictProR:::gp_cv_generate_plots_enabled(ctx))

  ctx$cv_generate_plots <- FALSE
  expect_false(PredictProR:::gp_cv_generate_plots_enabled(ctx))
})

test_that("trait prediction test summary is distinct from CV holdout size", {
  pheno <- data.frame(
    ID = letters[1:5],
    T1 = c(1, NA, 3, NA, 5),
    T2 = c(NA, 2, 3, 4, 5),
    stringsAsFactors = FALSE
  )
  pheno_clean <- list(
    pheno_clean_data = pheno,
    test_set = c("e", "z"),
    test_set_by_trait = list(T1 = c("b", "d"), T2 = "a")
  )

  summary <- PredictProR:::gp_trait_prediction_test_set_summary(
    pheno_clean = pheno_clean,
    response = c("T1", "T2"),
    gen_name = "ID"
  )

  expect_equal(summary$prediction_test_n[summary$trait == "T1"], 4)
  expect_equal(summary$missing_response_test_n[summary$trait == "T1"], 2)
  expect_equal(summary$global_test_n[summary$trait == "T1"], 2)
  expect_equal(summary$training_records_n[summary$trait == "T1"], 2)

  cv_results <- list(list(
    trait = "T1",
    model = "Kernel-GBLUP",
    ypred_cv_Reps_all = data.frame(
      y = c(1, 3, 5),
      yhat = c(NA, 3.1, NA),
      cv_role = c("train", "test", "train")
    )
  ))

  annotated <- PredictProR:::gp_attach_trait_prediction_test_info_to_cv_results(
    cv_results = cv_results,
    pheno_clean = pheno_clean,
    response = c("T1", "T2"),
    gen_name = "ID"
  )

  expect_equal(annotated[[1]]$cv_info$n_cv_test, 1)
  expect_equal(annotated[[1]]$cv_info$n_prediction_test, 4)
  expect_true(annotated[[1]]$cv_info$n_test_is_cv_holdout)
})

test_that("GP true prediction values include PEV-based reliability metadata", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, NA, 3, 5),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = pheno$GID,
    Predicted_value = c(1.1, 2.2, 3.1, 4.8),
    Prediction_SE_observed = c(0.30, 0.40, 0.35, 0.30),
    Prediction_Var_observed = c(0.09, 0.16, 0.1225, 0.09),
    Prediction_SE_latent = c(0.10, 0.20, 0.15, 0.10),
    Prediction_Var_latent = c(0.01, 0.04, 0.0225, 0.01),
    stringsAsFactors = FALSE
  )
  gp_result <- list(
    per_env = data.frame(
      Env = "ENV1",
      total_genetic_variance = 2,
      stringsAsFactors = FALSE
    )
  )

  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = predictions,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    test_set = "g2",
    gp_result = gp_result
  )

  expect_true(all(c(
    "Standard_error",
    "PEV",
    "Genetic_variance",
    "Reliability_reference_variance",
    "Reliability",
    "Reliability_remarks",
    "Reliability_percentage",
    "Reliability_basis",
    "Uncertainty_remarks"
  ) %in% names(out)))
  expect_equal(out$Standard_error, predictions$Prediction_SE_latent)
  expect_equal(out$PEV, predictions$Prediction_Var_latent)
  expect_equal(out$Genetic_variance, rep(2, 4))
  expect_equal(out$Reliability_reference_variance, rep(2, 4))
  expect_equal(out$Reliability[2], 1 - 0.04 / 2)
  expect_identical(out$Reliability_basis, rep("reml_genetic_variance", 4))
  expect_identical(
    out$Uncertainty_remarks,
    rep("model_based_prediction_interval", 4)
  )
})

test_that("GP reliability is NA (no phenotypic fallback) when REML genetic variance is not estimable", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3", "g4"),
    Yield = c(1, NA, 3, 5),
    stringsAsFactors = FALSE
  )
  predictions <- data.frame(
    GID = pheno$GID,
    Predicted_value = c(1.1, 2.2, 3.1, 4.8),
    Prediction_SE_observed = c(0.30, 0.40, 0.35, 0.30),
    Prediction_Var_observed = c(0.09, 0.16, 0.1225, 0.09),
    Prediction_SE_latent = c(0.10, 0.20, 0.15, 0.10),
    Prediction_Var_latent = c(0.01, 0.04, 0.0225, 0.01),
    stringsAsFactors = FALSE
  )
  gp_result_no_vc <- list()

  out <- PredictProR:::gp_model_execute_single_trait_predicted_values(
    predictions = predictions,
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    test_set = "g2",
    gp_result = gp_result_no_vc
  )

  expect_true(all(is.na(out$Genetic_variance)))
  expect_true(all(is.na(out$Reliability_reference_variance)))
  expect_true(all(is.na(out$Reliability)))
  expect_identical(
    out$Reliability_basis,
    rep("reml_genetic_variance_not_estimable", nrow(out))
  )
  expect_identical(
    out$Reliability_remarks,
    rep("sigma2_g not estimable", nrow(out))
  )

  train_y <- pheno$Yield[!is.na(pheno$Yield) & pheno$GID != "g2"]
  phenotypic_var <- stats::var(train_y)
  reliabilities_under_phenotypic_fallback <-
    1 - predictions$Prediction_Var_latent / phenotypic_var
  expect_false(isTRUE(all.equal(
    reliabilities_under_phenotypic_fallback,
    out$Reliability
  )))
})

test_that("GP final prediction after CV requests uncertainty output", {
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Yield = c(1, NA, 3),
    stringsAsFactors = FALSE
  )
  K <- diag(3)
  rownames(K) <- colnames(K) <- pheno$GID
  seen <- list()

  local_mocked_bindings(
    gp_backend_gaussian_model = function(gp_output_level = NULL,
                                         gp_prediction_output = NULL,
                                         test_set = NULL,
                                         test_set_source = NULL,
                                         ...) {
      seen$gp_output_level <<- gp_output_level
      seen$gp_prediction_output <<- gp_prediction_output
      seen$test_set <<- test_set
      seen$test_set_source <<- test_set_source
      list(
        predicted_values = data.frame(
          GID = pheno$GID,
          Predicted_value = c(1.1, 2.2, 3.3),
          Train_Test_Label = c("Train", "Test", "Train"),
          Standard_error = c(0.1, 0.2, 0.1),
          PEV = c(0.01, 0.04, 0.01),
          Observed_value = pheno$Yield,
          stringsAsFactors = FALSE
        ),
        model_parameters = data.frame(stat = "gp_model", summary = "KRR"),
        lowrank_feature_matrix = K,
        diagnostic_plots = list(predicted = "plot")
      )
    },
    summary_statistics_AI = function(predicted_object = NULL, ...) {
      testthat::expect_true(all(c("Standard_error", "PEV") %in% names(predicted_object)))
      list(summary_statistics = data.frame(stat = "ok", summary = "ok"))
    },
    .package = "PredictProR"
  )

  ctx <- list(
    response = "Yield",
    GS_model = "KRR",
    pheno_clean = list(
      pheno_clean_data = pheno,
      test_set_by_trait = list(Yield = "g2"),
      test_set_by_trait_source = "inferred_missing_response"
    ),
    ml_dat_res = list(),
    gen_name = "GID",
    gmatrix_model_ready = K,
    gmatrix_kernel_model_ready_list = list(gmatrix_model_ready = K),
    gp_output_level = "predict_only",
    gp_prediction_output = "test_only",
    cross_validation = TRUE,
    cv_evaluation_only = FALSE,
    eval_metrics = "mean_squared_error",
    response_family = "gaussian",
    system_database = FALSE,
    friendly_name_lookup = PredictProR:::gp_model_friendly_lookup(),
    AI_valid_models = character(),
    bayes_valid_models = c("BRR", "BayesA", "BayesB", "BayesC", "BL"),
    bayes_gblup_valid_models = c("GBLUP_BRR", "RKHS"),
    msg = ""
  )

  out <- PredictProR:::gp_run_best_model_task(ctx)

  expect_identical(seen$gp_output_level, "full_vc")
  expect_identical(seen$gp_prediction_output, "all")
  expect_identical(seen$test_set, "g2")
  expect_identical(seen$test_set_source, "inferred_missing_response")
  expect_true("predicted_values" %in% names(out$res_model_output))
})

test_that("model prediction outputs expose a consistent predicted_values table", {
  bayes_like <- list(
    bayes_result = list(
      Predicted_value = data.frame(
        NAME = c("g1", "g2"),
        BLUP = c(1.2, 2.4),
        SE = c(0.11, 0.22),
        Prediction_error_variance = c(0.0121, 0.0484),
        stringsAsFactors = FALSE
      )
    )
  )

  out <- PredictProR:::gp_standardize_model_prediction_outputs(
    bayes_like,
    gen_name = "NAME"
  )

  expect_true("predicted_values" %in% names(out))
  expect_true(all(c("GID", "Predicted_value", "Standard_error", "PEV") %in% names(out$predicted_values)))
  expect_identical(out$predicted_values$GID, c("g1", "g2"))
  expect_equal(out$predicted_values$Predicted_value, c(1.2, 2.4))
  expect_equal(out$predicted_values$Standard_error, c(0.11, 0.22))
  expect_equal(out$predicted_values$PEV, c(0.0121, 0.0484))
  expect_true("predicted_values" %in% names(out$bayes_result))
})

test_that("true-prediction task table reduces requested GP CV models to one best per trait", {
  # Pre-0.20.9 the GP family was unique: when all requested CV models were
  # in the GP-lowrank family, the true-prediction task table expanded into
  # a trait x model cross-product so EVERY GP model was re-fit on EVERY
  # trait. Every other family (Bayes / ASReml / ML / DL) already took only
  # the best model per trait. The GP cross-product turned a 4-trait /
  # 4-model run into 16 final fits instead of 4. This test pins the new
  # uniform contract: one best model per trait for every family,
  # including the GP family.
  cv_results_processed <- list(
    aggregated_data_list = list(
      aggregated_across_reps = data.frame(
        trait = rep(c("Yield", "Height"), each = 2),
        model = rep(c("Kernel-GBLUP", "Scalable-GBLUP"), times = 2),
        accuracy = c(0.5, 0.4, 0.6, 0.7),
        stringsAsFactors = FALSE
      )
    ),
    best_models_list = list(
      accuracy = data.frame(
        trait = c("Yield", "Height"),
        model = c("Kernel-GBLUP", "Scalable-GBLUP"),
        accuracy = c(0.5, 0.7),
        stringsAsFactors = FALSE
      )
    )
  )

  tasks <- PredictProR:::gp_build_true_prediction_task_table(
    best_models = cv_results_processed$best_models_list$accuracy,
    GS_model = NULL,
    GS_model_cv = c("Kernel-GBLUP", "Scalable-GBLUP"),
    response = c("Yield", "Height"),
    cv_results_processed = cv_results_processed,
    metric_for_ranking = "accuracy",
    cross_validation = TRUE
  )

  expect_equal(nrow(tasks), 2L)
  expect_identical(tasks$trait, c("Yield", "Height"))
  expect_identical(tasks$model, c("KRR", "LowRankGP"))
  if ("is_cv_best" %in% names(tasks)) {
    expect_true(all(tasks$is_cv_best))
  }
})
