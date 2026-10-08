make_met_classification_cv_entry <- function(model, rep, pooled_accuracy,
                                             environment_accuracy,
                                             method = "CV1") {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  if (identical(model, "ModelA")) {
    prob_yes <- c(0.2, 0.7, 0.6, 0.3)
  } else {
    prob_yes <- c(0.1, 0.8, 0.3, 0.7)
  }
  pred <- data.frame(
    row_id = seq_along(y),
    GID = paste0("g", seq_along(y)),
    Env = rep(c("E1", "E2"), each = 2L),
    y = y,
    yhat = prob_yes,
    cv_role = "test",
    stringsAsFactors = FALSE
  )
  prob <- data.frame(
    no = 1 - prob_yes,
    yes = prob_yes,
    check.names = FALSE
  )
  list(
    trait = "Trait",
    rep = as.integer(rep),
    model = model,
    model_canonical = model,
    response_family = "binary",
    class_levels = c("no", "yes"),
    positive_class = "yes",
    cv_info = list(method = method),
    eval_metrics_reps = data.frame(
      Rep = as.integer(rep),
      Env = c("E1", "E2", "ALL"),
      balanced_accuracy = c(
        environment_accuracy,
        environment_accuracy,
        pooled_accuracy
      ),
      log_loss = c(0.5, 0.5, 1 - pooled_accuracy),
      stringsAsFactors = FALSE
    ),
    ypred_cv_Reps_all = pred,
    yprob_cv_Reps_all = prob
  )
}

test_that("single-trait MET classification exposes only purpose-fit models", {
  predictive <- PredictProR:::gp_met_supported_models()
  expect_setequal(
    PredictProR:::gp_met_classification_supported_models("binary"),
    c(predictive, "RKHS")
  )
  expect_setequal(
    PredictProR:::gp_met_classification_supported_models("multiclass"),
    predictive
  )
  expect_setequal(
    PredictProR:::gp_met_classification_supported_models("ordinal"),
    "RKHS"
  )

  standard <- PredictProR::met_data_standard()
  expect_setequal(
    standard$current_scope$response_families,
    c("gaussian", "binary", "multiclass", "ordinal")
  )
  expect_match(
    standard$current_scope$classification_cv_selection$ordinal,
    "quadratic-weighted kappa"
  )
})

test_that("MET classification guardrails reject structurally invalid routes", {
  pheno <- data.frame(
    GID = rep(paste0("g", 1:4), each = 2L),
    Env = rep(c("E1", "E2"), times = 4L),
    Class = factor(rep(c("A", "B", "C", "A"), each = 2L)),
    stringsAsFactors = FALSE
  )
  base <- list(
    pheno_data = pheno,
    response = "Class",
    gen_name = "GID",
    heter_groups = "Env",
    GS_model_cv = character()
  )

  multiclass_rkhs <- c(base, list(
    response_family = "multiclass",
    GS_model = "RKHS"
  ))
  expect_error(
    PredictProR:::gp_validate_single_trait_met_model_scope(
      multiclass_rkhs, is_met = TRUE
    ),
    "Unsupported.*multiclass.*RKHS.*stopped before fitting"
  )

  ordinal_catboost <- c(base, list(
    response_family = "ordinal",
    GS_model = "CatBoost"
  ))
  expect_error(
    PredictProR:::gp_validate_single_trait_met_model_scope(
      ordinal_catboost, is_met = TRUE
    ),
    "Unsupported.*ordinal.*CatBoost.*stopped before fitting"
  )

  ordinal_duplicate_gblup <- c(base, list(
    response_family = "ordinal",
    GS_model = "GBLUP_BRR",
    fixed = ~ Env,
    random = ~ GID + GID:Env
  ))
  expect_error(
    PredictProR:::gp_validate_single_trait_met_model_scope(
      ordinal_duplicate_gblup, is_met = TRUE
    ),
    "Unsupported.*ordinal.*GBLUP_BRR.*stopped before fitting"
  )

  ordinal_missing_gxe <- c(base, list(
    response_family = "ordinal",
    GS_model = "RKHS",
    fixed = ~ Env,
    random = ~ GID
  ))
  expect_error(
    PredictProR:::gp_validate_single_trait_met_model_scope(
      ordinal_missing_gxe, is_met = TRUE
    ),
    "require an explicit environment main effect and GxE structure"
  )

  ordinal_valid <- c(base, list(
    response_family = "ordinal",
    GS_model = "RKHS",
    fixed = ~ Env,
    random = ~ GID + GID:Env
  ))
  expect_invisible(
    PredictProR:::gp_validate_single_trait_met_model_scope(
      ordinal_valid, is_met = TRUE
    )
  )
})

test_that("MET classification ranks pooled OOF metrics without double-counting ALL", {
  entries <- list(
    make_met_classification_cv_entry("ModelA", 1L, 0.40, 0.90),
    make_met_classification_cv_entry("ModelA", 2L, 0.40, 0.90),
    make_met_classification_cv_entry("ModelB", 1L, 0.80, 0.50),
    make_met_classification_cv_entry("ModelB", 2L, 0.80, 0.50)
  )

  out <- PredictProR:::cv1_cv2_and_across_env_result_plot_process(
    entries,
    eval_metrics = c("balanced_accuracy", "log_loss"),
    heter_groups = "Env",
    metric_for_ranking = "auto",
    positive_class = "yes"
  )

  candidates <- out$classification_model_selection$candidate_metrics
  candidates <- candidates[order(candidates$model), , drop = FALSE]
  expect_identical(candidates$model, c("ModelA", "ModelB"))
  expect_equal(candidates$balanced_accuracy, c(0.40, 0.80))
  expect_identical(
    out$classification_model_selection$selected_models$model,
    "ModelB"
  )

  policy <- out$classification_model_selection$policy
  expect_true(policy$multi_environment)
  expect_identical(policy$environment_column, "Env")
  expect_identical(policy$cv_scenario, "CV1")
  expect_identical(
    policy$selection_scope,
    "pooled_out_of_fold_across_environments"
  )
  expect_identical(
    policy$environment_aggregation,
    "ALL_rows_mean_across_replications"
  )

  environment_metrics <- out$aggregated_data_list$balanced_accuracy$aggregated_data
  expect_false(any(as.character(environment_metrics$Env) == "ALL"))
  reliability <- out$classification_reliability
  expect_true("Env" %in% names(reliability$out_of_fold_predictions))
  expect_setequal(
    unique(reliability$summary_by_environment$Env),
    c("E1", "E2")
  )
  confusion <- out$classification_metrics$confusion_matrices
  expect_true(is.data.frame(confusion$per_environment_aggregated))
  expect_equal(
    sum(confusion$per_environment_aggregated$n_predictions),
    nrow(reliability$out_of_fold_predictions)
  )
})

test_that("fast MET classification processing uses the same pooled contract", {
  entries <- list(
    make_met_classification_cv_entry("ModelA", 1L, 0.40, 0.90, "CV2"),
    make_met_classification_cv_entry("ModelB", 1L, 0.80, 0.50, "CV2")
  )
  out <- PredictProR:::gp_cv_process_results_fast(
    entries,
    eval_metrics = c("balanced_accuracy", "log_loss"),
    heter_groups = "Env",
    is_multi = TRUE,
    metric_for_ranking = "balanced_accuracy",
    ranking_tie_breakers = "log_loss",
    positive_class = "yes"
  )
  expect_identical(out$best_models_list$balanced_accuracy$model, "ModelB")
  expect_equal(
    out$classification_metrics$aggregated$balanced_accuracy[
      order(out$classification_metrics$aggregated$model)
    ],
    c(0.40, 0.80)
  )
  expect_identical(
    out$classification_model_selection$policy$cv_scenario,
    "CV2"
  )
  expect_true(is.data.frame(
    out$classification_probability_summaries$aggregated_probabilities
  ))
})

test_that("classification selection policy does not call one model a comparison", {
  entry <- make_met_classification_cv_entry(
    "RKHS", 1L, pooled_accuracy = 0.60, environment_accuracy = 0.55
  )
  out <- PredictProR:::cv1_cv2_and_across_env_result_plot_process(
    list(entry),
    eval_metrics = c("balanced_accuracy", "log_loss"),
    heter_groups = "Env",
    metric_for_ranking = "balanced_accuracy",
    positive_class = "yes"
  )
  expect_false(
    out$classification_model_selection$policy$cv_compares_models
  )
})
