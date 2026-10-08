make_variance_contract_gaussian_predictions <- function(include_trait = FALSE,
                                                          include_env = FALSE) {
  dimensions <- list(GID = paste0("g", seq_len(4L)))
  if (isTRUE(include_trait)) dimensions$Trait <- c("T1", "T2")
  if (isTRUE(include_env)) dimensions$Env <- c("E1", "E2")
  out <- do.call(
    expand.grid,
    c(dimensions, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  )
  within_group <- ave(seq_len(nrow(out)), interaction(
    if (isTRUE(include_trait)) out$Trait else "trait",
    if (isTRUE(include_env)) out$Env else "env",
    drop = TRUE
  ), FUN = seq_along)
  out$Predicted_value <- as.numeric(within_group) / 4 +
    if (isTRUE(include_trait)) as.numeric(factor(out$Trait)) else 0 +
    if (isTRUE(include_env)) 2 * as.numeric(factor(out$Env)) else 0
  out$Observed_value <- out$Predicted_value + rep(c(-0.1, 0.1), length.out = nrow(out))
  out$Train_Test_Label <- ifelse(within_group == 4L, "Test", "Train")
  out$PEV <- 0.04 + as.numeric(within_group) / 100
  out
}

make_variance_contract_classification_predictions <- function(include_trait = FALSE,
                                                               include_env = FALSE) {
  out <- make_variance_contract_gaussian_predictions(include_trait, include_env)
  out$Predicted_value <- NULL
  out$Observed_value <- NULL
  out$PEV <- NULL
  within_group <- ave(seq_len(nrow(out)), interaction(
    if (isTRUE(include_trait)) out$Trait else "trait",
    if (isTRUE(include_env)) out$Env else "env",
    drop = TRUE
  ), FUN = seq_along)
  out$Probability_A <- c(0.8, 0.65, 0.35, 0.2)[within_group]
  out$Probability_B <- 1 - out$Probability_A
  out$Predicted_class <- ifelse(out$Probability_A >= 0.5, "A", "B")
  out$Observed_class <- c("A", "A", "B", NA_character_)[within_group]
  out$Prediction_confidence <- pmax(out$Probability_A, out$Probability_B)
  out$Classification_uncertainty <- 1 - out$Prediction_confidence
  out$Reliability <- out$Prediction_confidence
  out$Reliability_remarks <- "predictive class confidence; not genetic reliability"
  out
}

expect_nonidentifiable_variance_contract <- function(vc, include_trait = FALSE,
                                                     include_env = FALSE,
                                                     info = NULL) {
  expect_true(is.data.frame(vc), info = info)
  expect_true(nrow(vc) > 0L, info = info)
  expect_true(all(c(
    "genetic_variance_not_identifiable",
    "residual_variance_not_identifiable",
    "heritability_not_identifiable"
  ) %in% vc$Component), info = info)
  expect_true(
    all(is.na(vc$Components[grepl("not_identifiable$", vc$Component)])),
    info = info
  )
  expect_identical("Trait" %in% names(vc), isTRUE(include_trait), info = info)
  expect_identical("Env" %in% names(vc), isTRUE(include_env), info = info)
}

finalize_variance_contract_branch <- function(model, payload, response_family,
                                               cv_selected) {
  result <- list(
    trait = "Yield",
    response = "Yield",
    GS_model = model,
    model = model,
    is_cv_best = isTRUE(cv_selected),
    res_model_output = payload,
    res_summary_stat = list(
      summary_statistics = data.frame(
        stat = c("trait", "requested_model", "response_family"),
        summary = c("Yield", model, response_family),
        stringsAsFactors = FALSE
      )
    )
  )
  out <- PredictProR:::gp_finalize_model_execute_results(
    results = list(result),
    GS_model = model,
    gen_name = "GID",
    heter_groups = NULL,
    pheno_data = NULL,
    response_family = response_family,
    best_models_ggplot_rep = NULL,
    best_models_ggplot_mean = NULL,
    cv_results_predicted_vs_observed = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = list(best_models_list = list()),
    cv_results = list(),
    run_metadata = list(
      task_unit = if (isTRUE(cv_selected)) "cv_best_model_trait" else "model_trait"
    ),
    system_database = TRUE,
    plot_extension = "pdf",
    plot_width = 7,
    plot_height = 5,
    plot_units = "in",
    plot_dpi = 72
  )
  out$model_results_by_model[[model]][["Yield"]]
}

test_that("variance capability registry covers every supported model and alias", {
  predictive <- unique(c(
    PredictProR:::gp_classical_ml_supported_models(),
    PredictProR:::gp_deep_learning_supported_models(),
    PredictProR:::gp_deep_learning_friendly_models()
  ))
  model_based <- unique(c(
    "GBLUP", "GBLUP_BRR", "RKHS", "BRR", "BayesA", "BayesB", "BayesC", "BL",
    PredictProR:::gp_lowrank_supported_models(),
    PredictProR:::gp_lowrank_friendly_models()
  ))

  expect_true(all(vapply(
    predictive,
    function(model) identical(PredictProR:::gp_model_variance_role(model), "predictive_only"),
    logical(1L)
  )))
  expect_true(all(vapply(
    model_based,
    function(model) identical(PredictProR:::gp_model_variance_role(model), "model_based"),
    logical(1L)
  )))
  expect_identical(PredictProR:::gp_model_variance_role("unregistered-model"), "unknown")
})

test_that("ML and DL Gaussian variance semantics are workflow invariant", {
  layouts <- list(
    single_trait_single_environment = c(FALSE, FALSE),
    multi_trait_single_environment = c(TRUE, FALSE),
    met_single_trait = c(FALSE, TRUE),
    multi_trait_met = c(TRUE, TRUE),
    multi_omics = c(FALSE, FALSE),
    feature_selection = c(FALSE, FALSE)
  )
  for (layout_name in names(layouts)) {
    layout <- layouts[[layout_name]]
    pred <- make_variance_contract_gaussian_predictions(layout[[1L]], layout[[2L]])
    for (model in c("RandomForest", "TabTransformer")) {
      vc <- PredictProR:::gp_public_variance_components(
        list(predicted_values = pred),
        model = model,
        response_family = "gaussian"
      )
      expect_nonidentifiable_variance_contract(
        vc,
        include_trait = layout[[1L]],
        include_env = layout[[2L]],
        info = paste(layout_name, model)
      )
      expect_true(all(c(
        "phenotypic_reference_variance", "prediction_variance",
        "estimated_prediction_mse"
      ) %in% vc$Component), info = paste(layout_name, model))
    }
  }
})

test_that("ML and DL classification variance limits are explicit in every table layout", {
  layouts <- list(
    single_trait_single_environment = c(FALSE, FALSE),
    multi_trait_single_environment = c(TRUE, FALSE),
    met_single_trait = c(FALSE, TRUE),
    multi_trait_met = c(TRUE, TRUE)
  )
  for (layout_name in names(layouts)) {
    layout <- layouts[[layout_name]]
    pred <- make_variance_contract_classification_predictions(layout[[1L]], layout[[2L]])
    for (model in c("CatBoost", "DenseNeuralNet")) {
      vc <- PredictProR:::gp_public_variance_components(
        list(predicted_values = pred),
        model = model,
        response_family = "binary"
      )
      expect_nonidentifiable_variance_contract(
        vc,
        include_trait = layout[[1L]],
        include_env = layout[[2L]],
        info = paste(layout_name, model)
      )
      expect_true(all(c(
        "mean_prediction_confidence", "prediction_confidence_variance",
        "mean_classification_uncertainty"
      ) %in% vc$Component), info = paste(layout_name, model))
      expect_true(any(is.finite(vc$Components[!grepl("not_identifiable$", vc$Component)])))
    }
  }
})

test_that("model-based outputs preserve fitted variance and never use ML/DL labels", {
  pred <- make_variance_contract_gaussian_predictions(include_trait = TRUE, include_env = TRUE)
  native <- data.frame(
    Trait = rep(c("T1", "T2"), each = 3L),
    Component = rep(c("genetic_variance", "residual_variance", "heritability"), 2L),
    Components = c(0.6, 0.4, 0.6, 0.8, 0.2, 0.8),
    Standard_error = c(0.05, 0.04, 0.03, 0.06, 0.03, 0.02),
    stringsAsFactors = FALSE
  )
  models <- PredictProR:::gp_model_based_variance_models()
  for (model in models) {
    vc <- PredictProR:::gp_public_variance_components(
      list(predicted_values = pred, variance_components = native),
      model = model,
      response_family = "gaussian"
    )
    expect_true(all(c("genetic_variance", "residual_variance", "heritability") %in% vc$Component))
    expect_false(any(grepl("not_identifiable$", vc$Component)), info = model)

    missing <- PredictProR:::gp_public_variance_components(
      list(predicted_values = pred),
      model = model,
      response_family = "gaussian"
    )
    expect_false(any(grepl("not_identifiable$", missing$Component)), info = model)
    expect_identical(missing$Component, "model_variance_components_not_returned")
  }
})

test_that("Bayesian classification liability variance remains authoritative", {
  pred <- make_variance_contract_classification_predictions()
  liability <- data.frame(
    Component = c(
      "liability_scale_genetic_variance",
      "liability_scale_residual_variance_fixed",
      "liability_scale_heritability"
    ),
    Components = c(1.5, 1, 0.6),
    Standard_error = c(0.2, 0, 0.04),
    stringsAsFactors = FALSE
  )
  for (model in c("GBLUP_BRR", "RKHS")) {
    vc <- PredictProR:::gp_public_variance_components(
      list(predicted_values = pred, variance_components = liability),
      model = model,
      response_family = "binary"
    )
    expect_identical(vc$Component, liability$Component)
    expect_equal(vc$Components, liability$Components)
    expect_false(any(grepl("not_identifiable$", vc$Component)))
  }
})

test_that("unknown model identity cannot trigger the ML/DL variance fallback", {
  pred <- make_variance_contract_gaussian_predictions()
  vc <- PredictProR:::gp_public_variance_components(list(predicted_values = pred))
  expect_true(is.data.frame(vc))
  expect_equal(nrow(vc), 0L)
})

test_that("variance capability is recovered from model metadata automatically", {
  gaussian_pred <- make_variance_contract_gaussian_predictions()
  ml_result <- list(
    predicted_values = gaussian_pred,
    model_parameters = data.frame(
      stat = c("requested_model", "response_family"),
      summary = c("Ridge_Regression", "gaussian"),
      stringsAsFactors = FALSE
    )
  )
  ml_vc <- PredictProR:::gp_public_variance_components(ml_result)
  expect_nonidentifiable_variance_contract(ml_vc)
  expect_true("prediction_variance" %in% ml_vc$Component)

  structural_result <- list(
    predicted_values = gaussian_pred,
    model_parameters = data.frame(
      stat = c("requested_model", "response_family"),
      summary = c("RKHS", "gaussian"),
      stringsAsFactors = FALSE
    )
  )
  structural_vc <- PredictProR:::gp_public_variance_components(structural_result)
  expect_identical(structural_vc$Component, "model_variance_components_not_returned")
  expect_false(any(grepl("not_identifiable$", structural_vc$Component)))

  classification_pred <- make_variance_contract_classification_predictions()
  classification_result <- list(
    predicted_values = classification_pred,
    model_name = "TabTransformer",
    model_parameters = data.frame(
      stat = "response_family",
      summary = "binary",
      stringsAsFactors = FALSE
    )
  )
  classification_vc <- PredictProR:::gp_public_variance_components(classification_result)
  expect_nonidentifiable_variance_contract(classification_vc)
  expect_true("mean_prediction_confidence" %in% classification_vc$Component)
})

test_that("CV-selected and direct true-prediction fits use the same variance rules", {
  gaussian_pred <- make_variance_contract_gaussian_predictions()
  classification_pred <- make_variance_contract_classification_predictions()
  structural_gaussian_vc <- data.frame(
    Component = c("genetic_variance", "residual_variance", "heritability"),
    Components = c(0.7, 0.3, 0.7),
    Standard_error = c(0.05, 0.04, 0.03),
    stringsAsFactors = FALSE
  )
  structural_class_vc <- data.frame(
    Component = c(
      "liability_scale_genetic_variance",
      "liability_scale_residual_variance_fixed",
      "liability_scale_heritability"
    ),
    Components = c(1.5, 1, 0.6),
    Standard_error = c(0.2, 0, 0.04),
    stringsAsFactors = FALSE
  )

  cases <- list(
    list(
      model = "RandomForest", family = "gaussian",
      payload = list(predicted_values = gaussian_pred), predictive_only = TRUE
    ),
    list(
      model = "RKHS", family = "gaussian",
      payload = list(
        predicted_values = gaussian_pred,
        variance_components = structural_gaussian_vc
      ), predictive_only = FALSE
    ),
    list(
      model = "TabTransformer", family = "binary",
      payload = list(predicted_values = classification_pred), predictive_only = TRUE
    ),
    list(
      model = "GBLUP_BRR", family = "binary",
      payload = list(
        predicted_values = classification_pred,
        variance_components = structural_class_vc
      ), predictive_only = FALSE
    )
  )

  for (case in cases) {
    direct <- finalize_variance_contract_branch(
      case$model, case$payload, case$family, cv_selected = FALSE
    )
    post_cv <- finalize_variance_contract_branch(
      case$model, case$payload, case$family, cv_selected = TRUE
    )
    expect_equal(post_cv$variance_components, direct$variance_components, info = case$model)
    if (isTRUE(case$predictive_only)) {
      expect_nonidentifiable_variance_contract(direct$variance_components, info = case$model)
    } else {
      expect_false(
        any(grepl("not_identifiable$", direct$variance_components$Component)),
        info = case$model
      )
    }
  }
})

test_that("hybrid variance fallback is restricted to ML and DL", {
  pred <- data.frame(
    Hybrid_ID = paste0("H", seq_len(4L)),
    Female = paste0("F", seq_len(4L)),
    Male = paste0("M", seq_len(4L)),
    Predicted_value = c(2, 3, 4, 5),
    Train_Test_Label = c("Train", "Train", "Train", "Test"),
    Observed_value = c(2.1, 2.9, 4.2, NA_real_),
    Female_additive_contribution = c(0.2, 0.3, 0.5, 0.7),
    Male_additive_contribution = c(0.1, 0.4, 0.2, 0.6),
    Hybrid_interaction_contribution = c(0.05, 0.1, 0.3, 0.4),
    stringsAsFactors = FALSE
  )
  predictive <- PredictProR:::gp_format_hybrid_variance_components(
    list(predicted_values = pred),
    model = "RandomForest"
  )
  expect_nonidentifiable_variance_contract(predictive)
  expect_true(any(grepl("prediction_variance$", predictive$Component)))

  structural_missing <- PredictProR:::gp_format_hybrid_variance_components(
    list(predicted_values = pred),
    model = "GBLUP"
  )
  expect_identical(structural_missing$Component, "model_variance_components_not_returned")
  expect_false(any(grepl("not_identifiable$", structural_missing$Component)))
  expect_false(any(grepl("prediction_variance$", structural_missing$Component)))

  native <- data.frame(
    Component = c("female_gca_genetic_variance", "residual_variance", "heritability"),
    Components = c(0.4, 0.6, 0.4),
    Standard_error = c(0.05, 0.06, 0.03),
    stringsAsFactors = FALSE
  )
  structural <- PredictProR:::gp_format_hybrid_variance_components(
    list(predicted_values = pred, variance_components = native),
    model = "GBLUP"
  )
  expect_true(all(native$Component %in% structural$Component))
  expect_false(any(grepl("not_identifiable$", structural$Component)))
})
