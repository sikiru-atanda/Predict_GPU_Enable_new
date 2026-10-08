.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all(".")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

model_routes <- list(
  hybrid_asreml = PredictProR:::gp_hybrid_asreml_supported_models(),
  hybrid_bayes = PredictProR:::gp_hybrid_bayes_supported_models(),
  hybrid_gp = PredictProR:::gp_hybrid_gp_supported_models(),
  hybrid_ml = PredictProR:::gp_hybrid_ml_supported_models(),
  hybrid_dl = PredictProR:::gp_hybrid_dl_supported_models()
)
route_audit <- do.call(rbind, lapply(names(model_routes), function(mode) {
  do.call(rbind, lapply(model_routes[[mode]], function(model) {
    PredictProR:::gp_validate_hybrid_model_selection(
      mode = mode,
      GS_model = model,
      cross_validation = FALSE
    )
    data.frame(mode = mode, model = model, gaussian_only = TRUE, route_pass = TRUE)
  }))
}))
rownames(route_audit) <- NULL
assert_true(nrow(route_audit) == 15L, "Expected 15 supported hybrid engine/model routes.")

shape_audit <- do.call(rbind, lapply(
  list(
    single_trait = c(FALSE, FALSE),
    multi_environment = c(TRUE, FALSE),
    multi_trait = c(FALSE, TRUE),
    multi_trait_multi_environment = c(TRUE, TRUE)
  ),
  function(flags) {
    cols <- PredictProR:::gp_hybrid_prediction_columns(
      include_env = flags[[1L]],
      include_trait = flags[[2L]]
    )
    data.frame(
      include_env = flags[[1L]],
      include_trait = flags[[2L]],
      columns = length(cols),
      exact_unique_order = identical(cols, unique(cols))
    )
  }
))
shape_audit$shape <- rownames(shape_audit)
rownames(shape_audit) <- NULL
shape_audit <- shape_audit[, c("shape", "include_env", "include_trait", "columns", "exact_unique_order")]
assert_true(all(shape_audit$exact_unique_order), "Hybrid column contract contains duplicates.")

classification_error <- tryCatch(
  {
    prediction_output_standard(response_family = "classification", task = "hybrid")
    NULL
  },
  error = identity
)
assert_true(
  inherits(classification_error, "error") &&
    grepl("supports gaussian responses only", conditionMessage(classification_error), fixed = TRUE),
  "Hybrid classification requests must fail with the gaussian-only scope message."
)

model_pred <- data.frame(
  HybridID = paste0("H", 1:4),
  Female = c("F1", "F1", "F2", "F2"),
  Male = c("M1", "M2", "M1", "M2"),
  Predicted_value = c(1.0, 1.4, 1.8, 2.1),
  Train_Test_Label = c("Train", "Train", "Train", "Test"),
  Observed_value = c(1.1, 1.3, 1.9, NA_real_),
  Standard_error = c(0.20, 0.25, 0.30, 0.35),
  PEV = c(0.04, 0.0625, 0.09, 0.1225),
  Reliability_reference_variance = rep(0.50, 4L),
  stringsAsFactors = FALSE
)
model_result <- PredictProR:::gp_standardize_hybrid_model_result(list(
  bayes_model = list(fitted = TRUE),
  predicted_values = model_pred,
  variance_components = data.frame(
    Component = c("female_gca_variance", "male_gca_variance", "sca_variance", "residual_variance"),
    Components = c(0.20, 0.15, 0.15, 0.30),
    Standard_error = c(0.02, 0.02, 0.01, 0.03)
  )
))
model_public <- model_result$predicted_values
assert_true(
  identical(names(model_public), PredictProR:::gp_hybrid_prediction_columns()),
  "Model-implied hybrid prediction schema mismatch."
)
assert_true(all(is.finite(model_public$PEV)), "Model-implied hybrid PEV must be finite.")
assert_true(
  isTRUE(all.equal(model_public$Standard_error^2, model_public$PEV, tolerance = 1e-12)),
  "Model-implied hybrid SE squared must equal PEV."
)
assert_true(
  isTRUE(all.equal(
    model_public$Reliability,
    pmax(0, pmin(1, 1 - model_public$PEV / model_public$Reliability_reference_variance)),
    tolerance = 1e-12
  )),
  "Model-implied hybrid reliability formula mismatch."
)

predictive_pred <- model_pred[, c(
  "HybridID", "Female", "Male", "Predicted_value", "Train_Test_Label", "Observed_value"
)]
predictive_pred$Female_additive_contribution <- c(-0.2, -0.2, 0.3, 0.3)
predictive_pred$Male_additive_contribution <- c(-0.1, 0.1, -0.1, 0.1)
predictive_pred$Hybrid_interaction_contribution <- c(0.0, 0.1, 0.2, -0.1)
predictive_result <- PredictProR:::gp_standardize_hybrid_model_result(list(
  predicted_values = predictive_pred
))
predictive_public <- predictive_result$predicted_values
assert_true(
  identical(names(predictive_public), PredictProR:::gp_hybrid_prediction_columns()),
  "Predictive hybrid ML/DL schema mismatch."
)
assert_true(all(is.na(predictive_public$PEV)), "Uncalibrated hybrid ML/DL true prediction must not fabricate PEV.")
assert_true(all(is.finite(predictive_public$Reliability)), "Hybrid ML/DL plotting reliability must be finite.")
assert_true(
  isTRUE(all.equal(
    predictive_public$Reliability,
    pmax(0, pmin(
      1,
      1 - predictive_public$Reliability_variance_input /
        predictive_public$Reliability_reference_variance
    )),
    tolerance = 1e-12
  )),
  "Hybrid ML/DL marker-adjustment reliability formula mismatch."
)
assert_true(
  all(grepl("not genetic reliability", predictive_public$Reliability_basis, fixed = TRUE)),
  "Hybrid ML/DL reliability basis must state that it is not genetic reliability."
)

cat("HYBRID MODEL ROUTE AUDIT\n")
print(route_audit, row.names = FALSE)
cat("\nHYBRID OUTPUT SHAPE AUDIT\n")
print(shape_audit, row.names = FALSE)
cat("\nHYBRID STATISTICAL CONTRACT AUDIT\n")
print(data.frame(
  contract = c(
    "model_implied_PEV_and_reliability",
    "ML_DL_plotting_reliability_without_fabricated_PEV",
    "classification_rejected"
  ),
  pass = TRUE
), row.names = FALSE)
cat("\nAUDIT SUMMARY: 15 model routes passed; 4 output shapes passed; 3 statistical contracts passed.\n")
cat("audit_hybrid_prediction_consistency: ok\n")
