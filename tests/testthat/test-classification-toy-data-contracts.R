mk_toy_classification_data <- function(n = 18L, p = 6L) {
  set.seed(20260517)
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  rownames(X) <- paste0("g", seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))

  binary <- factor(ifelse(X[, 1] - 0.4 * X[, 2] > 0, "yes", "no"), levels = c("no", "yes"))
  ordinal_score <- X[, 1] + X[, 2]
  ordinal <- ordered(
    cut(
      ordinal_score,
      breaks = stats::quantile(ordinal_score, probs = c(0, 1 / 3, 2 / 3, 1)),
      include.lowest = TRUE,
      labels = c("low", "medium", "high")
    ),
    levels = c("low", "medium", "high")
  )
  multiclass_score <- cbind(A = X[, 1], B = X[, 2], C = -X[, 1] - X[, 2])
  multiclass <- factor(colnames(multiclass_score)[max.col(multiclass_score)], levels = c("A", "B", "C"))

  binary[c(4, 15)] <- NA
  ordinal[c(5, 16)] <- NA
  multiclass[c(6, 17)] <- NA

  list(
    pheno = data.frame(
      GID = rownames(X),
      BinaryTrait = binary,
      OrdinalTrait = ordinal,
      MultiTrait = multiclass,
      stringsAsFactors = FALSE
    ),
    geno = X
  )
}

expect_classification_public_contract <- function(out, probability_columns) {
  expect_identical(
    names(out),
    c(
      "GID",
      "Predicted_class",
      "Train_Test_Label",
      "Observed_class",
      "Prediction_confidence",
      "Classification_uncertainty",
      "Reliability",
      "Reliability_remarks",
      probability_columns
    )
  )
  expect_false(any(c("Predicted_value", "Standard_error", "PEV") %in% names(out)))
  expect_true(all(is.finite(out$Prediction_confidence)))
  expect_true(all(is.finite(out$Classification_uncertainty)))
  expect_true(all(is.finite(out$Reliability)))
  expect_true(all(!is.na(out$Reliability_remarks) & nzchar(out$Reliability_remarks)))

  probability_matrix <- as.matrix(out[, probability_columns, drop = FALSE])
  expect_true(all(is.finite(probability_matrix)))
  expect_true(all(probability_matrix >= 0 & probability_matrix <= 1))
  expect_equal(rowSums(probability_matrix), rep(1, nrow(out)), tolerance = 1e-12)
  expected_confidence <- apply(probability_matrix, 1L, max)
  expect_equal(out$Prediction_confidence, expected_confidence, tolerance = 1e-12)
  expect_equal(out$Classification_uncertainty, 1 - expected_confidence, tolerance = 1e-12)
  expect_equal(out$Reliability, expected_confidence, tolerance = 1e-12)

  predicted_probability_column <- vapply(
    out$Predicted_class,
    PredictProR:::gp_classification_probability_name,
    character(1L)
  )
  predicted_probability <- vapply(seq_len(nrow(out)), function(i) {
    probability_matrix[i, match(predicted_probability_column[[i]], probability_columns)]
  }, numeric(1L))
  expect_equal(predicted_probability, expected_confidence, tolerance = 1e-12)
}

test_that("toy classification data resolves binary, ordinal, and multiclass families", {
  toy <- mk_toy_classification_data()

  expect_identical(
    PredictProR:::gp_resolve_response_family_for_responses(toy$pheno, "BinaryTrait", "auto"),
    "binary"
  )
  expect_identical(
    PredictProR:::gp_resolve_response_family_for_responses(toy$pheno, "OrdinalTrait", "auto"),
    "ordinal"
  )
  expect_identical(
    PredictProR:::gp_resolve_response_family_for_responses(toy$pheno, "MultiTrait", "auto"),
    "multiclass"
  )
  expect_identical(PredictProR:::gp_response_class_levels(toy$pheno$BinaryTrait, "binary"), c("no", "yes"))
  expect_identical(PredictProR:::gp_response_class_levels(toy$pheno$OrdinalTrait, "ordinal"), c("low", "medium", "high"))
  expect_identical(PredictProR:::gp_response_class_levels(toy$pheno$MultiTrait, "multiclass"), c("A", "B", "C"))
  expect_error(
    PredictProR:::gp_resolve_response_family_for_responses(
      toy$pheno,
      c("BinaryTrait", "OrdinalTrait"),
      "auto"
    ),
    "Multiple response families were inferred"
  )
})

test_that("documented ordinal and nominal aliases normalize consistently", {
  ordinal_aliases <- c("ordinal", "ordered", "ordered-categorical", "ordered_categorical")
  nominal_aliases <- c(
    "multiclass", "multi-class", "multi_class", "multinomial", "nominal",
    "categorical", "unordered", "unordered-categorical", "unordered_categorical"
  )

  expect_identical(
    vapply(ordinal_aliases, PredictProR:::gp_normalize_response_family, character(1L)),
    stats::setNames(rep("ordinal", length(ordinal_aliases)), ordinal_aliases)
  )
  expect_identical(
    vapply(nominal_aliases, PredictProR:::gp_normalize_response_family, character(1L)),
    stats::setNames(rep("multiclass", length(nominal_aliases)), nominal_aliases)
  )
  expect_identical(
    PredictProR:::gp_public_prediction_columns("nominal"),
    PredictProR:::gp_public_prediction_columns("multiclass")
  )
  expect_identical(
    PredictProR:::gp_public_prediction_columns("multi-class", include_env = TRUE, include_trait = TRUE),
    c(
      "GID", "Env", "Trait", "Predicted_class", "Train_Test_Label",
      "Observed_class", "Prediction_confidence", "Classification_uncertainty",
      "Reliability", "Reliability_remarks"
    )
  )
})

test_that("toy binary prediction tables use the classification public contract", {
  toy <- mk_toy_classification_data()
  predicted_class <- ifelse(seq_len(nrow(toy$pheno)) %% 2L == 0L, "yes", "no")
  probability_yes <- ifelse(predicted_class == "yes", 0.8, 0.2)
  raw <- data.frame(
    GID = toy$pheno$GID,
    Predicted_value = predicted_class,
    Observed_value = as.character(toy$pheno$BinaryTrait),
    Train_Test_Label = ifelse(is.na(toy$pheno$BinaryTrait), "Test", "Train"),
    Prediction_confidence = rep(0.8, nrow(toy$pheno)),
    Classification_uncertainty = rep(0.2, nrow(toy$pheno)),
    Prob_no = 1 - probability_yes,
    Prob_yes = probability_yes,
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")

  expect_classification_public_contract(out, c("Probability_no", "Probability_yes"))
  expect_identical(out$Predicted_class, raw$Predicted_value)
  expect_identical(out$Observed_class, raw$Observed_value)
  expect_equal(sum(out$Train_Test_Label == "Test"), sum(is.na(toy$pheno$BinaryTrait)))
})

test_that("toy ordinal prediction tables retain ordered class labels and probabilities", {
  toy <- mk_toy_classification_data()
  predicted_class <- rep(c("low", "medium", "high"), length.out = nrow(toy$pheno))
  probability_matrix <- t(vapply(predicted_class, function(level) {
    out <- stats::setNames(rep(0.15, 3L), c("low", "medium", "high"))
    out[[level]] <- 0.7
    out
  }, numeric(3L)))
  raw <- data.frame(
    GID = toy$pheno$GID,
    Predicted_class = predicted_class,
    Observed_class = as.character(toy$pheno$OrdinalTrait),
    Train_Test_Label = ifelse(is.na(toy$pheno$OrdinalTrait), "Test", "Train"),
    Prediction_confidence = rep(0.7, nrow(toy$pheno)),
    Classification_uncertainty = rep(0.3, nrow(toy$pheno)),
    Probability_low = probability_matrix[, "low"],
    Probability_medium = probability_matrix[, "medium"],
    Probability_high = probability_matrix[, "high"],
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")

  expect_classification_public_contract(
    out,
    c("Probability_low", "Probability_medium", "Probability_high")
  )
  expect_equal(
    PredictProR::evaluation_metrics(
      y_observed = ordered(c("low", "medium", "high"), levels = c("low", "medium", "high")),
      y_predicted = c("low", "high", "medium"),
      eval_metrics = "mean_absolute_error_class",
      response_family = "ordinal"
    ),
    2 / 3
  )
})

test_that("toy multiclass prediction tables retain all class probability columns", {
  toy <- mk_toy_classification_data()
  predicted_class <- rep(c("A", "B", "C"), length.out = nrow(toy$pheno))
  probability_matrix <- t(vapply(predicted_class, function(level) {
    out <- stats::setNames(rep(0.125, 3L), c("A", "B", "C"))
    out[[level]] <- 0.75
    out
  }, numeric(3L)))
  raw <- data.frame(
    GID = toy$pheno$GID,
    Predicted_value = predicted_class,
    Observed_value = as.character(toy$pheno$MultiTrait),
    Train_Test_Label = ifelse(is.na(toy$pheno$MultiTrait), "Test", "Train"),
    Prediction_confidence = rep(0.75, nrow(toy$pheno)),
    Classification_uncertainty = rep(0.25, nrow(toy$pheno)),
    Prob_A = probability_matrix[, "A"],
    Prob_B = probability_matrix[, "B"],
    Prob_C = probability_matrix[, "C"],
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")

  expect_classification_public_contract(out, c("Probability_A", "Probability_B", "Probability_C"))
  expect_identical(out$Predicted_class, raw$Predicted_value)
  expect_equal(sum(out$Train_Test_Label == "Test"), sum(is.na(toy$pheno$MultiTrait)))
})

test_that("classification confidence, uncertainty, and reliability are reconstructed from probabilities", {
  raw <- data.frame(
    GID = c("g1", "g2", "g3"),
    Predicted_class = c("A", "B", "C"),
    Train_Test_Label = c("Train", "Train", "Test"),
    Observed_class = c("A", "B", NA),
    Probability_A = c(0.8, 0.1, 0.1),
    Probability_B = c(0.1, 0.7, 0.2),
    Probability_C = c(0.1, 0.2, 0.7),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")
  expect_classification_public_contract(
    out,
    c("Probability_A", "Probability_B", "Probability_C")
  )
  expect_equal(out$Reliability, c(0.8, 0.7, 0.7))
})

test_that("multiclass probabilities are not fabricated from confidence alone", {
  raw <- data.frame(
    GID = c("g1", "g2", "g3"),
    Predicted_class = c("A", "B", "C"),
    Observed_class = c("A", "B", "C"),
    Train_Test_Label = rep("Train", 3L),
    Prediction_confidence = c(0.8, 0.7, 0.6),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")
  expect_length(PredictProR:::gp_classification_probability_columns(out), 0L)
  expect_equal(out$Reliability, raw$Prediction_confidence)
  expect_equal(out$Classification_uncertainty, 1 - raw$Prediction_confidence)
})

test_that("classification observed classes are restored from phenotype by GID", {
  pred <- data.frame(
    GID = c("g3", "g1", "g2"),
    Predicted_class = c("yes", "no", "yes"),
    Train_Test_Label = c("Test", "Train", "Train"),
    Observed_class = NA_character_,
    Prediction_confidence = c(0.6, 0.8, 0.7),
    Probability_no = c(0.4, 0.8, 0.3),
    Probability_yes = c(0.6, 0.2, 0.7),
    stringsAsFactors = FALSE
  )
  pheno <- data.frame(
    GID = c("g1", "g2", "g3"),
    Trait = factor(c("no", "yes", NA), levels = c("no", "yes")),
    stringsAsFactors = FALSE
  )

  out <- PredictProR:::gp_merge_classification_observed_values(
    pred,
    pheno_data = pheno,
    response = "Trait",
    gen_name = "GID"
  )

  expect_identical(out$Observed_class, c(NA_character_, "no", "yes"))
})

test_that("toy classification model guardrails route only compatible model families", {
  toy <- mk_toy_classification_data()
  expect_silent(PredictProR:::gp_validate_model_response_family("RandomForest", "binary"))
  expect_silent(PredictProR:::gp_validate_model_response_family("RandomForest", "multiclass"))
  expect_silent(PredictProR:::gp_validate_model_response_family("BayesA", "ordinal"))
  expect_silent(PredictProR:::gp_validate_model_response_family("RKHS", "binary"))
  expect_error(
    PredictProR:::gp_validate_model_response_family("KRR", "binary"),
    "does not support response_family='binary'"
  )
  expect_error(
    PredictProR:::gp_validate_model_response_family("CatBoost", "ordinal"),
    "does not support response_family='ordinal'"
  )
  expect_error(
    PredictProR:::gp_resolve_response_family_for_responses(
      toy$pheno,
      c("BinaryTrait", "MultiTrait"),
      "auto"
    ),
    "Multiple response families were inferred"
  )
})
