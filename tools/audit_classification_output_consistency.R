#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "tools/audit_classification_output_consistency.R"
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
local_lib <- file.path(root, ".r-lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The classification audit requires pkgload.", call. = FALSE)
}
pkgload::load_all(root, quiet = TRUE)

fail <- function(...) stop(..., call. = FALSE)
assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
}

ordinal_aliases <- c("ordinal", "ordered", "ordered-categorical", "ordered_categorical")
nominal_aliases <- c(
  "multiclass", "multi-class", "multi_class", "multinomial", "nominal",
  "categorical", "unordered", "unordered-categorical", "unordered_categorical"
)

alias_results <- data.frame(
  input = c(ordinal_aliases, nominal_aliases),
  expected = c(rep("ordinal", length(ordinal_aliases)), rep("multiclass", length(nominal_aliases))),
  stringsAsFactors = FALSE
)
alias_results$actual <- vapply(
  alias_results$input,
  PredictProR:::gp_normalize_response_family,
  character(1L)
)
alias_results$pass <- alias_results$actual == alias_results$expected
assert_true(all(alias_results$pass), "One or more documented classification aliases normalize incorrectly.")

make_raw_prediction <- function(class_levels) {
  n_class <- length(class_levels)
  probability <- matrix(
    (1 - 0.8) / (n_class - 1),
    nrow = n_class,
    ncol = n_class,
    dimnames = list(NULL, paste0("Probability_", make.names(class_levels, unique = TRUE)))
  )
  diag(probability) <- 0.8
  data.frame(
    GID = paste0("g", seq_len(n_class)),
    Predicted_class = class_levels,
    Train_Test_Label = c(rep("Train", n_class - 1L), "Test"),
    Observed_class = c(class_levels[-n_class], NA_character_),
    probability,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

check_contract <- function(label, class_levels) {
  raw <- make_raw_prediction(class_levels)
  out <- PredictProR:::gp_public_prediction_table(raw, gen_name = "GID")
  probability_columns <- paste0("Probability_", make.names(class_levels, unique = TRUE))
  expected_columns <- c(
    PredictProR:::gp_public_prediction_columns(label),
    probability_columns
  )
  assert_true(identical(names(out), expected_columns), paste(label, "column contract differs from the canonical classification schema."))
  prob <- as.matrix(out[, probability_columns, drop = FALSE])
  assert_true(all(is.finite(prob)), paste(label, "contains non-finite class probabilities."))
  assert_true(all(prob >= 0 & prob <= 1), paste(label, "contains probabilities outside [0, 1]."))
  assert_true(max(abs(rowSums(prob) - 1)) <= 1e-12, paste(label, "probability rows do not sum to one."))
  confidence <- apply(prob, 1L, max)
  assert_true(max(abs(out$Prediction_confidence - confidence)) <= 1e-12, paste(label, "confidence is not max class probability."))
  assert_true(max(abs(out$Classification_uncertainty - (1 - confidence))) <= 1e-12, paste(label, "uncertainty is not 1 - confidence."))
  assert_true(max(abs(out$Reliability - confidence)) <= 1e-12, paste(label, "reliability is not classification confidence."))
  assert_true(all(is.finite(out$Reliability)), paste(label, "contains non-finite reliability."))
  assert_true(identical(out$Predicted_class, class_levels[max.col(prob, ties.method = "first")]), paste(label, "predicted class is not the probability argmax."))
  data.frame(
    family = label,
    rows = nrow(out),
    probability_columns = length(probability_columns),
    reliability_definition = "max_class_probability",
    pass = TRUE,
    stringsAsFactors = FALSE
  )
}

contract_results <- do.call(rbind, list(
  check_contract("ordinal", c("low", "medium", "high")),
  check_contract("nominal", c("A", "B", "C")),
  check_contract("multi-class", c("A", "B", "C"))
))

bayesian_models <- c("BayesA", "BayesB", "BayesC", "BL", "BRR", "GBLUP_BRR", "RKHS")
ml_models <- c(
  "Xgboost", "RandomForest", "CatBoost", "LightGBM",
  "SupportVectorMachine", "K-NearestNeighbors"
)
dl_models <- PredictProR:::gp_deep_learning_supported_models()

route_rows <- list()
add_route <- function(model, supplied_family, expected_family, expected_support) {
  normalized <- PredictProR:::gp_normalize_response_family(supplied_family)
  error <- tryCatch({
    PredictProR:::gp_validate_model_response_family(model, supplied_family)
    NULL
  }, error = identity)
  supported <- is.null(error)
  pass <- identical(normalized, expected_family) && identical(supported, expected_support)
  route_rows[[length(route_rows) + 1L]] <<- data.frame(
    model = model,
    supplied_family = supplied_family,
    normalized_family = normalized,
    supported = supported,
    expected_support = expected_support,
    pass = pass,
    stringsAsFactors = FALSE
  )
}
for (model in bayesian_models) {
  add_route(model, "ordinal", "ordinal", TRUE)
  add_route(model, "nominal", "multiclass", FALSE)
}
for (model in c(ml_models, dl_models)) {
  add_route(model, "multi-class", "multiclass", TRUE)
  add_route(model, "ordinal", "ordinal", FALSE)
}
route_results <- do.call(rbind, route_rows)
assert_true(all(route_results$pass), "One or more model/family routing checks failed.")

cat("\nCLASSIFICATION ALIAS AUDIT\n")
print(alias_results, row.names = FALSE)
cat("\nCLASSIFICATION OUTPUT CONTRACT AUDIT\n")
print(contract_results, row.names = FALSE)
cat("\nMODEL/FAMILY ROUTING AUDIT\n")
print(route_results, row.names = FALSE)
cat(
  "\nAUDIT SUMMARY:",
  sum(alias_results$pass), "aliases passed;",
  sum(contract_results$pass), "output contracts passed;",
  sum(route_results$pass), "model routes passed.\n"
)
cat("audit_classification_output_consistency: ok\n")
