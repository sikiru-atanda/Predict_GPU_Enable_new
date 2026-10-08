# Class-probability columns are named after each trait's class labels. With
# two separately fitted traits that do not share labels, stacking every
# trait's probability table into one data frame failed with
# "names do not match previous names" after all CV fits had finished.
one_cv_result <- function(trait, model, labels, env = NULL) {
  n <- 12L
  y <- sample(labels, n, replace = TRUE)
  prob <- matrix(stats::runif(n * length(labels)), n, dimnames = list(NULL, labels))
  prob <- prob / rowSums(prob)
  pred <- data.frame(row_id = seq_len(n), GID = sprintf("L%02d", seq_len(n)),
                     y = y, yhat = sample(labels, n, replace = TRUE),
                     cv_role = "test", stringsAsFactors = FALSE)
  if (!is.null(env)) pred$Env <- rep(env, length.out = n)
  list(trait = trait, rep = 1L, model = model,
       ypred_cv_Reps_all = pred, yprob_cv_Reps_all = as.data.frame(prob))
}
two_label_sets <- function(env = NULL) {
  set.seed(7)
  list(
    one_cv_result("Grain_colour", "RandomForest", c("red", "white", "amber"), env),
    one_cv_result("Plant_type", "RandomForest", c("erect", "semi", "prostrate"), env),
    one_cv_result("Grain_colour", "Xgboost", c("red", "white", "amber"), env),
    one_cv_result("Plant_type", "Xgboost", c("erect", "semi", "prostrate"), env)
  )
}

test_that("CV diagnostics combine classification traits with different class labels", {
  results <- two_label_sets()
  out <- NULL
  expect_no_error(
    out <- suppressWarnings(PredictProR:::single_predicted_vs_observed_result_plots_process(
      results = results, pheno_data = NULL,
      abs_very_close_threshold = 0.1, abs_close_threshold = 0.2
    ))
  )
  expect_false(is.null(out))
})

test_that("MET CV diagnostics combine classification traits with different class labels", {
  results <- two_label_sets(env = c("E1", "E2"))
  out <- NULL
  expect_no_error(
    out <- suppressWarnings(PredictProR:::multi_env_predicted_vs_observed_result_plots_process(
      results = results, pheno_data = NULL, heter_groups = "Env",
      abs_very_close_threshold = 0.1, abs_close_threshold = 0.2
    ))
  )
  expect_setequal(names(out$predicted_vs_observed_plots), c("Grain_colour", "Plant_type"))
})

test_that("CV classification reliability combines traits with different class labels", {
  results <- lapply(two_label_sets(), function(r) {
    r$response_family <- "multiclass"
    r$class_levels <- colnames(r$yprob_cv_Reps_all)
    r
  })
  out <- NULL
  expect_no_error(out <- PredictProR:::aggregate_cv_classification_reliability(results))
  expect_setequal(unique(out$out_of_fold_predictions$trait), c("Grain_colour", "Plant_type"))
  expect_equal(nrow(out$summary), 4L)
  expect_true(all(is.finite(out$summary$ece)))
})

test_that("CV probability summaries keep every trait's class columns", {
  out <- NULL
  expect_no_error(out <- PredictProR:::aggregate_cv_probability_summaries(two_label_sets()))
  expect_setequal(out$probability_columns,
                  c("red", "white", "amber", "erect", "semi", "prostrate"))
  # a trait's own classes are filled; the other trait's classes are NA
  gc_rows <- out$combined_probabilities$trait == "Grain_colour"
  expect_false(anyNA(out$combined_probabilities[gc_rows, c("red", "white", "amber")]))
  expect_true(all(is.na(out$combined_probabilities[gc_rows, c("erect", "semi", "prostrate")])))
  expect_equal(nrow(out$aggregated_probabilities), 4L * 12L)
  agg <- out$aggregated_probabilities
  expect_false(any(vapply(agg[out$probability_columns], function(v) any(is.nan(v)), logical(1))))
})
