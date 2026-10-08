make_classification_cv_entry <- function() {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  pred <- data.frame(
    row_id = seq_along(y),
    GID = paste0("g", seq_along(y)),
    y = y,
    yhat = c(0.1, 0.9, 0.6, 0.4),
    cv_role = "test",
    stringsAsFactors = FALSE
  )
  prob <- data.frame(
    no = c(0.9, 0.1, 0.4, 0.6),
    yes = c(0.1, 0.9, 0.6, 0.4),
    check.names = FALSE
  )
  list(
    trait = "Trait",
    rep = 1L,
    model = "ToyClassifier",
    model_canonical = "ToyClassifier",
    response_family = "binary",
    class_levels = c("no", "yes"),
    positive_class = "yes",
    eval_metrics_reps = data.frame(
      Rep = 1L, accuracy = 0.5, log_loss = 0.8573992, ece = 0.35
    ),
    results_eval_metrics_reps = data.frame(
      Rep = 1L, accuracy = 0.5, log_loss = 0.8573992, ece = 0.35
    ),
    ypred_cv_Reps_all = pred,
    yprob_cv_Reps_all = prob
  )
}

test_that("classification metric fallback and ranking directions are correct", {
  expect_identical(PredictProR:::gp_eval_metric_direction("balanced_accuracy"), "maximize")
  expect_identical(PredictProR:::gp_eval_metric_direction("macro_f1"), "maximize")
  expect_identical(PredictProR:::gp_eval_metric_direction("log_loss"), "minimize")
  expect_identical(PredictProR:::gp_eval_metric_direction("ece"), "minimize")

  expect_lt(PredictProR:::metric_worst_value("balanced_accuracy"), 0)
  expect_gt(PredictProR:::metric_worst_value("log_loss"), 0)
  expect_true(is.na(PredictProR:::safe_metric_value(
    factor(c("no", "yes")), c(NA_real_, NA_real_),
    "balanced_accuracy", response_family = "binary"
  )))

  dat <- data.frame(
    trait = "Trait",
    model = c("good", "bad"),
    balanced_accuracy = c(0.8, 0.2),
    ece = c(0.1, 0.4),
    stringsAsFactors = FALSE
  )
  expect_identical(PredictProR:::rank_models(dat, "balanced_accuracy")$model, "good")
  expect_identical(PredictProR:::rank_models(dat, "ece")$model, "good")
  dat$mcc <- NA_real_
  expect_equal(nrow(PredictProR:::rank_models(dat, "mcc")), 0L)
})

test_that("classification CV defaults are family aware and include selection metrics", {
  binary <- PredictProR:::gp_prepare_cv_metric_selection(
    response_family = "binary"
  )
  multiclass <- PredictProR:::gp_prepare_cv_metric_selection(
    response_family = "multiclass"
  )
  ordinal <- PredictProR:::gp_prepare_cv_metric_selection(
    response_family = "ordinal"
  )

  expect_identical(binary$metric_for_ranking, "balanced_accuracy")
  expect_identical(binary$ranking_tie_breakers, "log_loss")
  expect_true(all(PredictProR:::gp_eval_metric_catalog()$binary %in% binary$eval_metrics))
  expect_identical(multiclass$metric_for_ranking, "macro_f1")
  expect_identical(multiclass$ranking_tie_breakers, "log_loss")
  expect_identical(ordinal$metric_for_ranking, "quadratic_weighted_kappa")
  expect_identical(
    ordinal$ranking_tie_breakers,
    c("mean_absolute_error_class", "log_loss")
  )
  expect_identical(formals(PredictProR::model_execute)$metric_for_ranking, "auto")
})

test_that("primary ties use declared secondary metrics and then model name", {
  dat <- data.frame(
    trait = "Trait",
    model = c("ModelB", "ModelA", "ModelC"),
    balanced_accuracy = c(0.8, 0.8, 0.7),
    log_loss = c(0.25, 0.35, 0.10),
    stringsAsFactors = FALSE
  )
  best <- PredictProR:::rank_models(
    dat, "balanced_accuracy", tie_breakers = "log_loss"
  )
  expect_identical(best$model, "ModelB")
  expect_true(all(c("balanced_accuracy", "log_loss") %in% names(best)))

  dat$log_loss[1:2] <- 0.25
  expect_identical(
    PredictProR:::rank_models(
      dat, "balanced_accuracy", tie_breakers = "log_loss"
    )$model,
    "ModelA"
  )
})

test_that("explicit binary positive class controls event metrics and probability columns", {
  y <- factor(c("yes", "yes", "no", "no"), levels = c("yes", "no"))
  p_yes <- c(0.9, 0.8, 0.7, 0.1)
  expect_equal(
    PredictProR::evaluation_metrics(
      y, p_yes, "precision", response_family = "binary",
      positive_class = "yes"
    ),
    2 / 3
  )
  expect_equal(
    PredictProR::evaluation_metrics(
      y, p_yes, "recall", response_family = "binary",
      positive_class = "yes"
    ),
    1
  )
  prob <- cbind(yes = p_yes, no = 1 - p_yes)
  expect_equal(
    PredictProR::evaluation_metrics(
      y, prob, "brier_score", response_family = "binary",
      positive_class = "yes"
    ),
    mean((as.numeric(y == "yes") - p_yes)^2)
  )
  expect_error(
    PredictProR::evaluation_metrics(
      factor(c("A", "B", "C")), diag(3), "accuracy",
      response_family = "multiclass", positive_class = "A"
    ),
    "only valid"
  )
})

test_that("binary ECE calibrates predicted-class confidence", {
  y <- factor(c("no", "yes", "no", "yes"), levels = c("no", "yes"))
  p_yes <- c(0.1, 0.9, 0.6, 0.4)
  expect_equal(
    PredictProR::evaluation_metrics(y, p_yes, "ece", response_family = "binary"),
    0.35,
    tolerance = 1e-12
  )
})

test_that("classification CV exposes metrics and out-of-fold reliability", {
  entry <- make_classification_cv_entry()
  out <- PredictProR:::cv_single_loc_result_plot_process(
    list(entry),
    eval_metrics = c("accuracy", "log_loss", "ece"),
    metric_for_ranking = "accuracy",
    render_interactive = FALSE
  )

  expect_true(is.list(out$classification_metrics))
  expect_true(all(c("per_replication", "aggregated", "metric_directions",
                    "response_families", "confusion_matrices") %in%
                  names(out$classification_metrics)))
  expect_identical(
    out$classification_metrics$metric_directions$direction,
    c("maximize", "minimize", "minimize")
  )

  reliability <- out$classification_reliability
  expect_true(is.list(reliability))
  expect_true(all(c("out_of_fold_predictions", "calibration_by_bin",
                    "summary", "definition") %in% names(reliability)))
  point <- reliability$out_of_fold_predictions
  expect_true(all(c(
    "GID", "Observed_class", "Predicted_class", "Correct",
    "Prediction_confidence", "Classification_uncertainty", "Reliability",
    "Probability_no", "Probability_yes"
  ) %in% names(point)))
  expect_identical(point$GID, paste0("g", 1:4))
  expect_equal(point$Prediction_confidence, c(0.9, 0.9, 0.6, 0.6))
  expect_equal(point$Reliability, point$Prediction_confidence)
  expect_equal(point$Classification_uncertainty, 1 - point$Reliability)
  expect_equal(reliability$summary$ece, 0.35, tolerance = 1e-12)
  expect_true(is.data.frame(
    out$classification_metrics$confusion_matrices$aggregated
  ))
})

test_that("classification CV selection compares models and records the policy", {
  model_a <- make_classification_cv_entry()
  model_a$model <- "ModelA"
  model_a$eval_metrics_reps <- data.frame(
    Rep = 1L, balanced_accuracy = 0.75, log_loss = 0.45
  )
  model_a$results_eval_metrics_reps <- model_a$eval_metrics_reps
  model_b <- make_classification_cv_entry()
  model_b$model <- "ModelB"
  model_b$eval_metrics_reps <- data.frame(
    Rep = 1L, balanced_accuracy = 0.75, log_loss = 0.30
  )
  model_b$results_eval_metrics_reps <- model_b$eval_metrics_reps

  out <- PredictProR:::cv_single_loc_result_plot_process(
    list(model_a, model_b),
    eval_metrics = c("balanced_accuracy", "log_loss"),
    metric_for_ranking = "auto",
    positive_class = "yes",
    render_interactive = FALSE
  )
  selection <- out$classification_model_selection
  expect_identical(selection$policy$primary_metric, "balanced_accuracy")
  expect_identical(selection$policy$tie_breakers, "log_loss")
  expect_identical(selection$policy$positive_class, "yes")
  expect_identical(selection$models_compared$model, c("ModelA", "ModelB"))
  expect_identical(selection$selected_models$model, "ModelB")
  expect_true(selection$policy$cv_compares_models)
  expect_false(selection$policy$true_prediction_performed)
})

test_that("one-column binary probabilities are completed without inventing multiclass probabilities", {
  entry <- make_classification_cv_entry()
  entry$yprob_cv_Reps_all <- entry$yprob_cv_Reps_all["yes"]
  reliability <- PredictProR:::aggregate_cv_classification_reliability(list(entry))
  point <- reliability$out_of_fold_predictions
  expect_true(all(c("Probability_no", "Probability_yes") %in% names(point)))
  expect_equal(point$Probability_no + point$Probability_yes, rep(1, nrow(point)))
})

test_that("undefined classification metrics remain visible as NA in processed CV output", {
  entry <- make_classification_cv_entry()
  entry$eval_metrics_reps$mcc <- NA_real_
  entry$results_eval_metrics_reps$mcc <- NA_real_
  out <- PredictProR:::cv_single_loc_result_plot_process(
    list(entry), eval_metrics = c("accuracy", "mcc"),
    metric_for_ranking = "accuracy", render_interactive = FALSE
  )
  expect_true("mcc" %in% names(out$classification_metrics$aggregated))
  expect_true(is.na(out$classification_metrics$aggregated$mcc[[1L]]))
  expect_equal(nrow(out$best_models_list$mcc), 0L)
})
