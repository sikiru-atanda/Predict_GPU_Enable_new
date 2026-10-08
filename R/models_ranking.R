# library(dplyr)
#
# # Assuming cv_results_processed$aggregated_data_list$aggregated_across_reps is your dataframe:
# aggregated_data_sik <- cv_results_processed$aggregated_data_list$aggregated_across_reps


# Function to rank models for each metric
rank_models <- function(data, metric, tie_breakers = NULL) {

  msg <- ""

  # List of metrics that need to be minimized (lower is better)
  minimize_metrics <- c("mean_squared_error", "bias",
                        "root_mean_squared_error", "relative_squared_error",
                        "mean_absolute_error", "mean_absolute_percent_error",
                        "log_loss", "brier_score", "mean_absolute_error_class",
                        "ece")

  maximize_metrics <- c("accuracy", "kendalls_tau",
                        "balanced_accuracy", "precision", "recall",
                        "specificity", "f1", "mcc",
                        "macro_precision", "macro_recall", "macro_f1",
                        "within_one_class_accuracy", "quadratic_weighted_kappa")

  metric_key <- tolower(as.character(metric)[1L])
  metric_column <- names(data)[tolower(names(data)) == metric_key]
  if (length(metric_column) != 1L) {
    stop(paste0(msg, "Metric column '", metric,
                "' is not present exactly once in the CV results."),
         call. = FALSE)
  }
  metric_column <- metric_column[[1L]]
  direction <- if (exists("gp_eval_metric_direction", mode = "function")) {
    gp_eval_metric_direction(metric_key)
  } else if (metric_key %in% minimize_metrics) {
    "minimize"
  } else if (metric_key %in% maximize_metrics) {
    "maximize"
  } else {
    NA_character_
  }
  metric_values <- suppressWarnings(as.numeric(data[[metric_column]]))
  data <- data[is.finite(metric_values), , drop = FALSE]

  tie_breakers <- unique(tolower(as.character(tie_breakers %||% character())))
  tie_breakers <- setdiff(tie_breakers[!is.na(tie_breakers) & nzchar(tie_breakers)], metric_key)
  tie_columns <- vapply(tie_breakers, function(tie_metric) {
    matched <- names(data)[tolower(names(data)) == tie_metric]
    if (length(matched) == 1L) matched[[1L]] else NA_character_
  }, character(1L))
  tie_columns <- tie_columns[!is.na(tie_columns)]

  if (!direction %in% c("minimize", "maximize")) {
    stop(paste0(msg,
                "Unknown evaluation metric: '", metric, "'. Please select from the following metrics: ",
                paste(c(minimize_metrics, maximize_metrics), collapse = ", ")),
         call. = FALSE)
  }

  # Use the declared secondary order only when primary values tie. A final
  # canonical model-name key makes the result deterministic even when all
  # reported metrics are equal.
  order_keys <- list(as.character(data$trait))
  metric_order_value <- function(column) {
    values <- suppressWarnings(as.numeric(data[[column]]))
    metric_direction <- gp_eval_metric_direction(tolower(column))
    if (identical(metric_direction, "maximize")) -values else values
  }
  order_keys[[length(order_keys) + 1L]] <- metric_order_value(metric_column)
  for (tie_column in tie_columns) {
    order_keys[[length(order_keys) + 1L]] <- metric_order_value(tie_column)
  }
  order_keys[[length(order_keys) + 1L]] <- tolower(as.character(data$model))
  ord <- do.call(order, c(order_keys, list(na.last = TRUE)))
  arrange_func <- data[ord, , drop = FALSE]

  # Apply the ranking
  rankings <- arrange_func |>
    dplyr::group_by(trait) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # Filter to keep only the top-ranked model for each trait
  best_models <- rankings |>
    dplyr::filter(rank == 1) |>
    dplyr::select(
      trait,
      model,
      dplyr::any_of(c("feature_k", "feature_scoring_model", "feature_scoring_cv")),
      dplyr::all_of(unique(c(metric_column, tie_columns)))
    )

  return(best_models)
}

# metrics <- c("accuracy","mean_squared_error","bias",
#              "root_mean_squared_error")
#
# best_models_list <- lapply(metrics, function(m) rank_models(aggregated_data_sik, m))
#
# names(best_models_list) <- metrics



