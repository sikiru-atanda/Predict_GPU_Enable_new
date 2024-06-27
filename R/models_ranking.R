# library(dplyr)
#
# # Assuming cv_results_processed$aggregated_data_list$aggregated_across_reps is your dataframe:
# aggregated_data_sik <- cv_results_processed$aggregated_data_list$aggregated_across_reps


# Function to rank models for each metric
rank_models <- function(data, metric) {

  # List of metrics that need to be minimized (lower is better)
  minimize_metrics <- c("mean_squared_error", "bias",
                        "root_mean_squared_error", "relative_squared_error",
                        "mean_absolute_error", "mean_absolute_percent_error")

  maximize_metrics <- c("accuracy", "kendalls_tau")

  # Decide the ranking order based on the type of metric
  arrange_func <- if (metric %in% minimize_metrics) {
    dplyr::arrange(data, trait, !!dplyr::sym(metric))
  } else if(metric %in% maximize_metrics) {
    dplyr::arrange(data, trait, dplyr::desc(!!dplyr::sym(metric)))
  } else {
    stop("unknown metric")
  }

  # Apply the ranking
  rankings <- arrange_func |>
    dplyr::group_by(trait) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # Filter to keep only the top-ranked model for each trait
  best_models <- rankings |>
    dplyr::filter(rank == 1) |>
    dplyr::select(trait, model, dplyr::all_of(metric))

  return(best_models)
}

# metrics <- c("accuracy","mean_squared_error","bias",
#              "root_mean_squared_error")
#
# best_models_list <- lapply(metrics, function(m) rank_models(aggregated_data_sik, m))
#
# names(best_models_list) <- metrics



