

#' Title
#'
#' @param cv_results_data
#' @param eval_metrics
#' @param metric_for_ranking
#'
#' @return
#' @export
#'
#' @examples
cv1_cv2_and_across_env_result_plot_process <- function(cv_results_data = NULL,
                                        eval_metrics = NULL,
                                        metric_for_ranking = "accuracy") {

  # Combine all data frames into one
  combined_df <- do.call(rbind, lapply(cv_results_data, function(x) {
    cbind(trait = x$trait, rep = x$rep, model = x$model, x$eval_metrics_reps)
  }))

  # Ensure 'Rep' is treated as a factor for proper aggregation later
  combined_df$Rep <- as.factor(combined_df$rep)

  aggregated_data_list <- lapply(eval_metrics, function(metric) {
    # Aggregate across reps and environments for each metric
    agg_formula <- as.formula(paste(metric, "~ trait + model + Env"))
    aggregated_data <- aggregate(agg_formula, data = combined_df, FUN = mean, na.rm = TRUE)

    # Further aggregate to get mean across environments for each metric
    agg_formula_env <- as.formula(paste(metric, "~ trait + model"))
    aggregated_data_env_mean <- aggregate(agg_formula_env, data = aggregated_data, FUN = mean, na.rm = TRUE)

    list(metric = metric,
         aggregated_data = aggregated_data,
         aggregated_data_env_mean = aggregated_data_env_mean)
  })
  # Assign names to the elements of the list based on metrics
  names(aggregated_data_list) <- eval_metrics
  # Create plots for each metric
  plot_reps_list <- lapply(eval_metrics, function(metric) {

    ggplot_boxplot_reps <- ggplot2::ggplot(combined_df, ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~paste(trait, Env), scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Model Performance by", toupper(metric), ", Trait, and Environment"),
           x = "Model",
           y = toupper(metric))

    plotly_boxplot_reps <- plotly::ggplotly(ggplot_boxplot_reps)

    list(metric = metric,
         ggplot_boxplot_reps = ggplot_boxplot_reps,
         plotly_boxplot_reps = plotly_boxplot_reps)
  })
  names(plot_reps_list) <-  eval_metrics
  ####
  plot_mean_list <- lapply(eval_metrics, function(metric) {
    aggregated_data <- aggregated_data_list[[metric]][["aggregated_data_env_mean"]]
    ggplot_barplot_mean <-  ggplot2::ggplot(aggregated_data, ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
      ggplot2::facet_wrap(~trait, scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Mean Model Performance by Trait",
         x = "Model",
         y = "Mean Accuracy")

  plotly_barplot_mean <- plotly::ggplotly(ggplot_barplot_mean)

  list(metric = metric,
       ggplot_barplot_mean = ggplot_barplot_mean,
       plotly_barplot_mean = plotly_barplot_mean)
       })
  names(plot_reps_list) <- eval_metrics
  ########## Select best models for each trait
  aggregated_data_across_env <- combined_df |>
    dplyr::group_by(trait, model) |>
    dplyr::summarise(mean_metric = mean(!!dplyr::sym(metric_for_ranking), na.rm = TRUE), .groups = "drop")

  # Step 2: Rank models for each trait based on the aggregated metric
  rankings <- aggregated_data_across_env |>
    dplyr::group_by(trait) |>
    dplyr::arrange(trait, dplyr::desc(mean_metric)) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()

  # Filter to keep only the top-ranked model for each trait
  best_models <- rankings |>
    dplyr::filter(rank == 1) |>
    dplyr::select(trait, model, mean_metric)
  colnames(best_models)[3] <- metric_for_ranking
  # Return the list of plots
  return(list(aggregated_data_list = aggregated_data_list,
              plot_reps_list = plot_reps_list,
              plot_mean_list = plot_mean_list,
              best_models = best_models))
}

# Example usage

#results <- cv1_cv2_and_across_env_result_plot_process(cv_results_data=TT, eval_metrics = eval_metrics)


