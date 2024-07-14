
aggregate_metrics <- function(cv_results_data) {
  # Convert the list to a data frame
  combined_df <- do.call(rbind, lapply(cv_results_data, function(x) {
    cbind(data.frame(trait = x$trait, model = x$model, rep = x$rep), x$eval_metrics_reps)
  }))

  combined_df <- combined_df |>
    dplyr::mutate(dplyr::across(-c(trait, model), as.numeric))
  # Calculate mean evaluation metrics across reps for each trait and model
  aggregated_across_reps <- stats::aggregate(. ~ trait + model, data = combined_df, FUN = mean)

  return(list(aggregated_across_reps = aggregated_across_reps,
              combined_df = combined_df))
}

# plot_model_metric_by_trait <- function(data, metric = "accuracy") {
#   # Validate metric
#   if (!metric %in% names(data)) {
#     stop("The specified metric is not found in the data.")
#   }
#
#   # Construct the plot
#   ggplot_boxplot_reps <- ggplot2::ggplot(data, ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
#     ggplot2::geom_violin(trim = FALSE) +
#     #ggplot2::geom_boxplot() +
#     ggplot2::facet_wrap(~ trait, scales = "free") +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
#     ggplot2::labs(title = paste("Boxplot of Model", metric, "by Trait"),
#                   x = "Model",
#                   y = metric)
#
#   plotly_boxplot_reps <- plotly::ggplotly(ggplot_boxplot_reps)
#
#   return(list(ggplot_boxplot_reps = ggplot_boxplot_reps,
#               plotly_boxplot_reps = plotly_boxplot_reps))
#   #return(plotly::ggplotly(p))
# }



cv_single_loc_result_plot_process <- function(cv_results_data = NULL,
                                              eval_metrics = NULL,
                                              metric_for_ranking = "accuracy") {


  aggregated_data_list <- aggregate_metrics(cv_results_data)

  plot_reps_list <- lapply(eval_metrics, function(metric) {

    ggplot_boxplot_reps <- ggplot2::ggplot(aggregated_data_list[["combined_df"]], ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
      #ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ trait, scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Boxplot of Model", metric, "by Trait"),
                    x = "Model",
                    y = metric)

    plotly_boxplot_reps <- plotly::ggplotly(ggplot_boxplot_reps)

    list(metric = metric,
         ggplot_boxplot_reps = ggplot_boxplot_reps,
         plotly_boxplot_reps = plotly_boxplot_reps)
  })
  names(plot_reps_list) <-  eval_metrics
  ###################
  plot_mean_list <- lapply(eval_metrics, function(metric) {

    # ggplot_barplot_mean <-  ggplot2::ggplot(aggregated_data_list[["aggregated_across_reps"]], ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
    #   ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    #   ggplot2::facet_wrap(~trait, scales = "free_x") +
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    #   ggplot2::labs(title = "Mean Model Performance by Trait",
    #                 x = "Model",
    #                 y = "Mean Accuracy")
    ######
    ggplot_lineplot_mean <-   ggplot2::ggplot(aggregated_data_list[["aggregated_across_reps"]], ggplot2::aes_string(x = "model", y = metric, group = "trait", color = "trait")) +
      ggplot2::geom_point(size = 3) + # Add points
      ggplot2::geom_line(ggplot2::aes_string(linetype = "trait"), size = 1) + # Add lines
      ggplot2::scale_color_viridis_d() + # Dynamic color assignment
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Comparison of Prediction Accuracy Across Models for Different Traits",
                    x = "Model",
                    y = "Accuracy") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
                     legend.position = "right", # Move the legend to the bottom
                     panel.grid.major.x = ggplot2::element_blank(), # Turn off x-axis major grid lines
                     panel.grid.minor.x = ggplot2::element_blank(), # Turn off x-axis minor grid lines
                     panel.grid.major.y = ggplot2::element_blank(), # Add grey lines for y-axis major grid lines
                     panel.grid.minor.y = ggplot2::element_blank(), # Turn off y-axis minor grid lines
                     axis.line = ggplot2::element_line(colour = "black"), # Add black line for x and y axis
                     axis.ticks = ggplot2::element_line(color = "black"), # Change axis tick color if needed
                     axis.ticks.length = ggplot2::unit(-0.1, "cm")) # Set axis ticks to point inward

    plotly_lineplot_mean <- plotly::ggplotly(ggplot_lineplot_mean)


    plotly_lineplot_mean <- plotly::ggplotly(ggplot_lineplot_mean)

    list(metric = metric,
         ggplot_lineplot_mean = ggplot_lineplot_mean,
         plotly_lineplot_mean = plotly_lineplot_mean)
  })
  names(plot_mean_list) <- eval_metrics
  ##### Initiate process for ranking of models per traits
  # aggregated_data_across_reps <- aggregated_data_list[["combined_df"]] |>
  #   dplyr::group_by(trait, model) |>
  #   dplyr::summarise(mean_metric = mean(!!rlang::sym(metric_for_ranking), na.rm = TRUE), .groups = "drop")

  # Step 2: Rank models for each trait based on the aggregated metric
 #  rankings <- aggregated_data_across_reps |>
 #    dplyr::group_by(trait) |>
 #    dplyr::arrange(trait, dplyr::desc(mean_metric)) |>
 #    dplyr::mutate(rank = dplyr::row_number()) |>
 #    dplyr::ungroup()
 #
 #  # Filter to keep only the top-ranked model for each trait
 #  best_models <- rankings |>
 #    dplyr::filter(rank == 1) |>
 #    dplyr::select(trait, model, mean_metric)
 #
 # colnames(best_models)[3] <- metric_for_ranking
  best_models_list <- lapply(eval_metrics, function(metric) rank_models(data = aggregated_data_list[["aggregated_across_reps"]], metric = metric))

  names(best_models_list) <- eval_metrics


 return(list(aggregated_data_list = aggregated_data_list,
             plot_reps_list = plot_reps_list,
             plot_mean_list = plot_mean_list,
             best_models_list = best_models_list))
}


 #MM = cv_single_loc_result_plot_process(TT, eval_metrics)

