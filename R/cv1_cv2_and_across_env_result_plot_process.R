

#' Title
#'
#' @param cv_results_data
#' @param eval_metrics
#' @param heter_groups Environment-grouping column used for the MET summary.
#' @param metric_for_ranking
#' @param ranking_tie_breakers Optional ordered secondary metrics used to
#'   resolve exact primary-metric ties.
#' @param positive_class Optional event class for binary classification
#'   summaries.
#'
#' @return
#' @export
#'
#' @examples
cv1_cv2_and_across_env_result_plot_process <- function(cv_results_data = NULL,
                                        eval_metrics = NULL,
                                        heter_groups = NULL,
                                        metric_for_ranking = "auto",
                                        ranking_tie_breakers = NULL,
                                        positive_class = NULL) {

  if (is.null(heter_groups)) stop(paste(heter_groups,  "is not defined. This is required for multi-environment GS."))
  gaussian_uncertainty_summaries <- aggregate_cv_gaussian_uncertainty_summaries(cv_results_data)
  gaussian_risk_summaries <- aggregate_cv_gaussian_risk_summaries(cv_results_data)
  # Combine all data frames into one
  combined_df <- do.call(rbind, lapply(cv_results_data, function(x) {
    cbind(trait = x$trait, rep = x$rep, model = x$model, x$eval_metrics_reps)
  }))

  # combined_df <- combined_df |>
  #   dplyr::mutate(dplyr::across(-c(trait, model, !!sym(heter_groups)), as.numeric))

  feature_group_cols <- intersect(
    c("feature_k", "feature_scoring_model", "feature_scoring_cv"),
    names(combined_df)
  )
  columns_to_exclude <- c("trait", "model", heter_groups, feature_group_cols)
  columns_to_convert <- setdiff(names(combined_df), columns_to_exclude)

  combined_df <- combined_df |>
    dplyr::mutate(dplyr::across(dplyr::all_of(columns_to_convert), as.numeric))

  # Ensure 'Rep' is treated as a factor for proper aggregation later
  combined_df$Rep <- as.factor(combined_df$rep)
  environment_value <- as.character(combined_df[[heter_groups]])
  pooled_df <- combined_df[environment_value == "ALL", , drop = FALSE]
  environment_df <- combined_df[environment_value != "ALL", , drop = FALSE]
  if (!nrow(environment_df)) {
    environment_df <- combined_df
  }
  pooled_group_cols <- c("trait", "model", feature_group_cols)
  classification_metric_base <- list(
    combined_df = combined_df,
    aggregated_across_reps = unique(combined_df[, pooled_group_cols, drop = FALSE])
  )

  aggregated_data_list <- lapply(eval_metrics, function(metric) {
    # Environment rows remain descriptive. The pre-computed ALL rows contain
    # the pooled out-of-fold metric and are intentionally excluded here so
    # they are not counted as an additional environment.
    env_group_cols <- c(pooled_group_cols, heter_groups)
    aggregated_data <- gp_cv_metric_mean_by(environment_df, env_group_cols, metric)

    # This equal-environment mean is a descriptive summary only. Model
    # selection below is based on pooled ALL rows across replications.
    aggregated_data_env_mean <- gp_cv_metric_mean_by(
      aggregated_data, pooled_group_cols, metric
    )

    list(metric = metric,
         aggregated_data = aggregated_data,
         aggregated_data_env_mean = aggregated_data_env_mean)
  })
  # Assign names to the elements of the list based on metrics
  names(aggregated_data_list) <- eval_metrics
  pooled_metric_frames <- lapply(eval_metrics, function(metric) {
    source_df <- if (nrow(pooled_df)) pooled_df else environment_df
    gp_cv_metric_mean_by(source_df, pooled_group_cols, metric)
  })
  classification_metric_base$aggregated_across_reps <- gp_cv_merge_metric_frames(
    pooled_metric_frames,
    pooled_group_cols
  )
  classification_metrics <- aggregate_cv_classification_metrics(
    cv_results_data, classification_metric_base, eval_metrics,
    heter_groups = heter_groups
  )
  classification_reliability <- aggregate_cv_classification_reliability(
    cv_results_data, heter_groups = heter_groups
  )
  classification_confusion_matrices <-
    aggregate_cv_classification_confusion_matrices(
      classification_reliability, heter_groups = heter_groups
    )
  if (!is.null(classification_metrics)) {
    classification_metrics$confusion_matrices <- classification_confusion_matrices
  }
  response_family <- unique(vapply(cv_results_data, function(x) {
    as.character(x$response_family %||% "gaussian")
  }, character(1L)))[[1L]]
  metric_for_ranking <- gp_resolve_metric_for_ranking(
    metric_for_ranking, response_family
  )
  ranking_tie_breakers <- gp_resolve_ranking_tie_breakers(
    ranking_tie_breakers, response_family, metric_for_ranking
  )
  # Create plots for each metric
  plot_reps_list <- lapply(eval_metrics, function(metric) {
    metric_sym <- rlang::sym(metric)

    ggplot_boxplot_reps <- ggplot2::ggplot(
      environment_df,
      ggplot2::aes(x = .data$model, y = !!metric_sym, fill = .data$model)
    ) +
      ggplot2::geom_boxplot() +
      #ggplot2::facet_wrap(~paste(trait, Env), scales = "free") +
      ggplot2::facet_wrap(as.formula(paste("~ trait +", heter_groups)), scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Model Performance by", toupper(metric), ", Trait, and Environment"),
           x = "Model",
           y = toupper(metric))

    #plotly_boxplot_reps <- plotly::ggplotly(ggplot_boxplot_reps)

    list(metric = metric,
         ggplot_boxplot_reps = ggplot_boxplot_reps
         #plotly_boxplot_reps = plotly_boxplot_reps
         )
  })
  names(plot_reps_list) <-  eval_metrics
  ####
  plot_mean_list <- lapply(eval_metrics, function(metric) {
    aggregated_data <- aggregated_data_list[[metric]][["aggregated_data_env_mean"]]
    metric_sym <- rlang::sym(metric)
    # ggplot_barplot_mean <-  ggplot2::ggplot(aggregated_data, ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
    #   ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    #   ggplot2::facet_wrap(~trait, scales = "free_x") +
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    #   ggplot2::labs(title = "Mean Model Performance by Trait",
    #      x = "Model",
    #      y = "Mean Accuracy")
    ####
    ggplot_lineplot_mean <- ggplot2::ggplot(
      aggregated_data,
      ggplot2::aes(x = .data$model, y = !!metric_sym, group = .data$trait, color = .data$trait)
    ) +
      ggplot2::geom_point(size = 3) + # Add points
      ggplot2::geom_line(ggplot2::aes(linetype = .data$trait), linewidth = 1) + # Add lines
      ggplot2::scale_color_discrete() + # Dynamic color assignment without optional viridisLite
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Comparison of", metric, "Across Models for Different Traits"),
                    x = "Model",
                    y = metric) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
                     legend.position = "right", # Move the legend to the bottom
                     panel.grid.major.x = ggplot2::element_blank(), # Turn off x-axis major grid lines
                     panel.grid.minor.x = ggplot2::element_blank(), # Turn off x-axis minor grid lines
                     panel.grid.major.y = ggplot2::element_blank(), # Add grey lines for y-axis major grid lines
                     panel.grid.minor.y = ggplot2::element_blank(), # Turn off y-axis minor grid lines
                     axis.line = ggplot2::element_line(colour = "black"), # Add black line for x and y axis
                     axis.ticks = ggplot2::element_line(color = "black"), # Change axis tick color if needed
                     axis.ticks.length = ggplot2::unit(-0.1, "cm")) # Set axis ticks to point inward

  #plotly_lineplot_mean <- plotly::ggplotly(ggplot_lineplot_mean)

  list(metric = metric,
       ggplot_lineplot_mean = ggplot_lineplot_mean
       #plotly_lineplot_mean = plotly_lineplot_mean
       )
       })
  names(plot_mean_list) <- eval_metrics
  ########## Select best models for each trait
  # aggregated_data_across_env <- combined_df |>
  #   dplyr::group_by(trait, model) |>
  #   dplyr::summarise(mean_metric = mean(!!rlang::sym(metric_for_ranking), na.rm = TRUE), .groups = "drop")

  # Step 2: Rank models for each trait based on the aggregated metric
  # rankings <- aggregated_data_across_env |>
  #   dplyr::group_by(trait) |>
  #   dplyr::arrange(trait, dplyr::desc(mean_metric)) |>
  #   dplyr::mutate(rank = dplyr::row_number()) |>
  #   dplyr::ungroup()
  #
  # # Filter to keep only the top-ranked model for each trait
  # best_models <- rankings |>
  #   dplyr::filter(rank == 1) |>
  #   dplyr::select(trait, model, mean_metric)
  # colnames(best_models)[3] <- metric_for_ranking
  #
  best_models_list <- lapply(eval_metrics, function(m) {
    rank_models(
      classification_metric_base[["aggregated_across_reps"]],
      m,
      tie_breakers = if (identical(m, metric_for_ranking)) ranking_tie_breakers else NULL
    )
  })

  names(best_models_list) <- eval_metrics

  # Return the list of plots
  return(list(aggregated_data_list = aggregated_data_list,
              plot_reps_list = plot_reps_list,
              plot_mean_list = plot_mean_list,
              best_models_list = best_models_list,
              classification_model_selection = gp_classification_cv_selection_summary(
                classification_metrics = classification_metrics,
                best_models_list = best_models_list,
                response_family = response_family,
                metric_for_ranking = metric_for_ranking,
                ranking_tie_breakers = ranking_tie_breakers,
                positive_class = positive_class,
                cv_results_data = cv_results_data,
                heter_groups = heter_groups,
                is_multi = TRUE
              ),
              classification_probability_summaries = aggregate_cv_probability_summaries(
                cv_results_data, heter_groups = heter_groups
              ),
              classification_metrics = classification_metrics,
              classification_reliability = classification_reliability,
              gaussian_uncertainty_summaries = gaussian_uncertainty_summaries,
              gaussian_risk_summaries = gaussian_risk_summaries,
              gaussian_risk_plots = build_cv_gaussian_risk_plots(gaussian_risk_summaries)))
}

# Example usage

#results <- cv1_cv2_and_across_env_result_plot_process(cv_results_data=TT, eval_metrics = eval_metrics)


