multi_env_predicted_vs_observed_result_plots_process <- function(results,
                                                                 pheno_data,
                                                                 heter_groups,
                                                                 abs_very_close_threshold,
                                                                 abs_close_threshold) {
  combined_results <- data.frame()
  # Probability columns are named after each trait's class labels, which differ
  # between traits, so keep one table per model x trait instead of stacking them.
  combined_prob_results <- list()
  prob_key <- function(model, trait) paste(model, trait, sep = "\r")

  for (i in seq_along(results)) {
    current_result <- results[[i]]
    ypred_data <- as.data.frame(current_result$ypred_cv_Reps_all)
    ypred_data$trait <- current_result$trait
    ypred_data$rep <- current_result$rep
    ypred_data$model <- current_result$model
    combined_results <- rbind(combined_results, ypred_data)

    if (!is.null(current_result$yprob_cv_Reps_all) &&
        is.data.frame(current_result$yprob_cv_Reps_all) &&
        nrow(current_result$yprob_cv_Reps_all) == nrow(ypred_data)) {
      yprob_data <- as.data.frame(current_result$yprob_cv_Reps_all)
      if (!"row_id" %in% names(yprob_data) && "row_id" %in% names(ypred_data)) {
        yprob_data$row_id <- ypred_data$row_id
      }
      yprob_data$trait <- current_result$trait
      yprob_data$rep <- current_result$rep
      yprob_data$model <- current_result$model
      if (heter_groups %in% names(ypred_data) && !heter_groups %in% names(yprob_data)) {
        yprob_data[[heter_groups]] <- ypred_data[[heter_groups]]
      }
      key <- prob_key(current_result$model, current_result$trait)
      combined_prob_results[[key]] <- rbind(combined_prob_results[[key]], yprob_data)
    }
  }

  all_traits <- unique(combined_results$trait)
  mod_res_per_trait_per_model <- list()
  plot_reps_list <- list()
  classification_diagnostic_summaries <- list()

  for (trait in all_traits) {
    dat_trait <- combined_results[combined_results$trait == trait, , drop = FALSE]
    trait_models <- unique(dat_trait$model)
    plot_all_models <- list()

    for (mod in trait_models) {
      dat_mod <- dat_trait[dat_trait$model == mod, , drop = FALSE]
      envs <- unique(dat_mod[[heter_groups]])

      for (env_i in envs) {
        dat_env <- dat_mod[dat_mod[[heter_groups]] == env_i, , drop = FALSE]
        prob_env <- combined_prob_results[[prob_key(mod, trait)]]
        if (!is.null(prob_env)) {
          prob_env <- prob_env[prob_env[[heter_groups]] == env_i, , drop = FALSE]
        }

        stats <- if (is.numeric(dat_env$y) || is.integer(dat_env$y)) {
          calculate_cv_prediction_statistics(dat_env)
        } else {
          calculate_cv_classification_statistics(dat_env, prob_env)
        }

        model_env_name <- paste(mod, env_i, sep = "_")
        if (is.null(mod_res_per_trait_per_model[[model_env_name]])) {
          mod_res_per_trait_per_model[[model_env_name]] <- list()
        }
        mod_res_per_trait_per_model[[model_env_name]][[trait]] <- stats
        if (!is.numeric(dat_env$y) && !is.integer(dat_env$y)) {
          classification_diagnostic_summaries[[paste(model_env_name, trait, sep = "__")]] <-
            cv_classification_diagnostic_summaries(
              stats_df = stats,
              trait = paste(trait, env_i, sep = " - "),
              model = mod
            )
        }

        is_regression_stats <- is.list(stats) && all(c("y_mean", "pred_mean") %in% names(stats))
        plot_all_models[[model_env_name]] <- if (is_regression_stats) {
          predicted_vs_observed_ranking_plot(
            observed_value = as.numeric(stats$y_mean),
            predicted_value = as.numeric(stats$pred_mean),
            abs_very_close_threshold = abs_very_close_threshold,
            abs_close_threshold = abs_close_threshold
          )
        } else {
          cv_classification_diagnostic_plot(
            stats_df = stats,
            trait = paste(trait, env_i, sep = " - "),
            model = mod
          )
        }
      }
    }

    plot_reps_list[[trait]] <- list(
      traits = trait,
      predicted_vs_observed_plots = plot_all_models
    )
  }

  list(
    predicted_vs_observed_plots = plot_reps_list,
    mod_res_per_trait_per_model = mod_res_per_trait_per_model,
    classification_diagnostic_summaries = classification_diagnostic_summaries
  )
}
