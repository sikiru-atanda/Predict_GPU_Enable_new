
calculate_cv_prediction_statistics <- function(dattt) {
  if (!"row_id" %in% names(dattt)) {
    dattt$row_id <- seq_len(nrow(dattt))
  }

  dattt <- dattt[is.finite(dattt$yhat) & is.finite(dattt$y), , drop = FALSE]
  if (nrow(dattt) == 0L) {
    return(list(pred_mean = numeric(0), y_mean = numeric(0)))
  }

  aggregated <- dattt |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      pred_mean = mean(yhat, na.rm = TRUE),
      y_mean = mean(y, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id)

  list(
    pred_mean = aggregated$pred_mean,
    y_mean = aggregated$y_mean
  )
}

calculate_cv_classification_statistics <- function(dattt, prob_dat = NULL) {
  if (!"row_id" %in% names(dattt)) {
    dattt$row_id <- seq_len(nrow(dattt))
  }

  dattt <- dattt[!is.na(dattt$y), , drop = FALSE]
  if (!nrow(dattt)) {
    return(NULL)
  }

  prob_cols <- character(0)
  if (!is.null(prob_dat) && is.data.frame(prob_dat)) {
    if (!"row_id" %in% names(prob_dat)) {
      prob_dat$row_id <- seq_len(nrow(prob_dat))
    }
    prob_cols <- grep("^Prob_", names(prob_dat), value = TRUE)
  }

  label_mode <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    tbl <- sort(table(as.character(x)), decreasing = TRUE)
    names(tbl)[1]
  }

  aggregated_labels <- dattt |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      observed = label_mode(y),
      predicted = label_mode(yhat),
      cv_role = label_mode(cv_role),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id)

  if (!length(prob_cols)) {
    return(aggregated_labels)
  }

  aggregated_prob <- prob_dat |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(prob_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(row_id)

  dplyr::left_join(aggregated_labels, aggregated_prob, by = "row_id")
}

cv_classification_diagnostic_summaries <- function(stats_df, trait, model) {
  if (is.null(stats_df) || !nrow(stats_df)) {
    return(NULL)
  }

  prob_cols <- grep("^Prob_", names(stats_df), value = TRUE)
  if (!length(prob_cols)) {
    return(NULL)
  }

  class_levels <- gsub("^Prob_", "", prob_cols)
  prob_mat <- as.matrix(stats_df[, prob_cols, drop = FALSE])
  colnames(prob_mat) <- class_levels
  pred_class <- class_levels[max.col(prob_mat, ties.method = "first")]
  confidence <- apply(prob_mat, 1, max, na.rm = TRUE)
  observed <- as.character(stats_df$observed)
  cv_role <- stats_df$cv_role %||% NA_character_

  confusion <- data.frame(
    trait = trait,
    model = model,
    cv_role = cv_role,
    Observed = factor(observed, levels = class_levels),
    Predicted = factor(pred_class, levels = class_levels),
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$trait, .data$model, .data$cv_role, .data$Observed, .data$Predicted) |>
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop")

  if (length(class_levels) == 2L) {
    positive_class <- class_levels[2]
    calibration_summary <- data.frame(
      trait = trait,
      model = model,
      cv_role = cv_role,
      observed = observed == positive_class,
      confidence = prob_mat[, 2],
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        bin = cut(.data$confidence, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
      ) |>
      dplyr::group_by(.data$trait, .data$model, .data$cv_role, .data$bin) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_confidence = mean(.data$confidence, na.rm = TRUE),
        observed_rate = mean(.data$observed, na.rm = TRUE),
        .groups = "drop"
      )

    confidence_summary <- data.frame(
      trait = trait,
      model = model,
      cv_role = cv_role,
      Correct = ifelse(pred_class == observed, "Correct", "Incorrect"),
      Confidence = confidence,
      stringsAsFactors = FALSE
    ) |>
      dplyr::group_by(.data$trait, .data$model, .data$cv_role, .data$Correct) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_confidence = mean(.data$Confidence, na.rm = TRUE),
        median_confidence = stats::median(.data$Confidence, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    calibration_summary <- NULL
    confidence_summary <- data.frame(
      trait = trait,
      model = model,
      cv_role = cv_role,
      Predicted = factor(pred_class, levels = class_levels),
      Confidence = confidence,
      stringsAsFactors = FALSE
    ) |>
      dplyr::group_by(.data$trait, .data$model, .data$cv_role, .data$Predicted) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_confidence = mean(.data$Confidence, na.rm = TRUE),
        median_confidence = stats::median(.data$Confidence, na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    confusion_matrix = confusion,
    calibration_summary = calibration_summary,
    confidence_summary = confidence_summary
  )
}

cv_classification_diagnostic_plot <- function(stats_df, trait, model) {
  prob_cols <- grep("^Prob_", names(stats_df), value = TRUE)
  if (!length(prob_cols)) {
    return(NULL)
  }

  class_levels <- gsub("^Prob_", "", prob_cols)
  prob_mat <- as.matrix(stats_df[, prob_cols, drop = FALSE])
  colnames(prob_mat) <- class_levels
  pred_class <- class_levels[max.col(prob_mat, ties.method = "first")]
  confidence <- apply(prob_mat, 1, max, na.rm = TRUE)
  observed <- as.character(stats_df$observed)

  if (length(class_levels) == 2L) {
    positive_class <- class_levels[2]
    p_pos <- prob_mat[, 2]
    calib_df <- data.frame(
      observed = observed == positive_class,
      prob = p_pos
    ) |>
      dplyr::mutate(bin = cut(prob, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) |>
      dplyr::group_by(bin) |>
      dplyr::summarise(
        mean_prob = mean(prob, na.rm = TRUE),
        observed_rate = mean(observed, na.rm = TRUE),
        .groups = "drop"
      )

    p1 <- ggplot2::ggplot(calib_df, ggplot2::aes(x = mean_prob, y = observed_rate)) +
      ggplot2::geom_point(color = "steelblue", size = 3) +
      ggplot2::geom_line(color = "steelblue") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("Binary Calibration:", trait, "-", model),
        x = paste("Mean predicted P(", positive_class, ")", sep = ""),
        y = paste("Observed rate of", positive_class)
      ) +
      ggplot2::theme_minimal()

    conf_df <- data.frame(
      Confidence = confidence,
      Correct = factor(ifelse(pred_class == observed, "Correct", "Incorrect"),
                       levels = c("Correct", "Incorrect"))
    )
    p2 <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Confidence, fill = Correct)) +
      ggplot2::geom_histogram(position = "identity", alpha = 0.6, bins = 10) +
      ggplot2::labs(
        title = paste("Prediction Confidence:", trait, "-", model),
        x = "Predicted class confidence",
        y = "Count"
      ) +
      ggplot2::scale_fill_manual(values = c("Correct" = "forestgreen", "Incorrect" = "firebrick")) +
      ggplot2::theme_minimal()
  } else {
    cm_df <- as.data.frame(table(
      Observed = factor(observed, levels = class_levels),
      Predicted = factor(pred_class, levels = class_levels)
    ), stringsAsFactors = FALSE)

    p1 <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = Observed, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = Freq), color = "white") +
      ggplot2::scale_fill_gradient(low = "grey90", high = "steelblue") +
      ggplot2::labs(
        title = paste("Confusion Matrix:", trait, "-", model),
        x = "Predicted class",
        y = "Observed class"
      ) +
      ggplot2::theme_minimal()

    conf_df <- data.frame(
      Predicted = factor(pred_class, levels = class_levels),
      Confidence = confidence
    )
    p2 <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Predicted, y = Confidence, fill = Predicted)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("Class Confidence:", trait, "-", model),
        x = "Predicted class",
        y = "Predicted class confidence"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }

  gp_arrange_grob_safely(p1, p2, ncol = 2)
}

single_predicted_vs_observed_result_plots_process <- function(results,
                                                              pheno_data,
                                                              abs_very_close_threshold ,
                                                              abs_close_threshold){
#browser()
  combined_results <- data.frame()
  # Probability columns are named after each trait's class labels, which differ
  # between traits, so keep one table per model x trait instead of stacking them.
  combined_prob_results <- list()
  prob_key <- function(model, trait) paste(model, trait, sep = "\r")

  # Loop through each element in the results list
  for (i in seq_along(results)) {
    # Extract the current element
    current_result <- results[[i]]

    # Extract the ypred_cv_Reps_all data frame
    ypred_data <- as.data.frame(current_result$ypred_cv_Reps_all)

    # Add the trait, rep, and model columns
    ypred_data$trait <- current_result$trait
    ypred_data$rep <- current_result$rep
    ypred_data$model <- current_result$model

    # Combine the current ypred_data with the combined_results data frame
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
      key <- prob_key(current_result$model, current_result$trait)
      combined_prob_results[[key]] <- rbind(combined_prob_results[[key]], yprob_data)
    }
  }

  models <- unique(combined_results$model)
  all_traits <- unique(combined_results$trait)

  # Initialize lists to store results
  mod_res_per_trait_per_model <- list()
  classification_diagnostic_summaries <- list()

  # Iterate through each model
  for (mod in models) {
    # Filter data for the current model
    dat <- combined_results[combined_results[["model"]] == mod, ]

    traits <- unique(dat$trait)
    # Initialize a list to store results per trait for the current model
    mod_res_per_trait <- list()

    # Iterate through each trait
    for (tt in traits) {
      # Filter data for the current trait and extract 'yhat'
      dattt <- dat[dat[["trait"]] == tt, ]
      prob_dattt <- combined_prob_results[[prob_key(mod, tt)]]
      stats <- if (is.numeric(dattt$y) || is.integer(dattt$y)) {
        calculate_cv_prediction_statistics(dattt)
      } else {
        calculate_cv_classification_statistics(dattt, prob_dattt)
      }


      # Store the statistics in the list for the current trait
      mod_res_per_trait[[tt]] <- stats
      if (!is.numeric(dattt$y) && !is.integer(dattt$y)) {
        classification_diagnostic_summaries[[paste(mod, tt, sep = "__")]] <-
          cv_classification_diagnostic_summaries(stats, trait = tt, model = mod)
      }
    }

    # Store the trait results list in the list for the current model
    mod_res_per_trait_per_model[[mod]] <- mod_res_per_trait
  }

  #return(mod_res_per_trait_per_model)
  plot_reps_list <- lapply(all_traits, function(trait) {

    plot_all_models <- list()
    trait_models <- names(mod_res_per_trait_per_model)[vapply(
      mod_res_per_trait_per_model,
      function(x) trait %in% names(x),
      logical(1)
    )]

    for (mod in trait_models) {

      prob_dattt <- combined_prob_results[[prob_key(mod, trait)]]

      trait_stats <- mod_res_per_trait_per_model[[mod]][[trait]]
      is_regression_stats <- is.list(trait_stats) &&
        all(c("y_mean", "pred_mean") %in% names(trait_stats))

      if (is_regression_stats) {
        sik <-  predicted_vs_observed_ranking_plot(
          observed_value = as.numeric(trait_stats$y_mean),
          predicted_value = as.numeric(trait_stats$pred_mean),
          abs_very_close_threshold = abs_very_close_threshold,
          abs_close_threshold =abs_close_threshold
        )
      } else {
        dattt <- combined_results[
          combined_results[["trait"]] == trait & combined_results[["model"]] == mod,
          ,
          drop = FALSE
        ]
        class_stats <- calculate_cv_classification_statistics(dattt, prob_dattt)
        sik <- cv_classification_diagnostic_plot(class_stats, trait = trait, model = mod)
      }


      plot_all_models[[mod]] <- sik
    }


    list(traits = trait,
         predicted_vs_observed_plots = plot_all_models
    )
  })

  names(plot_reps_list) <-  all_traits

  return(list(predicted_vs_observed_plots = plot_reps_list,
              mod_res_per_trait_per_model = mod_res_per_trait_per_model,
              classification_diagnostic_summaries = classification_diagnostic_summaries))

}




