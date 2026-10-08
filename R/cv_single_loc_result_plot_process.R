
aggregate_metrics <- function(cv_results_data) {
  # Convert the list to a data frame
  combined_df <- do.call(rbind, lapply(cv_results_data, function(x) {
    cbind(data.frame(trait = x$trait, model = x$model, rep = x$rep), x$eval_metrics_reps)
  }))

  feature_group_cols <- intersect(
    c("feature_k", "feature_scoring_model", "feature_scoring_cv"),
    names(combined_df)
  )
  group_cols <- c("trait", "model", feature_group_cols)
  metric_cols <- setdiff(names(combined_df), c(group_cols, "rep"))
  combined_df[metric_cols] <- lapply(combined_df[metric_cols], function(col) {
    suppressWarnings(as.numeric(col))
  })
  if (!length(metric_cols)) {
    return(list(
      aggregated_across_reps = combined_df[0, group_cols, drop = FALSE],
      combined_df = combined_df
    ))
  }

  # Calculate mean evaluation metrics across reps for each trait and model
  aggregated_across_reps <- stats::aggregate(
    combined_df[metric_cols],
    by = combined_df[group_cols],
    FUN = function(x) {
      out <- mean(x, na.rm = TRUE)
      if (is.nan(out)) NA_real_ else out
    }
  )

  return(list(aggregated_across_reps = aggregated_across_reps,
              combined_df = combined_df))
}

# Stack per trait x model CV tables whose class-probability columns are named
# after that trait's class labels. Traits fitted in one run can have different
# labels, so take the union of columns and fill NA for classes a trait lacks.
gp_rbind_class_probability_entries <- function(entries) {
  all_cols <- unique(unlist(lapply(entries, names), use.names = FALSE))
  entries <- lapply(entries, function(entry) {
    for (col in setdiff(all_cols, names(entry))) entry[[col]] <- NA_real_
    entry[all_cols]
  })
  do.call(rbind, entries)
}

aggregate_cv_probability_summaries <- function(cv_results_data,
                                               heter_groups = NULL) {
  prob_entries <- lapply(cv_results_data, function(x) {
    prob_df <- x$yprob_cv_Reps_all
    pred_df <- x$ypred_cv_Reps_all
    if (is.null(prob_df) || !is.data.frame(prob_df) || nrow(prob_df) == 0L) {
      return(NULL)
    }
    prob_df <- as.data.frame(prob_df, check.names = FALSE, stringsAsFactors = FALSE)
    prob_cols_current <- names(prob_df)[vapply(prob_df, is.numeric, logical(1))]
    prob_cols_current <- setdiff(
      prob_cols_current,
      c("trait", "model", "rep", "row_id", "y", "yhat", "cv_role", "Env", "env", heter_groups)
    )

    row_id <- if (!is.null(pred_df) && "row_id" %in% names(pred_df)) pred_df$row_id else seq_len(nrow(prob_df))
    cv_role <- if (!is.null(pred_df) && "cv_role" %in% names(pred_df)) pred_df$cv_role else NA_character_
    observed_y <- if (!is.null(pred_df) && "y" %in% names(pred_df)) pred_df$y else NA
    pred_yhat <- if (!is.null(pred_df) && "yhat" %in% names(pred_df)) pred_df$yhat else NA
    extra_cols <- NULL
    if (!is.null(pred_df)) {
      extra_names <- setdiff(names(pred_df), c("row_id", "y", "yhat", "cv_role"))
      if (length(extra_names)) {
        extra_cols <- pred_df[, extra_names, drop = FALSE]
      }
    }

    out <- data.frame(
      trait = x$trait,
      model = x$model,
      rep = x$rep,
      row_id = row_id,
      y = observed_y,
      yhat = pred_yhat,
      cv_role = cv_role,
      prob_df,
      stringsAsFactors = FALSE
    )
    if (!is.null(extra_cols)) {
      out <- cbind(out, extra_cols)
    }
    attr(out, "probability_columns") <- prob_cols_current
    out
  })

  prob_entries <- Filter(Negate(is.null), prob_entries)
  if (!length(prob_entries)) {
    return(NULL)
  }

  prob_cols <- unique(unlist(
    lapply(prob_entries, function(entry) attr(entry, "probability_columns", exact = TRUE)),
    use.names = FALSE
  ))
  combined_probabilities <- gp_rbind_class_probability_entries(prob_entries)
  prob_cols <- intersect(prob_cols, names(combined_probabilities))
  prob_cols <- prob_cols[vapply(combined_probabilities[prob_cols], is.numeric, logical(1))]

  if (!length(prob_cols)) {
    return(NULL)
  }

  group_cols <- c("trait", "model", "row_id")
  extra_group_cols <- intersect(
    setdiff(names(combined_probabilities), c(group_cols, "rep", prob_cols)),
    c("y", "yhat", "cv_role", "Env", "env", heter_groups)
  )
  grouped <- combined_probabilities |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, extra_group_cols)))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(prob_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$trait, .data$model, .data$row_id)
  grouped[prob_cols] <- lapply(grouped[prob_cols], function(v) replace(v, is.nan(v), NA_real_))

  list(
    combined_probabilities = combined_probabilities,
    aggregated_probabilities = grouped,
    probability_columns = prob_cols
  )
}

aggregate_cv_classification_metrics <- function(cv_results_data,
                                                aggregated_data_list,
                                                eval_metrics,
                                                heter_groups = NULL) {
  has_classification <- any(vapply(cv_results_data, function(x) {
    fam <- as.character(x$response_family %||% "")
    nzchar(fam) && !identical(fam, "gaussian") ||
      is.data.frame(x$yprob_cv_Reps_all)
  }, logical(1L)))
  if (!has_classification) return(NULL)

  metric_cols <- intersect(as.character(eval_metrics),
                           names(aggregated_data_list$combined_df))
  id_cols_rep <- intersect(
    c("trait", "model", "rep", "Rep", heter_groups, "feature_k",
      "feature_scoring_model", "feature_scoring_cv"),
    names(aggregated_data_list$combined_df)
  )
  id_cols_agg <- intersect(
    c("trait", "model", "feature_k", "feature_scoring_model",
      "feature_scoring_cv"),
    names(aggregated_data_list$aggregated_across_reps)
  )
  family_rows <- do.call(rbind, lapply(cv_results_data, function(x) {
    pred_y <- x$ypred_cv_Reps_all$y %||% NULL
    fam <- as.character(x$response_family %||% if (is.ordered(pred_y)) {
      "ordinal"
    } else if (length(x$class_levels %||% character()) == 2L) {
      "binary"
    } else {
      "multiclass"
    })
    data.frame(trait = as.character(x$trait), model = as.character(x$model),
               response_family = fam,
               positive_class = as.character(x$positive_class %||% NA_character_),
               stringsAsFactors = FALSE)
  }))

  list(
    per_replication = aggregated_data_list$combined_df[
      , unique(c(id_cols_rep, metric_cols)), drop = FALSE
    ],
    aggregated = aggregated_data_list$aggregated_across_reps[
      , unique(c(id_cols_agg, metric_cols)), drop = FALSE
    ],
    metric_directions = data.frame(
      metric = metric_cols,
      direction = vapply(metric_cols, gp_eval_metric_direction, character(1L)),
      stringsAsFactors = FALSE
    ),
    response_families = unique(family_rows)
  )
}

aggregate_cv_classification_confusion_matrices <- function(classification_reliability,
                                                           heter_groups = NULL) {
  predictions <- classification_reliability$out_of_fold_predictions %||% NULL
  if (!is.data.frame(predictions) || !nrow(predictions)) {
    return(NULL)
  }
  group_count <- function(group_cols) {
    keys <- interaction(predictions[group_cols], drop = TRUE, lex.order = TRUE)
    rows <- lapply(split(seq_len(nrow(predictions)), keys), function(idx) {
      out <- predictions[idx[[1L]], group_cols, drop = FALSE]
      out$n_predictions <- length(idx)
      out
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out[do.call(order, out[group_cols]), , drop = FALSE]
  }
  base_cols <- c(
    "trait", "model", "response_family", "Observed_class", "Predicted_class"
  )
  out <- list(
    per_replication = group_count(c(base_cols[1:3], "rep", base_cols[4:5])),
    aggregated = group_count(base_cols)
  )
  if (!is.null(heter_groups) && length(heter_groups) == 1L &&
      heter_groups %in% names(predictions)) {
    out$per_environment_per_replication <- group_count(
      c(base_cols[1:3], "rep", heter_groups, base_cols[4:5])
    )
    out$per_environment_aggregated <- group_count(
      c(base_cols[1:3], heter_groups, base_cols[4:5])
    )
  }
  out
}

aggregate_cv_classification_reliability <- function(cv_results_data,
                                                    heter_groups = NULL,
                                                    bins = 10L,
                                                    high_confidence = 0.9,
                                                    low_confidence = 0.5) {
  bins <- max(2L, as.integer(bins)[1L])
  entries <- lapply(cv_results_data, function(x) {
    pred <- x$ypred_cv_Reps_all
    prob <- x$yprob_cv_Reps_all
    if (!is.data.frame(pred) || !is.data.frame(prob) ||
        nrow(pred) != nrow(prob) || !nrow(pred)) return(NULL)

    fam <- as.character(x$response_family %||% if (is.ordered(pred$y)) {
      "ordinal"
    } else if (length(x$class_levels %||% character()) == 2L) {
      "binary"
    } else {
      "multiclass"
    })
    if (identical(fam, "gaussian")) return(NULL)
    class_levels <- as.character(x$class_levels %||%
      if (is.factor(pred$y)) levels(pred$y) else character())

    prob <- as.data.frame(prob, check.names = FALSE, stringsAsFactors = FALSE)
    prob_mat <- suppressWarnings(as.matrix(data.frame(
      lapply(prob, as.numeric), check.names = FALSE
    )))
    storage.mode(prob_mat) <- "double"
    prob_labels <- sub("^(Probability|Prob)_", "", colnames(prob_mat))
    if (!length(class_levels)) class_levels <- prob_labels

    if (identical(fam, "binary") && ncol(prob_mat) == 1L &&
        length(class_levels) == 2L) {
      positive <- prob_mat[, 1L]
      prob_mat <- cbind(1 - positive, positive)
      colnames(prob_mat) <- class_levels
    } else {
      matched <- match(make.names(class_levels), make.names(prob_labels))
      if (length(matched) == length(class_levels) && !anyNA(matched)) {
        prob_mat <- prob_mat[, matched, drop = FALSE]
        colnames(prob_mat) <- class_levels
      } else if (ncol(prob_mat) == length(class_levels)) {
        colnames(prob_mat) <- class_levels
      } else {
        return(NULL)
      }
    }

    observed <- as.character(pred$y)
    role <- as.character(pred$cv_role %||% rep("test", nrow(pred)))
    row_sum <- rowSums(prob_mat)
    keep <- role == "test" & !is.na(observed) &
      rowSums(is.finite(prob_mat)) == ncol(prob_mat) &
      is.finite(row_sum) & row_sum > 0
    if (!any(keep)) return(NULL)
    prob_mat <- prob_mat[keep, , drop = FALSE]
    prob_mat <- prob_mat / rowSums(prob_mat)
    pred_idx <- max.col(prob_mat, ties.method = "first")
    predicted <- colnames(prob_mat)[pred_idx]
    confidence <- prob_mat[cbind(seq_len(nrow(prob_mat)), pred_idx)]
    correct <- observed[keep] == predicted
    bin_id <- cut(confidence, breaks = seq(0, 1, length.out = bins + 1L),
                  include.lowest = TRUE, labels = FALSE)
    remarks <- ifelse(confidence >= high_confidence, "High Confidence",
                      ifelse(confidence <= low_confidence, "Low Confidence",
                             "Moderate Confidence"))

    probability_names <- make.unique(vapply(
      colnames(prob_mat), gp_classification_probability_name, character(1L)
    ))
    colnames(prob_mat) <- probability_names
    out <- data.frame(
      trait = as.character(x$trait),
      model = as.character(x$model),
      response_family = fam,
      rep = as.integer(x$rep),
      row_id = pred$row_id[keep] %||% which(keep),
      GID = as.character(pred$GID[keep] %||% pred$row_id[keep] %||% which(keep)),
      Observed_class = observed[keep],
      Predicted_class = predicted,
      Correct = as.logical(correct),
      Prediction_confidence = confidence,
      Classification_uncertainty = 1 - confidence,
      Reliability = confidence,
      Reliability_remarks = remarks,
      Calibration_bin = as.integer(bin_id),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!is.null(heter_groups) && length(heter_groups) == 1L &&
        heter_groups %in% names(pred)) {
      out[[heter_groups]] <- as.character(pred[[heter_groups]][keep])
    }
    cbind(out, as.data.frame(prob_mat, check.names = FALSE))
  })
  entries <- Filter(Negate(is.null), entries)
  if (!length(entries)) return(NULL)
  predictions <- gp_rbind_class_probability_entries(entries)
  rownames(predictions) <- NULL

  group_summary <- function(df, include_bin = FALSE) {
    group_cols <- c("trait", "model", if (include_bin) "Calibration_bin")
    key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
    rows <- lapply(split(seq_len(nrow(df)), key), function(idx) {
      z <- df[idx, , drop = FALSE]
      base <- z[1L, group_cols, drop = FALSE]
      base$n_predictions <- nrow(z)
      base$mean_confidence <- mean(z$Prediction_confidence)
      base$empirical_accuracy <- mean(z$Correct)
      base$mean_classification_uncertainty <- mean(z$Classification_uncertainty)
      base$calibration_gap <- base$empirical_accuracy - base$mean_confidence
      base$absolute_calibration_gap <- abs(base$calibration_gap)
      base
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  calibration <- group_summary(predictions, include_bin = TRUE)
  summary <- group_summary(predictions, include_bin = FALSE)
  summary$ece <- vapply(seq_len(nrow(summary)), function(i) {
    z <- calibration[
      calibration$trait == summary$trait[[i]] &
        calibration$model == summary$model[[i]], , drop = FALSE
    ]
    sum((z$n_predictions / sum(z$n_predictions)) * z$absolute_calibration_gap)
  }, numeric(1L))

  out <- list(
    out_of_fold_predictions = predictions,
    calibration_by_bin = calibration,
    summary = summary,
    definition = data.frame(
      quantity = c("Reliability", "Classification_uncertainty", "ece"),
      definition = c(
        "Maximum predicted class probability for the out-of-fold prediction; not PEV or genetic reliability.",
        "One minus maximum predicted class probability.",
        "Frequency-weighted absolute gap between mean confidence and empirical accuracy across calibration bins."
      ),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(heter_groups) && length(heter_groups) == 1L &&
      heter_groups %in% names(predictions)) {
    group_summary_environment <- function(df, include_bin = FALSE) {
      group_cols <- c(
        "trait", "model", heter_groups,
        if (include_bin) "Calibration_bin"
      )
      key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
      rows <- lapply(split(seq_len(nrow(df)), key), function(idx) {
        z <- df[idx, , drop = FALSE]
        base <- z[1L, group_cols, drop = FALSE]
        base$n_predictions <- nrow(z)
        base$mean_confidence <- mean(z$Prediction_confidence)
        base$empirical_accuracy <- mean(z$Correct)
        base$mean_classification_uncertainty <- mean(z$Classification_uncertainty)
        base$calibration_gap <- base$empirical_accuracy - base$mean_confidence
        base$absolute_calibration_gap <- abs(base$calibration_gap)
        base
      })
      env_out <- do.call(rbind, rows)
      rownames(env_out) <- NULL
      env_out
    }
    calibration_env <- group_summary_environment(predictions, include_bin = TRUE)
    summary_env <- group_summary_environment(predictions, include_bin = FALSE)
    summary_env$ece <- vapply(seq_len(nrow(summary_env)), function(i) {
      z <- calibration_env[
        calibration_env$trait == summary_env$trait[[i]] &
          calibration_env$model == summary_env$model[[i]] &
          calibration_env[[heter_groups]] == summary_env[[heter_groups]][[i]],
        , drop = FALSE
      ]
      sum((z$n_predictions / sum(z$n_predictions)) * z$absolute_calibration_gap)
    }, numeric(1L))
    out$calibration_by_environment_bin <- calibration_env
    out$summary_by_environment <- summary_env
  }
  out
}

aggregate_cv_gaussian_uncertainty_summaries <- function(cv_results_data) {
  gaussian_entries <- lapply(cv_results_data, function(x) {
    pred_df <- x$ypred_cv_Reps_all
    if (is.null(pred_df) || !is.data.frame(pred_df) || !all(c("y", "yhat") %in% names(pred_df))) {
      return(NULL)
    }

    y_num <- suppressWarnings(as.numeric(pred_df$y))
    yhat_num <- suppressWarnings(as.numeric(pred_df$yhat))
    ok_numeric <- is.finite(y_num) & is.finite(yhat_num)
    if (!any(ok_numeric)) {
      return(NULL)
    }

    row_id <- if ("row_id" %in% names(pred_df)) pred_df$row_id else seq_len(nrow(pred_df))
    cv_role <- if ("cv_role" %in% names(pred_df)) pred_df$cv_role else NA_character_
    env_col <- intersect(names(pred_df), c("Env", "env"))
    env_val <- if (length(env_col)) pred_df[[env_col[[1]]]] else NULL

    out <- data.frame(
      trait = x$trait,
      model = x$model,
      rep = x$rep,
      row_id = row_id,
      y = y_num,
      yhat = yhat_num,
      cv_role = cv_role,
      stringsAsFactors = FALSE
    )
    if (!is.null(env_val)) {
      out$Env <- env_val
    }

    out$residual <- out$yhat - out$y
    out$absolute_error <- abs(out$residual)
    out$squared_error <- out$residual^2

    if ("Standard_error" %in% names(pred_df)) {
      out$Standard_error <- suppressWarnings(as.numeric(pred_df$Standard_error))
    }
    if ("PEV" %in% names(pred_df)) {
      out$PEV <- suppressWarnings(as.numeric(pred_df$PEV))
    }
    if (all(c("lower_bound", "upper_bound") %in% names(pred_df))) {
      out$lower_bound <- suppressWarnings(as.numeric(pred_df$lower_bound))
      out$upper_bound <- suppressWarnings(as.numeric(pred_df$upper_bound))
      out$interval_width <- out$upper_bound - out$lower_bound
      out$covered <- out$y >= out$lower_bound & out$y <= out$upper_bound
    }

    out
  })

  gaussian_entries <- Filter(Negate(is.null), gaussian_entries)
  if (!length(gaussian_entries)) {
    return(NULL)
  }

  combined_uncertainty <- do.call(rbind, gaussian_entries)
  if ("cv_role" %in% names(combined_uncertainty) && any(!is.na(combined_uncertainty$cv_role))) {
    combined_uncertainty <- combined_uncertainty[combined_uncertainty$cv_role %in% c("test", NA), , drop = FALSE]
  }
  if (!nrow(combined_uncertainty)) {
    return(NULL)
  }

  group_cols <- c("trait", "model")
  if ("Env" %in% names(combined_uncertainty)) {
    group_cols <- c(group_cols, "Env")
  }

  aggregated_uncertainty <- combined_uncertainty |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_observed = mean(.data$y, na.rm = TRUE),
      mean_predicted = mean(.data$yhat, na.rm = TRUE),
      mean_absolute_error = mean(.data$absolute_error, na.rm = TRUE),
      mean_squared_error = mean(.data$squared_error, na.rm = TRUE),
      root_mean_squared_error = sqrt(mean(.data$squared_error, na.rm = TRUE)),
      residual_sd = stats::sd(.data$residual, na.rm = TRUE),
      mean_standard_error = if ("Standard_error" %in% names(combined_uncertainty)) mean(.data$Standard_error, na.rm = TRUE) else NA_real_,
      mean_pev = if ("PEV" %in% names(combined_uncertainty)) mean(.data$PEV, na.rm = TRUE) else NA_real_,
      mean_interval_width = if ("interval_width" %in% names(combined_uncertainty)) mean(.data$interval_width, na.rm = TRUE) else NA_real_,
      empirical_coverage = if ("covered" %in% names(combined_uncertainty)) mean(.data$covered, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))

  by_row_group_cols <- c("trait", "model", "row_id")
  if ("Env" %in% names(combined_uncertainty)) {
    by_row_group_cols <- c(by_row_group_cols, "Env")
  }
  aggregated_by_row <- combined_uncertainty |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by_row_group_cols))) |>
    dplyr::summarise(
      y = mean(.data$y, na.rm = TRUE),
      yhat = mean(.data$yhat, na.rm = TRUE),
      mean_absolute_error = mean(.data$absolute_error, na.rm = TRUE),
      mean_squared_error = mean(.data$squared_error, na.rm = TRUE),
      mean_standard_error = if ("Standard_error" %in% names(combined_uncertainty)) mean(.data$Standard_error, na.rm = TRUE) else NA_real_,
      mean_pev = if ("PEV" %in% names(combined_uncertainty)) mean(.data$PEV, na.rm = TRUE) else NA_real_,
      mean_interval_width = if ("interval_width" %in% names(combined_uncertainty)) mean(.data$interval_width, na.rm = TRUE) else NA_real_,
      empirical_coverage = if ("covered" %in% names(combined_uncertainty)) mean(.data$covered, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(by_row_group_cols)))

  list(
    combined_uncertainty = combined_uncertainty,
    aggregated_uncertainty = aggregated_uncertainty,
    aggregated_by_row = aggregated_by_row
  )
}

aggregate_cv_gaussian_risk_summaries <- function(cv_results_data) {
  risk_entries <- lapply(cv_results_data, function(x) {
    pred_df <- x$ypred_cv_Reps_all
    if (is.null(pred_df) || !is.data.frame(pred_df) || !all(c("yhat", "Prediction_risk_score") %in% names(pred_df))) {
      return(NULL)
    }

    row_id <- if ("row_id" %in% names(pred_df)) pred_df$row_id else seq_len(nrow(pred_df))
    cv_role <- if ("cv_role" %in% names(pred_df)) pred_df$cv_role else NA_character_
    env_col <- intersect(names(pred_df), c("Env", "env"))
    env_val <- if (length(env_col)) pred_df[[env_col[[1]]]] else NULL

    out <- data.frame(
      trait = x$trait,
      model = x$model,
      rep = x$rep,
      row_id = row_id,
      yhat = suppressWarnings(as.numeric(pred_df$yhat)),
      cv_role = cv_role,
      Prediction_confidence = if ("Prediction_confidence" %in% names(pred_df)) suppressWarnings(as.numeric(pred_df$Prediction_confidence)) else NA_real_,
      Prediction_unfamiliarity = if ("Prediction_unfamiliarity" %in% names(pred_df)) suppressWarnings(as.numeric(pred_df$Prediction_unfamiliarity)) else NA_real_,
      Prediction_unfamiliarity_remarks = if ("Prediction_unfamiliarity_remarks" %in% names(pred_df)) as.character(pred_df$Prediction_unfamiliarity_remarks) else NA_character_,
      Prediction_unfamiliarity_basis = if ("Prediction_unfamiliarity_basis" %in% names(pred_df)) as.character(pred_df$Prediction_unfamiliarity_basis) else NA_character_,
      Prediction_risk_score = suppressWarnings(as.numeric(pred_df$Prediction_risk_score)),
      Prediction_risk_remarks = if ("Prediction_risk_remarks" %in% names(pred_df)) as.character(pred_df$Prediction_risk_remarks) else NA_character_,
      Prediction_risk_basis = if ("Prediction_risk_basis" %in% names(pred_df)) as.character(pred_df$Prediction_risk_basis) else NA_character_,
      Rank_instability_risk = if ("Rank_instability_risk" %in% names(pred_df)) suppressWarnings(as.numeric(pred_df$Rank_instability_risk)) else suppressWarnings(as.numeric(pred_df$Prediction_risk_score)),
      Rank_instability_risk_remarks = if ("Rank_instability_risk_remarks" %in% names(pred_df)) as.character(pred_df$Rank_instability_risk_remarks) else if ("Prediction_risk_remarks" %in% names(pred_df)) as.character(pred_df$Prediction_risk_remarks) else NA_character_,
      Rank_instability_risk_basis = if ("Rank_instability_risk_basis" %in% names(pred_df)) as.character(pred_df$Rank_instability_risk_basis) else if ("Prediction_risk_basis" %in% names(pred_df)) as.character(pred_df$Prediction_risk_basis) else NA_character_,
      stringsAsFactors = FALSE
    )
    if (!is.null(env_val)) {
      out$Env <- env_val
    }
    out
  })

  risk_entries <- Filter(Negate(is.null), risk_entries)
  if (!length(risk_entries)) {
    return(NULL)
  }

  combined_risk <- do.call(rbind, risk_entries)
  if ("cv_role" %in% names(combined_risk) && any(!is.na(combined_risk$cv_role))) {
    combined_risk <- combined_risk[combined_risk$cv_role %in% c("test", NA), , drop = FALSE]
  }
  if (!nrow(combined_risk)) {
    return(NULL)
  }

  group_cols <- c("trait", "model")
  if ("Env" %in% names(combined_risk)) {
    group_cols <- c(group_cols, "Env")
  }

  aggregated_risk <- combined_risk |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_prediction_confidence = if ("Prediction_confidence" %in% names(combined_risk)) mean(.data$Prediction_confidence, na.rm = TRUE) else NA_real_,
      mean_prediction_unfamiliarity = if ("Prediction_unfamiliarity" %in% names(combined_risk)) mean(.data$Prediction_unfamiliarity, na.rm = TRUE) else NA_real_,
      mean_prediction_risk_score = mean(.data$Prediction_risk_score, na.rm = TRUE),
      mean_rank_instability_risk = mean(.data$Rank_instability_risk, na.rm = TRUE),
      low_risk_percentage = 100 * mean(.data$Prediction_risk_remarks == "Low Risk", na.rm = TRUE),
      moderate_risk_percentage = 100 * mean(.data$Prediction_risk_remarks == "Moderate Risk", na.rm = TRUE),
      high_risk_percentage = 100 * mean(.data$Prediction_risk_remarks == "High Risk", na.rm = TRUE),
      low_unfamiliarity_percentage = 100 * mean(.data$Prediction_unfamiliarity_remarks == "Familiar", na.rm = TRUE),
      moderate_unfamiliarity_percentage = 100 * mean(.data$Prediction_unfamiliarity_remarks == "Moderately Unfamiliar", na.rm = TRUE),
      high_unfamiliarity_percentage = 100 * mean(.data$Prediction_unfamiliarity_remarks == "Highly Unfamiliar", na.rm = TRUE),
      prediction_unfamiliarity_basis = paste(unique(stats::na.omit(.data$Prediction_unfamiliarity_basis)), collapse = "; "),
      prediction_risk_basis = paste(unique(stats::na.omit(.data$Prediction_risk_basis)), collapse = "; "),
      rank_instability_risk_basis = paste(unique(stats::na.omit(.data$Rank_instability_risk_basis)), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))

  aggregated_by_row <- combined_risk |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "row_id")))) |>
    dplyr::summarise(
      yhat = mean(.data$yhat, na.rm = TRUE),
      mean_prediction_confidence = if ("Prediction_confidence" %in% names(combined_risk)) mean(.data$Prediction_confidence, na.rm = TRUE) else NA_real_,
      mean_prediction_unfamiliarity = if ("Prediction_unfamiliarity" %in% names(combined_risk)) mean(.data$Prediction_unfamiliarity, na.rm = TRUE) else NA_real_,
      mean_prediction_risk_score = mean(.data$Prediction_risk_score, na.rm = TRUE),
      mean_rank_instability_risk = mean(.data$Rank_instability_risk, na.rm = TRUE),
      modal_prediction_risk_remark = names(sort(table(.data$Prediction_risk_remarks), decreasing = TRUE))[1],
      modal_prediction_unfamiliarity_remark = names(sort(table(.data$Prediction_unfamiliarity_remarks), decreasing = TRUE))[1],
      modal_rank_instability_risk_remark = names(sort(table(.data$Rank_instability_risk_remarks), decreasing = TRUE))[1],
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(c(group_cols, "row_id"))))

  list(
    combined_risk = combined_risk,
    aggregated_risk = aggregated_risk,
    aggregated_by_row = aggregated_by_row
  )
}

build_cv_gaussian_risk_plots <- function(gaussian_risk_summaries) {
  if (is.null(gaussian_risk_summaries) ||
      !is.list(gaussian_risk_summaries) ||
      is.null(gaussian_risk_summaries$aggregated_risk) ||
      !is.data.frame(gaussian_risk_summaries$aggregated_risk) ||
      !nrow(gaussian_risk_summaries$aggregated_risk)) {
    return(NULL)
  }

  agg <- gaussian_risk_summaries$aggregated_risk
  facet_formula <- if ("Env" %in% names(agg)) {
    stats::as.formula("~ trait + Env")
  } else {
    stats::as.formula("~ trait")
  }

  unfamiliarity_bar <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = .data$model, y = .data$mean_prediction_unfamiliarity, fill = .data$model)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(facet_formula, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Gaussian Mean Prediction Unfamiliarity by Model",
      subtitle = "Lower values indicate predictions that stay closer to the training feature distribution",
      x = "Model",
      y = "Mean Prediction Unfamiliarity"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  risk_bar <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = .data$model, y = .data$mean_rank_instability_risk, fill = .data$model)
  ) +
    ggplot2::geom_col() +
      ggplot2::facet_wrap(facet_formula, scales = "free_x") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
      title = "Gaussian Mean Rank Instability Risk by Model",
      subtitle = "Lower values indicate lower calibrated ranking risk on held-out predictions",
      x = "Model",
      y = "Mean Rank Instability Risk"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  id_cols <- intersect(c("trait", "model", "Env"), names(agg))
  risk_pct <- rbind(
    cbind(
      agg[, id_cols, drop = FALSE],
      risk_band = "Low Risk",
      percentage = agg$low_risk_percentage,
      stringsAsFactors = FALSE
    ),
    cbind(
      agg[, id_cols, drop = FALSE],
      risk_band = "Moderate Risk",
      percentage = agg$moderate_risk_percentage,
      stringsAsFactors = FALSE
    ),
    cbind(
      agg[, id_cols, drop = FALSE],
      risk_band = "High Risk",
      percentage = agg$high_risk_percentage,
      stringsAsFactors = FALSE
    )
  )
  risk_pct$risk_band <- factor(
    risk_pct$risk_band,
    levels = c("Low Risk", "Moderate Risk", "High Risk")
  )

  risk_stack <- ggplot2::ggplot(
    risk_pct,
    ggplot2::aes(x = .data$model, y = .data$percentage, fill = .data$risk_band)
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::facet_wrap(facet_formula, scales = "free_x") +
      ggplot2::scale_fill_manual(values = c("Low Risk" = "green", "Moderate Risk" = "lightblue", "High Risk" = "red")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
      title = "Gaussian Risk Remark Distribution by Model",
      subtitle = "Stacked Low/Moderate/High calibrated ranking-risk remarks",
      x = "Model",
      y = "Percentage",
      fill = "Risk Band"
      ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(
    mean_prediction_unfamiliarity = unfamiliarity_bar,
    mean_rank_instability_risk = risk_bar,
    risk_distribution = risk_stack
  )
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
                                              metric_for_ranking = "auto",
                                              ranking_tie_breakers = NULL,
                                              positive_class = NULL,
                                              render_interactive = TRUE) {


  aggregated_data_list <- aggregate_metrics(cv_results_data)
  classification_metrics <- aggregate_cv_classification_metrics(
    cv_results_data, aggregated_data_list, eval_metrics
  )
  classification_reliability <- aggregate_cv_classification_reliability(cv_results_data)
  classification_confusion_matrices <-
    aggregate_cv_classification_confusion_matrices(classification_reliability)
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
  gaussian_uncertainty_summaries <- aggregate_cv_gaussian_uncertainty_summaries(cv_results_data)
  gaussian_risk_summaries <- aggregate_cv_gaussian_risk_summaries(cv_results_data)
  plotly_available <- isTRUE(render_interactive) && requireNamespace("plotly", quietly = TRUE)

  plot_reps_list <- lapply(eval_metrics, function(metric) {
    metric_sym <- rlang::sym(metric)

    ggplot_boxplot_reps <- ggplot2::ggplot(
      aggregated_data_list[["combined_df"]],
      ggplot2::aes(x = .data$model, y = !!metric_sym, fill = .data$model)
    ) +
      #ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ trait, scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Boxplot of Model", metric, "by Trait"),
                    x = "Model",
                    y = metric)

    plotly_boxplot_reps <- if (isTRUE(plotly_available)) {
      tryCatch(
        plotly::ggplotly(ggplot_boxplot_reps),
        error = function(e) NULL
      )
    } else {
      NULL
    }

    list(metric = metric,
         ggplot_boxplot_reps = ggplot_boxplot_reps,
         plotly_boxplot_reps = plotly_boxplot_reps)
  })
  names(plot_reps_list) <-  eval_metrics
  ###################
  plot_mean_list <- lapply(eval_metrics, function(metric) {
    metric_sym <- rlang::sym(metric)

    # ggplot_barplot_mean <-  ggplot2::ggplot(aggregated_data_list[["aggregated_across_reps"]], ggplot2::aes_string(x = "model", y = metric, fill = "model")) +
    #   ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    #   ggplot2::facet_wrap(~trait, scales = "free_x") +
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    #   ggplot2::labs(title = "Mean Model Performance by Trait",
    #                 x = "Model",
    #                 y = "Mean Accuracy")
    ######
    ggplot_lineplot_mean <- ggplot2::ggplot(
      aggregated_data_list[["aggregated_across_reps"]],
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


    #plotly_lineplot_mean <- plotly::ggplotly(ggplot_lineplot_mean)

    list(metric = metric,
         ggplot_lineplot_mean = ggplot_lineplot_mean
         #plotly_lineplot_mean = plotly_lineplot_mean
         )
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
  best_models_list <- lapply(eval_metrics, function(metric) {
    rank_models(
      data = aggregated_data_list[["aggregated_across_reps"]],
      metric = metric,
      tie_breakers = if (identical(metric, metric_for_ranking)) ranking_tie_breakers else NULL
    )
  })

  names(best_models_list) <- eval_metrics


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
                is_multi = FALSE
              ),
              classification_probability_summaries = aggregate_cv_probability_summaries(cv_results_data),
              classification_metrics = classification_metrics,
              classification_reliability = classification_reliability,
              gaussian_uncertainty_summaries = gaussian_uncertainty_summaries,
              gaussian_risk_summaries = gaussian_risk_summaries,
              gaussian_risk_plots = build_cv_gaussian_risk_plots(gaussian_risk_summaries)))
}


 #MM = cv_single_loc_result_plot_process(TT, eval_metrics)

