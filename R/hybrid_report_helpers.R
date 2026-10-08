#' Null-coalesce helper for the hybrid-report code
#'
#' Returns `x` if it is non-NULL and not a single `NA`; otherwise returns the
#' fallback `y`. A small convenience used across the hybrid-report helpers.
#'
#' @param x Primary value (returned when non-NULL and not a single `NA`).
#' @param y Fallback value returned when `x` is `NULL` or a single `NA`.
#'
#' @return Either `x` or `y`.
#' @export
gp_hybrid_null_coalesce <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

#' Build Hybrid Benchmark Ranking Reports
#'
#' Create ranked hybrid benchmark summaries and comparison plots from a
#' standardized benchmark table with one row per model/split combination.
#'
#' @param summary_table Data frame with at least \code{family}, \code{model},
#'   \code{split}, \code{rmse}, and \code{mae} columns.
#' @param masked_split_name Split label used for direct masked true prediction.
#' @param cv_split_name Split label used for the hybrid CV ranking.
#' @param output_dir Optional directory for CSV/PDF export.
#' @param export Logical; when \code{TRUE}, write CSV and PDF outputs to
#'   \code{output_dir}.
#'
#' @return A list containing ranked tables, recommendations, and plots.
#' @export
hybrid_benchmark_report <- function(summary_table,
                                    masked_split_name = "masked_true_prediction",
                                    cv_split_name = "Hybrid_One_New_Parent",
                                    output_dir = NULL,
                                    export = FALSE) {
  tab <- as.data.frame(summary_table, stringsAsFactors = FALSE)
  required_cols <- c("family", "model", "split", "rmse", "mae")
  missing_cols <- setdiff(required_cols, names(tab))
  if (length(missing_cols)) {
    stop(
      paste(
        "summary_table is missing required columns:",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ranked_masked <- data.frame()
  ranked_cv <- data.frame()
  recommendations <- data.frame()

  masked_rows <- tab[tab$split == masked_split_name, , drop = FALSE]
  cv_rows <- tab[tab$split == cv_split_name, , drop = FALSE]

  if (nrow(masked_rows)) {
    ranked_masked <- masked_rows[order(masked_rows$rmse, masked_rows$mae), , drop = FALSE]
    ranked_masked$rank <- seq_len(nrow(ranked_masked))
    ranked_masked <- ranked_masked[, c("rank", "family", "model", "split", "rmse", "mae"), drop = FALSE]
  }
  if (nrow(cv_rows)) {
    ranked_cv <- cv_rows[order(cv_rows$rmse, cv_rows$mae), , drop = FALSE]
    ranked_cv$rank <- seq_len(nrow(ranked_cv))
    ranked_cv <- ranked_cv[, c("rank", "family", "model", "split", "rmse", "mae"), drop = FALSE]
  }

  rec_rows <- list()
  ridx <- 1L
  if (nrow(ranked_masked)) {
    rec_rows[[ridx]] <- data.frame(
      recommendation = "Best masked true-prediction model",
      family = ranked_masked$family[1],
      model = ranked_masked$model[1],
      basis = paste("lowest RMSE then MAE on", masked_split_name),
      stringsAsFactors = FALSE
    )
    ridx <- ridx + 1L
  }
  if (nrow(ranked_cv)) {
    rec_rows[[ridx]] <- data.frame(
      recommendation = paste("Best", cv_split_name, "model"),
      family = ranked_cv$family[1],
      model = ranked_cv$model[1],
      basis = paste("lowest RMSE then MAE on", cv_split_name),
      stringsAsFactors = FALSE
    )
  }
  if (length(rec_rows)) {
    recommendations <- do.call(rbind, rec_rows)
  }

  plot_tab <- tab
  plot_tab$model <- factor(plot_tab$model, levels = unique(plot_tab$model))
  rmse_plot <- ggplot2::ggplot(
    plot_tab,
    ggplot2::aes(x = model, y = rmse, fill = split)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Hybrid benchmark RMSE by model",
      x = "Model",
      y = "RMSE",
      fill = "Split"
    )

  mae_plot <- ggplot2::ggplot(
    plot_tab,
    ggplot2::aes(x = model, y = mae, fill = split)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Hybrid benchmark MAE by model",
      x = "Model",
      y = "MAE",
      fill = "Split"
    )

  out <- list(
    summary = tab,
    ranked_masked_true_prediction = ranked_masked,
    ranked_cv = ranked_cv,
    recommendations = recommendations,
    plots = list(
      rmse_by_model = rmse_plot,
      mae_by_model = mae_plot
    )
  )

  if (isTRUE(export)) {
    if (is.null(output_dir) || !nzchar(output_dir)) {
      stop("Provide output_dir when export = TRUE.", call. = FALSE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(tab, file.path(output_dir, "hybrid_shortlist_summary.csv"), row.names = FALSE)
    utils::write.csv(ranked_masked, file.path(output_dir, "hybrid_shortlist_ranked_masked_true_prediction.csv"), row.names = FALSE)
    utils::write.csv(ranked_cv, file.path(output_dir, "hybrid_shortlist_ranked_cv.csv"), row.names = FALSE)
    utils::write.csv(recommendations, file.path(output_dir, "hybrid_shortlist_recommendations.csv"), row.names = FALSE)
    ggplot2::ggsave(
      filename = file.path(output_dir, "hybrid_shortlist_rmse_by_model.pdf"),
      plot = rmse_plot,
      width = 10,
      height = 6,
      units = "in"
    )
    ggplot2::ggsave(
      filename = file.path(output_dir, "hybrid_shortlist_mae_by_model.pdf"),
      plot = mae_plot,
      width = 10,
      height = 6,
      units = "in"
    )
  }

  out
}

gp_hybrid_extract_mode_model <- function(run_object,
                                         fallback_family = NULL,
                                         fallback_model = NULL) {
  family <- fallback_family
  model <- fallback_model

  summ <- gp_hybrid_null_coalesce(run_object$summary_statistic, run_object$res_summary_stat)
  if (is.list(summ) && is.data.frame(summ$summary_statistics)) {
    stats_tab <- summ$summary_statistics
    stat_vals <- setNames(as.character(stats_tab$summary), as.character(stats_tab$stat))
    family <- gp_hybrid_null_coalesce(family, gp_hybrid_null_coalesce(stat_vals[["Mode"]], gp_hybrid_null_coalesce(stat_vals[["mode"]], NA_character_)))
    model <- gp_hybrid_null_coalesce(model, gp_hybrid_null_coalesce(stat_vals[["Model_Type"]], gp_hybrid_null_coalesce(stat_vals[["model_type"]], NA_character_)))
  }

  if (is.null(family) || !nzchar(family)) {
    ctx <- run_object$run_metadata
    if (is.data.frame(ctx) && all(c("key", "value") %in% names(ctx))) {
      ctx_map <- setNames(as.character(ctx$value), as.character(ctx$key))
      family <- gp_hybrid_null_coalesce(ctx_map[["context"]], family)
    }
  }

  list(
    family = as.character(gp_hybrid_null_coalesce(family, NA_character_)),
    model = as.character(gp_hybrid_null_coalesce(model, NA_character_))
  )
}

gp_hybrid_metrics_from_truth <- function(pred_df,
                                         truth_lookup,
                                         hybrid_id_col,
                                         truth_col,
                                         masked_split_name) {
  if (is.null(pred_df) || !nrow(pred_df) || is.null(truth_lookup)) {
    return(data.frame())
  }
  if (!all(c(hybrid_id_col, truth_col) %in% names(truth_lookup))) {
    stop(
      paste(
        "truth_lookup must contain columns:",
        paste(c(hybrid_id_col, truth_col), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!all(c(hybrid_id_col, "Predicted_value", "Train_Test_Label") %in% names(pred_df))) {
    return(data.frame())
  }

  test_df <- merge(
    pred_df[pred_df$Train_Test_Label == "Test", , drop = FALSE],
    truth_lookup[, c(hybrid_id_col, truth_col), drop = FALSE],
    by = hybrid_id_col,
    all.x = TRUE,
    sort = FALSE
  )
  test_df <- test_df[!is.na(test_df[[truth_col]]), , drop = FALSE]
  if (!nrow(test_df)) {
    return(data.frame())
  }

  data.frame(
    split = masked_split_name,
    rmse = sqrt(mean((test_df[[truth_col]] - test_df$Predicted_value)^2, na.rm = TRUE)),
    mae = mean(abs(test_df[[truth_col]] - test_df$Predicted_value), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_metrics_from_cv <- function(cv_object,
                                      cv_split_name = NULL) {
  proc <- gp_hybrid_null_coalesce(cv_object$cv_results_processed, NULL)
  if (is.list(proc) && is.data.frame(proc$hybrid_metric_summary) && nrow(proc$hybrid_metric_summary)) {
    tab <- proc$hybrid_metric_summary
    split_col <- if ("cv_scenario" %in% names(tab)) "cv_scenario" else NULL
    if (!is.null(cv_split_name) && !is.null(split_col)) {
      tab <- tab[tab[[split_col]] == cv_split_name, , drop = FALSE]
    }
    if (!nrow(tab)) {
      return(data.frame())
    }
    if (!all(c("root_mean_squared_error", "mean_absolute_error") %in% names(tab))) {
      return(data.frame())
    }
    split_vals <- if (!is.null(split_col)) as.character(tab[[split_col]]) else rep(gp_hybrid_null_coalesce(cv_split_name, "hybrid_cv"), nrow(tab))
    return(data.frame(
      split = split_vals,
      rmse = as.numeric(tab$root_mean_squared_error),
      mae = as.numeric(tab$mean_absolute_error),
      stringsAsFactors = FALSE
    ))
  }

  raw <- gp_hybrid_null_coalesce(cv_object$cv_results_raw, cv_object$cv_results)
  if (is.list(raw) && length(raw)) {
    rows <- list()
    idx <- 1L
    for (entry in raw) {
      met <- gp_hybrid_null_coalesce(entry$eval_metrics_reps, NULL)
      if (!is.data.frame(met) || !nrow(met)) {
        next
      }
      if (!all(c("root_mean_squared_error", "mean_absolute_error") %in% names(met))) {
        next
      }
      split_val <- gp_hybrid_null_coalesce(
        gp_hybrid_null_coalesce(entry$cv_info$scenario, entry$cv_scenario),
        gp_hybrid_null_coalesce(cv_split_name, "hybrid_cv")
      )
      rows[[idx]] <- data.frame(
        split = as.character(split_val),
        rmse = mean(met$root_mean_squared_error, na.rm = TRUE),
        mae = mean(met$mean_absolute_error, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    if (length(rows)) {
      out <- do.call(rbind, rows)
      if (!is.null(cv_split_name)) {
        out <- out[out$split == cv_split_name, , drop = FALSE]
      }
      return(out)
    }
  }

  data.frame()
}

#' Build a Hybrid Benchmark Table From Run Objects
#'
#' @param run_entries A list of hybrid run entries. Each entry can be a single
#'   `model_execute()` result object, or a list containing `true_prediction`,
#'   `cv`, `family`, `model`, and optional `cv_split_name`.
#' @param truth_lookup Optional lookup table used to score masked true
#'   prediction test rows.
#' @param hybrid_id_col Hybrid ID column shared by `truth_lookup` and hybrid
#'   prediction outputs.
#' @param truth_col Observed truth column in `truth_lookup`.
#' @param masked_split_name Label used for masked true-prediction rows.
#'
#' @return A standardized benchmark table with `family`, `model`, `split`,
#'   `rmse`, and `mae` columns.
#' @export
hybrid_benchmark_table_from_runs <- function(run_entries,
                                             truth_lookup = NULL,
                                             hybrid_id_col = "HybridID",
                                             truth_col = "Observed_truth",
                                             masked_split_name = "masked_true_prediction") {
  if (!is.list(run_entries) || !length(run_entries)) {
    stop("run_entries must be a non-empty list.", call. = FALSE)
  }

  rows <- list()
  idx <- 1L

  entry_names <- names(run_entries)
  if (is.null(entry_names)) {
    entry_names <- rep("", length(run_entries))
  }
  for (entry_name in entry_names) {
    entry <- run_entries[[if (nzchar(entry_name)) entry_name else idx]]
    if (is.null(entry)) {
      idx <- idx + 1L
      next
    }

    is_wrapped <- is.list(entry) && (("true_prediction" %in% names(entry)) || ("cv" %in% names(entry)))
    true_obj <- if (is_wrapped) gp_hybrid_null_coalesce(entry$true_prediction, NULL) else if (is.list(entry) && !is.null(entry$model_results)) entry else NULL
    cv_obj <- if (is_wrapped) gp_hybrid_null_coalesce(entry$cv, NULL) else if (is.list(entry) && (!is.null(entry$cv_results_processed) || !is.null(entry$cv_results_raw))) entry else NULL

    meta <- gp_hybrid_extract_mode_model(
      run_object = gp_hybrid_null_coalesce(true_obj, gp_hybrid_null_coalesce(cv_obj, list())),
      fallback_family = if (is_wrapped) gp_hybrid_null_coalesce(entry$family, NULL) else NULL,
      fallback_model = if (is_wrapped) gp_hybrid_null_coalesce(entry$model, NULL) else NULL
    )
    family <- meta$family
    model <- meta$model
    cv_split_name <- if (is_wrapped) gp_hybrid_null_coalesce(entry$cv_split_name, NULL) else NULL

    if (!is.null(true_obj)) {
      pred_df <- gp_hybrid_null_coalesce(true_obj$model_results$predicted_values, true_obj$predicted_values)
      true_metrics <- gp_hybrid_metrics_from_truth(
        pred_df = pred_df,
        truth_lookup = truth_lookup,
        hybrid_id_col = hybrid_id_col,
        truth_col = truth_col,
        masked_split_name = masked_split_name
      )
      if (nrow(true_metrics)) {
        rows[[length(rows) + 1L]] <- data.frame(
          family = family,
          model = model,
          split = true_metrics$split,
          rmse = true_metrics$rmse,
          mae = true_metrics$mae,
          stringsAsFactors = FALSE
        )
      }
    }

    if (!is.null(cv_obj)) {
      cv_metrics <- gp_hybrid_metrics_from_cv(
        cv_object = cv_obj,
        cv_split_name = cv_split_name
      )
      if (nrow(cv_metrics)) {
        rows[[length(rows) + 1L]] <- data.frame(
          family = rep(family, nrow(cv_metrics)),
          model = rep(model, nrow(cv_metrics)),
          split = cv_metrics$split,
          rmse = cv_metrics$rmse,
          mae = cv_metrics$mae,
          stringsAsFactors = FALSE
        )
      }
    }

    idx <- idx + 1L
  }

  if (!length(rows)) {
    return(data.frame(
      family = character(),
      model = character(),
      split = character(),
      rmse = numeric(),
      mae = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}

#' Build a Processed Hybrid Benchmark Bundle
#'
#' @param run_entries A list of hybrid run objects or wrapped hybrid entries.
#' @param truth_lookup Optional lookup table used to score masked true
#'   prediction rows.
#' @param hybrid_id_col Hybrid ID column shared by `truth_lookup` and hybrid
#'   prediction outputs.
#' @param truth_col Observed truth column in `truth_lookup`.
#' @param masked_split_name Split label used for direct masked true prediction.
#' @param cv_split_name Split label used for the primary hybrid CV ranking.
#' @param output_dir Optional directory for CSV/PDF export.
#' @param export Logical; when `TRUE`, export processed CSV/PDF artifacts.
#'
#' @return A processed hybrid benchmark bundle with the benchmark table,
#'   ranked summaries, recommendations, and comparison plots.
#' @export
hybrid_benchmark_process <- function(run_entries,
                                     truth_lookup = NULL,
                                     hybrid_id_col = "HybridID",
                                     truth_col = "Observed_truth",
                                     masked_split_name = "masked_true_prediction",
                                     cv_split_name = "Hybrid_One_New_Parent",
                                     output_dir = NULL,
                                     export = FALSE) {
  benchmark_table <- hybrid_benchmark_table_from_runs(
    run_entries = run_entries,
    truth_lookup = truth_lookup,
    hybrid_id_col = hybrid_id_col,
    truth_col = truth_col,
    masked_split_name = masked_split_name
  )

  report <- hybrid_benchmark_report(
    summary_table = benchmark_table,
    masked_split_name = masked_split_name,
    cv_split_name = cv_split_name,
    output_dir = output_dir,
    export = export
  )

  list(
    hybrid_benchmark_table = benchmark_table,
    hybrid_benchmark_summary = report$summary,
    hybrid_benchmark_ranked_masked_true_prediction = report$ranked_masked_true_prediction,
    hybrid_benchmark_ranked_cv = report$ranked_cv,
    hybrid_benchmark_recommendations = report$recommendations,
    hybrid_benchmark_plots = report$plots
  )
}
