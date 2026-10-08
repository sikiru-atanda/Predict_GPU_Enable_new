#' Title
#'
#' @param res_plot
#' @param plot_filename
#' @param plot_extension
#' @param plot_width
#' @param plot_height
#' @param plot_units
#' @param plot_dpi
#'
#' @return
#' @export
#'
#' @examples
save_ggplot <- function(res_plot,
                        plot_filename = "trait",
                        plot_extension = "pdf",
                        plot_width = 17,
                        plot_height = 12,
                        plot_units = "in",
                        plot_dpi = 300) {


  # save plot plot with user defined name and parameters
  if(is.null(plot_filename)) plot_filename <- "trait"
  plot_name <- paste(plot_filename,plot_extension, sep = ".")

       ggplot2::ggsave(filename = plot_name,
                       plot = res_plot,
                       width = plot_width,
                       height = plot_height,
                       units = plot_units,
                       dpi = plot_dpi)


}


#' Title
#'
#' @param pathout
#' @param geno_omic_files
#'
#' @return
#' @export
#'
#' @examples
zipMMatrixModelReady <- function(pathout,
                                 geno_omic_files) {
  is_geno <- grepl("Geno", geno_omic_files, ignore.case = TRUE)
  is_omic <- grepl("Omic", geno_omic_files, ignore.case = TRUE)

  zip_name <- if (any(is_geno) && any(is_omic)) {
    "Geno_Omics_Clean_data.zip"
  } else if (any(is_geno)) {
    "Geno_Clean_data.zip"
  } else if (any(is_omic)) {
    "Omic_Clean_data.zip"
  } else {
    "X_variables_Clean_data.zip"
  }

  zip(file.path(pathout, zip_name), files = geno_omic_files, flags = "-q")
}


# Keep exported run directories discoverable without exposing platform-invalid
# path characters. Public model spelling (including hyphens and underscores)
# is otherwise preserved so the directory identifies the model the user chose.
gp_output_path_component <- function(value, fallback = "output") {
  value <- as.character(value %||% character())
  value <- value[!is.na(value) & nzchar(trimws(value))]
  value <- if (length(value)) trimws(value[[1L]]) else fallback
  value <- gsub("[[:cntrl:]]", "_", value)
  value <- gsub("[<>:\"/\\\\|?*]", "_", value)
  value <- sub("[. ]+$", "", value)
  if (!nzchar(value)) fallback else value
}

gp_output_public_model_labels <- function(GS_model = NULL) {
  models <- as.character(GS_model %||% character())
  models <- models[!is.na(models) & nzchar(trimws(models))]
  if (!length(models)) {
    return(character())
  }
  labels <- if (exists("gp_public_model_label", mode = "function")) {
    tryCatch(gp_public_model_label(models), error = function(e) models)
  } else {
    models
  }
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(trimws(labels))]
  unique(vapply(
    labels,
    gp_output_path_component,
    character(1L),
    fallback = "model"
  ))
}

gp_output_directory_stem <- function(output_file_name = "output",
                                     GS_model = NULL) {
  stem <- gp_output_path_component(output_file_name, fallback = "output")
  labels <- gp_output_public_model_labels(GS_model)
  if (!length(labels)) {
    return(stem)
  }

  # Some legacy callers already put the exact public model label in
  # plot_filename. Do not repeat those labels; append only missing models.
  present <- vapply(
    labels,
    function(label) grepl(label, stem, fixed = TRUE),
    logical(1L)
  )
  missing <- labels[!present]
  if (!length(missing)) {
    return(stem)
  }
  marker <- if (length(labels) == 1L) "model-" else "models-"
  paste0(stem, "__", marker, paste(missing, collapse = "+"))
}

format_model_table_export_name <- function(prefix) {
  # Lower-case first: the character class would otherwise treat capitals as
  # separators ("Predicted_value" -> "redicted_value").
  key <- gsub("[^a-z0-9]+", "_", tolower(as.character(prefix %||% "")))
  key <- gsub("^_+|_+$", "", key)
  if (key %in% c("predicted_value", "predicted_values", "prediction_value", "prediction_values")) {
    return("Predicted_Value")
  }
  # One public file name whatever the engine's list key ("Variance_components"
  # or "variance_components"); file names are case-sensitive on Linux/macOS.
  if (identical(key, "variance_components")) {
    return("variance_components")
  }
  as.character(prefix)
}

is_predicted_value_export_name <- function(prefix) {
  identical(format_model_table_export_name(prefix), "Predicted_Value")
}

standard_predicted_value_export_table <- function(x) {
  df <- tryCatch(
    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    return(x)
  }
  # Hybrid prediction tables have a distinct public identifier contract.  They
  # must be recognized before the generic genomic formatter, which otherwise
  # replaces Hybrid_ID/Female/Male with a synthetic GID column.
  if (exists("gp_is_hybrid_prediction_table", mode = "function") &&
      isTRUE(gp_is_hybrid_prediction_table(df)) &&
      exists("gp_format_hybrid_prediction_table", mode = "function")) {
    return(gp_format_hybrid_prediction_table(df))
  }
  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(df)) &&
      exists("gp_format_classification_prediction_table", mode = "function")) {
    return(gp_format_classification_prediction_table(df))
  }
  if (exists("gp_standardize_prediction_table", mode = "function")) {
    df <- tryCatch(
      gp_standardize_prediction_table(df),
      error = function(e) df
    )
  }
  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(df)) &&
      exists("gp_format_classification_prediction_table", mode = "function")) {
    return(gp_format_classification_prediction_table(df))
  }
  if (exists("gp_format_gaussian_prediction_table", mode = "function")) {
    return(gp_format_gaussian_prediction_table(df))
  }
  common_cols <- c(
    "GID", "Predicted_value", "Train_Test_Label", "Observed_value",
    "Standard_error", "PEV", "lower_bound", "upper_bound",
    "Uncertainty", "Uncertainty_remarks", "Reliability", "Reliability_remarks"
  )
  keep <- intersect(common_cols, names(df))
  if (!length(keep)) df else df[, keep, drop = FALSE]
}

write_predicted_value_export_tables <- function(x, pathout) {
  common <- standard_predicted_value_export_table(x)
  write.csv(
    common,
    file.path(pathout, "Predicted_Value.csv"),
    row.names = FALSE
  )
  invisible(NULL)
}

prepare_table_for_export <- function(x, prefix) {
  export_name <- format_model_table_export_name(prefix)
  if (identical(export_name, "Predicted_Value") && is.data.frame(x)) {
    return(standard_predicted_value_export_table(x))
  }
  x
}

write_nested_csv_tables <- function(x, pathout, prefix) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  if (inherits(x, c("ggplot", "gtable", "plotly", "htmlwidget"))) {
    return(invisible(NULL))
  }
  export_prefix <- format_model_table_export_name(prefix)
  if (is.data.frame(x)) {
    if (is_predicted_value_export_name(prefix)) {
      write_predicted_value_export_tables(x, pathout)
      return(invisible(NULL))
    }
    x <- prepare_table_for_export(x, export_prefix)
    write.csv(
      x,
      file.path(pathout, paste0(export_prefix, ".csv")),
      row.names = FALSE
    )
    return(invisible(NULL))
  }
  if (is.matrix(x)) {
    write.csv(
      as.data.frame(x),
      file.path(pathout, paste0(export_prefix, ".csv")),
      row.names = TRUE
    )
    return(invisible(NULL))
  }
  if (is.atomic(x) && length(x) > 0L) {
    write.csv(
      data.frame(value = as.vector(x), stringsAsFactors = FALSE),
      file.path(pathout, paste0(export_prefix, ".csv")),
      row.names = FALSE
    )
    return(invisible(NULL))
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      nms <- paste0("item", seq_along(x))
    }
    for (i in seq_along(x)) {
      write_nested_csv_tables(
        x[[i]],
        pathout = pathout,
        prefix = paste(prefix, nms[[i]], sep = "_")
      )
    }
  }
  invisible(NULL)
}

format_processed_plot_export_name <- function(prefix) {
  plot_name_map <- c(
    "cv_results_processed_gaussian_risk_plots_mean_prediction_unfamiliarity" = "CV_gaussian_mean_prediction_unfamiliarity",
    "cv_results_processed_gaussian_risk_plots_mean_rank_instability_risk" = "CV_gaussian_mean_rank_instability_risk",
    "cv_results_processed_gaussian_risk_plots_risk_distribution" = "CV_gaussian_risk_distribution",
    "cv_results_processed_multitrait_cv_plots" = "CV_multitrait_observed_vs_predicted",
    "cv_results_processed_multitrait_rank_risk_plot" = "CV_multitrait_mean_rank_instability_risk",
    "cv_results_processed_hybrid_cv_plot" = "CV_hybrid_observed_vs_predicted",
    "cv_results_processed_hybrid_benchmark_plots_rmse_by_model" = "CV_hybrid_benchmark_rmse_by_model",
    "cv_results_processed_hybrid_benchmark_plots_mae_by_model" = "CV_hybrid_benchmark_mae_by_model"
  )
  if (prefix %in% names(plot_name_map)) {
    return(unname(plot_name_map[[prefix]]))
  }
  prefix
}


write_nested_plot_tables <- function(x,
                                      pathout,
                                      prefix,
                                      plot_extension = "pdf",
                                      plot_width = 17,
                                     plot_height = 12,
                                     plot_units = "in",
                                     plot_dpi = 300) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
    if (inherits(x, "ggplot") || inherits(x, "gtable")) {
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(pathout)
      export_name <- format_processed_plot_export_name(prefix)
      save_ggplot(
        res_plot = x,
        plot_filename = export_name,
        plot_extension = plot_extension,
        plot_width = plot_width,
        plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi
    )
    return(invisible(NULL))
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      nms <- paste0("item", seq_along(x))
    }
    for (i in seq_along(x)) {
      write_nested_plot_tables(
        x[[i]],
        pathout = pathout,
        prefix = paste(prefix, nms[[i]], sep = "_"),
        plot_extension = plot_extension,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_units = plot_units,
        plot_dpi = plot_dpi
      )
    }
  }
  invisible(NULL)
}


#' Title
#'
#' @param GS_model
#' @param res_model_output
#' @param res_summary_stat
#' @param res_plot
#' @param system_database
#'
#' @return
#' @export
#'
#' @examples



#' Title
#'
#' @param GS_model
#' @param res_model_output
#' @param res_summary_stat
#' @param res_plot
#' @param res_plot_mean
#' @param res_plot_result_diagnostic
#' @param test_diagonistic_plots
#' @param res_mod_results_cv_per_trait_model
#' @param res_plot_result_diagnostic_cv_only
#' @param cv_results_processed
#' @param cv_results_raw
#' @param geno_qc_stat
#' @param system_database
#' @param plot_filename
#' @param plot_extension
#' @param plot_width
#' @param plot_height
#' @param plot_units
#' @param plot_dpi
#' @param feature_selected
#' @param feature_score_metadata
#' @param run_metadata
#'
#' @return
#' @export
#'
#' @examples
results_handling <-  function(GS_model = NULL,
                              res_model_output = NULL,
                              res_summary_stat = NULL,
                              res_plot = NULL,
                              res_plot_mean = NULL,
                              res_plot_result_diagnostic = NULL,
                              test_diagonistic_plots = NULL,
                              res_mod_results_cv_per_trait_model = NULL,
                              res_plot_result_diagnostic_cv_only = NULL,
                              cv_results_processed = NULL,
                              cv_results_raw = NULL,
                              geno_qc_stat = NULL,
                              system_database = TRUE,
                              plot_filename = "trait",
                              #Plot_name_result_diagnostic = "result_diagnostic",
                              plot_extension = "pdf",
                              plot_width = 17,
                              plot_height = 12,
                              plot_units = "in",
                              plot_dpi = 300,
                              feature_selected = NULL,
                              feature_score_metadata = NULL,
                              run_metadata = NULL
                              ){

  if (exists("gp_standardize_hybrid_model_result", mode = "function")) {
    res_model_output <- tryCatch(
      gp_standardize_hybrid_model_result(
        res_model_output,
        model = GS_model
      ),
      error = function(e) res_model_output
    )
  }

  # Backend bridge payloads are implementation details, not public model
  # outputs. Keeping them here makes the generic recursive writer emit dozens
  # of misleading CSV files (and can expose large optimizer arrays). Public
  # predictions, variance components, parameters, correlations and diagnostics
  # are retained by the standard result contract before this point.
  internal_result_fields <- c(
    "raw_python_result", "raw_gp_result", "raw_gp_info",
    "gp_result", "gp_info", "gp_fit"
  )
  if (is.list(res_model_output) && length(names(res_model_output))) {
    res_model_output <- res_model_output[
      !names(res_model_output) %in% internal_result_fields
    ]
  }

  saveOutput <- function(res_model_output,
                         res_summary_stat,
                         pathout,
                         GS_model,
                         res_plot, res_plot_mean,
                         res_plot_result_diagnostic,
                         test_diagonistic_plots,
                         res_plot_result_diagnostic_cv_only,
                         feature_selected,
                         feature_score_metadata,
                         run_metadata,
                         cv_results_processed = NULL,
                         cv_results_raw = NULL) {
    #browser()
    if(!is.null(res_model_output)){
      # if ("bayes_model" %in% names(res_model_output)) {
      #   res_model_output <- res_model_output[names(res_model_output) != "bayes_model"]
      # }

      if ("bayes_result" %in% names(res_model_output)) {
        bayes_extra_output <- res_model_output[setdiff(names(res_model_output), c("bayes_result", "bayes_model"))]
        res_model_output <- utils::modifyList(res_model_output[["bayes_result"]], bayes_extra_output)
      }

      plot_slots <- c("diagnostic_plots", "diagnostic_tst_plot")
      if (any(plot_slots %in% names(res_model_output))) {
        res_model_output <- res_model_output[!names(res_model_output) %in% plot_slots]
      }
    #for (i in 1:length(res_model_output)) {

      for(i in seq_along(res_model_output)){
        output_value <- res_model_output[[i]]
        # NULL placeholders and empty public tables are not user results.  In
        # particular, exporting NULL trained-model placeholders produces tiny
        # unreadable-looking .RData files, while empty parameter tables produce
        # header-only CSVs.  Omit both directly at the package writer.
        if (is.null(output_value) ||
            (is.data.frame(output_value) && nrow(output_value) == 0L)) {
          next
        }
        ### check which output is a list
        if(length(class(output_value))==1 && !inherits(output_value, "list")){
          export_name <- format_model_table_export_name(names(res_model_output)[i])
          if (is_predicted_value_export_name(names(res_model_output)[i])) {
            write_predicted_value_export_tables(output_value, pathout)
            next
          }
          export_object <- prepare_table_for_export(output_value, export_name)
          if(tolower(names(res_model_output)[i])=="variance_components"){
            write.csv(export_object,
                      file.path(pathout,paste(export_name, "csv", sep = ".")),
                      row.names = FALSE)

          } else if(names(res_model_output)[i]=="bayes_model"){
            base::saveRDS(output_value,
                          paste(GS_model, "Bayes_model.RData", sep = "_"))

          } else if(names(res_model_output)[i]=="Asreml_model"){
            base::saveRDS(output_value,
                          "asreml_model.RData")

          } else if(names(res_model_output)[i]=="trained_model"){
            base::saveRDS(output_value,
                          paste(GS_model, "trained_model.RData", sep = "_"))

          }else{

            write.csv(export_object,
                      file.path(pathout,paste(export_name, "csv", sep = ".")),
                      row.names = FALSE)
          }
        } else if (length(class(output_value)) > 1 &&
                   "asreml.predict" %in% class(output_value) &&
                   "data.frame" %in% class(output_value)) {

          write.csv(output_value,
                    file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                    row.names = TRUE)


        } else if (length(class(output_value)) > 1 &&
                   "matrix" %in% class(output_value) &&
                   "array" %in% class(output_value)) {

          write.csv(output_value,
                    file.path(pathout,paste(names(res_model_output)[i], "csv", sep = ".")),
                    row.names = TRUE)

        } else {

          if(length(class(output_value))==1 && inherits(output_value,"list")){
            if(names(res_model_output)[i]=="M_matrix_model_ready"){
              for(j in seq_along(output_value)){
                write.csv(output_value[[j]],
                          file.path(pathout, paste(names(output_value[j]), "csv", sep = ".")),
                          row.names = TRUE)

              }
            } else if(names(res_model_output)[i]=="Covariance" | names(res_model_output)[i]=="Correlation"){
              for(j in seq_along(output_value)){
                write.csv(output_value[[j]],
                          file.path(pathout, paste(names(output_value[j]), "csv", sep = ".")),
                          row.names = TRUE)

              }

            } else if(names(res_model_output)[i]=="Asreml_model"){
              base::saveRDS(output_value,
                            "asreml_model.RData")

            } else if(names(res_model_output)[i]=="trained_model"){
              base::saveRDS(output_value,
                            "AI_trained_model.RData")

            } else {
              write_nested_csv_tables(
                output_value,
                pathout = pathout,
                prefix = names(res_model_output)[i]
              )

            }

          }

        }


      }

    #}

    }

    if(!is.null(res_summary_stat)){
      if (is.data.frame(res_summary_stat)) {
        if (nrow(res_summary_stat) > 0L) {
          write.csv(
            res_summary_stat,
            file.path(pathout, "summary_statistics.csv"),
            row.names = FALSE
          )
        }
      } else {
        for(s in seq_along(res_summary_stat)){
          if(!inherits(res_summary_stat[[s]], "list") &&
             !is.null(res_summary_stat[[s]]) &&
             !(is.data.frame(res_summary_stat[[s]]) && nrow(res_summary_stat[[s]]) == 0L)){
            write.csv(res_summary_stat[[s]],
                      file.path(pathout,paste(names(res_summary_stat)[s], "csv", sep = ".")),
                      row.names = FALSE)
          }
        }
      }
    }

    if (!is.null(run_metadata)) {
      write.csv(
        run_metadata,
        file.path(pathout, "Run_metadata.csv"),
        row.names = FALSE
      )
    }

    if (!is.null(cv_results_processed) && length(cv_results_processed) > 0) {
      for (nm in names(cv_results_processed)) {
        obj <- cv_results_processed[[nm]]
        if (inherits(obj, "ggplot") || inherits(obj, "gtable")) {
          write_nested_plot_tables(
            obj,
            pathout = pathout,
            prefix = paste0("cv_results_processed_", nm),
            plot_extension = plot_extension,
            plot_width = plot_width,
            plot_height = plot_height,
            plot_units = plot_units,
            plot_dpi = plot_dpi
          )
        } else if (is.data.frame(obj)) {
          write.csv(
            obj,
            file.path(pathout, paste0("cv_results_processed_", nm, ".csv")),
            row.names = FALSE
          )
        } else if (is.list(obj)) {
          write_nested_csv_tables(
            obj,
            pathout = pathout,
            prefix = paste0("cv_results_processed_", nm)
          )
          write_nested_plot_tables(
            obj,
            pathout = pathout,
            prefix = paste0("cv_results_processed_", nm),
            plot_extension = plot_extension,
            plot_width = plot_width,
            plot_height = plot_height,
            plot_units = plot_units,
            plot_dpi = plot_dpi
          )
        }
      }
    }

    if (!is.null(cv_results_raw) && length(cv_results_raw) > 0) {
      cv_dir <- file.path(pathout, "cv_results_raw")
      dir.create(cv_dir, showWarnings = FALSE)
      for (i in seq_along(cv_results_raw)) {
        entry <- cv_results_raw[[i]]
        if (is.null(entry) || !is.list(entry)) {
          next
        }
        trait <- entry[["trait"]] %||% entry[["response"]] %||% paste0("trait", i)
        model <- entry[["model"]] %||% "model"
        rep_i <- entry[["rep"]] %||% 1L
        prefix <- paste0(make.names(trait), "_", make.names(model), "_rep", rep_i)

        if (is.data.frame(entry[["eval_metrics_reps"]])) {
          write.csv(
            entry[["eval_metrics_reps"]],
            file.path(cv_dir, paste0(prefix, "_eval_metrics_reps.csv")),
            row.names = FALSE
          )
        }
        if (is.data.frame(entry[["ypred_cv_Reps_all"]])) {
          write.csv(
            entry[["ypred_cv_Reps_all"]],
            file.path(cv_dir, paste0(prefix, "_ypred_cv_Reps_all.csv")),
            row.names = FALSE
          )
        }
        if (is.data.frame(entry[["yprob_cv_Reps_all"]])) {
          write.csv(
            entry[["yprob_cv_Reps_all"]],
            file.path(cv_dir, paste0(prefix, "_yprob_cv_Reps_all.csv")),
            row.names = FALSE
          )
        }
      }
    }
  ####
    if(!is.null(res_plot)){
      save_ggplot(res_plot = res_plot,
                  plot_filename = plot_filename,
                  plot_extension = plot_extension,
                  plot_width = plot_width,
                  plot_height = plot_height,
                  plot_units = plot_units,
                  plot_dpi = plot_dpi)
    }

    if(!is.null(res_plot_mean)){
      save_ggplot(res_plot = res_plot_mean,
                  plot_filename = paste(plot_filename, "mean", sep = "_"),
                  plot_extension = plot_extension,
                  plot_width = plot_width,
                  plot_height = plot_height,
                  plot_units = plot_units,
                  plot_dpi = plot_dpi)
    }

    if(!is.null(test_diagonistic_plots)){
      diag_plot_name <- if (!is.null(res_model_output) && "met_long_data" %in% names(res_model_output)) {
        paste("MET_test_diagonistic_plots", "GS_model", sep = "_")
      } else {
        paste("test_diagonistic_plots", "GS_model", sep = "_")
      }
      if(inherits(test_diagonistic_plots, "gtable") || inherits(test_diagonistic_plots, "ggplot")){
        save_ggplot(res_plot = test_diagonistic_plots,
                    plot_filename = diag_plot_name,
                    plot_extension = plot_extension,
                    plot_width = plot_width,
                    plot_height = plot_height,
                    plot_units = plot_units,
                    plot_dpi = plot_dpi)
      } else if (is.list(test_diagonistic_plots)) {
        write_nested_plot_tables(
          test_diagonistic_plots,
          pathout = pathout,
          prefix = diag_plot_name,
          plot_extension = plot_extension,
          plot_width = plot_width,
          plot_height = plot_height,
          plot_units = plot_units,
          plot_dpi = plot_dpi
        )
      }
    }


    if(!is.null(res_plot_result_diagnostic)){

        models <-   names(res_plot_result_diagnostic$predicted_vs_observed_plots)

        for (mod in models) {

          combined_plot <- res_plot_result_diagnostic$predicted_vs_observed_plots[[mod]]

          name_plot <- paste("Cross_validation_diagonistic_plots", mod, sep = "_")

if(inherits(combined_plot, "gtable") || inherits(combined_plot, "ggplot")){
          save_ggplot(res_plot = combined_plot,
                      plot_filename = name_plot,
                      plot_extension = plot_extension,
                      plot_width = plot_width,
                      plot_height = plot_height,
                      plot_units = plot_units,
                      plot_dpi = plot_dpi)

}

        }


      }

    #######
    if(!is.null(res_plot_result_diagnostic_cv_only)){
      traits <- names(res_plot_result_diagnostic_cv_only)

      for (trait in traits) {

       res_plot_cv <-  res_plot_result_diagnostic_cv_only[[trait]]
      models <-   names(res_plot_cv$predicted_vs_observed_plots)

      for (mod in models) {

        combined_plot <- res_plot_cv$predicted_vs_observed_plots[[mod]]

        name_plot <- paste(paste("Cross_validation_diagonistic_plots",  trait, sep = "_"), mod, sep = "_")

        if(inherits(combined_plot, "gtable") || inherits(combined_plot, "ggplot")){
          save_ggplot(res_plot = combined_plot,
                      plot_filename = name_plot,
                      plot_extension = plot_extension,
                      plot_width = plot_width,
                      plot_height = plot_height,
                      plot_units = plot_units,
                      plot_dpi = plot_dpi)

        }

      }


    }

    }

    if(!is.null(feature_selected)){
      save(feature_selected, file = "feature_selected.RData")
    }
    if(!is.null(feature_score_metadata)){
      save(feature_score_metadata, file = "feature_score_metadata.RData")
    }

}

  processMMatrixModelReady <- function(pathout) {
    files_in_directory <- list.files()
    geno_omic_files <- grep("ready", files_in_directory, value = TRUE)

    if (length(geno_omic_files) > 0) {
      zipMMatrixModelReady(pathout, geno_omic_files)
      unlink(geno_omic_files)
    }
  }

  saveOutputAndZip <- function(res_model_output, res_summary_stat, output_file_name,
                               #pathout,
                               GS_model, res_plot, res_plot_mean,
                               test_diagonistic_plots, res_plot_result_diagnostic,
                               res_plot_result_diagnostic_cv_only,
                               feature_selected,
                               feature_score_metadata,
                               run_metadata,
                               cv_results_processed = NULL,
                               cv_results_raw = NULL) {
    msg <- ""
    # Check if any required object is NULL
    # if ((is.null(res_model_output) && is.null(res_summary_stat)) && (is.null(test_diagonistic_plots) && is.null(res_plot_result_diagnostic))) {
    #   cat(paste(msg,paste("The output", "from", GS_model,"is NULL.", "No output will be saved.\n")))
    #   return("Failed: Required objects are NULL")
    # }

    has_cv_outputs <- (!is.null(cv_results_processed) && length(cv_results_processed) > 0) ||
      (!is.null(cv_results_raw) && length(cv_results_raw) > 0) ||
      !is.null(res_plot_result_diagnostic_cv_only)

    if (is.null(res_model_output) && is.null(res_summary_stat) &&
        is.null(res_plot) && is.null(res_plot_mean) &&
        is.null(test_diagonistic_plots) && is.null(res_plot_result_diagnostic) &&
        !has_cv_outputs) {

      return("Failed: Required objects are NULL")
    }

    mainDir <- gp_output_root(create = TRUE)
    systime <- format(Sys.time(), "%m-%d-%Y_%I-%M%p")
    directory_stem <- gp_output_directory_stem(
      output_file_name = output_file_name,
      GS_model = GS_model
    )
    subDir <- paste(directory_stem, systime, sep = "_")
    pathout <- file.path(mainDir, subDir)
    run_index <- 1L
    while (dir.exists(pathout)) {
      run_index <- run_index + 1L
      pathout <- file.path(
        mainDir,
        paste0(subDir, "_run-", sprintf("%02d", run_index))
      )
    }
    if (!dir.create(pathout, recursive = FALSE, showWarnings = FALSE)) {
      stop("Could not create the model output directory: ", pathout, call. = FALSE)
    }
    pathout <- normalizePath(pathout, winslash = "/", mustWork = TRUE)
    # Restore the caller's working directory, not the output root: with
    # PredictProR.output_dir / PREDICTPRO_OUTPUT_DIR set, the session was left
    # in the output folder after every exported run.
    caller_wd <- getwd()
    on.exit(setwd(caller_wd), add = TRUE)
    setwd(pathout)

    saveOutput(res_model_output = res_model_output,
               res_summary_stat = res_summary_stat,
               pathout = pathout,
               GS_model = GS_model,
               res_plot = res_plot,
               res_plot_mean = res_plot_mean,
               test_diagonistic_plots = test_diagonistic_plots,
               res_plot_result_diagnostic = res_plot_result_diagnostic,
               res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
               feature_selected = feature_selected,
               feature_score_metadata = feature_score_metadata,
               run_metadata = run_metadata,
               cv_results_processed = cv_results_processed,
               cv_results_raw = cv_results_raw)
    processMMatrixModelReady(pathout)

    return(list(status = "Successful", directory = pathout))
  }

  processData <- function(GS_model,
                          res_model_output,
                          res_summary_stat,
                          res_plot,
                          #pathout,
                          res_plot_mean,
                          test_diagonistic_plots,
                          res_plot_result_diagnostic,
                          res_plot_result_diagnostic_cv_only,
                          res_mod_results_cv_per_trait_model,
                          cv_results_processed,
                          cv_results_raw,
                          system_database,
                          plot_filename,
                          plot_extension,
                          plot_width,
                          plot_height,
                          plot_units,
                          plot_dpi,
                          feature_selected,
                          feature_score_metadata,
                          run_metadata) {

      cv_results_predicted_vs_observed <- NULL
      if (!is.null(res_plot_result_diagnostic_cv_only) ||
          !is.null(res_mod_results_cv_per_trait_model)) {
        cv_results_predicted_vs_observed <- list(
          predicted_vs_observed_plots = res_plot_result_diagnostic_cv_only,
          mod_res_per_trait_per_model = res_mod_results_cv_per_trait_model
        )
      }

      output <- list(model_results = res_model_output,
                     summary_statistic = res_summary_stat,
                     run_metadata = run_metadata,
                     res_plot = res_plot,
                     cv_results_processed = cv_results_processed,
                     cv_results_raw = cv_results_raw,
                     res_plot_result_diagnostic = res_plot_result_diagnostic,
                     res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
                     cv_results_predicted_vs_observed = cv_results_predicted_vs_observed,
                     test_diagonistic_plots = test_diagonistic_plots,
                     res_mod_results_cv_per_trait_model = res_mod_results_cv_per_trait_model,
                     feature_selected = feature_selected,
                     feature_score_metadata = feature_score_metadata)

      if (isFALSE(system_database)) {
        export_result <- saveOutputAndZip(res_model_output = res_model_output,
                                                 res_summary_stat = res_summary_stat,
                                                 output_file_name = plot_filename,
                                                 #pathout = pathout,
                                                 GS_model = GS_model,
                                                 res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
                                                 res_plot = res_plot, res_plot_mean = res_plot_mean,
                                                 test_diagonistic_plots = test_diagonistic_plots,
                                                 res_plot_result_diagnostic = res_plot_result_diagnostic,
                                                 feature_selected = feature_selected,
                                                 feature_score_metadata = feature_score_metadata,
                                                 run_metadata = run_metadata,
                                                 cv_results_processed = cv_results_processed,
                                                 cv_results_raw = cv_results_raw)
        if (is.list(export_result) && !is.null(export_result[["status"]])) {
          output$export_status <- export_result[["status"]]
          output$export_directory <- export_result[["directory"]]
        } else {
          output$export_status <- export_result
        }
      }
      return(output)

  }



  return(processData(GS_model = GS_model,
                     res_model_output = res_model_output,
                     res_summary_stat = res_summary_stat,
                     res_plot = res_plot,
                     #pathout = pathout,
                     res_plot_mean = res_plot_mean,
                     test_diagonistic_plots = test_diagonistic_plots,
                     res_plot_result_diagnostic_cv_only = res_plot_result_diagnostic_cv_only,
                     res_plot_result_diagnostic = res_plot_result_diagnostic,
                     res_mod_results_cv_per_trait_model = res_mod_results_cv_per_trait_model,
                     cv_results_processed = cv_results_processed,
                     cv_results_raw = cv_results_raw,
                     system_database = system_database,
                     plot_filename = plot_filename,
                     #Plot_name_result_diagnostic = Plot_name_result_diagnostic,
                     plot_extension = plot_extension,
                     plot_width = plot_width,
                     plot_height = plot_height,
                     plot_units = plot_units,
                     plot_dpi = plot_dpi,
                     feature_selected = feature_selected,
                     feature_score_metadata = feature_score_metadata,
                     run_metadata = run_metadata))

} ## end of function

