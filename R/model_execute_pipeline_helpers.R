gp_drop_samples_by_id <- function(x, removed_ids) {
  if (is.null(x) || is.null(removed_ids)) {
    return(x)
  }

  x[!rownames(x) %in% removed_ids, , drop = FALSE]
}

gp_drop_kernel_samples_by_id <- function(x, removed_ids) {
  if (is.null(x) || is.null(removed_ids)) {
    return(x)
  }

  x[!rownames(x) %in% removed_ids, !colnames(x) %in% removed_ids, drop = FALSE]
}

gp_add_model_input <- function(model_inputs,
                               kernel_inputs,
                               result,
                               model_key,
                               kernel_key) {
  if (length(result) != 0 && all(c("kernel", "omic_model_ready") %in% names(result))) {
    kernel_inputs <- gp_add_kernel_inputs(kernel_inputs, kernel_key, result[["kernel"]])
    model_inputs[[model_key]] <- result[["omic_model_ready"]]
  } else if (length(result) != 0 && "omic_model_ready" %in% names(result)) {
    model_inputs[[model_key]] <- result[["omic_model_ready"]]
  } else if (length(result) != 0 && "kernel" %in% names(result)) {
    kernel_inputs <- gp_add_kernel_inputs(kernel_inputs, kernel_key, result[["kernel"]])
  }

  list(model_inputs = model_inputs, kernel_inputs = kernel_inputs)
}

gp_normalize_test_set <- function(test_set) {
  if (is.null(test_set)) {
    return(NULL)
  }
  if (is.list(test_set) && !is.data.frame(test_set)) {
    stop("Test_set cannot be a list. Provide it as a vector, single column dataframe or matrix.", call. = FALSE)
  }
  if (is.data.frame(test_set) || is.matrix(test_set)) {
    test_set <- test_set[, 1]
  }

  test_set <- unique(as.character(test_set))
  test_set <- test_set[!is.na(test_set)]
  if (!length(test_set)) NULL else test_set
}

gp_drop_ids_from_set <- function(ids, removed_ids) {
  ids <- gp_normalize_test_set(ids)
  removed_ids <- gp_normalize_test_set(removed_ids)
  if (is.null(ids) || is.null(removed_ids)) {
    return(ids)
  }
  ids <- setdiff(ids, removed_ids)
  if (!length(ids)) NULL else ids
}

gp_merge_test_sets <- function(...) {
  parts <- list(...)
  out <- unique(unlist(lapply(parts, gp_normalize_test_set), use.names = FALSE))
  out <- out[!is.na(out)]
  if (!length(out)) NULL else out
}

gp_test_set_overrides_inferred <- function(pheno_clean) {
  identical(pheno_clean[["test_set_source"]], "explicit")
}

gp_trait_test_set_source <- function(pheno_clean, response, gen_name) {
  if (isTRUE(gp_test_set_overrides_inferred(pheno_clean))) {
    return("explicit")
  }

  sources <- character()
  if (!is.null(gp_normalize_test_set(pheno_clean[["test_set"]]))) {
    sources <- c(sources, pheno_clean[["test_set_source"]] %||% "unknown")
  }

  test_set_by_trait <- pheno_clean[["test_set_by_trait"]]
  if (!is.null(test_set_by_trait) && !is.null(test_set_by_trait[[response]])) {
    sources <- c(sources, pheno_clean[["test_set_by_trait_source"]] %||% "inferred_missing_response")
  } else {
    pheno_data <- pheno_clean[["pheno_clean_data"]]
    if (!is.null(pheno_data) && response %in% names(pheno_data)) {
      missing_response <- is.na(pheno_data[[response]])
      if (any(missing_response, na.rm = TRUE)) {
        sources <- c(sources, "inferred_missing_response")
      }
    }
  }

  sources <- unique(as.character(sources))
  sources <- sources[!is.na(sources) & nzchar(sources)]
  if (!length(sources)) NULL else paste(sources, collapse = "+")
}

gp_trait_test_set <- function(pheno_clean, response, gen_name) {
  global_test_set <- gp_normalize_test_set(pheno_clean[["test_set"]])
  if (isTRUE(gp_test_set_overrides_inferred(pheno_clean))) {
    return(global_test_set)
  }

  test_set_by_trait <- pheno_clean[["test_set_by_trait"]]
  trait_test_set <- NULL
  if (!is.null(test_set_by_trait) && !is.null(test_set_by_trait[[response]])) {
    trait_test_set <- gp_normalize_test_set(test_set_by_trait[[response]])
  }

  if (is.null(trait_test_set)) {
    pheno_data <- pheno_clean[["pheno_clean_data"]]
    if (!is.null(pheno_data) && response %in% names(pheno_data)) {
      trait_test_set <- gp_normalize_test_set(pheno_data[[gen_name]][is.na(pheno_data[[response]])])
    }
  }

  gp_merge_test_sets(global_test_set, trait_test_set)
}

gp_trait_prediction_test_set_summary <- function(pheno_clean, response, gen_name) {
  pheno_data <- pheno_clean[["pheno_clean_data"]]
  global_test_set <- gp_normalize_test_set(pheno_clean[["test_set"]])
  test_set_by_trait <- pheno_clean[["test_set_by_trait"]]
  explicit_override <- isTRUE(gp_test_set_overrides_inferred(pheno_clean))

  rows <- lapply(response, function(trait) {
    trait_test_set <- gp_trait_test_set(pheno_clean, trait, gen_name)
    trait_specific_test_set <- NULL
    if (!is.null(test_set_by_trait) && !is.null(test_set_by_trait[[trait]])) {
      trait_specific_test_set <- gp_normalize_test_set(test_set_by_trait[[trait]])
    } else if (!is.null(pheno_data) && trait %in% names(pheno_data)) {
      trait_specific_test_set <- gp_normalize_test_set(pheno_data[[gen_name]][is.na(pheno_data[[trait]])])
    }

    missing_records <- 0L
    training_records <- 0L
    phenotype_records <- 0L
    if (!is.null(pheno_data) && trait %in% names(pheno_data) && gen_name %in% names(pheno_data)) {
      ids <- as.character(pheno_data[[gen_name]])
      phenotype_records <- length(ids)
      missing_records <- sum(is.na(pheno_data[[trait]]))
      test_ids <- gp_normalize_test_set(trait_test_set)
      training_records <- sum(!is.na(pheno_data[[trait]]) & !ids %in% (test_ids %||% character()))
    }

    data.frame(
      trait = trait,
      prediction_test_n = length(gp_normalize_test_set(trait_test_set) %||% character()),
      missing_response_test_n = length(gp_normalize_test_set(trait_specific_test_set) %||% character()),
      missing_response_records_n = as.integer(missing_records),
      global_test_n = length(global_test_set %||% character()),
      additional_global_or_genomic_test_n = length(setdiff(global_test_set %||% character(), trait_specific_test_set %||% character())),
      training_records_n = as.integer(training_records),
      phenotype_records_n = as.integer(phenotype_records),
      explicit_test_set_overrides_inferred = explicit_override,
      test_set_basis = if (isTRUE(explicit_override)) {
        "explicit_test_set_only"
      } else {
        "trait_missing_response_plus_global_or_genomic_test"
      },
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gp_attach_trait_prediction_test_info_to_cv_results <- function(cv_results,
                                                               pheno_clean,
                                                               response,
                                                               gen_name) {
  if (is.null(cv_results) || !length(cv_results)) {
    return(cv_results)
  }
  summary <- gp_trait_prediction_test_set_summary(
    pheno_clean = pheno_clean,
    response = response,
    gen_name = gen_name
  )
  if (is.null(summary) || !nrow(summary)) {
    return(cv_results)
  }

  # lapply() drops list attributes; keep gp_execution_policy so Run_metadata
  # can report the backend, workers and policy_decision_reason actually used.
  policy <- attr(cv_results, "gp_execution_policy", exact = TRUE)
  out <- lapply(cv_results, function(entry) {
    if (!is.list(entry)) {
      return(entry)
    }
    trait <- as.character(entry[["trait"]] %||% response[[1L]])
    row <- summary[summary$trait == trait, , drop = FALSE]
    if (!nrow(row)) {
      return(entry)
    }
    pred_df <- entry[["ypred_cv_Reps_all"]]
    cv_test_n <- NA_integer_
    cv_train_n <- NA_integer_
    if (is.data.frame(pred_df)) {
      if ("cv_role" %in% names(pred_df)) {
        cv_test_n <- sum(pred_df$cv_role == "test", na.rm = TRUE)
        cv_train_n <- sum(pred_df$cv_role == "train", na.rm = TRUE)
      } else if ("yhat" %in% names(pred_df)) {
        cv_test_n <- sum(!is.na(pred_df$yhat))
        cv_train_n <- sum(is.na(pred_df$yhat))
      }
    }

    cv_info <- entry[["cv_info"]]
    if (is.null(cv_info) || !is.list(cv_info)) {
      cv_info <- list()
    }
    cv_info[["n_cv_test"]] <- as.integer(cv_test_n)
    cv_info[["n_cv_train"]] <- as.integer(cv_train_n)
    cv_info[["n_test_is_cv_holdout"]] <- TRUE
    cv_info[["n_prediction_test"]] <- as.integer(row$prediction_test_n[[1L]])
    cv_info[["n_missing_response_test"]] <- as.integer(row$missing_response_test_n[[1L]])
    cv_info[["n_missing_response_records"]] <- as.integer(row$missing_response_records_n[[1L]])
    cv_info[["n_global_or_genomic_test"]] <- as.integer(row$additional_global_or_genomic_test_n[[1L]])
    cv_info[["prediction_test_set_basis"]] <- as.character(row$test_set_basis[[1L]])
    cv_info[["cv_holdout_count_basis"]] <- "cross_validation_holdout_rows_with_observed_response"
    entry[["cv_info"]] <- cv_info
    entry
  })
  attr(out, "gp_execution_policy") <- policy
  out
}

gp_ml_merged_feature_matrix <- function(ml_dat_res) {
  merged_all <- ml_dat_res[["merged_data"]][["merge_data"]]
  merged_test <- ml_dat_res[["merged_data_test"]]
  if (is.null(merged_test)) {
    return(merged_all)
  }
  missing_test_rows <- setdiff(rownames(merged_test), rownames(merged_all))
  if (length(missing_test_rows) > 0L) {
    merged_all <- rbind(merged_all, merged_test[missing_test_rows, , drop = FALSE])
  }
  merged_all
}

gp_prepare_ml_data_for_trait <- function(ml_dat_res, pheno_clean, response, gen_name) {
  if (is.null(ml_dat_res) || length(ml_dat_res) == 0) {
    return(ml_dat_res)
  }

  trait_test_set <- gp_trait_test_set(pheno_clean, response, gen_name)
  pheno_data <- pheno_clean[["pheno_clean_data"]]
  merged_all <- gp_ml_merged_feature_matrix(ml_dat_res)
  observed_train_ids <- rownames(merged_all)
  if (!is.null(pheno_data) && response %in% names(pheno_data)) {
    observed_train_ids <- unique(as.character(pheno_data[[gen_name]][!is.na(pheno_data[[response]])]))
  }

  if (is.null(trait_test_set) || !length(trait_test_set)) {
    train_ids <- intersect(rownames(merged_all), observed_train_ids)
    out <- ml_dat_res
    out[["pheno_clean_data"]] <- if (!is.null(pheno_data)) {
      pheno_data[as.character(pheno_data[[gen_name]]) %in% train_ids, , drop = FALSE]
    } else {
      pheno_data
    }
    out[["test_set"]] <- NULL
    out[["merged_data"]][["merge_data"]] <- merged_all[train_ids, , drop = FALSE]
    out[["merged_data_test"]] <- NULL
    return(out)
  }

  train_ids <- setdiff(intersect(rownames(merged_all), observed_train_ids), trait_test_set)
  pheno_train <- if (!is.null(pheno_data)) {
    pheno_data[
      as.character(pheno_data[[gen_name]]) %in% train_ids &
        !as.character(pheno_data[[gen_name]]) %in% trait_test_set,
      ,
      drop = FALSE
    ]
  } else {
    pheno_data
  }
  merged_train <- merged_all[rownames(merged_all) %in% train_ids, , drop = FALSE]
  merged_train <- merged_train[match(train_ids, rownames(merged_train)), , drop = FALSE]
  merged_test <- merged_all[rownames(merged_all) %in% trait_test_set, , drop = FALSE]
  merged_test <- merged_test[match(intersect(trait_test_set, rownames(merged_test)), rownames(merged_test)), , drop = FALSE]

  out <- ml_dat_res
  out[["pheno_clean_data"]] <- pheno_train
  out[["test_set"]] <- trait_test_set
  out[["merged_data"]] <- out[["merged_data"]]
  out[["merged_data"]][["merge_data"]] <- merged_train
  out[["merged_data_test"]] <- merged_test
  out
}

gp_prepare_pheno_for_trait <- function(pheno_clean, response, gen_name) {
  pheno_data <- pheno_clean[["pheno_clean_data"]]
  trait_test_set <- gp_trait_test_set(pheno_clean, response, gen_name)
  trait_test_set_source <- gp_trait_test_set_source(pheno_clean, response, gen_name)

  if (is.null(pheno_data) || !response %in% names(pheno_data)) {
    return(list(
      pheno_data = pheno_data,
      test_set = trait_test_set,
      test_set_source = trait_test_set_source
    ))
  }

  keep_rows <- !is.na(pheno_data[[response]])
  if (!is.null(trait_test_set) && length(trait_test_set) > 0L) {
    keep_rows <- keep_rows | as.character(pheno_data[[gen_name]]) %in% trait_test_set
  }

  list(
    pheno_data = pheno_data[keep_rows, , drop = FALSE],
    test_set = trait_test_set,
    test_set_source = trait_test_set_source
  )
}

gp_connectivity_summary_tables <- function(pheno_data) {
  connectivity_summary <- attr(pheno_data, "connectivity_summary")
  if (is.null(connectivity_summary) || !is.list(connectivity_summary)) {
    return(NULL)
  }

  overall_fields <- setdiff(names(connectivity_summary), "genotype_count_by_env")
  overall_values <- lapply(connectivity_summary[overall_fields], function(x) {
    if (length(x) == 0 || is.null(x)) {
      return(NA)
    }
    if (length(x) > 1) {
      return(paste(x, collapse = ","))
    }
    x
  })

  overall_table <- data.frame(
    metric = overall_fields,
    value = unlist(overall_values, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  env_counts <- connectivity_summary[["genotype_count_by_env"]]
  env_table <- NULL
  if (!is.null(env_counts) && length(env_counts) > 0) {
    env_table <- data.frame(
      env = names(env_counts),
      genotype_count = as.integer(env_counts),
      stringsAsFactors = FALSE
    )
  }

  list(
    Connectivity_summary = overall_table,
    Connectivity_genotype_count_by_env = env_table
  )
}

gp_is_multi_environment_pheno <- function(pheno_data,
                                          gen_name,
                                          heter_groups = NULL,
                                          response = NULL,
                                          response_family = NULL) {
  if (is.null(pheno_data) || is.null(gen_name) || !gen_name %in% names(pheno_data)) {
    return(FALSE)
  }

  gp_is_multi_environment_trait_panel(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response,
    response_family = response_family
  )
}

gp_normalize_single_environment_heter_controls <- function(pheno_data,
                                                           gen_name,
                                                           heter_groups = NULL,
                                                           heter_resid = FALSE,
                                                           var_cov_str = NULL,
                                                           response = NULL,
                                                           response_family = NULL,
                                                           preserve_heter_groups = FALSE,
                                                           preserve_var_cov_str = FALSE) {
  out <- list(
    heter_groups = heter_groups,
    heter_resid = isTRUE(heter_resid),
    var_cov_str = var_cov_str,
    single_environment = FALSE,
    changed = FALSE
  )
  if (is.null(pheno_data) || is.null(gen_name) || !gen_name %in% names(pheno_data)) {
    return(out)
  }

  single_environment <- !gp_is_multi_environment_trait_panel(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response,
    response_family = response_family
  )
  out$single_environment <- single_environment
  if (!isTRUE(single_environment)) {
    return(out)
  }

  out$heter_groups <- NULL
  out$heter_resid <- FALSE
  if (!isTRUE(preserve_var_cov_str)) {
    out$var_cov_str <- NULL
  }
  out$changed <- !is.null(heter_groups) || isTRUE(heter_resid) ||
    (!isTRUE(preserve_var_cov_str) && !is.null(var_cov_str))
  out
}

gp_attach_connectivity_outputs <- function(res_model_output, res_summary_stat, pheno_data) {
  connectivity_tables <- gp_connectivity_summary_tables(pheno_data)
  if (is.null(connectivity_tables)) {
    return(list(
      res_model_output = res_model_output,
      res_summary_stat = res_summary_stat
    ))
  }

  if (is.null(res_model_output)) {
    res_model_output <- list()
  }
  if (is.null(res_summary_stat)) {
    res_summary_stat <- list()
  }

  for (nm in names(connectivity_tables)) {
    table_obj <- connectivity_tables[[nm]]
    if (!is.null(table_obj)) {
      res_model_output[[nm]] <- table_obj
      res_summary_stat[[nm]] <- table_obj
    }
  }

  list(
    res_model_output = res_model_output,
    res_summary_stat = res_summary_stat
  )
}

gp_standard_prediction_columns <- function() {
  # `Genetic_variance` is reserved for engines that estimate a genetic random
  # effect. Predictive ML/DL routes keep the phenotypic scale in
  # `Prediction_stability_reference_variance`. The public `Reliability` field
  # is a plotting-compatible alias of that descriptive score for ML/DL; it is
  # not genetic reliability, and `Genetic_variance` remains unavailable.
  c(
    "GID",
    "Predicted_value",
    "Train_Test_Label",
    "Standard_error",
    "PEV",
    "PEV_basis",
    "Prediction_error_variance",
    "Prediction_uncertainty_source",
    "Prediction_interval_method",
    "Prediction_interval_nominal_coverage",
    "Prediction_interval_calibration_n",
    "Genetic_variance",
    "Reliability_reference_variance",
    "Reliability",
    "Reliability_remarks",
    "Reliability_variance_input",
    "Reliability_percentage",
    "Reliability_basis",
    "Prediction_stability",
    "Prediction_stability_remarks",
    "Prediction_stability_reference_variance",
    "Observed_value",
    "Train_Test"
  )
}

gp_standardize_prediction_table <- function(x, gen_name = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) {
    return(x)
  }

  gid_candidates <- unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
  gid_col <- gid_candidates[gid_candidates %in% names(df)][1L]
  if (!is.na(gid_col)) {
    df[["GID"]] <- as.character(df[[gid_col]])
  }

  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(df))) {
    if (!"Train_Test_Label" %in% names(df) && "Train_Test" %in% names(df)) {
      df[["Train_Test_Label"]] <- as.character(df[["Train_Test"]])
    }
    return(df)
  }

  pred_candidates <- c(
    "Predicted_value",
    "Prediction",
    "prediction",
    "yhat",
    "predicted.value",
    "predicted_value",
    "BLUP"
  )
  pred_col <- pred_candidates[pred_candidates %in% names(df)][1L]
  if (!is.na(pred_col)) {
    df[["Predicted_value"]] <- suppressWarnings(as.numeric(df[[pred_col]]))
  }

  se_candidates <- c("Standard_error", "std.error", "std_error", "SE", "Prediction_SE", "Prediction_SE_latent", "Prediction_SE_observed")
  se_col <- se_candidates[se_candidates %in% names(df)][1L]
  if (!is.na(se_col)) {
    df[["Standard_error"]] <- suppressWarnings(as.numeric(df[[se_col]]))
  }

  pev_candidates <- c("PEV", "Prediction_error_variance", "Prediction_Var", "Prediction_Var_latent", "Prediction_Var_observed")
  pev_col <- pev_candidates[pev_candidates %in% names(df)][1L]
  if (!is.na(pev_col)) {
    df[["PEV"]] <- suppressWarnings(as.numeric(df[[pev_col]]))
  }
  if (!"PEV" %in% names(df) && "Standard_error" %in% names(df)) {
    df[["PEV"]] <- suppressWarnings(as.numeric(df[["Standard_error"]]))^2
  }
  if ("PEV" %in% names(df) && !"Prediction_error_variance" %in% names(df)) {
    df[["Prediction_error_variance"]] <- df[["PEV"]]
  }
  if ("Prediction_error_variance" %in% names(df) && !"PEV" %in% names(df)) {
    df[["PEV"]] <- suppressWarnings(as.numeric(df[["Prediction_error_variance"]]))
  }

  if (!"Train_Test_Label" %in% names(df) && "Train_Test" %in% names(df)) {
    df[["Train_Test_Label"]] <- as.character(df[["Train_Test"]])
  }
  if ("Train_Test_Label" %in% names(df) && !"Train_Test" %in% names(df)) {
    df[["Train_Test"]] <- as.character(df[["Train_Test_Label"]])
  }

  standard_cols <- gp_standard_prediction_columns()
  for (nm in standard_cols) {
    if (!nm %in% names(df)) {
      df[[nm]] <- NA
    }
  }
  if (exists("gp_contract_recompute_flat_reliability", mode = "function")) {
    df <- gp_contract_recompute_flat_reliability(df)
  }
  df[, c(standard_cols, setdiff(names(df), standard_cols)), drop = FALSE]
}

gp_standardize_model_prediction_outputs <- function(res_model_output, gen_name = NULL) {
  if (is.null(res_model_output) || !is.list(res_model_output)) {
    return(res_model_output)
  }

  if (is.data.frame(res_model_output[["predicted_values"]])) {
    res_model_output[["predicted_values"]] <- gp_standardize_prediction_table(
      res_model_output[["predicted_values"]],
      gen_name = gen_name
    )
    if (is.null(res_model_output[["Predicted_value"]])) {
      res_model_output[["Predicted_value"]] <- res_model_output[["predicted_values"]]
    }
  }
  if (is.data.frame(res_model_output[["Predicted_value"]])) {
    normalized_pred <- gp_standardize_prediction_table(
      res_model_output[["Predicted_value"]],
      gen_name = gen_name
    )
    res_model_output[["Predicted_value"]] <- normalized_pred
    if (is.null(res_model_output[["predicted_values"]])) {
      res_model_output[["predicted_values"]] <- normalized_pred
    }
  }
  if (is.null(res_model_output[["diagnostic_plots"]]) &&
      !is.null(res_model_output[["diagnostic_tst_plot"]])) {
    res_model_output[["diagnostic_plots"]] <- res_model_output[["diagnostic_tst_plot"]]
  }

  nested_result_names <- c("bayes_result", "asreml_result", "gp_result")
  for (nm in nested_result_names) {
    nested <- res_model_output[[nm]]
    if (!is.list(nested)) {
      next
    }
    if (is.data.frame(nested[["predicted_values"]])) {
      nested[["predicted_values"]] <- gp_standardize_prediction_table(
        nested[["predicted_values"]],
        gen_name = gen_name
      )
    }
    if (is.data.frame(nested[["Predicted_value"]])) {
      normalized_nested <- gp_standardize_prediction_table(
        nested[["Predicted_value"]],
        gen_name = gen_name
      )
      if (is.null(nested[["predicted_values"]])) {
        nested[["predicted_values"]] <- normalized_nested
      }
      if (is.null(res_model_output[["predicted_values"]])) {
        res_model_output[["predicted_values"]] <- normalized_nested
      }
    }
    if (is.data.frame(nested[["predicted_values"]]) && is.null(nested[["Predicted_value"]])) {
      nested[["Predicted_value"]] <- nested[["predicted_values"]]
    }
    if (is.null(nested[["diagnostic_plots"]]) &&
        !is.null(nested[["diagnostic_tst_plot"]])) {
      nested[["diagnostic_plots"]] <- nested[["diagnostic_tst_plot"]]
    }
    if (is.null(res_model_output[["diagnostic_plots"]]) &&
        !is.null(nested[["diagnostic_plots"]])) {
      res_model_output[["diagnostic_plots"]] <- nested[["diagnostic_plots"]]
    }
    res_model_output[[nm]] <- nested
  }

  res_model_output
}

gp_public_prediction_columns <- function(response_family = "gaussian",
                                          include_env = FALSE,
                                          include_trait = FALSE) {
  supplied_family <- tolower(as.character(response_family %||% "gaussian")[1L])
  fam <- if (identical(supplied_family, "classification")) {
    "classification"
  } else {
    gp_normalize_response_family(supplied_family)
  }
  if (fam %in% c("binary", "multiclass", "ordinal") || identical(fam, "classification")) {
    return(gp_classification_prediction_columns(
      include_env = include_env,
      include_trait = include_trait
    ))
  }
  gp_gaussian_prediction_columns(include_env = include_env, include_trait = include_trait)
}

gp_public_prediction_table <- function(x, gen_name = NULL, heter_groups = NULL) {
  if (is.null(x)) {
    return(gp_format_gaussian_prediction_table(NULL, gen_name = gen_name, heter_groups = heter_groups))
  }
  if (gp_is_classification_prediction_table(x)) {
    return(gp_format_classification_prediction_table(x, gen_name = gen_name, heter_groups = heter_groups))
  }
  gp_format_gaussian_prediction_table(x, gen_name = gen_name, heter_groups = heter_groups)
}

gp_public_model_parameters <- function(x) {
  if (is.data.frame(x)) {
    out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    out <- data.frame(stat = character(), summary = character(), stringsAsFactors = FALSE)
  }
  if (!"stat" %in% names(out)) {
    out[["stat"]] <- character(nrow(out))
  }
  if (!"summary" %in% names(out)) {
    out[["summary"]] <- character(nrow(out))
  }
  out[, c("stat", "summary", setdiff(names(out), c("stat", "summary"))), drop = FALSE]
}

gp_find_first_data_frame <- function(x, candidates) {
  for (nm in candidates) {
    obj <- x[[nm]]
    if (is.data.frame(obj)) {
      return(obj)
    }
    if (is.list(obj)) {
      nested <- gp_find_first_data_frame(obj, candidates)
      if (is.data.frame(nested)) {
        return(nested)
      }
    }
  }
  NULL
}

gp_public_variance_components <- function(x,
                                          heter_groups = NULL,
                                          model = NULL,
                                          response_family = NULL) {
  resolved_model <- if (exists("gp_vc_resolve_model", mode = "function")) {
    gp_vc_resolve_model(x, model = model)
  } else {
    model
  }
  variance_role <- if (exists("gp_model_variance_role", mode = "function")) {
    gp_model_variance_role(resolved_model)
  } else {
    "unknown"
  }
  is_gp_result <- exists("gp_vc_is_gp_result", mode = "function") &&
    isTRUE(gp_vc_is_gp_result(x))
  gp_method <- if (isTRUE(is_gp_result) && exists("gp_vc_gp_estimation_method", mode = "function")) {
    gp_vc_gp_estimation_method(x)
  } else {
    "model_reported"
  }
  variance_status <- if (exists("gp_vc_find_scalar", mode = "function")) {
    gp_vc_find_scalar(x, c("variance_component_status"))
  } else {
    NULL
  }
  variance_status_token <- as.character(variance_status %||% character())
  if (length(variance_status_token) &&
      identical(variance_status_token[[1L]], "not_promoted_nonconverged")) {
    trait_labels <- if (exists("gp_find_named_object", mode = "function")) {
      gp_find_named_object(x, c("trait_levels", "traits", "response"))
    } else {
      NULL
    }
    trait_labels <- unique(as.character(trait_labels %||% character()))
    trait_labels <- trait_labels[!is.na(trait_labels) & nzchar(trait_labels)]
    if (length(trait_labels) >= 2L) {
      unavailable <- do.call(rbind, lapply(trait_labels, function(trait) {
        data.frame(
          Trait = trait,
          Component = c(
            "genetic_variance_not_promoted_nonconverged",
            "residual_variance_not_promoted_nonconverged",
            "heritability_not_promoted_nonconverged"
          ),
          Components = NA_real_,
          Standard_error = NA_real_,
          stringsAsFactors = FALSE
        )
      }))
      return(gp_format_variance_component_columns(unavailable, include_trait = TRUE))
    }
    return(gp_unavailable_variance_components(
      reason = "variance_components_not_promoted_nonconverged",
      estimation_method = gp_method
    ))
  }
  # Prefer an explicitly named public variance-component table before trying
  # GP-specific structural formatters.  ASReml and joint Bayesian routes place
  # their already-labelled Trait-level table directly in the result bundle;
  # treating that general list as a GP payload can collapse it to the first
  # trait and make reliability provenance inconsistent.
  if (is.list(x)) {
    explicit_vc <- x[["variance_components"]] %||% x[["Variance_components"]]
    if (is.null(explicit_vc)) {
      for (nm in c("bayes_result", "asreml_result")) {
        nested <- x[[nm]]
        if (is.list(nested)) {
          explicit_vc <- nested[["variance_components"]] %||%
            nested[["Variance_components"]]
          if (!is.null(explicit_vc)) break
        }
      }
    }
    # Only short-circuit for an already-labelled public component table. Raw
    # GP varcomp payloads still need the structural GP formatters below, while
    # ASReml/Bayesian tables with Component/Components are authoritative even
    # if a complex surrounding bundle happens to look GP-like.
    if (is.data.frame(explicit_vc) &&
        all(c("Component", "Components") %in% names(explicit_vc))) {
      explicit_vc <- gp_format_variance_components(
        explicit_vc,
        default_estimation_method = gp_method
      )
      if (is.data.frame(explicit_vc) && nrow(explicit_vc) > 0L) {
        return(explicit_vc)
      }
    }
  }
  # A joint multi-trait subprocess fit may expose response-scale Sigma_G and
  # Sigma_eps only in its private fit payload.  Handle that specific, labeled
  # covariance case before the single-environment formatter can emit its
  # strict "varcomp not returned" placeholder.
  if (isTRUE(is_gp_result) &&
      exists("gp_variance_components_from_covariances", mode = "function")) {
    joint_cov_vc <- gp_variance_components_from_covariances(x)
    joint_traits <- if (is.data.frame(joint_cov_vc) && "Trait" %in% names(joint_cov_vc)) {
      unique(as.character(joint_cov_vc[["Trait"]]))
    } else {
      character()
    }
    joint_traits <- joint_traits[!is.na(joint_traits) & nzchar(joint_traits)]
    joint_components <- if (is.data.frame(joint_cov_vc) && "Component" %in% names(joint_cov_vc)) {
      as.character(joint_cov_vc[["Component"]])
    } else {
      character()
    }
    if (length(joint_traits) >= 2L &&
        all(c("genetic_variance", "residual_variance") %in% joint_components)) {
      joint_cov_vc <- if (exists("gp_vc_apply_estimation_method", mode = "function")) {
        gp_vc_apply_estimation_method(joint_cov_vc, gp_method)
      } else {
        joint_cov_vc
      }
      return(gp_format_variance_component_columns(joint_cov_vc, include_trait = TRUE))
    }
  }
  if (exists("gp_format_gp_multi_env_variance_components", mode = "function")) {
    # Phase 3.27: thread heter_groups so MET detection works when the user
    # supplied a custom env col name (e.g. "YYY").
    gp_vc <- gp_format_gp_multi_env_variance_components(x, heter_groups = heter_groups)
    if (is.data.frame(gp_vc) && nrow(gp_vc) > 0L) {
      # Append per-env Vg + h2 for MET tables (factor-analytic and compound-
      # symmetry layouts) so all GP models report per-env (cross-engine consistency).
      # Pass fa_derived through so the FA branch can use its delta-method SEs.
      fa_derived <- if (exists("gp_vc_extract_fa_derived", mode = "function")) {
        gp_vc_extract_fa_derived(x)
      } else NULL
      if (exists("gp_met_varcomp_augment_per_env", mode = "function")) {
        gp_vc <- gp_met_varcomp_augment_per_env(gp_vc, fa_derived = fa_derived)
      }
      # Phase 3.29: the per-env rows appended above inherit NA in the optional
      # `Boundary` flag column (it is only set on the raw varcomp rows). Per-env
      # derived rows are not boundary parameters, so normalise their flag to "".
      if ("Boundary" %in% names(gp_vc)) {
        gp_vc[["Boundary"]][is.na(gp_vc[["Boundary"]])] <- ""
      }
      # Rewrite the internal "Env" token in varcomp labels to the user's
      # heter_groups name (must run AFTER the per-env augment, which keys off
      # the "Env_" prefix).
      if (exists("gp_vc_relabel_env_to_heter_groups", mode = "function")) {
        gp_vc <- gp_vc_relabel_env_to_heter_groups(gp_vc, heter_groups)
      }
      return(gp_vc)
    }
  }
  if (exists("gp_format_gp_single_env_variance_components", mode = "function")) {
    gp_vc <- gp_format_gp_single_env_variance_components(x, heter_groups = heter_groups)
    if (is.data.frame(gp_vc) && nrow(gp_vc) > 0L) {
      return(gp_vc)
    }
  }
  if (isTRUE(is_gp_result) && identical(gp_method, "REML")) {
    return(if (exists("gp_unavailable_variance_components", mode = "function")) {
      gp_unavailable_variance_components(
        reason = "reml_variance_components_not_returned",
        estimation_method = "REML"
      )
    } else {
      gp_empty_variance_components()
    })
  }
  vc <- gp_format_variance_components(x, default_estimation_method = gp_method)
  if (is.data.frame(vc) && nrow(vc) > 0L) {
    return(vc)
  }
  if (isTRUE(is_gp_result)) {
    return(if (exists("gp_unavailable_variance_components", mode = "function")) {
      gp_unavailable_variance_components(
        reason = "variance_components_not_estimated",
        estimation_method = "not_estimated"
      )
    } else {
      gp_empty_variance_components()
    })
  }
  if (identical(variance_role, "model_based")) {
    return(if (exists("gp_unavailable_variance_components", mode = "function")) {
      gp_unavailable_variance_components(
        reason = "model_variance_components_not_returned",
        estimation_method = "model_reported"
      )
    } else {
      gp_empty_variance_components()
    })
  }
  if (identical(variance_role, "predictive_only") &&
      exists("gp_predictive_model_variance_components", mode = "function")) {
    pred <- if (exists("gp_vc_extract_prediction_table", mode = "function")) {
      gp_vc_extract_prediction_table(x)
    } else if (is.list(x)) {
      x[["predicted_values"]] %||% x[["Predicted_value"]]
    } else {
      NULL
    }
    resolved_family <- if (exists("gp_vc_resolve_response_family", mode = "function")) {
      gp_vc_resolve_response_family(x, response_family = response_family)
    } else {
      response_family %||% "auto"
    }
    pred_vc <- gp_predictive_model_variance_components(
      pred,
      response_family = resolved_family
    )
    if (is.data.frame(pred_vc) && nrow(pred_vc) > 0L) {
      return(pred_vc)
    }
  }
  vc
}

gp_public_diagnostic_plots <- function(x) {
  empty <- list(
    predicted_vs_observed = NULL,
    predicted_vs_reliability = NULL,
    prediction_interval = NULL,
    reliability_summary = NULL
  )
  plots <- if (is.list(x)) x[["diagnostic_plots"]] %||% x[["diagnostic_tst_plot"]] else NULL
  if (is.null(plots) && is.list(x)) {
    for (nm in c("bayes_result", "asreml_result", "gp_result")) {
      nested <- x[[nm]]
      if (is.list(nested)) {
        plots <- nested[["diagnostic_plots"]] %||% nested[["diagnostic_tst_plot"]]
        if (!is.null(plots)) {
          break
        }
      }
    }
  }
  if (inherits(plots, "ggplot") || inherits(plots, "gtable")) {
    empty[["predicted_vs_observed"]] <- plots
    return(empty)
  }
  if (!is.list(plots)) {
    return(empty)
  }
  first_non_null <- function(...) {
    vals <- list(...)
    for (val in vals) {
      if (!is.null(val)) {
        return(val)
      }
    }
    NULL
  }
  empty["predicted_vs_observed"] <- list(first_non_null(
    plots[["predicted_vs_observed"]],
    plots[["predicted_observed"]],
    plots[["diagnostic_tst_plot"]],
    plots[["predicted"]]
  ))
  empty["predicted_vs_reliability"] <- list(first_non_null(
    plots[["predicted_vs_reliability"]],
    plots[["reliability_plot"]]
  ))
  empty["prediction_interval"] <- list(first_non_null(
    plots[["prediction_interval"]],
    plots[["prediction_inter_vs_reliability"]],
    plots[["prediction_intervals"]]
  ))
  empty["reliability_summary"] <- list(first_non_null(
    plots[["reliability_summary"]],
    plots[["predicted_vs_composite_reliability"]]
  ))
  if (is.null(empty[["predicted_vs_observed"]]) &&
      !is.null(plots[["prediction_inter_vs_reliability"]]) &&
      exists("gp_arrange_grob_safely", mode = "function") &&
      requireNamespace("gridExtra", quietly = TRUE)) {
    plot_pieces <- Filter(
      function(p) inherits(p, "ggplot") || inherits(p, "gtable") || inherits(p, "grob") || inherits(p, "gTree"),
      list(
        empty[["predicted_vs_reliability"]],
        empty[["prediction_interval"]],
        empty[["reliability_summary"]]
      )
    )
    if (length(plot_pieces) > 0L) {
      empty <- list(
        predicted_vs_observed = do.call(
          gp_arrange_grob_safely,
          c(plot_pieces, list(ncol = min(2L, length(plot_pieces))))
        ),
        predicted_vs_reliability = NULL,
        prediction_interval = NULL,
        reliability_summary = NULL
      )
    }
  }
  empty
}

gp_standardize_public_model_result <- function(res_model_output,
                                               gen_name = NULL,
                                               heter_groups = NULL,
                                               pheno_data = NULL,
                                               response = NULL,
                                               model = NULL,
                                               response_family = NULL) {
  if (is.null(res_model_output) || !is.list(res_model_output)) {
    res_model_output <- list()
  }
  res_model_output <- gp_standardize_model_prediction_outputs(res_model_output, gen_name = gen_name)
  pred <- res_model_output[["predicted_values"]] %||% res_model_output[["Predicted_value"]]
  if (is.null(pred)) {
    for (nm in c("bayes_result", "asreml_result", "gp_result")) {
      nested <- res_model_output[[nm]]
      if (is.list(nested)) {
        pred <- nested[["predicted_values"]] %||% nested[["Predicted_value"]]
        if (!is.null(pred)) {
          break
        }
      }
    }
  }
  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      exists("gp_merge_classification_observed_values", mode = "function") &&
      gp_is_classification_prediction_table(pred)) {
    residual <- res_model_output[["Residual_value"]]
    if (is.null(residual)) {
      for (nm in c("bayes_result", "asreml_result", "gp_result")) {
        nested <- res_model_output[[nm]]
        if (is.list(nested) && !is.null(nested[["Residual_value"]])) {
          residual <- nested[["Residual_value"]]
          break
        }
      }
    }
    pred <- gp_merge_classification_observed_values(
      pred,
      residual = residual,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  }
  vc <- gp_public_variance_components(
    res_model_output,
    heter_groups = heter_groups,
    model = model,
    response_family = response_family
  )
  # warn once per model_execute() call (this standardisation can run more than
  # once on the same result)
  if (is.data.frame(vc) && "Boundary" %in% names(vc) &&
      any(vc[["Boundary"]] %in% "not_estimable")) {
    where <- as.character(vc[["Component"]][vc[["Boundary"]] %in% "not_estimable" &
                                              grepl("residual", vc[["Component"]], ignore.case = TRUE)])
    warn_key <- paste(model %||% "", paste(unique(where), collapse = ","))
    warned <- get0("warned", envir = PredictProR_runtime_cache, inherits = FALSE) %||% character()
    if (!warn_key %in% warned) {
      assign("warned", c(warned, warn_key), envir = PredictProR_runtime_cache)
      warning(
        model %||% "Model", ": the residual variance was estimated at zero (",
        paste(unique(where), collapse = ", "), "). Genetic and residual variance cannot be separated there, ",
        "so genetic variance, residual variance, heritability and reliability are reported as not estimable ",
        "(Boundary = \"not_estimable\"). Predictions are returned, but their standard errors and intervals ",
        "are not reliable.",
        call. = FALSE
      )
    }
  }
  if (exists("gp_standardize_multitrait_prediction_contract", mode = "function") &&
      isTRUE(gp_is_joint_multitrait_gaussian_prediction_table(pred))) {
    pred <- gp_standardize_multitrait_prediction_contract(pred, vc = vc, gen_name = gen_name)
    res_model_output[["predicted_values"]] <- pred
    res_model_output[["Predicted_value"]] <- pred
  }
  pred <- gp_public_prediction_table(pred, gen_name = gen_name, heter_groups = heter_groups)
  if (exists("gp_is_public_gaussian_table", mode = "function") &&
      exists("gp_merge_gaussian_observed_values", mode = "function") &&
      gp_is_public_gaussian_table(pred)) {
    pred <- gp_merge_gaussian_observed_values(
      pred,
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  }
  total_pred <- res_model_output[["across_environment_predicted_values"]] %||%
    res_model_output[["Total_Predicted_value"]]
  if (is.null(total_pred)) {
    for (nm in c("bayes_result", "asreml_result", "gp_result")) {
      nested <- res_model_output[[nm]]
      nested_total <- if (is.list(nested)) {
        nested[["across_environment_predicted_values"]] %||%
          nested[["Total_Predicted_value"]]
      } else {
        NULL
      }
      if (!is.null(nested_total)) {
        total_pred <- nested_total
        break
      }
    }
  }
  total_pred <- if (!is.null(total_pred)) {
    if (exists("gp_standardize_multitrait_prediction_contract", mode = "function") &&
        isTRUE(gp_is_joint_multitrait_gaussian_prediction_table(total_pred))) {
      total_pred <- gp_standardize_multitrait_prediction_contract(total_pred, vc = vc, gen_name = gen_name)
    }
    gp_format_total_prediction_table(total_pred, gen_name = gen_name, heter_groups = heter_groups)
  } else {
    NULL
  }
  plots <- gp_public_diagnostic_plots(res_model_output)
  if (exists("gp_is_public_classification_table", mode = "function") &&
      exists("gp_classification_diagnostic_plots", mode = "function") &&
      gp_is_public_classification_table(pred) &&
      !any(vapply(plots, function(x) !is.null(x), logical(1L)))) {
    plots <- gp_classification_diagnostic_plots(pred)
  }
  if (exists("gp_is_public_gaussian_table", mode = "function") &&
      exists("gp_gaussian_diagnostic_plots", mode = "function") &&
      gp_is_public_gaussian_table(pred) &&
      !any(vapply(plots, function(x) !is.null(x), logical(1L)))) {
    plots <- gp_gaussian_diagnostic_plots(pred)
  }
  model_parameters <- res_model_output[["model_parameters"]]
  if (!is.data.frame(model_parameters)) {
    for (nm in c("bayes_result", "asreml_result", "gp_result")) {
      nested <- res_model_output[[nm]]
      if (is.list(nested) && is.data.frame(nested[["model_parameters"]])) {
        model_parameters <- nested[["model_parameters"]]
        break
      }
    }
  }
  public <- list(
    model_parameters = gp_public_model_parameters(model_parameters),
    predicted_values = pred,
    diagnostic_plots = plots,
    variance_components = vc
  )
  for (nm in c(
    "multitrait_prediction_wide", "multitrait_trait_counts",
    "met_feature_summary", "variance_component_intervals",
    "kernel_variance_components", "Genetic_covariance_by_kernel",
    "genetic_covariance_by_kernel_long", "kernel_combination",
    "variance_component_status", "fa_parameters",
    "asreml_convergence_trace", "model_notes",
    "Residual_variance_by_environment",
    "Residual_variance_by_observation",
    "dl_seed_manifest", "dl_computation_plan",
    "dl_seed_predictions", "dl_seed_variability"
  )) {
    value <- res_model_output[[nm]]
    if (is.null(value)) {
      for (parent in c("bayes_result", "asreml_result", "gp_result")) {
        nested <- res_model_output[[parent]]
        if (is.list(nested) && !is.null(nested[[nm]])) {
          value <- nested[[nm]]
          break
        }
      }
    }
    if (!is.null(value)) {
      public[[nm]] <- value
    }
  }
  for (nm in c("feature_selection_metadata", "feature_score_metadata")) {
    if (!is.null(res_model_output[[nm]])) {
      public[[nm]] <- res_model_output[[nm]]
    }
  }
  if (!is.null(total_pred)) {
    public[["across_environment_predicted_values"]] <- total_pred
  }
  public <- gp_attach_correlation_matrices(public, res_model_output)
  matrix_se_names <- c(
    "Genetic_covariance_traits_SE", "Genetic_correlation_traits_SE",
    "Residual_covariance_traits_SE", "Residual_correlation_traits_SE",
    "Total_genetic_covariance_traits_SE", "Total_genetic_correlation_traits_SE",
    "Genetic_covariance_environments_SE", "Genetic_correlation_environments_SE",
    "Residual_covariance_environments_SE", "Residual_correlation_environments_SE",
    "Prediction_covariance_traits_SE", "Prediction_correlation_traits_SE",
    "Prediction_error_covariance_traits_SE", "Prediction_error_correlation_traits_SE",
    "Prediction_covariance_environments_SE", "Prediction_correlation_environments_SE",
    "Prediction_error_covariance_environments_SE", "Prediction_error_correlation_environments_SE",
    "Genetic_covariance_trait_environment_SE",
    "Genetic_correlation_trait_environment_SE"
  )
  for (nm in matrix_se_names) {
    value <- gp_find_named_object(res_model_output, nm)
    if (is.matrix(value)) {
      public[[nm]] <- value
    }
  }
  if (!is.null(res_model_output[["prediction_uncertainty_source"]])) {
    public[["prediction_uncertainty_source"]] <- res_model_output[["prediction_uncertainty_source"]]
  }
  if (!is.null(res_model_output[["predictive_covariance_basis"]])) {
    public[["predictive_covariance_basis"]] <- res_model_output[["predictive_covariance_basis"]]
  }
  # Reconcile the prediction table once more against the final public
  # variance-component table. Complex Bayesian and ASReml bundles can resolve
  # their labelled components only after nested objects have been flattened;
  # this pass guarantees the file-writing object carries the same Trait-level
  # reliability reference and formula as the exported component table.
  if (exists("gp_standardize_multitrait_prediction_contract", mode = "function") &&
      isTRUE(gp_is_joint_multitrait_gaussian_prediction_table(public[["predicted_values"]]))) {
    public[["predicted_values"]] <- gp_standardize_multitrait_prediction_contract(
      public[["predicted_values"]],
      vc = public[["variance_components"]],
      gen_name = gen_name
    )
    public[["predicted_values"]] <- gp_public_prediction_table(
      public[["predicted_values"]],
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  }
  class(public) <- unique(c("predictpror_model_result", class(public)))
  validate_prediction_output(
    public,
    response_family = "auto",
    gen_name = gen_name,
    heter_groups = heter_groups,
    strict = TRUE
  )
  public
}

gp_public_index_model_result <- function(model_result,
                                         gen_name = NULL,
                                         heter_groups = NULL,
                                         model = NULL,
                                         response_family = NULL) {
  if (is.null(model_result) || !is.list(model_result)) {
    return(gp_standardize_public_model_result(
      model_result,
      gen_name = gen_name,
      heter_groups = heter_groups,
      model = model,
      response_family = response_family
    ))
  }
  pred <- model_result[["predicted_values"]] %||% model_result[["Predicted_value"]]
  if (exists("gp_is_hybrid_prediction_table", mode = "function") &&
      is.data.frame(pred) &&
      isTRUE(gp_is_hybrid_prediction_table(pred)) &&
      exists("gp_standardize_hybrid_model_result", mode = "function")) {
    return(gp_standardize_hybrid_model_result(
      model_result,
      model = model,
      response_family = response_family
    ))
  }
  gp_standardize_public_model_result(
    model_result,
    gen_name = gen_name,
    heter_groups = heter_groups,
    model = model,
    response_family = response_family
  )
}

gp_index_specialized_model_execute_result <- function(handled_result,
                                                      model_label,
                                                      trait_name,
                                                      gen_name = NULL,
                                                      heter_groups = NULL,
                                                      run_metadata = NULL,
                                                      feature_score_metadata = NULL) {
  model_label <- as.character(model_label %||% "model")[1L]
  if (is.na(model_label) || !nzchar(model_label)) {
    model_label <- "model"
  }
  trait_name <- as.character(trait_name %||% "trait")[1L]
  if (is.na(trait_name) || !nzchar(trait_name)) {
    trait_name <- "trait"
  }
  # Phase 3.26: store under the user-supplied model name only. R's [[..]] and
  # backtick-`$name` accessors handle dashes; Phase 3.25's dual-storage
  # duplicated every entry in names(by_model) and is the wrong fix.
  list_key <- as.character(model_label)
  model_key <- make.names(as.character(model_label))  # plot file names only
  model_result <- if (is.list(handled_result)) handled_result[["model_results"]] else NULL
  model_result <- gp_public_index_model_result(
    model_result,
    gen_name = gen_name,
    heter_groups = heter_groups,
    model = model_label
  )

  by_trait <- list()
  by_trait[[trait_name]] <- list()
  by_trait[[trait_name]][[list_key]] <- model_result

  by_model <- list()
  by_model[[list_key]] <- list()
  by_model[[list_key]][[trait_name]] <- model_result

  handled <- if (is.list(handled_result)) handled_result else list()
  out <- handled
  out[["cv_results_raw"]] <- handled[["cv_results_raw"]] %||% NULL
  out[["cv_results_processed"]] <- handled[["cv_results_processed"]] %||% NULL
  out[["cv_results_predicted_vs_observed"]] <-
    handled[["cv_results_predicted_vs_observed"]] %||%
    handled[["res_plot_result_diagnostic"]] %||%
    NULL
  out[["model_results_by_trait"]] <- by_trait
  out[["model_results_by_model"]] <- by_model
  out[["model_results"]] <- model_result
  out[["run_metadata"]] <- run_metadata %||% handled[["run_metadata"]] %||% NULL
  out[["feature_score_metadata"]] <-
    feature_score_metadata %||% handled[["feature_score_metadata"]] %||% NULL
  out
}

gp_prepare_model_input_objects <- function(ctx) {
  ctx <- gp_apply_model_input_guardrails(ctx)
  ctx <- gp_feature_apply_raw_selection_to_inputs(ctx)
  low_call_rate_inds_removed <- NULL
  explicit_test_set_override <- gp_test_set_overrides_inferred(ctx$pheno_clean)

  if (isFALSE(((is.null(ctx$geno_data) & is.null(ctx$train_geno_data)) & is.null(ctx$test_geno_data)))) {
    geno_res <- process_geno_data(
      geno_data = ctx$geno_data,
      train_geno_data = ctx$train_geno_data,
      test_geno_data = ctx$test_geno_data,
      test_set = ctx$test_set,
      pheno_clean_list = ctx$pheno_clean,
      train_set = ctx$train_set,
      gen_name = ctx$gen_name,
      kernel_method = ctx$kernel_method,
      gmatrix_method = ctx$gmatrix_method,
      scale = ctx$scale,
      map_data = ctx$map_data,
      maf_threshold = ctx$maf_threshold,
      het_threshold = ctx$het_threshold,
      ind_call_rate_threshold = ctx$ind_call_rate_threshold,
      snp_call_rate_threshold = ctx$snp_call_rate_threshold,
      impute = ctx$impute,
      imputation_method = ctx$imputation_method,
      impute_knn_k = ctx$impute_knn_k,
      ploidy = ctx$ploidy %||% "auto",
      # LD pruning of file input happens only here (model_execute rejects the
      # recoder's PLINK ld_pruning), so it runs at most once.
      ld_prunning_qc = ctx$ld_prunning_qc,
      # A VCF/HapMap matrix was already QC-filtered by its recoder: no second pass.
      qc_filtering = if (!is.null(ctx$geno_data_process)) FALSE else ctx$qc_filtering,
      message = ctx$message,
      heter_groups = ctx$heter_groups,
      test_set_overrides_inferred = explicit_test_set_override
    )

    low_call_rate_inds_removed <- geno_res[["low_call_rate_inds_removed"]]
    if (!is.null(low_call_rate_inds_removed)) {
      low_call_rate_inds_removed <- gp_normalize_test_set(low_call_rate_inds_removed)
      ctx$test_set <- gp_drop_ids_from_set(ctx$test_set, low_call_rate_inds_removed)
      ctx$train_set <- gp_drop_ids_from_set(ctx$train_set, low_call_rate_inds_removed)
      ctx$pheno_clean[["test_set"]] <- gp_drop_ids_from_set(ctx$pheno_clean[["test_set"]], low_call_rate_inds_removed)
      if (!is.null(ctx$pheno_clean[["train_set"]])) {
        ctx$pheno_clean[["train_set"]] <- gp_drop_ids_from_set(ctx$pheno_clean[["train_set"]], low_call_rate_inds_removed)
      }
      if (!is.null(ctx$pheno_clean[["test_set_by_trait"]])) {
        ctx$pheno_clean[["test_set_by_trait"]] <- lapply(
          ctx$pheno_clean[["test_set_by_trait"]],
          gp_drop_ids_from_set,
          removed_ids = low_call_rate_inds_removed
        )
      }
      ctx$pheno_clean[["pheno_clean_data"]] <- ctx$pheno_clean[["pheno_clean_data"]][
        !ctx$pheno_clean[["pheno_clean_data"]][[ctx$gen_name]] %in% low_call_rate_inds_removed,
        ,
        drop = FALSE
      ]
      ctx$pheno_clean[["pheno_clean_data"]] <- gp_order_pheno_for_model_guardrail(
        pheno_data = ctx$pheno_clean[["pheno_clean_data"]],
        gen_name = ctx$gen_name,
        heter_groups = ctx$heter_groups,
        order_for_asreml = gp_asreml_order_required(ctx)
      )
    }
  } else {
    geno_res <- list()
  }

  ctx$omic1_data <- gp_drop_samples_by_id(ctx$omic1_data, low_call_rate_inds_removed)
  ctx$train_omic1_data <- gp_drop_samples_by_id(ctx$train_omic1_data, low_call_rate_inds_removed)
  ctx$test_omic1_data <- gp_drop_samples_by_id(ctx$test_omic1_data, low_call_rate_inds_removed)
  ctx$omic2_data <- gp_drop_samples_by_id(ctx$omic2_data, low_call_rate_inds_removed)
  ctx$train_omic2_data <- gp_drop_samples_by_id(ctx$train_omic2_data, low_call_rate_inds_removed)
  ctx$test_omic2_data <- gp_drop_samples_by_id(ctx$test_omic2_data, low_call_rate_inds_removed)
  ctx$omic3_data <- gp_drop_samples_by_id(ctx$omic3_data, low_call_rate_inds_removed)
  ctx$train_omic3_data <- gp_drop_samples_by_id(ctx$train_omic3_data, low_call_rate_inds_removed)
  ctx$test_omic3_data <- gp_drop_samples_by_id(ctx$test_omic3_data, low_call_rate_inds_removed)

  omic1_res <- process_omic_data(
    omic_data = ctx$omic1_data,
    train_omic_data = ctx$train_omic1_data,
    test_omic_data = ctx$test_omic1_data,
    kernel_method = ctx$kernel_method,
    pheno_clean_list = ctx$pheno_clean,
    gen_name = ctx$gen_name,
    test_set = ctx$test_set,
    train_set = ctx$train_set,
    message = ctx$message,
    heter_groups = ctx$heter_groups,
    impute_omic = ctx$impute_omic,
    imputation_method = ctx$imputation_method,
    impute_knn_k = ctx$impute_knn_k,
    na_threshold = ctx$na_threshold,
    test_set_overrides_inferred = explicit_test_set_override
  )

  omic2_res <- process_omic_data(
    omic_data = ctx$omic2_data,
    train_omic_data = ctx$train_omic2_data,
    test_omic_data = ctx$test_omic2_data,
    kernel_method = ctx$kernel_method,
    pheno_clean_list = ctx$pheno_clean,
    gen_name = ctx$gen_name,
    test_set = ctx$test_set,
    train_set = ctx$train_set,
    message = ctx$message,
    heter_groups = ctx$heter_groups,
    impute_omic = ctx$impute_omic,
    imputation_method = ctx$imputation_method,
    impute_knn_k = ctx$impute_knn_k,
    na_threshold = ctx$na_threshold,
    test_set_overrides_inferred = explicit_test_set_override
  )

  omic3_res <- process_omic_data(
    omic_data = ctx$omic3_data,
    train_omic_data = ctx$train_omic3_data,
    test_omic_data = ctx$test_omic3_data,
    kernel_method = ctx$kernel_method,
    pheno_clean_list = ctx$pheno_clean,
    gen_name = ctx$gen_name,
    test_set = ctx$test_set,
    train_set = ctx$train_set,
    message = ctx$message,
    heter_groups = ctx$heter_groups,
    impute_omic = ctx$impute_omic,
    imputation_method = ctx$imputation_method,
    impute_knn_k = ctx$impute_knn_k,
    na_threshold = ctx$na_threshold,
    test_set_overrides_inferred = explicit_test_set_override
  )

  resolved_test_set <- gp_merge_test_sets(
    ctx$test_set,
    geno_res[["test_set"]],
    omic1_res[["test_set"]],
    omic2_res[["test_set"]],
    omic3_res[["test_set"]]
  )
  if (!is.null(resolved_test_set)) {
    if (isTRUE(explicit_test_set_override)) {
      ctx$pheno_clean[["test_set"]] <- gp_normalize_test_set(ctx$pheno_clean[["test_set"]])
    } else {
      ctx$pheno_clean[["test_set"]] <- gp_merge_test_sets(ctx$pheno_clean[["test_set"]], resolved_test_set)
      if (is.null(ctx$pheno_clean[["test_set_source"]])) {
        ctx$pheno_clean[["test_set_source"]] <- "inferred_geno_omic"
      } else if (!grepl("geno_omic", ctx$pheno_clean[["test_set_source"]], fixed = TRUE)) {
        ctx$pheno_clean[["test_set_source"]] <- paste(ctx$pheno_clean[["test_set_source"]], "geno_omic", sep = "+")
      }
    }
  }

  geno_omic_model_ready_list <- list()
  gmatrix_kernel_model_ready_list <- list()

  if (length(geno_res) != 0 && all(c("gmatrix", "geno_model_ready") %in% names(geno_res))) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "gmatrix_model_ready",
      geno_res[["gmatrix"]]
    )
    geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
  } else if (length(geno_res) != 0 && "geno_model_ready" %in% names(geno_res)) {
    geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
  }

  if (!is.null(ctx$gmatrix)) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "gmatrix_model_ready",
      gp_drop_test_kernel_bank(ctx$gmatrix, low_call_rate_inds_removed, "gmatrix")
    )
  }
  if (!is.null(ctx$gkernel)) {
    gkernel_ready <- gp_drop_test_kernel_bank(
      ctx$gkernel, low_call_rate_inds_removed, "gkernel"
    )
    if (gp_is_kernel_matrix(gkernel_ready) &&
        is.null(attr(gkernel_ready, "predictpror_source_kernel", exact = TRUE))) {
      attr(gkernel_ready, "predictpror_source_kernel") <- "gkernel"
    }
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      if (is.null(gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]])) {
        "gmatrix_model_ready"
      } else {
        "gkernel_model_ready"
      },
      gkernel_ready
    )
  }
  if (!is.null(ctx$kernel_list)) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "kernel_list_model_ready",
      gp_drop_test_kernel_bank(ctx$kernel_list, low_call_rate_inds_removed, "kernel_list")
    )
  }

  omic1_inputs <- gp_add_model_input(
    model_inputs = geno_omic_model_ready_list,
    kernel_inputs = gmatrix_kernel_model_ready_list,
    result = omic1_res,
    model_key = "omic1_model_ready",
    kernel_key = "omic1_kernel_model_ready"
  )
  geno_omic_model_ready_list <- omic1_inputs$model_inputs
  gmatrix_kernel_model_ready_list <- omic1_inputs$kernel_inputs

  omic2_inputs <- gp_add_model_input(
    model_inputs = geno_omic_model_ready_list,
    kernel_inputs = gmatrix_kernel_model_ready_list,
    result = omic2_res,
    model_key = "omic2_model_ready",
    kernel_key = "omic2_kernel_model_ready"
  )
  geno_omic_model_ready_list <- omic2_inputs$model_inputs
  gmatrix_kernel_model_ready_list <- omic2_inputs$kernel_inputs

  omic3_inputs <- gp_add_model_input(
    model_inputs = geno_omic_model_ready_list,
    kernel_inputs = gmatrix_kernel_model_ready_list,
    result = omic3_res,
    model_key = "omic3_model_ready",
    kernel_key = "omic3_kernel_model_ready"
  )
  geno_omic_model_ready_list <- omic3_inputs$model_inputs
  gmatrix_kernel_model_ready_list <- omic3_inputs$kernel_inputs

  geno_omic_model_ready_list <- gp_normalize_feature_bank_guardrail(
    geno_omic_model_ready_list,
    "model_ready_feature",
    id_col = ctx$gen_name
  )

  if (!is.null(ctx$omic1_kernel) && is.null(gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]])) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "omic1_kernel_model_ready",
      gp_drop_test_kernel_bank(ctx$omic1_kernel, low_call_rate_inds_removed, "omic1_kernel")
    )
  }
  if (!is.null(ctx$omic2_kernel) && is.null(gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]])) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "omic2_kernel_model_ready",
      gp_drop_test_kernel_bank(ctx$omic2_kernel, low_call_rate_inds_removed, "omic2_kernel")
    )
  }
  if (!is.null(ctx$omic3_kernel) && is.null(gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]])) {
    gmatrix_kernel_model_ready_list <- gp_add_kernel_inputs(
      gmatrix_kernel_model_ready_list,
      "omic3_kernel_model_ready",
      gp_drop_test_kernel_bank(ctx$omic3_kernel, low_call_rate_inds_removed, "omic3_kernel")
    )
  }

  gmatrix_kernel_model_ready_list <- gp_normalize_kernel_guardrail(
    gmatrix_kernel_model_ready_list,
    "model_ready_kernel"
  )

  kernel_vars <- names(gmatrix_kernel_model_ready_list)

  for (kernel_var in kernel_vars) {
    if (gp_is_kernel_matrix(gmatrix_kernel_model_ready_list[[kernel_var]])) {
      source_kernel_name <- gp_kernel_source_name(
        gmatrix_kernel_model_ready_list[[kernel_var]],
        gp_kernel_key_source_name(kernel_var)
      )
      gmatrix_kernel_model_ready_list[[kernel_var]] <- grm_kernel_precheck(
        grm_kernel_data = gmatrix_kernel_model_ready_list[[kernel_var]],
        pedigree_matrix = ctx$pedigree_matrix,
        bending = ctx$bending,
        bend_value = ctx$bend_value,
        blending = ctx$blending,
        blending_value = ctx$blending_value,
        high_diag_cut_off = ctx$high_diag_cut_off,
        low_diag_cut_off = ctx$low_diag_cut_off,
        duplicate_cut_off = ctx$duplicate_cut_off,
        rcn_cutoff = ctx$rcn_cutoff,
        optimize_diagonal = ctx$optimize_diagonal,
        optimize_duplicate = ctx$optimize_duplicate,
        show_message = ctx$message,
        kernel_check_level = ctx$kernel_check_level %||% "auto",
        kernel_large_n_threshold = ctx$kernel_large_n_threshold %||% 5000L,
        duplicate_scan = ctx$duplicate_scan %||% "auto",
        duplicate_sample_size = ctx$duplicate_sample_size %||% 2000L,
        duplicate_block_size = ctx$duplicate_block_size %||% 1024L,
        duplicate_max_pairs = ctx$duplicate_max_pairs %||% 10000L,
        kernel_fix_method = ctx$kernel_fix_method %||% "auto",
        kernel_repair_priority = ctx$kernel_repair_priority %||% "speed",
        kernel_rcn_check = ctx$kernel_rcn_check %||% "auto",
        kernel_nearpd_size_limit = ctx$kernel_nearpd_size_limit %||% 2500L,
        kernel_cpp_repair_size_limit = ctx$kernel_cpp_repair_size_limit %||% 3000L,
        kernel_cpp_keep_diag = ctx$kernel_cpp_keep_diag %||% TRUE,
        kernel_pd_check = ctx$kernel_pd_check %||% "auto",
        kernel_pd_sample_size = ctx$kernel_pd_sample_size %||% 500L,
        kernel_sanitize = ctx$kernel_sanitize %||% "auto",
        kernel_sanitize_value = ctx$kernel_sanitize_value %||% NULL
      )
      attr(
        gmatrix_kernel_model_ready_list[[kernel_var]],
        "predictpror_source_kernel"
      ) <- source_kernel_name
    }
  }

  test_set <- resolved_test_set
  if (length(gmatrix_kernel_model_ready_list) != 0) {
    for (checked_kernel_var_name in names(gmatrix_kernel_model_ready_list)) {
      if (gp_is_kernel_matrix(gmatrix_kernel_model_ready_list[[checked_kernel_var_name]])) {
        source_kernel_name <- gp_kernel_source_name(
          gmatrix_kernel_model_ready_list[[checked_kernel_var_name]],
          gp_kernel_key_source_name(checked_kernel_var_name)
        )
        match_result <- pheno_geno_match(
          object_pheno = ctx$pheno_clean[["pheno_clean_data"]],
          object_geno = gmatrix_kernel_model_ready_list[[checked_kernel_var_name]],
          gen_name = ctx$gen_name,
          test_set = test_set,
          train_set = ctx$train_set,
          test_set_overrides_inferred = explicit_test_set_override,
          message = ctx$message
        )

        gmatrix_kernel_model_ready_list[[checked_kernel_var_name]] <- match_result[[1]]
        attr(
          gmatrix_kernel_model_ready_list[[checked_kernel_var_name]],
          "predictpror_source_kernel"
        ) <- source_kernel_name
        if (length(match_result) > 1) {
          test_set <- gp_normalize_test_set(match_result[[2]])
        }
      }
    }
  }
  if (!is.null(test_set) && !isTRUE(explicit_test_set_override)) {
    ctx$pheno_clean[["test_set"]] <- gp_merge_test_sets(ctx$pheno_clean[["test_set"]], test_set)
    if (is.null(ctx$pheno_clean[["test_set_source"]])) {
      ctx$pheno_clean[["test_set_source"]] <- "inferred_geno_omic"
    } else if (!grepl("geno_omic", ctx$pheno_clean[["test_set_source"]], fixed = TRUE)) {
      ctx$pheno_clean[["test_set_source"]] <- paste(ctx$pheno_clean[["test_set_source"]], "geno_omic", sep = "+")
    }
  } else if (isTRUE(explicit_test_set_override)) {
    test_set <- gp_normalize_test_set(test_set)
  }

  # Composite dispatchers execute one protected route at a time, but their
  # first child must prepare the union of inputs needed by every requested
  # model.  Otherwise a kernel-first child omits ML/DL data, while an ML-first
  # child can omit the relationship matrix later GP/Bayesian/ASReml children
  # require.
  requested_preprocess_models <- gp_shared_preprocess_requested_models(ctx)
  model_check <- any(requested_preprocess_models %in% ctx$AI_valid_models)
  model_check <- isTRUE(model_check) || isTRUE(ctx$feature_scoring) ||
    !is.null(ctx$feature_score_metadata) ||
    !is.null(ctx$feature_k_grid) ||
    !is.null(ctx$feature_k) ||
    !is.null(ctx$feature_selected)

  use_met_ml_dl <- isTRUE(ctx$met_ml_dl) &&
    gp_is_multi_environment_pheno(
      ctx$pheno_clean[["pheno_clean_data"]],
      ctx$gen_name,
      heter_groups = ctx$heter_groups,
      response = ctx$response,
      response_family = ctx$response_family
    ) &&
    any(requested_preprocess_models %in% gp_met_supported_models())

  # Classical ML and DL do not estimate covariance components from a kernel,
  # but any user-supplied PSD kernel has an exact feature-map representation.
  # For non-MET workflows, expose every explicitly supplied kernel as a named
  # eigenfeature block and concatenate it with any raw genomic/omics blocks.
  # Do not add the automatically constructed GRM when the user supplied only
  # raw markers: that would duplicate the same information and silently change
  # longstanding marker-level ML/DL fits.  MET has its own equivalent builder.
  explicit_kernel_requested <- any(vapply(
    list(
      ctx$gmatrix, ctx$gkernel, ctx$omic1_kernel, ctx$omic2_kernel,
      ctx$omic3_kernel, ctx$kernel_list
    ),
    Negate(is.null),
    logical(1L)
  ))
  kernel_feature_summary <- data.frame()
  if (isTRUE(model_check) && !isTRUE(use_met_ml_dl) &&
      isTRUE(explicit_kernel_requested) &&
      length(gmatrix_kernel_model_ready_list)) {
    kernel_features <- gp_prepare_met_kernel_features(
      kernel_list = gmatrix_kernel_model_ready_list,
      omics_kernel_label = ctx$omics_kernel_label,
      var_explained = ctx$met_kernel_var_explained %||% 0.95,
      min_ev = ctx$met_kernel_min_ev %||% 1e-8,
      max_pcs = ctx$met_kernel_max_pcs
    )
    if (!is.null(kernel_features$feature_table) &&
        ncol(kernel_features$feature_table) > 0L) {
      geno_omic_model_ready_list[["kernel_features_model_ready"]] <-
        kernel_features$feature_table
      kernel_feature_summary <- kernel_features$feature_summary
    }
  }

  ml_dat_res <- if (isTRUE(model_check) && !use_met_ml_dl) {
    AI_process_ml_data_if_valid(
      model_check = model_check,
      geno_omic_model_ready_list,
      ctx$pheno_clean,
      ctx$response,
      ctx$gen_name
    )
  } else {
    list()
  }
  if (length(ml_dat_res) && nrow(kernel_feature_summary)) {
    ml_dat_res[["kernel_feature_summary"]] <- kernel_feature_summary
    ml_dat_res[["kernel_feature_strategy"]] <-
      "concatenated_named_psd_eigenfeature_blocks"
  }

  list(
    low_call_rate_inds_removed = low_call_rate_inds_removed,
    pheno_clean = ctx$pheno_clean,
    geno_res = geno_res,
    omic1_res = omic1_res,
    omic2_res = omic2_res,
    omic3_res = omic3_res,
    geno_omic_model_ready_list = geno_omic_model_ready_list,
    gmatrix_kernel_model_ready_list = gmatrix_kernel_model_ready_list,
    test_set = test_set,
    ml_dat_res = ml_dat_res
  )
}

gp_cv_generate_plots_enabled <- function(ctx) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  requested <- ctx$cv_generate_plots %||% NULL
  if (is.character(requested) && length(requested) && identical(tolower(requested[[1L]]), "auto")) {
    requested <- NULL
  }
  if (!is.null(requested)) {
    return(isTRUE(requested))
  }

  cv_models <- as.character(ctx$GS_model_cv %||% character())
  gp_models <- as.character(ctx$gp_valid_models %||% character())
  gp_only <- length(cv_models) > 0L && length(gp_models) > 0L && all(cv_models %in% gp_models)

  TRUE
}

gp_cv_empty_plot_lists <- function(eval_metrics) {
  eval_metrics <- as.character(eval_metrics %||% character())
  list(
    plot_reps_list = stats::setNames(
      lapply(eval_metrics, function(metric) {
        list(metric = metric, ggplot_boxplot_reps = NULL, plotly_boxplot_reps = NULL)
      }),
      eval_metrics
    ),
    plot_mean_list = stats::setNames(
      lapply(eval_metrics, function(metric) {
        list(metric = metric, ggplot_lineplot_mean = NULL, plotly_lineplot_mean = NULL)
      }),
      eval_metrics
    )
  )
}

gp_cv_metric_mean_by <- function(data, by_cols, metric) {
  if (!metric %in% names(data)) {
    return(NULL)
  }
  stats::aggregate(
    stats::as.formula(paste(metric, "~", paste(by_cols, collapse = " + "))),
    data = data,
    FUN = mean,
    na.rm = TRUE
  )
}

gp_cv_merge_metric_frames <- function(frames, by_cols) {
  frames <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, frames)
  if (!length(frames)) {
    return(data.frame())
  }
  Reduce(function(x, y) merge(x, y, by = by_cols, all = TRUE, sort = FALSE), frames)
}

gp_cv_rank_metric_frames <- function(metric_data,
                                     eval_metrics,
                                     metric_for_ranking = "accuracy",
                                     ranking_tie_breakers = NULL) {
  eval_metrics <- as.character(eval_metrics %||% character())
  out <- lapply(eval_metrics, function(metric) {
    if (!is.data.frame(metric_data) || !metric %in% names(metric_data)) {
      return(data.frame())
    }
    tie_breakers <- if (identical(metric, metric_for_ranking)) {
      ranking_tie_breakers
    } else {
      NULL
    }
    tryCatch(
      rank_models(metric_data, metric, tie_breakers = tie_breakers),
      error = function(e) data.frame()
    )
  })
  stats::setNames(out, eval_metrics)
}

gp_classification_cv_selection_summary <- function(classification_metrics,
                                                   best_models_list,
                                                   response_family,
                                                   metric_for_ranking,
                                                   ranking_tie_breakers = NULL,
                                                   positive_class = NULL,
                                                   cv_results_data = NULL,
                                                   heter_groups = NULL,
                                                   is_multi = FALSE) {
  fam <- gp_resolve_response_family(response_family)
  if (identical(fam, "gaussian")) {
    return(NULL)
  }
  primary <- gp_resolve_metric_for_ranking(metric_for_ranking, fam)
  tie_breakers <- gp_resolve_ranking_tie_breakers(
    ranking_tie_breakers = ranking_tie_breakers,
    response_family = fam,
    metric_for_ranking = primary
  )
  candidates <- classification_metrics$aggregated %||% data.frame()
  selected <- best_models_list[[primary]] %||% data.frame()
  models <- if (is.data.frame(candidates) && "model" %in% names(candidates)) {
    sort(unique(as.character(candidates$model)))
  } else {
    character()
  }
  cv_scenarios <- if (is.list(cv_results_data) && length(cv_results_data)) {
    unique(stats::na.omit(vapply(cv_results_data, function(x) {
      method <- as.character(x$cv_info$method %||% NA_character_)
      if (length(method)) method[[1L]] else NA_character_
    }, character(1L))))
  } else {
    character()
  }
  list(
    policy = data.frame(
      response_family = fam,
      primary_metric = primary,
      primary_direction = gp_eval_metric_direction(primary),
      tie_breakers = paste(tie_breakers, collapse = ","),
      positive_class = as.character(positive_class %||% NA_character_),
      cv_compares_models = length(models) > 1L,
      true_prediction_performed = FALSE,
      multi_environment = isTRUE(is_multi),
      environment_column = if (isTRUE(is_multi)) {
        as.character(heter_groups %||% NA_character_)[[1L]]
      } else {
        NA_character_
      },
      cv_scenario = if (length(cv_scenarios)) {
        paste(cv_scenarios, collapse = ",")
      } else {
        NA_character_
      },
      selection_scope = if (isTRUE(is_multi)) {
        "pooled_out_of_fold_across_environments"
      } else {
        "pooled_out_of_fold"
      },
      environment_aggregation = if (isTRUE(is_multi)) {
        "ALL_rows_mean_across_replications"
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    ),
    models_compared = data.frame(model = models, stringsAsFactors = FALSE),
    candidate_metrics = candidates,
    selected_models = selected
  )
}

gp_cv_process_results_fast <- function(cv_results_data,
                                       eval_metrics,
                                       heter_groups = NULL,
                                       is_multi = FALSE,
                                       metric_for_ranking = "accuracy",
                                       ranking_tie_breakers = NULL,
                                       positive_class = NULL) {
  eval_metrics <- as.character(eval_metrics %||% character())
  combined_df <- do.call(rbind, lapply(cv_results_data, function(x) {
    cbind(
      data.frame(trait = x$trait, rep = x$rep, model = x$model, stringsAsFactors = FALSE),
      x$eval_metrics_reps
    )
  }))
  combined_df <- as.data.frame(combined_df, stringsAsFactors = FALSE)

  non_numeric <- c("trait", "model", heter_groups, "feature_scoring_model", "feature_scoring_cv")
  for (nm in setdiff(names(combined_df), non_numeric)) {
    combined_df[[nm]] <- suppressWarnings(as.numeric(combined_df[[nm]]))
  }
  if ("rep" %in% names(combined_df)) {
    combined_df$Rep <- as.factor(combined_df$rep)
  }

  empty_plots <- gp_cv_empty_plot_lists(eval_metrics)
  if (isTRUE(is_multi)) {
    if (is.null(heter_groups) || !heter_groups %in% names(combined_df)) {
      stop("heter_groups is required for fast multi-environment CV processing.", call. = FALSE)
    }
    feature_cols <- intersect(c("feature_k", "feature_scoring_model", "feature_scoring_cv"), names(combined_df))
    environment_value <- as.character(combined_df[[heter_groups]])
    pooled_df <- combined_df[environment_value == "ALL", , drop = FALSE]
    environment_df <- combined_df[environment_value != "ALL", , drop = FALSE]
    if (!nrow(environment_df)) {
      environment_df <- combined_df
    }
    pooled_group_cols <- c("trait", "model", feature_cols)
    aggregated_data_list <- lapply(eval_metrics, function(metric) {
      grouped_cols <- c("trait", "model", feature_cols, heter_groups)
      aggregated_data <- gp_cv_metric_mean_by(environment_df, grouped_cols, metric)
      aggregated_data_env_mean <- gp_cv_metric_mean_by(
        aggregated_data, pooled_group_cols, metric
      )
      list(
        metric = metric,
        aggregated_data = aggregated_data %||% data.frame(),
        aggregated_data_env_mean = aggregated_data_env_mean %||% data.frame()
      )
    })
    names(aggregated_data_list) <- eval_metrics
    ranking_source <- if (nrow(pooled_df)) {
      gp_cv_merge_metric_frames(
        lapply(eval_metrics, function(metric) {
          gp_cv_metric_mean_by(pooled_df, pooled_group_cols, metric)
        }),
        pooled_group_cols
      )
    } else {
      gp_cv_merge_metric_frames(
        lapply(aggregated_data_list, `[[`, "aggregated_data_env_mean"),
        pooled_group_cols
      )
    }
    aggregated_data_list$combined_df <- combined_df
    aggregated_data_list$aggregated_across_reps <- ranking_source
  } else {
    feature_cols <- intersect(c("feature_k", "feature_scoring_model", "feature_scoring_cv"), names(combined_df))
    aggregated_across_reps <- gp_cv_merge_metric_frames(
      lapply(eval_metrics, function(metric) gp_cv_metric_mean_by(combined_df, c("trait", "model", feature_cols), metric)),
      c("trait", "model", feature_cols)
    )
    aggregated_data_list <- list(
      aggregated_across_reps = aggregated_across_reps,
      combined_df = combined_df
    )
    ranking_source <- aggregated_across_reps
  }

  ranking_for_selection <- ranking_source
  if (is.data.frame(ranking_for_selection) &&
      "feature_scoring_cv" %in% names(ranking_for_selection) &&
      all(c("fixed", "fold_internal") %in% unique(as.character(ranking_for_selection[["feature_scoring_cv"]])))) {
    ranking_for_selection <- ranking_for_selection[
      as.character(ranking_for_selection[["feature_scoring_cv"]]) == "fold_internal",
      ,
      drop = FALSE
    ]
  }

  best_models_list <- gp_cv_rank_metric_frames(
    ranking_for_selection,
    eval_metrics,
    metric_for_ranking = metric_for_ranking,
    ranking_tie_breakers = ranking_tie_breakers
  )
  response_family <- unique(vapply(cv_results_data, function(x) {
    as.character(x$response_family %||% "gaussian")
  }, character(1L)))
  response_family <- response_family[[1L]] %||% "gaussian"
  classification_metrics <- aggregate_cv_classification_metrics(
    cv_results_data, aggregated_data_list, eval_metrics,
    heter_groups = if (isTRUE(is_multi)) heter_groups else NULL
  )
  classification_reliability <- aggregate_cv_classification_reliability(
    cv_results_data,
    heter_groups = if (isTRUE(is_multi)) heter_groups else NULL
  )
  if (!is.null(classification_metrics)) {
    classification_metrics$confusion_matrices <-
      aggregate_cv_classification_confusion_matrices(
        classification_reliability,
        heter_groups = if (isTRUE(is_multi)) heter_groups else NULL
      )
  }

  list(
    aggregated_data_list = aggregated_data_list,
    plot_reps_list = empty_plots$plot_reps_list,
    plot_mean_list = empty_plots$plot_mean_list,
    best_models_list = best_models_list,
    classification_probability_summaries = aggregate_cv_probability_summaries(
      cv_results_data,
      heter_groups = if (isTRUE(is_multi)) heter_groups else NULL
    ),
    classification_metrics = classification_metrics,
    classification_reliability = classification_reliability,
    classification_model_selection = gp_classification_cv_selection_summary(
      classification_metrics = classification_metrics,
      best_models_list = best_models_list,
      response_family = response_family,
      metric_for_ranking = metric_for_ranking,
      ranking_tie_breakers = ranking_tie_breakers,
      positive_class = positive_class,
      cv_results_data = cv_results_data,
      heter_groups = heter_groups,
      is_multi = is_multi
    ),
    gaussian_uncertainty_summaries = list(),
    gaussian_risk_summaries = list(),
    gaussian_risk_plots = list()
  )
}

gp_cv_frontdoor_fast_enabled <- function() {
  tolower(Sys.getenv("PREDICTPRO_GP_CV_FRONTDOOR_FAST", unset = "true")) %in%
    c("1", "true", "yes", "y", "auto")
}

gp_cv_frontdoor_fast_info_keys <- function() {
  c(
    "package_cv_batch", "bridge_command", "cv_fast_path",
    "cv_fast_path_skip_reason", "cv_fast_path_approximate",
    "cv_fast_path_exact_refit", "cv_fast_path_solver",
    "cv_fast_path_shared_hyperparameters",
    "cv_fast_path_likelihood_noise_variance",
    "cv_fast_path_likelihood_noise_std", "cv_fast_path_self_checked",
    "cv_fast_path_self_check_reason", "cv_preprocessed_context",
    "cv_preprocessed_context_rows",
    "cv_preprocessed_context_tensor_cache_entries", "gp_exact_fast_cv"
  )
}

gp_cv_frontdoor_fast_metric_frame <- function(ypred_cv,
                                              rep_i,
                                              eval_metrics,
                                              response_family = "gaussian",
                                              heter_groups = NULL,
                                              is_multi = FALSE) {
  eval_metrics <- as.character(eval_metrics %||% character())
  if (!isTRUE(is_multi)) {
    oof <- which(ypred_cv$cv_role == "test" & !is.na(ypred_cv$yhat))
    out <- data.frame(Rep = as.integer(rep_i), check.names = FALSE)
    for (metric in eval_metrics) {
      out[[metric]] <- tryCatch(
        safe_metric_value(
          y_true = ypred_cv$y[oof],
          y_pred = ypred_cv$yhat[oof],
          metric = metric,
          response_family = response_family
        ),
        error = function(e) NA_real_
      )
    }
    return(out)
  }

  if (is.null(heter_groups) || !heter_groups %in% names(ypred_cv)) {
    stop("heter_groups is required for GP front-door multi-environment CV.", call. = FALSE)
  }
  test_envs <- sort(unique(ypred_cv[[heter_groups]][ypred_cv$cv_role == "test"]))
  by_env <- do.call(
    rbind,
    lapply(test_envs, function(env_i) {
      idx <- which(ypred_cv$cv_role == "test" & ypred_cv[[heter_groups]] == env_i)
      row <- list(Rep = as.integer(rep_i))
      row[[heter_groups]] <- env_i
      for (metric in eval_metrics) {
        ok <- gp_metric_valid_rows(ypred_cv$y[idx], ypred_cv$yhat[idx], response_family = response_family)
        row[[metric]] <- if (any(ok)) {
          safe_metric_value(
            y_true = ypred_cv$y[idx][ok],
            y_pred = ypred_cv$yhat[idx][ok],
            metric = metric,
            response_family = response_family
          )
        } else {
          NA_real_
        }
      }
      as.data.frame(row, check.names = FALSE)
    })
  )

  oof <- which(ypred_cv$cv_role == "test")
  ok <- gp_metric_valid_rows(ypred_cv$y[oof], ypred_cv$yhat[oof], response_family = response_family)
  overall <- data.frame(Rep = as.integer(rep_i), check.names = FALSE)
  overall[[heter_groups]] <- "ALL"
  for (metric in eval_metrics) {
    overall[[metric]] <- if (any(ok)) {
      safe_metric_value(
        y_true = ypred_cv$y[oof][ok],
        y_pred = ypred_cv$yhat[oof][ok],
        metric = metric,
        response_family = response_family
      )
    } else {
      NA_real_
    }
  }
  rbind(by_env, overall)
}

gp_cv_frontdoor_fast_splits <- function(pheno_data,
                                        trait,
                                        gen_name,
                                        heter_groups,
                                        cv_token,
                                        cross_validation_meth,
                                        sampling_method,
                                        nfolds,
                                        test_size = NULL,
                                        random_state,
                                        rep_i,
                                        response_family = "gaussian",
                                        verbose = FALSE) {
  trait_values <- pheno_data[[trait]]
  valid_idx <- if (identical(response_family, "gaussian")) {
    which(is.finite(as.double(trait_values)))
  } else {
    which(!is.na(trait_values))
  }
  if (!length(valid_idx)) {
    stop("Trait has no finite observations for GP front-door CV.", call. = FALSE)
  }
  pheno_data_split <- pheno_data[valid_idx, , drop = FALSE]
  base_seed <- if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1L) {
    as.integer(random_state)
  } else {
    123L
  }
  new_seed <- (base_seed + as.integer(rep_i) * 10000L) %% .Machine$integer.max

  holds_tokens <- c("hold_out", "stratified_hold_out", "repeated_hold_out", "repeated_stratified_hold_out")
  kfold_tokens <- c("k_folds", "stratified_k_folds", "repeated_k_folds", "repeated_stratified_k_folds")
  multi_tokens <- c("cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2")
  loo_tokens <- c("leave_one_out")

  is_holdout <- cv_token %in% holds_tokens
  is_kfold <- cv_token %in% kfold_tokens
  is_multi <- cv_token %in% multi_tokens
  is_loo <- cv_token %in% loo_tokens
  if (!(is_holdout || is_kfold || is_multi || is_loo)) {
    stop("Unsupported GP front-door CV method.", call. = FALSE)
  }

  effective_nfolds <- as.integer(nfolds %||% 5L)
  if (is_holdout) {
    test_set_val <- hold_out_stratified_and_un(
      pheno_data = pheno_data_split,
      gen_name = gen_name,
      response = trait,
      test_size = test_size %||% 0.2,
      random_state = new_seed,
      replication = as.integer(rep_i),
      sampling_method = sampling_method
    )
    group <- rep(NA_integer_, nrow(pheno_data))
    tst <- valid_idx[test_set_val[[as.integer(rep_i)]]]
    group[tst] <- 1L
    return(list(group = group, folds = list(tst), nfolds = 1L, is_multi = FALSE))
  }

  if (is_loo) {
    effective_nfolds <- nrow(pheno_data_split)
    group <- rep(NA_integer_, nrow(pheno_data))
    group[valid_idx] <- seq_len(effective_nfolds)
    folds <- lapply(seq_len(effective_nfolds), function(j) which(group == j))
    return(list(group = group, folds = folds, nfolds = effective_nfolds, is_multi = FALSE))
  }

  if (is_kfold) {
    test_set_val <- kfolds_stratified_un(
      pheno_data = pheno_data_split,
      gen_name = gen_name,
      response = trait,
      nfolds = effective_nfolds,
      random_state = new_seed,
      replication = as.integer(rep_i),
      sampling_method = sampling_method
    )
    group <- rep(NA_integer_, nrow(pheno_data))
    group[valid_idx] <- test_set_val[[as.integer(rep_i)]]
    folds <- lapply(seq_len(effective_nfolds), function(j) which(group == j))
    return(list(group = group, folds = folds, nfolds = effective_nfolds, is_multi = FALSE))
  }

  if (is.null(heter_groups)) {
    stop("heter_groups is required for GP front-door CV0/CV1/CV2.", call. = FALSE)
  }
  CVi <- switch(
    cv_token,
    "cv0" = 0L,
    "repeated_cv0" = 0L,
    "cv1" = 1L,
    "repeated_cv1" = 1L,
    "cv2" = 2L,
    "repeated_cv2" = 2L,
    stop("Unsupported GP front-door multi-environment CV method.", call. = FALSE)
  )
  test_set_val <- CV0_CV1_CV2_for_multi_environment(
    pheno_data = pheno_data_split,
    gen_name = gen_name,
    heter_groups = heter_groups,
    CV = CVi,
    nfolds = effective_nfolds,
    random_state = new_seed,
    replication = as.integer(rep_i),
    message = isTRUE(verbose)
  )
  group <- rep(NA_integer_, nrow(pheno_data))
  group[valid_idx] <- test_set_val[[as.integer(rep_i)]]
  fold_ids <- sort(unique(group))
  fold_ids <- fold_ids[!is.na(fold_ids) & fold_ids > 0L]
  folds <- lapply(fold_ids, function(j) which(group == j))
  list(group = group, folds = folds, nfolds = NA_integer_, is_multi = TRUE)
}

gp_try_cv_frontdoor_fast_lane <- function(ctx) {
  if (!isTRUE(gp_cv_frontdoor_fast_enabled())) {
    return(NULL)
  }
  if (!isTRUE(ctx$cross_validation) || !isTRUE(ctx$cv_evaluation_only)) {
    return(NULL)
  }
  if (!identical(ctx$response_family %||% "gaussian", "gaussian")) {
    return(NULL)
  }
  if (isTRUE(ctx$hybrid_dl) || isTRUE(ctx$hybrid_ml) ||
      isTRUE(ctx$multi_trait_ml) || isTRUE(ctx$multi_trait_dl) ||
      isTRUE(ctx$multi_trait_asreml) || isTRUE(ctx$hybrid_asreml) ||
      isTRUE(ctx$hybrid_bayes) || isTRUE(ctx$hybrid_gp)) {
    return(NULL)
  }

  response <- as.character(ctx$response %||% character())
  model <- as.character(ctx$GS_model_cv %||% character())
  model_label <- as.character(ctx$GS_model_cv_display %||% gp_display_supported_model_names(model))
  if (length(model_label) != length(model)) {
    model_label <- gp_display_supported_model_names(model)
  }
  gp_valid_models <- as.character(ctx$gp_valid_models %||% gp_lowrank_supported_models())
  if (length(response) != 1L || length(model) != 1L || !model %in% gp_valid_models) {
    return(NULL)
  }
  if (!exists("gp_backend_cv_predict_batch", mode = "function") ||
      !isTRUE(gp_backend_cv_batch_enabled())) {
    return(NULL)
  }
  if (exists("gp_bridge_direct_execution_enabled", mode = "function") &&
      !isTRUE(gp_bridge_direct_execution_enabled())) {
    return(NULL)
  }
  if ((!is.null(ctx$ml_dat_res) && length(ctx$ml_dat_res) > 0L) || isTRUE(ctx$met_ml_dl)) {
    return(NULL)
  }

  pheno_data <- ctx$pheno_clean[["pheno_clean_data"]]
  if (is.null(pheno_data) || !is.data.frame(pheno_data)) {
    return(NULL)
  }
  test_set <- ctx$test_set %||% NULL
  if (!is.null(test_set)) {
    test_set <- unique(as.character(test_set))
    pheno_data <- pheno_data[!as.character(pheno_data[[ctx$gen_name]]) %in% test_set, , drop = FALSE]
  }
  if (!nrow(pheno_data)) {
    return(NULL)
  }

  gmatrix_model_ready <- ctx$gmatrix_model_ready %||%
    ctx$gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] %||%
    ctx$gmatrix %||% NULL
  if (is.null(gmatrix_model_ready)) {
    return(NULL)
  }

  cv_token <- normalize_cv_token(ctx$cross_validation_meth)
  supported_tokens <- c(
    "hold_out", "stratified_hold_out", "repeated_hold_out", "repeated_stratified_hold_out",
    "k_folds", "stratified_k_folds", "repeated_k_folds", "repeated_stratified_k_folds",
    "cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2",
    "leave_one_out"
  )
  if (!cv_token %in% supported_tokens) {
    return(NULL)
  }
  is_multi_cv <- cv_token %in% c("cv0", "cv1", "cv2", "repeated_cv0", "repeated_cv1", "repeated_cv2")
  heter_groups <- ctx$heter_groups %||% NULL
  if (isTRUE(is_multi_cv) && (is.null(heter_groups) || !heter_groups %in% names(pheno_data))) {
    return(NULL)
  }
  if (!is.null(ctx$env_similarity) && !is.null(ctx$env_covariates)) {
    return(NULL)
  }

  replication <- as.integer(ctx$replication %||% 1L)
  if (!is.finite(replication) || replication < 1L) {
    return(NULL)
  }
  sampling_method <- gp_resolve_cv_sampling_method(ctx$cross_validation_meth, ctx$sampling_method)
  cv_profile <- gp_runtime_profile_new("gp_cv_frontdoor_fast_lane")
  trait <- response[[1L]]
  y <- as.double(pheno_data[[trait]])
  params <- list(
    pheno_data = pheno_data,
    response = trait,
    gen_name = ctx$gen_name,
    heter_groups = heter_groups,
    gmatrix_model_ready = gmatrix_model_ready,
    env_similarity = ctx$env_similarity %||% NULL,
    env_ids = ctx$env_ids %||% NULL,
    env_covariates = ctx$env_covariates %||% NULL,
    reaction_norm_feature_qc = ctx$reaction_norm_feature_qc %||% TRUE,
    kenv_kernel = ctx$kenv_kernel %||% "matern32",
    kenv_bandwidth = ctx$kenv_bandwidth %||% 1.0,
    kenv_kernel_kwargs = ctx$kenv_kernel_kwargs %||% NULL,
    gp_backend = ctx$gp_backend %||% "auto",
    gp_engine = ctx$gp_engine %||% "auto",
    gp_learn_scales = ctx$gp_learn_scales %||% NULL,
    gp_varcomp_mode = ctx$gp_varcomp_mode %||% "reml",
    gp_fa_rank = ctx$gp_fa_rank %||% 1L,
    gp_factor_cache = ctx$gp_factor_cache %||% NULL,
    gp_iters = ctx$gp_iters %||% NULL,
    gp_lr = ctx$gp_lr %||% NULL,
    gp_exact_fast_cv = ctx$gp_exact_fast_cv %||% NULL,
    random_state = ctx$random_state %||% 12345L,
    fixed_effects = ctx$fixed_effects %||% NULL,
    include_components = ctx$include_components %||% NULL,
    krr_lam = ctx$krr_lam %||% NULL,
    krr_lams = ctx$krr_lams %||% NULL,
    lam_select = ctx$lam_select %||% NULL,
    large_n_threshold = ctx$large_n_threshold %||% NULL,
    operator_tol = ctx$operator_tol %||% NULL,
    operator_max_iter = ctx$operator_max_iter %||% NULL,
    operator_dtype_compute = ctx$operator_dtype_compute %||% NULL,
    gp_dtype = ctx$gp_dtype %||% NULL,
    gp_tiered_dispatch = ctx$gp_tiered_dispatch %||% NULL,
    gp_tiered_verbose = ctx$gp_tiered_verbose %||% NULL,
    gp_se_hutchinson_probes = ctx$gp_se_hutchinson_probes %||% NULL
  )

  cv_results <- vector("list", replication)
  gp_info_keep <- gp_cv_frontdoor_fast_info_keys()
  for (rep_i in seq_len(replication)) {
    split <- gp_cv_frontdoor_fast_splits(
      pheno_data = pheno_data,
      trait = trait,
      gen_name = ctx$gen_name,
      heter_groups = heter_groups,
      cv_token = cv_token,
      cross_validation_meth = ctx$cross_validation_meth,
      sampling_method = sampling_method,
      nfolds = ctx$nfolds,
      test_size = ctx$test_size %||% NULL,
      random_state = ctx$random_state,
      rep_i = rep_i,
      response_family = ctx$response_family %||% "gaussian",
      verbose = ctx$verbose %||% FALSE
    )
    fold_list <- lapply(split$folds, function(idx) as.integer(idx[!is.na(idx)]))
    fold_list <- fold_list[lengths(fold_list) > 0L]
    if (!length(fold_list)) {
      stop("GP front-door CV produced no non-empty folds.", call. = FALSE)
    }

    batch_preds <- gp_backend_cv_predict_batch(
      model_name = model,
      y = y,
      folds = fold_list,
      additional_params = params
    )
    if (!is.list(batch_preds) || length(batch_preds) != length(fold_list)) {
      stop("GP front-door CV returned an unexpected fold count.", call. = FALSE)
    }

    ypred_cv <- data.frame(
      row_id = seq_len(nrow(pheno_data)),
      y = y,
      yhat = rep(NA_real_, nrow(pheno_data)),
      cv_role = "train",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (isTRUE(is_multi_cv)) {
      ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])
    }
    for (ii in seq_along(fold_list)) {
      tst <- fold_list[[ii]]
      pred <- as.double(batch_preds[[ii]])
      if (length(pred) != length(tst)) {
        names_pred <- names(batch_preds[[ii]])
        if (!is.null(names_pred) && !is.null(rownames(pheno_data))) {
          row_match <- match(names_pred, rownames(pheno_data))
          keep <- match(tst, row_match)
          pred <- pred[keep]
        }
      }
      if (length(pred) != length(tst)) {
        stop("GP front-door CV prediction length did not match fold size.", call. = FALSE)
      }
      ypred_cv$cv_role[tst] <- "test"
      ypred_cv$yhat[tst] <- pred
    }
    ypred_cv$yhat[ypred_cv$cv_role != "test"] <- NA_real_
    if (!any(ypred_cv$cv_role == "test")) {
      stop("GP front-door CV marked no test rows.", call. = FALSE)
    }

    info <- attr(batch_preds, "gp_info", exact = TRUE)
    gp_info <- if (is.list(info)) info[intersect(gp_info_keep, names(info))] else list()
    metric_frame <- gp_cv_frontdoor_fast_metric_frame(
      ypred_cv = ypred_cv,
      rep_i = rep_i,
      eval_metrics = ctx$eval_metrics,
      response_family = ctx$response_family %||% "gaussian",
      heter_groups = heter_groups,
      is_multi = is_multi_cv
    )
    cv_results[[rep_i]] <- list(
      trait = trait,
      rep = as.integer(rep_i),
      model = model_label,
      model_canonical = model,
      eval_metrics_reps = metric_frame,
      results_eval_metrics_reps = metric_frame,
      ypred_cv_Reps_all = ypred_cv,
      yprob_cv_Reps_all = NULL,
      cv_info = c(
        list(
          method = ctx$cross_validation_meth,
          token = cv_token,
          nfolds = if (cv_token %in% c("k_folds", "stratified_k_folds", "repeated_k_folds",
                                       "repeated_stratified_k_folds", "leave_one_out")) {
            as.integer(split$nfolds)
          } else {
            NA_integer_
          },
          replication = as.integer(rep_i),
          n_test = sum(ypred_cv$cv_role == "test"),
          n_train = sum(ypred_cv$cv_role == "train"),
          gp_cv_frontdoor_fast_lane = TRUE
        ),
        gp_info
      )
    )
  }
  cv_results <- gp_attach_trait_prediction_test_info_to_cv_results(
    cv_results = cv_results,
    pheno_clean = ctx$pheno_clean,
    response = response,
    gen_name = ctx$gen_name
  )
  attr(cv_results, "gp_execution_policy") <- list(
    backend = "sequential",
    workers = 1L,
    plan = NULL,
    chunk_size = NULL,
    backend_score = Inf,
    decision_reason = "gp_cv_frontdoor_fast_lane",
    num_gpus = NA_integer_,
    model_ids = model,
    user_sequential_models = unique(ctx$sequential_models %||% character()),
    internal_flags = TRUE,
    internal_sequential_count = length(cv_results),
    parallel_count = 0L
  )
  cv_profile <- gp_runtime_profile_mark(
    cv_profile,
    "frontdoor_cv_execution",
    detail = paste("folded_results", length(cv_results))
  )

  cv_generate_plots <- gp_cv_generate_plots_enabled(ctx)
  if (isTRUE(cv_generate_plots) && isTRUE(is_multi_cv)) {
    cv_results_predicted_vs_observed <- multi_env_predicted_vs_observed_result_plots_process(
      results = cv_results,
      pheno_data = pheno_data,
      heter_groups = heter_groups,
      abs_very_close_threshold = ctx$abs_very_close_threshold,
      abs_close_threshold = ctx$abs_close_threshold
    )
    cv_results_processed <- cv1_cv2_and_across_env_result_plot_process(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      heter_groups = heter_groups,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  } else if (isTRUE(cv_generate_plots)) {
    cv_results_predicted_vs_observed <- single_predicted_vs_observed_result_plots_process(
      results = cv_results,
      pheno_data = pheno_data,
      abs_very_close_threshold = ctx$abs_very_close_threshold,
      abs_close_threshold = ctx$abs_close_threshold
    )
    cv_results_processed <- cv_single_loc_result_plot_process(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  } else {
    cv_results_predicted_vs_observed <- list(
      predicted_vs_observed_plots = NULL,
      mod_res_per_trait_per_model = NULL,
      classification_diagnostic_summaries = NULL
    )
    cv_results_processed <- gp_cv_process_results_fast(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      heter_groups = heter_groups,
      is_multi = is_multi_cv,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  }
  cv_profile <- gp_runtime_profile_mark(
    cv_profile,
    "frontdoor_cv_postprocess",
    detail = paste("traits", length(response))
  )
  cv_results_processed[["trait_prediction_test_summary"]] <- gp_trait_prediction_test_set_summary(
    pheno_clean = ctx$pheno_clean,
    response = response,
    gen_name = ctx$gen_name
  )

  cv_info_values <- function(key) {
    vals <- unlist(
      lapply(cv_results %||% list(), function(res) {
        info <- if (is.list(res)) res$cv_info else NULL
        if (is.list(info)) info[[key]] else NULL
      }),
      use.names = FALSE
    )
    vals <- vals[!is.na(vals)]
    vals <- vals[nzchar(as.character(vals))]
    unique(as.character(vals))
  }
  cv_gp_runtime_fields <- list(gp_cv_frontdoor_fast_lane = TRUE)
  for (key in c(gp_info_keep, "gp_cv_frontdoor_fast_lane")) {
    vals <- cv_info_values(key)
    if (length(vals)) {
      cv_gp_runtime_fields[[paste0("gp_", key)]] <- paste(vals, collapse = ",")
    }
  }

  metadata_python <- tryCatch(gp_detect_python(), error = function(e) NULL)
  cv_run_metadata <- gp_runtime_metadata(
    execution_policy = attr(cv_results, "gp_execution_policy"),
    python_path = metadata_python,
    preferred_python = metadata_python,
    include_accelerator = FALSE,
    context = "cross_validation_gp_frontdoor_fast_lane",
    extra_fields = c(
      gp_runtime_profile_metadata_fields(cv_profile),
      cv_gp_runtime_fields,
      list(
        cv_generate_plots = cv_generate_plots,
        response_family = ctx$response_family %||% "gaussian",
        metric_for_ranking = ctx$metric_for_ranking,
        ranking_tie_breakers = paste(ctx$ranking_tie_breakers %||% character(), collapse = ","),
        positive_class = ctx$positive_class %||% NA_character_,
        task_unit = "trait_model_replication",
        task_models = 1L,
        task_traits = 1L,
        task_replications = as.integer(replication)
      )
    )
  )
  cv_results_processed[["run_metadata"]] <- cv_run_metadata
  best_models <- cv_results_processed[["best_models_list"]][[ctx$metric_for_ranking]]
  best_models_ggplot_rep <- cv_results_processed[["plot_reps_list"]][[ctx$metric_for_ranking]][["ggplot_boxplot_reps"]]
  best_models_ggplot_mean <- cv_results_processed[["plot_mean_list"]][[ctx$metric_for_ranking]][["ggplot_lineplot_mean"]]

  list(
    cv_results = cv_results,
    cv_results_processed = cv_results_processed,
    cv_results_predicted_vs_observed = cv_results_predicted_vs_observed,
    run_metadata = cv_run_metadata,
    run_profile = gp_runtime_profile_table(cv_profile),
    best_models = best_models,
    best_models_ggplot_rep = best_models_ggplot_rep,
    best_models_ggplot_mean = best_models_ggplot_mean
  )
}

gp_run_cross_validation_pipeline <- function(ctx) {
  cv_profile <- gp_runtime_profile_new("cross_validation_pipeline")
  cv_generate_plots <- gp_cv_generate_plots_enabled(ctx)
  gmatrix_model_ready <- ctx$gmatrix_model_ready %||% ctx$gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] %||% NULL
  omic1_kernel_model_ready <- ctx$omic1_kernel_model_ready %||% ctx$gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] %||% NULL
  omic2_kernel_model_ready <- ctx$omic2_kernel_model_ready %||% ctx$gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] %||% NULL
  omic3_kernel_model_ready <- ctx$omic3_kernel_model_ready %||% ctx$gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] %||% NULL
  ctx$gp_output_level <- "predict_only"
  ctx$gp_return_se <- FALSE
  ctx$gp_full_vc <- FALSE
  ctx$gp_prediction_output <- "test_only"
  ctx$gp_return_trait_correlations <- FALSE

  cv_results <- tryCatch({
    cv_params <- list(
      test_set = ctx$test_set,
      response = ctx$response,
      gen_name = ctx$gen_name,
      test_size = ctx$test_size,
      random_state = ctx$random_state,
      replication = ctx$replication,
      selected_raw = ctx$feature_selected,
      max_features = ctx$max_features,
      feature_score_metadata = ctx$feature_score_metadata,
      feature_k_grid = ctx$feature_k_grid,
      feature_scoring_model = ctx$feature_scoring_model,
      feature_scoring_cv = ctx$feature_scoring_cv,
      feature_scoring_seed = ctx$feature_scoring_seed,
      feature_ridge_lambda = ctx$feature_ridge_lambda,
      feature_bayes_nIter = ctx$feature_bayes_nIter,
      feature_bayes_burnIn = ctx$feature_bayes_burnIn,
      feature_bayes_thin = ctx$feature_bayes_thin,
      feature_source_map = ctx$feature_source_map,
      feature_source_matrices = list(
        geno_data = ctx$geno_omic_model_ready_list[["geno_model_ready"]] %||% NULL,
        omic1_data = ctx$geno_omic_model_ready_list[["omic1_model_ready"]] %||% NULL,
        omic2_data = ctx$geno_omic_model_ready_list[["omic2_model_ready"]] %||% NULL,
        omic3_data = ctx$geno_omic_model_ready_list[["omic3_model_ready"]] %||% NULL
      ),
      feature_precomputed_kernel_supplied = any(vapply(
        list(ctx$gmatrix, ctx$gkernel, ctx$omic1_kernel, ctx$omic2_kernel, ctx$omic3_kernel, ctx$kernel_list),
        Negate(is.null),
        logical(1L)
      )),
      gmatrix_method = ctx$gmatrix_method,
      ploidy = ctx$ploidy %||% "auto",
      kernel_method = ctx$kernel_method,
      fixed = ctx$fixed,
      random = ctx$random,
      cova = ctx$cova,
      weights = ctx$weights,
      fixed_term_model_bayesian = ctx$fixed_term_model_bayesian,
      rand_term_model_bayesian = ctx$rand_term_model_bayesian,
      nIter = ctx$nIter,
      burnIn = ctx$burnIn,
      thin = ctx$thin,
      omics_data_label = ctx$omics_data_label,
      heter_resid = ctx$heter_resid,
      bayes_kernel_heter_resid = ctx$bayes_kernel_heter_resid,
      var_cov_str = ctx$var_cov_str,
      inverse = ctx$inverse,
      epsilon = ctx$epsilon,
      workspace = ctx$workspace,
      pworkspace = ctx$pworkspace,
      maxit = ctx$maxit,
      heter_groups = ctx$heter_groups,
      cross_validation_meth = ctx$cross_validation_meth,
      nfolds = ctx$nfolds,
      sampling_method = ctx$sampling_method,
      model_prep_all_bayes_cv = ctx$model_prep_all_bayes_cv,
      asreml_models_prep_cv = ctx$asreml_models_prep_cv,
      engine = ctx$engine,
      ml_dat_res = ctx$ml_dat_res,
      met_ml_dl = ctx$met_ml_dl,
      gmatrix_model_ready = gmatrix_model_ready,
      omic1_kernel_model_ready = omic1_kernel_model_ready,
      omic2_kernel_model_ready = omic2_kernel_model_ready,
      omic3_kernel_model_ready = omic3_kernel_model_ready,
      kernel_list = ctx$gmatrix_kernel_model_ready_list,
      omics_kernel_label = ctx$omics_kernel_label,
      met_kernel_var_explained = ctx$met_kernel_var_explained,
      met_kernel_min_ev = ctx$met_kernel_min_ev,
      met_kernel_max_pcs = ctx$met_kernel_max_pcs,
      env_similarity = ctx$env_similarity,
      env_ids = ctx$env_ids,
      env_covariates = ctx$env_covariates,
      reaction_norm_feature_qc = ctx$reaction_norm_feature_qc,
      kenv_kernel = ctx$kenv_kernel,
      kenv_bandwidth = ctx$kenv_bandwidth,
      kenv_kernel_kwargs = ctx$kenv_kernel_kwargs,
      GS_model_cv = ctx$GS_model_cv,
      crossval = TRUE,
      num_cores = ctx$num_cores,
      eval_metrics = ctx$eval_metrics,
      response_family = ctx$response_family,
      positive_class = ctx$positive_class,
      scaling = ctx$scaling,
      centering = ctx$centering,
      eta = ctx$learning_rate,
      nrounds = ctx$iteration,
      max_depth = ctx$max_depth,
      xgb_gamma = ctx$xgb_gamma,
      subsample = ctx$subsample,
      colsample_bytree = ctx$colsample_bytree,
      xgb_alpha = ctx$xgb_alpha,
      xgb_lambda = ctx$xgb_lambda,
      min_child_weight = ctx$min_child_weight,
      early_stop_for_iteration_xgb = ctx$early_stop_for_iteration_xgb,
      xgb_booster = ctx$xgb_booster,
      ncomp = ctx$ncomp,
      ntree = ctx$ntree,
      mtry = ctx$mtry,
      maxnodes = ctx$maxnodes,
      nodesize = ctx$nodesize,
      rf_n_jobs = ctx$rf_n_jobs,
      k = ctx$k,
      svm_kernel = ctx$svm_kernel,
      sigma_value = ctx$sigma_value,
      C_value = ctx$C_value,
      degree_value = ctx$degree_value,
      scale_value = ctx$scale_value,
      offset_value = ctx$offset_value,
      early_stop = ctx$early_stop,
      optimizer_name = ctx$optimizer_name,
      use_amp = ctx$use_amp,
      max_grad_norm = ctx$max_grad_norm,
      auto_class_weights = ctx$auto_class_weights,
      cnn_neurons_per_layer = ctx$cnn_neurons_per_layer,
      cnn_kernel_size = ctx$cnn_kernel_size,
      cnn_dense_layers = ctx$cnn_dense_layers,
      cnn_use_max_pool = ctx$cnn_use_max_pool,
      cnn_pool_kernel = ctx$cnn_pool_kernel,
      cnn_pool_stride = ctx$cnn_pool_stride,
      cnn_pool_padding = ctx$cnn_pool_padding,
      cnn_learning_rate = ctx$cnn_learning_rate,
      cnn_separable = ctx$cnn_separable,
      cnn_dilations = ctx$cnn_dilations,
      cnn_use_se = ctx$cnn_use_se,
      cnn_norm_type = ctx$cnn_norm_type,
      cnn_pool_type = ctx$cnn_pool_type,
      cnn_use_global_pool = ctx$cnn_use_global_pool,
      resnet_neurons_per_block = ctx$resnet_neurons_per_block,
      resnet_blocks = ctx$resnet_blocks,
      resnet_learning_rate = ctx$resnet_learning_rate,
      ft_d_model = ctx$ft_d_model,
      ft_heads = ctx$ft_heads,
      ft_layers = ctx$ft_layers,
      ft_ff_mult = ctx$ft_ff_mult,
      ft_dropout = ctx$ft_dropout,
      ft_token_dropout = ctx$ft_token_dropout,
      ft_use_cls = ctx$ft_use_cls,
      saint_d_model = ctx$saint_d_model,
      saint_heads = ctx$saint_heads,
      saint_layers = ctx$saint_layers,
      saint_ff_mult = ctx$saint_ff_mult,
      saint_dropout = ctx$saint_dropout,
      saint_token_dropout = ctx$saint_token_dropout,
      saint_use_cls = ctx$saint_use_cls,
      use_grouping = ctx$use_grouping,
      group_trigger = ctx$group_trigger,
      group_method = ctx$group_method,
      init_group_size = ctx$init_group_size,
      max_tokens = ctx$max_tokens,
      kmeans_batch = ctx$kmeans_batch,
      kmeans_iter = ctx$kmeans_iter,
      tabnet_steps = ctx$tabnet_steps,
      tabnet_feature_dim = ctx$tabnet_feature_dim,
      tabnet_output_dim = ctx$tabnet_output_dim,
      tabnet_gamma = ctx$tabnet_gamma,
      tabnet_lambda_sparse = ctx$tabnet_lambda_sparse,
      node_trees = ctx$node_trees,
      node_depth = ctx$node_depth,
      deepfm_k = ctx$deepfm_k,
      deepfm_hidden = ctx$deepfm_hidden,
      dcn_layers = ctx$dcn_layers,
      dcn_hidden = ctx$dcn_hidden,
      nam_hidden = ctx$nam_hidden,
      nam_activation = ctx$nam_activation,
      nam_add_linear = ctx$nam_add_linear,
      nam_l1 = ctx$nam_l1,
      moe_n_experts = ctx$moe_n_experts,
      moe_expert_hidden = ctx$moe_expert_hidden,
      moe_gate_hidden = ctx$moe_gate_hidden,
      moe_temperature = ctx$moe_temperature,
      moe_sparse_topk = ctx$moe_sparse_topk,
      moe_entropy_reg = ctx$moe_entropy_reg,
      gp_use_variational = ctx$gp_use_variational,
      gp_num_inducing = ctx$gp_num_inducing,
      gp_feature_dim = ctx$gp_feature_dim,
      gp_kernel = ctx$gp_kernel,
      gp_ard = ctx$gp_ard,
      gp_lr_mult = ctx$gp_lr_mult,
      rff_features = ctx$rff_features,
      rff_lengthscale = ctx$rff_lengthscale,
      rff_deep_hidden = ctx$rff_deep_hidden,
      model_type = ctx$model_type,
      epochs = ctx$epochs,
      batch_size = ctx$batch_size,
      dropout = ctx$dropout,
      l2_weight_decay = ctx$l2_weight_decay,
      l2_regularizer_dp = ctx$l2_regularizer_dp,
      dropout_rate = ctx$dropout_rate,
      batch_norm = ctx$batch_norm,
      validation_split = ctx$validation_split,
      compile_model = ctx$compile_model,
      deterministic = ctx$deterministic,
      random_seed = ctx$random_seed,
      device = ctx$device,
      mlp_neurons_per_layer = ctx$mlp_neurons_per_layer,
      mlp_learning_rate = ctx$mlp_learning_rate,
      final_attention = ctx$final_attention,
      attention_across_multiple_layers = ctx$attention_across_multiple_layers,
      heteroscedastic = ctx$heteroscedastic,
      gp_return_se = ctx$gp_return_se,
      gp_full_vc = ctx$gp_full_vc,
      gp_backend = ctx$gp_backend,
      gp_engine = ctx$gp_engine %||% "auto",
      gp_learn_scales = ctx$gp_learn_scales %||% NULL,
      gp_output_level = ctx$gp_output_level,
      gp_varcomp_mode = ctx$gp_varcomp_mode,
      gp_fa_rank = ctx$gp_fa_rank,
      gp_prediction_output = ctx$gp_prediction_output,
      gp_factor_cache = ctx$gp_factor_cache,
      gp_iters = ctx$gp_iters,
      gp_lr = ctx$gp_lr,
      gp_exact_fast_cv = ctx$gp_exact_fast_cv,
      docker_nd_usage = ctx$docker_nd_usage,
      globals_max_GB = ctx$globals_max_GB,
      worker_memory_gb = ctx$worker_memory_gb,
      memory_budget_gb = ctx$memory_budget_gb,
      sequential_models = ctx$sequential_models,
      parallel_mode = ctx$parallel_mode,
      parallel_backend_prefer_fork = ctx$parallel_backend_prefer_fork,
      GS_model_cv_display = ctx$GS_model_cv_display
    )

    models_execute_crossval(
      pheno_data = ctx$pheno_data,
      params = cv_params,
      verbose = ctx$verbose %||% TRUE
    )
  }, error = function(e) {
    # Surface engine failures (e.g. a parallel backend refusing the job)
    # instead of returning NULL results with only a warning.
    stop(
      paste(ctx$msg, "Cross-validation failed:", conditionMessage(e),
            "\nIf this came from the parallel backend, rerun with parallel_mode = \"sequential\" or \"base_parallel\"."),
      call. = FALSE
    )
  })
  cv_profile <- gp_runtime_profile_mark(
    cv_profile,
    "cv_execution",
    detail = paste("results", length(cv_results %||% list()))
  )

  if (length(cv_results) == 0) {
    warning(
      paste(
        ctx$msg,
        paste(
          "Error in cross-validation.",
          "No result for cross-validation.",
          "Check the model and the data to ensure the data is correct and the model(s) is well specified."
        )
      ),
      call. = FALSE
    )
    return(NULL)
  }

  cv_results <- gp_attach_trait_prediction_test_info_to_cv_results(
    cv_results = cv_results,
    pheno_clean = ctx$pheno_clean,
    response = ctx$response,
    gen_name = ctx$gen_name
  )

  is_multi_cv <- ctx$cross_validation_meth %in% c("CV0", "CV1", "CV2", "Repeated_CV0", "Repeated_CV1", "Repeated_CV2")
  if (isTRUE(cv_generate_plots) && is_multi_cv) {
    cv_results_predicted_vs_observed <- multi_env_predicted_vs_observed_result_plots_process(
      results = cv_results,
      pheno_data = ctx$pheno_data,
      heter_groups = ctx$heter_groups,
      abs_very_close_threshold = ctx$abs_very_close_threshold,
      abs_close_threshold = ctx$abs_close_threshold
    )
    cv_results_processed <- cv1_cv2_and_across_env_result_plot_process(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      heter_groups = ctx$heter_groups,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  } else if (isTRUE(cv_generate_plots)) {
    cv_results_predicted_vs_observed <- single_predicted_vs_observed_result_plots_process(
      results = cv_results,
      pheno_data = ctx$pheno_data,
      abs_very_close_threshold = ctx$abs_very_close_threshold,
      abs_close_threshold = ctx$abs_close_threshold
    )
    cv_results_processed <- cv_single_loc_result_plot_process(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  } else {
    cv_results_predicted_vs_observed <- list(
      predicted_vs_observed_plots = NULL,
      mod_res_per_trait_per_model = NULL,
      classification_diagnostic_summaries = NULL
    )
    cv_results_processed <- gp_cv_process_results_fast(
      cv_results_data = cv_results,
      eval_metrics = ctx$eval_metrics,
      heter_groups = ctx$heter_groups,
      is_multi = is_multi_cv,
      metric_for_ranking = ctx$metric_for_ranking,
      ranking_tie_breakers = ctx$ranking_tie_breakers,
      positive_class = ctx$positive_class
    )
  }
  cv_profile <- gp_runtime_profile_mark(
    cv_profile,
    "cv_postprocess",
    detail = paste("traits", length(ctx$response %||% character()))
  )
  cv_results_processed[["trait_prediction_test_summary"]] <- gp_trait_prediction_test_set_summary(
    pheno_clean = ctx$pheno_clean,
    response = ctx$response,
    gen_name = ctx$gen_name
  )

  cv_gp_runtime_fields <- list()
  cv_info_values <- function(key) {
    vals <- unlist(
      lapply(cv_results %||% list(), function(res) {
        info <- if (is.list(res)) res$cv_info else NULL
        if (is.list(info)) info[[key]] else NULL
      }),
      use.names = FALSE
    )
    vals <- vals[!is.na(vals)]
    vals <- vals[nzchar(as.character(vals))]
    unique(as.character(vals))
  }
  cv_info_keys <- c(
    "package_cv_batch", "bridge_command", "cv_fast_path",
    "cv_fast_path_skip_reason", "cv_fast_path_approximate",
    "cv_fast_path_exact_refit", "cv_fast_path_solver",
    "cv_fast_path_shared_hyperparameters",
    "cv_fast_path_likelihood_noise_variance",
    "cv_fast_path_likelihood_noise_std", "cv_fast_path_self_checked",
    "cv_fast_path_self_check_reason", "cv_preprocessed_context",
    "cv_preprocessed_context_rows",
    "cv_preprocessed_context_tensor_cache_entries", "gp_exact_fast_cv"
  )
  for (key in cv_info_keys) {
    vals <- cv_info_values(key)
    if (length(vals)) {
      cv_gp_runtime_fields[[paste0("gp_", key)]] <- paste(vals, collapse = ",")
    }
  }

  metadata_python <- tryCatch(gp_detect_python(), error = function(e) NULL)
  cv_run_metadata <- gp_runtime_metadata(
    execution_policy = attr(cv_results, "gp_execution_policy"),
    python_path = metadata_python,
    preferred_python = metadata_python,
    include_accelerator = FALSE,
    context = "cross_validation",
    extra_fields = c(
      gp_runtime_profile_metadata_fields(cv_profile),
      cv_gp_runtime_fields,
      list(
        cv_generate_plots = cv_generate_plots,
        response_family = ctx$response_family %||% "gaussian",
        metric_for_ranking = ctx$metric_for_ranking,
        ranking_tie_breakers = paste(ctx$ranking_tie_breakers %||% character(), collapse = ","),
        positive_class = ctx$positive_class %||% NA_character_,
        task_unit = "trait_model_replication",
        task_models = length(unique(ctx$GS_model_cv %||% character())),
        task_traits = length(ctx$response %||% character()),
        task_replications = as.integer(ctx$replication %||% 1L)
      )
    )
  )
  if (!is.null(cv_results_predicted_vs_observed$classification_diagnostic_summaries)) {
    cv_results_processed[["classification_diagnostic_summaries"]] <-
      cv_results_predicted_vs_observed$classification_diagnostic_summaries
  }
  if (!is.null(ctx$feature_score_metadata)) {
    cv_results_processed[["feature_score_metadata"]] <- ctx$feature_score_metadata
  }
  if (!is.null(ctx$feature_score_metadata) || !is.null(ctx$feature_selected)) {
    cv_results_processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
      lapply(cv_results, function(x) x[["feature_selection_metadata"]])
    )
    cv_results_processed[["feature_k_grid"]] <- ctx$feature_k_grid
    cv_results_processed[["feature_scoring_model"]] <- if (is.null(ctx$feature_score_metadata)) {
      "user_supplied"
    } else {
      ctx$feature_scoring_model
    }
    cv_results_processed[["feature_scoring_cv"]] <- if (is.null(ctx$feature_score_metadata)) {
      "explicit"
    } else {
      ctx$feature_scoring_cv
    }
  }
  cv_results_processed[["run_metadata"]] <- cv_run_metadata

  best_models <- cv_results_processed[["best_models_list"]][[ctx$metric_for_ranking]]
  best_models_ggplot_rep <- cv_results_processed[["plot_reps_list"]][[ctx$metric_for_ranking]][["ggplot_boxplot_reps"]]
  best_models_ggplot_mean <- cv_results_processed[["plot_mean_list"]][[ctx$metric_for_ranking]][["ggplot_lineplot_mean"]]

  list(
    cv_results = cv_results,
    cv_results_processed = cv_results_processed,
    cv_results_predicted_vs_observed = cv_results_predicted_vs_observed,
    run_metadata = cv_run_metadata,
    run_profile = gp_runtime_profile_table(cv_profile),
    best_models = best_models,
    best_models_ggplot_rep = best_models_ggplot_rep,
    best_models_ggplot_mean = best_models_ggplot_mean
  )
}

gp_prepare_cross_validation_artifacts <- function(ctx) {
  pheno_data <- ctx$pheno_clean[["pheno_clean_data"]]
  if (!is.null(ctx$test_set)) {
    pheno_data <- pheno_data[!pheno_data[[ctx$gen_name]] %in% ctx$test_set, , drop = FALSE]
  }

  geno_model_ready <- ctx$geno_omic_model_ready_list[["geno_model_ready"]] %||% NULL
  omic1_model_ready <- ctx$geno_omic_model_ready_list[["omic1_model_ready"]] %||% NULL
  omic2_model_ready <- ctx$geno_omic_model_ready_list[["omic2_model_ready"]] %||% NULL
  omic3_model_ready <- ctx$geno_omic_model_ready_list[["omic3_model_ready"]] %||% NULL
  gmatrix_model_ready <- ctx$gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] %||% NULL
  omic1_kernel_model_ready <- ctx$gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] %||% NULL
  omic2_kernel_model_ready <- ctx$gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] %||% NULL
  omic3_kernel_model_ready <- ctx$gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] %||% NULL

  model_prep_all_bayes_cv <- NULL
  if (any(ctx$GS_model_cv %in% c(ctx$bayes_valid_models, ctx$bayes_gblup_valid_models))) {
    model_prep_all_bayes_cv <- model_prep_bayes_cv(
      fixed = ctx$fixed,
      random = ctx$random,
      GS_model_cv = ctx$GS_model_cv,
      response = ctx$response,
      gen_name = ctx$gen_name,
      pheno_data = pheno_data,
      test_set = ctx$test_set,
      weights = ctx$weights,
      fixed_term_model_bayesian = ctx$fixed_term_model_bayesian,
      rand_term_model_bayesian = ctx$rand_term_model_bayesian,
      nIter = ctx$nIter,
      burnIn = ctx$burnIn,
      thin = ctx$thin,
      geno_data = geno_model_ready,
      omic1_data = omic1_model_ready,
      omic2_data = omic2_model_ready,
      omic3_data = omic3_model_ready,
      omics_data_label = ctx$omics_data_label,
      gmatrix = gmatrix_model_ready,
      omic1_kernel = omic1_kernel_model_ready,
      omic2_kernel = omic2_kernel_model_ready,
      omic3_kernel = omic3_kernel_model_ready,
      kernel_list = ctx$gmatrix_kernel_model_ready_list,
      heter_groups = ctx$heter_groups,
      heter_resid = ctx$heter_resid,
      bayes_kernel_heter_resid = ctx$bayes_kernel_heter_resid,
      omics_kernel_label = ctx$omics_kernel_label,
      cross_validation = TRUE,
      response_family = ctx$response_family
    )
  }

  asreml_models_prep_cv <- NULL
  if (any(ctx$GS_model_cv %in% "GBLUP")) {
    asreml_models_prep_cv <- asreml_utilis_new(
      fixed = ctx$fixed,
      random = ctx$random,
      engine = ctx$engine,
      cova = ctx$cova,
      GS_model = "GBLUP",
      response = ctx$response,
      pheno_data = pheno_data,
      gmatrix = gmatrix_model_ready,
      omic1_kernel = omic1_kernel_model_ready,
      omic2_kernel = omic2_kernel_model_ready,
      omic3_kernel = omic3_kernel_model_ready,
      kernel_list = ctx$gmatrix_kernel_model_ready_list,
      inverse = ctx$inverse,
      epsilon = ctx$epsilon,
      gen_name = ctx$gen_name,
      heter_groups = ctx$heter_groups,
      heter_resid = ctx$heter_resid,
      var_cov_str = ctx[["var_cov_str"]],
      weights = ctx$weights,
      workspace = ctx$workspace,
      pworkspace = ctx$pworkspace,
      maxit = ctx$maxit,
      cross_validation = TRUE
    )
  }

  list(
    pheno_data = pheno_data,
    model_prep_all_bayes_cv = model_prep_all_bayes_cv,
    asreml_models_prep_cv = asreml_models_prep_cv
  )
}

gp_cv_metric_minimize_metrics <- function() {
  c(
    "mean_squared_error", "bias", "root_mean_squared_error",
    "relative_squared_error", "mean_absolute_error",
    "mean_absolute_percent_error", "log_loss", "brier_score",
    "mean_absolute_error_class", "ece"
  )
}

gp_cv_metric_source_frame <- function(cv_results_processed, metric_for_ranking = "accuracy") {
  aggregated <- cv_results_processed[["aggregated_data_list"]]
  if (!is.list(aggregated)) {
    return(data.frame())
  }
  if (is.data.frame(aggregated[["aggregated_across_reps"]])) {
    return(aggregated[["aggregated_across_reps"]])
  }
  metric_entry <- aggregated[[metric_for_ranking]]
  if (is.list(metric_entry)) {
    for (nm in c("aggregated_data_env_mean", "aggregated_across_reps", "aggregated_data")) {
      if (is.data.frame(metric_entry[[nm]])) {
        return(metric_entry[[nm]])
      }
    }
  }
  data_frames <- Filter(is.data.frame, aggregated)
  if (length(data_frames)) {
    return(data_frames[[1L]])
  }
  data.frame()
}

gp_best_cv_row_per_trait_model <- function(metric_data, metric_for_ranking = "accuracy") {
  if (!is.data.frame(metric_data) ||
      !all(c("trait", "model") %in% names(metric_data)) ||
      !nrow(metric_data)) {
    return(data.frame())
  }
  data <- as.data.frame(metric_data, stringsAsFactors = FALSE, check.names = FALSE)
  data[["trait"]] <- as.character(data[["trait"]])
  data[["model_canonical"]] <- gp_canonicalize_supported_model_names(data[["model"]])
  data[["model_canonical"]] <- as.character(data[["model_canonical"]])

  if (metric_for_ranking %in% names(data)) {
    metric_value <- suppressWarnings(as.numeric(data[[metric_for_ranking]]))
    metric_value[!is.finite(metric_value)] <- NA_real_
    if (metric_for_ranking %in% gp_cv_metric_minimize_metrics()) {
      ord <- order(data[["trait"]], data[["model_canonical"]], metric_value, na.last = TRUE)
    } else {
      ord <- order(data[["trait"]], data[["model_canonical"]], -metric_value, na.last = TRUE)
    }
    data <- data[ord, , drop = FALSE]
  } else {
    data <- data[order(data[["trait"]], data[["model_canonical"]]), , drop = FALSE]
  }

  group_key <- paste(data[["trait"]], data[["model_canonical"]], sep = "\r")
  data[!duplicated(group_key), , drop = FALSE]
}

gp_best_cv_row_per_trait <- function(metric_data, metric_for_ranking = "accuracy") {
  if (!is.data.frame(metric_data) ||
      !all(c("trait", "model") %in% names(metric_data)) ||
      !nrow(metric_data)) {
    return(data.frame())
  }
  data <- as.data.frame(metric_data, stringsAsFactors = FALSE, check.names = FALSE)
  data[["trait"]] <- as.character(data[["trait"]])
  data[["model_canonical"]] <- gp_canonicalize_supported_model_names(data[["model"]])
  data[["model_canonical"]] <- as.character(data[["model_canonical"]])

  if (metric_for_ranking %in% names(data)) {
    metric_value <- suppressWarnings(as.numeric(data[[metric_for_ranking]]))
    metric_value[!is.finite(metric_value)] <- NA_real_
    if (metric_for_ranking %in% gp_cv_metric_minimize_metrics()) {
      ord <- order(data[["trait"]], metric_value, data[["model_canonical"]], na.last = TRUE)
    } else {
      ord <- order(data[["trait"]], -metric_value, data[["model_canonical"]], na.last = TRUE)
    }
    data <- data[ord, , drop = FALSE]
  } else {
    data <- data[order(data[["trait"]], data[["model_canonical"]]), , drop = FALSE]
  }

  data[!duplicated(data[["trait"]]), , drop = FALSE]
}

gp_best_models_from_cv_processed <- function(cv_results_processed,
                                             metric_for_ranking = "accuracy") {
  best_list <- cv_results_processed[["best_models_list"]]
  if (is.list(best_list)) {
    metric_names <- unique(c(metric_for_ranking, names(best_list)))
    metric_names <- metric_names[nzchar(metric_names)]
    for (nm in metric_names) {
      candidate <- best_list[[nm]]
      if (is.data.frame(candidate) &&
          all(c("trait", "model") %in% names(candidate)) &&
          nrow(candidate)) {
        return(as.data.frame(candidate, stringsAsFactors = FALSE, check.names = FALSE))
      }
    }
  }

  metric_source <- gp_cv_metric_source_frame(cv_results_processed, metric_for_ranking)
  best <- gp_best_cv_row_per_trait(metric_source, metric_for_ranking)
  if (nrow(best)) {
    return(best)
  }
  NULL
}

gp_reduce_cv_best_models <- function(best_models,
                                     metric_for_ranking = "accuracy") {
  if (!is.data.frame(best_models) ||
      !all(c("trait", "model") %in% names(best_models)) ||
      !nrow(best_models)) {
    return(best_models)
  }

  data <- as.data.frame(best_models, stringsAsFactors = FALSE, check.names = FALSE)
  data[["trait"]] <- as.character(data[["trait"]])
  data[["model"]] <- as.character(data[["model"]])

  if (!anyDuplicated(data[["trait"]])) {
    return(data)
  }

  if ("rank" %in% names(data)) {
    rank_value <- suppressWarnings(as.numeric(data[["rank"]]))
    if (any(is.finite(rank_value))) {
      data[[".predictpror_rank_value"]] <- rank_value
      data <- data[order(data[["trait"]], data[[".predictpror_rank_value"]], na.last = TRUE), , drop = FALSE]
      data <- data[!duplicated(data[["trait"]]), , drop = FALSE]
      data[[".predictpror_rank_value"]] <- NULL
      rownames(data) <- NULL
      return(data)
    }
  }

  if (!is.null(metric_for_ranking) && metric_for_ranking %in% names(data)) {
    reduced <- gp_best_cv_row_per_trait(data, metric_for_ranking)
    if (is.data.frame(reduced) && nrow(reduced)) {
      rownames(reduced) <- NULL
      return(reduced)
    }
  }

  data <- data[!duplicated(data[["trait"]]), , drop = FALSE]
  rownames(data) <- NULL
  data
}

gp_build_true_prediction_task_table <- function(best_models = NULL,
                                                GS_model = NULL,
                                                GS_model_cv = NULL,
                                                response = NULL,
                                                cv_results_processed = NULL,
                                                metric_for_ranking = "accuracy",
                                                cross_validation = FALSE) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  response <- as.character(response %||% character())
  requested_cv_models <- if (isTRUE(cross_validation) && length(GS_model_cv %||% character()) > 0L) {
    gp_canonicalize_supported_model_names(GS_model_cv)
  } else {
    character()
  }
  requested_cv_models <- unique(as.character(requested_cv_models[!is.na(requested_cv_models) & nzchar(requested_cv_models)]))

  # Previously, when *all* requested CV models were in the GP-lowrank family
  # (Kernel-GBLUP / Gaussian-Process-GBLUP / FA-GBLUP / Scalable-GBLUP), this
  # function expanded the post-CV true-prediction task table into a
  # (trait x model) cross-product so every GP model was re-fit on every
  # trait. That made the GP family behave differently from every other
  # family (Bayes, ASReml, ML, DL all take only the best model per trait
  # after CV) and turned a 4-trait/4-GP-model run into 16 final fits
  # instead of 4. The cross-product branch has been removed so every
  # family follows the same "one best model per trait" contract.

  best_models <- if (is.data.frame(best_models)) {
    as.data.frame(best_models, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    NULL
  }
  if (is.null(best_models) && isTRUE(cross_validation)) {
    best_models <- gp_best_models_from_cv_processed(
      cv_results_processed = cv_results_processed,
      metric_for_ranking = metric_for_ranking
    )
  }
  if (!is.null(best_models) && isTRUE(cross_validation)) {
    best_models <- gp_reduce_cv_best_models(
      best_models = best_models,
      metric_for_ranking = metric_for_ranking
    )
  }

  if (!is.null(best_models) && all(c("trait", "model") %in% names(best_models))) {
    tasks <- best_models[, c("trait", "model"), drop = FALSE]
    tasks[["trait"]] <- as.character(tasks[["trait"]])
    tasks[["model"]] <- gp_canonicalize_supported_model_names(tasks[["model"]])
  } else if (isTRUE(cross_validation) && length(GS_model_cv %||% character()) > 0L) {
    requested_models <- requested_cv_models
    if (!length(response) && !is.null(best_models) && "trait" %in% names(best_models)) {
      response <- unique(as.character(best_models[["trait"]]))
    }
    if (length(requested_models) == 1L) {
      tasks <- data.frame(
        trait = response,
        model = requested_models[[1L]],
        stringsAsFactors = FALSE
      )
    } else {
      stop(
        paste0(
          "\n==================================================\n",
          "Unable to determine one CV-selected best model per trait. ",
          "True prediction after cross-validation requires a best_models table ",
          "or cv_results_processed$best_models_list[[metric_for_ranking]]."
        ),
        call. = FALSE
      )
    }
  } else {
    model <- gp_canonicalize_supported_model_names(GS_model)
    if (length(model) > 1L && length(response) > 1L && length(model) == length(response)) {
      tasks <- data.frame(trait = response, model = as.character(model), stringsAsFactors = FALSE)
    } else {
      task_rows <- lapply(response, function(trait_i) {
        data.frame(trait = trait_i, model = as.character(model), stringsAsFactors = FALSE)
      })
      tasks <- if (length(task_rows)) do.call(rbind, task_rows) else data.frame(trait = character(), model = character())
    }
  }

  if (!nrow(tasks)) {
    return(tasks)
  }

  tasks <- as.data.frame(tasks, stringsAsFactors = FALSE, check.names = FALSE)
  tasks[["trait"]] <- as.character(tasks[["trait"]])
  tasks[["model"]] <- gp_canonicalize_supported_model_names(tasks[["model"]])
  tasks[["model"]] <- as.character(tasks[["model"]])
  tasks <- tasks[!is.na(tasks[["trait"]]) & nzchar(tasks[["trait"]]) &
                   !is.na(tasks[["model"]]) & nzchar(tasks[["model"]]), , drop = FALSE]
  if (!nrow(tasks)) {
    return(tasks)
  }
  tasks <- unique(tasks)

  metric_source <- gp_best_cv_row_per_trait_model(
    gp_cv_metric_source_frame(cv_results_processed, metric_for_ranking),
    metric_for_ranking
  )
  if (nrow(metric_source)) {
    for (i in seq_len(nrow(tasks))) {
      idx <- which(metric_source[["trait"]] == tasks[["trait"]][[i]] &
                     metric_source[["model_canonical"]] == tasks[["model"]][[i]])
      if (!length(idx)) {
        next
      }
      source_row <- metric_source[idx[[1L]], , drop = FALSE]
      copy_cols <- setdiff(
        names(source_row),
        c("trait", "model", "model_canonical")
      )
      for (nm in copy_cols) {
        if (!nm %in% names(tasks)) {
          tasks[[nm]] <- NA
        }
        tasks[[nm]][[i]] <- source_row[[nm]][[1L]]
      }
      if (!"cv_model_label" %in% names(tasks)) {
        tasks[["cv_model_label"]] <- NA_character_
      }
      tasks[["cv_model_label"]][[i]] <- as.character(source_row[["model"]][[1L]])
    }
  }

  if (!is.null(best_models) && all(c("trait", "model") %in% names(best_models))) {
    # After the reduction above best_models is already one row per trait
    # (the CV winner), so it is itself the best-reference used to flag
    # is_cv_best on each task row. Pre-0.20.9 a separate cv_best_reference
    # existed only because the GP-family branch expanded into a
    # trait x model cross-product; that branch has been removed so the
    # reference and the task list coincide for every family.
    best_key_source <- best_models
  } else {
    best_key_source <- NULL
  }
  best_key <- if (!is.null(best_key_source) && all(c("trait", "model") %in% names(best_key_source))) {
    paste(
      as.character(best_key_source[["trait"]]),
      as.character(gp_canonicalize_supported_model_names(best_key_source[["model"]])),
      sep = "\r"
    )
  } else {
    character()
  }
  task_key <- paste(tasks[["trait"]], tasks[["model"]], sep = "\r")
  tasks[["is_cv_best"]] <- task_key %in% best_key
  tasks[["model_public"]] <- gp_public_model_label(tasks[["model"]])
  tasks[["model_key"]] <- make.names(tasks[["model_public"]])
  tasks[["task_name"]] <- make.unique(make.names(paste(tasks[["trait"]], tasks[["model_public"]], sep = "_")))
  rownames(tasks) <- NULL
  tasks
}

gp_execute_best_model_tasks <- function(best_models,
                                        run_one_best,
                                        sequential_models = NULL,
                                        canonical_names,
                                        friendly_names,
                                        policy_param_source = list(),
                                        parallel_mode,
                                        num_cores,
                                        globals_max_GB,
                                        verbose,
                                        parallel_backend_prefer_fork,
                                        pheno_clean,
                                        ml_dat_res,
                                        gmatrix_kernel_model_ready_list,
                                        geno_omic_model_ready_list,
                                        response,
                                        init_py,
                                        worker_memory_gb = NULL,
                                        memory_budget_gb = NULL) {
  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (nrow(best_models) == 0L) {
    results <- list()
    names(results) <- response
    return(results)
  }

  user_seq_opt <- getOption("gp.force.sequential.models", NULL)
  seq_models <- unique(c(sequential_models %||% character(), user_seq_opt %||% character()))
  if (length(seq_models) > 0) {
    seq_models <- gp_canonicalize_supported_model_names(seq_models)
  } else {
    seq_models <- character()
  }

  globals_for_size <- list(
    pheno_clean = pheno_clean,
    ml_dat_res = ml_dat_res,
    gmatrix_kernel_model_ready_list = gmatrix_kernel_model_ready_list,
    geno_omic_model_ready_list = geno_omic_model_ready_list
  )
  direct_cli_models <- unique(c(
    "Xgboost", "RandomForest", "CatBoost", "LightGBM", "PartialLeastSquare",
    "SupportVectorMachine", "K-NearestNeighbors", "Lasso", "Ridge_Regression",
    canonical_names,
    tryCatch(gp_lowrank_supported_models(), error = function(e) character()),
    # BGLR/ASReml run in R and must not trigger an embedded-Python init.
    gp_parallel_r_only_models()
  ))
  needs_initializer <- any(!tolower(as.character(best_models$model)) %in% tolower(direct_cli_models))

  results <- gp_execute_partitioned_tasks(
    model_ids = best_models$model,
    run_task = run_one_best,
    sequential_models = seq_models,
    decision_args = list(
      model_params_list = lapply(best_models$model, function(mod) {
        gp_parallel_extract_model_policy_params(mod, policy_param_source)
      }),
      globals = globals_for_size,
      user_mode = parallel_mode %||% "auto",
      num_cores = num_cores,
      globals_max_GB = globals_max_GB %||% 4,
      worker_memory_gb = worker_memory_gb,
      memory_budget_gb = memory_budget_gb,
      # One full-data fit per trait x model task.
      task_cost = tryCatch(list(
        n_train = NROW(pheno_clean[["pheno_clean_data"]]),
        n_markers = NCOL(ml_dat_res[["merged_data"]][["merge_data"]] %||%
                           gmatrix_kernel_model_ready_list[[1L]] %||% matrix(0, 1, 1000)),
        nIter = policy_param_source$nIter %||% 6000,
        ntree = policy_param_source$ntree %||% 500,
        iterations = policy_param_source$iteration %||% 300,
        folds_per_task = 1
      ), error = function(e) NULL),
      verbose = isTRUE(verbose) %||% TRUE,
      sys_name = Sys.info()[["sysname"]],
      prefer_fork = gp_parallel_safe_fork_preference(
        prefer_fork = parallel_backend_prefer_fork,
        needs_process_initializer = needs_initializer
      )
    ),
    packages = c("dplyr"),
    seed = TRUE,
    initializer = if (isTRUE(needs_initializer)) init_py else NULL
  )

  if (length(results) != nrow(best_models)) {
    stop(
      "Internal true-prediction task/result count mismatch; refusing to label results.",
      call. = FALSE
    )
  }
  for (i in seq_len(nrow(best_models))) {
    result_i <- results[[i]]
    if (is.null(result_i) || !is.list(result_i)) {
      next
    }
    has_task_identity <- any(c("trait", "response", "model", "GS_model") %in% names(result_i))
    if (!isTRUE(has_task_identity)) {
      # Test doubles and third-party wrappers may intentionally return an
      # opaque payload. Production run_one_best() results always carry both
      # identities and are checked below.
      next
    }
    expected_trait <- as.character(best_models[["trait"]][[i]])
    expected_model <- as.character(gp_canonicalize_supported_model_names(
      best_models[["model"]][[i]]
    ))
    actual_trait <- as.character(result_i[["trait"]] %||% result_i[["response"]] %||% "")
    actual_model <- as.character(gp_canonicalize_supported_model_names(
      result_i[["model"]] %||% result_i[["GS_model"]] %||% ""
    ))
    if (!identical(actual_trait, expected_trait) ||
        !identical(actual_model, expected_model)) {
      stop(
        sprintf(
          paste0(
            "Internal true-prediction task/result alignment failure at task %d: ",
            "expected %s/%s but received %s/%s. Results were not labelled or exported."
          ),
          i,
          expected_trait,
          expected_model,
          actual_trait,
          actual_model
        ),
        call. = FALSE
      )
    }
  }

  names(results) <- best_models[["trait"]]
  results
}

gp_merge_task_context <- function(parent_env, local_env) {
  utils::modifyList(
    as.list(parent_env),
    as.list(local_env),
    keep.null = TRUE
  )
}

gp_guard_asreml_cv0_fixed_environment <- function(cross_validation = FALSE,
                                                  cross_validation_meth = NULL,
                                                  GS_model_cv = NULL,
                                                  engine = NULL,
                                                  fixed = NULL,
                                                  heter_groups = NULL,
                                                  msg = "\n==================================================\n") {
  if (!isTRUE(cross_validation) ||
      is.null(cross_validation_meth) ||
      is.null(GS_model_cv) ||
      is.null(engine) ||
      is.null(heter_groups)) {
    return(invisible(TRUE))
  }

  cv_token <- normalize_cv_token(as.character(cross_validation_meth))
  if (!any(cv_token %in% c("cv0", "repeated_cv0")) ||
      !any(as.character(GS_model_cv) %in% "GBLUP") ||
      !identical(tolower(as.character(engine)[[1L]]), "asreml")) {
    return(invisible(TRUE))
  }

  fixed_terms <- tryCatch(
    asreml_rhs_terms(fixed, labels = c("FIXED", "fixed_term_model", "NULL")),
    error = function(e) {
      if (inherits(fixed, "formula")) unique(all.vars(fixed)) else character()
    }
  )
  if (any(as.character(heter_groups) %in% fixed_terms)) {
    stop(
      paste(
        msg,
        "ASReml-R GBLUP CV0 with the environment term in fixed effects is not currently supported.",
        "In CV0, whole environments are held out; their fixed environment means are not estimable,",
        "and the ASReml-R predictor returns no finite held-out predictions.",
        "Use ASReml-R GBLUP with CV1/CV2, remove the held-out environment fixed effect only if your model is statistically justified,",
        "or use a supported CV0 MET model such as GBLUP_BRR, RKHS, GP, or MET ML/DL."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

gp_run_best_model_task <- function(ctx) {
  response <- ctx$response
  GS_model <- ctx$GS_model
  trait_pheno_res <- gp_prepare_pheno_for_trait(
    pheno_clean = ctx$pheno_clean,
    response = response,
    gen_name = ctx$gen_name
  )
  trait_pheno_data <- trait_pheno_res[["pheno_data"]]
  trait_test_set <- trait_pheno_res[["test_set"]]
  trait_test_set_source <- trait_pheno_res[["test_set_source"]]
  trait_ml_dat_res <- gp_prepare_ml_data_for_trait(
    ml_dat_res = ctx$ml_dat_res,
    pheno_clean = ctx$pheno_clean,
    response = response,
    gen_name = ctx$gen_name
  )
  trait_ml_dat_res <- gp_feature_apply_explicit_to_ml_dat_res(
    ml_dat_res = trait_ml_dat_res,
    feature_selected = ctx$feature_selected,
    trait = response
  )
  explicit_feature_meta <- NULL
  if (!is.null(ctx$feature_selected)) {
    source_matrices <- list(
      geno_data = ctx$geno_model_ready,
      omic1_data = ctx$omic1_model_ready,
      omic2_data = ctx$omic2_model_ready,
      omic3_data = ctx$omic3_model_ready
    )
    selected_by_source <- gp_feature_explicit_selected_by_source(
      feature_selected = ctx$feature_selected,
      trait = response,
      source_matrices = source_matrices,
      context = "feature_selected"
    )
    if (length(selected_by_source)) {
      selected_sources <- gp_feature_subset_source_matrices(
        source_matrices = source_matrices,
        selected_by_source = selected_by_source,
        context = paste0("Trait '", response, "' feature_selected")
      )
      ctx$geno_model_ready <- selected_sources[["geno_data"]]
      ctx$omic1_model_ready <- selected_sources[["omic1_data"]]
      ctx$omic2_model_ready <- selected_sources[["omic2_data"]]
      ctx$omic3_model_ready <- selected_sources[["omic3_data"]]
      precomputed_kernel_supplied <- any(vapply(
        list(ctx$gmatrix, ctx$gkernel, ctx$omic1_kernel, ctx$omic2_kernel, ctx$omic3_kernel, ctx$kernel_list),
        Negate(is.null),
        logical(1L)
      ))
      if (gp_feature_model_needs_kernel(GS_model, met_ml_dl = ctx$met_ml_dl) &&
          !isTRUE(precomputed_kernel_supplied)) {
        selected_kernel_bank <- gp_feature_rebuild_true_kernel_bank(
          model = GS_model,
          trait = response,
          selected_by_source = selected_by_source,
          source_matrices = source_matrices,
          gmatrix_method = ctx$gmatrix_method,
          kernel_method = ctx$kernel_method,
          ploidy = ctx$ploidy %||% "auto",
          scaling = ctx$scaling %||% TRUE,
          centering = ctx$centering %||% FALSE,
          met_ml_dl = ctx$met_ml_dl,
          precomputed_kernel_supplied = FALSE
        )
        if (!is.null(selected_kernel_bank)) {
          ctx$gmatrix_kernel_model_ready_list <- selected_kernel_bank
          ctx$gmatrix_model_ready <- selected_kernel_bank[["gmatrix_model_ready"]] %||% NULL
          ctx$omic1_kernel_model_ready <- selected_kernel_bank[["omic1_kernel_model_ready"]] %||% NULL
          ctx$omic2_kernel_model_ready <- selected_kernel_bank[["omic2_kernel_model_ready"]] %||% NULL
          ctx$omic3_kernel_model_ready <- selected_kernel_bank[["omic3_kernel_model_ready"]] %||% NULL
        }
      }
    }
    metadata_sources <- selected_by_source
    if (!length(metadata_sources)) {
      selected_names <- trait_ml_dat_res[["feature_selected"]] %||% character()
      if (length(selected_names)) {
        metadata_sources <- list(merged = unique(selected_names))
      } else {
        metadata_sources <- gp_feature_declared_by_source(ctx$feature_selected, trait = response)
      }
    }
    explicit_feature_meta <- gp_feature_explicit_metadata(
      selected_by_source = metadata_sources,
      trait = response,
      prediction_model = as.character(GS_model),
      n_train = sum(!is.na(trait_pheno_data[[response]]))
    )
  }
  trait_ml_dat_res <- gp_feature_apply_to_ml_dat_res(
    ml_dat_res = trait_ml_dat_res,
    feature_score_metadata = ctx$feature_score_metadata,
    trait = response,
    k = ctx$feature_k
  )
  if (!is.null(ctx$feature_score_metadata) && !is.null(ctx$feature_k)) {
    ctx$geno_model_ready <- gp_feature_select_matrix_by_metadata(
      ctx$geno_model_ready, ctx$feature_score_metadata, response, ctx$feature_k,
      source = "geno_data"
    )
    ctx$omic1_model_ready <- gp_feature_select_matrix_by_metadata(
      ctx$omic1_model_ready, ctx$feature_score_metadata, response, ctx$feature_k,
      source = "omic1_data"
    )
    ctx$omic2_model_ready <- gp_feature_select_matrix_by_metadata(
      ctx$omic2_model_ready, ctx$feature_score_metadata, response, ctx$feature_k,
      source = "omic2_data"
    )
    ctx$omic3_model_ready <- gp_feature_select_matrix_by_metadata(
      ctx$omic3_model_ready, ctx$feature_score_metadata, response, ctx$feature_k,
      source = "omic3_data"
    )
    selected_by_source <- gp_feature_selected_by_source(
      ctx$feature_score_metadata,
      trait = response,
      k = ctx$feature_k
    )
    selected_kernel_bank <- gp_feature_rebuild_true_kernel_bank(
      model = GS_model,
      trait = response,
      selected_by_source = selected_by_source,
      source_matrices = list(
        geno_data = ctx$geno_model_ready,
        omic1_data = ctx$omic1_model_ready,
        omic2_data = ctx$omic2_model_ready,
        omic3_data = ctx$omic3_model_ready
      ),
      gmatrix_method = ctx$gmatrix_method,
      kernel_method = ctx$kernel_method,
      ploidy = ctx$ploidy %||% "auto",
      scaling = ctx$scaling %||% TRUE,
      centering = ctx$centering %||% FALSE,
      met_ml_dl = ctx$met_ml_dl,
      precomputed_kernel_supplied = any(vapply(
        list(ctx$gmatrix, ctx$gkernel, ctx$omic1_kernel, ctx$omic2_kernel, ctx$omic3_kernel, ctx$kernel_list),
        Negate(is.null),
        logical(1L)
      ))
    )
    if (!is.null(selected_kernel_bank)) {
      ctx$gmatrix_kernel_model_ready_list <- selected_kernel_bank
      ctx$gmatrix_model_ready <- selected_kernel_bank[["gmatrix_model_ready"]] %||% NULL
      ctx$omic1_kernel_model_ready <- selected_kernel_bank[["omic1_kernel_model_ready"]] %||% NULL
      ctx$omic2_kernel_model_ready <- selected_kernel_bank[["omic2_kernel_model_ready"]] %||% NULL
      ctx$omic3_kernel_model_ready <- selected_kernel_bank[["omic3_kernel_model_ready"]] %||% NULL
    }
  }
  res_model_output <- NULL
  res_summary_stat <- NULL
  is_multi_env_bayes <- gp_is_multi_environment_pheno(
    pheno_data = trait_pheno_data,
    gen_name = ctx$gen_name,
    heter_groups = ctx$heter_groups,
    response = response,
    response_family = ctx$response_family
  )

  if (isTRUE(is_multi_env_bayes) && !is.null(GS_model) && GS_model %in% ctx$bayes_valid_models) {
    warning(
      paste(
        "Bayesian marker-regression models",
        paste(ctx$bayes_valid_models, collapse = ", "),
        "are currently single-environment only in this package.",
        "Use RKHS or GBLUP_BRR for multi-environment Bayesian prediction."
      ),
      call. = FALSE
    )
  }

  if (is.null(res_model_output) &&
      length(GS_model) == 1L &&
      !is.na(GS_model) &&
      as.character(GS_model) %in% gp_lowrank_supported_models()) {
    gp_true_prediction_run <- !isTRUE(ctx$cross_validation) ||
      (isTRUE(ctx$cross_validation) && !isTRUE(ctx$cv_evaluation_only))
    requested_gp_output_level <- ctx$gp_output_level %||% "predict_only"
    gp_final_output_level <- if (isTRUE(gp_true_prediction_run) &&
                                 identical(requested_gp_output_level, "predict_only")) {
      "full_vc"
    } else {
      requested_gp_output_level
    }
    gp_final_prediction_output <- if (isTRUE(gp_true_prediction_run)) {
      "all"
    } else {
      ctx$gp_prediction_output %||% "all"
    }
    lowrank_fit <- tryCatch({
      gp_backend_gaussian_model(
        model_name = GS_model,
        pheno_data = trait_pheno_data,
        response = response,
        gen_name = ctx$gen_name,
        gmatrix = ctx$gmatrix_model_ready,
        omic1_kernel = ctx$omic1_kernel_model_ready,
        omic2_kernel = ctx$omic2_kernel_model_ready,
        omic3_kernel = ctx$omic3_kernel_model_ready,
        kernel_list = ctx$gmatrix_kernel_model_ready_list,
        heter_groups = ctx$heter_groups,
        test_set = trait_test_set,
        test_set_source = trait_test_set_source,
        lowrank_eps_trace = ctx$lowrank_eps_trace %||% 1e-6,
        lowrank_max_rank = ctx$lowrank_max_rank,
        lowrank_jitter = ctx$lowrank_jitter,
        lowrank_noise_grid = ctx$lowrank_noise_grid,
        lowrank_kernel_weights = ctx$lowrank_kernel_weights,
        system_database = ctx$system_database,
        gp_backend = ctx$gp_backend %||% "auto",
        gp_output_level = gp_final_output_level,
        gp_varcomp_mode = ctx$gp_varcomp_mode %||% "reml",
        gp_fa_rank = ctx$gp_fa_rank %||% 1L,
        gp_prediction_output = gp_final_prediction_output,
        gp_factor_cache = ctx$gp_factor_cache %||% NULL,
        gp_iters = ctx$gp_iters %||% NULL,
        gp_lr = ctx$gp_lr %||% NULL,
        env_similarity = ctx$env_similarity %||% NULL,
        env_ids = ctx$env_ids %||% NULL,
        env_covariates = ctx$env_covariates %||% NULL,
        reaction_norm_feature_qc = ctx$reaction_norm_feature_qc %||% TRUE,
        kenv_kernel = ctx$kenv_kernel %||% "matern32",
        kenv_bandwidth = ctx$kenv_bandwidth %||% 1.0,
        kenv_kernel_kwargs = ctx$kenv_kernel_kwargs %||% NULL,
        gp_engine = ctx$gp_engine %||% "auto",
        gp_learn_scales = ctx$gp_learn_scales %||% NULL,
        observation_weights = ctx$weights
      )
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(lowrank_fit)) {
      res_model_output <- lowrank_fit
      GS_model_public <- gp_public_model_label(GS_model)[1]
      res_model_output[["model_parameters"]] <- gp_relabel_public_model_parameters(
        res_model_output[["model_parameters"]],
        GS_model
      )
      res_summary_stat <- tryCatch({
        summary_statistics_AI(
          predicted_object = res_model_output[["predicted_values"]],
          pheno_object = trait_pheno_data[!is.na(trait_pheno_data[[response]]), , drop = FALSE],
          response = response,
          test_set = trait_test_set,
          geno_omic_object = res_model_output[["lowrank_feature_matrix"]],
          model_parameters = res_model_output[["model_parameters"]],
          eval_metrics = ctx$eval_metrics,
          response_family = ctx$response_family,
          GS_model = GS_model_public
        )
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        NULL
      })
    }
  }

  if (length(trait_pheno_data[, ctx$gen_name]) ==
      length(unique(trait_pheno_data[, ctx$gen_name]))) {
    if ((GS_model %in% ctx$bayes_valid_models && is.null(ctx$rand_term_model_bayesian)) ||
        (is.null(GS_model) && any(ctx$rand_term_model_bayesian %in% ctx$bayes_valid_models)) ||
        (!is.null(GS_model) && any(ctx$rand_term_model_bayesian %in% ctx$bayes_valid_models))) {
      bayes_A_B_C_BRR_mod_process <- tryCatch({
        bayes_finalize_A_B_C_BL_BRR(
          fixed = ctx$fixed,
          random = ctx$random,
          GS_model = GS_model,
          response = response,
          weights = ctx$weights,
          fixed_term_model_bayesian = ctx$fixed_term_model_bayesian,
          rand_term_model_bayesian = ctx$rand_term_model_bayesian,
          pheno_data = trait_pheno_data,
          geno_data = ctx$geno_model_ready,
          omic1_data = ctx$omic1_model_ready,
          omic2_data = ctx$omic2_model_ready,
          omic3_data = ctx$omic3_model_ready,
          gen_name = ctx$gen_name,
          nIter = ctx$nIter,
          burnIn = ctx$burnIn,
          thin = ctx$thin,
          omics_data_label = ctx$omics_data_label,
          scaling = ctx$scaling,
          CI_width_thresholds = ctx$CI_width_thresholds,
          confidence_level = ctx$confidence_level,
          high_reliability_thres = ctx$high_reliability_thres,
          low_reliability_thres = ctx$low_reliability_thres,
          n_components = ctx$n_components,
          threshold = ctx$threshold,
          target = "test_set",
          interval_width_high_threshold = ctx$interval_width_high_threshold,
          interval_width_low_threshold = ctx$interval_width_low_threshold,
          interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
          response_family = ctx$response_family
        )
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        NULL
      })

      if (!is.null(bayes_A_B_C_BRR_mod_process)) {
        res_model_output <- bayes_A_B_C_BRR_mod_process
        bayes_summary_stat_process <- tryCatch({
          summary_statistics_bayes(
            mod = res_model_output[["bayes_model"]],
            eval_metrics = ctx$eval_metrics,
            model_result = res_model_output[["bayes_result"]],
            GS_model = GS_model,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            CI_width_thresholds = ctx$CI_width_thresholds,
            confidence_level = ctx$confidence_level,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
            ,
            response_family = ctx$response_family
          )
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          NULL
        })

        if (!is.null(bayes_summary_stat_process)) {
          res_summary_stat <- bayes_summary_stat_process
          if ("diagnostic_tst_plot" %in% names(res_summary_stat)) {
            res_model_output[["diagnostic_plots"]] <- res_summary_stat[["diagnostic_tst_plot"]]
            res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
          }
        } else {
          cat(sprintf("Bayesian %s summary statistics failed.\n", GS_model))
        }
      } else {
        cat(sprintf("Bayesian %s model failed.\n", GS_model))
      }
    }
  }

  if (is.null(res_model_output) && (
    (GS_model %in% c("RKHS", "GBLUP_BRR", "GBLUP") && is.null(ctx$rand_term_model_bayesian)) ||
    (is.null(GS_model) && any(ctx$rand_term_model_bayesian %in% ctx$bayes_gblup_valid_models)) ||
    (!is.null(GS_model) && any(ctx$rand_term_model_bayesian %in% ctx$bayes_gblup_valid_models))
  )) {
    if (GS_model %in% c("BRR", "RKHS", "GBLUP_BRR")) {
      # Phase 3.21: resolve the kernel-Bayes-specific heter_resid override.
      bayes_kernel_heter_resid_resolved <- if (!is.null(ctx$bayes_kernel_heter_resid)) {
        isTRUE(ctx$bayes_kernel_heter_resid)
      } else {
        ctx$heter_resid
      }
      bayes_RKHS_GBLUP_BRR_mod_process <- tryCatch({
        bayes_finalize_RKHS_GBLUPBRR(
          fixed = ctx$fixed,
          random = ctx$random,
          GS_model = if (GS_model == "GBLUP_BRR") "BRR" else GS_model,
          response = response,
          weights = ctx$weights,
          fixed_term_model_bayesian = ctx$fixed_term_model_bayesian,
          rand_term_model_bayesian = ctx$rand_term_model_bayesian,
          pheno_data = trait_pheno_data,
          gmatrix = ctx$gmatrix_model_ready,
          omic1_kernel = ctx$omic1_kernel_model_ready,
          omic2_kernel = ctx$omic2_kernel_model_ready,
          omic3_kernel = ctx$omic3_kernel_model_ready,
          kernel_list = ctx$gmatrix_kernel_model_ready_list,
          gen_name = ctx$gen_name,
          nIter = ctx$nIter,
          burnIn = ctx$burnIn,
          thin = ctx$thin,
          heter_groups = ctx$heter_groups,
          heter_resid = bayes_kernel_heter_resid_resolved,
          omics_kernel_label = ctx$omics_kernel_label,
          CI_width_thresholds = ctx$CI_width_thresholds,
          confidence_level = ctx$confidence_level,
          high_reliability_thres = ctx$high_reliability_thres,
          low_reliability_thres = ctx$low_reliability_thres,
          n_components = ctx$n_components,
          threshold = ctx$threshold,
          target = "test_set",
          cross_validation = FALSE,
          interval_width_low_threshold = ctx$interval_width_low_threshold,
          interval_width_high_threshold = ctx$interval_width_high_threshold,
          interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
          response_family = ctx$response_family
        )
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        NULL
      })

      if (!is.null(bayes_RKHS_GBLUP_BRR_mod_process)) {
        res_model_output <- bayes_RKHS_GBLUP_BRR_mod_process
        bayes_GBLUP_summary_stat_process <- tryCatch({
          summary_statistics_bayes(
            mod = res_model_output[["bayes_model"]],
            eval_metrics = ctx$eval_metrics,
            model_result = res_model_output[["bayes_result"]],
            GS_model = GS_model,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            CI_width_thresholds = ctx$CI_width_thresholds,
            confidence_level = ctx$confidence_level,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
            ,
            response_family = ctx$response_family
          )
        }, error = function(e) {
          cat("Error:", conditionMessage(e), "\n")
          NULL
        })

        if (!is.null(bayes_GBLUP_summary_stat_process)) {
          res_summary_stat <- bayes_GBLUP_summary_stat_process
          if ("diagnostic_tst_plot" %in% names(res_summary_stat)) {
            res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
          }
        } else {
          cat(sprintf("Bayesian %s summary statistics failed.\n", GS_model))
        }
      }
    } else if (isTRUE(length(GS_model) == 1L && !is.na(GS_model) &&
                      identical(as.character(GS_model), "GBLUP") &&
                      length(ctx$engine) == 1L && !is.na(ctx$engine) &&
                      identical(as.character(ctx$engine), "asreml"))) {
      result_asreml_mod <- tryCatch({
        asreml_utilis_new(
          fixed = ctx$fixed,
          random = ctx$random,
          cova = ctx$cova,
          GS_model = GS_model,
          response = response,
          pheno_data = trait_pheno_data,
          gmatrix = ctx$gmatrix_model_ready,
          omic1_kernel = ctx$omic1_kernel_model_ready,
          omic2_kernel = ctx$omic2_kernel_model_ready,
          omic3_kernel = ctx$omic3_kernel_model_ready,
          kernel_list = ctx$gmatrix_kernel_model_ready_list,
          gen_name = ctx$gen_name,
          heter_groups = ctx$heter_groups,
          heter_resid = ctx$heter_resid,
          var_cov_str = ctx[["var_cov_str"]],
          weights = ctx$weights,
          pworkspace = ctx$pworkspace,
          workspace = ctx$workspace,
          maxit = ctx$maxit,
          inverse = ctx$inverse,
          epsilon = ctx$epsilon,
          cross_validation = FALSE,
          engine = ctx$engine
        )
      }, error = function(e) {
        stop("ASReml GBLUP fit failed: ", conditionMessage(e), call. = FALSE)
      })

      if (!is.null(result_asreml_mod)) {
        mod <- result_asreml_mod
        vc <- asreml_varcomp_table(mod$model)
        res_comp_checkk <- tryCatch({
          VAR_check_Pos <- which(vc$bound == "?" | vc$bound == "S")
          if (length(VAR_check_Pos) > 0) {
            unstable_components <- paste(
              "variance component for",
              paste(rownames(vc)[VAR_check_Pos], collapse = " and"),
              "is unstable, refix the model"
            )
            stop(unstable_components, call. = FALSE)
          }
          vc
        }, error = function(e) {
          stop("ASReml variance-component validation failed: ",
               conditionMessage(e), call. = FALSE)
        })

        if (!is.null(res_comp_checkk)) {
          asreml_mod_output_process <- tryCatch({
            asreml_mod_output_new(
              mod_asreml = mod,
              pheno_data = trait_pheno_data,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              heter_groups = ctx$heter_groups,
              gen_name = ctx$gen_name,
              response = response,
                var_cov_str = ctx[["var_cov_str"]],
              heter_resid = ctx$heter_resid,
              pworkspace = ctx$pworkspace,
              workspace = ctx$workspace,
              maxit = ctx$maxit
            )
          }, error = function(e) {
            stop("ASReml output processing failed: ",
                 conditionMessage(e), call. = FALSE)
          })

          if (!is.null(asreml_mod_output_process)) {
            res_model_output <- asreml_mod_output_process
            asreml_summary_stat_process <- tryCatch({
              summary_statistics_asreml(
                mod = res_model_output[["Asreml_model"]],
                response = response,
                pheno_data = trait_pheno_data,
                heter_groups = ctx$heter_groups,
                GID_names = res_model_output[["Predicted_value"]][ctx$gen_name],
                predicted_value = if ("Predicted_value" %in% colnames(res_model_output[["Predicted_value"]]))
                  res_model_output[["Predicted_value"]]["Predicted_value"] else res_model_output[["Predicted_value"]]["BLUP"],
                standard_errors = res_model_output[["Predicted_value"]]["Standard_error"],
                prediction_error_var = res_model_output[["Predicted_value"]]["Prediction_error_variance"],
                genetic_var = {
                  asreml_vc <- res_model_output[["Variance_components"]]
                  genetic_row <- if (
                    "mean_genetic_variance_across_environments" %in% rownames(asreml_vc)
                  ) {
                    "mean_genetic_variance_across_environments"
                  } else {
                    "genetic_variance"
                  }
                  sum(asreml_vc[genetic_row, "Components"])
                },
                pred_heter_groups = NULL,
                variance_components = res_model_output[["Variance_components"]],
                eval_metrics = ctx$eval_metrics,
                gen_name = ctx$gen_name,
                CI_width_thresholds = ctx$CI_width_thresholds,
                confidence_level = ctx$confidence_level,
                high_reliability_thres = ctx$high_reliability_thres,
                low_reliability_thres = ctx$low_reliability_thres,
                system_database = ctx$system_database,
                response_family = ctx$response_family
              )
            }, error = function(e) {
              stop("ASReml summary processing failed: ",
                   conditionMessage(e), call. = FALSE)
            })

            if (!is.null(asreml_summary_stat_process)) {
              res_summary_stat <- asreml_summary_stat_process
              if ("diagnostic_tst_plot" %in% names(res_summary_stat)) {
                res_model_output[["diagnostic_plots"]] <- res_summary_stat[["diagnostic_tst_plot"]]
                res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
              }
            } else {
              cat(paste(ctx$msg, "Error processing summary statistics for asreml result.\n"))
            }
          } else {
            cat(paste(ctx$msg, "Error processing the output of asreml result.\n"))
          }
        }
      }
    }
  }

  if (is.null(res_model_output) && GS_model %in% ctx$AI_valid_models) {
    if (isTRUE(is_multi_env_bayes) && isTRUE(ctx$met_ml_dl) && GS_model %in% c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "mlp", "ft_transformer", "saint", "tabnet", "moe")) {
      res_model_output <- tryCatch({
        switch(
          GS_model,
          CatBoost = AI_CatBoost_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            para_tunning = ctx$para_tunning,
            catboost_iterations = ctx$catboost_iterations,
            catboost_depth = ctx$catboost_depth,
            catboost_learning_rate = ctx$catboost_learning_rate,
            catboost_l2_leaf_reg = ctx$catboost_l2_leaf_reg,
            catboost_thread_count = ctx$catboost_thread_count,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds
          ),
          LightGBM = AI_LightGBM_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            para_tunning = ctx$para_tunning,
            lightgbm_paras_tunning = ctx$lightgbm_paras_tunning,
            lightgbm_nrounds = ctx$lightgbm_nrounds,
            lightgbm_learning_rate = ctx$lightgbm_learning_rate,
            lightgbm_num_leaves = ctx$lightgbm_num_leaves,
            lightgbm_feature_fraction = ctx$lightgbm_feature_fraction,
            lightgbm_bagging_fraction = ctx$lightgbm_bagging_fraction,
            lightgbm_min_data_in_leaf = ctx$lightgbm_min_data_in_leaf,
            lightgbm_lambda_l1 = ctx$lightgbm_lambda_l1,
            lightgbm_lambda_l2 = ctx$lightgbm_lambda_l2,
            lightgbm_nthread = ctx$lightgbm_nthread,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds
          ),
          Xgboost = AI_Xgb_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            para_tunning = ctx$para_tunning,
            xgb_paras_tunning = ctx$xgb_paras_tunning,
            learning_rate = ctx$learning_rate,
            max_depth = ctx$max_depth,
            subsample = ctx$subsample,
            xgb_booster = ctx$xgb_booster,
            xgb_alpha = ctx$xgb_alpha,
            xgb_lambda = ctx$xgb_lambda,
            xgb_gamma = ctx$xgb_gamma,
            min_child_weight = ctx$min_child_weight,
            iteration = ctx$iteration,
            xgb_nthread = ctx$xgb_nthread,
            colsample_bytree = ctx$colsample_bytree,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds
          ),
          RandomForest = AI_randomForest_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            para_tunning = ctx$para_tunning,
            rf_paras_tunning = ctx$rf_paras_tunning,
            ntree = ctx$ntree,
            mtry = ctx$mtry,
            maxnodes = ctx$maxnodes,
            nodesize = ctx$nodesize,
            rf_n_jobs = ctx$rf_n_jobs,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds
          ),
          mlp = AI_mlp_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            compile_model = ctx$compile_model,
            deterministic = ctx$deterministic,
            random_seed = ctx$random_seed,
            device = ctx$device,
            use_amp = ctx$use_amp,
            batch_norm = ctx$batch_norm,
            validation_split = ctx$validation_split,
            epochs = ctx$epochs,
            batch_size = ctx$batch_size,
            mlp_neurons_per_layer = ctx$mlp_neurons_per_layer,
            mlp_learning_rate = ctx$mlp_learning_rate,
            dropout = ctx$dropout,
            l2_weight_decay = ctx$l2_weight_decay,
            final_attention = ctx$final_attention,
            attention_across_multiple_layers = ctx$attention_across_multiple_layers,
            optimizer_name = ctx$optimizer_name,
            max_grad_norm = ctx$max_grad_norm,
            heteroscedastic = ctx$heteroscedastic,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds,
            dl_internal_calibration = ctx$dl_internal_calibration,
            CI_width_thresholds = ctx$CI_width_thresholds,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
          ),
          ft_transformer = AI_ft_transformer_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            compile_model = ctx$compile_model,
            deterministic = ctx$deterministic,
            random_seed = ctx$random_seed,
            device = ctx$device,
            use_amp = ctx$use_amp,
            batch_norm = ctx$batch_norm,
            validation_split = ctx$validation_split,
            epochs = ctx$epochs,
            batch_size = ctx$batch_size,
            ft_d_model = ctx$ft_d_model,
            ft_heads = ctx$ft_heads,
            ft_layers = ctx$ft_layers,
            ft_ff_mult = ctx$ft_ff_mult,
            ft_dropout = ctx$ft_dropout,
            ft_token_dropout = ctx$ft_token_dropout,
            ft_use_cls = ctx$ft_use_cls,
            dropout = ctx$dropout,
            l2_weight_decay = ctx$l2_weight_decay,
            optimizer_name = ctx$optimizer_name,
            max_grad_norm = ctx$max_grad_norm,
            heteroscedastic = ctx$heteroscedastic,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds,
            dl_internal_calibration = ctx$dl_internal_calibration,
            CI_width_thresholds = ctx$CI_width_thresholds,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
          ),
          saint = AI_saint_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            compile_model = ctx$compile_model,
            deterministic = ctx$deterministic,
            random_seed = ctx$random_seed,
            device = ctx$device,
            use_amp = ctx$use_amp,
            batch_norm = ctx$batch_norm,
            validation_split = ctx$validation_split,
            epochs = ctx$epochs,
            batch_size = ctx$batch_size,
            saint_d_model = ctx$saint_d_model,
            saint_heads = ctx$saint_heads,
            saint_layers = ctx$saint_layers,
            saint_ff_mult = ctx$saint_ff_mult,
            saint_dropout = ctx$saint_dropout,
            saint_token_dropout = ctx$saint_token_dropout,
            saint_use_cls = ctx$saint_use_cls,
            dropout = ctx$dropout,
            l2_weight_decay = ctx$l2_weight_decay,
            optimizer_name = ctx$optimizer_name,
            max_grad_norm = ctx$max_grad_norm,
            heteroscedastic = ctx$heteroscedastic,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds,
            dl_internal_calibration = ctx$dl_internal_calibration,
            CI_width_thresholds = ctx$CI_width_thresholds,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
          ),
          tabnet = AI_tabnet_MET(
            pheno_object = trait_pheno_data,
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            response_family = ctx$response_family,
            gmatrix = ctx$gmatrix_model_ready,
            omic1_kernel = ctx$omic1_kernel_model_ready,
            omic2_kernel = ctx$omic2_kernel_model_ready,
            omic3_kernel = ctx$omic3_kernel_model_ready,
            kernel_list = ctx$gmatrix_kernel_model_ready_list,
            omics_kernel_label = ctx$omics_kernel_label,
            met_kernel_var_explained = ctx$met_kernel_var_explained,
            met_kernel_min_ev = ctx$met_kernel_min_ev,
            met_kernel_max_pcs = ctx$met_kernel_max_pcs,
            compile_model = ctx$compile_model,
            deterministic = ctx$deterministic,
            random_seed = ctx$random_seed,
            device = ctx$device,
            use_amp = ctx$use_amp,
            batch_norm = ctx$batch_norm,
            validation_split = ctx$validation_split,
            epochs = ctx$epochs,
            batch_size = ctx$batch_size,
            tabnet_steps = ctx$tabnet_steps,
            tabnet_feature_dim = ctx$tabnet_feature_dim,
            tabnet_output_dim = ctx$tabnet_output_dim,
            tabnet_gamma = ctx$tabnet_gamma,
            tabnet_lambda_sparse = ctx$tabnet_lambda_sparse,
            dropout = ctx$dropout,
            l2_weight_decay = ctx$l2_weight_decay,
            optimizer_name = ctx$optimizer_name,
            max_grad_norm = ctx$max_grad_norm,
            heteroscedastic = ctx$heteroscedastic,
            n_bootstrap = ctx$n_bootstrap,
            internal_cv_nfolds = ctx$internal_cv_nfolds,
            dl_internal_calibration = ctx$dl_internal_calibration,
            CI_width_thresholds = ctx$CI_width_thresholds,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            system_database = ctx$system_database
          )
        )
      }, error = function(e) {
        cat(paste(ctx$msg, "Error in", GS_model, "MET model:", conditionMessage(e), "\n"))
        NULL
      })
    }

    run_xgboost <- function() tryCatch({
      AI_Xgb(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        para_tunning = ctx$para_tunning,
        xgb_paras_tunning = ctx$xgb_paras_tunning,
        resample_method_tune = ctx$resample_method_tune,
        number_of_fold_tune = ctx$number_of_fold_tune,
        learning_rate = ctx$learning_rate,
        xgb_gamma = ctx$xgb_gamma,
        xgb_lambda = ctx$xgb_lambda,
        xgb_alpha = ctx$xgb_alpha,
        max_depth = ctx$max_depth,
        subsample = ctx$subsample,
        xgb_booster = ctx$xgb_booster,
        colsample_bytree = ctx$colsample_bytree,
        alpha = ctx$xgb_alpha,
        lambda = ctx$xgb_lambda,
        iteration = ctx$iteration,
        xgb_rate_drop = ctx$xgb_rate_drop,
        xgb_skip_drop = ctx$xgb_skip_drop,
        xgb_objective = ctx$xgb_objective,
        xgb_sample_type = ctx$xgb_sample_type,
        xgb_normalize_type = ctx$xgb_normalize_type,
        xgb_nthread = ctx$xgb_nthread,
        early_stop_for_iteration_xgb = ctx$early_stop_for_iteration_xgb,
        N_feature_impo = ctx$N_feature_impo,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in Xgboost model:", conditionMessage(e), "\n"))
      NULL
    })

    run_random_forest <- function() tryCatch({
      AI_randomForest(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        para_tunning = ctx$para_tunning,
        rf_paras_tunning = ctx$rf_paras_tunning,
        ntree = ctx$ntree,
        mtry = ctx$mtry,
        maxnodes = ctx$maxnodes,
        nodesize = ctx$nodesize,
        rf_n_jobs = ctx$rf_n_jobs,
        importance = ctx$importance,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in RandomForest model:", conditionMessage(e), "\n"))
      NULL
    })

    run_catboost <- function() tryCatch({
      AI_CatBoost(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        para_tunning = ctx$para_tunning,
        catboost_iterations = ctx$catboost_iterations,
        catboost_depth = ctx$catboost_depth,
        catboost_learning_rate = ctx$catboost_learning_rate,
        catboost_l2_leaf_reg = ctx$catboost_l2_leaf_reg,
        catboost_thread_count = ctx$catboost_thread_count,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in CatBoost model:", conditionMessage(e), "\n"))
      NULL
    })

    run_lightgbm <- function() tryCatch({
      AI_LightGBM(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        para_tunning = ctx$para_tunning,
        lightgbm_nrounds = ctx$lightgbm_nrounds,
        lightgbm_learning_rate = ctx$lightgbm_learning_rate,
        lightgbm_num_leaves = ctx$lightgbm_num_leaves,
        lightgbm_feature_fraction = ctx$lightgbm_feature_fraction,
        lightgbm_bagging_fraction = ctx$lightgbm_bagging_fraction,
        lightgbm_min_data_in_leaf = ctx$lightgbm_min_data_in_leaf,
        lightgbm_lambda_l1 = ctx$lightgbm_lambda_l1,
        lightgbm_lambda_l2 = ctx$lightgbm_lambda_l2,
        lightgbm_nthread = ctx$lightgbm_nthread,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in LightGBM model:", conditionMessage(e), "\n"))
      NULL
    })

    run_pls <- function() tryCatch({
      AI_pls(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        para_tunning = ctx$para_tunning,
        ncomp = ctx$ncomp,
        pls_paras_tunning = ctx$pls_paras_tunning,
        resample_method_tune = ctx$resample_method_tune,
        N_feature_impo = ctx$N_feature_impo,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in PartialLeastSquare model:", conditionMessage(e), "\n"))
      NULL
    })

    run_svm <- function() tryCatch({
      AI_svm(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        para_tunning = ctx$para_tunning,
        svm_paras_tunning = ctx$svm_paras_tunning,
        svm_type = ctx$svm_type,
        svm_kernel = ctx$svm_kernel,
        sigma_value = ctx$sigma_value,
        C_value = ctx$C_value,
        degree_value = ctx$degree_value,
        scale_value = ctx$scale_value,
        offset_value = ctx$offset_value,
        gamma_value = ctx$gamma_value,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in SupportVectorMachine model:", conditionMessage(e), "\n"))
      NULL
    })

    run_knn <- function() tryCatch({
      AI_knn(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        response_family = ctx$response_family,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        message = ctx$message,
        gen_name = ctx$gen_name,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        para_tunning = ctx$para_tunning,
        knn_paras_tunning = ctx$knn_paras_tunning,
        k = ctx$k,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in K-NearestNeighbors model:", conditionMessage(e), "\n"))
      NULL
    })

    run_lasso <- function() tryCatch({
      AI_RidgeRegression_Lasso(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        gen_name = ctx$gen_name,
        para_tunning = ctx$para_tunning,
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        lasso_paras_tunning = ctx$lasso_paras_tunning,
        message = ctx$message,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        GS_model = GS_model,
        lambda_rr = ctx$lambda_rr,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in Lasso model:", conditionMessage(e), "\n"))
      NULL
    })

    run_ridge <- function() tryCatch({
      AI_RidgeRegression_Lasso(
        pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
        response = response,
        geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
        geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
        gen_name = ctx$gen_name,
        para_tunning = ctx$para_tunning,
        AI_cv_nfolds = ctx$AI_cv_nfolds,
        lasso_paras_tunning = ctx$rr_paras_tunning,
        message = ctx$message,
        scaling = ctx$scaling,
        centering = ctx$centering,
        omic_count = trait_ml_dat_res[["omic_count"]],
        GS_model = GS_model,
        lambda_rr = ctx$lambda_rr,
        CI_width_thresholds = ctx$CI_width_thresholds,
        high_reliability_thres = ctx$high_reliability_thres,
        low_reliability_thres = ctx$low_reliability_thres,
        n_components = ctx$n_components,
        threshold = ctx$threshold,
        target = "test_set",
        interval_width_low_threshold = ctx$interval_width_low_threshold,
        interval_width_high_threshold = ctx$interval_width_high_threshold,
        interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
        n_bootstrap = ctx$n_bootstrap
      )
    }, error = function(e) {
      cat(paste(ctx$msg, "Error in Ridge Regression model:", conditionMessage(e), "\n"))
      NULL
    })

    if (GS_model %in% ctx$canonical_names) {
      res_model_output <- tryCatch({
        if (isTRUE(is_multi_env_bayes) && isTRUE(ctx$met_ml_dl) && GS_model %in% c("mlp", "ft_transformer", "saint", "tabnet", "moe")) {
          switch(
            GS_model,
            mlp = AI_mlp_MET(
              pheno_object = trait_pheno_data,
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              response_family = ctx$response_family,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              omics_kernel_label = ctx$omics_kernel_label,
              met_kernel_var_explained = ctx$met_kernel_var_explained,
              met_kernel_min_ev = ctx$met_kernel_min_ev,
              met_kernel_max_pcs = ctx$met_kernel_max_pcs,
              compile_model = ctx$compile_model,
              deterministic = ctx$deterministic,
              random_seed = ctx$random_seed,
              device = ctx$device,
              use_amp = ctx$use_amp,
              batch_norm = ctx$batch_norm,
              validation_split = ctx$validation_split,
              epochs = ctx$epochs,
              batch_size = ctx$batch_size,
              mlp_neurons_per_layer = ctx$mlp_neurons_per_layer,
              mlp_learning_rate = ctx$mlp_learning_rate,
              dropout = ctx$dropout,
              l2_weight_decay = ctx$l2_weight_decay,
              final_attention = ctx$final_attention,
              attention_across_multiple_layers = ctx$attention_across_multiple_layers,
              optimizer_name = ctx$optimizer_name,
              max_grad_norm = ctx$max_grad_norm,
              heteroscedastic = ctx$heteroscedastic,
              n_bootstrap = ctx$n_bootstrap,
              internal_cv_nfolds = ctx$internal_cv_nfolds,
              dl_internal_calibration = ctx$dl_internal_calibration,
              CI_width_thresholds = ctx$CI_width_thresholds,
              high_reliability_thres = ctx$high_reliability_thres,
              low_reliability_thres = ctx$low_reliability_thres,
              system_database = ctx$system_database
            ),
            ft_transformer = AI_ft_transformer_MET(
              pheno_object = trait_pheno_data,
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              response_family = ctx$response_family,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              omics_kernel_label = ctx$omics_kernel_label,
              met_kernel_var_explained = ctx$met_kernel_var_explained,
              met_kernel_min_ev = ctx$met_kernel_min_ev,
              met_kernel_max_pcs = ctx$met_kernel_max_pcs,
              compile_model = ctx$compile_model,
              deterministic = ctx$deterministic,
              random_seed = ctx$random_seed,
              device = ctx$device,
              use_amp = ctx$use_amp,
              batch_norm = ctx$batch_norm,
              validation_split = ctx$validation_split,
              epochs = ctx$epochs,
              batch_size = ctx$batch_size,
              ft_d_model = ctx$ft_d_model,
              ft_heads = ctx$ft_heads,
              ft_layers = ctx$ft_layers,
              ft_ff_mult = ctx$ft_ff_mult,
              ft_dropout = ctx$ft_dropout,
              ft_token_dropout = ctx$ft_token_dropout,
              ft_use_cls = ctx$ft_use_cls,
              dropout = ctx$dropout,
              l2_weight_decay = ctx$l2_weight_decay,
              optimizer_name = ctx$optimizer_name,
              max_grad_norm = ctx$max_grad_norm,
              heteroscedastic = ctx$heteroscedastic,
              n_bootstrap = ctx$n_bootstrap,
              internal_cv_nfolds = ctx$internal_cv_nfolds,
              dl_internal_calibration = ctx$dl_internal_calibration,
              CI_width_thresholds = ctx$CI_width_thresholds,
              high_reliability_thres = ctx$high_reliability_thres,
              low_reliability_thres = ctx$low_reliability_thres,
              system_database = ctx$system_database
            ),
            saint = AI_saint_MET(
              pheno_object = trait_pheno_data,
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              response_family = ctx$response_family,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              omics_kernel_label = ctx$omics_kernel_label,
              met_kernel_var_explained = ctx$met_kernel_var_explained,
              met_kernel_min_ev = ctx$met_kernel_min_ev,
              met_kernel_max_pcs = ctx$met_kernel_max_pcs,
              compile_model = ctx$compile_model,
              deterministic = ctx$deterministic,
              random_seed = ctx$random_seed,
              device = ctx$device,
              use_amp = ctx$use_amp,
              batch_norm = ctx$batch_norm,
              validation_split = ctx$validation_split,
              epochs = ctx$epochs,
              batch_size = ctx$batch_size,
              saint_d_model = ctx$saint_d_model,
              saint_heads = ctx$saint_heads,
              saint_layers = ctx$saint_layers,
              saint_ff_mult = ctx$saint_ff_mult,
              saint_dropout = ctx$saint_dropout,
              saint_token_dropout = ctx$saint_token_dropout,
              saint_use_cls = ctx$saint_use_cls,
              dropout = ctx$dropout,
              l2_weight_decay = ctx$l2_weight_decay,
              optimizer_name = ctx$optimizer_name,
              max_grad_norm = ctx$max_grad_norm,
              heteroscedastic = ctx$heteroscedastic,
              n_bootstrap = ctx$n_bootstrap,
              internal_cv_nfolds = ctx$internal_cv_nfolds,
              dl_internal_calibration = ctx$dl_internal_calibration,
              CI_width_thresholds = ctx$CI_width_thresholds,
              high_reliability_thres = ctx$high_reliability_thres,
              low_reliability_thres = ctx$low_reliability_thres,
              system_database = ctx$system_database
            ),
            tabnet = AI_tabnet_MET(
              pheno_object = trait_pheno_data,
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              response_family = ctx$response_family,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              omics_kernel_label = ctx$omics_kernel_label,
              met_kernel_var_explained = ctx$met_kernel_var_explained,
              met_kernel_min_ev = ctx$met_kernel_min_ev,
              met_kernel_max_pcs = ctx$met_kernel_max_pcs,
              compile_model = ctx$compile_model,
              deterministic = ctx$deterministic,
              random_seed = ctx$random_seed,
              device = ctx$device,
              use_amp = ctx$use_amp,
              batch_norm = ctx$batch_norm,
              validation_split = ctx$validation_split,
              epochs = ctx$epochs,
              batch_size = ctx$batch_size,
              tabnet_steps = ctx$tabnet_steps,
              tabnet_feature_dim = ctx$tabnet_feature_dim,
              tabnet_output_dim = ctx$tabnet_output_dim,
              tabnet_gamma = ctx$tabnet_gamma,
              tabnet_lambda_sparse = ctx$tabnet_lambda_sparse,
              dropout = ctx$dropout,
              l2_weight_decay = ctx$l2_weight_decay,
              optimizer_name = ctx$optimizer_name,
              max_grad_norm = ctx$max_grad_norm,
              heteroscedastic = ctx$heteroscedastic,
              n_bootstrap = ctx$n_bootstrap,
              internal_cv_nfolds = ctx$internal_cv_nfolds,
              dl_internal_calibration = ctx$dl_internal_calibration,
              CI_width_thresholds = ctx$CI_width_thresholds,
              high_reliability_thres = ctx$high_reliability_thres,
              low_reliability_thres = ctx$low_reliability_thres,
              system_database = ctx$system_database
            ),
          moe = AI_moe_MET(
              pheno_object = trait_pheno_data,
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              response_family = ctx$response_family,
              gmatrix = ctx$gmatrix_model_ready,
              omic1_kernel = ctx$omic1_kernel_model_ready,
              omic2_kernel = ctx$omic2_kernel_model_ready,
              omic3_kernel = ctx$omic3_kernel_model_ready,
              kernel_list = ctx$gmatrix_kernel_model_ready_list,
              omics_kernel_label = ctx$omics_kernel_label,
              met_kernel_var_explained = ctx$met_kernel_var_explained,
              met_kernel_min_ev = ctx$met_kernel_min_ev,
              met_kernel_max_pcs = ctx$met_kernel_max_pcs,
              compile_model = ctx$compile_model,
              deterministic = ctx$deterministic,
              random_seed = ctx$random_seed,
              device = ctx$device,
              use_amp = ctx$use_amp,
              batch_norm = ctx$batch_norm,
              validation_split = ctx$validation_split,
              epochs = ctx$epochs,
              batch_size = ctx$batch_size,
              moe_n_experts = ctx$moe_n_experts,
              moe_expert_hidden = ctx$moe_expert_hidden,
              moe_gate_hidden = ctx$moe_gate_hidden,
              moe_temperature = ctx$moe_temperature,
              moe_sparse_topk = ctx$moe_sparse_topk,
              moe_entropy_reg = ctx$moe_entropy_reg,
              dropout = ctx$dropout,
              l2_weight_decay = ctx$l2_weight_decay,
              optimizer_name = ctx$optimizer_name,
              max_grad_norm = ctx$max_grad_norm,
              heteroscedastic = ctx$heteroscedastic,
              n_bootstrap = ctx$n_bootstrap,
              internal_cv_nfolds = ctx$internal_cv_nfolds,
              dl_internal_calibration = ctx$dl_internal_calibration,
              CI_width_thresholds = ctx$CI_width_thresholds,
              high_reliability_thres = ctx$high_reliability_thres,
              low_reliability_thres = ctx$low_reliability_thres,
              system_database = ctx$system_database
            )
          )
        } else {
          deep_learning_model(
            pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
            geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
            geno_omic_test_object = trait_ml_dat_res[["merged_data_test"]],
            response = response,
            gen_name = ctx$gen_name,
            response_family = ctx$response_family,
            model_type = GS_model,
            optimizer_name = ctx$optimizer_name,
            use_amp = ctx$use_amp,
            max_grad_norm = ctx$max_grad_norm,
            auto_class_weights = ctx$auto_class_weights,
            cnn_neurons_per_layer = ctx$cnn_neurons_per_layer,
            cnn_kernel_size = ctx$cnn_kernel_size,
            cnn_dense_layers = ctx$cnn_dense_layers,
            cnn_use_max_pool = ctx$cnn_use_max_pool,
            cnn_pool_kernel = ctx$cnn_pool_kernel,
            cnn_pool_stride = ctx$cnn_pool_stride,
            cnn_pool_padding = ctx$cnn_pool_padding,
            cnn_learning_rate = ctx$cnn_learning_rate,
            cnn_separable = ctx$cnn_separable,
            cnn_dilations = ctx$cnn_dilations,
            cnn_use_se = ctx$cnn_use_se,
            cnn_norm_type = ctx$cnn_norm_type,
            cnn_pool_type = ctx$cnn_pool_type,
            cnn_use_global_pool = ctx$cnn_use_global_pool,
            resnet_neurons_per_block = ctx$resnet_neurons_per_block,
            resnet_blocks = ctx$resnet_blocks,
            resnet_learning_rate = ctx$resnet_learning_rate,
            ft_d_model = ctx$ft_d_model,
            ft_heads = ctx$ft_heads,
            ft_layers = ctx$ft_layers,
            ft_ff_mult = ctx$ft_ff_mult,
            ft_dropout = ctx$ft_dropout,
            ft_token_dropout = ctx$ft_token_dropout,
            ft_use_cls = ctx$ft_use_cls,
            saint_d_model = ctx$saint_d_model,
            saint_heads = ctx$saint_heads,
            saint_layers = ctx$saint_layers,
            saint_ff_mult = ctx$saint_ff_mult,
            saint_dropout = ctx$saint_dropout,
            saint_token_dropout = ctx$saint_token_dropout,
            saint_use_cls = ctx$saint_use_cls,
            use_grouping = ctx$use_grouping,
            group_trigger = ctx$group_trigger,
            group_method = ctx$group_method,
            init_group_size = ctx$init_group_size,
            max_tokens = ctx$max_tokens,
            kmeans_batch = ctx$kmeans_batch,
            kmeans_iter = ctx$kmeans_iter,
            tabnet_steps = ctx$tabnet_steps,
            tabnet_feature_dim = ctx$tabnet_feature_dim,
            tabnet_output_dim = ctx$tabnet_output_dim,
            tabnet_gamma = ctx$tabnet_gamma,
            tabnet_lambda_sparse = ctx$tabnet_lambda_sparse,
            node_trees = ctx$node_trees,
            node_depth = ctx$node_depth,
            deepfm_k = ctx$deepfm_k,
            deepfm_hidden = ctx$deepfm_hidden,
            dcn_layers = ctx$dcn_layers,
            dcn_hidden = ctx$dcn_hidden,
            nam_hidden = ctx$nam_hidden,
            nam_activation = ctx$nam_activation,
            nam_add_linear = ctx$nam_add_linear,
            nam_l1 = ctx$nam_l1,
            moe_n_experts = ctx$moe_n_experts,
            moe_expert_hidden = ctx$moe_expert_hidden,
            moe_gate_hidden = ctx$moe_gate_hidden,
            moe_temperature = ctx$moe_temperature,
            moe_sparse_topk = ctx$moe_sparse_topk,
            moe_entropy_reg = ctx$moe_entropy_reg,
            gp_use_variational = ctx$gp_use_variational,
            gp_num_inducing = ctx$gp_num_inducing,
            gp_feature_dim = ctx$gp_feature_dim,
            gp_kernel = ctx$gp_kernel,
            gp_ard = ctx$gp_ard,
            gp_lr_mult = ctx$gp_lr_mult,
            rff_features = ctx$rff_features,
            rff_lengthscale = ctx$rff_lengthscale,
            rff_deep_hidden = ctx$rff_deep_hidden,
            epochs = ctx$epochs,
            batch_size = ctx$batch_size,
            dropout = ctx$dropout,
            l2_weight_decay = ctx$l2_weight_decay,
            l2_regularizer_dp = ctx$l2_regularizer_dp,
            dropout_rate = ctx$dropout_rate,
            batch_norm = ctx$batch_norm,
            validation_split = ctx$validation_split,
            compile_model = ctx$compile_model,
            deterministic = ctx$deterministic,
            random_seed = ctx$random_seed,
            dl_n_seeds = ctx$dl_n_seeds,
            dl_seeds = ctx$dl_seeds,
            dl_seed_aggregation = ctx$dl_seed_aggregation,
            device = ctx$device,
            mlp_neurons_per_layer = ctx$mlp_neurons_per_layer,
            mlp_learning_rate = ctx$mlp_learning_rate,
            final_attention = ctx$final_attention,
            attention_across_multiple_layers = ctx$attention_across_multiple_layers,
            heteroscedastic = ctx$heteroscedastic,
            message = ctx$message,
            scaling = ctx$scaling,
            centering = ctx$centering,
            omic_count = trait_ml_dat_res[["omic_count"]],
            para_tunning = if (ctx$cross_validation) FALSE else ctx$para_tunning,
            param_grid = if (ctx$cross_validation) NULL else ctx$dpl_paras_tunning,
            CI_width_thresholds = ctx$CI_width_thresholds,
            high_reliability_thres = ctx$high_reliability_thres,
            low_reliability_thres = ctx$low_reliability_thres,
            threshold = ctx$threshold,
            target = "test_set",
            interval_width_low_threshold = ctx$interval_width_low_threshold,
            interval_width_high_threshold = ctx$interval_width_high_threshold,
            interval_width_moderate_threshold = ctx$interval_width_moderate_threshold,
            n_bootstrap = ctx$n_bootstrap,
            dl_internal_calibration = ctx$dl_internal_calibration,
            crossval = FALSE
          )
        }
      }, error = function(e) {
        cat(paste(ctx$msg, "Error in deep learning model:", conditionMessage(e), "\n"))
        NULL
      })
    } else if (is.null(res_model_output)) {
      res_model_output <- switch(
        GS_model,
        Xgboost = run_xgboost(),
        RandomForest = run_random_forest(),
        CatBoost = run_catboost(),
        LightGBM = run_lightgbm(),
        PartialLeastSquare = run_pls(),
        SupportVectorMachine = run_svm(),
        `K-NearestNeighbors` = run_knn(),
        Lasso = run_lasso(),
        Ridge_Regression = run_ridge(),
        stop(paste(ctx$msg, "Unknown GS_model:", GS_model), call. = FALSE)
      )
    }

    is_met_ml_dl_result <- isTRUE(is_multi_env_bayes) &&
      isTRUE(ctx$met_ml_dl) &&
      GS_model %in% c("CatBoost", "LightGBM", "Xgboost", "RandomForest", "mlp", "ft_transformer", "saint", "tabnet", "moe")

    if (is.null(res_model_output)) {
      cat(paste(ctx$msg, GS_model, "model fitting failed.\n"))
      res_summary_stat <- NULL
    } else {
      GS_model_fr <- replace_with_friendly_name(GS_model, ctx$friendly_name_lookup)
      res_model_output[["model_parameters"]] <-
        gp_annotate_requested_model_parameters(
          res_model_output[["model_parameters"]],
          GS_model
        )
      AI_summary_stat_process <- tryCatch({
        if (isTRUE(is_met_ml_dl_result)) {
          met_summary <- gp_met_summary_statistics(
            predicted_values = res_model_output[["predicted_values"]],
            pheno_data = res_model_output[["met_long_data"]],
            response = response,
            gen_name = ctx$gen_name,
            heter_groups = ctx$heter_groups,
            model_parameters = res_model_output[["model_parameters"]],
            feature_summary = res_model_output[["met_feature_summary"]],
            GS_model = GS_model_fr,
            response_family = ctx$response_family
          )
          res_model_output[["MET_prediction_table"]] <- met_summary[["MET_prediction_table"]]
          if (!is.null(met_summary[["MET_prediction_table"]])) {
            res_model_output[["predicted_values"]] <- met_summary[["MET_prediction_table"]]
          }
          met_total_pred <- met_summary[["across_environment_predicted_values"]] %||%
            met_summary[["Total_Predicted_value"]]
          if (!is.null(met_total_pred)) {
            res_model_output[["across_environment_predicted_values"]] <- met_total_pred
            res_model_output[["Total_Predicted_value"]] <- met_total_pred
          }
          met_response_family <- gp_resolve_response_family(ctx$response_family)
          if (identical(met_response_family, "gaussian")) {
            uncertainty_table <- met_summary[["MET_prediction_table"]]
            uncertainty_source <- character()
            if (is.data.frame(uncertainty_table) &&
                "Standard_error" %in% names(uncertainty_table) &&
                any(is.finite(suppressWarnings(as.numeric(
                  uncertainty_table[["Standard_error"]]
                ))))) {
              if ("Prediction_uncertainty_source" %in% names(uncertainty_table)) {
                uncertainty_source <- unique(as.character(
                  uncertainty_table[["Prediction_uncertainty_source"]]
                ))
                uncertainty_source <- uncertainty_source[
                  !is.na(uncertainty_source) & nzchar(uncertainty_source)
                ]
              }
            }
            res_model_output[["prediction_uncertainty_source"]] <- if (
              length(uncertainty_source)
            ) {
              paste(uncertainty_source, collapse = ";")
            } else {
              "MET ML/DL uncertainty unavailable without held-out calibration"
            }
            matrix_names <- grep(
              "^(Prediction_|predictive_covariance_basis$)",
              names(met_summary),
              value = TRUE
            )
            for (matrix_name in matrix_names) {
              res_model_output[[matrix_name]] <- met_summary[[matrix_name]]
            }
          }
          res_model_output[["variance_components"]] <-
            gp_predictive_model_variance_components(
              met_summary[["MET_prediction_table"]],
              response_family = met_response_family
            )
          res_model_output[["Variance_components"]] <-
            res_model_output[["variance_components"]]
          if (is.null(res_model_output[["diagnostic_plots"]])) {
            res_model_output[["diagnostic_plots"]] <- gp_met_true_prediction_plot(
              predicted_values = res_model_output[["predicted_values"]],
              pheno_data = res_model_output[["met_long_data"]],
              response = response,
              gen_name = ctx$gen_name,
              heter_groups = ctx$heter_groups,
              model_label = GS_model_fr,
              response_family = ctx$response_family
            )
          }
          met_summary[setdiff(
            names(met_summary),
            c("MET_prediction_table", "across_environment_predicted_values", "Total_Predicted_value")
          )]
        } else {
          summary_statistics_AI(
            predicted_object = res_model_output[["predicted_values"]],
            pheno_object = trait_ml_dat_res[["pheno_clean_data"]],
            response = response,
            test_set = trait_ml_dat_res[["test_set"]],
            geno_omic_object = trait_ml_dat_res[["merged_data"]][["merge_data"]],
            eval_metrics = ctx$eval_metrics,
            model_parameters = res_model_output[["model_parameters"]],
            GS_model = GS_model_fr,
            response_family = ctx$response_family
          )
        }
      }, error = function(e) {
        cat(paste(ctx$msg, "Error:", conditionMessage(e), "\n"))
        NULL
      })
      if (!is.null(AI_summary_stat_process)) {
        res_summary_stat <- AI_summary_stat_process
      } else {
        cat(paste(ctx$msg, "Error processing machine learning model output.\n"))
      }
    }
  }

  GS_model_fr <- replace_with_friendly_name(GS_model, ctx$friendly_name_lookup)
  if (!is.null(res_model_output) &&
      (!is.null(ctx$feature_score_metadata) || !is.null(ctx$feature_selected))) {
    trait_feature_meta <- NULL
    if (!is.null(ctx$feature_score_metadata)) {
      trait_feature_meta <- ctx$feature_score_metadata[
        as.character(ctx$feature_score_metadata$trait) == response,
        ,
        drop = FALSE
      ]
    }
    selected_names <- trait_ml_dat_res[["feature_selected"]]
    if (!is.null(trait_feature_meta) && nrow(trait_feature_meta)) {
      res_model_output[["feature_score_metadata"]] <- trait_feature_meta
    }
    if (!is.null(explicit_feature_meta) && nrow(explicit_feature_meta)) {
      res_model_output[["feature_selection_metadata"]] <- explicit_feature_meta
    }
    if (!is.null(trait_feature_meta) && nrow(trait_feature_meta) ||
        !is.null(explicit_feature_meta) && nrow(explicit_feature_meta) ||
        length(selected_names)) {
      if (!length(selected_names) && !is.null(explicit_feature_meta)) {
        selected_names <- unique(as.character(explicit_feature_meta[["predictor"]]))
      }
      n_predictors_used <- if (!is.null(trait_ml_dat_res[["merged_data"]][["merge_data"]])) {
        ncol(trait_ml_dat_res[["merged_data"]][["merge_data"]])
      } else {
        NA_integer_
      }
      feature_param_rows <- data.frame(
        stat = c(
          if (!is.null(trait_feature_meta) && nrow(trait_feature_meta)) "feature_scoring_model" else character(),
          if (!is.null(ctx$feature_k)) "feature_k" else character(),
          if (length(selected_names)) "feature_selected" else character(),
          "n_predictors_used"
        ),
        summary = c(
          if (!is.null(trait_feature_meta) && nrow(trait_feature_meta)) {
            as.character(unique(trait_feature_meta$scoring_model))[1]
          } else {
            character()
          },
          if (!is.null(ctx$feature_k)) as.character(ctx$feature_k) else character(),
          if (length(selected_names)) paste(selected_names, collapse = ",") else character(),
          as.character(n_predictors_used)
        ),
        stringsAsFactors = FALSE
      )
      if (is.data.frame(res_model_output[["model_parameters"]])) {
        res_model_output[["model_parameters"]] <- rbind(
          res_model_output[["model_parameters"]],
          feature_param_rows
        )
      } else {
        res_model_output[["model_parameters"]] <- feature_param_rows
      }
    }
  }
  res_model_output <- gp_standardize_model_prediction_outputs(
    res_model_output = res_model_output,
    gen_name = ctx$gen_name
  )
  # Multi-environment: label each line x environment row as Train (observed),
  # Test (line never observed) or Unobserved (line observed elsewhere), the
  # same way for every engine.
  if (!is.null(ctx$heter_groups) && is.data.frame(trait_pheno_data) &&
      ctx$heter_groups %in% names(trait_pheno_data) &&
      anyDuplicated(as.character(trait_pheno_data[[ctx$gen_name]])) > 0L) {
    for (nm in c("predicted_values", "Predicted_value")) {
      if (is.data.frame(res_model_output[[nm]])) {
        res_model_output[[nm]] <- gp_met_relabel_predictions(
          res_model_output[[nm]], trait_pheno_data, ctx$gen_name, ctx$heter_groups, response,
          input_keys = ctx$met_input_keys
        )
      }
    }
  }
  connectivity_outputs <- gp_attach_connectivity_outputs(
    res_model_output = res_model_output,
    res_summary_stat = res_summary_stat,
    pheno_data = trait_pheno_data
  )
  list(
    GS_model = GS_model_fr,
    res_model_output = connectivity_outputs[["res_model_output"]],
    res_summary_stat = connectivity_outputs[["res_summary_stat"]]
  )
}

gp_finalize_model_execute_results <- function(results,
                                              GS_model,
                                              gen_name = NULL,
                                              heter_groups = NULL,
                                              pheno_data = NULL,
                                              response_family = NULL,
                                              best_models_ggplot_rep,
                                              best_models_ggplot_mean,
                                              cv_results_predicted_vs_observed,
                                              geno_qc_stat,
                                              cv_results_processed,
                                              cv_results,
                                              feature_selected = NULL,
                                              feature_score_metadata = NULL,
                                              run_metadata = NULL,
                                              system_database,
                                              plot_extension,
                                              plot_width,
                                              plot_height,
                                              plot_units,
                                              plot_dpi) {
  all_results <- list()
  mainDirt <- getwd()

  name_cv_entries <- function(x) {
    if (is.null(x) || !length(x) || !is.list(x)) {
      return(x)
    }
    nms <- vapply(seq_along(x), function(i) {
      entry <- x[[i]]
      if (!is.list(entry)) {
        return(paste0("result_", i))
      }
      trait <- as.character(entry[["trait"]] %||% entry[["response"]] %||% paste0("trait", i))[1L]
      model <- as.character(entry[["model"]] %||% entry[["model_canonical"]] %||% "model")[1L]
      rep_i <- as.character(entry[["rep"]] %||% entry[["replication"]] %||% i)[1L]
      make.names(paste(trait, model, paste0("rep", rep_i), sep = "_"))
    }, character(1))
    names(x) <- make.unique(nms)
    x
  }

  cv_results <- name_cv_entries(cv_results)

  safe_index <- function(x, idx) {
    if (is.null(x) || length(x) < idx) {
      return(NULL)
    }
    x[[idx]]
  }

  cv_plot_bundle_for_trait <- function(cv_results_predicted_vs_observed, trait_name, idx) {
    plot_bundle <- if (is.list(cv_results_predicted_vs_observed)) {
      cv_results_predicted_vs_observed$predicted_vs_observed_plots
    } else {
      NULL
    }
    if (is.null(plot_bundle) || !length(plot_bundle)) {
      return(NULL)
    }
    if (!is.null(names(plot_bundle)) && trait_name %in% names(plot_bundle)) {
      return(plot_bundle[[trait_name]])
    }
    safe_index(plot_bundle, idx)
  }

  model_results_by_trait <- list()
  model_results_by_model <- list()
  processed_results_by_trait <- list()
  processed_results_by_model <- list()
  trait_model_counts <- list()

  add_model_aliases <- function(processed_result, model_label, model_result) {
    model_aliases <- unique(c(as.character(model_label), make.names(as.character(model_label))))
    model_aliases <- model_aliases[nzchar(model_aliases)]
    if (!length(model_aliases)) {
      return(processed_result)
    }
    existing_models <- processed_result[["models"]]
    if (!is.list(existing_models)) {
      existing_models <- list()
    }
    for (alias in model_aliases) {
      existing_models[[alias]] <- model_result
    }
    processed_result[["models"]] <- existing_models
    processed_result
  }

  for (res in seq_along(results)) {
    tryCatch({
      model_label <- results[[res]][["GS_model"]] %||% GS_model
      model_label <- as.character(model_label)[1L]
      trait_name <- results[[res]][["trait"]] %||%
        results[[res]][["response"]] %||%
        names(results)[res] %||%
        paste("trait", res, sep = "_")
      trait_name <- as.character(trait_name)[1L]
      if (is.na(trait_name) || !nzchar(trait_name)) {
        trait_name <- paste("trait", res, sep = "_")
      }
      # Phase 3.26: store under the user's original model name only. R's
      # [[name]] and backtick-`$name` accessors handle dashes correctly,
      # so duplicating entries under both the original and make.names()
      # alias (as Phase 3.25 did) is wrong -- it produces double rows in
      # names(model_results_by_model). make.names() output is reserved for
      # filesystem-safe plot names only.
      list_key <- if (!is.na(model_label) && nzchar(model_label)) {
        as.character(model_label)
      } else {
        paste0("model_", res)
      }
      model_key <- if (!is.na(model_label) && nzchar(model_label)) {
        make.names(as.character(model_label))
      } else {
        paste0("model_", res)
      }
      plot_name <- if (!is.null(names(results)[res]) && nzchar(names(results)[res])) {
        names(results)[res]
      } else {
        make.names(paste(trait_name, model_key, sep = "_"))
      }
      public_model_output <- gp_standardize_public_model_result(
        if ("res_model_output" %in% names(results[[res]])) results[[res]][["res_model_output"]] else NULL,
        gen_name = gen_name,
        heter_groups = heter_groups,
        pheno_data = pheno_data,
        response = trait_name,
        model = model_label,
        response_family = response_family
      )

      processed_result <- results_handling(
        GS_model = model_label,
        res_model_output = public_model_output,
        res_summary_stat = if ("res_summary_stat" %in% names(results[[res]])) results[[res]][["res_summary_stat"]] else NULL,
        res_plot = best_models_ggplot_rep,
        res_plot_mean = best_models_ggplot_mean,
        res_plot_result_diagnostic = cv_plot_bundle_for_trait(cv_results_predicted_vs_observed, trait_name, res),
        test_diagonistic_plots = if (!is.null(results[[res]]$res_model_output)) results[[res]]$res_model_output$diagnostic_plots else NULL,
        res_mod_results_cv_per_trait_model = if (!is.null(cv_results_predicted_vs_observed)) cv_results_predicted_vs_observed$mod_res_per_trait_per_model else NULL,
        geno_qc_stat = geno_qc_stat,
        res_plot_result_diagnostic_cv_only = NULL,
        cv_results_processed = cv_results_processed,
        cv_results_raw = cv_results,
        feature_selected = feature_selected,
        feature_score_metadata = feature_score_metadata,
        run_metadata = run_metadata,
        system_database = system_database,
        plot_filename = plot_name,
        plot_extension = plot_extension,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_units = plot_units,
        plot_dpi = plot_dpi
      )

      if (!is.null(model_label) && length(model_label) == 1L && !is.na(model_label) && nzchar(model_label)) {
        processed_result[["model_name"]] <- as.character(model_label)
      }
      processed_result[["trait"]] <- trait_name
      processed_result[["model_key"]] <- model_key
      processed_result[["is_cv_best"]] <- isTRUE(results[[res]][["is_cv_best"]])

      model_result <- processed_result[["model_results"]]
      if (is.list(model_result)) {
        model_result <- gp_standardize_public_model_result(
          model_result,
          gen_name = gen_name,
          heter_groups = heter_groups,
          model = model_label,
          response_family = response_family
        )
        processed_result[["model_results"]] <- model_result
      }
      processed_result <- add_model_aliases(processed_result, model_label, model_result)

      if (is.null(model_results_by_trait[[trait_name]])) {
        model_results_by_trait[[trait_name]] <- list()
      }
      model_results_by_trait[[trait_name]][[list_key]] <- model_result

      if (is.null(model_results_by_model[[list_key]])) {
        model_results_by_model[[list_key]] <- list()
      }
      model_results_by_model[[list_key]][[trait_name]] <- model_result

      if (is.null(processed_results_by_trait[[trait_name]])) {
        processed_results_by_trait[[trait_name]] <- list()
      }
      processed_results_by_trait[[trait_name]][[list_key]] <- processed_result

      if (is.null(processed_results_by_model[[list_key]])) {
        processed_results_by_model[[list_key]] <- list()
      }
      processed_results_by_model[[list_key]][[trait_name]] <- processed_result

      if (is.null(trait_model_counts[[trait_name]])) {
        trait_model_counts[[trait_name]] <- 0L
      }
      trait_model_counts[[trait_name]] <- trait_model_counts[[trait_name]] + 1L

    }, error = function(e) {
      setwd(mainDirt)
      message(paste("Error processing result", res, ":", e$message))
    })
  }

  all_results[["cv_results_raw"]] <- cv_results
  all_results[["cv_results_processed"]] <- cv_results_processed
  all_results[["cv_results_predicted_vs_observed"]] <- cv_results_predicted_vs_observed
  all_results[["model_results_by_trait"]] <- model_results_by_trait
  all_results[["model_results_by_model"]] <- model_results_by_model
  if (length(model_results_by_model) == 1L) {
    single_model_block <- model_results_by_model[[1L]]
    if (is.list(single_model_block) && length(single_model_block) == 1L) {
      all_results[["model_results"]] <- single_model_block[[1L]]
    }
  }
  if (length(processed_results_by_model) == 1L) {
    single_processed_block <- processed_results_by_model[[1L]]
    if (is.list(single_processed_block) && length(single_processed_block) == 1L) {
      single_processed <- single_processed_block[[1L]]
      for (nm in c(
        "summary_statistic",
        "res_plot",
        "res_plot_result_diagnostic",
        "res_plot_result_diagnostic_cv_only",
        "cv_results_predicted_vs_observed",
        "test_diagonistic_plots",
        "res_mod_results_cv_per_trait_model",
        "feature_selected",
        "feature_score_metadata",
        "export_status",
        "model_name",
        "trait",
        "model_key",
        "is_cv_best"
      )) {
        if (!is.null(single_processed[[nm]])) {
          all_results[[nm]] <- single_processed[[nm]]
        }
      }
    }
  }
  all_results[["run_metadata"]] <- run_metadata
  all_results[["feature_selected"]] <- feature_selected
  all_results[["feature_score_metadata"]] <- feature_score_metadata
  all_results[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(results, function(x) {
      x[["res_model_output"]][["feature_selection_metadata"]] %||% NULL
    })
  )
  all_results
}
