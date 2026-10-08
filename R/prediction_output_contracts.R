#' Describe the canonical PredictProR output contract
#'
#' Model families share the same four top-level result fields. Prediction
#' tables then use a Gaussian, classification, or hybrid schema, with optional
#' environment and trait identifier columns.
#'
#' @param response_family One of `"gaussian"`, `"classification"`, or
#'   `"auto"`.
#' @param task One of `"single_trait"`, `"multi_trait"`,
#'   `"multi_environment"`, `"multi_trait_multi_environment"`, or
#'   the corresponding `"hybrid*"` task.
#' @param gen_name Output genotype identifier column name.
#' @param heter_groups Output environment identifier column name.
#'
#' @return A list describing required and optional result fields and ordered
#'   prediction columns.
#' @export
prediction_output_standard <- function(response_family = c("gaussian", "classification", "auto"),
                                       task = c(
                                         "single_trait", "multi_trait", "multi_environment",
                                         "multi_trait_multi_environment", "hybrid",
                                         "hybrid_multi_trait", "hybrid_multi_environment",
                                         "hybrid_multi_trait_multi_environment"
                                       ),
                                       gen_name = "GID",
                                       heter_groups = "Env") {
  response_family <- match.arg(response_family)
  task <- match.arg(task)
  include_env <- task %in% c(
    "multi_environment", "multi_trait_multi_environment",
    "hybrid_multi_environment", "hybrid_multi_trait_multi_environment"
  )
  include_trait <- task %in% c(
    "multi_trait", "multi_trait_multi_environment",
    "hybrid_multi_trait", "hybrid_multi_trait_multi_environment"
  )
  if (startsWith(task, "hybrid") && identical(response_family, "classification")) {
    stop("Hybrid prediction currently supports gaussian responses only.", call. = FALSE)
  }
  gen_name <- as.character(gen_name %||% "GID")[[1L]]
  heter_groups <- as.character(heter_groups %||% "Env")[[1L]]

  prediction_columns <- if (startsWith(task, "hybrid")) {
    gp_hybrid_prediction_columns(include_env = include_env, include_trait = include_trait)
  } else if (identical(response_family, "classification")) {
    c(
      gen_name,
      if (include_env) heter_groups,
      if (include_trait) "Trait",
      "Predicted_class", "Train_Test_Label", "Observed_class",
      "Prediction_confidence", "Classification_uncertainty",
      "Reliability", "Reliability_remarks", "Probability_<class>"
    )
  } else if (identical(response_family, "gaussian")) {
    gp_gaussian_prediction_columns(
      include_env = include_env,
      include_trait = include_trait,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  } else {
    list(
      gaussian = gp_gaussian_prediction_columns(
        include_env = include_env,
        include_trait = include_trait,
        gen_name = gen_name,
        heter_groups = heter_groups
      ),
      classification = c(
        gen_name,
        if (include_env) heter_groups,
        if (include_trait) "Trait",
        "Predicted_class", "Train_Test_Label", "Observed_class",
        "Prediction_confidence", "Classification_uncertainty",
        "Reliability", "Reliability_remarks", "Probability_<class>"
      )
    )
  }

  list(
    top_level_required = c(
      "model_parameters", "predicted_values", "diagnostic_plots",
      "variance_components"
    ),
    top_level_optional = c(
      "across_environment_predicted_values", "covariance_components",
      "correlation_components", "Genetic_covariance_traits",
       "Genetic_correlation_traits", "Residual_covariance_traits",
       "Residual_correlation_traits", "Genetic_covariance_environments",
       "Genetic_correlation_environments",
       "Residual_covariance_environments", "Residual_correlation_environments",
       "Prediction_covariance_traits", "Prediction_correlation_traits",
       "Prediction_error_covariance_traits", "Prediction_error_correlation_traits",
       "Prediction_covariance_environments", "Prediction_correlation_environments",
       "Prediction_error_covariance_environments", "Prediction_error_correlation_environments",
       "Total_genetic_covariance_traits", "Total_genetic_correlation_traits",
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
       "Genetic_correlation_trait_environment_SE",
      "Genetic_covariance_trait_environment",
      "Genetic_correlation_trait_environment",
      "variance_component_intervals",
      "prediction_uncertainty_source",
      "predictive_covariance_basis",
      "multitrait_prediction_wide", "multitrait_trait_counts"
    ),
    prediction_columns = prediction_columns,
    variance_component_columns = gp_variance_component_columns(
      include_trait = include_trait,
      include_env = include_env
    ),
    notes = c(
      "PEV equals Standard_error squared wherever both are finite.",
      "Reliability is on the [0, 1] scale; inspect Reliability_basis for its estimand.",
      paste(
        "For ML/DL, Reliability is available only with finite target-specific PEV and is the",
        "non-genetic predictive precision index V_y/(V_y + PEV_i). Prediction_stability is a",
        "separate descriptive bootstrap-resampling quantity and must not be substituted when PEV is unavailable."
      ),
      "ML/DL PEV is a held-out estimate of predictive MSE, not mixed-model breeding-value PEV.",
      "ML/DL prediction intervals target marginal coverage under exchangeability, not conditional coverage.",
      "ML/DL genetic/residual variance and heritability are unavailable unless an explicit random-effects model estimates them.",
      "Prediction covariance is descriptive and is not genetic covariance.",
      "Classification probability columns are named Probability_<class>.",
      "model_execute() indexes these model results by both model and trait."
    )
  )
}

gp_output_contract_collect_results <- function(x) {
  if (!is.list(x)) {
    return(list(input = x))
  }
  if (all(c("predicted_values", "model_parameters") %in% names(x))) {
    return(list(model_result = x))
  }
  by_model <- x[["model_results_by_model"]]
  if (is.list(by_model) && length(by_model)) {
    out <- list()
    for (model_name in names(by_model) %||% as.character(seq_along(by_model))) {
      traits <- by_model[[model_name]]
      if (!is.list(traits)) {
        next
      }
      for (trait_name in names(traits) %||% as.character(seq_along(traits))) {
        out[[paste(model_name, trait_name, sep = "/")]] <- traits[[trait_name]]
      }
    }
    if (length(out)) {
      return(out)
    }
  }
  if (is.list(x[["model_results"]])) {
    return(list(model_result = x[["model_results"]]))
  }
  list(input = x)
}

gp_output_contract_task <- function(pred, heter_groups = NULL) {
  if (all(c("Hybrid_ID", "Female", "Male") %in% names(pred))) {
    env_candidates <- unique(c(as.character(heter_groups %||% character()), "Env"))
    has_env <- any(env_candidates %in% names(pred))
    has_trait <- "Trait" %in% names(pred)
    if (has_env && has_trait) return("hybrid_multi_trait_multi_environment")
    if (has_env) return("hybrid_multi_environment")
    if (has_trait) return("hybrid_multi_trait")
    return("hybrid")
  }
  env_candidates <- unique(c(as.character(heter_groups %||% character()), "Env"))
  has_env <- any(env_candidates %in% names(pred))
  has_trait <- "Trait" %in% names(pred)
  if (has_env && has_trait) {
    return("multi_trait_multi_environment")
  }
  if (has_env) {
    return("multi_environment")
  }
  if (has_trait) {
    return("multi_trait")
  }
  "single_trait"
}

gp_output_contract_validate_one <- function(result,
                                            path,
                                            response_family = "auto",
                                            gen_name = NULL,
                                            heter_groups = NULL,
                                            tolerance = 1e-8) {
  issues <- character()
  required <- prediction_output_standard()$top_level_required
  if (!is.list(result)) {
    issues <- "model result must be a list"
    return(data.frame(
      path = path, valid = FALSE, response_family = NA_character_,
      task = NA_character_, n_predictions = NA_integer_,
      issues = paste(issues, collapse = "; "), stringsAsFactors = FALSE
    ))
  }
  missing_slots <- setdiff(required, names(result))
  if (length(missing_slots)) {
    issues <- c(issues, paste("missing result fields", paste(missing_slots, collapse = ", ")))
  }
  pred <- result[["predicted_values"]]
  if (!is.data.frame(pred)) {
    issues <- c(issues, "predicted_values must be a data frame")
    return(data.frame(
      path = path, valid = FALSE, response_family = NA_character_,
      task = NA_character_, n_predictions = NA_integer_,
      issues = paste(unique(issues), collapse = "; "), stringsAsFactors = FALSE
    ))
  }

  family <- response_family
  if (identical(family, "auto")) {
    family <- if (gp_is_classification_prediction_table(pred)) "classification" else "gaussian"
  }
  task <- gp_output_contract_task(pred, heter_groups = heter_groups)
  id_name <- if (startsWith(task, "hybrid")) {
    "Hybrid_ID"
  } else if (!is.null(gen_name) && nzchar(as.character(gen_name)[1L])) {
    as.character(gen_name)[1L]
  } else {
    names(pred)[[1L]] %||% "GID"
  }
  env_name <- if (!is.null(heter_groups) && nzchar(as.character(heter_groups)[1L])) {
    as.character(heter_groups)[1L]
  } else {
    "Env"
  }
  standard <- prediction_output_standard(
    response_family = family,
    task = task,
    gen_name = id_name,
    heter_groups = env_name
  )
  expected <- standard$prediction_columns
  if (identical(family, "classification")) {
    expected_base <- expected[expected != "Probability_<class>"]
    probability_columns <- grep("^Probability_", names(pred), value = TRUE)
    expected <- c(expected_base, probability_columns)
    if (!length(probability_columns)) {
      issues <- c(issues, "classification output requires Probability_<class> columns")
    }
  }
  if (!identical(names(pred), expected)) {
    issues <- c(issues, "prediction columns do not match the canonical order")
  }

  ids <- as.character(pred[[id_name]] %||% character())
  if (length(ids) != nrow(pred) || anyNA(ids) || any(!nzchar(ids))) {
    issues <- c(issues, paste0(id_name, " must be complete and non-empty"))
  }

  if (identical(family, "gaussian")) {
    numeric_required <- c(
      "Predicted_value", "Observed_value", "Standard_error", "PEV",
      "lower_bound", "upper_bound", "Uncertainty", "Reliability",
      "Reliability_variance_input", "Reliability_reference_variance"
    )
    if (!all(vapply(pred[intersect(numeric_required, names(pred))], is.numeric, logical(1)))) {
      issues <- c(issues, "Gaussian estimate and uncertainty columns must be numeric")
    }
    if (all(c("Standard_error", "PEV") %in% names(pred))) {
      idx <- is.finite(pred$Standard_error) & is.finite(pred$PEV)
      if (any(idx) && any(abs(pred$PEV[idx] - pred$Standard_error[idx]^2) >
                          tolerance * (1 + abs(pred$PEV[idx])))) {
        issues <- c(issues, "PEV must equal Standard_error squared")
      }
    }
    if (all(c("lower_bound", "Predicted_value", "upper_bound") %in% names(pred))) {
      idx <- is.finite(pred$lower_bound) & is.finite(pred$Predicted_value) & is.finite(pred$upper_bound)
      if (any(idx) && any(pred$lower_bound[idx] > pred$Predicted_value[idx] + tolerance |
                          pred$upper_bound[idx] < pred$Predicted_value[idx] - tolerance)) {
        issues <- c(issues, "prediction intervals must contain Predicted_value")
      }
    }
  } else {
    probability_columns <- grep("^Probability_", names(pred), value = TRUE)
    probabilities <- unlist(pred[probability_columns], use.names = FALSE)
    finite_probabilities <- probabilities[is.finite(probabilities)]
    if (length(finite_probabilities) &&
        any(finite_probabilities < -tolerance | finite_probabilities > 1 + tolerance)) {
      issues <- c(issues, "classification probabilities must be in [0, 1]")
    }
  }

  if ("Reliability" %in% names(pred)) {
    reliability <- pred$Reliability[is.finite(pred$Reliability)]
    if (length(reliability) && any(reliability < -tolerance | reliability > 1 + tolerance)) {
      issues <- c(issues, "Reliability must be in [0, 1]")
    }
  }
  if (identical(family, "gaussian") && all(c(
    "Reliability", "Reliability_variance_input",
    "Reliability_reference_variance", "Reliability_basis"
  ) %in% names(pred))) {
    predictive_alias <- gp_contract_predictive_stability_basis(pred[["Reliability_basis"]])
    if (any(predictive_alias, na.rm = TRUE)) {
      basis <- as.character(pred[["Reliability_basis"]])
      predictive_pev_basis <- predictive_alias & grepl(
        "predictive precision index|target-specific predictive PEV|Reliability unavailable",
        ifelse(is.na(basis), "", basis),
        ignore.case = TRUE
      )
      pev_unavailable <- predictive_pev_basis &
        (!is.finite(pred[["PEV"]]) | pred[["PEV"]] < 0)
      if (any(is.finite(pred[["Reliability"]][pev_unavailable]))) {
        issues <- c(issues, "ML/DL Reliability must be NA when target-specific PEV is unavailable")
      }

      formula_idx <- predictive_pev_basis &
        is.finite(pred[["PEV"]]) & pred[["PEV"]] >= 0 &
        is.finite(pred[["Reliability_reference_variance"]]) &
        pred[["Reliability_reference_variance"]] > 0
      if (any(!is.finite(pred[["Reliability"]][formula_idx]))) {
        issues <- c(issues, "ML/DL Reliability must be finite when target-specific PEV and its reference variance are finite")
      }
      if (any(formula_idx)) {
        input_matches_pev <- abs(
          pred[["Reliability_variance_input"]][formula_idx] - pred[["PEV"]][formula_idx]
        ) <= tolerance * (1 + abs(pred[["PEV"]][formula_idx]))
        if (any(!is.finite(input_matches_pev) | !input_matches_pev)) {
          issues <- c(issues, "ML/DL Reliability_variance_input must equal target-specific PEV")
        }
        expected_reliability <-
          pred[["Reliability_reference_variance"]][formula_idx] /
          (pred[["Reliability_reference_variance"]][formula_idx] +
             pred[["PEV"]][formula_idx])
        if (any(abs(pred[["Reliability"]][formula_idx] - expected_reliability) >
                tolerance * (1 + abs(expected_reliability)))) {
          issues <- c(issues, "ML/DL Reliability must equal V_y/(V_y + PEV_i)")
        }
      }

      legacy_alias <- predictive_alias & !predictive_pev_basis
      if (any(!is.finite(pred[["Reliability"]][legacy_alias]))) {
        issues <- c(issues, "Legacy ML/DL predictive-stability Reliability must be finite")
      }
      legacy_formula_idx <- legacy_alias &
        is.finite(pred[["Reliability"]]) &
        is.finite(pred[["Reliability_variance_input"]]) &
        pred[["Reliability_variance_input"]] >= 0 &
        is.finite(pred[["Reliability_reference_variance"]]) &
        pred[["Reliability_reference_variance"]] > 0
      if (any(legacy_formula_idx)) {
        expected_reliability <- 1 -
          pred[["Reliability_variance_input"]][legacy_formula_idx] /
          pred[["Reliability_reference_variance"]][legacy_formula_idx]
        expected_reliability <- pmax(0, pmin(1, expected_reliability))
        if (any(abs(pred[["Reliability"]][legacy_formula_idx] - expected_reliability) >
                tolerance * (1 + abs(expected_reliability)))) {
          issues <- c(issues, "Legacy ML/DL Reliability must equal clipped 1 - U_i/V_yhat")
        }
      }
      if ("Prediction_stability" %in% names(pred)) {
        alias_idx <- legacy_alias & is.finite(pred[["Prediction_stability"]]) &
          is.finite(pred[["Reliability"]])
        if (any(alias_idx) && any(abs(
          pred[["Reliability"]][alias_idx] - pred[["Prediction_stability"]][alias_idx]
        ) > tolerance)) {
          issues <- c(issues, "Legacy ML/DL Reliability must equal Prediction_stability")
        }
      }
      if (any(!grepl(
        "not genetic reliability",
        pred[["Reliability_basis"]][predictive_alias],
        fixed = TRUE
      ))) {
        issues <- c(issues, "ML/DL Reliability_basis must say not genetic reliability")
      }
    }
  }
  vc <- result[["variance_components"]]
  if (!is.data.frame(vc)) {
    issues <- c(issues, "variance_components must be a data frame")
  } else if (!all(c("Component", "Components", "Standard_error") %in% names(vc))) {
    issues <- c(issues, "variance_components is missing canonical columns")
  }

  data.frame(
    path = path,
    valid = length(issues) == 0L,
    response_family = family,
    task = task,
    n_predictions = nrow(pred),
    issues = paste(unique(issues), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

#' Validate PredictProR model output
#'
#' Accepts either one canonical model result or the indexed list returned by
#' `model_execute()`. The validator checks field names, prediction-column
#' order, uncertainty identities, reliability bounds, and variance-component
#' schema.
#'
#' @param x A PredictProR model result or complete `model_execute()` result.
#' @param response_family One of `"auto"`, `"gaussian"`, or
#'   `"classification"`.
#' @param gen_name Optional genotype identifier column name.
#' @param heter_groups Optional environment identifier column name.
#' @param tolerance Numerical comparison tolerance.
#' @param strict Logical; stop when any result violates the contract.
#'
#' @return Invisibly returns one validation row per model/trait result.
#' @rdname prediction_output_standard
#' @export
validate_prediction_output <- function(x,
                                       response_family = c("auto", "gaussian", "classification"),
                                       gen_name = NULL,
                                       heter_groups = NULL,
                                       tolerance = 1e-8,
                                       strict = TRUE) {
  response_family <- match.arg(response_family)
  results <- gp_output_contract_collect_results(x)
  report <- do.call(rbind, lapply(seq_along(results), function(i) {
    result_name <- names(results)[[i]] %||% paste0("result_", i)
    gp_output_contract_validate_one(
      result = results[[i]],
      path = result_name,
      response_family = response_family,
      gen_name = gen_name,
      heter_groups = heter_groups,
      tolerance = tolerance
    )
  }))
  rownames(report) <- NULL
  if (isTRUE(strict) && any(!report$valid)) {
    failures <- paste0(report$path[!report$valid], ": ", report$issues[!report$valid])
    stop(
      paste("PredictProR output contract validation failed:", paste(failures, collapse = " | ")),
      call. = FALSE
    )
  }
  invisible(report)
}
