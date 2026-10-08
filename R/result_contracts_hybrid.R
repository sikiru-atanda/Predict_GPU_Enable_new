gp_hybrid_prediction_id_columns <- function(include_env = FALSE, include_trait = FALSE) {
  cols <- c("Hybrid_ID", "Female", "Male")
  if (isTRUE(include_env)) {
    cols <- c(cols, "Env")
  }
  if (isTRUE(include_trait)) {
    cols <- c(cols, "Trait")
  }
  cols
}

gp_hybrid_prediction_columns <- function(include_env = FALSE, include_trait = FALSE) {
  c(
    gp_hybrid_prediction_id_columns(include_env = include_env, include_trait = include_trait),
    "Predicted_value",
    "Train_Test_Label",
    "Observed_value",
    "Standard_error",
    "PEV",
    "PEV_basis",
    "Prediction_uncertainty_source",
    "Prediction_interval_method",
    "Prediction_interval_nominal_coverage",
    "Prediction_interval_calibration_n",
    "lower_bound",
    "upper_bound",
    "Uncertainty",
    "Uncertainty_remarks",
    "Prediction_stability",
    "Prediction_stability_remarks",
    "Prediction_stability_reference_variance",
    "Reliability",
    "Reliability_remarks",
    "Reliability_variance_input",
    "Reliability_reference_variance",
    "Reliability_basis"
  )
}

gp_is_hybrid_prediction_table <- function(x) {
  df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    return(FALSE)
  }
  has_parent_cols <- any(c("Female", "female_parent", "FemaleParent") %in% names(df)) &&
    any(c("Male", "male_parent", "MaleParent") %in% names(df))
  has_hybrid_id <- any(c("Hybrid_ID", "HybridID", "HybridCross", "Cross", "GID") %in% names(df))
  has_parent_cols && has_hybrid_id
}

gp_format_hybrid_prediction_table <- function(x,
                                              hybrid_id = "Hybrid_ID",
                                              female_parent = "Female",
                                              male_parent = "Male",
                                              include_env = NULL,
                                              include_trait = NULL) {
  df <- gp_format_gaussian_prediction_table(
    x = x,
    gen_name = hybrid_id,
    include_env = include_env,
    include_trait = include_trait
  )
  if (!is.data.frame(df)) {
    return(df)
  }
  original <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(original)) {
    return(df)
  }
  id_candidates <- unique(stats::na.omit(c("Hybrid_ID", hybrid_id, "HybridID", "HybridCross", "Cross", "GID")))
  female_candidates <- unique(stats::na.omit(c("Female", female_parent, "female_parent", "FemaleParent")))
  male_candidates <- unique(stats::na.omit(c("Male", male_parent, "male_parent", "MaleParent")))
  id_col <- gp_contract_first_name(names(original), id_candidates)
  female_col <- gp_contract_first_name(names(original), female_candidates)
  male_col <- gp_contract_first_name(names(original), male_candidates)
  if (!is.na(id_col)) {
    df[["Hybrid_ID"]] <- as.character(original[[id_col]])
  } else {
    df[["Hybrid_ID"]] <- df[["GID"]]
  }
  df[["Female"]] <- if (!is.na(female_col)) as.character(original[[female_col]]) else NA_character_
  df[["Male"]] <- if (!is.na(male_col)) as.character(original[[male_col]]) else NA_character_
  include_env <- "Env" %in% names(df)
  include_trait <- "Trait" %in% names(df)
  ordered_cols <- gp_hybrid_prediction_columns(include_env = include_env, include_trait = include_trait)
  df[, ordered_cols, drop = FALSE]
}

gp_hybrid_variance_component_columns <- function(include_trait = FALSE, include_env = FALSE) {
  c(
    if (isTRUE(include_trait)) "Trait",
    if (isTRUE(include_env)) "Env",
    "Hybrid_component",
    "Component",
    "Components",
    "Standard_error"
  )
}

gp_empty_hybrid_variance_components <- function() {
  data.frame(
    Hybrid_component = character(),
    Component = character(),
    Components = numeric(),
    Standard_error = numeric(),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_default_variance_method <- function(x) {
  if (is.list(x)) {
    if (!is.null(x[["Asreml_model"]]) || !is.null(x[["asreml_model"]])) {
      return("ASReml_REML")
    }
    if (!is.null(x[["bayes_model"]]) || !is.null(x[["bayes_result"]])) {
      return("posterior_component_empirical")
    }
    if (!is.null(x[["gp_model"]]) || !is.null(x[["gp_result"]])) {
      return("kernel_prediction_component")
    }
    if (!is.null(x[["predicted_values"]]) || !is.null(x[["Predicted_value"]])) {
      return("prediction_component_empirical")
    }
  }
  "model_reported"
}

gp_hybrid_component_class <- function(component) {
  component <- tolower(as.character(component))
  out <- rep("hybrid_other", length(component))
  female <- grepl("female|dam|maternal", component)
  male <- !female & grepl("(^|[^[:alpha:]])male|sire|paternal", component)
  out[female & grepl("gca|additive|parent", component)] <- "female_gca"
  out[male & grepl("gca|additive|parent", component)] <- "male_gca"
  out[grepl("sca|specific|interaction|dominance|hybrid", component)] <- "sca"
  out[grepl("resid|residual|error|noise|units", component)] <- "residual"
  out[grepl("heritability", component)] <- "heritability"
  out[is.na(component) | !nzchar(component)] <- "hybrid_other"
  out
}

gp_format_hybrid_variance_component_columns <- function(out,
                                                        include_trait = NULL,
                                                        include_env = NULL) {
  if (!is.data.frame(out)) {
    return(gp_empty_hybrid_variance_components())
  }
  include_trait <- if (is.null(include_trait)) {
    "Trait" %in% names(out) && gp_vc_has_nonempty_labels(out[["Trait"]])
  } else {
    isTRUE(include_trait)
  }
  include_env <- if (is.null(include_env)) {
    "Env" %in% names(out) && gp_vc_has_nonempty_labels(out[["Env"]])
  } else {
    isTRUE(include_env)
  }
  if (isTRUE(include_trait) && !"Trait" %in% names(out)) {
    out[["Trait"]] <- NA_character_
  }
  if (isTRUE(include_env) && !"Env" %in% names(out)) {
    out[["Env"]] <- NA_character_
  }
  if (!"Component" %in% names(out)) {
    out[["Component"]] <- rownames(out)
  }
  if (!"Hybrid_component" %in% names(out)) {
    out[["Hybrid_component"]] <- gp_hybrid_component_class(out[["Component"]])
  }
  required <- c("Component", "Components", "Standard_error")
  for (nm in required) {
    if (!nm %in% names(out)) {
      out[[nm]] <- if (nm %in% c("Components", "Standard_error")) NA_real_ else NA_character_
    }
  }
  out <- out[, gp_hybrid_variance_component_columns(include_trait = include_trait, include_env = include_env), drop = FALSE]
  rownames(out) <- NULL
  out
}

gp_format_hybrid_variance_components <- function(x,
                                                 model = NULL,
                                                 response_family = "gaussian") {
  default_method <- gp_hybrid_default_variance_method(x)
  vc <- gp_format_variance_components(x, default_estimation_method = default_method)
  if (!is.data.frame(vc) || !nrow(vc)) {
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
    model_implied_engine <- is.list(x) && (
      !is.null(x[["bayes_model"]]) ||
      !is.null(x[["Asreml_model"]]) ||
      !is.null(x[["asreml_model"]]) ||
      !is.null(x[["gp_model"]]) ||
      !is.null(x[["gp_result"]])
    )
    if (identical(variance_role, "model_based") ||
        (identical(variance_role, "unknown") && isTRUE(model_implied_engine))) {
      return(gp_format_hybrid_variance_component_columns(data.frame(
        Hybrid_component = "hybrid_other",
        Component = "model_variance_components_not_returned",
        Components = NA_real_,
        Standard_error = NA_real_,
        stringsAsFactors = FALSE
      )))
    }
    if (!identical(variance_role, "predictive_only")) {
      return(gp_empty_hybrid_variance_components())
    }
    pred <- if (is.list(x)) x[["predicted_values"]] %||% x[["Predicted_value"]] else NULL
    pred <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    component_map <- c(
      female_gca_prediction_variance = "Female_additive_contribution",
      male_gca_prediction_variance = "Male_additive_contribution",
      sca_prediction_variance = "Hybrid_interaction_contribution",
      female_gca_effect_prediction_variance = "Female_GCA",
      male_gca_effect_prediction_variance = "Male_GCA",
      sca_effect_prediction_variance = "SCA_effect"
    )
    if (is.data.frame(pred) && nrow(pred)) {
      group_cols <- character()
      if ("Trait" %in% names(pred) &&
          any(!is.na(pred[["Trait"]]) & nzchar(as.character(pred[["Trait"]]))) &&
          gp_contract_multiple_nonempty_values(pred[["Trait"]])) {
        group_cols <- c(group_cols, "Trait")
      }
      if ("Env" %in% names(pred) &&
          any(!is.na(pred[["Env"]]) & nzchar(as.character(pred[["Env"]]))) &&
          gp_contract_multiple_nonempty_values(pred[["Env"]])) {
        group_cols <- c(group_cols, "Env")
      }
      pred_groups <- if (length(group_cols)) {
        split(pred, interaction(pred[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE), drop = TRUE)
      } else {
        list(pred)
      }
      pieces <- unlist(lapply(pred_groups, function(dat) {
        lapply(names(component_map), function(component) {
          col <- component_map[[component]]
          if (!col %in% names(dat)) {
            return(NULL)
          }
          vals <- suppressWarnings(as.numeric(dat[[col]]))
          vals <- vals[is.finite(vals)]
          if (length(vals) < 2L) {
            return(NULL)
          }
          out <- data.frame(
            Hybrid_component = gp_hybrid_component_class(component),
            Component = component,
            Components = stats::var(vals, na.rm = TRUE),
            Standard_error = NA_real_,
            Estimation_method = default_method,
            stringsAsFactors = FALSE
          )
          for (nm in group_cols) {
            out[[nm]] <- as.character(dat[[nm]][[1L]])
          }
          out
        })
      }), recursive = FALSE)
      pieces <- pieces[vapply(pieces, is.data.frame, logical(1L))]
      if (length(pieces)) {
        vc <- do.call(rbind, pieces)
        unavailable <- data.frame(
          Hybrid_component = "hybrid_other",
          Component = c(
            "genetic_variance_not_identifiable",
            "residual_variance_not_identifiable",
            "heritability_not_identifiable"
          ),
          Components = NA_real_,
          Standard_error = NA_real_,
          Estimation_method = "not_identifiable_from_predictive_model",
          stringsAsFactors = FALSE
        )
        for (nm in group_cols) {
          unavailable[[nm]] <- NA_character_
        }
        unavailable <- unavailable[, names(vc), drop = FALSE]
        vc <- rbind(vc, unavailable)
        return(gp_format_hybrid_variance_component_columns(
          vc,
          include_trait = "Trait" %in% group_cols,
          include_env = "Env" %in% group_cols
        ))
      }
    }
    return(gp_empty_hybrid_variance_components())
  }
  if (!"Component" %in% names(vc)) {
    vc[["Component"]] <- rownames(vc)
  }
  if (!"Estimation_method" %in% names(vc)) {
    vc[["Estimation_method"]] <- default_method
  }
  vc[["Component"]] <- as.character(vc[["Component"]])
  if ("Trait" %in% names(vc)) {
    vc[["Trait"]] <- as.character(vc[["Trait"]])
  }
  if ("Env" %in% names(vc)) {
    vc[["Env"]] <- as.character(vc[["Env"]])
  }
  vc[["Estimation_method"]] <- as.character(vc[["Estimation_method"]])
  missing_method <- is.na(vc[["Estimation_method"]]) | !nzchar(vc[["Estimation_method"]])
  vc[["Estimation_method"]][missing_method] <- default_method
  gp_format_hybrid_variance_component_columns(vc)
}

gp_hybrid_reference_component <- function(component) {
  component <- tolower(as.character(component))
  component <- ifelse(is.na(component), "", component)
  genetic_signal <- grepl(
    "gca|sca|additive|interaction|gxe|genetic|female|male|hybrid",
    component
  )
  nuisance <- grepl(
    "residual|error|noise|lambda|selected_lambda|heritability|correlation|covariance|intercept|fixed|prediction",
    component
  )
  genetic_signal & !nuisance
}

gp_hybrid_positive_variance_sum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x > 0]
  if (length(x)) sum(x, na.rm = TRUE) else NA_real_
}

gp_hybrid_add_reliability_reference_variance <- function(pred, vc) {
  out <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(out) || !nrow(out)) {
    return(pred)
  }

  ref <- if ("Reliability_reference_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Reliability_reference_variance"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  needs_ref <- !is.finite(ref) | ref <= 0
  if (!any(needs_ref, na.rm = TRUE)) {
    return(out)
  }

  vc <- tryCatch(as.data.frame(vc, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.data.frame(vc) && nrow(vc) &&
      all(c("Component", "Components") %in% names(vc))) {
    comp_ok <- gp_hybrid_reference_component(vc[["Component"]])
    vc <- vc[comp_ok & is.finite(suppressWarnings(as.numeric(vc[["Components"]]))) &
               suppressWarnings(as.numeric(vc[["Components"]])) > 0, , drop = FALSE]
    if (nrow(vc)) {
      vc[["Components"]] <- suppressWarnings(as.numeric(vc[["Components"]]))
      key_sets <- list(c("Trait", "Env"), "Trait", "Env")
      for (keys in key_sets) {
        if (!all(keys %in% names(out)) || !all(keys %in% names(vc))) {
          next
        }
        vc_key_data <- vc[, keys, drop = FALSE]
        vc_has_key <- apply(vc_key_data, 1L, function(z) all(!is.na(z) & nzchar(as.character(z))))
        if (!any(vc_has_key)) {
          next
        }
        pred_key <- do.call(paste, c(out[, keys, drop = FALSE], sep = "\r"))
        vc_key <- do.call(paste, c(vc[vc_has_key, keys, drop = FALSE], sep = "\r"))
        sums <- stats::aggregate(
          vc[["Components"]][vc_has_key],
          by = list(key = vc_key),
          FUN = sum,
          na.rm = TRUE
        )
        names(sums) <- c("key", "reference_variance")
        matched <- sums$reference_variance[match(pred_key, sums$key)]
        fill <- needs_ref & is.finite(matched) & matched > 0
        ref[fill] <- matched[fill]
        needs_ref <- !is.finite(ref) | ref <= 0
        if (!any(needs_ref, na.rm = TRUE)) {
          break
        }
      }

      global_ref <- gp_hybrid_positive_variance_sum(vc[["Components"]])
      fill <- needs_ref & is.finite(global_ref) & global_ref > 0
      ref[fill] <- global_ref
    }
  }

  needs_ref <- !is.finite(ref) | ref <= 0
  if (any(needs_ref, na.rm = TRUE) && "Observed_value" %in% names(out)) {
    obs_val <- suppressWarnings(as.numeric(out[["Observed_value"]]))
    fallback_ref <- gp_hybrid_positive_variance_sum(stats::var(obs_val, na.rm = TRUE))
    fill <- needs_ref & is.finite(fallback_ref) & fallback_ref > 0
    ref[fill] <- fallback_ref
  }
  needs_ref <- !is.finite(ref) | ref <= 0
  if (any(needs_ref, na.rm = TRUE) && "Predicted_value" %in% names(out)) {
    pred_val <- suppressWarnings(as.numeric(out[["Predicted_value"]]))
    train_label <- if ("Train_Test_Label" %in% names(out)) as.character(out[["Train_Test_Label"]]) else rep(NA_character_, nrow(out))
    train_row <- is.na(train_label) | train_label != "Test"
    fallback_ref <- gp_hybrid_positive_variance_sum(stats::var(pred_val[train_row], na.rm = TRUE))
    if (!is.finite(fallback_ref) || fallback_ref <= 0) {
      fallback_ref <- gp_hybrid_positive_variance_sum(stats::var(pred_val, na.rm = TRUE))
    }
    fill <- needs_ref & is.finite(fallback_ref) & fallback_ref > 0
    ref[fill] <- fallback_ref
  }

  out[["Reliability_reference_variance"]] <- ref
  out
}

gp_hybrid_add_true_prediction_uncertainty <- function(pred,
                                                      hybrid_id,
                                                      ml_dl_estimand = FALSE) {
  out <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(out) || !nrow(out)) {
    return(pred)
  }
  if (any(c("cv_scenario", "fold", "rep", "cv_role") %in% names(out))) {
    return(out)
  }
  if (!exists("gp_met_add_gaussian_uncertainty", mode = "function")) {
    return(out)
  }
  env_col <- gp_contract_first_name(
    names(out),
    c("Env", "env", "Environment", "environment", "Location", "location")
  )
  model_uncertainty <- ("Standard_error" %in% names(out) &&
                          any(is.finite(suppressWarnings(as.numeric(out[["Standard_error"]]))))) ||
    ("PEV" %in% names(out) && any(is.finite(suppressWarnings(as.numeric(out[["PEV"]]))))) ||
    ("Prediction_error_variance" %in% names(out) &&
       any(is.finite(suppressWarnings(as.numeric(out[["Prediction_error_variance"]])))))
  missing_unc_remarks <- if ("Uncertainty_remarks" %in% names(out)) {
    is.na(out[["Uncertainty_remarks"]]) | !nzchar(as.character(out[["Uncertainty_remarks"]]))
  } else {
    rep(TRUE, nrow(out))
  }
  out <- tryCatch(
    gp_met_add_gaussian_uncertainty(
      pred_obs = out,
      gen_name = hybrid_id,
      heter_groups = if (!is.na(env_col)) env_col else NULL,
      ml_dl_estimand = ml_dl_estimand
    ),
    error = function(e) out
  )
  if (isTRUE(model_uncertainty) && "Uncertainty_remarks" %in% names(out)) {
    out[["Uncertainty_remarks"]][missing_unc_remarks] <- "model_reported_prediction_error_variance"
  }
  out
}

gp_standardize_hybrid_model_result <- function(res_model_output,
                                               hybrid_id = NULL,
                                               female_parent = NULL,
                                               male_parent = NULL,
                                               model = NULL,
                                               response_family = "gaussian") {
  if (is.null(res_model_output) || !is.list(res_model_output)) {
    return(res_model_output)
  }
  pred <- res_model_output[["predicted_values"]] %||% res_model_output[["Predicted_value"]]
  if (!is.data.frame(pred) || !gp_is_hybrid_prediction_table(pred)) {
    return(res_model_output)
  }

  hybrid_id <- hybrid_id %||% gp_contract_first_name(
    names(pred),
    c("Hybrid_ID", "HybridID", "HybridCross", "Cross", "GID")
  )
  female_parent <- female_parent %||% gp_contract_first_name(
    names(pred),
    c("Female", "female_parent", "FemaleParent")
  )
  male_parent <- male_parent %||% gp_contract_first_name(
    names(pred),
    c("Male", "male_parent", "MaleParent")
  )
  if (length(hybrid_id) == 0L || is.na(hybrid_id) || !nzchar(hybrid_id)) {
    hybrid_id <- "Hybrid_ID"
  }
  if (length(female_parent) == 0L || is.na(female_parent) || !nzchar(female_parent)) {
    female_parent <- "Female"
  }
  if (length(male_parent) == 0L || is.na(male_parent) || !nzchar(male_parent)) {
    male_parent <- "Male"
  }

  vc <- gp_format_hybrid_variance_components(
    res_model_output,
    model = model,
    response_family = response_family
  )
  model_implied_engine <- !is.null(res_model_output[["bayes_model"]]) ||
    !is.null(res_model_output[["Asreml_model"]]) ||
    !is.null(res_model_output[["asreml_model"]]) ||
    !is.null(res_model_output[["gp_model"]]) ||
    !is.null(res_model_output[["gp_result"]])
  ml_dl_estimand <- !isTRUE(model_implied_engine) &&
    is.data.frame(vc) && "Component" %in% names(vc) && any(
      grepl("not_identifiable$", as.character(vc[["Component"]])),
      na.rm = TRUE
    )
  pred <- gp_hybrid_add_reliability_reference_variance(pred, vc)
  pred <- gp_hybrid_add_true_prediction_uncertainty(
    pred,
    hybrid_id = hybrid_id,
    ml_dl_estimand = ml_dl_estimand
  )

  pred_public <- gp_format_hybrid_prediction_table(
    pred,
    hybrid_id = hybrid_id,
    female_parent = female_parent,
    male_parent = male_parent
  )
  total_pred <- res_model_output[["Total_Predicted_value"]] %||% res_model_output[["total_predicted_values"]]
  if (is.data.frame(total_pred) && gp_is_hybrid_prediction_table(total_pred)) {
    total_pred <- gp_format_hybrid_prediction_table(
      total_pred,
      hybrid_id = hybrid_id,
      female_parent = female_parent,
      male_parent = male_parent,
      include_env = FALSE
    )
  } else {
    total_pred <- NULL
  }

  plots <- gp_public_diagnostic_plots(res_model_output)
  public <- list(
    model_parameters = gp_public_model_parameters(res_model_output[["model_parameters"]]),
    predicted_values = pred_public,
    diagnostic_plots = plots,
    variance_components = vc
  )
  if (!is.null(total_pred)) {
    public[["across_environment_predicted_values"]] <- total_pred
  }
  public <- gp_attach_correlation_matrices(public, res_model_output)
  predictive_se_names <- grep(
    "^Prediction_(error_)?(covariance|correlation)_(traits|environments)_SE$",
    names(res_model_output),
    value = TRUE
  )
  for (nm in predictive_se_names) {
    if (is.matrix(res_model_output[[nm]])) {
      public[[nm]] <- res_model_output[[nm]]
    }
  }
  for (nm in c("prediction_uncertainty_source", "predictive_covariance_basis")) {
    if (!is.null(res_model_output[[nm]])) {
      public[[nm]] <- res_model_output[[nm]]
    }
  }
  class(public) <- unique(c("predictpror_model_result", class(public)))
  validate_prediction_output(public, response_family = "gaussian", strict = TRUE)
  public
}
