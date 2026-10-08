gp_format_total_prediction_table <- function(x, gen_name = NULL, heter_groups = NULL, include_trait = NULL) {
  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      exists("gp_format_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(x))) {
    return(gp_format_classification_prediction_table(
      x = x,
      gen_name = gen_name,
      heter_groups = heter_groups,
      include_env = FALSE,
      include_trait = include_trait
    ))
  }
  gp_format_gaussian_prediction_table(
    x = x,
    gen_name = gen_name,
    heter_groups = heter_groups,
    include_env = FALSE,
    include_trait = include_trait
  )
}

gp_is_joint_multitrait_gaussian_prediction_table <- function(x) {
  df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df) ||
      !all(c("Trait", "Predicted_value", "Train_Test_Label") %in% names(df))) {
    return(FALSE)
  }
  if (exists("gp_is_classification_prediction_table", mode = "function") &&
      isTRUE(gp_is_classification_prediction_table(df))) {
    return(FALSE)
  }
  gp_contract_multiple_nonempty_values(df[["Trait"]])
}

gp_multitrait_has_cv_columns <- function(x) {
  any(c("cv_scenario", "fold", "rep", "cv_role", "Fold", "Rep") %in% names(x))
}

gp_multitrait_reference_component <- function(component) {
  component <- tolower(as.character(component))
  component[is.na(component)] <- ""
  signal <- grepl(
    "genetic|gmatrix|kernel|vm\\(|gid|geno|additive|gxe",
    component
  )
  nuisance <- grepl(
    "resid|residual|noise|heritability|correlation|covariance|environment",
    component
  )
  signal & !nuisance
}

gp_multitrait_positive_variance <- function(x, min_var = .Machine$double.eps) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) > 1L) {
    val <- stats::var(x, na.rm = TRUE)
    if (is.finite(val) && val > 0) {
      return(max(val, min_var))
    }
  }
  NA_real_
}

gp_multitrait_fill_reference_from_vc <- function(out, ref, vc, needs_ref) {
  vc <- tryCatch(as.data.frame(vc, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (!is.data.frame(vc) || !nrow(vc) || !all(c("Component", "Components") %in% names(vc))) {
    return(ref)
  }
  vc[["Components"]] <- suppressWarnings(as.numeric(vc[["Components"]]))
  keep <- gp_multitrait_reference_component(vc[["Component"]]) &
    is.finite(vc[["Components"]]) & vc[["Components"]] > 0
  vc <- vc[keep, , drop = FALSE]
  if (!nrow(vc)) {
    return(ref)
  }

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
      return(ref)
    }
  }

  global_ref <- sum(vc[["Components"]][is.finite(vc[["Components"]]) & vc[["Components"]] > 0], na.rm = TRUE)
  fill <- needs_ref & is.finite(global_ref) & global_ref > 0
  ref[fill] <- global_ref
  ref
}

gp_multitrait_fill_reference_empirically <- function(out,
                                                     ref,
                                                     needs_ref,
                                                     include_all_final_predictions = FALSE) {
  if (!"Predicted_value" %in% names(out)) {
    return(ref)
  }
  pred_val <- suppressWarnings(as.numeric(out[["Predicted_value"]]))
  obs_val <- if ("Observed_value" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Observed_value"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  train_label <- if ("Train_Test_Label" %in% names(out)) {
    as.character(out[["Train_Test_Label"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  train_row <- is.na(train_label) | train_label != "Test"
  if (!any(train_row, na.rm = TRUE)) {
    train_row <- rep(TRUE, nrow(out))
  }
  reference_row <- if (isTRUE(include_all_final_predictions)) {
    rep(TRUE, nrow(out))
  } else {
    train_row
  }

  scale_values <- c(pred_val, obs_val)
  scale_values <- scale_values[is.finite(scale_values)]
  scale_ref <- if (length(scale_values)) stats::median(abs(scale_values), na.rm = TRUE) else 1
  if (!is.finite(scale_ref) || scale_ref <= 0) {
    scale_ref <- 1
  }
  min_var <- max((scale_ref * 1e-4)^2, .Machine$double.eps)

  key_sets <- list(c("Trait", "Env"), "Trait", "Env")
  for (keys in key_sets) {
    if (!all(keys %in% names(out))) {
      next
    }
    key_data <- out[, keys, drop = FALSE]
    has_key <- apply(key_data, 1L, function(z) all(!is.na(z) & nzchar(as.character(z))))
    if (!any(has_key)) {
      next
    }
    row_key <- do.call(paste, c(key_data, sep = "\r"))
    keys_unique <- unique(row_key[has_key])
    ref_by_key <- stats::setNames(rep(NA_real_, length(keys_unique)), keys_unique)
    for (key in keys_unique) {
      idx <- row_key == key
      val <- gp_multitrait_positive_variance(pred_val[idx & reference_row], min_var = min_var)
      if (!is.finite(val)) {
        val <- gp_multitrait_positive_variance(pred_val[idx], min_var = min_var)
      }
      if (!is.finite(val)) {
        val <- gp_multitrait_positive_variance(obs_val[idx], min_var = min_var)
      }
      ref_by_key[[key]] <- val
    }
    matched <- ref_by_key[row_key]
    fill <- needs_ref & is.finite(matched) & matched > 0
    ref[fill] <- as.numeric(matched[fill])
    needs_ref <- !is.finite(ref) | ref <= 0
    if (!any(needs_ref, na.rm = TRUE)) {
      return(ref)
    }
  }

  fallback_ref <- gp_multitrait_positive_variance(pred_val[reference_row], min_var = min_var)
  if (!is.finite(fallback_ref)) {
    fallback_ref <- gp_multitrait_positive_variance(pred_val, min_var = min_var)
  }
  if (!is.finite(fallback_ref)) {
    fallback_ref <- gp_multitrait_positive_variance(obs_val, min_var = min_var)
  }
  fill <- needs_ref & is.finite(fallback_ref) & fallback_ref > 0
  ref[fill] <- fallback_ref
  ref
}

gp_multitrait_add_reliability_reference_variance <- function(pred,
                                                             vc = NULL,
                                                             ml_dl_estimand = FALSE) {
  out <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(out) || !nrow(out) || !"Trait" %in% names(out)) {
    return(pred)
  }

  ref <- if ("Reliability_reference_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Reliability_reference_variance"]]))
  } else if ("Genetic_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Genetic_variance"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  # For model-based multi-trait results, the labelled variance-component table
  # is authoritative.  Earlier generic formatting may have recycled the first
  # trait's variance across all rows; replace such provisional values with the
  # exact Trait/Env-matched model component before calculating reliability.
  if (!isTRUE(ml_dl_estimand)) {
    vc_ref <- gp_multitrait_fill_reference_from_vc(
      out,
      ref = rep(NA_real_, nrow(out)),
      vc = vc,
      needs_ref = rep(TRUE, nrow(out))
    )
    use_vc_ref <- is.finite(vc_ref) & vc_ref > 0
    ref[use_vc_ref] <- vc_ref[use_vc_ref]
  }
  needs_ref <- !is.finite(ref) | ref <= 0
  if (any(needs_ref, na.rm = TRUE)) {
    ref <- gp_multitrait_fill_reference_from_vc(out, ref, vc, needs_ref)
  }
  needs_ref <- !is.finite(ref) | ref <= 0
  if (any(needs_ref, na.rm = TRUE)) {
    ref <- gp_multitrait_fill_reference_empirically(
      out,
      ref,
      needs_ref,
      include_all_final_predictions = ml_dl_estimand
    )
  }

  out[["Reliability_reference_variance"]] <- ref
  out
}

gp_multitrait_add_true_prediction_uncertainty <- function(pred,
                                                          gen_name = NULL,
                                                          ml_dl_estimand = FALSE,
                                                          genetic_model = FALSE) {
  out <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(out) || !nrow(out) ||
      !gp_is_joint_multitrait_gaussian_prediction_table(out) ||
      gp_multitrait_has_cv_columns(out) ||
      !exists("gp_met_add_gaussian_uncertainty", mode = "function")) {
    return(pred)
  }

  model_uncertainty <- ("Standard_error" %in% names(out) &&
                          any(is.finite(suppressWarnings(as.numeric(out[["Standard_error"]]))))) ||
    ("PEV" %in% names(out) &&
       any(is.finite(suppressWarnings(as.numeric(out[["PEV"]]))))) ||
    ("Prediction_error_variance" %in% names(out) &&
       any(is.finite(suppressWarnings(as.numeric(out[["Prediction_error_variance"]])))))
  missing_unc_remarks <- if ("Uncertainty_remarks" %in% names(out)) {
    is.na(out[["Uncertainty_remarks"]]) | !nzchar(as.character(out[["Uncertainty_remarks"]]))
  } else {
    rep(TRUE, nrow(out))
  }

  had_env <- "Env" %in% names(out)
  group_col <- "Trait"
  temp_group_col <- NULL
  if (had_env && gp_contract_multiple_nonempty_values(out[["Env"]])) {
    temp_group_col <- ".PredictProR_trait_env_group"
    out[[temp_group_col]] <- paste(as.character(out[["Trait"]]), as.character(out[["Env"]]), sep = "\r")
    group_col <- temp_group_col
  }

  out <- tryCatch(
    gp_met_add_gaussian_uncertainty(
      pred_obs = out,
      gen_name = gen_name %||% "GID",
      heter_groups = group_col,
      ml_dl_estimand = ml_dl_estimand
    ),
    error = function(e) out
  )

  if (!is.null(temp_group_col) && temp_group_col %in% names(out)) {
    out[[temp_group_col]] <- NULL
  }
  if (!had_env && "Env" %in% names(out) &&
      identical(as.character(out[["Env"]]), as.character(out[["Trait"]]))) {
    out[["Env"]] <- NULL
  }
  if (isTRUE(model_uncertainty) && "Uncertainty_remarks" %in% names(out)) {
    out[["Uncertainty_remarks"]][missing_unc_remarks] <- "model_reported_prediction_error_variance"
  }
  # A genetic (GP / mixed) model without PEV, e.g. the Scalable-GBLUP
  # method-of-moments route: say so, instead of the ML/DL wording.
  if (isTRUE(genetic_model) && !isTRUE(ml_dl_estimand) && "Reliability_basis" %in% names(out)) {
    rel <- suppressWarnings(as.numeric(out[["Reliability"]] %||% rep(NA_real_, nrow(out))))
    no_pev <- grepl("^Reliability unavailable", as.character(out[["Reliability_basis"]])) & !is.finite(rel)
    if (any(no_pev)) {
      out[["Reliability_basis"]][no_pev] <- paste(
        "Reliability unavailable: the genetic model returned no prediction error",
        "variance (PEV), so 1 - PEV / genetic variance cannot be computed"
      )
      if ("Reliability_remarks" %in% names(out)) {
        out[["Reliability_remarks"]][no_pev] <- "PEV not returned by the model"
      }
    }
  }
  out
}

gp_standardize_multitrait_prediction_contract <- function(pred, vc = NULL, gen_name = NULL) {
  if (!gp_is_joint_multitrait_gaussian_prediction_table(pred)) {
    return(pred)
  }
  ml_dl_estimand <- is.data.frame(vc) && "Component" %in% names(vc) && any(
    grepl("not_identifiable$", as.character(vc[["Component"]])),
    na.rm = TRUE
  )
  pred <- gp_multitrait_add_reliability_reference_variance(
    pred,
    vc = vc,
    ml_dl_estimand = ml_dl_estimand
  )
  # a genetic model: its variance table reports a finite genetic variance
  genetic_model <- !isTRUE(ml_dl_estimand) && is.data.frame(vc) &&
    all(c("Component", "Components") %in% names(vc)) &&
    any(grepl("^genetic_variance", as.character(vc[["Component"]])) &
          is.finite(suppressWarnings(as.numeric(vc[["Components"]]))))
  gp_multitrait_add_true_prediction_uncertainty(
    pred,
    gen_name = gen_name,
    ml_dl_estimand = ml_dl_estimand,
    genetic_model = genetic_model
  )
}

gp_reconcile_multitrait_model_output <- function(model_output,
                                                 gen_name = NULL,
                                                 heter_groups = NULL,
                                                 model = NULL,
                                                 response_family = NULL) {
  if (!is.list(model_output)) {
    return(model_output)
  }
  pred <- model_output[["predicted_values"]] %||% model_output[["Predicted_value"]]
  if (!isTRUE(gp_is_joint_multitrait_gaussian_prediction_table(pred))) {
    return(model_output)
  }

  vc <- gp_public_variance_components(
    model_output,
    heter_groups = heter_groups,
    model = model,
    response_family = response_family
  )
  pred <- gp_standardize_multitrait_prediction_contract(
    pred,
    vc = vc,
    gen_name = gen_name
  )
  pred <- gp_public_prediction_table(
    pred,
    gen_name = gen_name,
    heter_groups = heter_groups
  )

  # Update only the public tables. Fitted model objects and route-specific
  # diagnostics remain in the bundle for normal package export/return.
  model_output[["predicted_values"]] <- pred
  if ("Predicted_value" %in% names(model_output)) {
    model_output[["Predicted_value"]] <- pred
  }
  model_output[["variance_components"]] <- vc
  if ("Variance_components" %in% names(model_output)) {
    model_output[["Variance_components"]] <- vc
  }
  model_output
}

gp_multitrait_across_environment_prediction <- function(pred,
                                                        gen_name = NULL,
                                                        heter_groups = NULL,
                                                        confidence_level = 0.95) {
  if (!gp_is_joint_multitrait_gaussian_prediction_table(pred)) {
    return(NULL)
  }
  pred <- as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE)
  # Phase 3.26 fix: respect user-supplied heter_groups -- the public
  # multi-trait formatter renames "Env" to the user's heter_groups col
  # (e.g. "YYY") so hardcoded "Env" returns NULL.
  env_candidates <- unique(stats::na.omit(c(
    as.character(heter_groups),
    "Env", "env", "Environment", "environment", "Loc", "loc", "Location", "location"
  )))
  env_resolved <- env_candidates[nzchar(env_candidates) & env_candidates %in% names(pred)][1L]
  if (is.na(env_resolved) || !nzchar(env_resolved)) {
    return(NULL)
  }
  if (length(unique(as.character(pred[[env_resolved]]))) <= 1L) {
    return(NULL)
  }
  gid_col <- gp_contract_first_name(
    names(pred),
    unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
  )
  if (is.na(gid_col)) {
    return(NULL)
  }

  key <- paste(as.character(pred[[gid_col]]), as.character(pred[["Trait"]]), sep = "\r")
  split_rows <- split(seq_len(nrow(pred)), key, drop = TRUE)
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  rows <- lapply(split_rows, function(idx) {
    gid <- as.character(pred[[gid_col]][idx][[1L]])
    trait <- as.character(pred[["Trait"]][idx][[1L]])
    n_env <- length(unique(as.character(pred[[env_resolved]][idx])))
    pred_vals <- suppressWarnings(as.numeric(pred[["Predicted_value"]][idx]))
    obs_vals <- if ("Observed_value" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Observed_value"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    pev_vals <- if ("PEV" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["PEV"]][idx]))
    } else if ("Prediction_error_variance" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Prediction_error_variance"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    rel_vals <- if ("Reliability" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Reliability"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    rel_variance_vals <- if ("Reliability_variance_input" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Reliability_variance_input"]][idx]))
    } else {
      pev_vals
    }
    ref_vals <- if ("Reliability_reference_variance" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Reliability_reference_variance"]][idx]))
    } else if ("Genetic_variance" %in% names(pred)) {
      suppressWarnings(as.numeric(pred[["Genetic_variance"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    basis_vals <- if ("Reliability_basis" %in% names(pred)) {
      as.character(pred[["Reliability_basis"]][idx])
    } else {
      rep(NA_character_, length(idx))
    }
    predictive_basis <- any(gp_contract_predictive_stability_basis(basis_vals), na.rm = TRUE)

    pred_mean <- mean(pred_vals, na.rm = TRUE)
    if (!is.finite(pred_mean)) pred_mean <- NA_real_
    obs_mean <- mean(obs_vals, na.rm = TRUE)
    if (!is.finite(obs_mean)) obs_mean <- NA_real_
    finite_pev <- pev_vals[is.finite(pev_vals) & pev_vals >= 0]
    pev <- if (length(finite_pev) && n_env > 0L) {
      sum(finite_pev, na.rm = TRUE) / (n_env^2)
    } else {
      NA_real_
    }
    se <- if (is.finite(pev)) sqrt(pmax(pev, 0)) else NA_real_
    ref <- mean(ref_vals[is.finite(ref_vals) & ref_vals > 0], na.rm = TRUE)
    if (!is.finite(ref) || ref <= 0) {
      ref <- NA_real_
    }
    finite_rel_variance <- rel_variance_vals[
      is.finite(rel_variance_vals) & rel_variance_vals >= 0
    ]
    rel_variance <- if (length(finite_rel_variance) && n_env > 0L) {
      sum(finite_rel_variance, na.rm = TRUE) / (n_env^2)
    } else {
      NA_real_
    }
    stability <- gp_ml_gaussian_stability(rel_variance, ref)
    predictive_precision <- gp_ml_predictive_reliability(pev, ref)
    rel <- if (predictive_basis) {
      predictive_precision$reliability
    } else if (is.finite(pev) && is.finite(ref) && ref > 0 &&
               exists("gp_met_empirical_reliability", mode = "function")) {
      gp_met_empirical_reliability(pev, ref)
    } else {
      rel <- mean(rel_vals[is.finite(rel_vals)], na.rm = TRUE)
      if (is.finite(rel)) pmax(0, pmin(1, rel)) else NA_real_
    }
    label <- if (any(as.character(pred[["Train_Test_Label"]][idx]) == "Test", na.rm = TRUE)) {
      "Test"
    } else {
      "Train"
    }
    lower <- if (is.finite(pred_mean) && is.finite(se)) pred_mean - z * se else NA_real_
    upper <- if (is.finite(pred_mean) && is.finite(se)) pred_mean + z * se else NA_real_
    data.frame(
      GID = gid,
      Trait = trait,
      Predicted_value = pred_mean,
      Train_Test_Label = label,
      Observed_value = obs_mean,
      Standard_error = se,
      PEV = pev,
      lower_bound = lower,
      upper_bound = upper,
      Uncertainty = if (is.finite(lower) && is.finite(upper)) upper - lower else NA_real_,
      Uncertainty_remarks = "across_environment_aggregate",
      Prediction_stability = if (predictive_basis) stability$stability else NA_real_,
      Prediction_stability_remarks = if (predictive_basis) {
        stability$remarks
      } else {
        NA_character_
      },
      Prediction_stability_reference_variance = if (predictive_basis) ref else NA_real_,
      Reliability_variance_input = if (predictive_basis) predictive_precision$variance_input else pev,
      Reliability_reference_variance = if (predictive_basis) predictive_precision$reference_variance else ref,
      Reliability = rel,
      Reliability_remarks = if (predictive_basis) {
        predictive_precision$remarks
      } else {
        gp_contract_reliability_remarks(rel)
      },
      Reliability_basis = if (predictive_basis) {
        predictive_precision$basis
      } else {
        NA_character_
      },
      Prediction_uncertainty_source = if (predictive_basis) {
        "across_environment_independence_approximation"
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[vapply(rows, is.data.frame, logical(1L))]
  if (!length(rows)) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gp_attach_correlation_matrices <- function(public_result, source_result) {
  mats <- gp_extract_correlation_matrices(source_result)
  for (nm in names(mats)) {
    mat <- mats[[nm]]
    # Reticulate can preserve pandas column labels while replacing its row
    # index with "1", "2", ... . Public square Cov/Cor matrices must name
    # the same variables on both axes, so recover those row labels before
    # the generic CSV writer sees the matrix.
    if (is.matrix(mat) && nrow(mat) == ncol(mat)) {
      rows <- rownames(mat)
      cols <- colnames(mat)
      rows_are_integer_strings <- !is.null(rows) && length(rows) == nrow(mat) &&
        all(nzchar(rows)) && all(grepl("^[0-9]+$", rows))
      cols_are_informative <- !is.null(cols) && length(cols) == ncol(mat) &&
        all(nzchar(cols)) && !all(grepl("^[0-9]+$", cols))
      if (rows_are_integer_strings && cols_are_informative) {
        rownames(mat) <- cols
      }
    }
    public_result[[nm]] <- mat
  }
  # Phase 3.27: env labels to relabel numpy-array matrices (e.g. Python's
  # `env_covariance_fa` is a plain np.ndarray with no axis labels). Pull
  # them from the public VC table's per-env Vg rows so the public Cov/Cor
  # tables carry env names on both axes consistently.
  env_labels_from_vc <- NULL
  vc_for_labels <- public_result[["variance_components"]]
  if (is.data.frame(vc_for_labels) && nrow(vc_for_labels) > 0L &&
      "Component" %in% names(vc_for_labels)) {
    hits <- grep("^genetic_variance_", as.character(vc_for_labels[["Component"]]))
    if (length(hits)) {
      env_labels_from_vc <- sub("^genetic_variance_", "",
                                as.character(vc_for_labels[["Component"]][hits]))
    }
  }
  if (exists("gp_extract_matrix_component_table", mode = "function")) {
    covariance_components <- gp_extract_matrix_component_table(
      source_result,
      value_type = "covariance",
      default_estimation_method = "model_reported",
      env_labels = env_labels_from_vc
    )
    correlation_components <- gp_extract_matrix_component_table(
      source_result,
      value_type = "correlation",
      default_estimation_method = "model_reported",
      env_labels = env_labels_from_vc
    )
    if (is.data.frame(covariance_components) && nrow(covariance_components)) {
      public_result[["covariance_components"]] <- covariance_components
    }
    if (is.data.frame(correlation_components) && nrow(correlation_components)) {
      public_result[["correlation_components"]] <- correlation_components
    }

    # Phase 3.27: ensure every GP MET model emits at least the CS-implied
    # env-by-env Cov_g / Cor_g, even single-kernel models (Kernel-GBLUP /
    # GP-GBLUP / Scalable-GBLUP) where Python correctly emits a near-zero
    # interaction matrix because the model carries no env-interaction
    # kernel. Compound symmetry implies Cov_g(i,j) = sqrt(Vg_i * Vg_j) and
    # Cor_g(i,j) = 1 between any two envs -- this is what the model
    # actually assumes; reporting it is more useful than reporting zeros.
    vc_for_cs <- public_result[["variance_components"]]
    if (is.data.frame(vc_for_cs) && nrow(vc_for_cs) > 0L) {
      cs_pair <- gp_cs_env_cov_corr_from_per_env_vg(vc_for_cs)
      if (!is.null(cs_pair) && !is.null(cs_pair$cov) && nrow(cs_pair$cov) >= 2L) {
        cov_tbl <- public_result[["covariance_components"]]
        if (gp_cov_table_is_degenerate(cov_tbl)) {
          public_result[["covariance_components"]] <- gp_matrix_to_component_table(
            cs_pair$cov,
            matrix_name = "Genetic_covariance_environments_CS",
            component = "environment_covariance",
            estimation_method = "compound_symmetry_from_per_env_Vg"
          )
        }
        cor_tbl <- public_result[["correlation_components"]]
        if (gp_cov_table_is_degenerate(cor_tbl)) {
          public_result[["correlation_components"]] <- gp_matrix_to_component_table(
            cs_pair$cor,
            matrix_name = "Genetic_correlation_environments_CS",
            component = "environment_correlation",
            estimation_method = "compound_symmetry_from_per_env_Vg"
          )
        }
      }
    }
  }
  public_result
}

# Phase 3.27: utility - is a cov/corr component table effectively empty?
# A degenerate table is one whose every sub-matrix carries no model info.
# Sub-matrices are identified by `Matrix`; each one is degenerate when:
#   - it has no rows, OR
#   - no finite values, OR
#   - all finite values are exactly zero (e.g. degenerate Cov), OR
#   - off-diagonal entries are all zero AND on-diagonal entries are all 1
#     (an identity matrix shape -- for Cor that's the no-GxE case, which
#     contradicts what compound-symmetry single-kernel models actually
#     assume; backfilling with CS Cor=J is the faithful interpretation).
# A table that mixes a zero matrix and an identity matrix is also
# degenerate (every constituent matrix carries no information), so the
# check is run per-Matrix and the overall table is degenerate iff every
# constituent is degenerate.
gp_cov_table_is_degenerate <- function(tbl) {
  if (!is.data.frame(tbl) || !nrow(tbl)) {
    return(TRUE)
  }
  if (!"Value" %in% names(tbl)) {
    return(TRUE)
  }
  vals <- suppressWarnings(as.numeric(tbl[["Value"]]))
  if (!any(is.finite(vals))) {
    return(TRUE)
  }
  rows <- as.character(tbl[["Row"]])
  cols <- as.character(tbl[["Column"]])
  group_keys <- if ("Matrix" %in% names(tbl)) as.character(tbl[["Matrix"]]) else rep("_", nrow(tbl))
  is_block_degenerate <- function(idx) {
    block_vals <- vals[idx]
    finite_block <- block_vals[is.finite(block_vals)]
    if (!length(finite_block)) {
      return(TRUE)
    }
    if (all(finite_block == 0)) {
      return(TRUE)
    }
    if (length(rows) == length(vals) && length(cols) == length(vals)) {
      diag_mask <- rows[idx] == cols[idx]
      on <- block_vals[diag_mask]
      off <- block_vals[!diag_mask]
      on <- on[is.finite(on)]
      off <- off[is.finite(off)]
      if (length(on) && length(off) &&
          all(off == 0) &&
          all(abs(on - 1) < 1e-12)) {
        return(TRUE)
      }
    }
    FALSE
  }
  for (key in unique(group_keys)) {
    idx <- which(group_keys == key)
    if (!is_block_degenerate(idx)) {
      return(FALSE)
    }
  }
  TRUE
}

# Phase 3.27: build the compound-symmetry env-by-env Cov_g + Cor_g from
# per-environment genetic-variance rows in a public VC table. Looks for
# rows named genetic_variance_<env> (the layout
# gp_variance_components_add_per_env appends across all GP MET models).
# Returns NULL when not enough per-env Vg rows are present.
gp_cs_env_cov_corr_from_per_env_vg <- function(vc) {
  if (!is.data.frame(vc) || !nrow(vc)) {
    return(NULL)
  }
  comp_col <- if ("Component" %in% names(vc)) "Component" else NULL
  if (is.null(comp_col)) {
    return(NULL)
  }
  comp <- as.character(vc[[comp_col]])
  hits <- grepl("^genetic_variance_", comp)
  if (sum(hits) < 2L) {
    return(NULL)
  }
  env_labels <- sub("^genetic_variance_", "", comp[hits])
  vg_vals <- suppressWarnings(as.numeric(vc[["Components"]][hits]))
  if (!any(is.finite(vg_vals) & vg_vals > 0)) {
    return(NULL)
  }
  vg_vals[!is.finite(vg_vals) | vg_vals < 0] <- NA_real_
  sd_vals <- sqrt(vg_vals)
  cov_mat <- outer(sd_vals, sd_vals)
  cor_mat <- matrix(1, nrow = length(env_labels), ncol = length(env_labels))
  dimnames(cov_mat) <- list(env_labels, env_labels)
  dimnames(cor_mat) <- list(env_labels, env_labels)
  list(cov = cov_mat, cor = cor_mat, env_labels = env_labels)
}
