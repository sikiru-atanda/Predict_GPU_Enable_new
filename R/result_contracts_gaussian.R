gp_gaussian_prediction_columns <- function(include_env = FALSE, include_trait = FALSE,
                                           gen_name = NULL, heter_groups = NULL) {
  gid_col <- if (!is.null(gen_name) && length(gen_name) == 1L && nzchar(as.character(gen_name))) {
    as.character(gen_name)
  } else {
    "GID"
  }
  env_col <- if (!is.null(heter_groups) && length(heter_groups) == 1L && nzchar(as.character(heter_groups))) {
    as.character(heter_groups)
  } else {
    "Env"
  }
  cols <- c(gid_col)
  if (isTRUE(include_env)) {
    cols <- c(cols, env_col)
  }
  if (isTRUE(include_trait)) {
    cols <- c(cols, "Trait")
  }
  c(
    cols,
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

gp_contract_predictive_stability_basis <- function(x) {
  grepl(
    "predictive_stability|Prediction_stability|not genetic reliability|not identifiable for ML/DL",
    ifelse(is.na(x), "", as.character(x)),
    ignore.case = TRUE
  )
}

gp_contract_nonempty_column <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  x <- as.character(x)
  any(!is.na(x) & nzchar(x))
}

gp_contract_multiple_nonempty_values <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  length(unique(x)) > 1L
}

gp_contract_first_name <- function(names_x, candidates) {
  hit <- candidates[candidates %in% names_x]
  if (length(hit)) hit[[1L]] else NA_character_
}

gp_contract_as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

gp_merge_gaussian_observed_values <- function(pred,
                                               pheno_data = NULL,
                                               response = NULL,
                                               gen_name = NULL,
                                               heter_groups = NULL) {
  pred_df <- tryCatch(
    as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  pheno_df <- tryCatch(
    as.data.frame(pheno_data, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  response <- as.character(response %||% character())
  response <- response[!is.na(response) & nzchar(response)]
  if (is.null(pred_df) || is.null(pheno_df) || !nrow(pred_df) ||
      !nrow(pheno_df) || length(response) != 1L ||
      !response %in% names(pheno_df)) {
    return(pred)
  }

  train_col <- gp_contract_first_name(
    names(pred_df),
    c("Train_Test_Label", "Train_Test", "set", "Set")
  )
  obs_col <- gp_contract_first_name(
    names(pred_df),
    c("Observed_value", "Observed", "observed", "y")
  )
  pred_id_col <- gp_contract_first_name(
    names(pred_df),
    unique(stats::na.omit(c(gen_name, "GID", "gid", "Name", "name", "ID", "id")))
  )
  pheno_id_col <- gp_contract_first_name(
    names(pheno_df),
    unique(stats::na.omit(c(gen_name, "GID", "gid", "Name", "name", "ID", "id")))
  )
  if (is.na(train_col) || is.na(pred_id_col) || is.na(pheno_id_col)) {
    return(pred_df)
  }

  pred_key_cols <- pred_id_col
  pheno_key_cols <- pheno_id_col
  pred_env_col <- gp_contract_first_name(
    names(pred_df),
    unique(stats::na.omit(c(heter_groups, "Env", "Environment", "env", "environment")))
  )
  pheno_env_col <- gp_contract_first_name(
    names(pheno_df),
    unique(stats::na.omit(c(heter_groups, "Env", "Environment", "env", "environment")))
  )
  if (!is.na(pred_env_col) && !is.na(pheno_env_col)) {
    pred_key_cols <- c(pred_key_cols, pred_env_col)
    pheno_key_cols <- c(pheno_key_cols, pheno_env_col)
  }

  make_key <- function(df, cols) {
    do.call(paste, c(lapply(cols, function(nm) as.character(df[[nm]])), sep = "\r"))
  }
  pheno_key <- make_key(pheno_df, pheno_key_cols)
  pheno_value <- gp_contract_as_numeric(pheno_df[[response]])
  finite_pheno <- is.finite(pheno_value) & !is.na(pheno_key) & nzchar(pheno_key)
  if (!any(finite_pheno)) {
    return(pred_df)
  }

  # Replicated Gaussian observations are reduced to the phenotype mean at the
  # prediction-table key. This matches the genotype or genotype-by-environment
  # estimand exposed by the public table while retaining only supplied values.
  observed_by_key <- stats::aggregate(
    pheno_value[finite_pheno],
    by = list(key = pheno_key[finite_pheno]),
    FUN = mean,
    na.rm = TRUE
  )
  names(observed_by_key) <- c("key", "observed")
  matched <- observed_by_key$observed[
    match(make_key(pred_df, pred_key_cols), observed_by_key$key)
  ]
  train <- !is.na(pred_df[[train_col]]) &
    tolower(trimws(as.character(pred_df[[train_col]]))) == "train"
  can_fill <- train & is.finite(matched)
  if (is.na(obs_col)) {
    pred_df[["Observed_value"]] <- rep(NA_real_, nrow(pred_df))
    obs_col <- "Observed_value"
  }
  current <- gp_contract_as_numeric(pred_df[[obs_col]])
  fill <- can_fill & !is.finite(current)
  current[fill] <- matched[fill]
  pred_df[[obs_col]] <- current
  pred_df
}

gp_contract_first_numeric_name <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) {
    return(NA_character_)
  }
  finite_hit <- hit[vapply(
    hit,
    function(nm) any(is.finite(gp_contract_as_numeric(df[[nm]]))),
    logical(1L)
  )]
  if (length(finite_hit)) finite_hit[[1L]] else hit[[1L]]
}

gp_contract_first_nonempty_name <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) {
    return(NA_character_)
  }
  nonempty_hit <- hit[vapply(
    hit,
    function(nm) {
      val <- as.character(df[[nm]])
      any(!is.na(val) & nzchar(trimws(val)))
    },
    logical(1L)
  )]
  if (length(nonempty_hit)) nonempty_hit[[1L]] else hit[[1L]]
}

gp_contract_fill_numeric <- function(n) {
  rep(NA_real_, n)
}

gp_contract_reliability_remarks <- function(reliability,
                                            high_reliability_thres = 0.9,
                                            low_reliability_thres = 0.5) {
  reliability <- gp_contract_as_numeric(reliability)
  ifelse(
    !is.finite(reliability),
    NA_character_,
    ifelse(
      reliability >= high_reliability_thres,
      "Reliable",
      ifelse(reliability >= low_reliability_thres, "Acceptable", "Unreliable")
    )
  )
}

gp_contract_positive_variance <- function(x, min_var = 0) {
  x <- gp_contract_as_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  v <- stats::var(x, na.rm = TRUE)
  if (is.finite(v) && v > min_var) v else NA_real_
}

gp_contract_recompute_flat_reliability <- function(out, df = out) {
  if (!is.data.frame(out) || !nrow(out) ||
      !"Reliability" %in% names(out) || !"PEV" %in% names(out)) {
    return(out)
  }
  names_x <- names(df)
  reliability_basis <- if ("Reliability_basis" %in% names(out)) {
    as.character(out[["Reliability_basis"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  protected_nonidentifiable <- gp_contract_predictive_stability_basis(reliability_basis)
  rel_ref_col <- gp_contract_first_numeric_name(
    df,
    c(
      "Reliability_reference_variance", "Genetic_variance", "genetic_variance",
      "Reference_variance", "reference_variance", "Prediction_reference_variance"
    )
  )
  ref_var <- if (!is.na(rel_ref_col)) {
    gp_contract_as_numeric(df[[rel_ref_col]])
  } else if ("Reliability_reference_variance" %in% names(out)) {
    gp_contract_as_numeric(out[["Reliability_reference_variance"]])
  } else if ("Genetic_variance" %in% names(out)) {
    gp_contract_as_numeric(out[["Genetic_variance"]])
  } else {
    rep(NA_real_, nrow(out))
  }
  if (length(ref_var) != nrow(out)) {
    ref_var <- rep(NA_real_, nrow(out))
  }

  rel_group_col <- gp_contract_first_name(
    names_x,
    c("Env", "env", "Environment", "environment", "Location", "location")
  )
  rel_group <- if (!is.na(rel_group_col)) {
    as.character(df[[rel_group_col]])
  } else if ("Env" %in% names(out)) {
    as.character(out[["Env"]])
  } else {
    rep("All", nrow(out))
  }
  rel_group[is.na(rel_group) | !nzchar(rel_group)] <- "All"

  pred <- if ("Predicted_value" %in% names(out)) {
    gp_contract_as_numeric(out[["Predicted_value"]])
  } else {
    rep(NA_real_, nrow(out))
  }
  obs <- if ("Observed_value" %in% names(out)) {
    gp_contract_as_numeric(out[["Observed_value"]])
  } else {
    rep(NA_real_, nrow(out))
  }
  train_label <- if ("Train_Test_Label" %in% names(out)) {
    as.character(out[["Train_Test_Label"]])
  } else {
    rep(NA_character_, nrow(out))
  }

  recompute_rel <- !is.finite(gp_contract_as_numeric(out[["Reliability"]])) &
    !protected_nonidentifiable
  pev <- gp_contract_as_numeric(out[["PEV"]])
  mark_flat_reliability_group <- function(group_idx) {
    if (length(group_idx) < 3L) {
      return(invisible(NULL))
    }
    finite_rel <- gp_contract_as_numeric(out[["Reliability"]][group_idx])
    finite_rel <- finite_rel[is.finite(finite_rel)]
    finite_pev <- pev[group_idx]
    finite_pev <- finite_pev[is.finite(finite_pev)]
    flat_rel_with_varying_pev <- length(finite_rel) >= 3L &&
      length(unique(round(finite_rel, 12L))) <= 1L &&
      length(unique(round(finite_pev, 12L))) > 1L
    if (!isTRUE(flat_rel_with_varying_pev)) {
      return(invisible(NULL))
    }

    group_ref <- ref_var[group_idx]
    if (!any(is.finite(group_ref) & group_ref > 0, na.rm = TRUE)) {
      train_idx <- group_idx[is.na(train_label[group_idx]) | train_label[group_idx] != "Test"]
      if (!length(train_idx)) {
        train_idx <- group_idx
      }
      scale_values <- c(pred[group_idx], obs[group_idx])
      scale_values <- scale_values[is.finite(scale_values)]
      scale_ref <- if (length(scale_values)) stats::median(abs(scale_values), na.rm = TRUE) else 1
      if (!is.finite(scale_ref) || scale_ref <= 0) {
        scale_ref <- 1
      }
      min_var <- max((scale_ref * 1e-4)^2, .Machine$double.eps)
      fallback_ref <- gp_contract_positive_variance(pred[train_idx], min_var = min_var)
      if (!is.finite(fallback_ref)) {
        fallback_ref <- gp_contract_positive_variance(pred[group_idx], min_var = min_var)
      }
      if (!is.finite(fallback_ref)) {
        fallback_ref <- gp_contract_positive_variance(obs[group_idx], min_var = min_var)
      }
      if (!is.finite(fallback_ref) && length(finite_pev)) {
        fallback_ref <- max(mean(finite_pev, na.rm = TRUE), min_var)
      }
      ref_var[group_idx] <<- fallback_ref
    }
    recompute_rel[group_idx] <<- TRUE
    invisible(NULL)
  }

  for (lev in unique(rel_group)) {
    mark_flat_reliability_group(which(rel_group == lev))
  }
  label_group <- paste(rel_group, train_label, sep = "\r")
  label_group[is.na(label_group) | !nzchar(label_group)] <- rel_group[is.na(label_group) | !nzchar(label_group)]
  for (lev in unique(label_group)) {
    mark_flat_reliability_group(which(label_group == lev))
  }

  can_rel <- recompute_rel & !protected_nonidentifiable & is.finite(pev) & pev >= 0 &
    is.finite(ref_var) & ref_var > 0
  if (any(can_rel, na.rm = TRUE)) {
    rel_vals <- gp_met_empirical_reliability(
      pev[can_rel],
      ref_var[can_rel]
    )
    out[["Reliability"]][can_rel] <- rel_vals
    out <- gp_contract_mark_ratio_fallback(out, which(can_rel), attr(rel_vals, "ratio_fallback"))
    if ("Reliability_reference_variance" %in% names(out)) {
      out[["Reliability_reference_variance"]][can_rel] <- ref_var[can_rel]
    }
    if ("Genetic_variance" %in% names(out)) {
      out[["Genetic_variance"]][can_rel] <- ref_var[can_rel]
    }
  }
  out
}

gp_contract_mark_ratio_fallback <- function(out, idx, fb) {
  # Append a transparent note to Reliability_remarks for rows where the
  # reliability ratio-form fallback replaced the standard 1 - PEV/sigma2_g, so
  # the mixed definition is never silent. `idx` are the row positions passed to
  # gp_met_empirical_reliability; `fb` is its `ratio_fallback` attribute.
  if (is.null(fb) || !length(idx)) {
    return(out)
  }
  hit <- idx[which(as.logical(fb))]
  if (!length(hit)) {
    return(out)
  }
  if (!"Reliability_remarks" %in% names(out)) {
    out[["Reliability_remarks"]] <- rep(NA_character_, nrow(out))
  }
  note <- "ratio-form fallback used"
  cur <- as.character(out[["Reliability_remarks"]][hit])
  cur[is.na(cur)] <- ""
  already <- grepl(note, cur, fixed = TRUE)
  sep <- ifelse(nzchar(cur) & !already, "; ", "")
  out[["Reliability_remarks"]][hit] <- ifelse(already, cur, paste0(cur, sep, note))
  out
}

gp_format_gaussian_prediction_table <- function(x,
                                                gen_name = NULL,
                                                heter_groups = NULL,
                                                include_env = NULL,
                                                include_trait = NULL,
                                                confidence_level = 0.95) {
  # `heter_groups` (when supplied by the caller from the user's input) is
  # treated as the AUTHORITATIVE env-column name. Hardcoded fallbacks are
  # kept only for legacy paths that don't thread it through.
  if (is.null(x)) {
    return(data.frame(
      GID = character(),
      Predicted_value = numeric(),
      Train_Test_Label = character(),
      Observed_value = numeric(),
      Standard_error = numeric(),
      PEV = numeric(),
      Prediction_uncertainty_source = character(),
      lower_bound = numeric(),
      upper_bound = numeric(),
      Uncertainty = numeric(),
      Uncertainty_remarks = character(),
      Reliability = numeric(),
      Reliability_remarks = character(),
      Reliability_variance_input = numeric(),
      Reliability_reference_variance = numeric(),
      Reliability_basis = character(),
      stringsAsFactors = FALSE
    ))
  }
  df <- tryCatch(
    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    return(x)
  }

  n <- nrow(df)
  names_x <- names(df)
  gid_col <- gp_contract_first_name(
    names_x,
    unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
  )
  pred_col <- gp_contract_first_numeric_name(
    df,
    c("Predicted_value", "Prediction", "prediction", "yhat", "predicted.value",
      "predicted_value", "BLUP", "Estimated_breeding_value")
  )
  train_col <- gp_contract_first_name(names_x, c("Train_Test_Label", "Train_Test", "set", "Set"))
  obs_col <- gp_contract_first_name(names_x, c("Observed_value", "Observed", "observed", "y", "response"))
  se_col <- gp_contract_first_numeric_name(
    df,
    c("Standard_error", "std.error", "std_error", "SE", "Prediction_SE",
      "Prediction_SE_latent", "Prediction_SE_observed")
  )
  pev_col <- gp_contract_first_numeric_name(
    df,
    c("PEV", "Prediction_error_variance", "Prediction_Var",
      "Prediction_Var_latent", "Prediction_Var_observed", "prediction_error_var")
  )
  pev_basis_col <- gp_contract_first_nonempty_name(df, c("PEV_basis", "Prediction_error_estimand"))
  interval_method_col <- gp_contract_first_nonempty_name(df, c("Prediction_interval_method"))
  interval_coverage_col <- gp_contract_first_numeric_name(df, c("Prediction_interval_nominal_coverage"))
  interval_n_col <- gp_contract_first_numeric_name(df, c("Prediction_interval_calibration_n"))
  stability_col <- gp_contract_first_numeric_name(df, c("Prediction_stability"))
  stability_remarks_col <- gp_contract_first_nonempty_name(df, c("Prediction_stability_remarks"))
  stability_ref_col <- gp_contract_first_numeric_name(df, c("Prediction_stability_reference_variance"))
  lower_col <- gp_contract_first_numeric_name(df, c("lower_bound", "lower", "lower_conf", "Lower"))
  upper_col <- gp_contract_first_numeric_name(df, c("upper_bound", "upper", "upper_conf", "Upper"))
  unc_col <- gp_contract_first_numeric_name(df, c("Uncertainty", "MPIW", "interval_width"))
  unc_rem_col <- gp_contract_first_nonempty_name(
    df,
    c("Uncertainty_remarks", "reliability_remarks", "Prediction_stability_remarks",
      "Prediction_confidence_remarks")
  )
  rel_col <- gp_contract_first_numeric_name(
    df,
    c("Reliability")
  )
  rel_ref_col <- gp_contract_first_numeric_name(
    df,
    c(
      "Reliability_reference_variance", "Genetic_variance", "genetic_variance",
      "Reference_variance", "reference_variance", "Prediction_reference_variance"
    )
  )
  rel_variance_col <- gp_contract_first_numeric_name(
    df,
    c("Reliability_variance_input", "Prediction_stability_variance_input")
  )
  rel_rem_col <- gp_contract_first_nonempty_name(
    df,
    c("Reliability_remarks")
  )
  reliability_basis_col <- gp_contract_first_nonempty_name(
    df,
    c("Reliability_basis", "reliability_basis")
  )
  uncertainty_source_col <- gp_contract_first_nonempty_name(
    df,
    c(
      "Prediction_uncertainty_source", "prediction_uncertainty_source",
      "Uncertainty_source", "uncertainty_source"
    )
  )

  # Honor the user-supplied gen_name for the output column name (e.g. "Line"
  # or "Sample_ID"). Fallback is "GID" only when gen_name is not provided.
  gid_out_name <- if (!is.null(gen_name) &&
                      length(gen_name) == 1L &&
                      nzchar(as.character(gen_name))) {
    as.character(gen_name)
  } else {
    "GID"
  }
  out <- data.frame(
    .gid_placeholder = if (!is.na(gid_col)) as.character(df[[gid_col]]) else rep(NA_character_, n),
    Predicted_value = if (!is.na(pred_col)) gp_contract_as_numeric(df[[pred_col]]) else gp_contract_fill_numeric(n),
    Train_Test_Label = if (!is.na(train_col)) as.character(df[[train_col]]) else rep(NA_character_, n),
    Observed_value = if (!is.na(obs_col)) gp_contract_as_numeric(df[[obs_col]]) else gp_contract_fill_numeric(n),
    Standard_error = if (!is.na(se_col)) gp_contract_as_numeric(df[[se_col]]) else gp_contract_fill_numeric(n),
    PEV = if (!is.na(pev_col)) gp_contract_as_numeric(df[[pev_col]]) else gp_contract_fill_numeric(n),
    PEV_basis = if (!is.na(pev_basis_col)) as.character(df[[pev_basis_col]]) else rep(NA_character_, n),
    Prediction_uncertainty_source = if (!is.na(uncertainty_source_col)) as.character(df[[uncertainty_source_col]]) else rep(NA_character_, n),
    Prediction_interval_method = if (!is.na(interval_method_col)) as.character(df[[interval_method_col]]) else rep(NA_character_, n),
    Prediction_interval_nominal_coverage = if (!is.na(interval_coverage_col)) gp_contract_as_numeric(df[[interval_coverage_col]]) else gp_contract_fill_numeric(n),
    Prediction_interval_calibration_n = if (!is.na(interval_n_col)) gp_contract_as_numeric(df[[interval_n_col]]) else gp_contract_fill_numeric(n),
    lower_bound = if (!is.na(lower_col)) gp_contract_as_numeric(df[[lower_col]]) else gp_contract_fill_numeric(n),
    upper_bound = if (!is.na(upper_col)) gp_contract_as_numeric(df[[upper_col]]) else gp_contract_fill_numeric(n),
    Uncertainty = if (!is.na(unc_col)) gp_contract_as_numeric(df[[unc_col]]) else gp_contract_fill_numeric(n),
    Uncertainty_remarks = if (!is.na(unc_rem_col)) as.character(df[[unc_rem_col]]) else rep(NA_character_, n),
    Prediction_stability = if (!is.na(stability_col)) gp_contract_as_numeric(df[[stability_col]]) else gp_contract_fill_numeric(n),
    Prediction_stability_remarks = if (!is.na(stability_remarks_col)) as.character(df[[stability_remarks_col]]) else rep(NA_character_, n),
    Prediction_stability_reference_variance = if (!is.na(stability_ref_col)) gp_contract_as_numeric(df[[stability_ref_col]]) else gp_contract_fill_numeric(n),
    Reliability = if (!is.na(rel_col)) gp_contract_as_numeric(df[[rel_col]]) else gp_contract_fill_numeric(n),
    Reliability_remarks = if (!is.na(rel_rem_col)) as.character(df[[rel_rem_col]]) else rep(NA_character_, n),
    Reliability_variance_input = if (!is.na(rel_variance_col)) gp_contract_as_numeric(df[[rel_variance_col]]) else gp_contract_fill_numeric(n),
    Reliability_reference_variance = if (!is.na(rel_ref_col)) gp_contract_as_numeric(df[[rel_ref_col]]) else gp_contract_fill_numeric(n),
    Reliability_basis = if (!is.na(reliability_basis_col)) as.character(df[[reliability_basis_col]]) else rep(NA_character_, n),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out)[names(out) == ".gid_placeholder"] <- gid_out_name

  pev_missing_before_derivation <- !is.finite(out[["PEV"]])
  if (all(pev_missing_before_derivation) && any(is.finite(out[["Standard_error"]]))) {
    out[["PEV"]] <- out[["Standard_error"]]^2
  }
  if (all(!is.finite(out[["Standard_error"]])) && any(is.finite(out[["PEV"]]))) {
    out[["Standard_error"]] <- sqrt(pmax(out[["PEV"]], 0))
  }

  predictive_alias <- !is.finite(out[["Reliability"]]) &
    is.finite(out[["Prediction_stability"]]) &
    is.finite(out[["PEV"]]) & out[["PEV"]] >= 0
  if (any(predictive_alias, na.rm = TRUE)) {
    out[["Reliability"]][predictive_alias] <- out[["Prediction_stability"]][predictive_alias]
    missing_rel_rem <- predictive_alias &
      (is.na(out[["Reliability_remarks"]]) | !nzchar(out[["Reliability_remarks"]]))
    out[["Reliability_remarks"]][missing_rel_rem] <-
      out[["Prediction_stability_remarks"]][missing_rel_rem]
    missing_ref <- predictive_alias & !is.finite(out[["Reliability_reference_variance"]]) &
      is.finite(out[["Prediction_stability_reference_variance"]])
    out[["Reliability_reference_variance"]][missing_ref] <-
      out[["Prediction_stability_reference_variance"]][missing_ref]
    missing_variance_input <- predictive_alias & !is.finite(out[["Reliability_variance_input"]]) &
      is.finite(out[["PEV"]]) & out[["PEV"]] >= 0
    out[["Reliability_variance_input"]][missing_variance_input] <-
      out[["PEV"]][missing_variance_input]
    missing_basis <- predictive_alias & (
      is.na(out[["Reliability_basis"]]) |
        !nzchar(out[["Reliability_basis"]]) |
        grepl("not identifiable for ML/DL", out[["Reliability_basis"]], fixed = TRUE)
    )
    out[["Reliability_basis"]][missing_basis] <- paste0(
      "Reliability equals the supplied Prediction_stability value; ",
      "compatibility alias for plotting; original variance inputs unavailable; ",
      "descriptive, not genetic reliability"
    )
  }

  if (!is.na(rel_ref_col)) {
    ref_var <- gp_contract_as_numeric(df[[rel_ref_col]])
    protected_nonidentifiable <- gp_contract_predictive_stability_basis(
      out[["Reliability_basis"]]
    )
    recompute_rel <- !is.finite(out[["Reliability"]]) & !protected_nonidentifiable
    rel_group_col <- gp_contract_first_name(
      names_x,
      c("Env", "env", "Environment", "environment", "Location", "location")
    )
    rel_group <- if (!is.na(rel_group_col)) {
      as.character(df[[rel_group_col]])
    } else {
      rep("All", n)
    }
    rel_group[is.na(rel_group) | !nzchar(rel_group)] <- "All"
    for (lev in unique(rel_group)) {
      group_idx <- which(rel_group == lev)
      finite_rel <- out[["Reliability"]][group_idx]
      finite_rel <- finite_rel[is.finite(finite_rel)]
      finite_pev <- out[["PEV"]][group_idx]
      finite_pev <- finite_pev[is.finite(finite_pev)]
      if (length(finite_rel) >= 3L &&
          length(unique(round(finite_rel, 12L))) <= 1L &&
          length(unique(round(finite_pev, 12L))) > 1L) {
        recompute_rel[group_idx] <- TRUE
      }
    }
    can_rel <- recompute_rel & !protected_nonidentifiable &
      is.finite(out[["PEV"]]) & out[["PEV"]] >= 0 &
      is.finite(ref_var) & ref_var > 0
    if (any(can_rel, na.rm = TRUE)) {
      rel_vals <- gp_met_empirical_reliability(
        out[["PEV"]][can_rel],
        ref_var[can_rel]
      )
      out[["Reliability"]][can_rel] <- rel_vals
      out <- gp_contract_mark_ratio_fallback(out, which(can_rel), attr(rel_vals, "ratio_fallback"))
    }
  }
  out <- gp_contract_recompute_flat_reliability(out, df = df)

  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  missing_bounds <- !is.finite(out[["lower_bound"]]) | !is.finite(out[["upper_bound"]])
  can_bound <- is.finite(out[["Predicted_value"]]) & is.finite(out[["Standard_error"]])
  idx <- missing_bounds & can_bound
  out[["lower_bound"]][idx] <- out[["Predicted_value"]][idx] - z * out[["Standard_error"]][idx]
  out[["upper_bound"]][idx] <- out[["Predicted_value"]][idx] + z * out[["Standard_error"]][idx]

  # A finite model-reported SE/PEV and Gaussian interval must carry its
  # estimand and interval provenance in the same public table.  The formatter
  # derives missing bounds by a normal approximation above, so leaving these
  # metadata fields blank makes an otherwise valid interval uninterpretable.
  # Existing backend-specific metadata (for example posterior quantiles or
  # held-out calibration) is always preserved.
  model_uncertainty <- is.finite(out[["Standard_error"]]) &
    out[["Standard_error"]] >= 0 &
    is.finite(out[["PEV"]]) & out[["PEV"]] >= 0
  blank_character <- function(value) {
    is.na(value) | !nzchar(trimws(ifelse(is.na(value), "", as.character(value))))
  }
  missing_pev_basis <- blank_character(out[["PEV_basis"]])
  from_se <- model_uncertainty & pev_missing_before_derivation
  out[["PEV_basis"]][missing_pev_basis & model_uncertainty] <-
    "model-reported prediction-error variance"
  out[["PEV_basis"]][missing_pev_basis & from_se] <-
    "square of model-reported prediction standard error"

  missing_source <- blank_character(out[["Prediction_uncertainty_source"]])
  out[["Prediction_uncertainty_source"]][missing_source & model_uncertainty] <-
    "model_reported_prediction_uncertainty"

  expected_lower <- out[["Predicted_value"]] - z * out[["Standard_error"]]
  expected_upper <- out[["Predicted_value"]] + z * out[["Standard_error"]]
  interval_tolerance <- 1e-8 * pmax(
    1,
    abs(expected_lower), abs(expected_upper),
    abs(out[["lower_bound"]]), abs(out[["upper_bound"]]),
    na.rm = TRUE
  )
  normal_interval <- model_uncertainty &
    is.finite(out[["lower_bound"]]) & is.finite(out[["upper_bound"]]) &
    abs(out[["lower_bound"]] - expected_lower) <= interval_tolerance &
    abs(out[["upper_bound"]] - expected_upper) <= interval_tolerance
  missing_method <- blank_character(out[["Prediction_interval_method"]])
  out[["Prediction_interval_method"]][missing_method & normal_interval] <-
    "Gaussian approximation from model-reported standard error"
  missing_coverage <- !is.finite(out[["Prediction_interval_nominal_coverage"]])
  out[["Prediction_interval_nominal_coverage"]][missing_coverage & normal_interval] <-
    confidence_level

  missing_unc <- !is.finite(out[["Uncertainty"]])
  can_unc <- is.finite(out[["lower_bound"]]) & is.finite(out[["upper_bound"]])
  out[["Uncertainty"]][missing_unc & can_unc] <-
    out[["upper_bound"]][missing_unc & can_unc] - out[["lower_bound"]][missing_unc & can_unc]

  reliability_remarks <- as.character(out[["Reliability_remarks"]])
  missing_rel_remarks <- is.na(reliability_remarks) |
    !nzchar(trimws(ifelse(is.na(reliability_remarks), "", reliability_remarks)))
  derived_rel_remarks <- gp_contract_reliability_remarks(out[["Reliability"]])
  fill_rel_remarks <- missing_rel_remarks & !is.na(derived_rel_remarks)
  reliability_remarks[fill_rel_remarks] <- derived_rel_remarks[fill_rel_remarks]
  out[["Reliability_remarks"]] <- reliability_remarks

  # Resolve the env column. Authoritative: the user-supplied `heter_groups`
  # column from model_execute() (can be ANY name -- "Loc", "yearxLoc", etc.).
  # Fallbacks (common names) only apply when heter_groups was not threaded.
  env_col <- if (!is.null(heter_groups) &&
                 length(heter_groups) == 1L &&
                 nzchar(as.character(heter_groups)) &&
                 as.character(heter_groups) %in% names_x) {
    as.character(heter_groups)
  } else {
    gp_contract_first_name(names_x, c("Env", "env", "Environment", "environment", "Location", "location"))
  }
  trait_col <- gp_contract_first_name(names_x, c("Trait", "trait", "Response", "response_name"))
  if (is.null(include_env)) {
    include_env <- !is.na(env_col) && gp_contract_multiple_nonempty_values(df[[env_col]])
  }
  if (is.null(include_trait)) {
    include_trait <- !is.na(trait_col) && gp_contract_multiple_nonempty_values(df[[trait_col]])
  }
  # Preserve the user's chosen env-column name in the output (if heter_groups
  # was supplied), otherwise use the legacy "Env" alias for back-compatibility.
  env_out_name <- if (!is.null(heter_groups) &&
                      length(heter_groups) == 1L &&
                      nzchar(as.character(heter_groups))) {
    as.character(heter_groups)
  } else {
    "Env"
  }
  if (isTRUE(include_env)) {
    out[[env_out_name]] <- if (!is.na(env_col)) as.character(df[[env_col]]) else rep(NA_character_, n)
  }
  if (isTRUE(include_trait)) {
    out[["Trait"]] <- if (!is.na(trait_col)) as.character(df[[trait_col]]) else rep(NA_character_, n)
  }

  keep_cols <- gp_gaussian_prediction_columns(
    include_env = include_env, include_trait = include_trait,
    gen_name = gid_out_name, heter_groups = env_out_name
  )
  out[, keep_cols, drop = FALSE]
}

gp_is_public_gaussian_table <- function(x) {
  is.data.frame(x) &&
    all(c(
      "Predicted_value", "Train_Test_Label", "Observed_value",
      "Standard_error", "PEV", "lower_bound", "upper_bound",
      "Uncertainty", "Uncertainty_remarks", "Reliability", "Reliability_remarks"
    ) %in% names(x))
}

gp_gaussian_diagnostic_plots <- function(pred) {
  if (!gp_is_public_gaussian_table(pred)) {
    return(NULL)
  }

  obs <- pred[is.finite(pred[["Observed_value"]]) & is.finite(pred[["Predicted_value"]]), , drop = FALSE]
  p_obs <- NULL
  if (nrow(obs)) {
    p_obs <- ggplot2::ggplot(obs, ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
      ggplot2::labs(x = "Observed value", y = "Predicted value") +
      ggplot2::theme_minimal()
  }

  finite_reliability <- is.finite(pred[["Reliability"]])
  use_prediction_stability <- !any(finite_reliability, na.rm = TRUE) &&
    "Prediction_stability" %in% names(pred) &&
    any(is.finite(pred[["Prediction_stability"]]), na.rm = TRUE)
  predictive_reliability <- "Reliability_basis" %in% names(pred) && any(
    gp_contract_predictive_stability_basis(pred[["Reliability_basis"]]),
    na.rm = TRUE
  )
  reliability_axis <- if (isTRUE(use_prediction_stability)) {
    "Prediction stability (bootstrap descriptive index; non-genetic)"
  } else if (isTRUE(predictive_reliability)) {
    "Predictive precision index (non-genetic reliability)"
  } else {
    "Reliability"
  }
  metric_values <- if (isTRUE(use_prediction_stability)) {
    pred[["Prediction_stability"]]
  } else {
    pred[["Reliability"]]
  }
  rel <- pred[is.finite(metric_values) & is.finite(pred[["Predicted_value"]]), , drop = FALSE]
  rel[[".diagnostic_metric"]] <- metric_values[is.finite(metric_values) & is.finite(pred[["Predicted_value"]])]
  p_rel <- NULL
  if (nrow(rel)) {
    .diagnostic_metric <- NULL
    p_rel <- ggplot2::ggplot(rel, ggplot2::aes(x = .diagnostic_metric, y = Predicted_value, color = Train_Test_Label)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(x = reliability_axis, y = "Predicted value") +
      ggplot2::theme_minimal()
  }

  int <- pred[
    is.finite(pred[["Predicted_value"]]) &
      is.finite(pred[["lower_bound"]]) &
      is.finite(pred[["upper_bound"]]),
    ,
    drop = FALSE
  ]
  p_int <- NULL
  if (nrow(int)) {
    int[["row_id"]] <- seq_len(nrow(int))
    p_int <- ggplot2::ggplot(int, ggplot2::aes(x = row_id, y = Predicted_value, ymin = lower_bound, ymax = upper_bound, color = Train_Test_Label)) +
      ggplot2::geom_errorbar(alpha = 0.55, width = 0) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::labs(x = "Prediction row", y = "Predicted value") +
      ggplot2::theme_minimal()
  }

  p_summary <- NULL
  summary_remarks_col <- if (isTRUE(use_prediction_stability)) {
    "Prediction_stability_remarks"
  } else {
    "Reliability_remarks"
  }
  if (summary_remarks_col %in% names(pred)) {
    summ <- pred[
      !is.na(pred[[summary_remarks_col]]) & nzchar(as.character(pred[[summary_remarks_col]])),
      ,
      drop = FALSE
    ]
    if (nrow(summ)) {
      tab <- as.data.frame(table(Reliability_remarks = summ[[summary_remarks_col]]), stringsAsFactors = FALSE)
      p_summary <- ggplot2::ggplot(tab, ggplot2::aes(x = Reliability_remarks, y = Freq, fill = Reliability_remarks)) +
        ggplot2::geom_col(alpha = 0.85) +
        ggplot2::labs(x = NULL, y = "Prediction count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    }
  }

  list(
    predicted_vs_observed = p_obs,
    predicted_vs_reliability = p_rel,
    prediction_interval = p_int,
    reliability_summary = p_summary
  )
}
