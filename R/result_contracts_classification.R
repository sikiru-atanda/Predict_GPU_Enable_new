gp_classification_probability_columns <- function(df) {
  if (is.null(df) || !length(names(df))) {
    return(character())
  }
  grep("^(Probability|Prob)_", names(df), value = TRUE)
}

gp_classification_probability_name <- function(class_label) {
  label <- as.character(class_label %||% "")
  label <- if (is.na(label) || !nzchar(label)) "Unknown" else label
  paste0("Probability_", make.names(label))
}

gp_classification_prediction_columns <- function(probability_columns = character(),
                                                 include_env = FALSE,
                                                 include_trait = FALSE,
                                                 id_column = "GID",
                                                 env_column = "Env") {
  leading <- as.character(id_column)[1L]
  if (isTRUE(include_env)) {
    leading <- c(leading, as.character(env_column)[1L])
  }
  if (isTRUE(include_trait)) {
    leading <- c(leading, "Trait")
  }
  c(
    leading,
    "Predicted_class", "Train_Test_Label", "Observed_class",
    "Prediction_confidence", "Classification_uncertainty",
    "Reliability", "Reliability_remarks",
    unique(as.character(probability_columns))
  )
}

gp_fill_classification_probabilities_from_confidence <- function(out) {
  if (!is.data.frame(out) || !"Predicted_class" %in% names(out)) {
    return(out)
  }
  existing <- gp_classification_probability_columns(out)
  if (length(existing)) {
    return(out)
  }
  classes <- unique(c(as.character(out$Observed_class), as.character(out$Predicted_class)))
  classes <- classes[!is.na(classes) & nzchar(classes)]
  # A single confidence value identifies the complementary probability only
  # for a binary outcome.  For three or more classes, inventing the remaining
  # probabilities would not be statistically identified.
  if (length(classes) != 2L) {
    return(out)
  }
  classes <- sort(classes)
  prob_names <- vapply(classes, gp_classification_probability_name, character(1))
  prob_names <- make.unique(prob_names)
  for (nm in prob_names) {
    out[[nm]] <- NA_real_
  }
  confidence <- if ("Prediction_confidence" %in% names(out)) {
    suppressWarnings(as.numeric(out$Prediction_confidence))
  } else {
    rep(NA_real_, nrow(out))
  }
  confidence <- pmax(pmin(confidence, 1), 0)
  pred <- as.character(out$Predicted_class)
  for (i in seq_len(nrow(out))) {
    if (is.na(pred[[i]]) || !nzchar(pred[[i]]) || is.na(confidence[[i]])) {
      next
    }
    hit <- match(pred[[i]], classes)
    if (is.na(hit)) {
      next
    }
    out[[prob_names[[hit]]]][[i]] <- confidence[[i]]
    other <- setdiff(seq_along(classes), hit)
    out[[prob_names[[other]]]][[i]] <- 1 - confidence[[i]]
  }
  out
}

gp_classification_confidence_from_probabilities <- function(out,
                                                            tolerance = 1e-6) {
  prob_cols <- gp_classification_probability_columns(out)
  if (!length(prob_cols) || !nrow(out)) {
    return(rep(NA_real_, nrow(out)))
  }
  prob <- suppressWarnings(as.matrix(data.frame(
    lapply(out[, prob_cols, drop = FALSE], as.numeric),
    check.names = FALSE
  )))
  storage.mode(prob) <- "double"
  row_sums <- rowSums(prob)
  valid <- rowSums(is.finite(prob)) == ncol(prob) &
    rowSums(prob < -tolerance | prob > 1 + tolerance) == 0L &
    is.finite(row_sums) &
    abs(row_sums - 1) <= tolerance
  confidence <- rep(NA_real_, nrow(prob))
  if (any(valid)) {
    confidence[valid] <- apply(prob[valid, , drop = FALSE], 1L, max)
    confidence[valid] <- pmax(0, pmin(1, confidence[valid]))
  }
  confidence
}

gp_complete_classification_confidence_fields <- function(out,
                                                         high_confidence = 0.9,
                                                         low_confidence = 0.5) {
  if (!is.data.frame(out) || !nrow(out)) {
    return(out)
  }
  probability_confidence <- gp_classification_confidence_from_probabilities(out)
  supplied_confidence <- suppressWarnings(as.numeric(out[["Prediction_confidence"]]))
  supplied_valid <- is.finite(supplied_confidence) &
    supplied_confidence >= 0 & supplied_confidence <= 1
  confidence <- ifelse(is.finite(probability_confidence), probability_confidence,
                       ifelse(supplied_valid, supplied_confidence, NA_real_))

  # Classification reliability is predictive class confidence, not the
  # Gaussian PEV-based reliability and not a genetic heritability quantity.
  out[["Prediction_confidence"]] <- confidence
  out[["Classification_uncertainty"]] <- ifelse(
    is.finite(confidence), 1 - confidence, NA_real_
  )
  out[["Reliability"]] <- confidence

  remarks <- as.character(out[["Reliability_remarks"]])
  missing_remarks <- is.na(remarks) | !nzchar(remarks)
  derived_remarks <- ifelse(
    confidence >= high_confidence, "High Confidence",
    ifelse(confidence <= low_confidence, "Low Confidence", "Moderate Confidence")
  )
  remarks[missing_remarks & is.finite(confidence)] <-
    derived_remarks[missing_remarks & is.finite(confidence)]
  out[["Reliability_remarks"]] <- remarks
  out
}

gp_is_classification_prediction_table <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    return(FALSE)
  }
  if (length(gp_classification_probability_columns(df)) > 0L) {
    return(TRUE)
  }
  class_cols <- intersect(c("Predicted_class", "Observed_class"), names(df))
  if (length(class_cols) &&
      any(vapply(class_cols, function(nm) gp_contract_nonempty_column(df[[nm]]), logical(1)))) {
    return(TRUE)
  }
  value_cols <- intersect(c("Predicted_value", "Observed_value"), names(df))
  if (!length(value_cols)) {
    return(FALSE)
  }
  any(vapply(value_cols, function(nm) {
    vals <- as.character(df[[nm]])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) {
      return(FALSE)
    }
    numeric_vals <- suppressWarnings(as.numeric(vals))
    any(is.na(numeric_vals))
  }, logical(1)))
}

gp_format_classification_prediction_table <- function(x,
                                                      gen_name = NULL,
                                                      heter_groups = NULL,
                                                      include_env = NULL,
                                                      include_trait = NULL) {
  df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    return(x)
  }
  n <- nrow(df)
  names_x <- names(df)
  gid_col <- gp_contract_first_name(
    names_x,
    unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
  )
  pred_col <- gp_contract_first_name(names_x, c("Predicted_class", "Predicted_value", "class", "prediction"))
  obs_col <- gp_contract_first_name(names_x, c("Observed_class", "Observed_value", "Observed", "observed", "y"))
  train_col <- gp_contract_first_name(names_x, c("Train_Test_Label", "Train_Test", "set", "Set"))
  conf_col <- gp_contract_first_name(names_x, c("Prediction_confidence", "Predicted_probability", "confidence"))
  unc_col <- gp_contract_first_name(names_x, c("Classification_uncertainty", "Uncertainty"))
  rel_col <- gp_contract_first_name(names_x, c("Reliability", "Prediction_confidence"))
  rel_rem_col <- gp_contract_first_name(names_x, c("Reliability_remarks", "Prediction_confidence_remarks"))
  # Honor user-supplied heter_groups as the authoritative env-column name.
  env_col <- if (!is.null(heter_groups) &&
                 length(heter_groups) == 1L &&
                 nzchar(as.character(heter_groups)) &&
                 as.character(heter_groups) %in% names_x) {
    as.character(heter_groups)
  } else {
    gp_contract_first_name(names_x, c("Env", "env", "Environment", "environment", "Location", "location"))
  }
  trait_col <- gp_contract_first_name(names_x, c("Trait", "trait", "Response", "response_name"))

  gid_out_name <- if (!is.null(gen_name) && length(gen_name) == 1L && nzchar(as.character(gen_name))) {
    as.character(gen_name)
  } else {
    "GID"
  }
  env_out_name <- if (!is.null(heter_groups) && length(heter_groups) == 1L && nzchar(as.character(heter_groups))) {
    as.character(heter_groups)
  } else {
    "Env"
  }

  out <- data.frame(
    .gid = if (!is.na(gid_col)) as.character(df[[gid_col]]) else rep(NA_character_, n),
    Predicted_class = if (!is.na(pred_col)) as.character(df[[pred_col]]) else rep(NA_character_, n),
    Train_Test_Label = if (!is.na(train_col)) as.character(df[[train_col]]) else rep(NA_character_, n),
    Observed_class = if (!is.na(obs_col)) as.character(df[[obs_col]]) else rep(NA_character_, n),
    Prediction_confidence = if (!is.na(conf_col)) gp_contract_as_numeric(df[[conf_col]]) else gp_contract_fill_numeric(n),
    Classification_uncertainty = if (!is.na(unc_col)) gp_contract_as_numeric(df[[unc_col]]) else gp_contract_fill_numeric(n),
    Reliability = if (!is.na(rel_col)) gp_contract_as_numeric(df[[rel_col]]) else gp_contract_fill_numeric(n),
    Reliability_remarks = if (!is.na(rel_rem_col)) as.character(df[[rel_rem_col]]) else rep(NA_character_, n),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out)[names(out) == ".gid"] <- gid_out_name
  prob_cols <- gp_classification_probability_columns(df)
  for (nm in prob_cols) {
    public_nm <- sub("^Prob_", "Probability_", nm)
    out[[public_nm]] <- gp_contract_as_numeric(df[[nm]])
  }
  out <- gp_fill_classification_probabilities_from_confidence(out)
  out <- gp_complete_classification_confidence_fields(out)
  prob_cols <- sub("^Prob_", "Probability_", prob_cols)
  prob_cols <- unique(c(prob_cols, gp_classification_probability_columns(out)))

  if (is.null(include_env)) {
    include_env <- !is.na(env_col) && gp_contract_multiple_nonempty_values(df[[env_col]])
  }
  if (is.null(include_trait)) {
    include_trait <- !is.na(trait_col) && gp_contract_multiple_nonempty_values(df[[trait_col]])
  }
  leading_cols <- gid_out_name
  if (isTRUE(include_env)) {
    out[[env_out_name]] <- if (!is.na(env_col)) as.character(df[[env_col]]) else rep(NA_character_, n)
    leading_cols <- c(leading_cols, env_out_name)
  }
  if (isTRUE(include_trait)) {
    out[["Trait"]] <- if (!is.na(trait_col)) as.character(df[[trait_col]]) else rep(NA_character_, n)
    leading_cols <- c(leading_cols, "Trait")
  }

  out[, gp_classification_prediction_columns(
    probability_columns = prob_cols,
    include_env = isTRUE(include_env),
    include_trait = isTRUE(include_trait),
    id_column = gid_out_name,
    env_column = env_out_name
  ), drop = FALSE]
}

gp_merge_classification_observed_values <- function(pred,
                                                     residual = NULL,
                                                     pheno_data = NULL,
                                                     response = NULL,
                                                     gen_name = NULL,
                                                     heter_groups = NULL) {
  pred_df <- tryCatch(as.data.frame(pred, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  resid_df <- tryCatch(as.data.frame(residual, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  pheno_df <- tryCatch(as.data.frame(pheno_data, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(pred_df) || !nrow(pred_df)) {
    return(pred)
  }

  pred_obs_col <- gp_contract_first_name(names(pred_df), c("Observed_class", "Observed_value", "Observed", "observed", "y"))
  if (is.na(pred_obs_col)) {
    pred_obs_col <- "Observed_class"
    pred_df[[pred_obs_col]] <- rep(NA_character_, nrow(pred_df))
  }
  pred_id_col <- gp_contract_first_name(
    names(pred_df),
    unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
  )
  if (is.na(pred_id_col)) {
    return(pred_df)
  }

  make_key <- function(df, cols) {
    do.call(paste, c(lapply(cols, function(nm) as.character(df[[nm]])), sep = "\r"))
  }

  fill_from_source <- function(out, source, source_obs_candidates) {
    if (is.null(source) || !nrow(source)) {
      return(out)
    }
    source_obs_col <- gp_contract_first_name(names(source), source_obs_candidates)
    source_id_col <- gp_contract_first_name(
      names(source),
      unique(stats::na.omit(c("GID", "gid", gen_name, "Name", "name", "ID", "id")))
    )
    if (is.na(source_obs_col) || is.na(source_id_col)) {
      return(out)
    }
    pred_key_cols <- pred_id_col
    source_key_cols <- source_id_col
    pred_env_col <- gp_contract_first_name(
      names(out),
      unique(stats::na.omit(c(heter_groups, "Env", "env", "Environment", "environment")))
    )
    source_env_col <- gp_contract_first_name(
      names(source),
      unique(stats::na.omit(c(heter_groups, "Env", "env", "Environment", "environment")))
    )
    if (!is.na(pred_env_col) && !is.na(source_env_col)) {
      pred_key_cols <- c(pred_key_cols, pred_env_col)
      source_key_cols <- c(source_key_cols, source_env_col)
    }
    idx <- match(make_key(out, pred_key_cols), make_key(source, source_key_cols))
    candidate <- rep(NA_character_, nrow(out))
    matched <- !is.na(idx)
    candidate[matched] <- as.character(source[[source_obs_col]][idx[matched]])
    current <- as.character(out[[pred_obs_col]])
    missing <- is.na(current) | !nzchar(trimws(current)) | toupper(trimws(current)) == "NA"
    usable <- !is.na(candidate) & nzchar(trimws(candidate)) &
      toupper(trimws(candidate)) != "NA"
    current[missing & usable] <- candidate[missing & usable]
    out[[pred_obs_col]] <- current
    out
  }

  pred_df <- fill_from_source(
    pred_df,
    resid_df,
    c("Observed_class", "Observed_value", "Observed", "observed", "y")
  )
  response_candidates <- unique(stats::na.omit(c(
    as.character(response %||% character()),
    "Observed_class", "Observed_value", "Observed", "observed", "y"
  )))
  pred_df <- fill_from_source(pred_df, pheno_df, response_candidates)
  pred_df
}

gp_is_public_classification_table <- function(x) {
  is.data.frame(x) &&
    all(c(
      "GID", "Predicted_class", "Train_Test_Label", "Observed_class",
      "Prediction_confidence", "Classification_uncertainty"
    ) %in% names(x))
}

gp_classification_diagnostic_plots <- function(pred) {
  if (!gp_is_public_classification_table(pred)) {
    return(NULL)
  }

  prob_cols <- gp_classification_probability_columns(pred)
  observed_rows <- pred[!is.na(pred[["Observed_class"]]) & nzchar(as.character(pred[["Observed_class"]])), , drop = FALSE]

  p_observed <- NULL
  if (nrow(observed_rows)) {
    cm <- as.data.frame(
      as.table(table(
        Observed_class = observed_rows[["Observed_class"]],
        Predicted_class = observed_rows[["Predicted_class"]]
      )),
      stringsAsFactors = FALSE
    )
    p_observed <- ggplot2::ggplot(cm, ggplot2::aes(x = Observed_class, y = Predicted_class, fill = Freq)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = Freq), size = 3) +
      ggplot2::scale_fill_gradient(low = "#eff6ff", high = "#2563eb") +
      ggplot2::labs(x = "Observed class", y = "Predicted class", fill = "Count") +
      ggplot2::theme_minimal()
  }

  conf_df <- pred[is.finite(pred[["Prediction_confidence"]]), , drop = FALSE]
  p_conf <- NULL
  if (nrow(conf_df)) {
    p_conf <- ggplot2::ggplot(
      conf_df,
      ggplot2::aes(x = Train_Test_Label, y = Prediction_confidence, fill = Train_Test_Label)
    ) +
      ggplot2::geom_boxplot(alpha = 0.75, outlier.alpha = 0.5) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(x = NULL, y = "Prediction confidence") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }

  p_prob <- NULL
  if (length(prob_cols)) {
    prob_long <- stats::reshape(
      pred[, c("GID", "Train_Test_Label", prob_cols), drop = FALSE],
      varying = prob_cols,
      v.names = "Probability",
      timevar = "Class",
      times = sub("^Probability_", "", prob_cols),
      direction = "long"
    )
    rownames(prob_long) <- NULL
    prob_long[["Probability"]] <- suppressWarnings(as.numeric(prob_long[["Probability"]]))
    prob_long <- prob_long[is.finite(prob_long[["Probability"]]), , drop = FALSE]
    if (nrow(prob_long)) {
      p_prob <- ggplot2::ggplot(prob_long, ggplot2::aes(x = Class, y = Probability, fill = Class)) +
        ggplot2::geom_boxplot(alpha = 0.75, outlier.alpha = 0.4) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = "Class", y = "Posterior/probability") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    }
  }

  p_rel <- NULL
  if ("Reliability_remarks" %in% names(pred)) {
    rel <- pred[!is.na(pred[["Reliability_remarks"]]) & nzchar(as.character(pred[["Reliability_remarks"]])), , drop = FALSE]
    if (nrow(rel)) {
      rel_tab <- as.data.frame(table(Reliability_remarks = rel[["Reliability_remarks"]]), stringsAsFactors = FALSE)
      p_rel <- ggplot2::ggplot(rel_tab, ggplot2::aes(x = Reliability_remarks, y = Freq, fill = Reliability_remarks)) +
        ggplot2::geom_col(alpha = 0.85) +
        ggplot2::labs(x = NULL, y = "Prediction count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    }
  }

  list(
    predicted_vs_observed = p_observed,
    predicted_vs_reliability = p_conf,
    prediction_interval = p_prob,
    reliability_summary = p_rel
  )
}
