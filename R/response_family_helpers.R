# Response-family helpers

`%||%` <- function(a, b) if (!is.null(a)) a else b

gp_supported_response_families <- function() {
  c("gaussian", "binary", "ordinal", "multiclass")
}

gp_normalize_response_family <- function(response_family = NULL) {
  if (is.null(response_family) || length(response_family) == 0L) {
    return("auto")
  }

  fam <- trimws(as.character(response_family)[1L])
  if (is.na(fam) || !nzchar(fam)) {
    return("auto")
  }

  fam <- tolower(fam)
  aliases <- c(
    continuous = "gaussian",
    numeric = "gaussian",
    regression = "gaussian",
    gaussian = "gaussian",
    binomial = "binary",
    binary = "binary",
    ordered = "ordinal",
    `ordered-categorical` = "ordinal",
    ordered_categorical = "ordinal",
    ordinal = "ordinal",
    multiclass = "multiclass",
    `multi-class` = "multiclass",
    multi_class = "multiclass",
    multinomial = "multiclass",
    nominal = "multiclass",
    unordered = "multiclass",
    `unordered-categorical` = "multiclass",
    unordered_categorical = "multiclass",
    categorical = "multiclass",
    auto = "auto"
  )

  alias <- unname(aliases[fam])
  fam <- if (!is.na(alias)) alias else fam
  valid <- c("auto", gp_supported_response_families())
  if (!fam %in% valid) {
    stop(
      "Unsupported response_family. Choose from: ",
      paste(valid, collapse = ", "),
      call. = FALSE
    )
  }
  fam
}

gp_infer_response_family <- function(y) {
  if (is.null(y)) {
    return("gaussian")
  }

  y_no_na <- y[!is.na(y)]
  if (!length(y_no_na)) {
    return("gaussian")
  }

  if (is.ordered(y)) {
    return("ordinal")
  }
  if (is.factor(y) || is.character(y) || is.logical(y)) {
    nlev <- length(unique(as.character(y_no_na)))
    return(if (nlev <= 2L) "binary" else "multiclass")
  }

  if (is.numeric(y_no_na) || is.integer(y_no_na)) {
    uniq <- sort(unique(y_no_na))
    integer_like <- all(is.finite(uniq)) &&
      all(abs(uniq - round(uniq)) < sqrt(.Machine$double.eps))
    if (length(uniq) == 2L &&
        integer_like &&
        (all(uniq %in% c(0, 1)) || all(uniq %in% c(1, 2)))) {
      return("binary")
    }
    return("gaussian")
  }

  "gaussian"
}

gp_resolve_response_family <- function(response_family = NULL, y = NULL) {
  fam <- gp_normalize_response_family(response_family)
  if (identical(fam, "auto")) {
    return(gp_infer_response_family(y))
  }
  fam
}

gp_resolve_response_family_for_responses <- function(pheno_data, response, response_family = NULL) {
  fam <- gp_normalize_response_family(response_family)
  if (!identical(fam, "auto")) {
    return(fam)
  }
  if (is.null(pheno_data) || is.null(response) || !length(response)) {
    return("gaussian")
  }
  inferred <- unique(vapply(response, function(col) gp_infer_response_family(pheno_data[[col]]), character(1)))
  if (length(inferred) > 1L) {
    stop(
      "Multiple response families were inferred across the requested traits: ",
      paste(inferred, collapse = ", "),
      ". Please provide a single explicit response_family for this run.",
      call. = FALSE
    )
  }
  inferred[[1]]
}

gp_eval_metric_catalog <- function() {
  list(
    gaussian = c(
      "accuracy",
      "mean_squared_error",
      "bias",
      "root_mean_squared_error",
      "relative_squared_error",
      "mean_absolute_error",
      "mean_absolute_percent_error",
      "kendalls_tau"
    ),
    binary = c(
      "accuracy",
      "balanced_accuracy",
      "precision",
      "recall",
      "specificity",
      "f1",
      "log_loss",
      "brier_score",
      "mcc",
      "ece"
    ),
    ordinal = c(
      "accuracy",
      "balanced_accuracy",
      "mean_absolute_error_class",
      "within_one_class_accuracy",
      "quadratic_weighted_kappa",
      "log_loss",
      "ece"
    ),
    multiclass = c(
      "accuracy",
      "balanced_accuracy",
      "macro_precision",
      "macro_recall",
      "macro_f1",
      "log_loss",
      "ece"
    )
  )
}

gp_cv_selection_defaults <- function(response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  switch(
    fam,
    gaussian = list(primary_metric = "accuracy", tie_breakers = character()),
    binary = list(primary_metric = "balanced_accuracy", tie_breakers = "log_loss"),
    multiclass = list(primary_metric = "macro_f1", tie_breakers = "log_loss"),
    ordinal = list(
      primary_metric = "quadratic_weighted_kappa",
      tie_breakers = c("mean_absolute_error_class", "log_loss")
    )
  )
}

gp_resolve_metric_for_ranking <- function(metric_for_ranking = "auto",
                                          response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  metric <- as.character(metric_for_ranking %||% "auto")
  if (!length(metric) || is.na(metric[[1L]]) || !nzchar(trimws(metric[[1L]]))) {
    metric <- "auto"
  }
  if (length(metric) != 1L) {
    stop("metric_for_ranking must be one metric or 'auto'.", call. = FALSE)
  }
  metric <- tolower(trimws(metric[[1L]]))
  if (identical(metric, "auto")) {
    metric <- gp_cv_selection_defaults(fam)$primary_metric
  }
  allowed <- gp_eval_metric_catalog()[[fam]]
  if (!metric %in% allowed) {
    stop(
      "metric_for_ranking='", metric, "' is not valid for response_family='",
      fam, "'. Allowed metrics: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  metric
}

gp_resolve_ranking_tie_breakers <- function(ranking_tie_breakers = NULL,
                                            response_family = "gaussian",
                                            metric_for_ranking = "auto") {
  fam <- gp_resolve_response_family(response_family)
  primary <- gp_resolve_metric_for_ranking(metric_for_ranking, fam)
  tie_breakers <- if (is.null(ranking_tie_breakers)) {
    gp_cv_selection_defaults(fam)$tie_breakers
  } else {
    tolower(trimws(as.character(ranking_tie_breakers)))
  }
  if (length(tie_breakers) == 1L && identical(tie_breakers, "none")) {
    return(character())
  }
  tie_breakers <- unique(tie_breakers[!is.na(tie_breakers) & nzchar(tie_breakers)])
  tie_breakers <- setdiff(tie_breakers, primary)
  allowed <- gp_eval_metric_catalog()[[fam]]
  bad <- setdiff(tie_breakers, allowed)
  if (length(bad)) {
    stop(
      "Unsupported ranking_tie_breakers for response_family='", fam, "': ",
      paste(bad, collapse = ", "), ". Allowed metrics: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  tie_breakers
}

gp_prepare_cv_metric_selection <- function(eval_metrics = NULL,
                                           metric_for_ranking = "auto",
                                           ranking_tie_breakers = NULL,
                                           response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  primary <- gp_resolve_metric_for_ranking(metric_for_ranking, fam)
  tie_breakers <- gp_resolve_ranking_tie_breakers(
    ranking_tie_breakers = ranking_tie_breakers,
    response_family = fam,
    metric_for_ranking = primary
  )
  metrics <- if (is.null(eval_metrics)) {
    if (identical(fam, "gaussian")) primary else gp_eval_metric_catalog()[[fam]]
  } else {
    tolower(as.character(eval_metrics))
  }
  metrics <- unique(c(metrics, primary, tie_breakers))
  gp_validate_eval_metrics(metrics, fam)
  list(
    response_family = fam,
    eval_metrics = metrics,
    metric_for_ranking = primary,
    ranking_tie_breakers = tie_breakers
  )
}

gp_eval_metric_direction <- function(metric) {
  metric <- tolower(as.character(metric)[1L])
  minimize <- c(
    "mean_squared_error", "bias", "root_mean_squared_error",
    "relative_squared_error", "mean_absolute_error",
    "mean_absolute_percent_error", "log_loss", "brier_score",
    "mean_absolute_error_class", "ece"
  )
  maximize <- c(
    "accuracy", "kendalls_tau", "spearman", "pearson",
    "correlation", "rsq", "r2", "balanced_accuracy", "precision",
    "recall", "specificity", "f1", "mcc", "macro_precision",
    "macro_recall", "macro_f1", "within_one_class_accuracy",
    "quadratic_weighted_kappa"
  )
  if (metric %in% minimize) return("minimize")
  if (metric %in% maximize) return("maximize")
  NA_character_
}

gp_validate_eval_metrics <- function(eval_metrics = NULL, response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  if (is.null(eval_metrics)) {
    return(invisible(fam))
  }
  metrics <- tolower(as.character(eval_metrics))
  allowed <- gp_eval_metric_catalog()[[fam]]
  bad <- setdiff(metrics, allowed)
  if (length(bad) > 0L) {
    stop(
      "Unsupported eval_metrics for response_family='", fam, "': ",
      paste(bad, collapse = ", "),
      ". Allowed metrics: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(fam)
}

gp_validate_model_response_family <- function(model = NULL, response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  if (is.null(model) || length(model) == 0L || is.na(model)) {
    return(invisible(fam))
  }

  model_label <- as.character(model)[1]
  # display names (e.g. "Kernel-GBLUP") map to the internal names used below
  model <- gp_canonicalize_supported_model_names(model_label)[1]
  if (model %in% gp_deep_learning_supported_models()) {
    allowed <- c("gaussian", "binary", "multiclass")
  } else {
    allowed <- switch(
    model,
    BayesA = c("gaussian", "binary", "ordinal"),
    BayesB = c("gaussian", "binary", "ordinal"),
    BayesC = c("gaussian", "binary", "ordinal"),
    BL = c("gaussian", "binary", "ordinal"),
    BRR = c("gaussian", "binary", "ordinal"),
    GBLUP_BRR = c("gaussian", "binary", "ordinal"),
    RKHS = c("gaussian", "binary", "ordinal"),
    KRR = "gaussian",
    GP = "gaussian",
    GP_FA = "gaussian",
    LowRankGP = "gaussian",
    GBLUP = "gaussian",
    PartialLeastSquare = "gaussian",
    Lasso = "gaussian",
    Ridge_Regression = "gaussian",
    Xgboost = c("gaussian", "binary", "multiclass"),
    RandomForest = c("gaussian", "binary", "multiclass"),
    CatBoost = c("gaussian", "binary", "multiclass"),
    LightGBM = c("gaussian", "binary", "multiclass"),
    SupportVectorMachine = c("gaussian", "binary", "multiclass"),
    `K-NearestNeighbors` = c("gaussian", "binary", "multiclass"),
    c("gaussian", "binary", "multiclass")
    )
  }

  if (!fam %in% allowed) {
    stop(
      "Model '", model_label, "' does not support response_family='", fam,
      "'. Allowed response families: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(fam)
}

gp_to_class_labels <- function(y) {
  if (is.factor(y)) return(as.character(y))
  if (is.logical(y)) return(as.character(y))
  if (is.character(y)) return(y)
  as.character(y)
}

gp_response_class_levels <- function(y, response_family = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  if (identical(fam, "gaussian")) {
    return(NULL)
  }
  if (is.factor(y) || is.ordered(y)) {
    return(as.character(levels(y)))
  }
  vals <- unique(as.character(stats::na.omit(y)))
  if (identical(fam, "ordinal")) {
    return(vals)
  }
  sort(vals)
}

gp_resolve_positive_class <- function(y,
                                      response_family = "binary",
                                      positive_class = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y)
  if (!identical(fam, "binary")) {
    if (!is.null(positive_class)) {
      stop("positive_class is only valid for response_family='binary'.", call. = FALSE)
    }
    return(NULL)
  }
  levs <- gp_response_class_levels(y, fam)
  if (length(levs) != 2L) {
    stop("Binary responses must contain exactly two class levels.", call. = FALSE)
  }
  if (is.null(positive_class)) {
    return(levs[[2L]])
  }
  positive_class <- as.character(positive_class)
  if (length(positive_class) != 1L || is.na(positive_class) || !nzchar(positive_class)) {
    stop("positive_class must name exactly one non-missing binary class.", call. = FALSE)
  }
  if (!positive_class %in% levs) {
    stop(
      "positive_class='", positive_class, "' is not present in the binary response. Available classes: ",
      paste(levs, collapse = ", "),
      call. = FALSE
    )
  }
  positive_class
}

gp_apply_binary_positive_class <- function(y, positive_class = NULL) {
  positive_class <- gp_resolve_positive_class(y, "binary", positive_class)
  levs <- gp_response_class_levels(y, "binary")
  ordered_levels <- c(setdiff(levs, positive_class), positive_class)
  factor(as.character(y), levels = ordered_levels, ordered = is.ordered(y))
}

gp_configure_positive_class_responses <- function(pheno_data,
                                                  response,
                                                  response_family,
                                                  positive_class = NULL) {
  fam <- gp_resolve_response_family(response_family)
  if (!identical(fam, "binary")) {
    if (!is.null(positive_class)) {
      stop("positive_class is only valid for response_family='binary'.", call. = FALSE)
    }
    return(list(pheno_data = pheno_data, positive_class = NULL))
  }
  response <- as.character(response)
  resolved <- vapply(response, function(trait) {
    gp_resolve_positive_class(pheno_data[[trait]], fam, positive_class)
  }, character(1L))
  if (length(unique(resolved)) != 1L) {
    stop(
      "All binary response traits in one run must use the same positive_class.",
      call. = FALSE
    )
  }
  for (trait in response) {
    pheno_data[[trait]] <- gp_apply_binary_positive_class(
      pheno_data[[trait]], resolved[[1L]]
    )
  }
  list(pheno_data = pheno_data, positive_class = resolved[[1L]])
}

gp_binary_positive_probability <- function(y_predicted,
                                           class_levels,
                                           positive_class) {
  if (!is.matrix(y_predicted) && !is.data.frame(y_predicted)) {
    return(as.numeric(y_predicted))
  }
  prob <- as.matrix(y_predicted)
  if (ncol(prob) == 1L) {
    return(as.numeric(prob[, 1L]))
  }
  prob_names <- colnames(prob)
  if (!is.null(prob_names)) {
    cleaned <- sub("^(Probability|Prob)_", "", prob_names)
    idx <- match(make.names(positive_class), make.names(cleaned))
    if (!is.na(idx)) {
      return(as.numeric(prob[, idx]))
    }
  }
  idx <- match(positive_class, class_levels)
  if (is.na(idx) || idx > ncol(prob)) {
    stop("Unable to identify the positive-class probability column.", call. = FALSE)
  }
  as.numeric(prob[, idx])
}

gp_pred_to_class_labels <- function(y_predicted, positive_class = NULL, class_levels = NULL) {
  if (is.factor(y_predicted) || is.character(y_predicted) || is.logical(y_predicted)) {
    return(as.character(y_predicted))
  }
  class_levels <- if (!is.null(class_levels) && length(class_levels)) as.character(class_levels) else NULL
  if (is.data.frame(y_predicted)) {
    y_predicted <- as.matrix(y_predicted)
  }
  if (is.matrix(y_predicted)) {
    if (ncol(y_predicted) == 1L) {
      levs <- class_levels %||% c("0", positive_class %||% "1")
      pos <- positive_class %||% levs[min(2L, length(levs))]
      neg <- setdiff(levs, pos)[[1L]]
      return(as.character(ifelse(y_predicted[, 1] >= 0.5, pos, neg)))
    }
    cls <- if (!is.null(class_levels) && length(class_levels) == ncol(y_predicted)) {
      class_levels
    } else {
      colnames(y_predicted) %||% as.character(seq_len(ncol(y_predicted)))
    }
    return(as.character(cls[max.col(y_predicted, ties.method = "first")]))
  }
  levs <- class_levels %||% c("0", positive_class %||% "1")
  pos <- positive_class %||% levs[min(2L, length(levs))]
  neg <- setdiff(levs, pos)[[1L]]
  ifelse(as.numeric(y_predicted) >= 0.5, pos, neg)
}

gp_safe_divide <- function(num, den) {
  ifelse(is.finite(den) & den != 0, num / den, NA_real_)
}

gp_metric_valid_rows <- function(y_true, y_pred, response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family, y = y_true)
  if (fam == "gaussian") {
    return(is.finite(y_true) & is.finite(y_pred))
  }
  if (is.matrix(y_pred) || is.data.frame(y_pred)) {
    return(!is.na(y_true))
  }
  !is.na(y_true) & !is.na(y_pred)
}

gp_confusion_from_labels <- function(y_true, y_pred, levels = NULL) {
  if (is.null(levels)) {
    levels <- sort(unique(c(y_true, y_pred)))
  }
  table(
    factor(y_true, levels = levels),
    factor(y_pred, levels = levels)
  )
}

gp_multiclass_precision_recall_f1 <- function(cm) {
  tp <- diag(cm)
  fp <- colSums(cm) - tp
  fn <- rowSums(cm) - tp
  precision <- gp_safe_divide(tp, tp + fp)
  recall <- gp_safe_divide(tp, tp + fn)
  f1 <- gp_safe_divide(2 * precision * recall, precision + recall)
  list(
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1 = mean(f1, na.rm = TRUE)
  )
}

gp_balanced_accuracy <- function(cm) {
  sens <- gp_safe_divide(diag(cm), rowSums(cm))
  mean(sens, na.rm = TRUE)
}

gp_log_loss <- function(y_true, y_predicted, response_family, positive_class = NULL) {
  eps <- 1e-15
  fam <- gp_resolve_response_family(response_family, y = y_true)
  if (fam == "binary") {
    levs <- gp_response_class_levels(y_true, fam)
    positive_class <- gp_resolve_positive_class(y_true, fam, positive_class)
    p <- gp_binary_positive_probability(y_predicted, levs, positive_class)
    p <- pmin(pmax(p, eps), 1 - eps)
    labs <- gp_to_class_labels(y_true)
    y01 <- as.numeric(labs == positive_class)
    return(-mean(y01 * log(p) + (1 - y01) * log(1 - p)))
  }

  prob <- as.matrix(y_predicted)
  prob <- pmin(pmax(prob, eps), 1 - eps)
  prob <- prob / rowSums(prob)
  labs <- gp_to_class_labels(y_true)
  cls <- colnames(prob) %||% gp_response_class_levels(y_true, fam) %||% sort(unique(labs))
  if (!is.null(cls)) {
    cls <- sub("^Prob_", "", cls)
    cls <- make.names(cls, unique = FALSE)
    labs <- make.names(labs, unique = FALSE)
  }
  idx <- match(labs, cls)
  if (anyNA(idx)) {
    return(NA_real_)
  }
  return(-mean(log(prob[cbind(seq_along(idx), idx)])))
}

gp_expected_calibration_error <- function(y_true,
                                          y_predicted,
                                          response_family = "binary",
                                          bins = 10L,
                                          positive_class = NULL) {
  fam <- gp_resolve_response_family(response_family, y = y_true)
  labs <- gp_to_class_labels(y_true)

  if (fam == "binary") {
    levs <- gp_response_class_levels(y_true, fam)
    positive_class <- gp_resolve_positive_class(y_true, fam, positive_class)
    p <- gp_binary_positive_probability(y_predicted, levs, positive_class)
    ok <- is.finite(p) & !is.na(labs)
    p <- p[ok]
    labs <- labs[ok]
    if (!length(p)) {
      return(NA_real_)
    }
    p <- pmin(pmax(p, 0), 1)
    negative_class <- setdiff(levs, positive_class)[[1L]]
    predicted_class <- ifelse(p >= 0.5, positive_class, negative_class)
    confidence <- ifelse(p >= 0.5, p, 1 - p)
    correct <- as.numeric(labs == predicted_class)
  } else {
    prob <- as.matrix(y_predicted)
    if (!nrow(prob)) {
      return(NA_real_)
    }
    row_ok <- rowSums(is.finite(prob)) == ncol(prob) & !is.na(labs)
    prob <- prob[row_ok, , drop = FALSE]
    labs <- labs[row_ok]
    if (!nrow(prob)) {
      return(NA_real_)
    }
    prob <- prob / rowSums(prob)
    cls <- colnames(prob) %||% gp_response_class_levels(y_true, fam) %||% sort(unique(labs))
    pred_idx <- max.col(prob, ties.method = "first")
    pred_class <- cls[pred_idx]
    confidence <- prob[cbind(seq_len(nrow(prob)), pred_idx)]
    correct <- as.numeric(pred_class == labs)
  }

  bin_ids <- cut(
    confidence,
    breaks = seq(0, 1, length.out = bins + 1L),
    include.lowest = TRUE,
    labels = FALSE
  )
  calib <- data.frame(confidence = confidence, correct = correct, bin = bin_ids)
  calib <- calib[!is.na(calib$bin), , drop = FALSE]
  if (!nrow(calib)) {
    return(NA_real_)
  }

  by_bin <- calib |>
    dplyr::group_by(.data$bin) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_confidence = mean(.data$confidence, na.rm = TRUE),
      observed_accuracy = mean(.data$correct, na.rm = TRUE),
      .groups = "drop"
    )

  sum((by_bin$n / sum(by_bin$n)) * abs(by_bin$observed_accuracy - by_bin$mean_confidence))
}

gp_quadratic_weighted_kappa <- function(y_true, y_pred, levels = NULL) {
  # Distance weights need the ordinal class order; sorting labels would order
  # them alphabetically (low < medium < high became high, low, medium).
  labs <- if (length(levels)) {
    as.character(levels)
  } else {
    sort(unique(c(y_true, y_pred)))
  }
  cm <- gp_confusion_from_labels(y_true, y_pred, levels = labs)
  n <- sum(cm)
  if (n == 0) return(NA_real_)
  k <- length(labs)
  w <- outer(seq_len(k), seq_len(k), function(i, j) ((i - j)^2) / ((k - 1)^2))
  obs <- cm / n
  exp <- outer(rowSums(cm), colSums(cm)) / (n^2)
  1 - sum(w * obs) / sum(w * exp)
}

gp_classification_prediction_summary <- function(prob, class_levels = NULL,
                                                 high_confidence = 0.9,
                                                 low_confidence = 0.5) {
  if (is.data.frame(prob)) {
    prob <- as.matrix(prob)
  }

  if (!is.matrix(prob) && !is.null(class_levels) && length(class_levels) > 2L) {
    p <- as.numeric(prob)
    n_classes <- length(class_levels)
    if (length(p) %% n_classes == 0L) {
      prob <- matrix(p, ncol = n_classes, byrow = FALSE)
    }
  }

  if (is.matrix(prob)) {
    levs <- class_levels %||% colnames(prob) %||% as.character(seq_len(ncol(prob)))
    idx <- max.col(prob, ties.method = "first")
    pred_class <- levs[idx]
    confidence <- prob[cbind(seq_len(nrow(prob)), idx)]
  } else {
    p <- as.numeric(prob)
    levs <- class_levels %||% c("0", "1")
    pred_class <- ifelse(p >= 0.5, levs[min(2, length(levs))], levs[1])
    confidence <- pmax(p, 1 - p)
  }

  remarks <- ifelse(
    confidence >= high_confidence, "High Confidence",
    ifelse(confidence <= low_confidence, "Low Confidence", "Moderate Confidence")
  )

  data.frame(
    Predicted_class = as.character(pred_class),
    Prediction_confidence = as.numeric(confidence),
    Classification_uncertainty = as.numeric(1 - confidence),
    Prediction_confidence_remarks = remarks,
    Prediction_confidence_percentage = mean(confidence, na.rm = TRUE) * 100,
    stringsAsFactors = FALSE
  )
}
