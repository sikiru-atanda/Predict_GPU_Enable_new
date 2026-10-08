#' Evaluate Prediction Performance Metrics
#'
#' This function computes various evaluation metrics to assess the performance of prediction models.
#'
#' @param y_observed A numeric vector of observed values.
#' @param y_predicted Predicted values. For classification this can be class labels,
#' class probabilities for the positive class, or a probability matrix with one
#' column per class.
#' @param eval_metrics A character vector specifying which evaluation metrics to compute.
#' @param response_family Response family. One of `gaussian`, `binary`,
#' `ordinal`, or `multiclass`.
#' @param positive_class For binary responses, the class represented by a
#'   one-column probability vector and used for precision, recall, specificity,
#'   F1, MCC, Brier score, and binary log loss. When omitted, the second factor
#'   level is used.
#' @return The function returns a numeric value representing the computed evaluation metric specified by `eval_metrics`.
#' @examples
#' \dontrun{
#' observed <- c(1, 2, 3, 4, 5)
#' predicted <- c(1.1, 1.9, 3.1, 3.9, 4.8)
#' mse <- evaluation_metrics(y_observed = observed,
#'                           y_predicted = predicted,
#'                           eval_metrics = "Mean_Squared_Error")
#' print(mse)
#'
#' mae <- evaluation_metrics(y_observed = observed,
#'                           y_predicted = predicted,
#'                           eval_metrics = "Mean_Absolute_Error")
#' print(mae)
#'}
#' @export
#'
evaluation_metrics <- function(y_observed = NULL,
                               y_predicted = NULL,
                               eval_metrics = NULL,
                               response_family = "gaussian",
                               positive_class = NULL) {

  msg <- ""

  if (is.null(y_observed) || is.null(y_predicted)) {
    stop("Both observed and predicted values must be provided.")
  }
  fam <- gp_resolve_response_family(response_family, y = y_observed)
  if (!identical(fam, "binary") && !is.null(positive_class)) {
    stop("positive_class is only valid for response_family='binary'.", call. = FALSE)
  }
  gp_validate_eval_metrics(eval_metrics = eval_metrics, response_family = fam)

  if (!is.matrix(y_predicted) && !is.data.frame(y_predicted) && fam %in% c("multiclass", "ordinal")) {
    n_obs <- length(y_observed)
    pred_num <- suppressWarnings(as.numeric(y_predicted))
    cls_levels <- sort(unique(as.character(stats::na.omit(y_observed))))
    if (length(pred_num) > n_obs &&
        n_obs > 0L &&
        length(cls_levels) > 1L &&
        length(pred_num) %% n_obs == 0L &&
        (length(pred_num) / n_obs) == length(cls_levels)) {
      y_predicted <- matrix(pred_num, nrow = n_obs, ncol = length(cls_levels), byrow = FALSE)
      colnames(y_predicted) <- cls_levels
    }
  }

  if (!is.matrix(y_predicted) && !is.data.frame(y_predicted) &&
      length(y_observed) != length(y_predicted)) {
    stop("Observed and predicted values must have the same length.")
  }
  if (fam == "gaussian") {
    metric_functions <- list(
      accuracy = function(y, y_hat) cor(y, y_hat),
      mean_squared_error = function(y, y_hat) mean((y - y_hat)^2),
      bias = function(y, y_hat) abs(mean(y - y_hat)),
      root_mean_squared_error = function(y, y_hat) sqrt(mean((y - y_hat)^2)),
      relative_squared_error = function(y, y_hat) sum((y - y_hat)^2) / sum((y - mean(y))^2),
      mean_absolute_error = function(y, y_hat) mean(abs(y - y_hat)),
      mean_absolute_percent_error = function(y, y_hat) mean(abs((y - y_hat) / y)) * 100,
      kendalls_tau = kendalls_tau
    )
  } else {
    y_true <- gp_to_class_labels(y_observed)
    observed_levels <- gp_response_class_levels(y_observed, fam)
    positive_class <- if (identical(fam, "binary")) {
      gp_resolve_positive_class(y_observed, fam, positive_class)
    } else {
      NULL
    }
    y_pred_labels <- gp_pred_to_class_labels(
      y_predicted = y_predicted,
      positive_class = positive_class,
      class_levels = observed_levels
    )
    cls_levels <- observed_levels %||% sort(unique(c(y_true, y_pred_labels)))
    cm <- gp_confusion_from_labels(y_true, y_pred_labels, levels = cls_levels)
    macro <- gp_multiclass_precision_recall_f1(cm)

    metric_functions <- list(
      accuracy = function(y, y_hat) mean(y_true == y_pred_labels),
      balanced_accuracy = function(y, y_hat) gp_balanced_accuracy(cm),
      precision = function(y, y_hat) {
        if (fam != "binary") stop(msg, "precision is only valid for binary response_family.", call. = FALSE)
        pos <- positive_class
        tp <- cm[pos, pos]
        fp <- sum(cm[, pos]) - tp
        gp_safe_divide(tp, tp + fp)
      },
      recall = function(y, y_hat) {
        if (fam != "binary") stop(msg, "recall is only valid for binary response_family.", call. = FALSE)
        pos <- positive_class
        tp <- cm[pos, pos]
        fn <- sum(cm[pos, ]) - tp
        gp_safe_divide(tp, tp + fn)
      },
      specificity = function(y, y_hat) {
        if (fam != "binary") stop(msg, "specificity is only valid for binary response_family.", call. = FALSE)
        pos <- positive_class
        neg <- setdiff(cls_levels, pos)
        tn <- sum(cm[neg, neg, drop = FALSE])
        fp <- sum(cm[neg, pos, drop = FALSE])
        gp_safe_divide(tn, tn + fp)
      },
      f1 = function(y, y_hat) {
        if (fam != "binary") stop(msg, "f1 is only valid for binary response_family.", call. = FALSE)
        pr <- metric_functions$precision(y, y_hat)
        rc <- metric_functions$recall(y, y_hat)
        gp_safe_divide(2 * pr * rc, pr + rc)
      },
      mcc = function(y, y_hat) {
        if (fam != "binary") stop(msg, "mcc is only valid for binary response_family.", call. = FALSE)
        pos <- positive_class
        neg <- setdiff(cls_levels, pos)
        # table() counts are integers; the four-way product below overflows
        # R's integer range for a few hundred predictions, so use doubles.
        tp <- as.numeric(cm[pos, pos])
        tn <- as.numeric(sum(cm[neg, neg, drop = FALSE]))
        fp <- as.numeric(sum(cm[neg, pos, drop = FALSE]))
        fn <- as.numeric(sum(cm[pos, neg, drop = FALSE]))
        den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        gp_safe_divide((tp * tn) - (fp * fn), den)
      },
      brier_score = function(y, y_hat) {
        if (fam != "binary") stop(msg, "brier_score is only valid for binary response_family.", call. = FALSE)
        p <- gp_binary_positive_probability(y_predicted, cls_levels, positive_class)
        y01 <- as.numeric(y_true == positive_class)
        mean((y01 - p)^2)
      },
      log_loss = function(y, y_hat) gp_log_loss(
        y_observed, y_predicted, fam, positive_class = positive_class
      ),
      ece = function(y, y_hat) gp_expected_calibration_error(
        y_observed, y_predicted, fam, positive_class = positive_class
      ),
      macro_precision = function(y, y_hat) {
        if (fam == "binary") stop(msg, "macro_precision is not used for binary response_family.", call. = FALSE)
        macro$precision
      },
      macro_recall = function(y, y_hat) {
        if (fam == "binary") stop(msg, "macro_recall is not used for binary response_family.", call. = FALSE)
        macro$recall
      },
      macro_f1 = function(y, y_hat) {
        if (fam == "binary") stop(msg, "macro_f1 is not used for binary response_family.", call. = FALSE)
        macro$f1
      },
      mean_absolute_error_class = function(y, y_hat) {
        if (fam != "ordinal") stop(msg, "mean_absolute_error_class is only valid for ordinal response_family.", call. = FALSE)
        yt <- match(y_true, cls_levels)
        yp <- match(y_pred_labels, cls_levels)
        mean(abs(yt - yp))
      },
      within_one_class_accuracy = function(y, y_hat) {
        if (fam != "ordinal") stop(msg, "within_one_class_accuracy is only valid for ordinal response_family.", call. = FALSE)
        yt <- match(y_true, cls_levels)
        yp <- match(y_pred_labels, cls_levels)
        mean(abs(yt - yp) <= 1)
      },
      quadratic_weighted_kappa = function(y, y_hat) {
        if (fam != "ordinal") stop(msg, "quadratic_weighted_kappa is only valid for ordinal response_family.", call. = FALSE)
        gp_quadratic_weighted_kappa(y_true, y_pred_labels, levels = cls_levels)
      }
    )
  }

  results <- NULL
  for (metric in tolower(eval_metrics)) {
    if (!metric %in% names(metric_functions)) {
      stop(msg, "Unknown evaluation metric for response_family='", fam, "'.", call. = FALSE)
    }
    results <- metric_functions[[metric]](y_observed, y_predicted)
  }
  results
}

