#' Validate predictive uncertainty on held-out observations
#'
#' Audits prediction-error calibration and prediction-interval coverage using
#' observations that were not used to fit the evaluated predictions. For
#' machine- and deep-learning outputs, `PEV` is treated as an estimated mean
#' squared prediction error, not as mixed-model breeding-value PEV. Genetic
#' reliability cannot be validated from phenotype predictions alone.
#'
#' @param predictions A prediction data frame containing observed and predicted
#'   values. By default, rows labelled `Test`, `Validation`, `Holdout`, or
#'   `Heldout` in `Train_Test_Label` are evaluated.
#' @param observed_col,predicted_col,pev_col,lower_col,upper_col Column names.
#' @param confidence_level Nominal prediction-interval coverage.
#' @param group_cols Optional columns for subgroup diagnostics.
#' @param evaluation_rows Optional logical vector or integer row indices defining
#'   an independently held-out evaluation set.
#' @param require_heldout If `TRUE`, fail when held-out rows cannot be identified.
#'
#' @return A list with `overall`, optional `by_group`, and evaluated row indices.
#'   `pev_to_mse_ratio` and `mean_squared_standardized_error` should be near one
#'   for a well-calibrated prediction-error estimate. Interval coverage should be
#'   interpreted with its Wilson binomial confidence limits.
#' @export
validate_predictive_uncertainty <- function(
    predictions,
    observed_col = "Observed_value",
    predicted_col = "Predicted_value",
    pev_col = "PEV",
    lower_col = "lower_bound",
    upper_col = "upper_bound",
    confidence_level = 0.95,
    group_cols = NULL,
    evaluation_rows = NULL,
    require_heldout = TRUE) {
  if (!is.data.frame(predictions)) {
    stop("`predictions` must be a data frame.", call. = FALSE)
  }
  required <- c(observed_col, predicted_col)
  missing_required <- setdiff(required, names(predictions))
  if (length(missing_required)) {
    stop(
      "Missing required prediction columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  confidence_level <- suppressWarnings(as.numeric(confidence_level))[1L]
  if (!is.finite(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    stop("`confidence_level` must be strictly between zero and one.", call. = FALSE)
  }

  n_rows <- nrow(predictions)
  if (is.null(evaluation_rows)) {
    if ("Train_Test_Label" %in% names(predictions)) {
      label <- tolower(trimws(as.character(predictions[["Train_Test_Label"]])))
      evaluation_rows <- label %in% c(
        "test", "validation", "validate", "holdout", "heldout", "held-out"
      )
    } else {
      evaluation_rows <- rep(FALSE, n_rows)
    }
  } else if (is.numeric(evaluation_rows)) {
    idx <- suppressWarnings(as.integer(evaluation_rows))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= n_rows]
    keep <- rep(FALSE, n_rows)
    keep[idx] <- TRUE
    evaluation_rows <- keep
  } else {
    evaluation_rows <- as.logical(evaluation_rows)
    if (length(evaluation_rows) != n_rows) {
      stop("Logical `evaluation_rows` must have one value per prediction row.", call. = FALSE)
    }
    evaluation_rows[is.na(evaluation_rows)] <- FALSE
  }

  if (!any(evaluation_rows)) {
    if (isTRUE(require_heldout)) {
      stop(
        "No held-out evaluation rows were identified. Supply `evaluation_rows`, ",
        "or label independent rows as Test/Validation/Heldout.",
        call. = FALSE
      )
    }
    evaluation_rows <- rep(TRUE, n_rows)
  }

  group_cols <- unique(as.character(group_cols %||% character()))
  missing_groups <- setdiff(group_cols, names(predictions))
  if (length(missing_groups)) {
    stop("Missing grouping columns: ", paste(missing_groups, collapse = ", "), call. = FALSE)
  }

  dat <- predictions[evaluation_rows, , drop = FALSE]
  dat[[".observed"]] <- suppressWarnings(as.numeric(dat[[observed_col]]))
  dat[[".predicted"]] <- suppressWarnings(as.numeric(dat[[predicted_col]]))
  dat <- dat[is.finite(dat[[".observed"]]) & is.finite(dat[[".predicted"]]), , drop = FALSE]
  if (!nrow(dat)) {
    stop("Held-out rows contain no finite observed/predicted pairs.", call. = FALSE)
  }

  wilson_limits <- function(successes, n, level = 0.95) {
    if (!is.finite(n) || n < 1L) return(c(NA_real_, NA_real_))
    z <- stats::qnorm(1 - (1 - level) / 2)
    p <- successes / n
    denominator <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denominator
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
    pmax(0, pmin(1, c(center - half, center + half)))
  }

  summarize_one <- function(x) {
    error <- x[[".observed"]] - x[[".predicted"]]
    squared_error <- error^2
    mse <- mean(squared_error)
    rmse <- sqrt(mse)

    pev <- if (pev_col %in% names(x)) suppressWarnings(as.numeric(x[[pev_col]])) else rep(NA_real_, nrow(x))
    valid_pev <- is.finite(pev) & pev > 0
    mean_pev <- if (any(valid_pev)) mean(pev[valid_pev]) else NA_real_
    pev_to_mse <- if (is.finite(mean_pev) && is.finite(mse) && mse > 0) mean_pev / mse else NA_real_
    standardized <- if (any(valid_pev)) mean(squared_error[valid_pev] / pev[valid_pev]) else NA_real_

    lower <- if (lower_col %in% names(x)) suppressWarnings(as.numeric(x[[lower_col]])) else rep(NA_real_, nrow(x))
    upper <- if (upper_col %in% names(x)) suppressWarnings(as.numeric(x[[upper_col]])) else rep(NA_real_, nrow(x))
    valid_interval <- is.finite(lower) & is.finite(upper) & lower <= upper
    interval_n <- sum(valid_interval)
    covered <- valid_interval & x[[".observed"]] >= lower & x[[".observed"]] <= upper
    coverage <- if (interval_n) sum(covered) / interval_n else NA_real_
    coverage_ci <- if (interval_n) wilson_limits(sum(covered), interval_n) else c(NA_real_, NA_real_)
    mean_width <- if (interval_n) mean(upper[valid_interval] - lower[valid_interval]) else NA_real_

    stability_correlation <- NA_real_
    if ("Prediction_stability" %in% names(x)) {
      stability <- suppressWarnings(as.numeric(x[["Prediction_stability"]]))
      valid_stability <- is.finite(stability) & is.finite(squared_error)
      if (sum(valid_stability) >= 3L && length(unique(stability[valid_stability])) > 1L) {
        stability_correlation <- suppressWarnings(stats::cor(
          stability[valid_stability],
          squared_error[valid_stability],
          method = "spearman"
        ))
      }
    }

    data.frame(
      n = nrow(x),
      mean_error = mean(error),
      mse = mse,
      rmse = rmse,
      mean_pev = mean_pev,
      pev_to_mse_ratio = pev_to_mse,
      mean_squared_standardized_error = standardized,
      interval_n = interval_n,
      nominal_coverage = confidence_level,
      empirical_coverage = coverage,
      coverage_gap = coverage - confidence_level,
      coverage_ci_lower = coverage_ci[[1L]],
      coverage_ci_upper = coverage_ci[[2L]],
      mean_interval_width = mean_width,
      stability_error_spearman = stability_correlation,
      stringsAsFactors = FALSE
    )
  }

  overall <- summarize_one(dat)
  by_group <- NULL
  if (length(group_cols)) {
    key <- interaction(dat[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    pieces <- lapply(split(dat, key, drop = TRUE), function(x) {
      ans <- summarize_one(x)
      for (nm in group_cols) ans[[nm]] <- as.character(x[[nm]][[1L]])
      ans[, c(group_cols, setdiff(names(ans), group_cols)), drop = FALSE]
    })
    by_group <- do.call(rbind, pieces)
    rownames(by_group) <- NULL
  }

  structure(
    list(
      overall = overall,
      by_group = by_group,
      evaluation_rows = which(evaluation_rows),
      claim_boundary = paste(
        "ML/DL PEV is assessed as predictive MSE.",
        "This function does not validate genetic variance or breeding-value reliability."
      )
    ),
    class = "PredictProR_uncertainty_validation"
  )
}
