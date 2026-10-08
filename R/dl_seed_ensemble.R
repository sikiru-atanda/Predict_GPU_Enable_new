#' Plan deep-learning seed ensembles
#'
#' Builds the deterministic training-seed schedule used by PredictProR deep-
#' learning true-prediction workflows and reports the corresponding model-fit
#' count before any model is trained. Bootstrap resampling and training seeds
#' are separate dimensions: each bootstrap sample is fitted once per training
#' seed, and predictions are averaged across seeds within that sample.
#'
#' @param n_bootstrap Positive integer number of bootstrap samples used for
#'   true-prediction uncertainty.
#' @param dl_n_seeds Optional positive integer number of training seeds. When
#'   `dl_seeds` is not supplied, the seeds are deterministically derived from
#'   `random_seed`.
#' @param dl_seeds Optional explicit vector of unique non-negative integer training
#'   seeds. If both `dl_n_seeds` and `dl_seeds` are supplied, their lengths must
#'   agree.
#' @param random_seed Non-negative integer master seed. It controls bootstrap
#'   resampling and is the first derived training seed.
#' @param calibration_folds Non-negative integer number of internal held-out
#'   calibration folds. Use zero when the response family does not run the
#'   Gaussian risk-calibration step.
#' @param dl_seed_aggregation Seed aggregation method. Currently only the
#'   arithmetic mean is supported; no seed is selected using test outcomes.
#'
#' @return A list containing `seed_manifest`, `fit_counts`, and the aggregation
#'   rule.
#' @export
dl_seed_plan <- function(n_bootstrap = 30L,
                         dl_n_seeds = NULL,
                         dl_seeds = NULL,
                         random_seed = 123L,
                         calibration_folds = 5L,
                         dl_seed_aggregation = "mean") {
  n_bootstrap <- gp_dl_positive_count(n_bootstrap, "n_bootstrap")
  calibration_folds <- gp_dl_nonnegative_count(
    calibration_folds,
    "calibration_folds"
  )
  aggregation <- gp_dl_seed_aggregation(dl_seed_aggregation)
  manifest <- gp_dl_seed_manifest(
    dl_n_seeds = dl_n_seeds,
    dl_seeds = dl_seeds,
    random_seed = random_seed
  )
  n_seeds <- nrow(manifest)
  bootstrap_fits <- as.double(n_bootstrap) * as.double(n_seeds)
  calibration_fits <- as.double(calibration_folds) * as.double(n_seeds)

  list(
    seed_manifest = manifest,
    fit_counts = data.frame(
      component = c(
        "bootstrap_true_prediction",
        "heldout_risk_calibration",
        "total"
      ),
      data_resamples_or_folds = c(
        n_bootstrap,
        calibration_folds,
        n_bootstrap + calibration_folds
      ),
      training_seeds_per_resample = rep(n_seeds, 3L),
      model_fits = c(
        bootstrap_fits,
        calibration_fits,
        bootstrap_fits + calibration_fits
      ),
      stringsAsFactors = FALSE
    ),
    aggregation = aggregation,
    selection_rule = "all_requested_seeds_are_averaged; no_test_outcome_selection"
  )
}

gp_dl_integerish_scalar <- function(x, name, allow_zero = FALSE) {
  value <- suppressWarnings(as.double(x))
  lower <- if (isTRUE(allow_zero)) 0 else 1
  if (length(value) != 1L || !is.finite(value) || value < lower ||
      value != floor(value) || value > .Machine$integer.max) {
    qualifier <- if (isTRUE(allow_zero)) "a non-negative" else "a positive"
    stop("`", name, "` must be ", qualifier, " integer.", call. = FALSE)
  }
  as.integer(value)
}

gp_dl_positive_count <- function(x, name) {
  gp_dl_integerish_scalar(x, name, allow_zero = FALSE)
}

gp_dl_calibration_enabled <- function(x) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`dl_internal_calibration` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(x)
}

gp_dl_nonnegative_count <- function(x, name) {
  gp_dl_integerish_scalar(x, name, allow_zero = TRUE)
}

gp_dl_seed_aggregation <- function(x) {
  value <- tolower(trimws(as.character(x %||% "mean")[[1L]]))
  if (!identical(value, "mean")) {
    stop(
      "`dl_seed_aggregation` currently supports only 'mean'; PredictProR does not select a best seed from test outcomes.",
      call. = FALSE
    )
  }
  value
}

gp_dl_bootstrap_response_scale <- function(boot_matrix,
                                           response_family,
                                           y_scaler = NULL) {
  out <- as.matrix(boot_matrix)
  storage.mode(out) <- "double"
  fam <- gp_dl_bridge_response_family(response_family)
  if (identical(fam, "gaussian") && !is.null(y_scaler)) {
    out <- out * as.numeric(y_scaler$std) + as.numeric(y_scaler$mean)
  } else if (identical(fam, "binary")) {
    # Scalar-first pmin()/pmax() drops the matrix dimensions in base R. Keep
    # the bootstrap-by-target layout because the downstream interval and
    # variance summaries operate over target columns.
    out[] <- pmax(0, pmin(1, as.numeric(out)))
  }
  out
}

gp_dl_validate_seed_vector <- function(x, name = "dl_seeds") {
  values <- suppressWarnings(as.double(x))
  if (!length(values) || any(!is.finite(values)) || any(values < 0) ||
      any(values != floor(values)) || any(values > .Machine$integer.max)) {
    stop("`", name, "` must contain non-negative integer seeds.", call. = FALSE)
  }
  values <- as.integer(values)
  if (anyDuplicated(values)) {
    stop("`", name, "` must contain unique seeds.", call. = FALSE)
  }
  values
}

gp_dl_seed_manifest <- function(dl_n_seeds = NULL,
                                dl_seeds = NULL,
                                random_seed = 123L) {
  base_seed <- gp_dl_integerish_scalar(
    random_seed,
    "random_seed",
    allow_zero = TRUE
  )
  if (!is.null(dl_seeds)) {
    seeds <- gp_dl_validate_seed_vector(dl_seeds)
    if (!is.null(dl_n_seeds)) {
      n_requested <- gp_dl_positive_count(dl_n_seeds, "dl_n_seeds")
      if (n_requested != length(seeds)) {
        stop(
          "`dl_n_seeds` must equal length(`dl_seeds`) when both are supplied.",
          call. = FALSE
        )
      }
    }
    source <- "explicit"
  } else {
    n_requested <- if (is.null(dl_n_seeds)) {
      1L
    } else {
      gp_dl_positive_count(dl_n_seeds, "dl_n_seeds")
    }
    # Keep n = 1 exactly backward compatible. Additional seeds use a fixed
    # prime stride and are reported in the manifest, so the derivation is
    # deterministic and inspectable rather than hidden in the backend.
    modulus <- as.double(.Machine$integer.max) + 1
    offsets <- as.double(seq_len(n_requested) - 1L) * 104729
    seeds <- as.integer((as.double(base_seed) + offsets) %% modulus)
    if (anyDuplicated(seeds)) {
      stop(
        "The derived DL seed schedule contains duplicates; provide explicit `dl_seeds`.",
        call. = FALSE
      )
    }
    source <- "derived_from_random_seed"
  }

  data.frame(
    seed_index = seq_along(seeds),
    training_seed = seeds,
    source = rep(source, length(seeds)),
    stringsAsFactors = FALSE
  )
}

gp_dl_seed_requested <- function(dl_n_seeds = NULL, dl_seeds = NULL) {
  (!is.null(dl_n_seeds) && suppressWarnings(as.double(dl_n_seeds[[1L]])) > 1) ||
    (!is.null(dl_seeds) && length(dl_seeds) > 1L)
}

gp_dl_set_job_training_seed <- function(job, seed) {
  job$dl_args <- job$dl_args %||% list()
  job$dl_args$random_seed <- as.integer(seed)
  job
}

gp_dl_aggregate_formatted_predictions <- function(predictions,
                                                   response_family,
                                                   class_levels = NULL) {
  if (!length(predictions)) {
    stop("No DL seed predictions were supplied for aggregation.", call. = FALSE)
  }
  fam <- gp_dl_bridge_response_family(response_family)
  if (length(predictions) == 1L) {
    return(predictions[[1L]])
  }

  if (identical(fam, "gaussian")) {
    values <- do.call(rbind, lapply(predictions, as.numeric))
    return(colMeans(values))
  }
  if (identical(fam, "multitask_regression")) {
    arrays <- lapply(predictions, as.matrix)
    total <- Reduce(`+`, arrays)
    return(total / length(arrays))
  }

  probabilities <- lapply(predictions, function(x) attr(x, "probabilities"))
  if (any(vapply(probabilities, is.null, logical(1)))) {
    stop(
      "Classification DL seed aggregation requires per-seed probabilities.",
      call. = FALSE
    )
  }
  mean_probability <- Reduce(`+`, lapply(probabilities, as.matrix)) /
    length(probabilities)
  row_total <- rowSums(mean_probability)
  row_total[!is.finite(row_total) | row_total <= 0] <- 1
  mean_probability <- mean_probability / row_total

  if (identical(fam, "binary")) {
    out <- as.numeric(mean_probability[, ncol(mean_probability), drop = TRUE])
    attr(out, "probabilities") <- mean_probability
    return(out)
  }

  levels <- as.character(
    class_levels %||% colnames(mean_probability) %||%
      paste0("class_", seq_len(ncol(mean_probability)))
  )
  if (length(levels) != ncol(mean_probability)) {
    stop("DL seed probability columns do not match the class levels.", call. = FALSE)
  }
  colnames(mean_probability) <- levels
  out <- levels[max.col(mean_probability, ties.method = "first")]
  attr(out, "probabilities") <- mean_probability
  out
}

gp_dl_bridge_fit_predict_seed_batch <- function(jobs, seed_manifest) {
  if (is.null(jobs) || !length(jobs)) {
    return(list(predictions = list(), per_seed = list()))
  }
  seeds <- as.integer(seed_manifest$training_seed)
  expanded <- unlist(lapply(seq_along(jobs), function(job_index) {
    lapply(seq_along(seeds), function(seed_index) {
      job <- gp_dl_set_job_training_seed(jobs[[job_index]], seeds[[seed_index]])
      job$id <- paste0(job$id %||% job_index, "__seed_", seed_index)
      job
    })
  }), recursive = FALSE)
  raw_predictions <- gp_dl_bridge_fit_predict_batch(expanded)
  grouped <- split(
    raw_predictions,
    rep(seq_along(jobs), each = length(seeds))
  )
  aggregated <- Map(function(seed_predictions, job) {
    gp_dl_aggregate_formatted_predictions(
      predictions = seed_predictions,
      response_family = job$response_family %||% "gaussian",
      class_levels = job$class_levels %||% NULL
    )
  }, grouped, jobs)
  list(predictions = unname(aggregated), per_seed = unname(grouped))
}

gp_dl_bridge_fit_predict_seed_serial <- function(model_type,
                                                 X_train,
                                                 y_train,
                                                 X_test,
                                                 response_family,
                                                 dl_args,
                                                 class_levels,
                                                 seed_manifest) {
  per_seed <- lapply(seed_manifest$training_seed, function(seed) {
    seeded_args <- dl_args %||% list()
    seeded_args$random_seed <- as.integer(seed)
    gp_dl_bridge_fit_predict(
      model_type = model_type,
      X_train = X_train,
      y_train = y_train,
      X_test = X_test,
      response_family = response_family,
      dl_args = seeded_args,
      class_levels = class_levels
    )
  })
  list(
    prediction = gp_dl_aggregate_formatted_predictions(
      predictions = per_seed,
      response_family = response_family,
      class_levels = class_levels
    ),
    per_seed = per_seed
  )
}

gp_dl_seed_output_tables <- function(boot_payload,
                                     ids,
                                     response_family,
                                     y_scaler = NULL,
                                     class_levels = NULL,
                                     gen_name = "GID") {
  manifest <- boot_payload$seed_manifest %||% data.frame()
  seed_matrix <- boot_payload$seed_predictions %||% NULL
  variability <- boot_payload$seed_variability %||% NULL
  if (is.null(seed_matrix) || !nrow(seed_matrix) || !length(ids)) {
    return(list(
      manifest = manifest,
      predictions = data.frame(),
      variability = data.frame()
    ))
  }

  seed_matrix <- as.matrix(seed_matrix)
  storage.mode(seed_matrix) <- "double"
  fam <- gp_dl_bridge_response_family(response_family)
  if (identical(fam, "gaussian") && !is.null(y_scaler)) {
    seed_matrix <- seed_matrix * y_scaler$std + y_scaler$mean
  }
  seed_values <- if (nrow(manifest) == nrow(seed_matrix)) {
    manifest$training_seed
  } else {
    seq_len(nrow(seed_matrix))
  }
  if (identical(fam, "multiclass")) {
    levels <- as.character(class_levels %||% character())
    if (!length(levels) || ncol(seed_matrix) != length(ids) * length(levels)) {
      stop(
        "Multiclass DL per-seed probabilities do not align with IDs and class levels.",
        call. = FALSE
      )
    }
    prediction_parts <- lapply(seq_len(nrow(seed_matrix)), function(i) {
      probability <- matrix(
        seed_matrix[i, ],
        nrow = length(ids),
        ncol = length(levels)
      )
      colnames(probability) <- paste0(
        "Probability_",
        make.names(levels, unique = TRUE)
      )
      data.frame(
        seed_index = i,
        training_seed = seed_values[[i]],
        id = as.character(ids),
        Predicted_class = levels[max.col(probability, ties.method = "first")],
        as.data.frame(probability, check.names = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    prediction_table <- do.call(rbind, prediction_parts)
  } else {
    if (ncol(seed_matrix) != length(ids)) {
      stop("DL per-seed predictions do not align with the prediction IDs.", call. = FALSE)
    }
    prediction_table <- data.frame(
      seed_index = rep(seq_len(nrow(seed_matrix)), each = ncol(seed_matrix)),
      training_seed = rep(seed_values, each = ncol(seed_matrix)),
      id = rep(as.character(ids), times = nrow(seed_matrix)),
      Predicted_value = as.vector(t(seed_matrix)),
      stringsAsFactors = FALSE
    )
  }
  names(prediction_table)[names(prediction_table) == "id"] <- gen_name

  if (is.null(variability) || !nrow(variability)) {
    variability_table <- data.frame(
      id = as.character(ids),
      DL_seed_prediction_SD = rep(NA_real_, length(ids)),
      DL_seed_prediction_range = rep(NA_real_, length(ids)),
      stringsAsFactors = FALSE
    )
  } else {
    variability <- as.data.frame(variability, stringsAsFactors = FALSE)
    scale_sd <- if (identical(fam, "gaussian") && !is.null(y_scaler)) {
      abs(y_scaler$std)
    } else {
      1
    }
    if (identical(fam, "multiclass")) {
      levels <- as.character(class_levels %||% character())
      if (!length(levels) || nrow(variability) != length(ids) * length(levels)) {
        stop(
          "Multiclass DL seed variability does not align with IDs and class levels.",
          call. = FALSE
        )
      }
      seed_sd <- matrix(
        sqrt(pmax(0, variability$seed_variance)),
        nrow = length(ids),
        ncol = length(levels)
      )
      seed_range <- matrix(
        variability$seed_range,
        nrow = length(ids),
        ncol = length(levels)
      )
      colnames(seed_sd) <- paste0(
        "DL_seed_probability_SD_",
        make.names(levels, unique = TRUE)
      )
      colnames(seed_range) <- paste0(
        "DL_seed_probability_range_",
        make.names(levels, unique = TRUE)
      )
      row_max_or_na <- function(x) {
        apply(x, 1L, function(values) {
          values <- values[is.finite(values)]
          if (length(values)) max(values) else NA_real_
        })
      }
      variability_table <- data.frame(
        id = as.character(ids),
        DL_seed_probability_SD = row_max_or_na(seed_sd),
        DL_seed_probability_range = row_max_or_na(seed_range),
        as.data.frame(seed_sd, check.names = FALSE),
        as.data.frame(seed_range, check.names = FALSE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      if (nrow(variability) != length(ids)) {
        stop("DL seed variability does not align with the prediction IDs.", call. = FALSE)
      }
      variability_table <- data.frame(
        id = as.character(ids),
        DL_seed_prediction_SD = sqrt(pmax(0, variability$seed_variance)) * scale_sd,
        DL_seed_prediction_range = variability$seed_range * scale_sd,
        stringsAsFactors = FALSE
      )
    }
  }
  names(variability_table)[names(variability_table) == "id"] <- gen_name

  list(
    manifest = manifest,
    predictions = prediction_table,
    variability = variability_table
  )
}

# Deep-learning true prediction: one network per training seed, each trained
# on ALL training lines, averaged. Returns the ensemble on the model scale:
# a numeric vector (gaussian: scaled response; binary: positive-class
# probability) or an n x K probability matrix for multiclass.
gp_dl_full_data_seed_ensemble <- function(model_type, X_train, y_train, X_pred,
                                          response_family, dl_args = list(),
                                          seeds = 123L, class_levels = NULL) {
  seeds <- unique(as.integer(seeds %||% 123L))
  jobs <- lapply(seq_along(seeds), function(i) {
    list(
      id = paste0("seed_", seeds[[i]]),
      model_type = model_type,
      X_train = X_train,
      y_train = y_train,
      X_test = X_pred,
      response_family = response_family,
      dl_args = utils::modifyList(dl_args %||% list(), list(random_seed = seeds[[i]])),
      class_levels = class_levels
    )
  })
  results <- gp_dl_bridge_fit_predict_batch_raw(jobs)
  fam <- gp_dl_bridge_response_family(response_family)
  if (identical(fam, "multiclass")) {
    probs <- lapply(results, function(r) {
      p <- as.matrix(r$probabilities)
      cls <- as.character(r$classes %||% class_levels %||% colnames(p))
      if (length(cls) == ncol(p)) colnames(p) <- cls
      p[, as.character(class_levels), drop = FALSE]
    })
    return(Reduce(`+`, probs) / length(probs))
  }
  preds <- vapply(results, function(r) as.numeric(r$predictions), numeric(NROW(X_pred)))
  rowMeans(matrix(preds, nrow = NROW(X_pred)))
}