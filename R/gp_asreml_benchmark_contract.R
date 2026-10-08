gp_asreml_benchmark_scenarios <- function() {
  data.frame(
    scenario = c(
      "single_trait_single_environment",
      "single_trait_met_single_kernel",
      "single_trait_met_multi_kernel",
      "multi_trait_single_environment",
      "multi_trait_multi_environment"
    ),
    model_family = c("GP", "GP_FA", "GP_FA", "MT_GP", "MT_MET_GP"),
    n_traits = c(1L, 1L, 1L, 2L, 2L),
    n_environments = c(1L, 3L, 3L, 1L, 3L),
    n_kernels = c(1L, 1L, 2L, 1L, 1L),
    parity_class = c(
      "matched",
      "matched",
      "nested_asreml_reference",
      "matched",
      "prediction_reference"
    ),
    stringsAsFactors = FALSE
  )
}

gp_asreml_benchmark_prediction_columns <- function() {
  c(
    "scenario", "dataset", "engine", "model_family",
    "GID", "Env", "Trait", "Observed_value", "Predicted_value",
    "Standard_error", "PEV", "Train_Test_Label"
  )
}

gp_asreml_benchmark_pick_column <- function(data, candidates, required = FALSE) {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) {
    if (isTRUE(required)) {
      stop(
        "Benchmark prediction data are missing a required column. Expected one of: ",
        paste(candidates, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(NULL)
  }
  hit[[1L]]
}

gp_asreml_benchmark_standardize_predictions <- function(data,
                                                         scenario,
                                                         dataset,
                                                         engine,
                                                         model_family,
                                                         observed = NULL,
                                                         default_env = NA_character_,
                                                         default_trait = NA_character_) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!nrow(data)) {
    stop("Benchmark prediction data must contain at least one row.", call. = FALSE)
  }

  gid_col <- gp_asreml_benchmark_pick_column(
    data, c("GID", "gid", "Name", "genotype"), required = TRUE
  )
  env_col <- gp_asreml_benchmark_pick_column(data, c("Env", "env", "Environment"))
  trait_col <- gp_asreml_benchmark_pick_column(data, c("Trait", "trait"))
  pred_col <- gp_asreml_benchmark_pick_column(
    data, c("Predicted_value", "Prediction", "predicted.value"), required = TRUE
  )
  obs_col <- gp_asreml_benchmark_pick_column(
    data, c("Observed_value", "y_true", "y", "observed")
  )
  se_col <- gp_asreml_benchmark_pick_column(
    data,
    c("Standard_error", "SE", "Prediction_SE_observed", "std.error")
  )
  pev_col <- gp_asreml_benchmark_pick_column(
    data,
    c(
      "Prediction_Var_observed", "PEV", "Prediction_error_variance",
      "Prediction_Var_latent"
    )
  )

  n <- nrow(data)
  obs <- if (!is.null(observed)) {
    as.numeric(observed)
  } else if (!is.null(obs_col)) {
    suppressWarnings(as.numeric(data[[obs_col]]))
  } else {
    rep(NA_real_, n)
  }
  if (length(obs) != n) {
    stop("`observed` must have one value per prediction row.", call. = FALSE)
  }
  se <- if (is.null(se_col)) rep(NA_real_, n) else suppressWarnings(as.numeric(data[[se_col]]))
  pev <- if (is.null(pev_col)) se^2 else suppressWarnings(as.numeric(data[[pev_col]]))

  out <- data.frame(
    scenario = rep(as.character(scenario)[1L], n),
    dataset = rep(as.character(dataset)[1L], n),
    engine = rep(as.character(engine)[1L], n),
    model_family = rep(as.character(model_family)[1L], n),
    GID = as.character(data[[gid_col]]),
    Env = if (is.null(env_col)) rep(as.character(default_env)[1L], n) else as.character(data[[env_col]]),
    Trait = if (is.null(trait_col)) rep(as.character(default_trait)[1L], n) else as.character(data[[trait_col]]),
    Observed_value = obs,
    Predicted_value = suppressWarnings(as.numeric(data[[pred_col]])),
    Standard_error = se,
    PEV = pev,
    Train_Test_Label = rep("Test", n),
    stringsAsFactors = FALSE
  )
  out[gp_asreml_benchmark_prediction_columns()]
}

gp_validate_asreml_benchmark_predictions <- function(data) {
  errors <- character()
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  expected <- gp_asreml_benchmark_prediction_columns()
  if (!identical(names(data), expected)) {
    errors <- c(errors, paste0("Columns must be exactly: ", paste(expected, collapse = ", "), "."))
  }
  if (nrow(data)) {
    required_text <- intersect(c("scenario", "dataset", "engine", "model_family", "GID"), names(data))
    for (nm in required_text) {
      values <- as.character(data[[nm]])
      if (any(is.na(values) | !nzchar(values))) {
        errors <- c(errors, paste0("`", nm, "` contains missing or empty values."))
      }
    }
    if ("engine" %in% names(data) && any(!data$engine %in% c("PredictProR_GP", "ASReml_R"))) {
      errors <- c(errors, "`engine` must contain only PredictProR_GP or ASReml_R.")
    }
    if ("Predicted_value" %in% names(data) && any(!is.finite(data$Predicted_value))) {
      errors <- c(errors, "`Predicted_value` contains non-finite values.")
    }
    if ("Standard_error" %in% names(data) && any(is.finite(data$Standard_error) & data$Standard_error < 0)) {
      errors <- c(errors, "`Standard_error` contains negative values.")
    }
    if ("PEV" %in% names(data) && any(is.finite(data$PEV) & data$PEV < 0)) {
      errors <- c(errors, "`PEV` contains negative values.")
    }
    if (all(c("Standard_error", "PEV") %in% names(data))) {
      finite_uncertainty <- is.finite(data$Standard_error) & is.finite(data$PEV)
      expected_pev <- data$Standard_error[finite_uncertainty]^2
      tolerance <- 1e-8 * pmax(1, abs(expected_pev), abs(data$PEV[finite_uncertainty]))
      if (any(abs(data$PEV[finite_uncertainty] - expected_pev) > tolerance)) {
        errors <- c(errors, "`PEV` must equal `Standard_error^2` within numerical tolerance.")
      }
    }
    key_cols <- intersect(c("scenario", "engine", "GID", "Env", "Trait"), names(data))
    if (length(key_cols) == 5L && anyDuplicated(data[key_cols])) {
      errors <- c(errors, "Prediction keys are duplicated within scenario and engine.")
    }
  }
  list(valid = !length(errors), errors = unique(errors))
}

gp_asreml_benchmark_covariance_diagnostics <- function(estimate,
                                                       truth = NULL,
                                                       tolerance = 1e-8) {
  estimate <- as.matrix(estimate)
  if (!is.numeric(estimate) || nrow(estimate) != ncol(estimate) || !nrow(estimate)) {
    stop("`estimate` must be a non-empty numeric square matrix.", call. = FALSE)
  }
  symmetry_error <- max(abs(estimate - t(estimate)), na.rm = TRUE)
  eigenvalues <- eigen((estimate + t(estimate)) / 2, symmetric = TRUE, only.values = TRUE)$values
  truth_error <- NA_real_
  if (!is.null(truth)) {
    truth <- as.matrix(truth)
    if (!identical(dim(truth), dim(estimate))) {
      stop("`truth` and `estimate` must have the same dimensions.", call. = FALSE)
    }
    denom <- sqrt(sum(truth^2))
    truth_error <- sqrt(sum((estimate - truth)^2)) / max(denom, .Machine$double.eps)
  }
  data.frame(
    symmetry_error = symmetry_error,
    min_eigenvalue = min(eigenvalues),
    positive_semidefinite = min(eigenvalues) >= -abs(tolerance),
    relative_frobenius_truth_error = truth_error,
    stringsAsFactors = FALSE
  )
}

gp_asreml_benchmark_metrics_by <- function(predictions, group_by = character()) {
  predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
  validation <- gp_validate_asreml_benchmark_predictions(predictions)
  if (!isTRUE(validation$valid)) {
    stop(paste(validation$errors, collapse = " "), call. = FALSE)
  }

  group_by <- unique(as.character(group_by))
  allowed_groups <- c("Env", "Trait")
  unknown_groups <- setdiff(group_by, allowed_groups)
  if (length(unknown_groups)) {
    stop(
      "`group_by` must contain only Env and/or Trait.",
      call. = FALSE
    )
  }
  if (length(group_by)) {
    missing_group <- vapply(group_by, function(nm) {
      values <- as.character(predictions[[nm]])
      any(is.na(values) | !nzchar(values))
    }, logical(1L))
    if (any(missing_group)) {
      stop(
        "Stratified metrics require non-missing values for: ",
        paste(group_by[missing_group], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  identity_columns <- c("scenario", "dataset", "engine", group_by)
  if (!nrow(predictions)) {
    out <- predictions[FALSE, identity_columns, drop = FALSE]
    out$n <- integer()
    out$RMSE <- numeric()
    out$MAE <- numeric()
    out$Pearson <- numeric()
    out$Spearman <- numeric()
    return(out)
  }

  split_key <- do.call(
    interaction,
    c(
      unname(predictions[identity_columns]),
      list(drop = TRUE, lex.order = TRUE)
    )
  )
  rows <- lapply(split(predictions, split_key), function(one) {
    ok <- is.finite(one$Observed_value) & is.finite(one$Predicted_value)
    obs <- one$Observed_value[ok]
    pred <- one$Predicted_value[ok]
    cbind(
      one[1L, identity_columns, drop = FALSE],
      data.frame(
        n = length(obs),
        RMSE = if (length(obs)) sqrt(mean((obs - pred)^2)) else NA_real_,
        MAE = if (length(obs)) mean(abs(obs - pred)) else NA_real_,
        Pearson = if (length(obs) >= 3L) stats::cor(obs, pred, method = "pearson") else NA_real_,
        Spearman = if (length(obs) >= 3L) stats::cor(obs, pred, method = "spearman") else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  ordering <- do.call(order, unname(out[identity_columns]))
  out[ordering, c(identity_columns, "n", "RMSE", "MAE", "Pearson", "Spearman"), drop = FALSE]
}

gp_asreml_benchmark_metrics <- function(predictions) {
  gp_asreml_benchmark_metrics_by(predictions, group_by = character())
}

gp_asreml_benchmark_metric_tables <- function(
    predictions,
    scenarios = gp_asreml_benchmark_scenarios()) {
  predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
  scenarios <- as.data.frame(scenarios, stringsAsFactors = FALSE)
  required_scenario_columns <- c("scenario", "n_traits", "n_environments")
  if (!all(required_scenario_columns %in% names(scenarios))) {
    stop(
      "`scenarios` must contain scenario, n_traits, and n_environments.",
      call. = FALSE
    )
  }
  if (anyDuplicated(scenarios$scenario)) {
    stop("`scenarios$scenario` must be unique.", call. = FALSE)
  }

  validation <- gp_validate_asreml_benchmark_predictions(predictions)
  if (!isTRUE(validation$valid)) {
    stop(paste(validation$errors, collapse = " "), call. = FALSE)
  }
  unknown_scenarios <- setdiff(unique(predictions$scenario), scenarios$scenario)
  if (length(unknown_scenarios)) {
    stop(
      "Prediction scenarios are absent from `scenarios`: ",
      paste(unknown_scenarios, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  environment_scenarios <- scenarios$scenario[
    scenarios$n_environments > 1L & scenarios$n_traits == 1L
  ]
  trait_scenarios <- scenarios$scenario[scenarios$n_traits > 1L]
  environment_trait_scenarios <- scenarios$scenario[
    scenarios$n_environments > 1L & scenarios$n_traits > 1L
  ]

  overall <- gp_asreml_benchmark_metrics(predictions)
  by_environment <- gp_asreml_benchmark_metrics_by(
    predictions[predictions$scenario %in% environment_scenarios, , drop = FALSE],
    "Env"
  )
  by_trait <- gp_asreml_benchmark_metrics_by(
    predictions[predictions$scenario %in% trait_scenarios, , drop = FALSE],
    "Trait"
  )
  by_environment_trait <- gp_asreml_benchmark_metrics_by(
    predictions[predictions$scenario %in% environment_trait_scenarios, , drop = FALSE],
    c("Env", "Trait")
  )

  add_scope <- function(data, metric_scope) {
    if (!"Env" %in% names(data)) data$Env <- rep(NA_character_, nrow(data))
    if (!"Trait" %in% names(data)) data$Trait <- rep(NA_character_, nrow(data))
    data$metric_scope <- rep(metric_scope, nrow(data))
    data[c(
      "scenario", "dataset", "engine", "metric_scope", "Env", "Trait",
      "n", "RMSE", "MAE", "Pearson", "Spearman"
    )]
  }
  stratified <- do.call(rbind, list(
    add_scope(overall, "overall"),
    add_scope(by_environment, "environment"),
    add_scope(by_trait, "trait"),
    add_scope(by_environment_trait, "environment_trait")
  ))
  rownames(stratified) <- NULL

  list(
    overall = overall,
    by_environment = by_environment,
    by_trait = by_trait,
    by_environment_trait = by_environment_trait,
    stratified = stratified
  )
}
