gp_tuning_default_weight_grid <- function() {
  data.frame(
    w_g = c(0.70, 0.80, 0.50, 1.00),
    w_ge = c(0.20, 0.10, 0.35, 0.00),
    w_e = c(0.30, 0.20, 0.30, 0.00),
    stringsAsFactors = FALSE
  )
}

gp_tuning_as_weight_grid <- function(tuning_weight_grid = NULL) {
  if (is.null(tuning_weight_grid)) {
    return(gp_tuning_default_weight_grid())
  }
  tuning_weight_grid <- as.data.frame(tuning_weight_grid, stringsAsFactors = FALSE)
  needed <- c("w_g", "w_ge", "w_e")
  if (!all(needed %in% names(tuning_weight_grid))) {
    stop(
      "`tuning_weight_grid` must contain columns w_g, w_ge, and w_e.",
      call. = FALSE
    )
  }
  tuning_weight_grid <- tuning_weight_grid[needed]
  tuning_weight_grid[] <- lapply(tuning_weight_grid, as.numeric)
  ok <- stats::complete.cases(tuning_weight_grid)
  tuning_weight_grid <- tuning_weight_grid[ok, , drop = FALSE]
  if (!nrow(tuning_weight_grid)) {
    stop("`tuning_weight_grid` does not contain any finite rows.", call. = FALSE)
  }
  tuning_weight_grid
}

gp_tuning_expand_grid <- function(tuning_krr_lams = c(0.03, 0.1, 0.3),
                                  tuning_kenv_kernels = c("matern32", "matern52", "rbf"),
                                  tuning_weight_grid = NULL,
                                  mean_modes = c("none", "location_offset")) {
  lambdas <- unique(as.numeric(tuning_krr_lams))
  lambdas <- lambdas[is.finite(lambdas) & lambdas > 0]
  if (!length(lambdas)) {
    stop("`tuning_krr_lams` must contain at least one positive finite value.", call. = FALSE)
  }
  kernels <- unique(tolower(as.character(tuning_kenv_kernels)))
  kernels <- kernels[nzchar(kernels)]
  valid_kernels <- c("linear", "rbf", "matern32", "matern52")
  bad <- setdiff(kernels, valid_kernels)
  if (length(bad)) {
    stop(
      "`tuning_kenv_kernels` contains unsupported values: ",
      paste(bad, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  mean_modes <- unique(tolower(as.character(mean_modes)))
  mean_modes <- mean_modes[mean_modes %in% c("none", "location_offset")]
  if (!length(mean_modes)) {
    mean_modes <- "none"
  }
  weights <- gp_tuning_as_weight_grid(tuning_weight_grid)
  grid <- expand.grid(
    lambda = lambdas,
    kenv_kernel = kernels,
    mean_mode = mean_modes,
    weight_row = seq_len(nrow(weights)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(
    grid[c("lambda", "kenv_kernel", "mean_mode")],
    weights[grid$weight_row, , drop = FALSE]
  )
  rownames(out) <- NULL
  out
}

gp_tuning_effective_reaction_norm_grid <- function(tune_grid,
                                                   w_ge_default = 0.05,
                                                   w_e_default = 0.5) {
  if (is.null(tune_grid) || !nrow(tune_grid)) {
    return(tune_grid)
  }
  out <- as.data.frame(tune_grid, stringsAsFactors = FALSE)
  if ("w_ge" %in% names(out)) {
    w_ge <- suppressWarnings(as.numeric(out$w_ge))
    zero_ge <- is.finite(w_ge) & abs(w_ge) <= .Machine$double.eps
    out$w_ge[zero_ge] <- as.numeric(w_ge_default)[1L]
  }
  if ("w_e" %in% names(out)) {
    w_e <- suppressWarnings(as.numeric(out$w_e))
    zero_e <- is.finite(w_e) & abs(w_e) <= .Machine$double.eps
    out$w_e[zero_e] <- as.numeric(w_e_default)[1L]
  }
  out
}

gp_tuning_candidate_complexity <- function(candidate) {
  kernel_score <- c(linear = 1, matern32 = 2, matern52 = 3, rbf = 4)
  kernel <- as.character(candidate$kenv_kernel)[1L]
  mode <- as.character(candidate$mean_mode)[1L]
  unname(kernel_score[[kernel]] %||% 5) +
    if (identical(mode, "location_offset")) 1 else 0
}

gp_tuning_normalize_objective <- function(tuning_objective = "rmse") {
  objective <- tolower(as.character(tuning_objective %||% "rmse")[1L])
  objective <- gsub("[^a-z0-9]+", "_", objective, perl = TRUE)
  objective <- gsub("^_+|_+$", "", objective, perl = TRUE)
  aliases <- c(
    raw_rmse = "rmse",
    overall_rmse = "rmse",
    validation_rmse = "rmse",
    centered = "centered_rmse",
    centered_rmse = "centered_rmse",
    within_env_centered_rmse = "centered_rmse",
    within_environment_centered_rmse = "centered_rmse",
    centered_cor = "centered_cor",
    centered_correlation = "centered_cor",
    within_env_centered_cor = "centered_cor",
    rank = "within_env_spearman",
    spearman = "within_env_spearman",
    rank_spearman = "within_env_spearman",
    within_env_rank = "within_env_spearman",
    within_env_spearman = "within_env_spearman",
    within_environment_spearman = "within_env_spearman"
  )
  if (objective %in% names(aliases)) {
    objective <- unname(aliases[[objective]])
  }
  choices <- c("rmse", "centered_rmse", "within_env_spearman", "centered_cor")
  if (!objective %in% choices) {
    stop(
      "`tuning_objective` must be one of: ",
      paste(choices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  objective
}

gp_tuning_objective_spec <- function(tuning_objective = "rmse") {
  objective <- gp_tuning_normalize_objective(tuning_objective)
  switch(
    objective,
    rmse = list(
      name = "rmse",
      column = "validation_rmse",
      se_column = "validation_rmse_se",
      maximize = FALSE
    ),
    centered_rmse = list(
      name = "centered_rmse",
      column = "validation_centered_rmse",
      se_column = "validation_centered_rmse_se",
      maximize = FALSE
    ),
    within_env_spearman = list(
      name = "within_env_spearman",
      column = "validation_within_env_spearman",
      se_column = NA_character_,
      maximize = TRUE
    ),
    centered_cor = list(
      name = "centered_cor",
      column = "validation_centered_cor",
      se_column = NA_character_,
      maximize = TRUE
    )
  )
}

gp_tuning_safe_cor <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L || stats::sd(x) <= 0 || stats::sd(y) <= 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method))
}

gp_tuning_score_vectors <- function(prediction, observed, env = NULL) {
  prediction <- as.numeric(prediction)
  observed <- as.numeric(observed)
  keep <- is.finite(prediction) & is.finite(observed)
  if (!is.null(env)) {
    env <- as.character(env)[keep]
  }
  prediction <- prediction[keep]
  observed <- observed[keep]
  if (!length(prediction)) {
    return(list(
      rmse = NA_real_,
      mae = NA_real_,
      cor = NA_real_,
      centered_rmse = NA_real_,
      centered_mae = NA_real_,
      centered_cor = NA_real_,
      within_env_spearman = NA_real_
    ))
  }
  centered_rmse <- NA_real_
  centered_mae <- NA_real_
  centered_cor <- NA_real_
  within_env_spearman <- NA_real_
  if (!is.null(env) && length(env) == length(prediction)) {
    by_env <- split(seq_along(prediction), env)
    env_n <- vapply(by_env, length, integer(1L))
    pred_centered <- prediction - ave(prediction, env, FUN = function(x) mean(x, na.rm = TRUE))
    obs_centered <- observed - ave(observed, env, FUN = function(x) mean(x, na.rm = TRUE))
    centered_eligible <- env_n[env] >= 2L
    if (any(centered_eligible)) {
      centered_rmse <- sqrt(mean((pred_centered[centered_eligible] - obs_centered[centered_eligible])^2))
      centered_mae <- mean(abs(pred_centered[centered_eligible] - obs_centered[centered_eligible]))
      centered_cor <- gp_tuning_safe_cor(pred_centered[centered_eligible], obs_centered[centered_eligible])
    }
    env_spearman <- vapply(by_env, function(idx) {
      if (length(idx) < 2L) {
        return(NA_real_)
      }
      gp_tuning_safe_cor(prediction[idx], observed[idx], method = "spearman")
    }, numeric(1L))
    ok <- is.finite(env_spearman)
    if (any(ok)) {
      within_env_spearman <- stats::weighted.mean(env_spearman[ok], env_n[ok])
    }
  }
  list(
    rmse = sqrt(mean((prediction - observed)^2)),
    mae = mean(abs(prediction - observed)),
    cor = gp_tuning_safe_cor(prediction, observed),
    centered_rmse = centered_rmse,
    centered_mae = centered_mae,
    centered_cor = centered_cor,
    within_env_spearman = within_env_spearman
  )
}

gp_tuning_blend_objective_values <- function(alpha_grid,
                                             gp_prediction,
                                             baseline_prediction,
                                             observed,
                                             env = NULL,
                                             objective = "rmse") {
  objective <- gp_tuning_normalize_objective(objective)
  alpha_grid <- as.numeric(alpha_grid)
  out <- rep(NA_real_, length(alpha_grid))
  if (!length(alpha_grid)) {
    return(out)
  }

  gp_prediction <- as.numeric(gp_prediction)
  baseline_prediction <- as.numeric(baseline_prediction)
  observed <- as.numeric(observed)
  if (!length(observed)) {
    return(out)
  }

  finite_alpha <- is.finite(alpha_grid)
  if (!any(finite_alpha)) {
    return(out)
  }

  if (identical(objective, "rmse")) {
    out[finite_alpha] <- vapply(alpha_grid[finite_alpha], function(alpha) {
      pred <- alpha * gp_prediction + (1 - alpha) * baseline_prediction
      sqrt(mean((pred - observed)^2))
    }, numeric(1L))
    return(out)
  }

  if (is.null(env) || length(env) != length(observed)) {
    return(out)
  }
  env <- as.character(env)
  by_env <- split(seq_along(observed), env)
  env_n <- vapply(by_env, length, integer(1L))
  centered_eligible <- env_n[env] >= 2L
  centered_eligible[is.na(centered_eligible)] <- FALSE
  if (!any(centered_eligible)) {
    return(out)
  }

  gp_centered <- gp_prediction -
    ave(gp_prediction, env, FUN = function(x) mean(x, na.rm = TRUE))
  baseline_centered <- baseline_prediction -
    ave(baseline_prediction, env, FUN = function(x) mean(x, na.rm = TRUE))
  observed_centered <- observed -
    ave(observed, env, FUN = function(x) mean(x, na.rm = TRUE))

  if (identical(objective, "centered_rmse")) {
    idx <- centered_eligible
    out[finite_alpha] <- vapply(alpha_grid[finite_alpha], function(alpha) {
      pred <- alpha * gp_centered[idx] + (1 - alpha) * baseline_centered[idx]
      sqrt(mean((pred - observed_centered[idx])^2))
    }, numeric(1L))
    return(out)
  }

  if (identical(objective, "centered_cor")) {
    idx <- centered_eligible
    out[finite_alpha] <- vapply(alpha_grid[finite_alpha], function(alpha) {
      pred <- alpha * gp_centered[idx] + (1 - alpha) * baseline_centered[idx]
      gp_tuning_safe_cor(pred, observed_centered[idx])
    }, numeric(1L))
    return(out)
  }

  if (identical(objective, "within_env_spearman")) {
    out[finite_alpha] <- vapply(alpha_grid[finite_alpha], function(alpha) {
      pred <- alpha * gp_prediction + (1 - alpha) * baseline_prediction
      gp_tuning_score_vectors(pred, observed, env)$within_env_spearman
    }, numeric(1L))
  }

  out
}

gp_tuning_rows_by_candidate_id <- function(x, id_col = "candidate_id") {
  if (is.null(x) || !is.data.frame(x) || !nrow(x) || !id_col %in% names(x)) {
    return(list())
  }
  ids <- suppressWarnings(as.integer(x[[id_col]]))
  keep <- is.finite(ids) & !is.na(ids)
  if (!any(keep)) {
    return(list())
  }
  split(x[keep, , drop = FALSE], as.character(ids[keep]), drop = TRUE)
}

gp_tuning_blend_alpha_grid <- function() {
  raw <- Sys.getenv("PREDICTPRO_GP_BLEND_ALPHA_GRID", unset = "")
  if (nzchar(raw)) {
    vals <- suppressWarnings(as.numeric(strsplit(raw, ",", fixed = TRUE)[[1L]]))
    vals <- vals[is.finite(vals)]
  } else {
    vals <- seq(0, 1, by = 0.1)
  }
  vals <- sort(unique(pmin(1, pmax(0, vals))))
  if (length(vals)) vals else seq(0, 1, by = 0.1)
}

gp_prediction_blend_auto_limits <- function() {
  min_improvement <- suppressWarnings(as.numeric(Sys.getenv(
    "PREDICTPRO_GP_BLEND_MIN_IMPROVEMENT",
    unset = "0.05"
  ))[1L])
  if (!is.finite(min_improvement) || is.na(min_improvement) || min_improvement < 0) {
    min_improvement <- 0.05
  }
  max_auto_alpha <- suppressWarnings(as.numeric(Sys.getenv(
    "PREDICTPRO_GP_BLEND_MAX_AUTO_ALPHA",
    unset = "0.10"
  ))[1L])
  if (!is.finite(max_auto_alpha) || is.na(max_auto_alpha)) {
    max_auto_alpha <- 0.10
  }
  list(
    min_relative_improvement = min_improvement,
    max_auto_alpha = pmin(1, pmax(0, max_auto_alpha))
  )
}

gp_prediction_blend_auto_decision <- function(prediction_blend,
                                              blend_alpha,
                                              relative_improvement,
                                              allow_auto = TRUE) {
  limits <- gp_prediction_blend_auto_limits()
  prediction_blend <- as.character(prediction_blend %||% "auto")[1L]
  if (identical(prediction_blend, "validation")) {
    return(c(
      list(apply = TRUE, reason = "validation_blend_forced"),
      limits
    ))
  }
  if (!isTRUE(allow_auto)) {
    return(c(
      list(apply = FALSE, reason = "validation_blend_auto_disabled_for_known_environment"),
      limits
    ))
  }
  alpha_ok <- is.finite(blend_alpha) &&
    blend_alpha <= limits$max_auto_alpha + sqrt(.Machine$double.eps)
  improvement_ok <- is.finite(relative_improvement) &&
    relative_improvement >= limits$min_relative_improvement
  if (alpha_ok && improvement_ok) {
    return(c(
      list(apply = TRUE, reason = "validation_blend_auto_baseline_dominant"),
      limits
    ))
  }
  reason <- if (!improvement_ok) {
    "validation_blend_improvement_below_threshold"
  } else {
    "validation_blend_alpha_above_auto_limit"
  }
  c(
    list(apply = FALSE, reason = reason),
    limits
  )
}

gp_tuning_select_blend <- function(gp_prediction,
                                   baseline_prediction,
                                   observed,
                                   env = NULL,
                                   objective = "rmse",
                                   alpha_grid = gp_tuning_blend_alpha_grid()) {
  objective <- gp_tuning_normalize_objective(objective)
  gp_prediction <- as.numeric(gp_prediction)
  baseline_prediction <- as.numeric(baseline_prediction)
  observed <- as.numeric(observed)
  keep <- is.finite(gp_prediction) & is.finite(baseline_prediction) & is.finite(observed)
  if (!any(keep)) {
    gp_metrics <- gp_tuning_score_vectors(gp_prediction, observed, env)
    return(list(
      blend_alpha = 1,
      blend_objective = objective,
      blend_objective_value = NA_real_,
      blend_rmse = gp_metrics$rmse,
      blend_cor = gp_metrics$cor,
      blend_centered_rmse = gp_metrics$centered_rmse,
      blend_centered_cor = gp_metrics$centered_cor,
      blend_within_env_spearman = gp_metrics$within_env_spearman,
      baseline_rmse = NA_real_,
      baseline_cor = NA_real_,
      baseline_centered_rmse = NA_real_,
      baseline_centered_cor = NA_real_,
      baseline_within_env_spearman = NA_real_
    ))
  }
  env_keep <- if (!is.null(env)) as.character(env)[keep] else NULL
  gp_prediction <- gp_prediction[keep]
  baseline_prediction <- baseline_prediction[keep]
  observed <- observed[keep]
  gp_metrics <- gp_tuning_score_vectors(gp_prediction, observed, env_keep)
  baseline_metrics <- gp_tuning_score_vectors(baseline_prediction, observed, env_keep)
  spec <- gp_tuning_objective_spec(objective)
  values <- gp_tuning_blend_objective_values(
    alpha_grid = alpha_grid,
    gp_prediction = gp_prediction,
    baseline_prediction = baseline_prediction,
    observed = observed,
    env = env_keep,
    objective = objective
  )
  ok <- is.finite(values)
  if (!any(ok)) {
    alpha <- 1
    metrics <- gp_metrics
    value <- switch(
      objective,
      rmse = metrics$rmse,
      centered_rmse = metrics$centered_rmse,
      within_env_spearman = metrics$within_env_spearman,
      centered_cor = metrics$centered_cor
    )
  } else {
    candidate_idx <- which(ok)
    best_local <- if (isTRUE(spec$maximize)) {
      candidate_idx[which.max(values[ok])]
    } else {
      candidate_idx[which.min(values[ok])]
    }
    alpha <- as.numeric(alpha_grid[[best_local]])
    pred <- alpha * gp_prediction + (1 - alpha) * baseline_prediction
    metrics <- gp_tuning_score_vectors(pred, observed, env_keep)
    value <- values[[best_local]]
  }
  list(
    blend_alpha = alpha,
    blend_objective = objective,
    blend_objective_value = as.numeric(value),
    blend_rmse = metrics$rmse,
    blend_cor = metrics$cor,
    blend_centered_rmse = metrics$centered_rmse,
    blend_centered_cor = metrics$centered_cor,
    blend_within_env_spearman = metrics$within_env_spearman,
    baseline_rmse = baseline_metrics$rmse,
    baseline_cor = baseline_metrics$cor,
    baseline_centered_rmse = baseline_metrics$centered_rmse,
    baseline_centered_cor = baseline_metrics$centered_cor,
    baseline_within_env_spearman = baseline_metrics$within_env_spearman
  )
}

gp_tuning_select_candidate <- function(tuning_detail, tuning_objective = "rmse") {
  requested_spec <- gp_tuning_objective_spec(tuning_objective)
  if (!"failed" %in% names(tuning_detail)) {
    tuning_detail$failed <- FALSE
  }
  failed <- as.logical(tuning_detail$failed)
  failed[is.na(failed)] <- FALSE
  objective_chain <- switch(
    requested_spec$name,
    rmse = c("rmse", "centered_rmse", "within_env_spearman", "centered_cor"),
    centered_rmse = c("centered_rmse", "rmse", "within_env_spearman", "centered_cor"),
    within_env_spearman = c("within_env_spearman", "centered_cor", "centered_rmse", "rmse"),
    centered_cor = c("centered_cor", "within_env_spearman", "centered_rmse", "rmse")
  )
  missing_requested <- !requested_spec$column %in% names(tuning_detail)
  if (isTRUE(missing_requested)) {
    stop("Tuning detail is missing objective column: ", requested_spec$column, call. = FALSE)
  }

  last_detail <- tuning_detail
  for (objective in objective_chain) {
    spec <- gp_tuning_objective_spec(objective)
    if (!spec$column %in% names(tuning_detail)) {
      next
    }
    detail <- tuning_detail
    detail$requested_tuning_objective <- requested_spec$name
    detail$tuning_objective <- spec$name
    detail$validation_objective <- as.numeric(detail[[spec$column]])
    detail$validation_objective_loss <- if (isTRUE(spec$maximize)) {
      -detail$validation_objective
    } else {
      detail$validation_objective
    }
    last_detail <- detail
    ok <- detail[
      is.finite(detail$validation_objective) & !failed,
      ,
      drop = FALSE
    ]
    if (!nrow(ok)) {
      next
    }

    objective_values <- ok$validation_objective
    best_idx <- if (isTRUE(spec$maximize)) {
      which.max(objective_values)
    } else {
      which.min(objective_values)
    }
    best_value <- objective_values[[best_idx]]
    se <- NA_real_
    if (!is.na(spec$se_column) && spec$se_column %in% names(ok)) {
      se <- as.numeric(ok[[spec$se_column]])[[best_idx]]
    }
    tolerance <- if (is.finite(se) && se > 0) {
      se
    } else if (isTRUE(spec$maximize)) {
      .Machine$double.eps
    } else {
      max(abs(best_value) * 0.01, .Machine$double.eps)
    }
    eligible <- if (isTRUE(spec$maximize)) {
      ok[ok$validation_objective >= best_value - tolerance, , drop = FALSE]
    } else {
      ok[ok$validation_objective <= best_value + tolerance, , drop = FALSE]
    }
    objective_order <- if (isTRUE(spec$maximize)) {
      -eligible$validation_objective
    } else {
      eligible$validation_objective
    }
    complexity_order <- if ("complexity" %in% names(eligible)) {
      suppressWarnings(as.numeric(eligible$complexity))
    } else {
      rep(0, nrow(eligible))
    }
    complexity_order[!is.finite(complexity_order)] <- Inf
    rmse_order <- if ("validation_rmse" %in% names(eligible)) {
      suppressWarnings(as.numeric(eligible$validation_rmse))
    } else {
      rep(Inf, nrow(eligible))
    }
    rmse_order[!is.finite(rmse_order)] <- Inf
    eligible <- eligible[
      order(complexity_order, objective_order, rmse_order),
      ,
      drop = FALSE
    ]
    fallback_reason <- if (identical(spec$name, requested_spec$name)) {
      NA_character_
    } else {
      paste0(
        "requested objective ", requested_spec$name,
        " had no finite candidates; selected by ", spec$name
      )
    }
    return(list(
      selected = eligible[1L, , drop = FALSE],
      tuning_detail = detail,
      objective = spec,
      requested_objective = requested_spec,
      best_value = best_value,
      tolerance = tolerance,
      reason = fallback_reason
    ))
  }

  list(
    selected = NULL,
    tuning_detail = last_detail,
    objective = requested_spec,
    reason = paste0("no tuning candidates had finite ", requested_spec$name, " objective")
  )
}

gp_tuning_parse_year <- function(env) {
  env <- as.character(env)
  hits <- gregexpr("(19|20)[0-9]{2}", env, perl = TRUE)
  out <- rep(NA_integer_, length(env))
  for (i in seq_along(env)) {
    m <- regmatches(env[i], hits[i])[[1L]]
    if (length(m) && !identical(m, "-1")) {
      out[[i]] <- suppressWarnings(as.integer(tail(m, 1L)))
    }
  }
  out
}

gp_tuning_infer_location <- function(env) {
  env <- as.character(env)
  loc <- gsub("(^|[_ .-])((19|20)[0-9]{2})([_ .-]|$)", "\\1", env, perl = TRUE)
  loc <- gsub("[_ .-]+", "_", loc, perl = TRUE)
  loc <- gsub("^[_ .-]+|[_ .-]+$", "", loc, perl = TRUE)
  loc[!nzchar(loc)] <- env[!nzchar(loc)]
  loc
}

gp_tuning_env_annotation <- function(value,
                                     env_levels,
                                     pheno_df = NULL,
                                     env_col = "Env",
                                     value_name = NULL) {
  env_levels <- as.character(env_levels)
  if (is.null(value)) {
    return(NULL)
  }
  if (is.character(value) && length(value) == 1L &&
      !is.null(pheno_df) && value %in% names(pheno_df)) {
    ann <- data.frame(
      Env = as.character(pheno_df[[env_col]]),
      value = pheno_df[[value]],
      stringsAsFactors = FALSE
    )
    ann <- ann[!is.na(ann$value) & !duplicated(ann$Env), , drop = FALSE]
    return(ann$value[match(env_levels, ann$Env)])
  }
  if (is.data.frame(value)) {
    env_name <- if (env_col %in% names(value)) {
      env_col
    } else if ("Env" %in% names(value)) {
      "Env"
    } else {
      names(value)[[1L]]
    }
    candidate_names <- setdiff(names(value), env_name)
    if (!is.null(value_name) && value_name %in% names(value)) {
      val_name <- value_name
    } else if (length(candidate_names)) {
      val_name <- candidate_names[[1L]]
    } else {
      stop("Environment annotation data frame must contain a value column.", call. = FALSE)
    }
    ann <- data.frame(
      Env = as.character(value[[env_name]]),
      value = value[[val_name]],
      stringsAsFactors = FALSE
    )
    ann <- ann[!is.na(ann$value) & !duplicated(ann$Env), , drop = FALSE]
    return(ann$value[match(env_levels, ann$Env)])
  }
  if (!is.null(names(value)) && any(nzchar(names(value)))) {
    return(unname(value[match(env_levels, names(value))]))
  }
  if (length(value) == length(env_levels)) {
    return(value)
  }
  stop(
    "Environment annotation must be a named vector, an env-level data frame, ",
    "a phenotype column name, or a vector aligned to the environment levels.",
    call. = FALSE
  )
}

gp_tuning_resolve_pheno_annotation <- function(value,
                                               pheno_data,
                                               pheno_env,
                                               target_name) {
  if (is.null(value) || is.null(pheno_data) || is.null(pheno_env)) {
    return(value)
  }
  if (!is.character(value) || length(value) != 1L || !value %in% names(pheno_data)) {
    return(value)
  }
  ann <- data.frame(
    Env = as.character(pheno_env),
    value = pheno_data[[value]],
    stringsAsFactors = FALSE
  )
  ann <- ann[!is.na(ann$value) & !duplicated(ann$Env), , drop = FALSE]
  names(ann)[names(ann) == "value"] <- target_name
  ann
}

gp_tuning_env_metadata <- function(env_levels,
                                   pheno_df = NULL,
                                   env_col = "Env",
                                   env_location = NULL,
                                   env_year = NULL) {
  env_levels <- unique(as.character(env_levels))
  location <- gp_tuning_env_annotation(
    env_location,
    env_levels,
    pheno_df = pheno_df,
    env_col = env_col,
    value_name = "location"
  )
  if (is.null(location)) {
    location <- gp_tuning_infer_location(env_levels)
  }
  location <- as.character(location)
  missing_location <- is.na(location) | !nzchar(location)
  location[missing_location] <- env_levels[missing_location]

  year <- gp_tuning_env_annotation(
    env_year,
    env_levels,
    pheno_df = pheno_df,
    env_col = env_col,
    value_name = "year"
  )
  if (is.null(year)) {
    year <- gp_tuning_parse_year(env_levels)
  }
  year <- suppressWarnings(as.integer(year))
  data.frame(
    Env = env_levels,
    location = location,
    year = year,
    env_order = seq_along(env_levels),
    stringsAsFactors = FALSE
  )
}

gp_tuning_latest_env <- function(meta) {
  ord <- order(
    ifelse(is.na(meta$year), -Inf, meta$year),
    meta$env_order,
    decreasing = TRUE
  )
  meta$Env[ord[[1L]]]
}

gp_tuning_blocked_env_split <- function(train_envs,
                                        pheno_df = NULL,
                                        env_col = "Env",
                                        env_location = NULL,
                                        env_year = NULL,
                                        validation_fraction = 0.25) {
  train_envs <- unique(as.character(train_envs))
  meta <- gp_tuning_env_metadata(
    train_envs,
    pheno_df = pheno_df,
    env_col = env_col,
    env_location = env_location,
    env_year = env_year
  )
  by_location <- split(meta, meta$location)
  validation_envs <- unname(vapply(
    by_location,
    function(x) if (nrow(x) >= 2L) gp_tuning_latest_env(x) else NA_character_,
    character(1L)
  ))
  validation_envs <- unique(validation_envs[!is.na(validation_envs)])
  calibration_envs <- setdiff(train_envs, validation_envs)

  if (!length(validation_envs) || !length(calibration_envs)) {
    n_hold <- max(1L, floor(length(train_envs) * validation_fraction))
    n_hold <- min(n_hold, length(train_envs) - 1L)
    ord <- order(
      ifelse(is.na(meta$year), -Inf, meta$year),
      meta$env_order,
      decreasing = TRUE
    )
    validation_envs <- meta$Env[head(ord, n_hold)]
    calibration_envs <- setdiff(train_envs, validation_envs)
  }

  list(
    validation_envs = validation_envs,
    calibration_envs = calibration_envs,
    env_metadata = meta
  )
}

gp_tuning_make_blocked_split <- function(pheno_df,
                                         train_mask,
                                         env_location = NULL,
                                         env_year = NULL,
                                         validation_fraction = 0.25,
                                         seed = 12345L) {
  if (!all(c("GID", "Env", "y") %in% names(pheno_df))) {
    stop("`pheno_df` must contain GID, Env, and y columns.", call. = FALSE)
  }
  validation_fraction <- as.numeric(validation_fraction)[1L]
  if (!is.finite(validation_fraction) || validation_fraction <= 0 || validation_fraction >= 1) {
    validation_fraction <- 0.25
  }
  finite_train <- train_mask & is.finite(as.numeric(pheno_df$y))
  if (sum(finite_train) < 4L) {
    return(list(status = "skipped", reason = "fewer than four finite training rows"))
  }
  train_envs <- unique(as.character(pheno_df$Env[finite_train]))
  train_ids <- unique(as.character(pheno_df$GID[finite_train]))
  if (length(train_envs) < 2L || length(train_ids) < 2L) {
    return(list(status = "skipped", reason = "blocked tuning needs at least two environments and two genotypes"))
  }

  env_split <- gp_tuning_blocked_env_split(
    train_envs,
    pheno_df = pheno_df[finite_train, , drop = FALSE],
    env_location = env_location,
    env_year = env_year,
    validation_fraction = validation_fraction
  )
  if (!length(env_split$validation_envs) || !length(env_split$calibration_envs)) {
    return(list(status = "skipped", reason = "could not form calibration and validation environment blocks"))
  }

  calibration_pool <- finite_train & pheno_df$Env %in% env_split$calibration_envs
  validation_pool <- finite_train & pheno_df$Env %in% env_split$validation_envs
  eligible_validation_ids <- unique(as.character(pheno_df$GID[validation_pool]))
  if (!length(eligible_validation_ids)) {
    return(list(status = "skipped", reason = "no finite rows in validation environments"))
  }

  gp_set_seed(as.integer(seed)[1L])
  eligible_validation_ids <- sample(eligible_validation_ids)
  target_n <- max(1L, ceiling(length(train_ids) * validation_fraction))
  target_n <- min(target_n, length(eligible_validation_ids), length(train_ids) - 1L)
  if (target_n < 1L) {
    return(list(status = "skipped", reason = "no genotype can be held out without emptying calibration"))
  }

  for (n_val in rev(seq_len(target_n))) {
    validation_ids <- head(eligible_validation_ids, n_val)
    calibration_rows <- which(calibration_pool & !pheno_df$GID %in% validation_ids)
    validation_rows <- which(validation_pool & pheno_df$GID %in% validation_ids)
    if (length(calibration_rows) >= 3L && length(validation_rows) >= 1L) {
      return(list(
        status = "ok",
        calibration_rows = calibration_rows,
        validation_rows = validation_rows,
        calibration_envs = env_split$calibration_envs,
        validation_envs = env_split$validation_envs,
        calibration_ids = unique(as.character(pheno_df$GID[calibration_rows])),
        validation_ids = unique(as.character(pheno_df$GID[validation_rows])),
        env_metadata = env_split$env_metadata
      ))
    }
  }

  list(status = "skipped", reason = "blocked split produced too few calibration or validation rows")
}

gp_tuning_make_genotype_split <- function(pheno_df,
                                          train_mask,
                                          env_location = NULL,
                                          env_year = NULL,
                                          validation_fraction = 0.25,
                                          seed = 12345L,
                                          allow_single_env = FALSE) {
  if (!all(c("GID", "Env", "y") %in% names(pheno_df))) {
    stop("`pheno_df` must contain GID, Env, and y columns.", call. = FALSE)
  }
  validation_fraction <- as.numeric(validation_fraction)[1L]
  if (!is.finite(validation_fraction) || validation_fraction <= 0 || validation_fraction >= 1) {
    validation_fraction <- 0.25
  }
  finite_train <- train_mask & is.finite(as.numeric(pheno_df$y))
  if (sum(finite_train) < 4L) {
    return(list(status = "skipped", reason = "fewer than four finite training rows"))
  }
  train_envs <- unique(as.character(pheno_df$Env[finite_train]))
  train_ids <- unique(as.character(pheno_df$GID[finite_train]))
  min_env <- if (isTRUE(allow_single_env)) 1L else 2L
  if (length(train_envs) < min_env || length(train_ids) < 3L) {
    env_text <- if (isTRUE(allow_single_env)) "one environment" else "two environments"
    return(list(
      status = "skipped",
      reason = paste("genotype-blocked tuning needs at least", env_text, "and three genotypes")
    ))
  }

  gp_set_seed(as.integer(seed)[1L])
  train_ids <- sample(train_ids)
  target_n <- max(1L, ceiling(length(train_ids) * validation_fraction))
  target_n <- min(target_n, length(train_ids) - 1L)
  if (target_n < 1L) {
    return(list(status = "skipped", reason = "no genotype can be held out without emptying calibration"))
  }

  for (n_val in rev(seq_len(target_n))) {
    validation_ids <- head(train_ids, n_val)
    calibration_rows <- which(finite_train & !pheno_df$GID %in% validation_ids)
    validation_rows <- which(finite_train & pheno_df$GID %in% validation_ids)
    calibration_envs <- unique(as.character(pheno_df$Env[calibration_rows]))
    validation_envs <- unique(as.character(pheno_df$Env[validation_rows]))
    if (length(calibration_rows) >= 3L && length(validation_rows) >= 1L &&
        length(calibration_envs) >= min_env && length(validation_envs) >= 1L) {
      return(list(
        status = "ok",
        calibration_rows = calibration_rows,
        validation_rows = validation_rows,
        calibration_envs = calibration_envs,
        validation_envs = validation_envs,
        calibration_ids = unique(as.character(pheno_df$GID[calibration_rows])),
        validation_ids = unique(as.character(pheno_df$GID[validation_rows])),
        env_metadata = gp_tuning_env_metadata(
          unique(c(calibration_envs, validation_envs)),
          pheno_df = pheno_df[finite_train, , drop = FALSE],
          env_location = env_location,
          env_year = env_year
        )
      ))
    }
  }

  list(status = "skipped", reason = "genotype-blocked split produced too few calibration or validation rows")
}

gp_tuning_max_or_inf <- function(x) {
  val <- suppressWarnings(as.numeric(x)[1L])
  if (!is.finite(val) || is.na(val) || val <= 0) {
    return(Inf)
  }
  floor(val)
}

gp_tuning_sample_rows_by_env <- function(rows, pheno_df, max_rows, seed = 12345L) {
  rows <- as.integer(rows)
  max_rows <- gp_tuning_max_or_inf(max_rows)
  if (!length(rows) || length(rows) <= max_rows) {
    return(rows)
  }
  by_env <- split(rows, as.character(pheno_df$Env[rows]))
  gp_set_seed(as.integer(seed)[1L])
  if (length(by_env) > max_rows) {
    env_keep <- sample(names(by_env), as.integer(max_rows))
    by_env <- by_env[env_keep]
  }
  per_env <- max(1L, floor(max_rows / max(1L, length(by_env))))
  picked <- unlist(
    lapply(by_env, function(x) sample(x, min(length(x), per_env))),
    use.names = FALSE
  )
  if (length(picked) < max_rows) {
    remaining <- setdiff(rows, picked)
    need <- min(length(remaining), as.integer(max_rows) - length(picked))
    if (need > 0L) {
      picked <- c(picked, sample(remaining, need))
    }
  }
  sort(as.integer(picked))
}

gp_tuning_limit_genotypes <- function(calibration_rows,
                                      validation_rows,
                                      pheno_df,
                                      max_genotypes,
                                      seed = 12345L) {
  max_genotypes <- gp_tuning_max_or_inf(max_genotypes)
  rows <- c(calibration_rows, validation_rows)
  ids <- unique(as.character(pheno_df$GID[rows]))
  if (!length(ids) || length(ids) <= max_genotypes) {
    return(list(
      calibration_rows = calibration_rows,
      validation_rows = validation_rows,
      limited = FALSE,
      n_ids_before = length(ids),
      n_ids_after = length(ids)
    ))
  }
  gp_set_seed(as.integer(seed)[1L] + 11L)
  validation_ids <- unique(as.character(pheno_df$GID[validation_rows]))
  calibration_ids <- unique(as.character(pheno_df$GID[calibration_rows]))
  n_validation_target <- min(
    length(validation_ids),
    max(1L, floor(as.integer(max_genotypes) * 0.35))
  )
  keep_validation <- if (n_validation_target > 0L) {
    sample(validation_ids, n_validation_target)
  } else {
    character()
  }
  remaining_n <- as.integer(max_genotypes) - length(keep_validation)
  calibration_pool <- setdiff(calibration_ids, keep_validation)
  keep_calibration <- if (remaining_n > 0L && length(calibration_pool)) {
    sample(calibration_pool, min(length(calibration_pool), remaining_n))
  } else {
    character()
  }
  keep_ids <- unique(c(keep_validation, keep_calibration))
  if (length(keep_ids) < max_genotypes) {
    rest <- setdiff(ids, keep_ids)
    need <- min(length(rest), as.integer(max_genotypes) - length(keep_ids))
    if (need > 0L) {
      keep_ids <- c(keep_ids, sample(rest, need))
    }
  }
  calibration_rows2 <- calibration_rows[as.character(pheno_df$GID[calibration_rows]) %in% keep_ids]
  validation_rows2 <- validation_rows[as.character(pheno_df$GID[validation_rows]) %in% keep_ids]
  list(
    calibration_rows = as.integer(calibration_rows2),
    validation_rows = as.integer(validation_rows2),
    limited = TRUE,
    n_ids_before = length(ids),
    n_ids_after = length(unique(as.character(pheno_df$GID[c(calibration_rows2, validation_rows2)])))
  )
}

gp_tuning_subset_blocked_rows <- function(pheno_df,
                                          calibration_rows,
                                          validation_rows,
                                          max_calibration_rows = NULL,
                                          max_validation_rows = NULL,
                                          max_genotypes = NULL,
                                          seed = 12345L) {
  n_cal_before <- length(calibration_rows)
  n_val_before <- length(validation_rows)
  ids_before <- length(unique(as.character(pheno_df$GID[c(calibration_rows, validation_rows)])))
  calibration_rows <- gp_tuning_sample_rows_by_env(
    calibration_rows,
    pheno_df = pheno_df,
    max_rows = max_calibration_rows,
    seed = seed
  )
  validation_rows <- gp_tuning_sample_rows_by_env(
    validation_rows,
    pheno_df = pheno_df,
    max_rows = max_validation_rows,
    seed = as.integer(seed)[1L] + 1L
  )
  limited_ids <- gp_tuning_limit_genotypes(
    calibration_rows,
    validation_rows,
    pheno_df = pheno_df,
    max_genotypes = max_genotypes,
    seed = seed
  )
  calibration_rows <- limited_ids$calibration_rows
  validation_rows <- limited_ids$validation_rows
  list(
    status = if (length(calibration_rows) >= 3L && length(validation_rows) >= 1L) "ok" else "skipped",
    reason = if (length(calibration_rows) >= 3L && length(validation_rows) >= 1L) {
      NA_character_
    } else {
      "bounded tuning subset produced too few calibration or validation rows"
    },
    calibration_rows = calibration_rows,
    validation_rows = validation_rows,
    sample_info = list(
      n_calibration_before = as.integer(n_cal_before),
      n_validation_before = as.integer(n_val_before),
      n_genotypes_before = as.integer(ids_before),
      n_calibration_after = as.integer(length(calibration_rows)),
      n_validation_after = as.integer(length(validation_rows)),
      n_genotypes_after = as.integer(length(unique(as.character(pheno_df$GID[c(calibration_rows, validation_rows)])))),
      max_calibration_rows = gp_tuning_max_or_inf(max_calibration_rows),
      max_validation_rows = gp_tuning_max_or_inf(max_validation_rows),
      max_genotypes = gp_tuning_max_or_inf(max_genotypes),
      genotype_limited = isTRUE(limited_ids$limited)
    )
  )
}

gp_tuning_dense_kernel_for_ids <- function(kernel,
                                           ids,
                                           gp_factor_cache = NULL,
                                           geno_ids = NULL,
                                           python_bin = NULL,
                                           project_root = NULL) {
  ids <- unique(as.character(ids))
  if (!length(ids)) {
    stop("No genotype IDs available for GP tuning.", call. = FALSE)
  }
  if (!is.null(kernel$K) && all(ids %in% rownames(kernel$K))) {
    return(kernel$K[ids, ids, drop = FALSE])
  }
  if (!is.null(gp_factor_cache)) {
    return(gp_public_factor_cache_kernel(
      gp_factor_cache = gp_factor_cache,
      ids = ids,
      geno_ids = geno_ids %||% kernel$geno_ids,
      python_bin = python_bin,
      project_root = project_root,
      dtype = "float64"
    ))
  }
  stop(
    "GP tuning needs either a dense `gmatrix` or a `gp_factor_cache` with `geno_ids`.",
    call. = FALSE
  )
}

gp_tuning_location_offsets <- function(pheno_df,
                                       training_rows,
                                       mode = "none",
                                       env_location = NULL,
                                       env_year = NULL) {
  mode <- tolower(as.character(mode)[1L])
  if (!identical(mode, "location_offset")) {
    return(rep(0, nrow(pheno_df)))
  }
  y <- as.numeric(pheno_df$y)
  training_rows <- training_rows[is.finite(y[training_rows])]
  if (!length(training_rows)) {
    return(rep(0, nrow(pheno_df)))
  }
  meta <- gp_tuning_env_metadata(
    unique(as.character(pheno_df$Env)),
    pheno_df = pheno_df,
    env_location = env_location,
    env_year = env_year
  )
  loc_by_env <- stats::setNames(meta$location, meta$Env)
  train_locs <- loc_by_env[as.character(pheno_df$Env[training_rows])]
  loc_means <- tapply(y[training_rows], train_locs, mean, na.rm = TRUE)
  global_mean <- mean(y[training_rows], na.rm = TRUE)
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  offsets <- as.numeric(loc_means[loc_by_env[as.character(pheno_df$Env)]])
  offsets[!is.finite(offsets)] <- global_mean
  offsets
}

gp_tuning_apply_mean_adjustment <- function(pheno_df,
                                            training_rows,
                                            mode = "none",
                                            env_location = NULL,
                                            env_year = NULL) {
  offsets <- gp_tuning_location_offsets(
    pheno_df,
    training_rows = training_rows,
    mode = mode,
    env_location = env_location,
    env_year = env_year
  )
  out <- pheno_df
  finite_y <- is.finite(as.numeric(out$y))
  out$y[finite_y] <- as.numeric(out$y[finite_y]) - offsets[finite_y]
  list(pheno_df = out, offsets = offsets, mode = tolower(as.character(mode)[1L]))
}

gp_tuning_envmean_baseline <- function(training_df,
                                       target_df,
                                       env_location = NULL,
                                       env_year = NULL) {
  if (!all(c("Env", "y") %in% names(training_df)) || !"Env" %in% names(target_df)) {
    return(rep(NA_real_, nrow(target_df)))
  }
  training_df <- as.data.frame(training_df, stringsAsFactors = FALSE)
  target_df <- as.data.frame(target_df, stringsAsFactors = FALSE)
  training_df$Env <- as.character(training_df$Env)
  target_df$Env <- as.character(target_df$Env)
  y <- as.numeric(training_df$y)
  finite <- is.finite(y)
  if (!any(finite)) {
    return(rep(NA_real_, nrow(target_df)))
  }
  global_mean <- mean(y[finite], na.rm = TRUE)
  env_mean <- stats::setNames(
    as.numeric(tapply(y[finite], training_df$Env[finite], mean, na.rm = TRUE)),
    names(tapply(y[finite], training_df$Env[finite], mean, na.rm = TRUE))
  )
  baseline <- as.numeric(env_mean[target_df$Env])
  missing <- !is.finite(baseline)
  if (any(missing)) {
    env_levels <- unique(c(training_df$Env, target_df$Env))
    meta_df <- training_df
    target_meta <- target_df
    all_cols <- union(names(meta_df), names(target_meta))
    for (nm in setdiff(all_cols, names(meta_df))) {
      meta_df[[nm]] <- NA
    }
    for (nm in setdiff(all_cols, names(target_meta))) {
      target_meta[[nm]] <- NA
    }
    meta_df <- rbind(meta_df[all_cols], target_meta[all_cols])
    meta <- gp_tuning_env_metadata(
      env_levels,
      pheno_df = meta_df,
      env_location = env_location,
      env_year = env_year
    )
    loc_map <- stats::setNames(as.character(meta$location), as.character(meta$Env))
    train_loc <- loc_map[training_df$Env]
    target_loc <- loc_map[target_df$Env]
    loc_mean <- stats::setNames(
      as.numeric(tapply(y[finite], train_loc[finite], mean, na.rm = TRUE)),
      names(tapply(y[finite], train_loc[finite], mean, na.rm = TRUE))
    )
    baseline[missing] <- as.numeric(loc_mean[target_loc[missing]])
  }
  baseline[!is.finite(baseline)] <- global_mean
  as.numeric(baseline)
}

gp_plain_data_frame <- function(x) {
  py_vector <- function(col) {
    if (!inherits(col, "python.builtin.object")) {
      return(col)
    }
    out <- tryCatch({
      if (requireNamespace("reticulate", quietly = TRUE) &&
          reticulate::py_has_attr(col, "to_numpy")) {
        reticulate::py_to_r(col$to_numpy())
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
    out <- tryCatch({
      if (requireNamespace("reticulate", quietly = TRUE) &&
          reticulate::py_has_attr(col, "tolist")) {
        reticulate::py_to_r(col$tolist())
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
    tryCatch(reticulate::py_to_r(col), error = function(e) col)
  }

  if (inherits(x, "python.builtin.object")) {
    x <- py_vector(x)
  }
  if (is.data.frame(x) || is.list(x)) {
    x_names <- names(x)
    out <- lapply(seq_along(x), function(i) {
      col <- py_vector(x[[i]])
      if (inherits(col, "python.builtin.object")) {
        col <- as.character(col)
      }
      if (is.data.frame(col) && ncol(col) == 1L) {
        col <- col[[1L]]
      }
      if (is.matrix(col) || (is.array(col) && length(dim(col)) <= 2L)) {
        col <- as.vector(col)
      }
      if (is.list(col) && !is.data.frame(col)) {
        flat <- tryCatch(unlist(col, recursive = FALSE, use.names = FALSE), error = function(e) col)
        if (!is.list(flat) || length(flat) == length(col)) {
          col <- flat
        }
      }
      col
    })
    out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    names(out) <- x_names
    return(out)
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  out <- lapply(x, function(col) {
    col <- py_vector(col)
    if (is.data.frame(col) && ncol(col) == 1L) {
      col <- col[[1L]]
    }
    if (is.matrix(col) || (is.array(col) && length(dim(col)) <= 2L)) {
      col <- as.vector(col)
    }
    if (is.list(col) && !is.data.frame(col)) {
      flat <- tryCatch(unlist(col, recursive = FALSE, use.names = FALSE), error = function(e) col)
      if (!is.list(flat) || length(flat) == length(col)) {
        col <- flat
      }
    }
    col
  })
  out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- names(x)
  out
}

gp_single_trait_prediction_keys <- function(pred) {
  pred <- gp_plain_data_frame(pred)
  if (!"GID" %in% names(pred) && "gid" %in% names(pred)) {
    names(pred)[names(pred) == "gid"] <- "GID"
  }
  if (!"GID" %in% names(pred) && "Hybrid" %in% names(pred)) {
    names(pred)[names(pred) == "Hybrid"] <- "GID"
  }
  if (!"GID" %in% names(pred) && "Genotype" %in% names(pred)) {
    names(pred)[names(pred) == "Genotype"] <- "GID"
  }
  if (!"GID" %in% names(pred) && "Name" %in% names(pred)) {
    names(pred)[names(pred) == "Name"] <- "GID"
  }
  if (!"Env" %in% names(pred) && "env" %in% names(pred)) {
    names(pred)[names(pred) == "env"] <- "Env"
  }
  pred
}

gp_apply_prediction_offset <- function(result, pheno_df, split, offsets, label = "location_offset") {
  if (is.null(result$predictions) || !"Prediction" %in% names(result$predictions)) {
    return(result)
  }
  pred <- as.data.frame(result$predictions, stringsAsFactors = FALSE)
  original_names <- names(pred)
  normalized <- gp_single_trait_prediction_keys(pred)
  applied <- rep(NA_real_, nrow(normalized))

  if (all(c("GID", "Env") %in% names(normalized))) {
    offset_df <- data.frame(
      GID = as.character(pheno_df$GID),
      Env = as.character(pheno_df$Env),
      offset = as.numeric(offsets),
      stringsAsFactors = FALSE
    )
    offset_df <- stats::aggregate(offset ~ GID + Env, offset_df, mean)
    key <- paste(normalized$GID, normalized$Env, sep = "\r")
    offset_key <- paste(offset_df$GID, offset_df$Env, sep = "\r")
    applied <- offset_df$offset[match(key, offset_key)]
  }

  if ((!length(applied) || all(!is.finite(applied))) &&
      nrow(pred) == sum(split$test_mask)) {
    applied <- as.numeric(offsets[split$test_mask])
  }

  applied[!is.finite(applied)] <- 0
  pred$Prediction <- as.numeric(pred$Prediction) + applied
  pred$Prediction_Mean_Offset <- applied
  pred$Prediction_Mean_Adjustment <- label
  names(pred)[seq_along(original_names)] <- original_names
  result$predictions <- pred
  if (!is.null(result$result) && !is.null(result$result$predictions)) {
    result$result$predictions <- pred
  }
  result
}

gp_apply_prediction_blend <- function(result,
                                      pheno_df,
                                      split,
                                      alpha,
                                      env_location = NULL,
                                      env_year = NULL,
                                      label = "envmean_validation_blend") {
  if (is.null(result$predictions) || !"Prediction" %in% names(result$predictions)) {
    return(list(result = result, info = list(applied = FALSE, reason = "prediction column unavailable")))
  }
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  if (!is.finite(alpha)) {
    return(list(result = result, info = list(applied = FALSE, reason = "blend alpha unavailable")))
  }
  alpha <- pmin(1, pmax(0, alpha))
  target_rows <- which(split$test_mask)
  if (!length(target_rows)) {
    return(list(result = result, info = list(applied = FALSE, reason = "no prediction rows")))
  }
  baseline <- gp_tuning_envmean_baseline(
    training_df = pheno_df[split$train_mask, c("GID", "Env", "y"), drop = FALSE],
    target_df = pheno_df[target_rows, c("GID", "Env", "y"), drop = FALSE],
    env_location = env_location,
    env_year = env_year
  )
  baseline_df <- data.frame(
    GID = as.character(pheno_df$GID[target_rows]),
    Env = as.character(pheno_df$Env[target_rows]),
    Prediction_EnvMean = as.numeric(baseline),
    stringsAsFactors = FALSE
  )
  baseline_df <- stats::aggregate(Prediction_EnvMean ~ GID + Env, baseline_df, mean)
  pred <- as.data.frame(result$predictions, stringsAsFactors = FALSE)
  original_names <- names(pred)
  normalized <- gp_single_trait_prediction_keys(pred)
  if (!all(c("GID", "Env") %in% names(normalized))) {
    if (nrow(normalized) == length(target_rows)) {
      normalized$GID <- as.character(pheno_df$GID[target_rows])
      normalized$Env <- as.character(pheno_df$Env[target_rows])
    } else {
      return(list(result = result, info = list(applied = FALSE, reason = "prediction keys unavailable")))
    }
  }
  key <- paste(as.character(normalized$GID), as.character(normalized$Env), sep = "\r")
  baseline_key <- paste(baseline_df$GID, baseline_df$Env, sep = "\r")
  baseline_match <- baseline_df$Prediction_EnvMean[match(key, baseline_key)]
  apply_rows <- is.finite(baseline_match) & is.finite(as.numeric(normalized$Prediction))
  if (!any(apply_rows)) {
    return(list(result = result, info = list(applied = FALSE, reason = "no predictions matched blend baseline")))
  }
  gp_pred <- as.numeric(normalized$Prediction)
  blended <- gp_pred
  blended[apply_rows] <- alpha * gp_pred[apply_rows] + (1 - alpha) * baseline_match[apply_rows]
  normalized$Prediction_GP <- gp_pred
  normalized$Prediction_EnvMean <- baseline_match
  normalized$Prediction <- blended
  normalized$Prediction_Blend_Alpha <- ifelse(apply_rows, alpha, NA_real_)
  normalized$Prediction_Blend_Method <- ifelse(apply_rows, label, NA_character_)
  var_cols <- intersect(c("Prediction_Var_observed", "Prediction_Var_latent", "Prediction_variance", "PEV"), names(normalized))
  for (nm in var_cols) {
    vals <- suppressWarnings(as.numeric(normalized[[nm]]))
    vals[apply_rows] <- vals[apply_rows] * alpha^2
    normalized[[nm]] <- vals
  }
  se_cols <- intersect(c("Prediction_SE_observed", "Prediction_SE_latent", "Prediction_SE", "SE"), names(normalized))
  for (nm in se_cols) {
    vals <- suppressWarnings(as.numeric(normalized[[nm]]))
    vals[apply_rows] <- vals[apply_rows] * abs(alpha)
    normalized[[nm]] <- vals
  }
  names(normalized)[seq_along(original_names)] <- original_names
  result$predictions <- normalized
  if (!is.null(result$result) && !is.null(result$result$predictions)) {
    result$result$predictions <- normalized
  }
  list(
    result = result,
    info = list(
      applied = TRUE,
      reason = NA_character_,
      alpha = alpha,
      n_blended = sum(apply_rows),
      baseline = "env_or_location_mean",
      se_variance_note = "Prediction SE/variance scaled by alpha and treats the baseline as fixed."
    )
  )
}

gp_apply_uncertainty_calibration <- function(result,
                                             calibration,
                                             label = "blocked_validation") {
  scale <- suppressWarnings(as.numeric(calibration$scale %||% NA_real_)[1L])
  if (!identical(calibration$status %||% NA_character_, "ok") ||
      !is.finite(scale) || scale <= 0) {
    return(list(
      result = result,
      info = c(
        list(requested = TRUE, applied = FALSE),
        calibration,
        list(reason = calibration$reason %||% "uncertainty calibration unavailable")
      )
    ))
  }
  if (is.null(result$predictions) || !nrow(gp_plain_data_frame(result$predictions))) {
    return(list(
      result = result,
      info = c(
        list(requested = TRUE, applied = FALSE),
        calibration,
        list(reason = "prediction table unavailable")
      )
    ))
  }
  pred <- gp_plain_data_frame(result$predictions)
  original_names <- names(pred)
  normalized <- gp_single_trait_prediction_keys(pred)
  var_cols <- intersect(
    c("Prediction_Var_observed", "Prediction_Var_latent", "Prediction_variance", "PEV"),
    names(normalized)
  )
  se_cols <- intersect(
    c("Prediction_SE_observed", "Prediction_SE_latent", "Prediction_SE", "SE", "SE_observed", "SE_latent"),
    names(normalized)
  )
  if (!length(var_cols) && !length(se_cols)) {
    return(list(
      result = result,
      info = c(
        list(requested = TRUE, applied = FALSE),
        calibration,
        list(reason = "prediction SE/variance columns unavailable")
      )
    ))
  }
  for (nm in var_cols) {
    vals <- suppressWarnings(as.numeric(normalized[[nm]]))
    normalized[[nm]] <- vals * scale^2
  }
  for (nm in se_cols) {
    vals <- suppressWarnings(as.numeric(normalized[[nm]]))
    normalized[[nm]] <- vals * scale
  }
  normalized$Prediction_SE_Calibration_Scale <- scale
  normalized$Prediction_SE_Calibration_Method <- label
  names(normalized)[seq_along(original_names)] <- original_names
  result$predictions <- normalized
  if (!is.null(result$result) && !is.null(result$result$predictions)) {
    result$result$predictions <- normalized
  }
  list(
    result = result,
    info = c(
      list(requested = TRUE, applied = TRUE),
      calibration,
      list(
        scale = scale,
        variance_scale = scale^2,
        label = label,
        scaled_se_columns = paste(se_cols, collapse = ","),
        scaled_variance_columns = paste(var_cols, collapse = ",")
      )
    )
  )
}

gp_tuning_score_predictions <- function(fit,
                                        truth,
                                        offsets,
                                        baseline = NULL,
                                        blend_objective = "rmse") {
  pred <- gp_single_trait_prediction_keys(fit$predictions)
  if (!all(c("GID", "Env", "Prediction") %in% names(pred))) {
    stop("Tuning predictions must contain genotype, environment, and Prediction columns.", call. = FALSE)
  }
  offset_df <- data.frame(
    GID = truth$GID,
    Env = truth$Env,
    Prediction_Offset = as.numeric(offsets),
    stringsAsFactors = FALSE
  )
  joined <- merge(pred, truth, by = c("GID", "Env"), sort = FALSE)
  joined <- merge(joined, offset_df, by = c("GID", "Env"), sort = FALSE)
  if (!nrow(joined)) {
    stop("Tuning candidate produced no validation predictions.", call. = FALSE)
  }
  joined$Prediction <- as.numeric(joined$Prediction) + as.numeric(joined$Prediction_Offset)
  if (!is.null(baseline)) {
    baseline_df <- data.frame(
      GID = truth$GID,
      Env = truth$Env,
      Prediction_Baseline = as.numeric(baseline),
      stringsAsFactors = FALSE
    )
    joined <- merge(joined, baseline_df, by = c("GID", "Env"), sort = FALSE)
  } else {
    joined$Prediction_Baseline <- NA_real_
  }
  keep <- is.finite(joined$Prediction) & is.finite(joined$Observed)
  joined <- joined[keep, , drop = FALSE]
  if (!nrow(joined)) {
    stop("Tuning candidate produced no finite validation predictions.", call. = FALSE)
  }
  by_env <- split(joined, joined$Env)
  env_rmse <- vapply(
    by_env,
    function(x) sqrt(mean((as.numeric(x$Prediction) - as.numeric(x$Observed))^2)),
    numeric(1L)
  )
  env_n <- vapply(by_env, nrow, integer(1L))
  pred_centered <- as.numeric(joined$Prediction) -
    ave(as.numeric(joined$Prediction), joined$Env, FUN = function(x) mean(x, na.rm = TRUE))
  obs_centered <- as.numeric(joined$Observed) -
    ave(as.numeric(joined$Observed), joined$Env, FUN = function(x) mean(x, na.rm = TRUE))
  centered_eligible <- env_n[as.character(joined$Env)] >= 2L
  centered_rmse <- if (any(centered_eligible)) {
    sqrt(mean((pred_centered[centered_eligible] - obs_centered[centered_eligible])^2))
  } else {
    NA_real_
  }
  centered_mae <- if (any(centered_eligible)) {
    mean(abs(pred_centered[centered_eligible] - obs_centered[centered_eligible]))
  } else {
    NA_real_
  }
  env_centered_rmse <- vapply(
    by_env,
    function(x) {
      if (nrow(x) < 2L) {
        return(NA_real_)
      }
      xp <- as.numeric(x$Prediction) - mean(as.numeric(x$Prediction), na.rm = TRUE)
      xo <- as.numeric(x$Observed) - mean(as.numeric(x$Observed), na.rm = TRUE)
      sqrt(mean((xp - xo)^2))
    },
    numeric(1L)
  )
  env_spearman <- vapply(
    by_env,
    function(x) {
      if (nrow(x) < 2L) {
        return(NA_real_)
      }
      gp_tuning_safe_cor(x$Prediction, x$Observed, method = "spearman")
    },
    numeric(1L)
  )
  env_centered_ok <- is.finite(env_centered_rmse)
  env_spearman_ok <- is.finite(env_spearman)
  blend <- gp_tuning_select_blend(
    gp_prediction = joined$Prediction,
    baseline_prediction = joined$Prediction_Baseline,
    observed = joined$Observed,
    env = joined$Env,
    objective = blend_objective
  )
  list(
    rmse = sqrt(mean((as.numeric(joined$Prediction) - as.numeric(joined$Observed))^2)),
    mae = mean(abs(as.numeric(joined$Prediction) - as.numeric(joined$Observed))),
    cor = gp_tuning_safe_cor(joined$Prediction, joined$Observed),
    centered_rmse = centered_rmse,
    centered_mae = centered_mae,
    centered_cor = if (any(centered_eligible)) {
      gp_tuning_safe_cor(pred_centered[centered_eligible], obs_centered[centered_eligible])
    } else {
      NA_real_
    },
    within_env_spearman = if (any(env_spearman_ok)) {
      stats::weighted.mean(env_spearman[env_spearman_ok], env_n[env_spearman_ok])
    } else {
      NA_real_
    },
    rmse_se = if (length(env_rmse) > 1L) stats::sd(env_rmse) / sqrt(length(env_rmse)) else NA_real_,
    centered_rmse_se = if (sum(env_centered_ok) > 1L) {
      stats::sd(env_centered_rmse[env_centered_ok]) / sqrt(sum(env_centered_ok))
    } else {
      NA_real_
    },
    within_env_spearman_se = if (sum(env_spearman_ok) > 1L) {
      stats::sd(env_spearman[env_spearman_ok]) / sqrt(sum(env_spearman_ok))
    } else {
      NA_real_
    },
    blend_alpha = as.numeric(blend$blend_alpha),
    blend_objective = as.character(blend$blend_objective),
    blend_objective_value = as.numeric(blend$blend_objective_value),
    blend_rmse = as.numeric(blend$blend_rmse),
    blend_cor = as.numeric(blend$blend_cor),
    blend_centered_rmse = as.numeric(blend$blend_centered_rmse),
    blend_centered_cor = as.numeric(blend$blend_centered_cor),
    blend_within_env_spearman = as.numeric(blend$blend_within_env_spearman),
    baseline_rmse = as.numeric(blend$baseline_rmse),
    baseline_cor = as.numeric(blend$baseline_cor),
    baseline_centered_rmse = as.numeric(blend$baseline_centered_rmse),
    baseline_centered_cor = as.numeric(blend$baseline_centered_cor),
    baseline_within_env_spearman = as.numeric(blend$baseline_within_env_spearman),
    n_validation = nrow(joined),
    env_detail = data.frame(
      env = names(env_rmse),
      n_validation = as.integer(env_n),
      rmse = as.numeric(env_rmse),
      centered_rmse = as.numeric(env_centered_rmse),
      spearman = as.numeric(env_spearman),
      stringsAsFactors = FALSE
    )
  )
}

gp_tuning_pick_numeric_column_info <- function(x, choices) {
  hits <- intersect(as.character(choices), names(x))
  if (!length(hits)) {
    return(list(values = rep(NA_real_, nrow(x)), column = NA_character_))
  }
  fallback <- rep(NA_real_, nrow(x))
  fallback_column <- NA_character_
  for (hit in hits) {
    vals <- suppressWarnings(as.numeric(x[[hit]]))
    if (any(is.finite(vals))) {
      return(list(values = vals, column = hit))
    }
    fallback <- vals
    fallback_column <- hit
  }
  list(values = fallback, column = fallback_column)
}

gp_tuning_pick_numeric_column <- function(x, choices) {
  gp_tuning_pick_numeric_column_info(x, choices)$values
}

gp_tuning_uncertainty_column_target <- function(column) {
  column <- as.character(column %||% NA_character_)[1L]
  if (!nzchar(column) || is.na(column)) {
    return(NA_character_)
  }
  if (column %in% c(
    "Prediction_Var_observed", "Prediction_variance",
    "Prediction_SE_observed", "Prediction_SE", "SE_observed", "SE"
  )) {
    return("observed")
  }
  if (column %in% c(
    "Prediction_Var_latent", "PEV",
    "Prediction_SE_latent", "SE_latent"
  )) {
    return("latent")
  }
  "unknown"
}

gp_tuning_uncertainty_calibration_limits <- function(min_scale = NULL,
                                                     max_scale = NULL,
                                                     prior_n = NULL,
                                                     prior_env = NULL) {
  if (is.null(min_scale)) {
    min_scale <- suppressWarnings(as.numeric(Sys.getenv(
      "PREDICTPRO_GP_SE_CALIBRATION_MIN_SCALE",
      unset = "0.5"
    ))[1L])
  }
  if (!is.finite(min_scale) || min_scale <= 0) {
    min_scale <- 0.5
  }
  if (is.null(max_scale)) {
    max_scale <- suppressWarnings(as.numeric(Sys.getenv(
      "PREDICTPRO_GP_SE_CALIBRATION_MAX_SCALE",
      unset = "2.5"
    ))[1L])
  }
  if (!is.finite(max_scale) || max_scale < min_scale) {
    max_scale <- max(min_scale, 2.5)
  }
  if (is.null(prior_n)) {
    prior_n <- suppressWarnings(as.numeric(Sys.getenv(
      "PREDICTPRO_GP_SE_CALIBRATION_PRIOR_N",
      unset = "50"
    ))[1L])
  }
  if (!is.finite(prior_n) || prior_n < 0) {
    prior_n <- 50
  }
  if (is.null(prior_env)) {
    prior_env <- suppressWarnings(as.numeric(Sys.getenv(
      "PREDICTPRO_GP_SE_CALIBRATION_PRIOR_ENV",
      unset = "3"
    ))[1L])
  }
  if (!is.finite(prior_env) || prior_env < 1) {
    prior_env <- 3
  }
  list(
    min_scale = as.numeric(min_scale),
    max_scale = as.numeric(max_scale),
    prior_n = as.numeric(prior_n),
    prior_env = as.numeric(prior_env)
  )
}

gp_tuning_regularize_uncertainty_scale <- function(raw_scale,
                                                   n_validation,
                                                   n_env,
                                                   min_scale = NULL,
                                                   max_scale = NULL,
                                                   prior_n = NULL,
                                                   prior_env = NULL) {
  limits <- gp_tuning_uncertainty_calibration_limits(
    min_scale = min_scale,
    max_scale = max_scale,
    prior_n = prior_n,
    prior_env = prior_env
  )
  raw_scale <- as.numeric(raw_scale)[1L]
  if (!is.finite(raw_scale) || raw_scale <= 0) {
    return(c(
      list(scale = NA_real_, shrinkage_weight = NA_real_),
      limits
    ))
  }
  n_validation <- max(0, as.numeric(n_validation)[1L] %||% 0)
  n_env <- max(1, as.numeric(n_env)[1L] %||% 1)
  n_weight <- n_validation / (n_validation + limits$prior_n)
  env_weight <- min(1, n_env / limits$prior_env)
  shrinkage_weight <- max(0, min(1, n_weight * env_weight))
  scale <- exp(log(raw_scale) * shrinkage_weight)
  scale <- min(limits$max_scale, max(limits$min_scale, scale))
  c(
    list(scale = as.numeric(scale), shrinkage_weight = as.numeric(shrinkage_weight)),
    limits
  )
}

gp_tuning_estimate_uncertainty_calibration <- function(fit,
                                                       truth,
                                                       offsets,
                                                       min_n = 20L,
                                                       min_scale = NULL,
                                                       max_scale = NULL,
                                                       prior_n = NULL,
                                                       prior_env = NULL) {
  pred <- gp_single_trait_prediction_keys(fit$predictions %||% data.frame())
  if (!all(c("GID", "Env", "Prediction") %in% names(pred))) {
    return(list(status = "skipped", reason = "prediction keys unavailable"))
  }
  offset_df <- data.frame(
    GID = truth$GID,
    Env = truth$Env,
    Prediction_Offset = as.numeric(offsets),
    stringsAsFactors = FALSE
  )
  joined <- merge(pred, truth, by = c("GID", "Env"), sort = FALSE)
  joined <- merge(joined, offset_df, by = c("GID", "Env"), sort = FALSE)
  if (!nrow(joined)) {
    return(list(status = "skipped", reason = "no validation predictions matched truth"))
  }
  joined$Prediction <- as.numeric(joined$Prediction) + as.numeric(joined$Prediction_Offset)
  pred_var_info <- gp_tuning_pick_numeric_column_info(
    joined,
    c("Prediction_Var_observed", "Prediction_variance", "Prediction_Var_latent", "PEV")
  )
  pred_var <- pred_var_info$values
  pred_se_info <- gp_tuning_pick_numeric_column_info(
    joined,
    c("Prediction_SE_observed", "SE", "SE_observed", "Prediction_SE", "Prediction_SE_latent", "SE_latent")
  )
  pred_se <- pred_se_info$values
  variance_from_se <- FALSE
  if (all(!is.finite(pred_var)) && any(is.finite(pred_se))) {
    pred_var <- pred_se^2
    pred_var_info$column <- pred_se_info$column
    variance_from_se <- TRUE
  }
  residual <- as.numeric(joined$Prediction) - as.numeric(joined$Observed)
  keep <- is.finite(residual) & is.finite(pred_var) & pred_var > 0
  min_n <- as.integer(min_n %||% 20L)[1L]
  if (!is.finite(min_n) || is.na(min_n) || min_n < 2L) {
    min_n <- 20L
  }
  if (sum(keep) < min_n) {
    return(list(
      status = "skipped",
      reason = paste0("too few validation rows with finite prediction variance (n=", sum(keep), ")"),
      n_validation_se = as.integer(sum(keep))
    ))
  }
  residual <- residual[keep]
  pred_var <- pred_var[keep]
  env <- as.character(joined$Env[keep])
  n_env <- length(unique(env))
  empirical_mse <- mean(residual^2)
  model_mse <- mean(pred_var)
  if (!is.finite(empirical_mse) || !is.finite(model_mse) || model_mse <= 0) {
    return(list(status = "skipped", reason = "non-finite validation uncertainty moments"))
  }
  global_scale <- sqrt(empirical_mse / model_mse)
  by_env <- split(seq_along(residual), env)
  env_scale <- vapply(by_env, function(idx) {
    if (length(idx) < 2L) {
      return(NA_real_)
    }
    denom <- mean(pred_var[idx])
    if (!is.finite(denom) || denom <= 0) {
      return(NA_real_)
    }
    sqrt(mean(residual[idx]^2) / denom)
  }, numeric(1L))
  env_scale <- env_scale[is.finite(env_scale) & env_scale > 0]
  robust_scale <- if (length(env_scale) >= 2L) stats::median(env_scale) else NA_real_
  raw_scale <- if (is.finite(robust_scale)) {
    exp(mean(log(c(global_scale, robust_scale))))
  } else {
    global_scale
  }
  regularized <- gp_tuning_regularize_uncertainty_scale(
    raw_scale,
    n_validation = length(residual),
    n_env = n_env,
    min_scale = min_scale,
    max_scale = max_scale,
    prior_n = prior_n,
    prior_env = prior_env
  )
  if (!is.finite(regularized$scale)) {
    return(list(status = "skipped", reason = "validation uncertainty scale was not finite"))
  }
  list(
    status = "ok",
    reason = NA_character_,
    scale = as.numeric(regularized$scale),
    variance_scale = as.numeric(regularized$scale)^2,
    raw_scale = as.numeric(raw_scale),
    global_scale = as.numeric(global_scale),
    env_median_scale = as.numeric(robust_scale),
    shrinkage_weight = as.numeric(regularized$shrinkage_weight),
    min_scale = as.numeric(regularized$min_scale),
    max_scale = as.numeric(regularized$max_scale),
    prior_n = as.numeric(regularized$prior_n),
    prior_env = as.numeric(regularized$prior_env),
    n_validation_se = as.integer(length(residual)),
    n_validation_env_se = as.integer(n_env),
    validation_rmse_for_se = sqrt(empirical_mse),
    validation_mean_model_se = sqrt(model_mse),
    variance_column = as.character(pred_var_info$column %||% NA_character_)[1L],
    se_column = as.character(pred_se_info$column %||% NA_character_)[1L],
    uncertainty_target = gp_tuning_uncertainty_column_target(pred_var_info$column),
    variance_from_se = isTRUE(variance_from_se),
    method = "blocked_validation_residual_mse_over_mean_prediction_variance"
  )
}

gp_tuning_fit_uncertainty_calibration <- function(selected,
                                                  mean_cache,
                                                  truth,
                                                  env_covariates = NULL,
                                                  env_kernel_cache,
                                                  include_components,
                                                  fixed_effects,
                                                  gp_backend,
                                                  reaction_norm_feature_qc,
                                                  kenv_bandwidth,
                                                  kenv_kernel_kwargs,
                                                  seed,
                                                  gp_se_hutchinson_probes = 32L,
                                                  python_bin = NULL,
                                                  project_root = NULL) {
  if (is.null(selected) || !nrow(selected)) {
    return(list(status = "skipped", reason = "no selected tuning candidate"))
  }
  candidate <- selected[1L, , drop = FALSE]
  mode <- as.character(candidate$mean_mode)[1L]
  prepped <- mean_cache[[mode]]
  if (is.null(prepped)) {
    return(list(status = "skipped", reason = "selected mean-adjustment cache unavailable"))
  }
  kernel_name <- as.character(candidate$kenv_kernel)[1L]
  cached_env <- env_kernel_cache$kernels[[kernel_name]]
  gp_se_hutchinson_probes <- as.integer(gp_se_hutchinson_probes %||% 32L)[1L]
  if (!is.finite(gp_se_hutchinson_probes) || is.na(gp_se_hutchinson_probes) ||
      gp_se_hutchinson_probes < 1L) {
    gp_se_hutchinson_probes <- 32L
  }
  fit_args <- list(
    pheno_data = prepped$fit_df,
    gmatrix = prepped$geno_kernel,
    response = "y",
    gen_name = "GID",
    heter_groups = "Env",
    include_components = include_components,
    fixed_effects = fixed_effects,
    gp_backend = gp_backend,
    gp_output_level = "predict_with_se",
    gp_return_se = TRUE,
    gp_prediction_output = "test_only",
    krr_lam = as.numeric(candidate$lambda)[1L],
    krr_lams = "none",
    lam_select = "fixed",
    w_g = as.numeric(candidate$w_g)[1L],
    w_ge = as.numeric(candidate$w_ge)[1L],
    w_e = as.numeric(candidate$w_e)[1L],
    seed = as.integer(seed),
    python_bin = python_bin,
    project_root = project_root,
    tune_gp = FALSE,
    mean_adjustment = "none",
    prediction_blend = "none",
    gp_auto_policy = FALSE,
    gp_se_hutchinson_probes = gp_se_hutchinson_probes
  )
  if (!is.null(cached_env)) {
    fit_args$env_similarity <- cached_env
  } else if (!is.null(env_covariates)) {
    fit_args$env_covariates <- env_covariates
    fit_args$reaction_norm_feature_qc <- isTRUE(reaction_norm_feature_qc)
    fit_args$kenv_kernel <- kernel_name
    fit_args$kenv_bandwidth <- as.numeric(kenv_bandwidth)[1L]
    if (!is.null(kenv_kernel_kwargs)) {
      fit_args$kenv_kernel_kwargs <- kenv_kernel_kwargs
    }
  }
  fit <- tryCatch(
    suppressWarnings(do.call(gp_single_trait_model, fit_args)),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      status = "skipped",
      reason = paste("validation SE fit failed:", gp_tuning_error_text(fit))
    ))
  }
  out <- gp_tuning_estimate_uncertainty_calibration(
    fit = fit,
    truth = truth,
    offsets = prepped$validation_offsets
  )
  rm(fit)
  if (requireNamespace("reticulate", quietly = TRUE)) {
    try({
      py_gc <- reticulate::import("gc", delay_load = FALSE)
      py_gc$collect()
    }, silent = TRUE)
  }
  invisible(gc(verbose = FALSE))
  out$candidate_id <- suppressWarnings(as.integer(candidate$candidate_id[[1L]] %||% NA_integer_))
  out
}

gp_tuning_env_kernel_detail_row <- function(kernel,
                                            status,
                                            info = list(),
                                            error = NA_character_) {
  data.frame(
    kenv_kernel = as.character(kernel)[1L],
    status = as.character(status)[1L],
    n_env = as.integer(info$n_env %||% NA_integer_),
    n_features_input = as.integer(info$n_features_input %||% NA_integer_),
    n_features_after_qc = as.integer(info$n_features_after_qc %||% NA_integer_),
    n_features_after_variance_filter = as.integer(info$n_features_after_variance_filter %||% NA_integer_),
    g2f_schema_detected = isTRUE(info$g2f_schema_detected),
    n_envs_imputed = as.integer(info$n_envs_imputed %||% NA_integer_),
    error = as.character(error %||% NA_character_)[1L],
    stringsAsFactors = FALSE
  )
}

gp_tuning_env_kernel_cache_work <- function(env_covariates, env_levels, kenv_kernels) {
  n_features <- NA_integer_
  if (is.data.frame(env_covariates)) {
    numeric_cols <- vapply(env_covariates, is.numeric, logical(1L))
    n_features <- sum(numeric_cols)
  } else if (is.matrix(env_covariates) || is.array(env_covariates)) {
    n_features <- ncol(env_covariates)
  }
  if (!is.finite(n_features) || is.na(n_features)) {
    return(Inf)
  }
  as.numeric(length(env_levels)) * as.numeric(n_features) * as.numeric(length(kenv_kernels))
}

gp_tuning_prepare_env_kernel_cache <- function(env_covariates,
                                               env_levels,
                                               kenv_kernels,
                                               reaction_norm_feature_qc = TRUE,
                                               kenv_bandwidth = 1.0,
                                               kenv_kernel_kwargs = NULL,
                                               python_bin = NULL,
                                               project_root = NULL) {
  empty <- list(
    available = FALSE,
    kernels = list(),
    detail = data.frame(),
    reason = NA_character_
  )
  if (is.null(env_covariates) || !length(env_levels)) {
    empty$reason <- "env_covariates or env_levels unavailable"
    return(empty)
  }
  if (exists("gp_bridge_direct_execution_enabled", mode = "function") &&
      isTRUE(gp_bridge_direct_execution_enabled())) {
    empty$reason <- "direct Python subprocess execution builds environment kernels inside fit calls"
    return(empty)
  }

  kernels <- unique(tolower(as.character(kenv_kernels)))
  kernels <- kernels[nzchar(kernels)]
  if (!length(kernels)) {
    empty$reason <- "no environment kernels requested"
    return(empty)
  }
  min_cache_work <- suppressWarnings(as.numeric(
    Sys.getenv("PREDICTPRO_GP_ENV_KERNEL_CACHE_MIN_WORK", unset = "500")
  )[1L])
  if (!is.finite(min_cache_work)) {
    min_cache_work <- 500
  }
  cache_work <- gp_tuning_env_kernel_cache_work(env_covariates, env_levels, kernels)
  if (is.finite(cache_work) && cache_work < min_cache_work) {
    empty$reason <- paste0(
      "small environment-covariate kernel workload left uncached (work=",
      round(cache_work, 2),
      ", threshold=",
      round(min_cache_work, 2),
      ")"
    )
    return(empty)
  }

  fw <- tryCatch(
    gp_import_gp_module("gp_framework", python_bin = python_bin, project_root = project_root),
    error = function(e) e
  )
  if (inherits(fw, "error")) {
    empty$reason <- conditionMessage(fw)
    return(empty)
  }
  has_helper <- tryCatch(
    reticulate::py_has_attr(fw, "prepare_env_similarity_from_covariates"),
    error = function(e) FALSE
  )
  if (!isTRUE(has_helper)) {
    empty$reason <- "Python helper prepare_env_similarity_from_covariates is unavailable"
    return(empty)
  }

  out <- list()
  detail <- vector("list", length(kernels))
  names(detail) <- kernels
  for (kernel_name in kernels) {
    built <- tryCatch(
      reticulate::py_to_r(fw$prepare_env_similarity_from_covariates(
        env_covariates = env_covariates,
        env_levels = as.list(as.character(env_levels)),
        env_col = "Env",
        feature_qc = isTRUE(reaction_norm_feature_qc),
        kenv_kernel = kernel_name,
        kenv_bandwidth = as.numeric(kenv_bandwidth)[1L],
        kenv_kernel_kwargs = kenv_kernel_kwargs
      )),
      error = function(e) e
    )
    if (inherits(built, "error")) {
      detail[[kernel_name]] <- gp_tuning_env_kernel_detail_row(
        kernel_name,
        status = "failed",
        error = conditionMessage(built)
      )
      next
    }
    info <- built$info %||% list()
    K_env <- built$env_similarity
    if (is.null(K_env)) {
      detail[[kernel_name]] <- gp_tuning_env_kernel_detail_row(
        kernel_name,
        status = "empty",
        info = info
      )
      next
    }
    K_env <- as.matrix(K_env)
    storage.mode(K_env) <- "double"
    rownames(K_env) <- colnames(K_env) <- as.character(env_levels)
    out[[kernel_name]] <- K_env
    detail[[kernel_name]] <- gp_tuning_env_kernel_detail_row(
      kernel_name,
      status = "ok",
      info = info
    )
  }

  list(
    available = length(out) > 0L,
    kernels = out,
    detail = if (length(detail)) do.call(rbind, detail) else data.frame(),
    reason = if (length(out)) NA_character_ else "no environment kernels were cached"
  )
}

gp_tuning_error_text <- function(e) {
  msg <- conditionMessage(e)
  if (requireNamespace("reticulate", quietly = TRUE)) {
    py <- tryCatch(reticulate::py_last_error(), error = function(err) NULL)
    tb <- py$traceback %||% py$message %||% NULL
    if (!is.null(tb) && length(tb) && nzchar(as.character(tb)[1L])) {
      msg <- paste(msg, paste(as.character(tb), collapse = "\n"), sep = "\n")
    }
  }
  msg
}

gp_tuning_empty_score <- function(error = NA_character_) {
  list(
    rmse = Inf,
    mae = Inf,
    cor = NA_real_,
    centered_rmse = NA_real_,
    centered_mae = NA_real_,
    centered_cor = NA_real_,
    within_env_spearman = NA_real_,
    baseline_rmse = NA_real_,
    baseline_cor = NA_real_,
    baseline_centered_rmse = NA_real_,
    baseline_centered_cor = NA_real_,
    baseline_within_env_spearman = NA_real_,
    blend_alpha = NA_real_,
    blend_objective = NA_character_,
    blend_objective_value = NA_real_,
    blend_rmse = NA_real_,
    blend_cor = NA_real_,
    blend_centered_rmse = NA_real_,
    blend_centered_cor = NA_real_,
    blend_within_env_spearman = NA_real_,
    rmse_se = NA_real_,
    centered_rmse_se = NA_real_,
    within_env_spearman_se = NA_real_,
    validation_uncertainty_scale = NA_real_,
    validation_uncertainty_raw_scale = NA_real_,
    validation_uncertainty_n = 0L,
    validation_uncertainty_n_env = 0L,
    n_validation = 0L,
    env_detail = NULL,
    failed = TRUE,
    error = as.character(error %||% NA_character_)[1L]
  )
}

gp_tuning_batch_candidates_direct <- function(prepped,
                                              candidate_df,
                                              env_covariates = NULL,
                                              env_kernel_cache = NULL,
                                              include_components,
                                              gp_backend,
                                              reaction_norm_feature_qc,
                                              kenv_bandwidth,
                                              kenv_kernel_kwargs,
                                              seed,
                                              python_bin = NULL,
                                              project_root = NULL) {
  python_bin <- python_bin %||% gp_detect_gp_python()
  project_root <- project_root %||% gp_project_root()
  dir_path <- tempfile("predictpror_gp_tune_batch_")
  input_dir <- file.path(dir_path, "input")
  out_dir <- file.path(dir_path, "out")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pheno_csv <- file.path(input_dir, "pheno.csv")
  kernel_bin <- file.path(input_dir, "geno_kernel.bin")
  kernel_meta <- file.path(input_dir, "geno_kernel_shape.txt")
  geno_ids_path <- file.path(input_dir, "geno_ids.txt")
  train_idx_path <- file.path(input_dir, "train_idx.txt")
  test_idx_path <- file.path(input_dir, "test_idx.txt")
  candidates_csv <- file.path(input_dir, "candidates.csv")
  utils::write.csv(prepped$fit_df, pheno_csv, row.names = FALSE, quote = TRUE)
  gp_bridge_write_matrix_bin(prepped$geno_kernel, kernel_bin, kernel_meta)
  gp_bridge_write_lines(prepped$ids, geno_ids_path)
  gp_bridge_write_lines(as.integer(prepped$training_idx), train_idx_path)
  gp_bridge_write_lines(as.integer(prepped$validation_idx), test_idx_path)
  utils::write.csv(candidate_df, candidates_csv, row.names = FALSE, quote = TRUE)

  spec <- list(
    project_root = normalizePath(project_root, winslash = "/", mustWork = TRUE),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    pheno_csv = normalizePath(pheno_csv, winslash = "/", mustWork = TRUE),
    kernel_bin = normalizePath(kernel_bin, winslash = "/", mustWork = TRUE),
    kernel_meta = normalizePath(kernel_meta, winslash = "/", mustWork = TRUE),
    geno_ids = normalizePath(geno_ids_path, winslash = "/", mustWork = TRUE),
    train_idx = normalizePath(train_idx_path, winslash = "/", mustWork = TRUE),
    test_idx = normalizePath(test_idx_path, winslash = "/", mustWork = TRUE),
    candidates_csv = normalizePath(candidates_csv, winslash = "/", mustWork = TRUE),
    gid_col = "GID",
    env_col = "Env",
    y_col = "y",
    include_components = as.list(include_components),
    fixed_effects = list(),
    backend = as.character(gp_backend)[1L],
    reaction_norm_feature_qc = isTRUE(reaction_norm_feature_qc),
    kenv_bandwidth = as.numeric(kenv_bandwidth)[1L],
    kenv_kernel_kwargs = gp_bridge_jsonable_value(kenv_kernel_kwargs),
    dtype = "float64",
    seed = as.integer(seed)[1L],
    standardize = "global",
    prediction_output = "test_only"
  )
  if (!is.null(env_covariates)) {
    env_covariates_csv <- file.path(input_dir, "env_covariates.csv")
    utils::write.csv(env_covariates, env_covariates_csv, row.names = FALSE, quote = TRUE)
    spec$env_covariates_csv <- normalizePath(env_covariates_csv, winslash = "/", mustWork = TRUE)
  }
  kernels <- env_kernel_cache$kernels %||% NULL
  if (length(kernels)) {
    env_specs <- list()
    for (kernel_name in names(kernels)) {
      if (!nzchar(as.character(kernel_name))) {
        next
      }
      env_bin <- file.path(input_dir, paste0("env_kernel_", make.names(kernel_name), ".bin"))
      env_meta <- file.path(input_dir, paste0("env_kernel_", make.names(kernel_name), "_shape.txt"))
      gp_bridge_write_matrix_bin(kernels[[kernel_name]], env_bin, env_meta)
      env_specs[[as.character(kernel_name)]] <- list(
        bin = normalizePath(env_bin, winslash = "/", mustWork = TRUE),
        meta = normalizePath(env_meta, winslash = "/", mustWork = TRUE)
      )
    }
    if (length(env_specs)) {
      spec$env_similarity_by_kernel <- env_specs
    }
  }

  spec_json <- file.path(input_dir, "batch_candidates_spec.json")
  gp_bridge_write_json(spec, spec_json)
  worker_ok <- FALSE
  if (gp_bridge_worker_use_enabled()) {
    worker_ok <- tryCatch({
      gp_bridge_worker_request(
        "batch-fit-candidates-spec",
        list(
          spec = normalizePath(spec_json, winslash = "/", mustWork = TRUE),
          out_dir = NULL,
          project_root = project_root
        ),
        python_bin = python_bin,
        project_root = project_root
      )
      TRUE
    }, error = function(e) {
      warning("Warm GP worker failed for batch tuning candidates; falling back to subprocess: ",
              conditionMessage(e), call. = FALSE)
      FALSE
    })
  }
  if (!isTRUE(worker_ok)) {
    gp_bridge_run_cli(
      "batch-fit-candidates-spec",
      c("--spec", normalizePath(spec_json, winslash = "/", mustWork = TRUE)),
      python_bin = python_bin
    )
  }

  pred_path <- file.path(out_dir, "predictions.csv")
  detail_path <- file.path(out_dir, "detail.csv")
  env_info_path <- file.path(out_dir, "env_kernel_cache_info.csv")
  if (!file.exists(pred_path) || !file.exists(detail_path)) {
    stop("GP batch tuning bridge did not write predictions/detail outputs.", call. = FALSE)
  }
  list(
    predictions = utils::read.csv(pred_path, stringsAsFactors = FALSE, check.names = FALSE),
    detail = utils::read.csv(detail_path, stringsAsFactors = FALSE, check.names = FALSE),
    env_kernel_cache_info = if (file.exists(env_info_path) &&
        is.finite(file.info(env_info_path)$size) &&
        file.info(env_info_path)$size > 1L) {
      tryCatch(
        utils::read.csv(env_info_path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) data.frame()
      )
    } else {
      data.frame()
    },
    info = gp_bridge_read_json(file.path(out_dir, "meta.json")),
    bridge_dir = normalizePath(dir_path, winslash = "/", mustWork = FALSE)
  )
}

gp_tuning_try_batch_candidate_scores <- function(tune_grid,
                                                 mean_cache,
                                                 truth,
                                                 baseline = NULL,
                                                 blend_objective = "rmse",
                                                 env_covariates = NULL,
                                                 env_kernel_cache,
                                                 include_components,
                                                 gp_backend,
                                                 reaction_norm_feature_qc,
                                                 kenv_bandwidth,
                                                 kenv_kernel_kwargs,
                                                 seed,
                                                 python_bin = NULL,
                                                 project_root = NULL) {
  use_direct <- exists("gp_bridge_direct_execution_enabled", mode = "function") &&
    isTRUE(gp_bridge_direct_execution_enabled())
  has_env_kernel_cache <- length(env_kernel_cache$kernels %||% list()) > 0L
  if (is.null(env_covariates) && !has_env_kernel_cache) {
    return(NULL)
  }
  fw <- NULL
  if (!isTRUE(use_direct)) {
    fw <- tryCatch(
      gp_import_gp_module("gp_framework", python_bin = python_bin, project_root = project_root),
      error = function(e) e
    )
    if (inherits(fw, "error")) {
      return(NULL)
    }
    has_helper <- tryCatch(
      reticulate::py_has_attr(fw, "batch_fit_mixed_model_candidates"),
      error = function(e) FALSE
    )
    if (!isTRUE(has_helper)) {
      return(NULL)
    }
  }

  details <- vector("list", nrow(tune_grid))
  env_details <- list()
  env_cache_info <- list()
  for (mode in names(mean_cache)) {
    idx <- which(as.character(tune_grid$mean_mode) == mode)
    if (!length(idx)) {
      next
    }
    prepped <- mean_cache[[mode]]
    candidate_df <- data.frame(
      candidate_id = idx,
      tune_grid[idx, , drop = FALSE],
      stringsAsFactors = FALSE
    )
    py_out <- tryCatch({
      if (isTRUE(use_direct)) {
        gp_tuning_batch_candidates_direct(
          prepped = prepped,
          candidate_df = candidate_df,
          env_covariates = env_covariates,
          env_kernel_cache = env_kernel_cache,
          include_components = include_components,
          gp_backend = gp_backend,
          reaction_norm_feature_qc = reaction_norm_feature_qc,
          kenv_bandwidth = kenv_bandwidth,
          kenv_kernel_kwargs = kenv_kernel_kwargs,
          seed = seed,
          python_bin = python_bin,
          project_root = project_root
        )
      } else {
        reticulate::py_to_r(fw$batch_fit_mixed_model_candidates(
          pheno_df = prepped$fit_df,
          gid_col = "GID",
          env_col = "Env",
          y_col = "y",
          geno_kernel = unname(prepped$geno_kernel),
          geno_ids = as.list(prepped$ids),
          train_idx = as.integer(prepped$training_idx),
          test_idx = as.integer(prepped$validation_idx),
          candidates = candidate_df,
          env_covariates = env_covariates,
          env_similarity_by_kernel = env_kernel_cache$kernels,
          include_components = as.list(include_components),
          fixed_effects = list(),
          backend = gp_backend,
          reaction_norm_feature_qc = isTRUE(reaction_norm_feature_qc),
          kenv_bandwidth = as.numeric(kenv_bandwidth)[1L],
          kenv_kernel_kwargs = kenv_kernel_kwargs,
          dtype = "float64",
          seed = as.integer(seed),
          standardize = "global",
          prediction_output = "test_only"
        ))
      }
    }, error = function(e) e)
    if (inherits(py_out, "error")) {
      return(NULL)
    }
    if (!is.null(py_out$env_kernel_cache_info) &&
        nrow(gp_plain_data_frame(py_out$env_kernel_cache_info))) {
      env_cache_info[[length(env_cache_info) + 1L]] <- gp_plain_data_frame(
        py_out$env_kernel_cache_info
      )
    }
    pred_all <- gp_plain_data_frame(py_out$predictions %||% data.frame())
    py_detail <- gp_plain_data_frame(py_out$detail %||% data.frame())
    pred_by_candidate <- gp_tuning_rows_by_candidate_id(pred_all)
    detail_by_candidate <- gp_tuning_rows_by_candidate_id(py_detail)
    for (candidate_id in idx) {
      candidate <- tune_grid[candidate_id, , drop = FALSE]
      candidate_key <- as.character(candidate_id)
      row_detail <- detail_by_candidate[[candidate_key]]
      if (is.null(row_detail)) {
        row_detail <- py_detail[FALSE, , drop = FALSE]
      }
      failed <- if (nrow(row_detail)) isTRUE(as.logical(row_detail$failed[[1L]])) else TRUE
      error <- if (nrow(row_detail)) as.character(row_detail$error[[1L]] %||% NA_character_) else "batch candidate missing"
      elapsed <- if (nrow(row_detail)) suppressWarnings(as.numeric(row_detail$elapsed_sec[[1L]])) else NA_real_
      score <- if (failed) {
        gp_tuning_empty_score(error)
      } else {
        pred <- pred_by_candidate[[candidate_key]]
        if (is.null(pred)) {
          pred <- pred_all[FALSE, , drop = FALSE]
        }
        if ("candidate_id" %in% names(pred)) {
          pred$candidate_id <- NULL
        }
        tryCatch({
          out <- gp_tuning_score_predictions(
            list(predictions = pred),
            truth,
            prepped$validation_offsets,
            baseline = baseline,
            blend_objective = blend_objective
          )
          out$failed <- FALSE
          out$error <- NA_character_
          out
        }, error = function(e) gp_tuning_empty_score(gp_tuning_error_text(e)))
      }
      if (!is.null(score$env_detail)) {
        env_detail <- score$env_detail
        env_detail$candidate_id <- candidate_id
        env_details[[length(env_details) + 1L]] <- env_detail
      }
      details[[candidate_id]] <- data.frame(
        candidate_id = candidate_id,
        candidate,
        validation_rmse = score$rmse,
        validation_mae = score$mae,
        validation_cor = score$cor,
        validation_centered_rmse = score$centered_rmse,
        validation_centered_mae = score$centered_mae,
        validation_centered_cor = score$centered_cor,
        validation_within_env_spearman = score$within_env_spearman,
        validation_envmean_rmse = score$baseline_rmse,
        validation_envmean_cor = score$baseline_cor,
        validation_envmean_centered_rmse = score$baseline_centered_rmse,
        validation_envmean_centered_cor = score$baseline_centered_cor,
        validation_envmean_within_env_spearman = score$baseline_within_env_spearman,
        validation_blend_alpha = score$blend_alpha,
        validation_blend_objective = score$blend_objective,
        validation_blend_objective_value = score$blend_objective_value,
        validation_blend_rmse = score$blend_rmse,
        validation_blend_cor = score$blend_cor,
        validation_blend_centered_rmse = score$blend_centered_rmse,
        validation_blend_centered_cor = score$blend_centered_cor,
        validation_blend_within_env_spearman = score$blend_within_env_spearman,
        validation_rmse_se = score$rmse_se,
        validation_centered_rmse_se = score$centered_rmse_se,
        validation_within_env_spearman_se = score$within_env_spearman_se,
        validation_uncertainty_scale = score$validation_uncertainty_scale %||% NA_real_,
        validation_uncertainty_raw_scale = score$validation_uncertainty_raw_scale %||% NA_real_,
        validation_uncertainty_n = score$validation_uncertainty_n %||% 0L,
        validation_uncertainty_n_env = score$validation_uncertainty_n_env %||% 0L,
        n_validation = score$n_validation,
        complexity = gp_tuning_candidate_complexity(candidate),
        elapsed_sec = if (is.finite(elapsed)) elapsed else NA_real_,
        failed = score$failed,
        error = score$error,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!all(vapply(details, is.data.frame, logical(1L)))) {
    return(NULL)
  }
  list(
    tuning_detail = do.call(rbind, details),
    tuning_env_detail = if (length(env_details)) do.call(rbind, env_details) else data.frame(),
    env_kernel_cache_detail = if (length(env_cache_info)) {
      unique(do.call(rbind, env_cache_info))
    } else {
      env_kernel_cache$detail
    },
    batch_used = TRUE,
    batch_execution = if (isTRUE(use_direct)) "python_subprocess" else "reticulate"
  )
}

gp_tune_single_trait_met <- function(pheno_df,
                                     split,
                                     kernel,
                                     gp_factor_cache = NULL,
                                     geno_ids = NULL,
                                     env_covariates = NULL,
                                     reaction_norm_feature_qc = TRUE,
                                     kenv_bandwidth = 1.0,
                                     kenv_kernel_kwargs = NULL,
                                     include_components = c("g", "ge", "e"),
                                     gp_backend = "auto",
                                     tuning_krr_lams = c(0.03, 0.1, 0.3),
                                     tuning_kenv_kernels = c("matern32", "matern52", "rbf"),
                                     tuning_weight_grid = NULL,
                                     tuning_objective = "rmse",
                                     mean_adjustment = "auto",
                                     fixed_effects = character(0),
                                     observation_weights = NULL,
                                     env_location = NULL,
                                     env_year = NULL,
                                     tuning_validation_fraction = 0.25,
                                     tuning_strategy = "blocked_environment",
                                     tuning_max_calibration_rows = NULL,
                                     tuning_max_validation_rows = NULL,
                                     tuning_max_genotypes = NULL,
                                     calibrate_uncertainty = FALSE,
                                     gp_se_hutchinson_probes = 32L,
                                     seed = 12345L,
                                     python_bin = NULL,
                                     project_root = NULL) {
  if (!is.null(observation_weights)) {
    observation_weights <- suppressWarnings(as.numeric(observation_weights))
    if (length(observation_weights) != nrow(pheno_df) ||
        any(!is.finite(observation_weights)) || any(observation_weights <= 0)) {
      stop(
        "GP tuning Stage 2 observation weights must contain one finite positive precision per phenotype row.",
        call. = FALSE
      )
    }
  }
  tuning_strategy <- gp_match_gp_choice(
    tuning_strategy,
    choices = c("blocked_environment", "blocked_genotype"),
    name = "tuning_strategy",
    default = "blocked_environment"
  )
  blocked <- if (identical(tuning_strategy, "blocked_genotype")) {
    gp_tuning_make_genotype_split(
      pheno_df,
      train_mask = split$train_mask,
      env_location = env_location,
      env_year = env_year,
      validation_fraction = tuning_validation_fraction,
      seed = seed,
      allow_single_env = length(unique(as.character(pheno_df$Env[split$train_mask]))) == 1L
    )
  } else {
    gp_tuning_make_blocked_split(
      pheno_df,
      train_mask = split$train_mask,
      env_location = env_location,
      env_year = env_year,
      validation_fraction = tuning_validation_fraction,
      seed = seed
    )
  }
  if (!identical(blocked$status, "ok")) {
    return(list(
      status = "skipped",
      reason = blocked$reason %||% "blocked tuning split was unavailable",
      tuning_objective = gp_tuning_normalize_objective(tuning_objective),
      tuning_strategy = tuning_strategy,
      selected = NULL,
      tuning_detail = data.frame(),
      tuning_env_detail = data.frame(),
      uncertainty_calibration = list(status = "skipped", reason = blocked$reason %||% "blocked tuning split was unavailable")
    ))
  }

  bounded <- gp_tuning_subset_blocked_rows(
    pheno_df,
    calibration_rows = blocked$calibration_rows,
    validation_rows = blocked$validation_rows,
    max_calibration_rows = tuning_max_calibration_rows,
    max_validation_rows = tuning_max_validation_rows,
    max_genotypes = tuning_max_genotypes,
    seed = seed
  )
  if (!identical(bounded$status, "ok")) {
    return(list(
      status = "skipped",
      reason = bounded$reason %||% "bounded tuning subset was unavailable",
      tuning_objective = gp_tuning_normalize_objective(tuning_objective),
      tuning_strategy = tuning_strategy,
      selected = NULL,
      tuning_detail = data.frame(),
      tuning_env_detail = data.frame(),
      uncertainty_calibration = list(status = "skipped", reason = bounded$reason %||% "bounded tuning subset was unavailable"),
      tuning_sample = bounded$sample_info
    ))
  }
  blocked$calibration_rows <- bounded$calibration_rows
  blocked$validation_rows <- bounded$validation_rows
  blocked$calibration_ids <- unique(as.character(pheno_df$GID[blocked$calibration_rows]))
  blocked$validation_ids <- unique(as.character(pheno_df$GID[blocked$validation_rows]))
  blocked$calibration_envs <- unique(as.character(pheno_df$Env[blocked$calibration_rows]))
  blocked$validation_envs <- unique(as.character(pheno_df$Env[blocked$validation_rows]))

  tuning_objective <- gp_tuning_normalize_objective(tuning_objective)

  mean_adjustment <- gp_match_gp_choice(
    mean_adjustment,
    choices = c("auto", "none", "location_offset"),
    name = "mean_adjustment",
    default = "auto"
  )
  mean_modes <- switch(
    mean_adjustment,
    auto = c("none", "location_offset"),
    none = "none",
    location_offset = "location_offset"
  )
  has_env_covariates <- !is.null(env_covariates)
  tuning_kenv_kernels_eff <- if (isTRUE(has_env_covariates)) {
    tuning_kenv_kernels
  } else {
    "linear"
  }
  tuning_krr_lams_eff <- if (isTRUE(has_env_covariates)) {
    tuning_krr_lams
  } else {
    unique(c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, as.numeric(tuning_krr_lams)))
  }
  tuning_weight_grid_eff <- if (isTRUE(has_env_covariates)) {
    tuning_weight_grid
  } else {
    data.frame(w_g = 1, w_ge = 0, w_e = 0, stringsAsFactors = FALSE)
  }
  include_components_eff <- if (isTRUE(has_env_covariates)) {
    include_components
  } else {
    "g"
  }
  fixed_effects_eff <- if (isTRUE(has_env_covariates)) {
    character(0)
  } else {
    as.character(fixed_effects %||% character(0))
  }
  tune_grid <- gp_tuning_expand_grid(
    tuning_krr_lams = tuning_krr_lams_eff,
    tuning_kenv_kernels = tuning_kenv_kernels_eff,
    tuning_weight_grid = tuning_weight_grid_eff,
    mean_modes = mean_modes
  )
  if (isTRUE(has_env_covariates)) {
    tune_grid <- gp_tuning_effective_reaction_norm_grid(tune_grid)
  }

  train_df <- pheno_df[blocked$calibration_rows, c("GID", "Env", "y"), drop = FALSE]
  validation_df <- pheno_df[blocked$validation_rows, c("GID", "Env", "y"), drop = FALSE]
  validation_baseline <- gp_tuning_envmean_baseline(
    training_df = train_df,
    target_df = validation_df,
    env_location = env_location,
    env_year = env_year
  )
  truth <- data.frame(
    GID = as.character(validation_df$GID),
    Env = as.character(validation_df$Env),
    Observed = as.numeric(validation_df$y),
    stringsAsFactors = FALSE
  )
  fit_df0 <- rbind(train_df, validation_df)
  fit_weights0 <- if (is.null(observation_weights)) {
    NULL
  } else {
    observation_weights[c(blocked$calibration_rows, blocked$validation_rows)]
  }
  training_rows <- seq_len(nrow(train_df))
  validation_rows <- seq.int(nrow(train_df) + 1L, nrow(fit_df0))
  mean_mode_values <- unique(as.character(tune_grid$mean_mode))
  mean_cache <- stats::setNames(vector("list", length(mean_mode_values)), mean_mode_values)
  for (mode in mean_mode_values) {
    adjusted <- gp_tuning_apply_mean_adjustment(
      fit_df0,
      training_rows = training_rows,
      mode = mode,
      env_location = env_location,
      env_year = env_year
    )
    fit_df <- adjusted$pheno_df
    fit_df$y[validation_rows] <- NA_real_
    ids <- unique(as.character(fit_df$GID))
    mean_cache[[mode]] <- list(
      fit_df = fit_df,
      validation_offsets = adjusted$offsets[validation_rows],
      ids = ids,
      geno_kernel = gp_tuning_dense_kernel_for_ids(
        kernel,
        ids = ids,
        gp_factor_cache = gp_factor_cache,
        geno_ids = geno_ids,
        python_bin = python_bin,
        project_root = project_root
      ),
      observation_weights = fit_weights0,
      training_idx = as.integer(training_rows - 1L),
      validation_idx = as.integer(validation_rows - 1L)
    )
  }
  env_levels <- sort(unique(as.character(fit_df0$Env)))
  env_kernel_cache <- gp_tuning_prepare_env_kernel_cache(
    env_covariates = env_covariates,
    env_levels = env_levels,
    kenv_kernels = unique(as.character(tune_grid$kenv_kernel)),
    reaction_norm_feature_qc = reaction_norm_feature_qc,
    kenv_bandwidth = kenv_bandwidth,
    kenv_kernel_kwargs = kenv_kernel_kwargs,
    python_bin = python_bin,
    project_root = project_root
  )

  # The batched tuning bridge does not yet carry observation-level residual
  # precisions. Weighted Stage 2 tuning therefore uses the serial candidate
  # path below so weights cannot be silently dropped.
  batch_scores <- if (is.null(observation_weights)) {
    gp_tuning_try_batch_candidate_scores(
      tune_grid = tune_grid,
      mean_cache = mean_cache,
      truth = truth,
      baseline = validation_baseline,
      blend_objective = tuning_objective,
      env_covariates = env_covariates,
      env_kernel_cache = env_kernel_cache,
      include_components = include_components_eff,
      gp_backend = gp_backend,
      reaction_norm_feature_qc = reaction_norm_feature_qc,
      kenv_bandwidth = kenv_bandwidth,
      kenv_kernel_kwargs = kenv_kernel_kwargs,
      seed = seed,
      python_bin = python_bin,
      project_root = project_root
    )
  } else {
    NULL
  }

  batch_candidate_used <- !is.null(batch_scores)
  batch_candidate_execution <- if (batch_candidate_used) {
    as.character(batch_scores$batch_execution %||% NA_character_)[1L]
  } else {
    NA_character_
  }

  if (!is.null(batch_scores)) {
    tuning_detail <- batch_scores$tuning_detail
    tuning_env_detail <- batch_scores$tuning_env_detail
    env_kernel_cache_detail <- batch_scores$env_kernel_cache_detail
  } else {
    details <- vector("list", nrow(tune_grid))
    env_details <- list()
    for (i in seq_len(nrow(tune_grid))) {
      candidate <- tune_grid[i, , drop = FALSE]
      t0 <- proc.time()[["elapsed"]]
      score <- tryCatch({
        mode <- as.character(candidate$mean_mode)[1L]
        kernel_name <- as.character(candidate$kenv_kernel)[1L]
        prepped <- mean_cache[[mode]]
        cached_env <- env_kernel_cache$kernels[[kernel_name]]
        fit_args <- list(
          pheno_data = prepped$fit_df,
          gmatrix = prepped$geno_kernel,
          response = "y",
          gen_name = "GID",
          heter_groups = "Env",
        include_components = include_components_eff,
        fixed_effects = fixed_effects_eff,
          gp_backend = gp_backend,
          gp_output_level = "predict_only",
          gp_prediction_output = "test_only",
          krr_lam = as.numeric(candidate$lambda)[1L],
          krr_lams = "none",
          lam_select = "fixed",
          w_g = as.numeric(candidate$w_g)[1L],
          w_ge = as.numeric(candidate$w_ge)[1L],
          w_e = as.numeric(candidate$w_e)[1L],
          seed = as.integer(seed),
          python_bin = python_bin,
          project_root = project_root,
          tune_gp = FALSE,
          mean_adjustment = "none",
          observation_weights = prepped$observation_weights,
          gp_auto_policy = FALSE
        )
        if (!is.null(cached_env)) {
          fit_args$env_similarity <- cached_env
        } else if (!is.null(env_covariates)) {
          fit_args$env_covariates <- env_covariates
          fit_args$reaction_norm_feature_qc <- isTRUE(reaction_norm_feature_qc)
          fit_args$kenv_kernel <- kernel_name
          fit_args$kenv_bandwidth <- as.numeric(kenv_bandwidth)[1L]
          if (!is.null(kenv_kernel_kwargs)) {
            fit_args$kenv_kernel_kwargs <- kenv_kernel_kwargs
          }
        }
        fit <- suppressWarnings(do.call(gp_single_trait_model, fit_args))
        out <- gp_tuning_score_predictions(
          fit,
          truth,
          prepped$validation_offsets,
          baseline = validation_baseline,
          blend_objective = tuning_objective
        )
        out$failed <- FALSE
        out$error <- NA_character_
        out
      }, error = function(e) gp_tuning_empty_score(gp_tuning_error_text(e)))

      if (!is.null(score$env_detail)) {
        env_detail <- score$env_detail
        env_detail$candidate_id <- i
        env_details[[length(env_details) + 1L]] <- env_detail
      }
      details[[i]] <- data.frame(
        candidate_id = i,
        candidate,
        validation_rmse = score$rmse,
        validation_mae = score$mae,
        validation_cor = score$cor,
        validation_centered_rmse = score$centered_rmse,
        validation_centered_mae = score$centered_mae,
        validation_centered_cor = score$centered_cor,
        validation_within_env_spearman = score$within_env_spearman,
        validation_envmean_rmse = score$baseline_rmse,
        validation_envmean_cor = score$baseline_cor,
        validation_envmean_centered_rmse = score$baseline_centered_rmse,
        validation_envmean_centered_cor = score$baseline_centered_cor,
        validation_envmean_within_env_spearman = score$baseline_within_env_spearman,
        validation_blend_alpha = score$blend_alpha,
        validation_blend_objective = score$blend_objective,
        validation_blend_objective_value = score$blend_objective_value,
        validation_blend_rmse = score$blend_rmse,
        validation_blend_cor = score$blend_cor,
        validation_blend_centered_rmse = score$blend_centered_rmse,
        validation_blend_centered_cor = score$blend_centered_cor,
        validation_blend_within_env_spearman = score$blend_within_env_spearman,
        validation_rmse_se = score$rmse_se,
        validation_centered_rmse_se = score$centered_rmse_se,
        validation_within_env_spearman_se = score$within_env_spearman_se,
        validation_uncertainty_scale = score$validation_uncertainty_scale %||% NA_real_,
        validation_uncertainty_raw_scale = score$validation_uncertainty_raw_scale %||% NA_real_,
        validation_uncertainty_n = score$validation_uncertainty_n %||% 0L,
        validation_uncertainty_n_env = score$validation_uncertainty_n_env %||% 0L,
        n_validation = score$n_validation,
        complexity = gp_tuning_candidate_complexity(candidate),
        elapsed_sec = proc.time()[["elapsed"]] - t0,
        failed = score$failed,
        error = score$error,
        stringsAsFactors = FALSE
      )
    }

    tuning_detail <- do.call(rbind, details)
    tuning_env_detail <- if (length(env_details)) do.call(rbind, env_details) else data.frame()
    env_kernel_cache_detail <- env_kernel_cache$detail
  }
  selected_info <- gp_tuning_select_candidate(tuning_detail, tuning_objective = tuning_objective)
  tuning_detail <- selected_info$tuning_detail
  if (is.null(selected_info$selected) || !nrow(selected_info$selected)) {
    return(list(
      status = "skipped",
      reason = selected_info$reason %||% "all tuning candidates failed",
      tuning_objective = tuning_objective,
      tuning_strategy = tuning_strategy,
      selected = NULL,
      tuning_detail = tuning_detail,
      tuning_env_detail = tuning_env_detail,
      env_kernel_cache_detail = env_kernel_cache_detail,
      tuning_sample = bounded$sample_info,
      validation_envs = blocked$validation_envs,
      calibration_envs = blocked$calibration_envs,
      validation_ids = blocked$validation_ids,
      batch_candidate_used = batch_candidate_used,
      batch_candidate_execution = batch_candidate_execution,
      uncertainty_calibration = list(status = "skipped", reason = selected_info$reason %||% "all tuning candidates failed")
    ))
  }

  selected <- selected_info$selected
  uncertainty_calibration <- if (isTRUE(calibrate_uncertainty)) {
    gp_tuning_fit_uncertainty_calibration(
      selected = selected,
      mean_cache = mean_cache,
      truth = truth,
      env_covariates = env_covariates,
      env_kernel_cache = env_kernel_cache,
      include_components = include_components_eff,
      fixed_effects = fixed_effects_eff,
      gp_backend = gp_backend,
      reaction_norm_feature_qc = reaction_norm_feature_qc,
      kenv_bandwidth = kenv_bandwidth,
      kenv_kernel_kwargs = kenv_kernel_kwargs,
      seed = seed,
      gp_se_hutchinson_probes = gp_se_hutchinson_probes,
      python_bin = python_bin,
      project_root = project_root
    )
  } else {
    list(status = "skipped", reason = "uncertainty calibration was not requested")
  }
  if (identical(uncertainty_calibration$status, "ok")) {
    selected$validation_uncertainty_scale <- as.numeric(uncertainty_calibration$scale)
    selected$validation_uncertainty_raw_scale <- as.numeric(uncertainty_calibration$raw_scale)
    selected$validation_uncertainty_n <- as.integer(uncertainty_calibration$n_validation_se)
    selected$validation_uncertainty_n_env <- as.integer(uncertainty_calibration$n_validation_env_se)
  }

  list(
    status = "ok",
    reason = NA_character_,
    tuning_objective = tuning_objective,
    tuning_strategy = tuning_strategy,
    tuning_objective_direction = if (isTRUE(selected_info$objective$maximize)) "maximize" else "minimize",
    tuning_objective_column = selected_info$objective$column,
    tuning_objective_best_value = selected_info$best_value,
    tuning_objective_tolerance = selected_info$tolerance,
    selected = selected,
    tuning_detail = tuning_detail,
    tuning_env_detail = tuning_env_detail,
    env_kernel_cache_detail = env_kernel_cache_detail,
    batch_candidate_used = batch_candidate_used,
    batch_candidate_execution = batch_candidate_execution,
    uncertainty_calibration = uncertainty_calibration,
    tuning_sample = bounded$sample_info,
    validation_envs = blocked$validation_envs,
    calibration_envs = blocked$calibration_envs,
    validation_ids = blocked$validation_ids,
    calibration_ids = blocked$calibration_ids,
    env_metadata = blocked$env_metadata
  )
}
