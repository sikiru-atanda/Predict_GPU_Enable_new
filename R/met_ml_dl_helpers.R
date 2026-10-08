gp_met_supported_models <- function() {
  c(
    "CatBoost", "LightGBM", "Xgboost", "RandomForest",
    "mlp", "ft_transformer", "saint", "tabnet", "moe"
  )
}

gp_is_met_ml_dl_model <- function(model) {
  !is.null(model) && any(as.character(model) %in% gp_met_supported_models())
}

gp_met_python_model_params <- function(model, params = list()) {
  model <- as.character(model %||% "")
  switch(
    model,
    CatBoost = list(
      catboost_iterations = params$catboost_iterations %||% 500,
      catboost_depth = params$catboost_depth %||% 6,
      catboost_learning_rate = params$catboost_learning_rate %||% 0.03,
      catboost_l2_leaf_reg = params$catboost_l2_leaf_reg %||% 3,
      catboost_thread_count = as.integer(params$catboost_thread_count %||% 1L)
    ),
    LightGBM = list(
      lightgbm_nrounds = params$lightgbm_nrounds %||% 100,
      lightgbm_learning_rate = params$lightgbm_learning_rate %||% 0.05,
      lightgbm_num_leaves = params$lightgbm_num_leaves %||% 31,
      lightgbm_feature_fraction = params$lightgbm_feature_fraction %||% 1.0,
      lightgbm_bagging_fraction = params$lightgbm_bagging_fraction %||% 1.0,
      lightgbm_min_data_in_leaf = params$lightgbm_min_data_in_leaf %||% 20,
      lightgbm_lambda_l1 = params$lightgbm_lambda_l1 %||% 0,
      lightgbm_lambda_l2 = params$lightgbm_lambda_l2 %||% 0,
      lightgbm_nthread = as.integer(params$lightgbm_nthread %||% 1L)
    ),
    Xgboost = list(
      nrounds = params$iteration %||% 100,
      eta = params$learning_rate %||% 0.01,
      max_depth = params$max_depth %||% 6,
      subsample = params$subsample %||% 0.7,
      xgb_gamma = params$xgb_gamma %||% 0.01,
      colsample_bytree = params$colsample_bytree %||% 0.7,
      min_child_weight = params$min_child_weight %||% 1,
      xgb_alpha = params$xgb_alpha %||% 0.001,
      xgb_lambda = params$xgb_lambda %||% 1.0,
      xgb_booster = params$xgb_booster %||% "gbtree",
      xgb_nthread = as.integer(params$xgb_nthread %||% 1L)
    ),
    RandomForest = list(
      ntree = params$ntree %||% 500,
      mtry = params$mtry %||% NULL,
      maxnodes = params$maxnodes %||% NULL,
      nodesize = params$nodesize %||% NULL,
      rf_n_jobs = as.integer(params$rf_n_jobs %||% 1L)
    ),
    stop("Unsupported MET ML/DL python model: ", model, call. = FALSE)
  )
}

gp_met_dl_model_params <- function(model, params = list()) {
  model <- as.character(model %||% "")
  fam <- params$response_family %||% "gaussian"
  switch(
    model,
    mlp = list(
      model_type = "mlp",
      response_family = fam,
      compile_model = params$compile_model %||% FALSE,
      deterministic = params$deterministic %||% TRUE,
      random_seed = params$random_seed %||% 123L,
      device = params$device %||% NULL,
      use_amp = params$use_amp %||% FALSE,
      batch_norm = params$batch_norm %||% FALSE,
      validation_split = params$validation_split %||% 0.2,
      epochs = params$epochs %||% 30L,
      batch_size = params$batch_size %||% 64L,
      mlp_neurons_per_layer = params$mlp_neurons_per_layer %||% as.integer(c(128L, 64L)),
      mlp_learning_rate = params$mlp_learning_rate %||% 1e-3,
      dropout = params$dropout %||% 0.2,
      dropout_rate = params$dropout_rate %||% (params$dropout %||% 0.2),
      l2_weight_decay = params$l2_weight_decay %||% 1e-4,
      l2_regularizer_dp = params$l2_regularizer_dp %||% (params$l2_weight_decay %||% 1e-4),
      final_attention = params$final_attention %||% FALSE,
      attention_across_multiple_layers = params$attention_across_multiple_layers %||% FALSE,
      optimizer_name = params$optimizer_name %||% "adam",
      max_grad_norm = params$max_grad_norm %||% NULL,
      heteroscedastic = params$heteroscedastic %||% FALSE
    ),
    ft_transformer = list(
      model_type = "ft_transformer",
      response_family = fam,
      compile_model = params$compile_model %||% FALSE,
      deterministic = params$deterministic %||% TRUE,
      random_seed = params$random_seed %||% 123L,
      device = params$device %||% NULL,
      use_amp = params$use_amp %||% FALSE,
      batch_norm = params$batch_norm %||% FALSE,
      validation_split = params$validation_split %||% 0.2,
      epochs = params$epochs %||% 30L,
      batch_size = params$batch_size %||% 64L,
      ft_d_model = params$ft_d_model %||% 64L,
      ft_heads = params$ft_heads %||% 4L,
      ft_layers = params$ft_layers %||% 2L,
      ft_ff_mult = params$ft_ff_mult %||% 2L,
      ft_dropout = params$ft_dropout %||% 0.1,
      ft_token_dropout = params$ft_token_dropout %||% 0.1,
      ft_use_cls = params$ft_use_cls %||% TRUE,
      dropout = params$dropout %||% 0.1,
      dropout_rate = params$dropout_rate %||% (params$dropout %||% 0.1),
      l2_weight_decay = params$l2_weight_decay %||% 1e-4,
      l2_regularizer_dp = params$l2_regularizer_dp %||% (params$l2_weight_decay %||% 1e-4),
      optimizer_name = params$optimizer_name %||% "adam",
      max_grad_norm = params$max_grad_norm %||% NULL,
      heteroscedastic = params$heteroscedastic %||% FALSE
    ),
    saint = list(
      model_type = "saint",
      response_family = fam,
      compile_model = params$compile_model %||% FALSE,
      deterministic = params$deterministic %||% TRUE,
      random_seed = params$random_seed %||% 123L,
      device = params$device %||% NULL,
      use_amp = params$use_amp %||% FALSE,
      batch_norm = params$batch_norm %||% FALSE,
      validation_split = params$validation_split %||% 0.2,
      epochs = params$epochs %||% 30L,
      batch_size = params$batch_size %||% 64L,
      saint_d_model = params$saint_d_model %||% 64L,
      saint_heads = params$saint_heads %||% 4L,
      saint_layers = params$saint_layers %||% 2L,
      saint_ff_mult = params$saint_ff_mult %||% 2L,
      saint_dropout = params$saint_dropout %||% 0.1,
      saint_token_dropout = params$saint_token_dropout %||% 0.1,
      saint_use_cls = params$saint_use_cls %||% TRUE,
      dropout = params$dropout %||% 0.1,
      dropout_rate = params$dropout_rate %||% (params$dropout %||% 0.1),
      l2_weight_decay = params$l2_weight_decay %||% 1e-4,
      l2_regularizer_dp = params$l2_regularizer_dp %||% (params$l2_weight_decay %||% 1e-4),
      optimizer_name = params$optimizer_name %||% "adam",
      max_grad_norm = params$max_grad_norm %||% NULL,
      heteroscedastic = params$heteroscedastic %||% FALSE
    ),
    tabnet = list(
      model_type = "tabnet",
      response_family = fam,
      compile_model = params$compile_model %||% FALSE,
      deterministic = params$deterministic %||% TRUE,
      random_seed = params$random_seed %||% 123L,
      device = params$device %||% NULL,
      use_amp = params$use_amp %||% FALSE,
      batch_norm = params$batch_norm %||% FALSE,
      validation_split = params$validation_split %||% 0.2,
      epochs = params$epochs %||% 30L,
      batch_size = params$batch_size %||% 64L,
      tabnet_steps = params$tabnet_steps %||% 3L,
      tabnet_feature_dim = params$tabnet_feature_dim %||% 16L,
      tabnet_output_dim = params$tabnet_output_dim %||% 16L,
      tabnet_gamma = params$tabnet_gamma %||% 1.3,
      tabnet_lambda_sparse = params$tabnet_lambda_sparse %||% 1e-4,
      dropout = params$dropout %||% 0.1,
      dropout_rate = params$dropout_rate %||% (params$dropout %||% 0.1),
      l2_weight_decay = params$l2_weight_decay %||% 1e-4,
      l2_regularizer_dp = params$l2_regularizer_dp %||% (params$l2_weight_decay %||% 1e-4),
      optimizer_name = params$optimizer_name %||% "adam",
      max_grad_norm = params$max_grad_norm %||% NULL,
      heteroscedastic = params$heteroscedastic %||% FALSE
    ),
    moe = list(
      model_type = "moe",
      response_family = fam,
      compile_model = params$compile_model %||% FALSE,
      deterministic = params$deterministic %||% TRUE,
      random_seed = params$random_seed %||% 123L,
      device = params$device %||% NULL,
      use_amp = params$use_amp %||% FALSE,
      batch_norm = params$batch_norm %||% FALSE,
      validation_split = params$validation_split %||% 0.2,
      epochs = params$epochs %||% 30L,
      batch_size = params$batch_size %||% 64L,
      moe_n_experts = params$moe_n_experts %||% 4L,
      moe_expert_hidden = params$moe_expert_hidden %||% as.integer(c(128L, 64L)),
      moe_gate_hidden = params$moe_gate_hidden %||% 128L,
      moe_temperature = params$moe_temperature %||% 1.0,
      moe_sparse_topk = params$moe_sparse_topk %||% NA,
      moe_entropy_reg = params$moe_entropy_reg %||% 0.0,
      dropout = params$dropout %||% 0.1,
      dropout_rate = params$dropout_rate %||% (params$dropout %||% 0.1),
      l2_weight_decay = params$l2_weight_decay %||% 1e-4,
      l2_regularizer_dp = params$l2_regularizer_dp %||% (params$l2_weight_decay %||% 1e-4),
      optimizer_name = params$optimizer_name %||% "adam",
      max_grad_norm = params$max_grad_norm %||% NULL,
      heteroscedastic = params$heteroscedastic %||% FALSE
    ),
    stop("Unsupported MET DL model: ", model, call. = FALSE)
  )
}

gp_met_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method, use = "complete.obs"))
}

gp_met_grouped_crossfit_calibration <- function(y,
                                                observed_idx,
                                                group_id,
                                                environment,
                                                n_targets,
                                                predictor_matrix = NULL,
                                                row_id = NULL,
                                                fit_predict,
                                                nfolds = 5L,
                                                seed = 123L) {
  y <- suppressWarnings(as.numeric(y))
  observed_idx <- as.integer(observed_idx)
  observed_idx <- observed_idx[
    is.finite(observed_idx) & observed_idx >= 1L & observed_idx <= length(y) &
      is.finite(y[observed_idx])
  ]
  group_id <- as.character(group_id)
  environment <- as.character(environment)
  if (is.null(row_id)) row_id <- as.character(seq_along(y))
  row_id <- as.character(row_id)
  predictor_payload_valid <- !is.null(predictor_matrix) &&
    nrow(as.matrix(predictor_matrix)) == length(y) &&
    length(row_id) == length(y) &&
    as.integer(n_targets) == length(y)
  if (!length(observed_idx) || length(group_id) != length(y) ||
      length(environment) != length(y)) {
    return(NULL)
  }
  observed_groups <- unique(group_id[observed_idx])
  observed_groups <- observed_groups[!is.na(observed_groups) & nzchar(observed_groups)]
  if (length(observed_groups) < 4L) {
    return(NULL)
  }
  nfolds <- suppressWarnings(as.integer(nfolds[[1L]]))
  if (!is.finite(nfolds) || nfolds < 2L) nfolds <- 2L
  nfolds <- min(nfolds, length(observed_groups))

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit(
    gp_clear_invalid_random_seed(old_seed, seed_existed),
    add = TRUE
  )
  gp_set_seed(seed)
  shuffled_groups <- sample(observed_groups, length(observed_groups), replace = FALSE)
  fold_groups <- split(
    shuffled_groups,
    rep(seq_len(nfolds), length.out = length(shuffled_groups))
  )

  cv_prediction <- rep(NA_real_, length(y))
  fold_id <- rep(NA_integer_, length(y))
  target_fold_predictions <- matrix(
    NA_real_, nrow = length(fold_groups), ncol = as.integer(n_targets)
  )
  for (i in seq_along(fold_groups)) {
    heldout_idx <- observed_idx[group_id[observed_idx] %in% fold_groups[[i]]]
    fit_idx <- setdiff(observed_idx, heldout_idx)
    if (!length(heldout_idx) || length(fit_idx) < 2L) next
    prediction_idx <- c(heldout_idx, seq_len(n_targets))
    fold_prediction <- suppressWarnings(as.numeric(
      fit_predict(fit_idx, prediction_idx, i)
    ))
    if (length(fold_prediction) != length(prediction_idx)) {
      stop("MET grouped cross-fit returned an unexpected prediction length.",
           call. = FALSE)
    }
    cv_prediction[heldout_idx] <- fold_prediction[seq_along(heldout_idx)]
    fold_id[heldout_idx] <- i
    target_fold_predictions[i, ] <- fold_prediction[
      length(heldout_idx) + seq_len(n_targets)
    ]
  }

  gp_ml_cv_rank_risk_calibration(
    predicted_cv = cv_prediction,
    observed_y = y,
    groups = environment,
    fold_id = fold_id,
    target_fold_predictions = target_fold_predictions,
    heldout_predictors = if (predictor_payload_valid) predictor_matrix else NULL,
    target_predictors = if (predictor_payload_valid) predictor_matrix else NULL,
    heldout_ids = if (predictor_payload_valid) row_id else NULL,
    target_ids = if (predictor_payload_valid) row_id else NULL
  )
}

gp_met_bootstrap_uncertainty_columns <- function(predicted_value,
                                                 bootstrap_matrix,
                                                 observed_y,
                                                 train_test_label,
                                                 heldout_calibration,
                                                 n_bootstrap,
                                                 confidence_level = 0.95) {
  predicted_value <- suppressWarnings(as.numeric(predicted_value))
  bootstrap_matrix <- as.matrix(bootstrap_matrix)
  if (ncol(bootstrap_matrix) != length(predicted_value) &&
      nrow(bootstrap_matrix) == length(predicted_value)) {
    bootstrap_matrix <- t(bootstrap_matrix)
  }
  if (ncol(bootstrap_matrix) != length(predicted_value) ||
      nrow(bootstrap_matrix) < 2L) {
    stop("MET bootstrap predictions do not align with the prediction rows.",
         call. = FALSE)
  }
  repeated_prediction_variance <- apply(
    bootstrap_matrix, 2L, stats::var, na.rm = TRUE
  )
  repeated_prediction_variance[
    !is.finite(repeated_prediction_variance) |
      repeated_prediction_variance < 0
  ] <- NA_real_
  calibrated <- gp_ml_calibrate_gaussian_uncertainty(
    predicted_value = predicted_value,
    pred_se = sqrt(repeated_prediction_variance),
    pred_variances = repeated_prediction_variance,
    observed_y = observed_y,
    train_test_label = train_test_label,
    boot_results = list(
      t = bootstrap_matrix,
      R = as.integer(n_bootstrap)
    ),
    heldout_calibration = heldout_calibration,
    confidence_level = confidence_level,
    require_heldout_calibration = TRUE
  )
  data.frame(
    Standard_error = calibrated$calibrated_standard_error,
    PEV = calibrated$calibrated_prediction_error_var,
    PEV_basis = calibrated$prediction_error_estimand,
    Prediction_uncertainty_source = calibrated$calibration_source,
    Prediction_interval_method = calibrated$prediction_interval_method,
    Prediction_interval_nominal_coverage =
      calibrated$prediction_interval_nominal_coverage,
    Prediction_interval_calibration_n =
      calibrated$prediction_interval_calibration_n,
    lower_bound = calibrated$calibrated_lower_bound,
    upper_bound = calibrated$calibrated_upper_bound,
    Repeated_prediction_variance = repeated_prediction_variance,
    stringsAsFactors = FALSE
  )
}

gp_met_attach_observed_predictions <- function(predicted_values,
                                               pheno_data,
                                               response,
                                               gen_name,
                                               heter_groups) {
  if (is.null(predicted_values)) {
    return(NULL)
  }
  if ("Observed_value" %in% names(predicted_values)) {
    out <- as.data.frame(predicted_values, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"Residual" %in% names(out) &&
        is.numeric(out$Observed_value) &&
        is.numeric(out$Predicted_value)) {
      out$Residual <- ifelse(
        is.na(out$Observed_value),
        NA_real_,
        out$Observed_value - out$Predicted_value
      )
    }
    return(out)
  }
  if (is.null(pheno_data)) {
    return(NULL)
  }
  join_cols <- c(gen_name, heter_groups)
  observed <- pheno_data[, c(join_cols, response), drop = FALSE]
  names(observed)[names(observed) == response] <- "Observed_value"
  out <- merge(
    predicted_values,
    observed,
    by = join_cols,
    all.x = TRUE,
    sort = FALSE
  )
  if (is.numeric(out$Observed_value) && is.numeric(out$Predicted_value)) {
    out$Residual <- ifelse(
      is.na(out$Observed_value),
      NA_real_,
      out$Observed_value - out$Predicted_value
    )
  } else {
    out$Residual <- NA_real_
  }
  out
}

gp_met_positive_variance <- function(x, min_var = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  val <- stats::var(x, na.rm = TRUE)
  if (is.finite(val) && val > min_var) val else NA_real_
}

gp_met_positive_mean_square <- function(x, min_var = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  val <- mean(x^2, na.rm = TRUE)
  if (is.finite(val) && val > min_var) val else NA_real_
}

gp_met_reliability_remarks <- function(reliability,
                                       high_reliability_thres = 0.9,
                                       low_reliability_thres = 0.5) {
  ifelse(
    is.na(reliability),
    NA_character_,
    ifelse(
      reliability >= high_reliability_thres,
      "Reliable",
      ifelse(reliability >= low_reliability_thres, "Acceptable", "Unreliable")
    )
  )
}

gp_met_empirical_reliability <- function(pev, reference_variance) {
  pev <- suppressWarnings(as.numeric(pev))
  reference_variance <- suppressWarnings(as.numeric(reference_variance))
  if (length(reference_variance) == 1L && length(pev) > 1L) {
    reference_variance <- rep(reference_variance, length(pev))
  }
  reliability <- rep(NA_real_, length(pev))
  fallback_used <- rep(FALSE, length(pev))
  valid <- is.finite(pev) & pev >= 0 &
    is.finite(reference_variance) & reference_variance > 0
  reliability[valid] <- pmax(0, pmin(1, 1 - pev[valid] / reference_variance[valid]))

  finite_rel <- reliability[valid & is.finite(reliability)]
  finite_pev <- pev[valid & is.finite(pev)]
  rel_is_flat <- length(unique(round(finite_rel, 12L))) <= 1L
  pev_varies <- length(unique(round(finite_pev, 12L))) > 1L
  single_saturated <- length(finite_rel) == 1L &&
    length(finite_pev) == 1L &&
    isTRUE(all.equal(finite_rel[[1L]], 0)) &&
    finite_pev[[1L]] > reference_variance[valid][[1L]]
  if ((isTRUE(rel_is_flat) && isTRUE(pev_varies)) || isTRUE(single_saturated)) {
    reliability[valid] <- reference_variance[valid] /
      (reference_variance[valid] + pev[valid])
    reliability[valid] <- pmax(0, pmin(1, reliability[valid]))
    # Mark which rows used the ratio-form fallback instead of the standard
    # r2 = 1 - PEV/sigma2_g, so callers can flag the mixed definition.
    fallback_used[valid] <- TRUE
  }

  attr(reliability, "ratio_fallback") <- fallback_used
  reliability
}

gp_met_add_gaussian_uncertainty <- function(pred_obs,
                                            gen_name,
                                            heter_groups,
                                            confidence_level = 0.95,
                                            high_reliability_thres = 0.9,
                                            low_reliability_thres = 0.5,
                                            ml_dl_estimand = FALSE) {
  if (is.null(pred_obs) || !is.data.frame(pred_obs) || !nrow(pred_obs) ||
      !"Predicted_value" %in% names(pred_obs)) {
    return(pred_obs)
  }
  out <- as.data.frame(pred_obs, stringsAsFactors = FALSE, check.names = FALSE)
  env_col <- if (length(heter_groups)) heter_groups[[1L]] else "Env"
  if (!env_col %in% names(out) && "Env" %in% names(out)) {
    env_col <- "Env"
  }
  if (!"Env" %in% names(out) && env_col %in% names(out)) {
    out[["Env"]] <- out[[env_col]]
  }

  pred <- suppressWarnings(as.numeric(out[["Predicted_value"]]))
  obs <- if ("Observed_value" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Observed_value"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  resid <- obs - pred
  train_label <- if ("Train_Test_Label" %in% names(out)) {
    as.character(out[["Train_Test_Label"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  train_row <- is.na(train_label) | train_label != "Test"
  if (!any(train_row, na.rm = TRUE)) {
    train_row <- rep(TRUE, nrow(out))
  }
  env <- if (env_col %in% names(out)) {
    as.character(out[[env_col]])
  } else {
    rep("All", nrow(out))
  }

  scale_values <- c(pred, obs)
  scale_values <- scale_values[is.finite(scale_values)]
  scale_ref <- if (length(scale_values)) stats::median(abs(scale_values), na.rm = TRUE) else 1
  if (!is.finite(scale_ref) || scale_ref <= 0) {
    scale_ref <- 1
  }
  min_var <- max((scale_ref * 1e-4)^2, .Machine$double.eps)

  global_pev <- gp_met_positive_mean_square(resid, min_var = min_var)
  if (!is.finite(global_pev)) {
    global_pev <- gp_met_positive_variance(resid, min_var = min_var)
  }
  if (!is.finite(global_pev)) {
    obs_var <- gp_met_positive_variance(obs, min_var = min_var)
    global_pev <- if (is.finite(obs_var)) obs_var * 0.25 else NA_real_
  }
  if (!is.finite(global_pev)) {
    pred_var <- gp_met_positive_variance(pred, min_var = min_var)
    global_pev <- if (is.finite(pred_var)) pred_var * 0.25 else min_var
  }
  global_pev <- max(global_pev, min_var)

  global_ref_var <- gp_met_positive_variance(obs[train_row], min_var = min_var)
  if (!is.finite(global_ref_var)) {
    global_ref_var <- gp_met_positive_variance(obs, min_var = min_var)
  }
  if (!is.finite(global_ref_var)) {
    global_ref_var <- gp_met_positive_variance(pred[train_row], min_var = min_var)
  }
  if (!is.finite(global_ref_var)) {
    global_ref_var <- max(global_pev, min_var)
  }

  env_levels <- unique(env[!is.na(env) & nzchar(env)])
  if (!length(env_levels)) {
    env_levels <- "All"
    env <- rep("All", nrow(out))
  }
  env_pev <- stats::setNames(rep(global_pev, length(env_levels)), env_levels)
  env_ref_var <- stats::setNames(rep(global_ref_var, length(env_levels)), env_levels)
  env_pred_var <- stats::setNames(rep(NA_real_, length(env_levels)), env_levels)
  env_pred_center <- stats::setNames(rep(NA_real_, length(env_levels)), env_levels)
  global_pred_var <- gp_met_positive_variance(pred, min_var = min_var)
  if (!is.finite(global_pred_var)) {
    global_pred_var <- max(global_ref_var, global_pev, min_var)
  }
  for (lev in env_levels) {
    idx <- env == lev
    lev_pev <- gp_met_positive_mean_square(resid[idx], min_var = min_var)
    if (!is.finite(lev_pev)) {
      lev_pev <- global_pev
    }
    lev_train <- idx & train_row
    if (!any(lev_train, na.rm = TRUE)) {
      lev_train <- idx
    }
    lev_ref <- gp_met_positive_variance(obs[lev_train], min_var = min_var)
    if (!is.finite(lev_ref)) {
      lev_ref <- global_ref_var
    }
    env_pev[[lev]] <- max(lev_pev, min_var)
    env_ref_var[[lev]] <- max(lev_ref, min_var)

    # The legacy marker-adjustment denominator used the complete final
    # prediction vector. Train_Test_Label identifies rows but does not subset
    # the fitted-prediction variance used for ML/DL reliability.
    lev_pred <- pred[idx]
    lev_pred <- lev_pred[is.finite(lev_pred)]
    env_pred_center[[lev]] <- if (length(lev_pred)) {
      stats::median(lev_pred, na.rm = TRUE)
    } else {
      NA_real_
    }
    lev_pred_var <- gp_met_positive_variance(lev_pred, min_var = min_var)
    env_pred_var[[lev]] <- if (is.finite(lev_pred_var)) lev_pred_var else global_pred_var
  }

  row_pev <- as.numeric(env_pev[env])
  row_pev[!is.finite(row_pev)] <- global_pev
  row_ref_var <- as.numeric(env_ref_var[env])
  row_ref_var[!is.finite(row_ref_var)] <- global_ref_var
  existing_ref_var <- if ("Reliability_reference_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Reliability_reference_variance"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  existing_reliability_variance <- if ("Reliability_variance_input" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Reliability_variance_input"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  existing_ref_ok <- is.finite(existing_ref_var) & existing_ref_var > 0
  row_ref_var[existing_ref_ok] <- existing_ref_var[existing_ref_ok]
  row_pred_center <- as.numeric(env_pred_center[env])
  row_pred_center[!is.finite(row_pred_center)] <- stats::median(pred[is.finite(pred)], na.rm = TRUE)
  row_pred_var <- as.numeric(env_pred_var[env])
  row_pred_var[!is.finite(row_pred_var) | row_pred_var <= 0] <- global_pred_var

  pred_distance_sq <- (pred - row_pred_center)^2
  dispersion_factor <- rep(1, nrow(out))
  dispersion_ok <- is.finite(pred_distance_sq) & is.finite(row_pred_var) & row_pred_var > 0
  dispersion_factor[dispersion_ok] <- 1 + pred_distance_sq[dispersion_ok] /
    (row_pred_var[dispersion_ok] + min_var)
  dispersion_factor <- pmax(0.5, pmin(3, dispersion_factor))

  empirical_row_pev <- row_pev * dispersion_factor
  resid_sq <- resid^2
  resid_ok <- is.finite(resid_sq)
  empirical_row_pev[resid_ok] <- 0.5 * empirical_row_pev[resid_ok] +
    0.5 * pmax(resid_sq[resid_ok], min_var)
  empirical_row_pev <- pmax(empirical_row_pev, min_var)

  existing_se <- if ("Standard_error" %in% names(out)) suppressWarnings(as.numeric(out[["Standard_error"]])) else rep(NA_real_, nrow(out))
  existing_pev <- if ("PEV" %in% names(out)) suppressWarnings(as.numeric(out[["PEV"]])) else rep(NA_real_, nrow(out))
  existing_prediction_var <- if ("Prediction_error_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Prediction_error_variance"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  missing_existing_pev <- !is.finite(existing_pev) | existing_pev < 0
  existing_pev[missing_existing_pev] <- existing_prediction_var[missing_existing_pev]
  missing_existing_pev <- !is.finite(existing_pev) | existing_pev < 0
  existing_pev[missing_existing_pev & is.finite(existing_se) & existing_se >= 0] <-
    existing_se[missing_existing_pev & is.finite(existing_se) & existing_se >= 0]^2
  validated_or_model_based <- is.finite(existing_pev) & existing_pev >= 0
  pev <- rep(NA_real_, nrow(out))
  pev[validated_or_model_based] <- pmax(existing_pev[validated_or_model_based], min_var)
  se <- rep(NA_real_, nrow(out))
  se[validated_or_model_based] <- sqrt(pmax(pev[validated_or_model_based], 0))
  keep_existing_se <- validated_or_model_based & is.finite(existing_se) & existing_se >= 0
  se[keep_existing_se] <- existing_se[keep_existing_se]

  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  lower <- if ("lower_bound" %in% names(out)) suppressWarnings(as.numeric(out[["lower_bound"]])) else rep(NA_real_, nrow(out))
  upper <- if ("upper_bound" %in% names(out)) suppressWarnings(as.numeric(out[["upper_bound"]])) else rep(NA_real_, nrow(out))
  lower[!validated_or_model_based] <- NA_real_
  upper[!validated_or_model_based] <- NA_real_
  bound_idx <- validated_or_model_based & (!is.finite(lower) | !is.finite(upper)) &
    is.finite(pred) & is.finite(se)
  lower[bound_idx] <- pred[bound_idx] - z * se[bound_idx]
  upper[bound_idx] <- pred[bound_idx] + z * se[bound_idx]

  uncertainty <- if ("Uncertainty" %in% names(out)) suppressWarnings(as.numeric(out[["Uncertainty"]])) else rep(NA_real_, nrow(out))
  uncertainty[!validated_or_model_based] <- NA_real_
  unc_idx <- !is.finite(uncertainty) & is.finite(lower) & is.finite(upper)
  uncertainty[unc_idx] <- upper[unc_idx] - lower[unc_idx]

  repeated_prediction_variance <- if ("Repeated_prediction_variance" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Repeated_prediction_variance"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  stability_variance_input <- if (isTRUE(ml_dl_estimand)) {
    empirical_row_pev
  } else {
    pev
  }
  use_repeated_prediction_variance <- isTRUE(ml_dl_estimand) &
    is.finite(repeated_prediction_variance) & repeated_prediction_variance >= 0
  stability_variance_input[use_repeated_prediction_variance] <-
    repeated_prediction_variance[use_repeated_prediction_variance]
  use_legacy_stability <- isTRUE(ml_dl_estimand) &
    (!is.finite(stability_variance_input) | stability_variance_input < 0)
  stability_variance_input[use_legacy_stability] <- empirical_row_pev[use_legacy_stability]
  stability_reference_variance <- if (isTRUE(ml_dl_estimand)) row_pred_var else row_ref_var
  stability_variance_source <- ifelse(
    use_repeated_prediction_variance,
    "repeated_prediction_se_squared_marker_adjustment",
    "legacy_in_sample_residual_dispersion_marker_adjustment"
  )
  stability <- gp_ml_gaussian_stability(
    prediction_error_var = stability_variance_input,
    reference_variance = stability_reference_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  predictive_precision <- gp_ml_predictive_reliability(
    prediction_error_var = pev,
    reference_variance = row_ref_var,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )
  existing_reliability <- if ("Reliability" %in% names(out)) {
    suppressWarnings(as.numeric(out[["Reliability"]]))
  } else {
    rep(NA_real_, nrow(out))
  }
  existing_reliability_basis <- if ("Reliability_basis" %in% names(out)) {
    as.character(out[["Reliability_basis"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  existing_pev_basis <- if ("PEV_basis" %in% names(out)) {
    as.character(out[["PEV_basis"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  existing_source <- if ("Prediction_uncertainty_source" %in% names(out)) {
    as.character(out[["Prediction_uncertainty_source"]])
  } else {
    rep(NA_character_, nrow(out))
  }
  ml_predictive_estimand <- isTRUE(ml_dl_estimand) | grepl(
    "not mixed-model PEV|cross-fitted|heldout",
    paste(
      ifelse(is.na(existing_pev_basis), "", existing_pev_basis),
      ifelse(is.na(existing_source), "", existing_source),
      ifelse(is.na(existing_reliability_basis), "", existing_reliability_basis)
    ),
    ignore.case = TRUE
  )
  marker_adjustment_surrogate <- is.finite(existing_reliability) & grepl(
    "marker-adjustment reliability surrogate",
    ifelse(is.na(existing_reliability_basis), "", existing_reliability_basis),
    ignore.case = TRUE
  )
  model_based_formula <- !ml_predictive_estimand & validated_or_model_based &
    existing_ref_ok
  model_based_reliability <- rep(NA_real_, nrow(out))
  model_based_reliability[model_based_formula] <- pmax(
    0,
    pmin(
      1,
      1 - pev[model_based_formula] / existing_ref_var[model_based_formula]
    )
  )
  preserve_reliability <- is.finite(existing_reliability) &
    !ml_predictive_estimand & !model_based_formula
  # Where the engine already reported this same value, keep its own basis text
  # (as PEV_basis and Prediction_uncertainty_source are kept); the generic
  # sentence is only right for rows whose reliability was actually recomputed.
  keep_engine_basis <- model_based_formula &
    is.finite(existing_reliability) & is.finite(model_based_reliability) &
    abs(existing_reliability - model_based_reliability) <= 1e-8 &
    !is.na(existing_reliability_basis) & nzchar(trimws(existing_reliability_basis))
  predictive_reliability <- ml_predictive_estimand

  out[["Standard_error"]] <- se
  out[["PEV"]] <- pev
  out[["Prediction_error_variance"]] <- pev
  out[["lower_bound"]] <- lower
  out[["upper_bound"]] <- upper
  out[["Uncertainty"]] <- uncertainty
  if (!"Uncertainty_remarks" %in% names(out)) {
    out[["Uncertainty_remarks"]] <- rep(NA_character_, nrow(out))
  }
  missing_unc_rem <- is.na(out[["Uncertainty_remarks"]]) | !nzchar(as.character(out[["Uncertainty_remarks"]]))
  out[["Uncertainty_remarks"]][missing_unc_rem & validated_or_model_based] <-
    "model_reported_prediction_uncertainty; inspect model assumptions"
  out[["Uncertainty_remarks"]][missing_unc_rem & !validated_or_model_based] <-
    "unavailable without held-out or model-based uncertainty"
  out[["Prediction_stability"]] <- ifelse(
    marker_adjustment_surrogate,
    existing_reliability,
    stability$stability
  )
  out[["Prediction_stability_remarks"]] <- ifelse(
    marker_adjustment_surrogate,
    ifelse(
      existing_reliability >= high_reliability_thres,
      "Stable",
      ifelse(existing_reliability <= low_reliability_thres, "Unstable", "Moderately Stable")
    ),
    stability$remarks
  )
  out[["Prediction_stability_reference_variance"]] <- ifelse(
    marker_adjustment_surrogate,
    existing_ref_var,
    stability_reference_variance
  )
  out[["Reliability_variance_input"]] <- ifelse(
    predictive_reliability,
    predictive_precision$variance_input,
    ifelse(model_based_formula, pev, ifelse(preserve_reliability, pev, NA_real_))
  )
  out[["Reliability_reference_variance"]] <- ifelse(
    predictive_reliability,
    predictive_precision$reference_variance,
    existing_ref_var
  )
  out[["Reliability"]] <- ifelse(
    predictive_reliability,
    predictive_precision$reliability,
    ifelse(
      model_based_formula,
      model_based_reliability,
      ifelse(preserve_reliability, existing_reliability, NA_real_)
    )
  )
  model_reliability_remarks <- if ("Reliability_remarks" %in% names(pred_obs)) {
    as.character(pred_obs[["Reliability_remarks"]])
  } else {
    rep("Model-based reliability", nrow(out))
  }
  out[["Reliability_remarks"]] <- ifelse(
    predictive_reliability,
    predictive_precision$remarks,
    ifelse(
      model_based_formula,
      gp_contract_reliability_remarks(model_based_reliability),
      ifelse(
        preserve_reliability,
        model_reliability_remarks,
        "Reliability unavailable until a model-specific estimand is supplied"
      )
    )
  )
  out[["Reliability_percentage"]] <- out[["Reliability"]] * 100
  out[["Reliability_basis"]] <- ifelse(
    predictive_reliability,
    predictive_precision$basis,
    ifelse(
      model_based_formula,
      ifelse(
        keep_engine_basis,
        existing_reliability_basis,
        paste(
          "1 - model-reported prediction-error variance / model-specific",
          "genetic variance; clipped to [0,1]"
        )
      ),
      ifelse(
        preserve_reliability,
        existing_reliability_basis,
        "unresolved; requires model-specific genetic variance and prediction-error estimand"
      )
    )
  )
  out[["Prediction_uncertainty_source"]] <- ifelse(
    validated_or_model_based,
    ifelse(is.na(existing_source) | !nzchar(existing_source),
           "model_reported_prediction_uncertainty", existing_source),
    "unavailable_no_heldout_or_model_based_uncertainty"
  )
  missing_pev_basis <- is.na(existing_pev_basis) | !nzchar(existing_pev_basis)
  existing_pev_basis[missing_pev_basis & validated_or_model_based] <-
    "model-reported prediction uncertainty; inspect model-specific assumptions"
  existing_pev_basis[missing_pev_basis & !validated_or_model_based] <-
    "unavailable without held-out or model-based uncertainty"
  out[["PEV_basis"]] <- existing_pev_basis
  out
}

gp_met_across_environment_prediction <- function(pred_obs,
                                                 gen_name,
                                                 heter_groups,
                                                 confidence_level = 0.95,
                                                 high_reliability_thres = 0.9,
                                                 low_reliability_thres = 0.5) {
  if (is.null(pred_obs) || !is.data.frame(pred_obs) || !nrow(pred_obs) ||
      !gen_name %in% names(pred_obs) || !"Predicted_value" %in% names(pred_obs)) {
    return(NULL)
  }
  env_col <- if (length(heter_groups)) heter_groups[[1L]] else "Env"
  if (!env_col %in% names(pred_obs) && "Env" %in% names(pred_obs)) {
    env_col <- "Env"
  }
  if (!env_col %in% names(pred_obs)) {
    return(NULL)
  }
  if (length(unique(as.character(pred_obs[[env_col]]))) <= 1L) {
    return(NULL)
  }

  split_rows <- split(seq_len(nrow(pred_obs)), as.character(pred_obs[[gen_name]]), drop = TRUE)
  rows <- lapply(names(split_rows), function(gid) {
    idx <- split_rows[[gid]]
    n_env <- length(unique(as.character(pred_obs[[env_col]][idx])))
    pred_vals <- suppressWarnings(as.numeric(pred_obs[["Predicted_value"]][idx]))
    obs_vals <- if ("Observed_value" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["Observed_value"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    pev_vals <- if ("PEV" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["PEV"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    ref_vals <- if ("Reliability_reference_variance" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["Reliability_reference_variance"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    stability_ref_vals <- if ("Prediction_stability_reference_variance" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["Prediction_stability_reference_variance"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    rel_vals <- if ("Reliability" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["Reliability"]][idx]))
    } else {
      rep(NA_real_, length(idx))
    }
    rel_variance_vals <- if ("Reliability_variance_input" %in% names(pred_obs)) {
      suppressWarnings(as.numeric(pred_obs[["Reliability_variance_input"]][idx]))
    } else {
      pev_vals
    }
    pred <- mean(pred_vals, na.rm = TRUE)
    if (!is.finite(pred)) pred <- NA_real_
    obs <- mean(obs_vals, na.rm = TRUE)
    if (!is.finite(obs)) obs <- NA_real_
    finite_pev <- pev_vals[is.finite(pev_vals) & pev_vals >= 0]
    pev <- if (length(finite_pev) && n_env > 0L) sum(finite_pev, na.rm = TRUE) / (n_env^2) else NA_real_
    se <- if (is.finite(pev)) sqrt(pmax(pev, 0)) else NA_real_
    ref_var <- mean(ref_vals[is.finite(ref_vals) & ref_vals > 0], na.rm = TRUE)
    if (!is.finite(ref_var)) ref_var <- NA_real_
    stability_ref_var <- mean(
      stability_ref_vals[is.finite(stability_ref_vals) & stability_ref_vals > 0],
      na.rm = TRUE
    )
    if (!is.finite(stability_ref_var)) stability_ref_var <- NA_real_
    finite_rel_variance <- rel_variance_vals[is.finite(rel_variance_vals) & rel_variance_vals >= 0]
    rel_variance <- if (length(finite_rel_variance) && n_env > 0L) {
      sum(finite_rel_variance, na.rm = TRUE) / (n_env^2)
    } else {
      NA_real_
    }
    stability <- gp_ml_gaussian_stability(
      rel_variance,
      stability_ref_var,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )
    predictive_precision <- gp_ml_predictive_reliability(
      prediction_error_var = pev,
      reference_variance = ref_var,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres
    )
    rel <- predictive_precision$reliability
    label <- if ("Train_Test_Label" %in% names(pred_obs) &&
                 any(as.character(pred_obs[["Train_Test_Label"]][idx]) == "Test", na.rm = TRUE)) {
      "Test"
    } else {
      "Train"
    }
    row <- data.frame(
      GID = gid,
      Predicted_value = pred,
      Train_Test_Label = label,
      Observed_value = obs,
      Standard_error = se,
      PEV = pev,
      PEV_basis = if (is.finite(pev)) {
        "across-environment predictive MSE under zero error covariance; not mixed-model PEV"
      } else {
        "unavailable without held-out or model-based uncertainty"
      },
      lower_bound = NA_real_,
      upper_bound = NA_real_,
      Uncertainty = NA_real_,
      Uncertainty_remarks = "across-environment interval unavailable without joint held-out calibration",
      Prediction_stability = stability$stability,
      Prediction_stability_remarks = stability$remarks,
      Prediction_stability_reference_variance = stability_ref_var,
      Reliability_variance_input = predictive_precision$variance_input,
      Reliability_reference_variance = predictive_precision$reference_variance,
      Reliability = rel,
      Reliability_remarks = predictive_precision$remarks,
      Reliability_percentage = rel * 100,
      Reliability_basis = predictive_precision$basis,
      Prediction_uncertainty_source = if (is.finite(pev)) {
        "across_environment_independence_approximation"
      } else {
        "unavailable_no_heldout_or_model_based_uncertainty"
      },
      Prediction_interval_method = "unavailable without joint held-out calibration",
      Prediction_interval_nominal_coverage = confidence_level,
      Prediction_interval_calibration_n = NA_real_,
      Prediction_target = "Across_environment_average",
      stringsAsFactors = FALSE
    )
    names(row)[names(row) == "GID"] <- gen_name
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

gp_met_majority_class <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1L]]
}

gp_met_normalize_probability_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) {
    return(x)
  }
  x <- pmax(0, pmin(1, x))
  total <- sum(x, na.rm = TRUE)
  if (is.finite(total) && total > 0) {
    x <- x / total
  }
  x
}

gp_met_classification_levels <- function(pred_obs, class_levels = NULL) {
  levs <- as.character(class_levels %||% character())
  if (length(levs)) {
    return(levs)
  }
  prob_cols <- gp_classification_probability_columns(pred_obs)
  if (length(prob_cols)) {
    levs <- sub("^(Probability|Prob)_", "", prob_cols)
    levs <- levs[!is.na(levs) & nzchar(levs)]
  }
  observed <- if ("Observed_value" %in% names(pred_obs)) {
    as.character(stats::na.omit(pred_obs[["Observed_value"]]))
  } else {
    character()
  }
  predicted <- unique(c(
    if ("Predicted_class" %in% names(pred_obs)) as.character(stats::na.omit(pred_obs[["Predicted_class"]])) else character(),
    if ("Predicted_value" %in% names(pred_obs)) as.character(stats::na.omit(pred_obs[["Predicted_value"]])) else character()
  ))
  levs <- unique(c(levs, observed, predicted))
  levs[!is.na(levs) & nzchar(levs)]
}

gp_met_across_classification_prediction <- function(pred_obs,
                                                    gen_name,
                                                    heter_groups,
                                                    class_levels = NULL,
                                                    response_family = "binary",
                                                    high_confidence = 0.9,
                                                    low_confidence = 0.5) {
  if (is.null(pred_obs) || !is.data.frame(pred_obs) || !nrow(pred_obs) ||
      !gen_name %in% names(pred_obs)) {
    return(NULL)
  }
  env_col <- if (length(heter_groups)) heter_groups[[1L]] else "Env"
  if (!env_col %in% names(pred_obs) && "Env" %in% names(pred_obs)) {
    env_col <- "Env"
  }
  if (!env_col %in% names(pred_obs)) {
    return(NULL)
  }
  if (length(unique(as.character(pred_obs[[env_col]]))) <= 1L) {
    return(NULL)
  }

  fam <- gp_resolve_response_family(response_family, y = pred_obs[["Observed_value"]] %||% NULL)
  if (!fam %in% c("binary", "multiclass")) {
    return(NULL)
  }
  levs <- gp_met_classification_levels(pred_obs, class_levels = class_levels)
  if (identical(fam, "binary") && length(levs) != 2L) {
    levs <- gp_ml_binary_levels(pred_obs[["Observed_value"]])
  }
  if (!length(levs)) {
    return(NULL)
  }

  prob_mat <- if (identical(fam, "binary")) {
    gp_met_binary_probability_matrix(pred_obs, class_levels = levs)
  } else {
    gp_met_multiclass_probability_matrix(pred_obs, class_levels = levs)
  }
  if (!is.matrix(prob_mat) || !nrow(prob_mat)) {
    return(NULL)
  }
  colnames(prob_mat) <- levs

  split_rows <- split(seq_len(nrow(pred_obs)), as.character(pred_obs[[gen_name]]), drop = TRUE)
  rows <- lapply(names(split_rows), function(gid) {
    idx <- split_rows[[gid]]
    mean_prob <- colMeans(prob_mat[idx, , drop = FALSE], na.rm = TRUE)
    mean_prob <- gp_met_normalize_probability_vector(mean_prob)
    if (all(is.na(mean_prob))) {
      mean_prob[] <- 1 / length(mean_prob)
    }
    cls <- gp_classification_prediction_summary(
      prob = matrix(mean_prob, nrow = 1L, dimnames = list(NULL, levs)),
      class_levels = levs,
      high_confidence = high_confidence,
      low_confidence = low_confidence
    )
    label <- if ("Train_Test_Label" %in% names(pred_obs) &&
                 any(as.character(pred_obs[["Train_Test_Label"]][idx]) == "Test", na.rm = TRUE)) {
      "Test"
    } else {
      "Train"
    }
    obs <- if ("Observed_value" %in% names(pred_obs)) {
      gp_met_majority_class(pred_obs[["Observed_value"]][idx])
    } else {
      NA_character_
    }
    row <- data.frame(
      GID = gid,
      Predicted_class = cls[["Predicted_class"]],
      Train_Test_Label = label,
      Observed_class = obs,
      Prediction_confidence = cls[["Prediction_confidence"]],
      Classification_uncertainty = cls[["Classification_uncertainty"]],
      Reliability = cls[["Prediction_confidence"]],
      Reliability_remarks = cls[["Prediction_confidence_remarks"]],
      stringsAsFactors = FALSE
    )
    prob_df <- as.data.frame(as.list(mean_prob), stringsAsFactors = FALSE)
    names(prob_df) <- vapply(levs, gp_classification_probability_name, character(1))
    out <- cbind(row, prob_df)
    names(out)[names(out) == "GID"] <- gen_name
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  gp_format_classification_prediction_table(
    out,
    gen_name = gen_name,
    include_env = FALSE
  )
}

gp_met_binary_probability_matrix <- function(pred_obs, class_levels = NULL) {
  prob_cols <- gp_classification_probability_columns(pred_obs)
  levs <- as.character(class_levels %||% gp_ml_binary_levels(pred_obs$Observed_value))

  if (length(prob_cols) >= 2L) {
    prob_mat <- as.matrix(pred_obs[, prob_cols, drop = FALSE])
    colnames(prob_mat) <- sub("^(Probability|Prob)_", "", prob_cols)
    missing_levels <- setdiff(levs, colnames(prob_mat))
    if (length(missing_levels)) {
      fill <- matrix(NA_real_, nrow = nrow(prob_mat), ncol = length(missing_levels))
      colnames(fill) <- missing_levels
      prob_mat <- cbind(prob_mat, fill)
    }
    return(prob_mat[, levs, drop = FALSE])
  }

  if ("Predicted_value" %in% names(pred_obs)) {
    p_pos <- as.numeric(pred_obs$Predicted_value)
    prob_mat <- cbind(1 - p_pos, p_pos)
    colnames(prob_mat) <- levs
    return(prob_mat)
  }

  matrix(NA_real_, nrow = nrow(pred_obs), ncol = length(levs), dimnames = list(NULL, levs))
}

gp_met_multiclass_probability_matrix <- function(pred_obs, class_levels = NULL) {
  prob_cols <- gp_classification_probability_columns(pred_obs)
  levs <- as.character(class_levels %||% sort(unique(c(
    as.character(stats::na.omit(pred_obs$Observed_value)),
    as.character(stats::na.omit(pred_obs$Predicted_class)),
    as.character(stats::na.omit(pred_obs$Predicted_value))
  ))))

  if (!length(prob_cols)) {
    return(matrix(NA_real_, nrow = nrow(pred_obs), ncol = length(levs), dimnames = list(NULL, levs)))
  }

  prob_mat <- as.matrix(pred_obs[, prob_cols, drop = FALSE])
  colnames(prob_mat) <- sub("^(Probability|Prob)_", "", prob_cols)
  missing_levels <- setdiff(levs, colnames(prob_mat))
  if (length(missing_levels)) {
    fill <- matrix(NA_real_, nrow = nrow(prob_mat), ncol = length(missing_levels))
    colnames(fill) <- missing_levels
    prob_mat <- cbind(prob_mat, fill)
  }
  prob_mat[, levs, drop = FALSE]
}

gp_met_summary_statistics <- function(predicted_values,
                                      pheno_data,
                                      response,
                                      gen_name,
                                      heter_groups,
                                      model_parameters = NULL,
                                      feature_summary = NULL,
                                      GS_model = NULL,
                                      response_family = "gaussian") {
  pred_obs <- gp_met_attach_observed_predictions(
    predicted_values = predicted_values,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  if (is.null(pred_obs) || !nrow(pred_obs)) {
    return(list(summary_statistics = NULL))
  }

  fam <- gp_resolve_response_family(response_family, y = pred_obs$Observed_value)
  across_env_predicted_values <- NULL
  if (identical(fam, "gaussian")) {
    pred_obs <- gp_met_add_gaussian_uncertainty(
      pred_obs = pred_obs,
      gen_name = gen_name,
      heter_groups = heter_groups,
      ml_dl_estimand = TRUE
    )
    across_env_predicted_values <- gp_met_across_environment_prediction(
      pred_obs = pred_obs,
      gen_name = gen_name,
      heter_groups = heter_groups
    )
  } else if (fam %in% c("binary", "multiclass")) {
    across_env_predicted_values <- gp_met_across_classification_prediction(
      pred_obs = pred_obs,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response_family = fam
    )
  }
  predictive_matrices <- if (identical(fam, "gaussian")) {
    gp_ml_predictive_matrix_bundle(
      pred = pred_obs,
      id_col = gen_name,
      dimension_col = heter_groups,
      error_pred = pred_obs
    )
  } else {
    list()
  }
  observed_rows <- pred_obs[!is.na(pred_obs$Observed_value), , drop = FALSE]
  n_train <- sum(pred_obs$Train_Test_Label == "Train", na.rm = TRUE)
  n_test <- sum(pred_obs$Train_Test_Label == "Test", na.rm = TRUE)
  n_env <- length(unique(pred_obs[[heter_groups]]))
  n_gid <- length(unique(pred_obs[[gen_name]]))
  n_blocks <- if (is.null(feature_summary)) 0L else nrow(feature_summary)
  n_pcs <- if (is.null(feature_summary) || !nrow(feature_summary)) 0L else sum(feature_summary$n_pcs)

  if (identical(fam, "binary")) {
    class_levels <- gp_ml_binary_levels(observed_rows$Observed_value)
    pooled_prob <- gp_met_binary_probability_matrix(observed_rows, class_levels = class_levels)
    pooled_metrics <- data.frame(
      metric = c("accuracy", "balanced_accuracy", "log_loss", "brier_score", "ece", "n_observed_rows"),
      summary = c(
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob[, 2], "accuracy", response_family = "binary") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob[, 2], "balanced_accuracy", response_family = "binary") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob[, 2], "log_loss", response_family = "binary") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob[, 2], "brier_score", response_family = "binary") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob[, 2], "ece", response_family = "binary") else NA_real_,
        nrow(observed_rows)
      ),
      stringsAsFactors = FALSE
    )

    env_metrics <- if (nrow(observed_rows)) {
      dplyr::bind_rows(lapply(split(observed_rows, observed_rows[[heter_groups]]), function(df) {
        prob_mat <- gp_met_binary_probability_matrix(df, class_levels = class_levels)
        data.frame(
          Env = unique(df[[heter_groups]])[1],
          n_observed_rows = nrow(df),
          accuracy = evaluation_metrics(df$Observed_value, prob_mat[, 2], "accuracy", response_family = "binary"),
          balanced_accuracy = evaluation_metrics(df$Observed_value, prob_mat[, 2], "balanced_accuracy", response_family = "binary"),
          log_loss = evaluation_metrics(df$Observed_value, prob_mat[, 2], "log_loss", response_family = "binary"),
          brier_score = evaluation_metrics(df$Observed_value, prob_mat[, 2], "brier_score", response_family = "binary"),
          ece = evaluation_metrics(df$Observed_value, prob_mat[, 2], "ece", response_family = "binary"),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame(
        Env = character(0),
        n_observed_rows = integer(0),
        accuracy = numeric(0),
        balanced_accuracy = numeric(0),
        log_loss = numeric(0),
        brier_score = numeric(0),
        ece = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  } else if (identical(fam, "multiclass")) {
    class_levels <- sort(unique(as.character(stats::na.omit(observed_rows$Observed_value))))
    pooled_prob <- gp_met_multiclass_probability_matrix(observed_rows, class_levels = class_levels)
    pooled_metrics <- data.frame(
      metric = c("accuracy", "balanced_accuracy", "macro_precision", "macro_recall", "macro_f1", "log_loss", "ece", "n_observed_rows"),
      summary = c(
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, observed_rows$Predicted_value, "accuracy", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, observed_rows$Predicted_value, "balanced_accuracy", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, observed_rows$Predicted_value, "macro_precision", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, observed_rows$Predicted_value, "macro_recall", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, observed_rows$Predicted_value, "macro_f1", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob, "log_loss", response_family = "multiclass") else NA_real_,
        if (nrow(observed_rows)) evaluation_metrics(observed_rows$Observed_value, pooled_prob, "ece", response_family = "multiclass") else NA_real_,
        nrow(observed_rows)
      ),
      stringsAsFactors = FALSE
    )

    env_metrics <- if (nrow(observed_rows)) {
      dplyr::bind_rows(lapply(split(observed_rows, observed_rows[[heter_groups]]), function(df) {
        prob_mat <- gp_met_multiclass_probability_matrix(df, class_levels = class_levels)
        data.frame(
          Env = unique(df[[heter_groups]])[1],
          n_observed_rows = nrow(df),
          accuracy = evaluation_metrics(df$Observed_value, df$Predicted_value, "accuracy", response_family = "multiclass"),
          balanced_accuracy = evaluation_metrics(df$Observed_value, df$Predicted_value, "balanced_accuracy", response_family = "multiclass"),
          macro_precision = evaluation_metrics(df$Observed_value, df$Predicted_value, "macro_precision", response_family = "multiclass"),
          macro_recall = evaluation_metrics(df$Observed_value, df$Predicted_value, "macro_recall", response_family = "multiclass"),
          macro_f1 = evaluation_metrics(df$Observed_value, df$Predicted_value, "macro_f1", response_family = "multiclass"),
          log_loss = evaluation_metrics(df$Observed_value, prob_mat, "log_loss", response_family = "multiclass"),
          ece = evaluation_metrics(df$Observed_value, prob_mat, "ece", response_family = "multiclass"),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame(
        Env = character(0),
        n_observed_rows = integer(0),
        accuracy = numeric(0),
        balanced_accuracy = numeric(0),
        macro_precision = numeric(0),
        macro_recall = numeric(0),
        macro_f1 = numeric(0),
        log_loss = numeric(0),
        ece = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  } else {
    pooled_metrics <- data.frame(
      metric = c(
        "root_mean_squared_error",
        "mean_absolute_error",
        "pearsons_correlation",
        "kendalls_tau",
        "n_observed_rows"
      ),
      summary = c(
        if (nrow(observed_rows)) sqrt(mean((observed_rows$Observed_value - observed_rows$Predicted_value)^2, na.rm = TRUE)) else NA_real_,
        if (nrow(observed_rows)) mean(abs(observed_rows$Observed_value - observed_rows$Predicted_value), na.rm = TRUE) else NA_real_,
        gp_met_safe_cor(observed_rows$Observed_value, observed_rows$Predicted_value, method = "pearson"),
        gp_met_safe_cor(observed_rows$Observed_value, observed_rows$Predicted_value, method = "kendall"),
        nrow(observed_rows)
      ),
      stringsAsFactors = FALSE
    )

    env_metrics <- if (nrow(observed_rows)) {
      dplyr::bind_rows(lapply(split(observed_rows, observed_rows[[heter_groups]]), function(df) {
        data.frame(
          Env = unique(df[[heter_groups]])[1],
          n_observed_rows = nrow(df),
          root_mean_squared_error = sqrt(mean((df$Observed_value - df$Predicted_value)^2, na.rm = TRUE)),
          mean_absolute_error = mean(abs(df$Observed_value - df$Predicted_value), na.rm = TRUE),
          pearsons_correlation = gp_met_safe_cor(df$Observed_value, df$Predicted_value, method = "pearson"),
          kendalls_tau = gp_met_safe_cor(df$Observed_value, df$Predicted_value, method = "kendall"),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame(
        Env = character(0),
        n_observed_rows = integer(0),
        root_mean_squared_error = numeric(0),
        mean_absolute_error = numeric(0),
        pearsons_correlation = numeric(0),
        kendalls_tau = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }

  counts_by_env <- pred_obs |>
    dplyr::group_by(.data[[heter_groups]]) |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_observed_rows = sum(!is.na(.data$Observed_value)),
      n_train = sum(.data$Train_Test_Label == "Train", na.rm = TRUE),
      n_test = sum(.data$Train_Test_Label == "Test", na.rm = TRUE),
      .groups = "drop"
    )
  names(counts_by_env)[1] <- "Env"

  summary_statistics <- data.frame(
    stat = c(
      "Min",
      "Max",
      "Phenotype_Variance",
      "Number_TrainingSet",
      "Number_TestingSet",
      "Number_Environments",
      "Number_Genotypes",
      "Number_KernelFeatureBlocks",
      "Number_KernelPCs",
      "MET_Output_Scope",
      "Response_Family",
      "GS_model"
    ),
    summary = c(
      if (identical(fam, "gaussian") && nrow(observed_rows)) round(min(observed_rows$Observed_value, na.rm = TRUE), 3) else NA,
      if (identical(fam, "gaussian") && nrow(observed_rows)) round(max(observed_rows$Observed_value, na.rm = TRUE), 3) else NA,
      if (identical(fam, "gaussian") && nrow(observed_rows)) round(stats::var(observed_rows$Observed_value, na.rm = TRUE), 3) else NA,
      n_train,
      n_test,
      n_env,
      n_gid,
      n_blocks,
      n_pcs,
      if (!is.null(across_env_predicted_values)) {
        "Per-environment and across-environment prediction"
      } else {
        "Per-environment prediction only"
      },
      fam,
      GS_model %||% NA_character_
    ),
    stringsAsFactors = FALSE
  )

  if (!is.null(model_parameters) && nrow(model_parameters)) {
    summary_statistics <- rbind(
      summary_statistics,
      data.frame(
        stat = paste0("model_parameter__", model_parameters$stat),
        summary = as.character(model_parameters$summary),
        stringsAsFactors = FALSE
      )
    )
  }

  c(list(
    summary_statistics = summary_statistics,
    MET_pooled_metrics = pooled_metrics,
    MET_metrics_by_environment = env_metrics,
    MET_prediction_counts_by_environment = counts_by_env,
    MET_prediction_table = pred_obs,
    across_environment_predicted_values = across_env_predicted_values,
    Total_Predicted_value = across_env_predicted_values,
    predictive_covariance_basis = if (identical(fam, "gaussian")) {
      "descriptive predictions and in-sample prediction errors; not genetic or residual covariance"
    } else {
      NULL
    }
  ), predictive_matrices)
}

gp_met_true_prediction_plot <- function(predicted_values,
                                        pheno_data,
                                        response,
                                        gen_name,
                                        heter_groups,
                                        model_label = NULL,
                                        response_family = "gaussian") {
  pred_obs <- gp_met_attach_observed_predictions(
    predicted_values = predicted_values,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups
  )
  if (is.null(pred_obs) || !nrow(pred_obs)) {
    return(NULL)
  }

  fam <- gp_resolve_response_family(response_family, y = pred_obs$Observed_value)
  if (identical(fam, "gaussian")) {
    pred_obs <- gp_met_add_gaussian_uncertainty(
      pred_obs = pred_obs,
      gen_name = gen_name,
      heter_groups = heter_groups,
      ml_dl_estimand = TRUE
    )
    if (exists("gp_format_gaussian_prediction_table", mode = "function") &&
        exists("gp_gaussian_diagnostic_plots", mode = "function")) {
      public_pred <- gp_format_gaussian_prediction_table(
        pred_obs,
        gen_name = gen_name,
        include_env = TRUE
      )
      gaussian_plots <- gp_gaussian_diagnostic_plots(public_pred)
      if (is.list(gaussian_plots) &&
          any(vapply(gaussian_plots, function(x) !is.null(x), logical(1L)))) {
        return(gaussian_plots)
      }
    }
  }
  observed_rows <- pred_obs[!is.na(pred_obs$Observed_value), , drop = FALSE]
  title_suffix <- paste(response, model_label %||% "", sep = " - ")

  if (identical(fam, "binary")) {
    class_levels <- gp_ml_binary_levels(observed_rows$Observed_value)
    prob_mat <- gp_met_binary_probability_matrix(observed_rows, class_levels = class_levels)
    p_pos <- prob_mat[, 2]
    calib_df <- data.frame(
      Env = observed_rows[[heter_groups]],
      observed = as.character(observed_rows$Observed_value) == class_levels[2],
      prob = p_pos,
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(bin = cut(.data$prob, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) |>
      dplyr::group_by(.data$Env, .data$bin) |>
      dplyr::summarise(
        mean_prob = mean(.data$prob, na.rm = TRUE),
        observed_rate = mean(.data$observed, na.rm = TRUE),
        .groups = "drop"
      )

    p1 <- ggplot2::ggplot(calib_df, ggplot2::aes(x = mean_prob, y = observed_rate, color = Env)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("MET Binary Calibration", title_suffix),
        x = paste("Mean predicted P(", class_levels[2], ")", sep = ""),
        y = paste("Observed rate of", class_levels[2]),
        color = heter_groups
      ) +
      ggplot2::theme_minimal()

    conf_df <- data.frame(
      Env = observed_rows[[heter_groups]],
      Confidence = pmax(p_pos, 1 - p_pos),
      Correct = ifelse(ifelse(p_pos >= 0.5, class_levels[2], class_levels[1]) == as.character(observed_rows$Observed_value), "Correct", "Incorrect"),
      stringsAsFactors = FALSE
    )
    p2 <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Env, y = Confidence, fill = Correct)) +
      ggplot2::geom_boxplot(alpha = 0.75) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("MET Binary Prediction Confidence", title_suffix),
        x = heter_groups,
        y = "Prediction confidence",
        fill = "Outcome"
      ) +
      ggplot2::theme_minimal()
  } else if (identical(fam, "multiclass")) {
    class_levels <- sort(unique(as.character(stats::na.omit(observed_rows$Observed_value))))
    cm_df <- as.data.frame(
      table(
        Observed = factor(observed_rows$Observed_value, levels = class_levels),
        Predicted = factor(observed_rows$Predicted_value, levels = class_levels)
      ),
      stringsAsFactors = FALSE
    )
    p1 <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Predicted, y = Observed, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = Freq), color = "white", size = 4) +
      ggplot2::labs(
        title = paste("MET Multiclass Confusion Matrix", title_suffix),
        x = "Predicted class",
        y = "Observed class",
        fill = "Count"
      ) +
      ggplot2::theme_minimal()

    conf_df <- data.frame(
      Env = observed_rows[[heter_groups]],
      Confidence = observed_rows$Prediction_confidence,
      Predicted_class = observed_rows$Predicted_value,
      stringsAsFactors = FALSE
    )
    p2 <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Env, y = Confidence, fill = Predicted_class)) +
      ggplot2::geom_boxplot(alpha = 0.75) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("MET Multiclass Prediction Confidence", title_suffix),
        x = heter_groups,
        y = "Prediction confidence",
        fill = "Predicted class"
      ) +
      ggplot2::theme_minimal()
    return(gridExtra::grid.arrange(p1, p2, ncol = 2))
  } else if (nrow(observed_rows) >= 2L) {
    p1 <- ggplot2::ggplot(
      observed_rows,
      ggplot2::aes(x = Observed_value, y = Predicted_value, color = .data[[heter_groups]])
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(
        title = paste("MET Observed vs Predicted", title_suffix),
        x = "Observed value",
        y = "Predicted value",
        color = heter_groups
      ) +
      ggplot2::theme_minimal()
  } else {
    p1 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "Not enough observed rows for scatter plot") +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste("MET Observed vs Predicted", title_suffix))
  }

  p2 <- ggplot2::ggplot(
    pred_obs,
    ggplot2::aes(x = .data[[heter_groups]], y = Predicted_value, fill = Train_Test_Label)
  ) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.alpha = 0.5) +
    ggplot2::labs(
      title = paste("MET Predicted Value by Environment", title_suffix),
      x = heter_groups,
      y = "Predicted value",
      fill = "Subset"
    ) +
    ggplot2::theme_minimal()

  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

gp_kernel_to_eigenfeatures <- function(kernel,
                                       prefix,
                                       var_explained = 0.95,
                                       min_ev = 1e-8,
                                       max_pcs = NULL) {
  if (is.null(kernel)) {
    return(NULL)
  }

  kernel <- as.matrix(kernel)
  if (is.null(rownames(kernel)) || is.null(colnames(kernel))) {
    stop("Kernel must have row and column names.", call. = FALSE)
  }
  if (!identical(sort(rownames(kernel)), sort(colnames(kernel)))) {
    stop("Kernel row/column names must represent the same genotype set.", call. = FALSE)
  }

  common_ids <- intersect(rownames(kernel), colnames(kernel))
  kernel <- kernel[common_ids, common_ids, drop = FALSE]
  if (!isTRUE(all.equal(kernel, t(kernel), tolerance = 1e-8))) {
    kernel <- (kernel + t(kernel)) / 2
  }

  evd <- eigen(kernel, symmetric = TRUE)
  keep <- which(is.finite(evd$values) & evd$values > min_ev)
  if (!length(keep)) {
    return(data.frame(row.names = common_ids))
  }

  vals <- evd$values[keep]
  vecs <- evd$vectors[, keep, drop = FALSE]
  total <- sum(vals)
  if (!is.finite(total) || total <= 0) {
    npc <- length(vals)
  } else {
    cumprop <- cumsum(vals) / total
    npc <- which(cumprop >= var_explained)[1]
    if (is.na(npc) || npc < 1L) npc <- length(vals)
  }
  if (!is.null(max_pcs)) {
    npc <- min(npc, as.integer(max_pcs))
  }
  npc <- max(1L, min(npc, length(vals)))

  vals <- vals[seq_len(npc)]
  vecs <- vecs[, seq_len(npc), drop = FALSE]
  feat <- sweep(vecs, 2, sqrt(pmax(vals, 0)), `*`)
  feat <- as.data.frame(feat, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(feat) <- paste0(prefix, "_PC", seq_len(ncol(feat)))
  rownames(feat) <- common_ids
  attr(feat, "eigen_values") <- vals
  attr(feat, "var_explained") <- if (total > 0) sum(vals) / total else NA_real_
  feat
}

gp_prepare_met_kernel_features <- function(gmatrix = NULL,
                                           omic1_kernel = NULL,
                                           omic2_kernel = NULL,
                                           omic3_kernel = NULL,
                                           kernel_list = NULL,
                                           omics_kernel_label = list(
                                             omic1_kernel = NULL,
                                             omic2_kernel = NULL,
                                             omic3_kernel = NULL
                                           ),
                                           var_explained = 0.95,
                                           min_ev = 1e-8,
                                           max_pcs = NULL) {
  kernels <- list(
    GRM = gmatrix,
    omic1 = omic1_kernel,
    omic2 = omic2_kernel,
    omic3 = omic3_kernel
  )
  if (gp_is_kernel_list(kernel_list)) {
    kernel_names <- names(kernel_list)
    if (is.null(kernel_names) || any(!nzchar(kernel_names))) {
      kernel_names <- paste0("kernel", seq_along(kernel_list))
    }
    canonical_keys <- c(
      gmatrix_model_ready = "GRM",
      omic1_kernel_model_ready = "omic1",
      omic2_kernel_model_ready = "omic2",
      omic3_kernel_model_ready = "omic3"
    )
    for (i in seq_along(kernel_list)) {
      source_name <- gp_kernel_source_name(kernel_list[[i]], kernel_names[[i]])
      canonical_name <- unname(canonical_keys[kernel_names[[i]]])
      if (length(canonical_name) && !is.na(canonical_name) && nzchar(canonical_name)) {
        # These prepared entries mirror the explicit primary arguments.
        if (is.null(kernels[[canonical_name]])) {
          kernels[[canonical_name]] <- kernel_list[[i]]
        }
      } else {
        # A source label is not an identity: append first, then disambiguate.
        kernels <- c(kernels, stats::setNames(
          list(kernel_list[[i]]), gp_sanitize_kernel_name(source_name)
        ))
      }
    }
  }
  kernels <- kernels[!vapply(kernels, is.null, logical(1L))]
  names(kernels) <- make.unique(names(kernels), sep = "_")

  omic_labels <- c(
    omic1 = omics_kernel_label[["omic1_kernel"]] %||% "omic1",
    omic2 = omics_kernel_label[["omic2_kernel"]] %||% "omic2",
    omic3 = omics_kernel_label[["omic3_kernel"]] %||% "omic3"
  )

  prefixes <- make.unique(vapply(names(kernels), function(nm) {
    if (identical(nm, "GRM")) {
      "GRM"
    } else if (!nm %in% names(omic_labels)) {
      make.names(nm, unique = TRUE)
    } else {
      make.names(omic_labels[[nm]], unique = TRUE)
    }
  }, character(1L)), sep = "_")

  blocks <- list()
  for (i in seq_along(kernels)) {
    nm <- names(kernels)[[i]]
    K <- kernels[[i]]
    block <- gp_kernel_to_eigenfeatures(
      kernel = K,
      prefix = prefixes[[i]],
      var_explained = var_explained,
      min_ev = min_ev,
      max_pcs = max_pcs
    )
    if (!is.null(block) && ncol(block) > 0) {
      blocks[[nm]] <- block
    }
  }

  if (!length(blocks)) {
    return(list(
      feature_table = NULL,
      feature_blocks = list(),
      feature_summary = data.frame()
    ))
  }

  common_ids <- Reduce(intersect, lapply(blocks, rownames))
  blocks <- lapply(blocks, function(x) x[common_ids, , drop = FALSE])
  feature_table <- do.call(cbind, unname(blocks))
  feature_summary <- do.call(
    rbind,
    lapply(names(blocks), function(nm) {
      blk <- blocks[[nm]]
      data.frame(
        source = nm,
        n_pcs = ncol(blk),
        stringsAsFactors = FALSE
      )
    })
  )

  list(
    feature_table = feature_table,
    feature_blocks = blocks,
    feature_summary = feature_summary
  )
}

gp_build_met_long_data <- function(pheno_data,
                                   response,
                                   gen_name,
                                   heter_groups,
                                   kernel_feature_table,
                                   env_covariates = NULL) {
  if (is.null(pheno_data) || !is.data.frame(pheno_data)) {
    stop("pheno_data must be a data.frame.", call. = FALSE)
  }
  if (is.null(kernel_feature_table) || !nrow(kernel_feature_table)) {
    stop("kernel_feature_table must contain at least one genotype.", call. = FALSE)
  }
  req <- c(response, gen_name, heter_groups)
  miss <- setdiff(req, names(pheno_data))
  if (length(miss)) {
    stop("Missing required columns in pheno_data: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  dat <- pheno_data[, req, drop = FALSE]
  dat[[gen_name]] <- as.character(dat[[gen_name]])
  dat[[heter_groups]] <- as.character(dat[[heter_groups]])
  dat <- dat[dat[[gen_name]] %in% rownames(kernel_feature_table), , drop = FALSE]

  feat <- kernel_feature_table[match(dat[[gen_name]], rownames(kernel_feature_table)), , drop = FALSE]
  rownames(feat) <- NULL
  long_dat <- cbind(dat, feat)

  if (!is.null(env_covariates)) {
    env_covariates <- as.data.frame(env_covariates, stringsAsFactors = FALSE)
    if (!heter_groups %in% names(env_covariates)) {
      stop("env_covariates must contain the heter_groups column.", call. = FALSE)
    }
    env_covariates[[heter_groups]] <- as.character(env_covariates[[heter_groups]])
    long_dat <- merge(
      long_dat,
      env_covariates,
      by = heter_groups,
      all.x = TRUE,
      sort = FALSE
    )
    long_dat <- long_dat[match(paste(dat[[gen_name]], dat[[heter_groups]], sep = "\r"),
                               paste(long_dat[[gen_name]], long_dat[[heter_groups]], sep = "\r")), , drop = FALSE]
  }

  long_dat$.met_row_id <- paste(long_dat[[gen_name]], long_dat[[heter_groups]], sep = "__")
  long_dat
}

gp_met_long_test_labels <- function(long_data,
                                    gen_name,
                                    heter_groups,
                                    test_rows = NULL,
                                    test_gids = NULL,
                                    test_envs = NULL) {
  n <- nrow(long_data)
  label <- rep("Train", n)
  if (!is.null(test_rows) && length(test_rows)) {
    label[test_rows] <- "Test"
  }
  if (!is.null(test_gids) && length(test_gids)) {
    label[as.character(long_data[[gen_name]]) %in% as.character(test_gids)] <- "Test"
  }
  if (!is.null(test_envs) && length(test_envs)) {
    label[as.character(long_data[[heter_groups]]) %in% as.character(test_envs)] <- "Test"
  }
  label
}

gp_met_cv_assignments <- function(pheno_data,
                                  gen_name,
                                  heter_groups,
                                  cross_validation_meth,
                                  nfolds = 5L,
                                  random_state = NULL,
                                  replication = 1L,
                                  message = FALSE) {
  cv_token <- normalize_cv_token(cross_validation_meth)
  CVi <- switch(
    cv_token,
    "cv0" = 0L,
    "repeated_cv0" = 0L,
    "cv1" = 1L,
    "repeated_cv1" = 1L,
    "cv2" = 2L,
    "repeated_cv2" = 2L,
    NULL
  )
  if (is.null(CVi)) {
    stop("Multi-environment cross-validation supports only CV0/CV1/CV2 scenario names.", call. = FALSE)
  }

  CV0_CV1_CV2_for_multi_environment(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    CV = CVi,
    nfolds = nfolds,
    random_state = random_state,
    replication = replication,
    message = message
  )
}

gp_met_python_cv_predict <- function(pheno_data,
                                     response,
                                     gen_name,
                                     heter_groups,
                                     response_family = "gaussian",
                                     model_type,
                                     gmatrix = NULL,
                                     omic1_kernel = NULL,
                                     omic2_kernel = NULL,
                                     omic3_kernel = NULL,
                                     kernel_list = NULL,
                                     omics_kernel_label = list(
                                       omic1_kernel = NULL,
                                       omic2_kernel = NULL,
                                       omic3_kernel = NULL
                                     ),
                                     tst,
                                     model_params = list(),
                                     met_kernel_var_explained = 0.95,
                                     met_kernel_min_ev = 1e-8,
                                     met_kernel_max_pcs = NULL) {
  feat_res <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )
  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feat_res$feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  fam <- gp_resolve_response_family(response_family, y = long_dat[[response]])
  y_all <- if (identical(fam, "gaussian")) as.numeric(long_dat[[response]]) else long_dat[[response]]
  train_idx <- setdiff(which(!is.na(y_all)), tst)
  if (!length(train_idx)) {
    stop("No observed training rows available for MET ", model_type, " CV.", call. = FALSE)
  }
  gp_py_ml_fit_predict(
    model_type = tolower(model_type),
    X_train = x_all[train_idx, , drop = FALSE],
    y_train = y_all[train_idx],
    X_test = x_all[tst, , drop = FALSE],
    response_family = fam,
    model_params = model_params,
    prefer_gpu = FALSE
  )
}

gp_met_python_cv_predict_batch <- function(pheno_data,
                                           response,
                                           gen_name,
                                           heter_groups,
                                           response_family = "gaussian",
                                           model_type,
                                           gmatrix = NULL,
                                           omic1_kernel = NULL,
                                           omic2_kernel = NULL,
                                           omic3_kernel = NULL,
                                           kernel_list = NULL,
                                           omics_kernel_label = list(
                                             omic1_kernel = NULL,
                                             omic2_kernel = NULL,
                                             omic3_kernel = NULL
                                           ),
                                           folds,
                                           model_params = list(),
                                           met_kernel_var_explained = 0.95,
                                           met_kernel_min_ev = 1e-8,
                                           met_kernel_max_pcs = NULL) {
  if (is.null(folds) || !length(folds)) {
    return(list())
  }
  feat_res <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )
  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feat_res$feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  fam <- gp_resolve_response_family(response_family, y = long_dat[[response]])
  y_all <- if (identical(fam, "gaussian")) as.numeric(long_dat[[response]]) else long_dat[[response]]
  observed_idx <- which(!is.na(y_all))
  jobs <- lapply(folds, function(tst) {
    tst <- as.integer(tst)
    train_idx <- setdiff(observed_idx, tst)
    if (!length(train_idx)) {
      stop("No observed training rows available for MET ", model_type, " CV.", call. = FALSE)
    }
    list(
      model_type = tolower(model_type),
      y_train = y_all[train_idx],
      build = function() list(              # built per job by the batch writer
        X_train = x_all[train_idx, , drop = FALSE],
        X_test = x_all[tst, , drop = FALSE]
      ),
      response_family = fam,
      model_params = model_params,
      prefer_gpu = FALSE
    )
  })
  gp_py_ml_fit_predict_batch(jobs)
}

gp_met_catboost_cv_predict <- function(...) {
  gp_met_python_cv_predict(model_type = "CatBoost", ...)
}
