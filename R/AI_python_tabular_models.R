gp_py_boosted_tree_local_uncertainty_models <- function() {
  sort(c("xgboost", "catboost", "lightgbm"))
}

gp_py_uses_local_target_uncertainty <- function(model_type) {
  model_type <- tolower(trimws(as.character(model_type %||% "")[1L]))
  nzchar(model_type) &&
    model_type %in% gp_py_boosted_tree_local_uncertainty_models()
}

gp_python_tabular_train_predict <- function(data_label_geno,
                                            indices,
                                            test_geno,
                                            model_type,
                                            model_params = list(),
                                            response_family = "gaussian",
                                            class_levels = NULL) {
  gp_py_ml_bootstrap_predict(
    data_label_geno = data_label_geno,
    indices = indices,
    test_geno = test_geno,
    model_type = model_type,
    response_family = response_family,
    model_params = model_params,
    class_levels = class_levels
  )
}

gp_python_tabular_model <- function(pheno_object = NULL,
                                    geno_omic_object = NULL,
                                    geno_omic_test_object = NULL,
                                    response = NULL,
                                    gen_name = NULL,
                                    response_family = "gaussian",
                                    model_type,
                                    model_label,
                                    model_params = list(),
                                    para_tunning = FALSE,
                                    tune_param_grid = NULL,
                                    tune_folds = 5L,
                                    tune_metric = NULL,
                                    scaling = TRUE,
                                    centering = FALSE,
                                    CI_width_thresholds = c(0.33, 0.66),
                                    high_reliability_thres = 0.9,
                                    low_reliability_thres = 0.5,
                                    n_bootstrap = 30,
                                    system_database = FALSE) {
  rng_seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng_seed_before <- if (rng_seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    gp_clear_invalid_random_seed(rng_seed_before, rng_seed_existed),
    add = TRUE
  )

  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  is_binary <- identical(fam, "binary")
  is_multiclass <- identical(fam, "multiclass")
  is_classification <- fam %in% c("binary", "multiclass")

  if (isTRUE(para_tunning) && (is.null(tune_param_grid) || !length(tune_param_grid))) {
    stop(model_label, " tuning requires a non-empty tune_param_grid.", call. = FALSE)
  }

  if (is.null(geno_omic_object) || is.null(pheno_object)) {
    stop("Training set missing.", call. = FALSE)
  }

  prep <- gp_ml_preprocess_predictors(geno_omic_object, scaling = scaling, centering = centering)
  geno_omic_object <- prep$data
  GID <- rownames(geno_omic_object)

  test_label <- NULL
  if (!is.null(geno_omic_test_object)) {
    test_label <- rownames(geno_omic_test_object)
    geno_omic_test_object <- gp_ml_apply_preprocessor(geno_omic_test_object, prep)
    geno_omic_test_object <- rbind(geno_omic_object, geno_omic_test_object)
    GID <- rownames(geno_omic_test_object)
  }

  y_train <- pheno_object[, response]
  if (is_classification) {
    class_levels <- gp_py_ml_class_levels(y_train, fam)
    y_train_processed <- as.character(y_train)
    y_scaler <- NULL
  } else {
    y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
    y_train_processed <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
    class_levels <- NULL
  }

  tuning_summary <- NULL
  if (isTRUE(para_tunning)) {
    if (is.null(tune_param_grid) || !length(tune_param_grid)) {
      stop(model_label, " tuning requires a non-empty tune_param_grid.", call. = FALSE)
    }
    tune_res <- gp_grid_tune_cv(
      y = y_train,
      response_family = fam,
      param_grid = tune_param_grid,
      tune_metric = tune_metric,
      nfolds = tune_folds,
      random_state = gp_ml_random_state(),
      predict_fun = function(tst, params) {
        pred <- gp_py_ml_fit_predict(
          model_type = model_type,
          X_train = geno_omic_object[-tst, , drop = FALSE],
          y_train = y_train[-tst],
          X_test = geno_omic_object[tst, , drop = FALSE],
          response_family = fam,
          model_params = utils::modifyList(model_params, params),
          prefer_gpu = FALSE
        )
        list(pred = pred, prob = attr(pred, "probabilities"))
      },
      batch_predict_fun = function(specs) {
        bridge_jobs <- lapply(specs, function(spec) {
          tst <- spec$tst
          list(
            id = paste(spec$param_index, spec$fold_index, sep = "_"),
            model_type = model_type,
            y_train = y_train[-tst],
            # built per job by the batch writer (grid x folds copies otherwise)
            build = function() list(
              X_train = geno_omic_object[-tst, , drop = FALSE],
              X_test = geno_omic_object[tst, , drop = FALSE]
            ),
            response_family = fam,
            model_params = utils::modifyList(model_params, spec$params),
            class_levels = class_levels
          )
        })
        preds <- gp_py_ml_fit_predict_batch(bridge_jobs)
        lapply(preds, function(pred) list(pred = pred, prob = attr(pred, "probabilities")))
      }
    )
    model_params <- utils::modifyList(model_params, tune_res$best_params)
    tuning_summary <- tune_res
  }

  cv_risk_calibration <- NULL
  if (identical(fam, "gaussian") && length(y_train) >= 6L) {
    seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    restore_seed <- seed_exists && is.integer(old_seed) && length(old_seed) > 1L
    gp_set_seed(gp_ml_random_state())
    cv_folds <- gp_tuning_folds(
      y = y_train,
      response_family = fam,
      nfolds = tune_folds,
      random_state = gp_ml_random_state()
    )
    cv_pred <- rep(NA_real_, length(y_train))
    cv_fold_id <- rep(NA_integer_, length(y_train))
    uncertainty_target_x <- geno_omic_test_object %||% geno_omic_object
    cv_target_predictions <- matrix(
      NA_real_,
      nrow = length(cv_folds),
      ncol = nrow(uncertainty_target_x)
    )
    # Fold matrices are built lazily, one job at a time, by the batch writer:
    # materialising all folds up front held ~5 full copies of the predictor
    # matrix at once (~1.4 GB at 1,814 x 10,346 markers).
    cv_jobs <- lapply(seq_along(cv_folds), function(i) {
      tst <- cv_folds[[i]]
      cv_fold_id[tst] <<- i
      list(
        id = paste0("risk_", i),
        model_type = model_type,
        y_train = y_train[-tst],
        build = function() list(
          X_train = geno_omic_object[-tst, , drop = FALSE],
          X_test = rbind(geno_omic_object[tst, , drop = FALSE], uncertainty_target_x)
        ),
        response_family = fam,
        model_params = model_params
      )
    })
    cv_pred_list <- tryCatch(
      gp_py_ml_fit_predict_batch(cv_jobs),
      error = function(e) {
        warning("Batched ML risk calibration failed; falling back to serial folds: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    rm(cv_jobs)
    if (!is.null(cv_pred_list)) {
      for (i in seq_along(cv_folds)) {
        values <- as.numeric(cv_pred_list[[i]])
        n_holdout <- length(cv_folds[[i]])
        cv_pred[cv_folds[[i]]] <- values[seq_len(n_holdout)]
        cv_target_predictions[i, ] <- values[n_holdout + seq_len(nrow(uncertainty_target_x))]
      }
    } else {
      for (i in seq_along(cv_folds)) {
        tst <- cv_folds[[i]]
        fold_pred <- gp_py_ml_fit_predict(
          model_type = model_type,
          X_train = geno_omic_object[-tst, , drop = FALSE],
          y_train = y_train[-tst],
          X_test = rbind(
            geno_omic_object[tst, , drop = FALSE],
            uncertainty_target_x
          ),
          response_family = fam,
          model_params = model_params,
          prefer_gpu = FALSE
        )
        values <- as.numeric(fold_pred)
        n_holdout <- length(tst)
        cv_pred[tst] <- values[seq_len(n_holdout)]
        cv_target_predictions[i, ] <- values[n_holdout + seq_len(nrow(uncertainty_target_x))]
      }
    }
    use_local_target_uncertainty <-
      gp_py_uses_local_target_uncertainty(model_type)
    cv_risk_calibration <- gp_ml_cv_rank_risk_calibration(
      predicted_cv = cv_pred,
      observed_y = y_train,
      fold_id = cv_fold_id,
      target_fold_predictions = cv_target_predictions,
      heldout_predictors = if (use_local_target_uncertainty) {
        geno_omic_object
      } else NULL,
      target_predictors = if (use_local_target_uncertainty) {
        uncertainty_target_x
      } else NULL,
      heldout_ids = if (use_local_target_uncertainty) {
        rownames(geno_omic_object)
      } else NULL,
      target_ids = if (use_local_target_uncertainty) {
        rownames(uncertainty_target_x)
      } else NULL
    )
    if (restore_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  # The published prediction is the model fitted once on ALL training lines;
  # the bootstrap refits only describe uncertainty (they used to be averaged
  # into the prediction, i.e. an extra bagging layer).
  full_fit <- gp_py_ml_fit_predict(
    model_type = model_type,
    X_train = geno_omic_object,
    y_train = y_train_processed,
    X_test = geno_omic_test_object %||% geno_omic_object,
    response_family = fam,
    model_params = model_params,
    class_levels = class_levels,
    prefer_gpu = identical(model_type, "xgboost")
  )
  full_fit_prob <- attr(full_fit, "probabilities", exact = TRUE)
  # An auto-tuned Lasso searches its penalty with LassoCV (minutes per search
  # at ~1-2k lines x 10k markers). The calibration folds above keep their own
  # search (no penalty chosen with their held-out rows); the bootstrap refits
  # reuse the penalty chosen on all training lines -- ~0.1 s per refit instead
  # of a fresh search for each of the n_bootstrap replicates.
  boot_params <- model_params
  if (identical(model_type, "lasso") && is.null(model_params$lasso_alpha) &&
      is.finite(attr(full_fit, "selected_alpha", exact = TRUE) %||% NA_real_)) {
    boot_params$lasso_alpha <- attr(full_fit, "selected_alpha", exact = TRUE)
  }
  boot_matrix <- gp_py_ml_bootstrap_matrix(
    model_type = model_type,
    X_train = geno_omic_object,
    y_train = y_train_processed,
    X_pred = geno_omic_test_object %||% geno_omic_object,
    response_family = fam,
    model_params = boot_params,
    class_levels = class_levels,
    n_bootstrap = n_bootstrap,
    seed = gp_ml_random_state(),
    prefer_gpu = identical(model_type, "xgboost")
  )
  boot_results <- list(t = boot_matrix, R = as.integer(n_bootstrap))

  if (!is_classification) {
    revert_scaling_ml <- function(scaled_values, scaler_mean, scaler_sd) {
      scaled_values * scaler_sd + scaler_mean
    }
    boot_results$t <- apply(
      boot_results$t,
      2,
      function(col_vec) revert_scaling_ml(col_vec, y_scaler$mean, y_scaler$std)
    )
  } else if (is_binary) {
    boot_results$t <- apply(boot_results$t, 2, function(col_vec) pmax(0, pmin(1, col_vec)))
  }

  train_test_label <- if (!is.null(geno_omic_test_object)) {
    ifelse(rownames(geno_omic_test_object) %in% test_label, "Test", "Train")
  } else {
    rep("Train", nrow(geno_omic_object))
  }

  if (is_multiclass) {
    ml_metrics <- gp_multiclass_bootstrap_metrics(
      boot_matrix = boot_results$t,
      n_obs = length(train_test_label),
      class_levels = class_levels,
      train_test_label = train_test_label
    )
    pred_variances <- ml_metrics$prediction_error_var
    pred_SE <- ml_metrics$standard_error
    # A classifier does not estimate a random-effect covariance component.
    # Keep the legacy prediction-table field explicitly unavailable; the
    # variance contract reports predictive probability dispersion separately.
    genetic_var <- NA_real_
    prob_mean <- ml_metrics$mean_prob
    if (!is.null(full_fit_prob) && all(class_levels %in% colnames(full_fit_prob)) &&
        nrow(full_fit_prob) == nrow(prob_mean)) {
      prob_mean <- full_fit_prob[, class_levels, drop = FALSE]
    }
    AI_pred_reverted <- class_levels[max.col(prob_mean, ties.method = "first")]
  } else {
    ml_metrics <- gp_ml_bootstrap_target_metrics(
      boot_matrix = boot_results$t,
      train_test_label = train_test_label,
      observed_y = if (identical(fam, "gaussian")) y_train else NULL
    )
    pred_variances <- ml_metrics$prediction_error_var
    pred_SE <- ml_metrics$standard_error
    AI_pred_reverted <- ml_metrics$predicted_mean
    if (identical(fam, "gaussian") && length(full_fit) == length(AI_pred_reverted)) {
      AI_pred_reverted <- as.numeric(full_fit) * y_scaler$std + y_scaler$mean
    } else if (is_binary && !is.null(full_fit_prob)) {
      positive <- as.character(class_levels[[min(2L, length(class_levels))]])
      if (positive %in% colnames(full_fit_prob) && nrow(full_fit_prob) == length(AI_pred_reverted)) {
        AI_pred_reverted <- pmax(0, pmin(1, as.numeric(full_fit_prob[, positive])))
      }
    }
    reference_variance <- ml_metrics$reference_variance
    genetic_var <- if (is_classification) NA_real_ else reference_variance
  }

  result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
    boot_results = boot_results,
    CI_width_thresholds = CI_width_thresholds
  )

  if (is_binary) {
    AI_preds <- gp_ml_classification_prediction_output(
      ids = GID,
      predicted_prob = AI_pred_reverted,
      train_test_label = train_test_label,
      pred_se = pred_SE,
      pred_variances = pred_variances,
      genetic_var = genetic_var,
      lower_bound = result_rel_MPIW$lower_bound,
      upper_bound = result_rel_MPIW$upper_bound,
      uncertainty = result_rel_MPIW$Uncertainty,
      uncertainty_remarks = result_rel_MPIW$reliability_remarks,
      class_levels = class_levels,
      gen_name = gen_name,
      high_confidence = high_reliability_thres,
      low_confidence = low_reliability_thres
    )
    diagnostic_plots <- diagnostic_plot_true_prediction(
      boot_results = boot_results,
      GID_names = GID,
      CI_width_thresholds = CI_width_thresholds,
      predictions = AI_pred_reverted,
      standard_errors = pred_SE,
      prediction_error_var = pred_variances,
      genetic_var = genetic_var,
      confidence_level = 0.95,
      model_for_CI_cal = "ML",
      response_family = fam,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      system_database = system_database
    )
  } else if (is_multiclass) {
    AI_preds <- gp_multiclass_prediction_output(
      ids = GID,
      prob_mean = prob_mean,
      train_test_label = train_test_label,
      pred_se = pred_SE,
      pred_variances = pred_variances,
      genetic_var = genetic_var,
      lower_bound = ml_metrics$lower_conf,
      upper_bound = ml_metrics$upper_conf,
      uncertainty = result_rel_MPIW$Uncertainty,
      uncertainty_remarks = result_rel_MPIW$reliability_remarks,
      class_levels = class_levels,
      gen_name = gen_name,
      high_confidence = high_reliability_thres,
      low_confidence = low_reliability_thres
    )
    diagnostic_plots <- NULL
  } else {
    AI_preds <- gp_ml_gaussian_prediction_output(
      ids = GID,
      predicted_value = AI_pred_reverted,
      train_test_label = train_test_label,
      pred_se = pred_SE,
      pred_variances = pred_variances,
      lower_bound = result_rel_MPIW$lower_bound,
      upper_bound = result_rel_MPIW$upper_bound,
      uncertainty = result_rel_MPIW$Uncertainty,
      uncertainty_remarks = result_rel_MPIW$reliability_remarks,
      reference_variance = reference_variance,
      observed_y = y_train,
      boot_results = boot_results,
      risk_calibration = cv_risk_calibration,
      predictor_matrix = geno_omic_test_object %||% geno_omic_object,
      gen_name = gen_name,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      confidence_level = 0.95,
      require_heldout_calibration = TRUE
    )
    diagnostic_plots <- if (!is.null(gp_ml_heldout_residual_summary(
      cv_risk_calibration,
      confidence_level = 0.95
    ))) {
      diagnostic_plot_true_prediction(
        boot_results = boot_results,
        GID_names = GID,
        CI_width_thresholds = CI_width_thresholds,
        predictions = AI_pred_reverted,
        standard_errors = pred_SE,
        prediction_error_var = pred_variances,
        reference_variance = reference_variance,
        lower_bound = result_rel_MPIW$lower_bound,
        upper_bound = result_rel_MPIW$upper_bound,
        observed_y = y_train,
        train_test_label = train_test_label,
        confidence_level = 0.95,
        risk_calibration = cv_risk_calibration,
        require_heldout_calibration = TRUE,
        model_for_CI_cal = "ML",
        high_reliability_thres = high_reliability_thres,
        low_reliability_thres = low_reliability_thres,
        system_database = system_database
      )
    } else {
      NULL
    }
  }

  model_para <- data.frame(
    stat = names(model_params),
    summary = vapply(model_params, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  if (!is.null(tuning_summary)) {
    model_para <- rbind(
      model_para,
      data.frame(stat = "tune_metric", summary = tuning_summary$metric, stringsAsFactors = FALSE),
      data.frame(stat = "tune_best_score", summary = as.character(tuning_summary$best_score), stringsAsFactors = FALSE)
    )
  }
  if (!nrow(model_para)) {
    model_para <- data.frame(stat = "model", summary = model_label, stringsAsFactors = FALSE)
  }

  variance_components <- if (identical(fam, "gaussian")) {
    gp_ml_gaussian_variance_components(
      boot_matrix = boot_results$t,
      train_test_label = train_test_label,
      observed_y = y_train,
      predicted_value = AI_pred_reverted,
      reference_variance = reference_variance,
      prediction_error_var = AI_preds[["PEV"]],
      marginal_prediction_mse = attr(
        AI_preds,
        "uncertainty_calibration"
      )$calibration_mean_pev
    )
  } else {
    gp_predictive_model_variance_components(
      AI_preds,
      response_family = fam
    )
  }

  list(
    model_parameters = model_para,
    predicted_values = AI_preds,
    diagnostic_plots = diagnostic_plots,
    variance_components = variance_components,
    Variance_components = variance_components
  )
}

gp_one_hot_env_block <- function(env_vec, prefix = "Env") {
  env_df <- data.frame(env = as.factor(env_vec), stringsAsFactors = FALSE)
  mm <- stats::model.matrix(~ env - 1, data = env_df)
  mm <- as.data.frame(mm, check.names = FALSE, stringsAsFactors = FALSE)
  names(mm) <- paste0(prefix, "_", make.names(sub("^env", "", names(mm))))
  mm
}

gp_python_tabular_met_model <- function(pheno_object = NULL,
                                        response = NULL,
                                        gen_name = NULL,
                                        heter_groups = NULL,
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
                                        response_family = "gaussian",
                                        model_type,
                                        model_label,
                                        model_params = list(),
                                        para_tunning = FALSE,
                                        tune_param_grid = NULL,
                                        tune_folds = 5L,
                                        tune_metric = NULL,
                                        met_kernel_var_explained = 0.95,
                                        met_kernel_min_ev = 1e-8,
                                        met_kernel_max_pcs = NULL,
                                        n_bootstrap = 30L,
                                        internal_cv_nfolds = 5L,
                                        confidence_level = 0.95) {
  rng_seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng_seed_before <- if (rng_seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    gp_clear_invalid_random_seed(rng_seed_before, rng_seed_existed),
    add = TRUE
  )

  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  if (!fam %in% c("gaussian", "binary", "multiclass")) {
    stop(model_label, " MET path currently supports gaussian, binary, and multiclass traits only.", call. = FALSE)
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
  feature_table <- feat_res$feature_table
  if (is.null(feature_table) || !nrow(feature_table)) {
    stop(model_label, " MET path requires at least one GRM/kernel source.", call. = FALSE)
  }

  long_dat <- gp_build_met_long_data(
    pheno_data = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernel_feature_table = feature_table
  )
  env_block <- gp_one_hot_env_block(long_dat[[heter_groups]])
  x_all <- cbind(
    long_dat[, setdiff(names(long_dat), c(response, gen_name, heter_groups, ".met_row_id")), drop = FALSE],
    env_block
  )
  y_all <- long_dat[[response]]
  train_idx <- which(!is.na(y_all))
  test_idx <- which(is.na(y_all))
  if (!length(train_idx)) {
    stop(model_label, " MET path requires non-missing training observations.", call. = FALSE)
  }
  class_levels <- gp_py_ml_class_levels(y_all[train_idx], fam)
  if (identical(fam, "gaussian")) {
    n_bootstrap <- suppressWarnings(as.integer(n_bootstrap[[1L]]))
    if (!is.finite(n_bootstrap) || n_bootstrap < 2L) {
      stop(model_label, " MET Gaussian uncertainty requires n_bootstrap >= 2.",
           call. = FALSE)
    }
  }

  tuning_summary <- NULL
  if (isTRUE(para_tunning)) {
    tune_res <- gp_grid_tune_cv(
      y = y_all[train_idx],
      response_family = fam,
      param_grid = tune_param_grid,
      tune_metric = tune_metric,
      nfolds = tune_folds,
      random_state = gp_ml_random_state(),
      predict_fun = function(tst, params) {
        full_train <- train_idx
        fit_idx <- full_train[-tst]
        pred_idx <- full_train[tst]
        pred <- gp_py_ml_fit_predict(
          model_type = model_type,
          X_train = x_all[fit_idx, , drop = FALSE],
          y_train = y_all[fit_idx],
          X_test = x_all[pred_idx, , drop = FALSE],
          response_family = fam,
          model_params = utils::modifyList(model_params, params),
          class_levels = class_levels,
          prefer_gpu = FALSE
        )
        list(pred = pred, prob = attr(pred, "probabilities"))
      },
      batch_predict_fun = function(specs) {
        full_train <- train_idx
        bridge_jobs <- lapply(specs, function(spec) {
          fit_idx <- full_train[-spec$tst]
          pred_idx <- full_train[spec$tst]
          list(
            id = paste(spec$param_index, spec$fold_index, sep = "_"),
            model_type = model_type,
            y_train = y_all[fit_idx],
            build = function() list(              # built per job by the batch writer
              X_train = x_all[fit_idx, , drop = FALSE],
              X_test = x_all[pred_idx, , drop = FALSE]
            ),
            response_family = fam,
            model_params = utils::modifyList(model_params, spec$params),
            class_levels = class_levels
          )
        })
        preds <- gp_py_ml_fit_predict_batch(bridge_jobs)
        lapply(preds, function(pred) list(pred = pred, prob = attr(pred, "probabilities")))
      }
    )
    model_params <- utils::modifyList(model_params, tune_res$best_params)
    tuning_summary <- tune_res
  }

  pred_all <- gp_py_ml_fit_predict(
    model_type = model_type,
    X_train = x_all[train_idx, , drop = FALSE],
    y_train = y_all[train_idx],
    X_test = x_all,
    response_family = fam,
    model_params = model_params,
    class_levels = class_levels,
    prefer_gpu = FALSE
  )
  train_test_label <- ifelse(seq_len(nrow(long_dat)) %in% test_idx, "Test", "Train")
  if (identical(fam, "binary")) {
    pred_prob <- as.numeric(pred_all)
    prob_mat <- attr(pred_all, "probabilities", exact = TRUE)
    if (is.null(prob_mat)) {
      prob_mat <- cbind(1 - pred_prob, pred_prob)
      colnames(prob_mat) <- class_levels
    } else {
      prob_mat <- as.matrix(prob_mat)
      if (ncol(prob_mat) == 1L) {
        prob_mat <- cbind(1 - pred_prob, pred_prob)
        colnames(prob_mat) <- class_levels
      } else {
        colnames(prob_mat) <- class_levels
      }
    }
    cls <- gp_classification_prediction_summary(
      prob = pred_prob,
      class_levels = class_levels
    )
    pred_df <- data.frame(
      GID = long_dat[[gen_name]],
      Env = long_dat[[heter_groups]],
      Predicted_value = pred_prob,
      Train_Test_Label = train_test_label,
      cls,
      stringsAsFactors = FALSE
    )
    prob_df <- as.data.frame(prob_mat, stringsAsFactors = FALSE)
    names(prob_df) <- gp_prob_column_names(class_levels)
    pred_df <- cbind(pred_df, prob_df)
  } else if (identical(fam, "multiclass")) {
    prob_mat <- attr(pred_all, "probabilities", exact = TRUE)
    if (is.null(prob_mat)) {
      prob_mat <- as.matrix(pred_all)
    } else {
      prob_mat <- as.matrix(prob_mat)
    }
    colnames(prob_mat) <- class_levels
    cls <- gp_classification_prediction_summary(
      prob = prob_mat,
      class_levels = class_levels
    )
    pred_df <- data.frame(
      GID = long_dat[[gen_name]],
      Env = long_dat[[heter_groups]],
      Predicted_value = cls$Predicted_class,
      Train_Test_Label = train_test_label,
      cls,
      stringsAsFactors = FALSE
    )
    prob_df <- as.data.frame(prob_mat, stringsAsFactors = FALSE)
    names(prob_df) <- gp_prob_column_names(class_levels)
    pred_df <- cbind(pred_df, prob_df)
  } else {
    uncertainty_seed <- gp_ml_random_state()
    bootstrap_matrix <- gp_py_ml_bootstrap_matrix(
      model_type = model_type,
      X_train = x_all[train_idx, , drop = FALSE],
      y_train = y_all[train_idx],
      X_pred = x_all,
      response_family = fam,
      model_params = model_params,
      class_levels = class_levels,
      n_bootstrap = n_bootstrap,
      seed = uncertainty_seed,
      prefer_gpu = FALSE
    )
    heldout_calibration <- gp_met_grouped_crossfit_calibration(
      y = y_all,
      observed_idx = train_idx,
      group_id = long_dat[[gen_name]],
      environment = long_dat[[heter_groups]],
      n_targets = nrow(x_all),
      predictor_matrix = x_all,
      row_id = long_dat$.met_row_id,
      fit_predict = function(fit_idx, prediction_idx, fold_index) {
        gp_py_ml_fit_predict(
          model_type = model_type,
          X_train = x_all[fit_idx, , drop = FALSE],
          y_train = y_all[fit_idx],
          X_test = x_all[prediction_idx, , drop = FALSE],
          response_family = fam,
          model_params = model_params,
          class_levels = class_levels,
          prefer_gpu = FALSE
        )
      },
      nfolds = internal_cv_nfolds,
      seed = uncertainty_seed
    )
    uncertainty <- gp_met_bootstrap_uncertainty_columns(
      predicted_value = as.numeric(pred_all),
      bootstrap_matrix = bootstrap_matrix,
      observed_y = y_all,
      train_test_label = train_test_label,
      heldout_calibration = heldout_calibration,
      n_bootstrap = n_bootstrap,
      confidence_level = confidence_level
    )
    pred_df <- data.frame(
      GID = long_dat[[gen_name]],
      Env = long_dat[[heter_groups]],
      Predicted_value = as.numeric(pred_all),
      Train_Test_Label = train_test_label,
      Observed_value = as.numeric(y_all),
      stringsAsFactors = FALSE
    )
    pred_df <- cbind(pred_df, uncertainty)
  }
  names(pred_df)[1] <- gen_name

  model_para <- data.frame(
    stat = names(model_params),
    summary = vapply(model_params, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  model_para <- rbind(
    model_para,
    data.frame(stat = "met_feature_source", summary = "kernel_eigen", stringsAsFactors = FALSE),
    data.frame(stat = "met_kernel_var_explained", summary = as.character(met_kernel_var_explained), stringsAsFactors = FALSE),
    data.frame(stat = "met_num_feature_blocks", summary = as.character(nrow(feat_res$feature_summary)), stringsAsFactors = FALSE),
    data.frame(stat = "met_total_kernel_pcs", summary = as.character(ncol(feature_table)), stringsAsFactors = FALSE),
    data.frame(stat = "met_bootstrap", summary = as.character(identical(fam, "gaussian")), stringsAsFactors = FALSE),
    data.frame(stat = "met_bootstrap_replicates", summary = if (identical(fam, "gaussian")) as.character(n_bootstrap) else "0", stringsAsFactors = FALSE),
    data.frame(stat = "met_uncertainty_calibration", summary = if (identical(fam, "gaussian")) "CV1_by_GID_grouped_internal_crossfit" else "not_applicable", stringsAsFactors = FALSE),
    data.frame(stat = "response_family", summary = fam, stringsAsFactors = FALSE)
  )
  met_kernel_names <- as.character(feat_res$feature_summary$source)
  met_kernel_names[met_kernel_names == "GRM"] <- "gmatrix"
  model_para <- rbind(
    model_para,
    gp_multi_kernel_parameter_rows(
      kernel_names = met_kernel_names,
      strategy = "concatenated_kernel_eigenfeature_blocks"
    )
  )
  if (!is.null(tuning_summary)) {
    model_para <- rbind(
      model_para,
      data.frame(stat = "tune_metric", summary = tuning_summary$metric, stringsAsFactors = FALSE),
      data.frame(stat = "tune_best_score", summary = as.character(tuning_summary$best_score), stringsAsFactors = FALSE)
    )
  }

  list(
    model_parameters = model_para,
    predicted_values = pred_df,
    diagnostic_plots = NULL,
    met_long_data = long_dat,
    met_feature_summary = feat_res$feature_summary
  )
}

gp_python_tabular_met_wrapper <- function(model_type,
                                          model_label,
                                          pheno_object = NULL,
                                          response = NULL,
                                          gen_name = NULL,
                                          heter_groups = NULL,
                                          response_family = "gaussian",
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
                                          model_params = list(),
                                          para_tunning = FALSE,
                                          tune_param_grid = NULL,
                                          met_kernel_var_explained = 0.95,
                                          met_kernel_min_ev = 1e-8,
                                          met_kernel_max_pcs = NULL,
                                          n_bootstrap = 30L,
                                          internal_cv_nfolds = 5L,
                                          confidence_level = 0.95) {
  gp_python_tabular_met_model(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response_family = response_family,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    model_type = model_type,
    model_label = model_label,
    model_params = model_params,
    para_tunning = para_tunning,
    tune_param_grid = tune_param_grid,
    met_kernel_var_explained = met_kernel_var_explained,
    met_kernel_min_ev = met_kernel_min_ev,
    met_kernel_max_pcs = met_kernel_max_pcs,
    n_bootstrap = n_bootstrap,
    internal_cv_nfolds = internal_cv_nfolds,
    confidence_level = confidence_level
  )
}

AI_CatBoost <- function(pheno_object = NULL,
                        geno_omic_object = NULL,
                        geno_omic_test_object = NULL,
                        response = NULL,
                        gen_name = NULL,
                        response_family = "gaussian",
                        message = TRUE,
                        scaling = TRUE,
                        centering = FALSE,
                        para_tunning = FALSE,
                        catboost_paras_tunning = list(
                          catboost_iterations = c(300, 500),
                          catboost_depth = c(4, 6, 8),
                          catboost_learning_rate = c(0.03, 0.1),
                          catboost_l2_leaf_reg = c(1, 3, 5)
                        ),
                        catboost_iterations = 500,
                        catboost_depth = 6,
                        catboost_learning_rate = 0.03,
                        catboost_l2_leaf_reg = 3,
                        catboost_thread_count = 1L,
                        CI_width_thresholds = c(0.33, 0.66),
                        high_reliability_thres = 0.9,
                        low_reliability_thres = 0.5,
                        n_bootstrap = 30,
                        system_database = FALSE,
                        ...) {
  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = response_family,
    model_type = "catboost",
    model_label = "CatBoost",
    model_params = list(
      catboost_iterations = catboost_iterations,
      catboost_depth = catboost_depth,
      catboost_learning_rate = catboost_learning_rate,
      catboost_l2_leaf_reg = catboost_l2_leaf_reg,
      catboost_thread_count = as.integer(catboost_thread_count %||% 1L)
    ),
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) catboost_paras_tunning else NULL,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}

AI_CatBoost_MET <- function(pheno_object = NULL,
                            response = NULL,
                            gen_name = NULL,
                            heter_groups = NULL,
                            response_family = "gaussian",
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
                            para_tunning = FALSE,
                            catboost_paras_tunning = list(
                              catboost_iterations = c(300, 500),
                              catboost_depth = c(4, 6, 8),
                              catboost_learning_rate = c(0.03, 0.1),
                              catboost_l2_leaf_reg = c(1, 3, 5)
                            ),
                            catboost_iterations = 500,
                            catboost_depth = 6,
                            catboost_learning_rate = 0.03,
                            catboost_l2_leaf_reg = 3,
                            catboost_thread_count = 1L,
                            met_kernel_var_explained = 0.95,
                            met_kernel_min_ev = 1e-8,
                            met_kernel_max_pcs = NULL,
                            n_bootstrap = 30L,
                            internal_cv_nfolds = 5L,
                            confidence_level = 0.95) {
  gp_python_tabular_met_wrapper(
    model_type = "catboost",
    model_label = "CatBoost",
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response_family = response_family,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    model_params = list(
      catboost_iterations = catboost_iterations,
      catboost_depth = catboost_depth,
      catboost_learning_rate = catboost_learning_rate,
      catboost_l2_leaf_reg = catboost_l2_leaf_reg,
      catboost_thread_count = as.integer(catboost_thread_count %||% 1L)
    ),
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) catboost_paras_tunning else NULL,
    met_kernel_var_explained = met_kernel_var_explained,
    met_kernel_min_ev = met_kernel_min_ev,
    met_kernel_max_pcs = met_kernel_max_pcs,
    n_bootstrap = n_bootstrap,
    internal_cv_nfolds = internal_cv_nfolds,
    confidence_level = confidence_level
  )
}

AI_LightGBM_MET <- function(pheno_object = NULL,
                            response = NULL,
                            gen_name = NULL,
                            heter_groups = NULL,
                            response_family = "gaussian",
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
                            para_tunning = FALSE,
                            lightgbm_paras_tunning = list(
                              lightgbm_nrounds = c(100, 300),
                              lightgbm_learning_rate = c(0.03, 0.1),
                              lightgbm_num_leaves = c(15, 31, 63),
                              lightgbm_feature_fraction = c(0.8, 1.0),
                              lightgbm_bagging_fraction = c(0.8, 1.0)
                            ),
                            lightgbm_nrounds = 100,
                            lightgbm_learning_rate = 0.05,
                            lightgbm_num_leaves = 31,
                            lightgbm_feature_fraction = 1.0,
                            lightgbm_bagging_fraction = 1.0,
                            lightgbm_min_data_in_leaf = 20,
                            lightgbm_lambda_l1 = 0,
                            lightgbm_lambda_l2 = 0,
                            lightgbm_nthread = 1L,
                            met_kernel_var_explained = 0.95,
                            met_kernel_min_ev = 1e-8,
                            met_kernel_max_pcs = NULL,
                            n_bootstrap = 30L,
                            internal_cv_nfolds = 5L,
                            confidence_level = 0.95) {
  gp_python_tabular_met_wrapper(
    model_type = "lightgbm",
    model_label = "LightGBM",
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response_family = response_family,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    model_params = list(
      lightgbm_nrounds = lightgbm_nrounds,
      lightgbm_learning_rate = lightgbm_learning_rate,
      lightgbm_num_leaves = lightgbm_num_leaves,
      lightgbm_feature_fraction = lightgbm_feature_fraction,
      lightgbm_bagging_fraction = lightgbm_bagging_fraction,
      lightgbm_min_data_in_leaf = lightgbm_min_data_in_leaf,
      lightgbm_lambda_l1 = lightgbm_lambda_l1,
      lightgbm_lambda_l2 = lightgbm_lambda_l2,
      lightgbm_nthread = as.integer(lightgbm_nthread %||% 1L)
    ),
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) lightgbm_paras_tunning else NULL,
    met_kernel_var_explained = met_kernel_var_explained,
    met_kernel_min_ev = met_kernel_min_ev,
    met_kernel_max_pcs = met_kernel_max_pcs,
    n_bootstrap = n_bootstrap,
    internal_cv_nfolds = internal_cv_nfolds,
    confidence_level = confidence_level
  )
}

AI_Xgb_MET <- function(pheno_object = NULL,
                       response = NULL,
                       gen_name = NULL,
                       heter_groups = NULL,
                       response_family = "gaussian",
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
                       para_tunning = FALSE,
                       xgb_paras_tunning = list(
                         Iter_tune = seq(100, 500, 100),
                         learning_rate_tune = seq(0.01, 0.1, 0.01),
                         max_depth = seq(3, 15, 2),
                         xgb_gamma = seq(0, 1, 0.01),
                         colsample_bytree = seq(0.1, 1, 0.1),
                         min_child_weight = seq(1, 10, 2),
                         subsample = seq(0.2, 1, 0.1),
                         L2_tune = seq(0, 1, 0.01),
                         L1_tune = seq(0, 1, 0.01)
                       ),
                       learning_rate = 0.01,
                       max_depth = 6,
                       subsample = 0.7,
                       xgb_booster = "gbtree",
                       xgb_alpha = 0.001,
                       xgb_lambda = 1.0,
                       xgb_gamma = 0.01,
                       min_child_weight = 1,
                       iteration = 100,
                       xgb_nthread = 1L,
                       colsample_bytree = 0.7,
                       met_kernel_var_explained = 0.95,
                       met_kernel_min_ev = 1e-8,
                       met_kernel_max_pcs = NULL,
                       n_bootstrap = 30L,
                       internal_cv_nfolds = 5L,
                       confidence_level = 0.95) {
  tune_grid <- if (isTRUE(para_tunning)) {
    list(
      nrounds = xgb_paras_tunning$Iter_tune %||% c(iteration),
      eta = xgb_paras_tunning$learning_rate_tune %||% c(learning_rate),
      max_depth = xgb_paras_tunning$max_depth %||% c(max_depth),
      xgb_gamma = xgb_paras_tunning$xgb_gamma %||% c(xgb_gamma),
      colsample_bytree = xgb_paras_tunning$colsample_bytree %||% c(colsample_bytree),
      min_child_weight = xgb_paras_tunning$min_child_weight %||% c(min_child_weight),
      subsample = xgb_paras_tunning$subsample %||% c(subsample),
      xgb_alpha = xgb_paras_tunning$L1_tune %||% c(xgb_alpha),
      xgb_lambda = xgb_paras_tunning$L2_tune %||% c(xgb_lambda)
    )
  } else {
    NULL
  }
  gp_python_tabular_met_wrapper(
    model_type = "xgboost",
    model_label = "Xgboost",
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response_family = response_family,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    model_params = list(
      nrounds = iteration,
      eta = learning_rate,
      max_depth = max_depth,
      subsample = subsample,
      xgb_gamma = xgb_gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      xgb_alpha = xgb_alpha,
      xgb_lambda = xgb_lambda,
      xgb_booster = xgb_booster,
      xgb_nthread = as.integer(xgb_nthread %||% 1L)
    ),
    para_tunning = para_tunning,
    tune_param_grid = tune_grid,
    met_kernel_var_explained = met_kernel_var_explained,
    met_kernel_min_ev = met_kernel_min_ev,
    met_kernel_max_pcs = met_kernel_max_pcs,
    n_bootstrap = n_bootstrap,
    internal_cv_nfolds = internal_cv_nfolds,
    confidence_level = confidence_level
  )
}

AI_randomForest_MET <- function(pheno_object = NULL,
                                response = NULL,
                                gen_name = NULL,
                                heter_groups = NULL,
                                response_family = "gaussian",
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
                                para_tunning = FALSE,
                                rf_paras_tunning = list(
                                  mtry = TRUE,
                                  ntree = c(500, 1000, 1500),
                                  nodesize = c(1, 5, 10),
                                  maxnodes = c(30, 50, NULL)
                                ),
                                ntree = 500,
                                mtry = NULL,
                                maxnodes = NULL,
                                nodesize = NULL,
                                rf_n_jobs = 1L,
                                met_kernel_var_explained = 0.95,
                                met_kernel_min_ev = 1e-8,
                                met_kernel_max_pcs = NULL,
                                n_bootstrap = 30L,
                                internal_cv_nfolds = 5L,
                                confidence_level = 0.95) {
  feature_table <- gp_prepare_met_kernel_features(
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    var_explained = met_kernel_var_explained,
    min_ev = met_kernel_min_ev,
    max_pcs = met_kernel_max_pcs
  )$feature_table
  default_mtry <- mtry %||% floor(sqrt(ncol(feature_table %||% matrix(0, nrow = 1, ncol = 1))))
  tune_grid <- if (isTRUE(para_tunning)) {
    list(
      ntree = rf_paras_tunning$ntree %||% c(ntree),
      mtry = if (isTRUE(rf_paras_tunning$mtry)) {
        unique(as.integer(c(sqrt(ncol(feature_table)), sqrt(ncol(feature_table)) / 2, ncol(feature_table) / 3)))
      } else {
        rf_paras_tunning$mtry %||% default_mtry
      },
      nodesize = rf_paras_tunning$nodesize %||% c(nodesize %||% 1L),
      maxnodes = rf_paras_tunning$maxnodes %||% c(maxnodes)
    )
  } else {
    NULL
  }
  gp_python_tabular_met_wrapper(
    model_type = "randomforest",
    model_label = "RandomForest",
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response_family = response_family,
    gmatrix = gmatrix,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list,
    omics_kernel_label = omics_kernel_label,
    model_params = list(
      ntree = ntree,
      mtry = default_mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    ),
    para_tunning = para_tunning,
    tune_param_grid = tune_grid,
    met_kernel_var_explained = met_kernel_var_explained,
    met_kernel_min_ev = met_kernel_min_ev,
    met_kernel_max_pcs = met_kernel_max_pcs,
    n_bootstrap = n_bootstrap,
    internal_cv_nfolds = internal_cv_nfolds,
    confidence_level = confidence_level
  )
}

AI_LightGBM <- function(pheno_object = NULL,
                        geno_omic_object = NULL,
                        geno_omic_test_object = NULL,
                        response = NULL,
                        gen_name = NULL,
                        response_family = "gaussian",
                        message = TRUE,
                        scaling = TRUE,
                        centering = FALSE,
                        para_tunning = FALSE,
                        lightgbm_paras_tunning = list(
                          lightgbm_nrounds = c(100, 300),
                          lightgbm_learning_rate = c(0.03, 0.1),
                          lightgbm_num_leaves = c(15, 31, 63),
                          lightgbm_feature_fraction = c(0.8, 1.0),
                          lightgbm_bagging_fraction = c(0.8, 1.0)
                        ),
                        lightgbm_nrounds = 100,
                        lightgbm_learning_rate = 0.05,
                        lightgbm_num_leaves = 31,
                        lightgbm_feature_fraction = 1.0,
                        lightgbm_bagging_fraction = 1.0,
                        lightgbm_min_data_in_leaf = 20,
                        lightgbm_lambda_l1 = 0,
                        lightgbm_lambda_l2 = 0,
                        lightgbm_nthread = 1L,
                        CI_width_thresholds = c(0.33, 0.66),
                        high_reliability_thres = 0.9,
                        low_reliability_thres = 0.5,
                        n_bootstrap = 30,
                        system_database = FALSE,
                        ...) {
  gp_python_tabular_model(
    pheno_object = pheno_object,
    geno_omic_object = geno_omic_object,
    geno_omic_test_object = geno_omic_test_object,
    response = response,
    gen_name = gen_name,
    response_family = response_family,
    model_type = "lightgbm",
    model_label = "LightGBM",
    model_params = list(
      lightgbm_nrounds = lightgbm_nrounds,
      lightgbm_learning_rate = lightgbm_learning_rate,
      lightgbm_num_leaves = lightgbm_num_leaves,
      lightgbm_feature_fraction = lightgbm_feature_fraction,
      lightgbm_bagging_fraction = lightgbm_bagging_fraction,
      lightgbm_min_data_in_leaf = lightgbm_min_data_in_leaf,
      lightgbm_lambda_l1 = lightgbm_lambda_l1,
      lightgbm_lambda_l2 = lightgbm_lambda_l2,
      lightgbm_nthread = as.integer(lightgbm_nthread %||% 1L)
    ),
    para_tunning = para_tunning,
    tune_param_grid = if (isTRUE(para_tunning)) lightgbm_paras_tunning else NULL,
    scaling = scaling,
    centering = centering,
    CI_width_thresholds = CI_width_thresholds,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    n_bootstrap = n_bootstrap,
    system_database = system_database
  )
}
