gp_multitrait_ml_supported_models <- function() {
  c(
    "RandomForest",
    "Ridge_Regression",
    "PartialLeastSquare",
    "Lasso",
    "SupportVectorMachine",
    "Xgboost",
    "CatBoost"
  )
}

gp_multitrait_ml_python_model_key <- function(model_type) {
  switch(
    as.character(model_type[[1]]),
    "Ridge_Regression" = "ridge",
    "PartialLeastSquare" = "pls",
    "Lasso" = "lasso",
    "SupportVectorMachine" = "svm",
    "RandomForest" = "randomforest",
    "Xgboost" = "xgboost",
    "CatBoost" = "catboost",
    stop("Unsupported multi-trait ML python model: ", model_type, call. = FALSE)
  )
}

gp_multitrait_ml_fit_trait <- function(model_type,
                                       X_train,
                                       y_train,
                                       X_test,
                                       ntree = 500,
                                       mtry = NULL,
                                       maxnodes = NULL,
                                       nodesize = NULL,
                                       lambda_rr = NULL,
                                       ncomp = 3L,
                                       svm_type = "eps-regression",
                                       svm_kernel = "Gaussian",
                                       sigma_value = NULL,
                                       C_value = 1,
                                       degree_value = 3,
                                       scale_value = 1,
                                       offset_value = 0,
                                       gamma_value = NULL,
                                       iteration = 100,
                                       learning_rate = 0.01,
                                       max_depth_xgb = 6,
                                       subsample_xgb = 0.7,
                                       xgb_booster = "gbtree",
                                       xgb_alpha = 0.001,
                                       xgb_lambda = 1.0,
                                       xgb_gamma = 0.01,
                                       min_child_weight = 1,
                                       xgb_nthread = 1L,
                                       colsample_bytree = 0.7,
                                       catboost_iterations = 500,
                                       catboost_depth = 6,
                                       catboost_learning_rate = 0.03,
                                       catboost_l2_leaf_reg = 3,
                                       catboost_thread_count = 1L,
                                       rf_n_jobs = 1L,
                                       return_job = FALSE) {
  model_type <- as.character(model_type[[1]])

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y_train)), method = c("center", "scale"))
  y_fit <- stats::predict(y_scaler, as.data.frame(as.matrix(y_train)))[, 1]
  model_key <- gp_multitrait_ml_python_model_key(model_type)
  model_params <- switch(
    model_type,
    "RandomForest" = list(
      ntree = ntree,
      mtry = mtry %||% max(1L, floor(sqrt(ncol(X_train)))),
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    ),
    "Ridge_Regression" = if (!is.null(lambda_rr)) list(ridge_alpha = lambda_rr) else list(),
    "Lasso" = list(),
    "SupportVectorMachine" = list(
      svm_type = svm_type,
      svm_kernel = svm_kernel,
      sigma_value = sigma_value,
      C_value = C_value,
      degree_value = degree_value,
      scale_value = scale_value,
      offset_value = offset_value,
      gamma_value = gamma_value
    ),
    "Xgboost" = list(
      nrounds = iteration,
      eta = learning_rate,
      max_depth = max_depth_xgb,
      subsample = subsample_xgb,
      xgb_gamma = xgb_gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      xgb_alpha = xgb_alpha,
      xgb_lambda = xgb_lambda,
      xgb_booster = xgb_booster,
      xgb_nthread = as.integer(xgb_nthread %||% 1L)
    ),
    "CatBoost" = list(
      catboost_iterations = catboost_iterations,
      catboost_depth = catboost_depth,
      catboost_learning_rate = catboost_learning_rate,
      catboost_l2_leaf_reg = catboost_l2_leaf_reg,
      catboost_thread_count = as.integer(catboost_thread_count %||% 1L)
    ),
    "PartialLeastSquare" = {
      max_comp <- max(1L, min(as.integer(ncomp %||% 3L), ncol(X_train), nrow(X_train) - 1L))
      list(ncomp = max_comp)
    },
    list()
  )

  # return_job = TRUE: hand back the fit as a batch job (predicting X_test) so
  # many fold fits can run in one Python process (gp_py_ml_fit_predict_batch);
  # finish it with gp_multitrait_ml_job_predictions().
  if (isTRUE(return_job)) {
    return(list(
      job = list(model_type = model_key, X_train = X_train, y_train = y_fit, X_test = X_test,
                 response_family = "gaussian", model_params = model_params),
      scaler = y_scaler
    ))
  }

  # One fit predicts the training and the target rows (two separate fits of
  # the same model gave the same predictions at twice the cost).
  n_train <- nrow(X_train)
  pred <- revert_scaling(as.numeric(gp_py_ml_fit_predict(
    model_type = model_key,
    X_train = X_train,
    y_train = y_fit,
    X_test = rbind(X_train, X_test),
    response_family = "gaussian",
    model_params = model_params,
    prefer_gpu = identical(model_key, "xgboost")
  )), y_scaler)

  list(
    train_pred = pred[seq_len(n_train)],
    test_pred = pred[n_train + seq_len(nrow(X_test))]
  )
}

gp_multitrait_ml_job_predictions <- function(result, scaler) {
  revert_scaling(as.numeric(result), scaler)
}

gp_multitrait_ml_summary_statistics <- function(pred_long, response, model_type) {
  counts <- do.call(rbind, lapply(response, function(tr) {
    dat <- pred_long[pred_long$Trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      observed_count = sum(dat$Train_Test_Label == "Train", na.rm = TRUE),
      predicted_count = sum(dat$Train_Test_Label == "Test", na.rm = TRUE),
      total_count = nrow(dat),
      stringsAsFactors = FALSE
    )
  }))
  data.frame(
    stat = c(
      "mode",
      "response_family",
      "model_type",
      "n_traits",
      "traits",
      "scope",
      paste0("observed_count_", counts$trait),
      paste0("predicted_count_", counts$trait)
    ),
    summary = c(
      "multi_trait_ml",
      "gaussian",
      model_type,
      length(response),
      paste(response, collapse = ", "),
      "unbalanced gaussian true prediction via auxiliary-trait multi-trait machine learning",
      counts$observed_count,
      counts$predicted_count
    ),
    stringsAsFactors = FALSE
  )
}

gp_multitrait_ml_diagnostic_plot <- function(pred_long, model_label = "RandomForest") {
  obs <- pred_long[!is.na(pred_long$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~Trait, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Multi-trait Gaussian", model_label, "observed vs predicted"),
      subtitle = "Each target trait borrows information from markers plus the other traits with training-only imputation and missingness indicators",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_multitrait_ml_prepare_aux_features <- function(aux_df, train_idx) {
  aux_df <- as.data.frame(aux_df, stringsAsFactors = FALSE)
  if (!ncol(aux_df)) {
    return(NULL)
  }

  aux_imp <- aux_df
  aux_ind <- aux_df
  for (nm in names(aux_df)) {
    miss <- is.na(aux_df[[nm]])
    train_mean <- mean(aux_df[[nm]][train_idx], na.rm = TRUE)
    if (!is.finite(train_mean)) {
      train_mean <- mean(aux_df[[nm]], na.rm = TRUE)
    }
    if (!is.finite(train_mean)) {
      train_mean <- 0
    }
    aux_imp[[nm]] <- as.numeric(aux_df[[nm]])
    aux_imp[[nm]][miss] <- train_mean
    aux_ind[[paste0(nm, "_missing")]] <- as.numeric(miss)
  }
  aux_ind <- aux_ind[, grep("_missing$", names(aux_ind)), drop = FALSE]
  out <- cbind(aux_imp, aux_ind)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- make.names(names(out), unique = TRUE)
  out
}

gp_multitrait_ml_prediction_wide <- function(pred_long, gen_name) {
  out <- stats::reshape(
    pred_long[, c(gen_name, "Trait", "Predicted_value"), drop = FALSE],
    idvar = gen_name,
    timevar = "Trait",
    direction = "wide"
  )
  names(out) <- sub("^Predicted_value\\.", "", names(out))
  out
}

gp_multitrait_ml_trait_counts <- function(pred_long, response) {
  do.call(rbind, lapply(response, function(tr) {
    dat <- pred_long[pred_long$Trait == tr, , drop = FALSE]
    data.frame(
      trait = tr,
      observed_count = sum(dat$Train_Test_Label == "Train", na.rm = TRUE),
      predicted_count = sum(dat$Train_Test_Label == "Test", na.rm = TRUE),
      total_count = nrow(dat),
      stringsAsFactors = FALSE
    )
  }))
}

gp_multitrait_ml_apply_risk <- function(pred_long,
                                        cv_pred_long,
                                        gen_name) {
  if (is.null(pred_long) || !nrow(pred_long) || is.null(cv_pred_long) || !nrow(cv_pred_long)) {
    return(pred_long)
  }
  pred_tmp <- pred_long
  cv_tmp <- cv_pred_long
  names(pred_tmp)[names(pred_tmp) == gen_name] <- "GID"
  names(cv_tmp)[names(cv_tmp) == gen_name] <- "GID"
  out <- gp_multitrait_dl_true_prediction_risk(
    pred_long = pred_tmp,
    cv_pred_long = cv_tmp,
    gen_name = "GID"
  )
  names(out)[names(out) == "GID"] <- gen_name
  out
}

gp_multitrait_ml_gaussian_model <- function(model_type = c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
                                            pheno_object,
                                            response,
                                            geno_omic_object,
                                            gen_name,
                                            ntree = 500,
                                            mtry = NULL,
                                            maxnodes = NULL,
                                            nodesize = NULL,
                                            lambda_rr = NULL,
                                            ncomp = 3L,
                                            svm_type = "eps-regression",
                                            svm_kernel = "Gaussian",
                                            sigma_value = NULL,
                                            C_value = 1,
                                            degree_value = 3,
                                            scale_value = 1,
                                            offset_value = 0,
                                            gamma_value = NULL,
                                            iteration = 100,
                                            learning_rate = 0.01,
                                            max_depth_xgb = 6,
                                            subsample_xgb = 0.7,
                                            xgb_booster = "gbtree",
                                            xgb_alpha = 0.001,
                                            xgb_lambda = 1.0,
                                            xgb_gamma = 0.01,
                                            min_child_weight = 1,
                                            xgb_nthread = 1L,
                                            colsample_bytree = 0.7,
                                            catboost_iterations = 500,
                                            catboost_depth = 6,
                                            catboost_learning_rate = 0.03,
                                            catboost_l2_leaf_reg = 3,
                                            catboost_thread_count = 1L,
                                            rf_n_jobs = 1L,
                                            internal_cv_nfolds = 5L,
                                            internal_cv_replication = 3L,
                                            random_seed = 123L) {
  model_type <- match.arg(model_type)

  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  if (!gen_name %in% names(ph)) {
    stop("The genotype ID column was not found in the phenotype data.", call. = FALSE)
  }
  if (anyDuplicated(ph[[gen_name]])) {
    stop("Multi-trait ML currently requires one phenotype row per genotype.", call. = FALSE)
  }

  ph_ids <- as.character(ph[[gen_name]])
  geno_ids <- rownames(geno_omic_object)
  matched_ph <- ph[ph_ids %in% geno_ids, c(gen_name, response), drop = FALSE]
  matched_ph <- matched_ph[!duplicated(as.character(matched_ph[[gen_name]])), , drop = FALSE]
  if (!nrow(matched_ph)) {
    stop("No overlapping genotype IDs between phenotypes and processed features.", call. = FALSE)
  }

  extra_geno_ids <- setdiff(geno_ids, as.character(matched_ph[[gen_name]]))
  if (length(extra_geno_ids)) {
    extra_ph <- as.data.frame(
      matrix(NA_real_, nrow = length(extra_geno_ids), ncol = length(response)),
      stringsAsFactors = FALSE
    )
    names(extra_ph) <- response
    extra_ph[[gen_name]] <- extra_geno_ids
    extra_ph <- extra_ph[, c(gen_name, response), drop = FALSE]
    ph <- rbind(matched_ph, extra_ph)
  } else {
    ph <- matched_ph
  }

  ids <- as.character(ph[[gen_name]])
  x_base <- as.matrix(geno_omic_object[ids, , drop = FALSE])
  storage.mode(x_base) <- "double"

  pred_rows <- vector("list", length(response))
  for (i in seq_along(response)) {
    tr <- response[[i]]
    y <- as.numeric(ph[[tr]])
    train_idx <- which(!is.na(y))
    if (!length(train_idx)) {
      stop("Trait '", tr, "' has no observed values for multi-trait ML.", call. = FALSE)
    }

    aux_traits <- setdiff(response, tr)
    aux_df <- ph[, aux_traits, drop = FALSE]
    aux_ready <- gp_multitrait_ml_prepare_aux_features(aux_df, train_idx = train_idx)
    x_all <- if (is.null(aux_ready)) {
      x_base
    } else {
      as.matrix(cbind(x_base, aux_ready))
    }
    storage.mode(x_all) <- "double"

    fit_res <- gp_multitrait_ml_fit_trait(
      model_type = model_type,
      X_train = x_all[train_idx, , drop = FALSE],
      y_train = y[train_idx],
      X_test = x_all[setdiff(seq_len(nrow(ph)), train_idx), , drop = FALSE],
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      lambda_rr = lambda_rr,
      ncomp = ncomp,
      svm_type = svm_type,
      svm_kernel = svm_kernel,
      sigma_value = sigma_value,
      C_value = C_value,
      degree_value = degree_value,
      scale_value = scale_value,
      offset_value = offset_value,
      gamma_value = gamma_value,
      iteration = iteration,
      learning_rate = learning_rate,
      max_depth_xgb = max_depth_xgb,
      subsample_xgb = subsample_xgb,
      xgb_booster = xgb_booster,
      xgb_alpha = xgb_alpha,
      xgb_lambda = xgb_lambda,
      xgb_gamma = xgb_gamma,
      min_child_weight = min_child_weight,
      xgb_nthread = xgb_nthread,
      colsample_bytree = colsample_bytree,
      catboost_iterations = catboost_iterations,
      catboost_depth = catboost_depth,
      catboost_learning_rate = catboost_learning_rate,
      catboost_l2_leaf_reg = catboost_l2_leaf_reg,
      catboost_thread_count = catboost_thread_count,
      rf_n_jobs = rf_n_jobs
    )

    pred_all <- rep(NA_real_, nrow(ph))
    pred_all[train_idx] <- fit_res$train_pred
    test_idx <- setdiff(seq_len(nrow(ph)), train_idx)
    if (length(test_idx)) {
      pred_all[test_idx] <- fit_res$test_pred
    }

    pred_rows[[i]] <- data.frame(
      stringsAsFactors = FALSE,
      id_value = ids,
      Trait = tr,
      Observed_value = y,
      Predicted_value = pred_all,
      Train_Test_Label = ifelse(is.na(y), "Test", "Train")
    )
    names(pred_rows[[i]])[1] <- gen_name
  }

  pred_long <- do.call(rbind, pred_rows)
  rownames(pred_long) <- NULL

  internal_cv <- gp_multitrait_ml_gaussian_cv(
    model_type = model_type,
    pheno_object = ph,
    response = response,
    geno_omic_object = geno_omic_object[ids, , drop = FALSE],
    gen_name = gen_name,
    cross_validation_meth = "K-Folds",
    nfolds = internal_cv_nfolds,
    replication = internal_cv_replication,
    eval_metrics = c("root_mean_squared_error", "mean_absolute_error"),
    ntree = ntree,
    mtry = mtry,
    maxnodes = maxnodes,
    nodesize = nodesize,
    lambda_rr = lambda_rr,
    ncomp = ncomp,
    svm_type = svm_type,
    svm_kernel = svm_kernel,
    sigma_value = sigma_value,
    C_value = C_value,
    degree_value = degree_value,
    scale_value = scale_value,
    offset_value = offset_value,
    gamma_value = gamma_value,
    iteration = iteration,
    learning_rate = learning_rate,
    max_depth_xgb = max_depth_xgb,
    subsample_xgb = subsample_xgb,
    xgb_booster = xgb_booster,
    xgb_alpha = xgb_alpha,
    xgb_lambda = xgb_lambda,
    xgb_gamma = xgb_gamma,
    min_child_weight = min_child_weight,
    xgb_nthread = xgb_nthread,
    colsample_bytree = colsample_bytree,
    catboost_iterations = catboost_iterations,
    catboost_depth = catboost_depth,
    catboost_learning_rate = catboost_learning_rate,
    catboost_l2_leaf_reg = catboost_l2_leaf_reg,
    catboost_thread_count = catboost_thread_count,
    rf_n_jobs = rf_n_jobs,
    random_seed = random_seed
  )
  cv_pred_long <- do.call(rbind, lapply(internal_cv$cv_results, function(x) x$ypred_cv_Reps_all))
  target_pred_long <- do.call(rbind, lapply(
    internal_cv$cv_results,
    function(x) x$target_prediction_Reps_all
  ))
  pred_long <- gp_multitrait_ml_apply_risk(
    pred_long = pred_long,
    cv_pred_long = cv_pred_long,
    gen_name = gen_name
  )
  pred_long <- gp_ml_apply_cv_uncertainty(
    pred = pred_long,
    cv_pred = cv_pred_long,
    id_col = gen_name,
    group_col = "Trait",
    target_pred = target_pred_long
  )
  variance_components <- gp_predictive_model_variance_components(
    pred_long,
    response_family = "gaussian"
  )
  predictive_matrices <- gp_ml_predictive_matrix_bundle(
    pred = pred_long,
    id_col = gen_name,
    dimension_col = "Trait",
    error_pred = cv_pred_long,
    random_seed = random_seed
  )

  c(list(
    predicted_values = pred_long,
    multitrait_prediction_wide = gp_multitrait_ml_prediction_wide(pred_long, gen_name = gen_name),
    multitrait_trait_counts = gp_multitrait_ml_trait_counts(pred_long, response),
    diagnostic_plots = gp_multitrait_ml_diagnostic_plot(pred_long, model_label = model_type),
    variance_components = variance_components,
    Variance_components = variance_components,
    prediction_uncertainty_source = "heldout_internal_cv",
    predictive_covariance_basis =
      "descriptive predictions and heldout prediction errors; not genetic or residual covariance"
  ), predictive_matrices)
}

gp_multitrait_ml_metric_table <- function(pred_long,
                                          eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                          model_type = "RandomForest",
                                          rep = 1L) {
  out <- list()
  idx <- 1L
  for (tr in unique(pred_long$Trait)) {
    dat_tr <- pred_long[pred_long$Trait == tr & !is.na(pred_long$Observed_value), , drop = FALSE]
    if (!nrow(dat_tr)) next
    for (fd in sort(unique(dat_tr$fold))) {
      dat_fd <- dat_tr[dat_tr$fold == fd, , drop = FALSE]
      if (!nrow(dat_fd)) next
      for (metric in eval_metrics) {
        val <- tryCatch(
          evaluation_metrics(
            y_observed = dat_fd$Observed_value,
            y_predicted = dat_fd$Predicted_value,
            eval_metrics = metric,
            response_family = "gaussian"
          ),
          error = function(e) NA_real_
        )
        out[[idx]] <- data.frame(
          trait = tr,
          model = model_type,
          rep = rep,
          fold = fd,
          metric = metric,
          value = as.numeric(val),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  if (!length(out)) {
    return(data.frame(
      trait = character(),
      model = character(),
      rep = integer(),
      fold = integer(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

gp_multitrait_ml_cv_process <- function(cv_results) {
  pred_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, function(x) x$ypred_cv_Reps_all))
  if (length(pred_parts)) {
    pred_all <- do.call(rbind, pred_parts)
    pred_all <- as.data.frame(pred_all, stringsAsFactors = FALSE)
    if (!"GID" %in% names(pred_all)) {
      non_id_cols <- c(
        "Trait", "model", "Observed_value", "Predicted_value",
        "Train_Test_Label", "fold", "rep", "cv_role",
        "Prediction_SE", "Prediction_SD", "Prediction_Lower", "Prediction_Upper"
      )
      id_candidates <- setdiff(names(pred_all), non_id_cols)
      if (length(id_candidates)) {
        names(pred_all)[names(pred_all) == id_candidates[[1L]]] <- "GID"
      }
    }
    counts <- stats::aggregate(
      Predicted_value ~ Trait + Train_Test_Label,
      data = pred_all,
      FUN = length
    )
    names(counts) <- c("Trait", "Train_Test_Label", "count")
  } else {
    pred_all <- data.frame(
      Trait = character(),
      Train_Test_Label = character(),
      Predicted_value = numeric(),
      stringsAsFactors = FALSE
    )
    counts <- data.frame(
      Trait = character(),
      Train_Test_Label = character(),
      count = integer(),
      stringsAsFactors = FALSE
    )
  }
  metric_parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, function(x) x$eval_metrics_reps))
  if (length(metric_parts)) {
    metrics_all <- do.call(rbind, metric_parts)
    if (!is.null(metrics_all) && nrow(metrics_all) > 0) {
      metric_summary <- tryCatch(
        stats::aggregate(
          value ~ trait + model + metric,
          data = metrics_all,
          FUN = function(v) mean(v, na.rm = TRUE)
        ),
        error = function(e) {
          data.frame(
            trait = character(),
            model = character(),
            metric = character(),
            value = numeric(),
            stringsAsFactors = FALSE
          )
        }
      )
    } else {
      metrics_all <- data.frame(
        trait = character(),
        model = character(),
        rep = integer(),
        fold = integer(),
        metric = character(),
        value = numeric(),
        stringsAsFactors = FALSE
      )
      metric_summary <- data.frame(
        trait = character(),
        model = character(),
        metric = character(),
        value = numeric(),
        stringsAsFactors = FALSE
      )
    }
  } else {
    metric_summary <- data.frame(
      trait = character(),
      model = character(),
      metric = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  }

  rank_risk_summary <- tryCatch(
    gp_multitrait_dl_rank_risk_summary(pred_all),
    error = function(e) {
      list(
        combined_risk = data.frame(
          GID = character(),
          Trait = character(),
          fold = integer(),
          rep = integer(),
          model = character(),
          rank_shift = numeric(),
          rank_instability_risk = numeric(),
          rank_instability_risk_remarks = character(),
          stringsAsFactors = FALSE
        ),
        aggregated_risk = data.frame(
          trait = character(),
          model = character(),
          mean_rank_shift = numeric(),
          mean_rank_instability_risk = numeric(),
          low_risk_percentage = numeric(),
          moderate_risk_percentage = numeric(),
          high_risk_percentage = numeric(),
          stringsAsFactors = FALSE
        )
      )
    }
  )

  list(
    multitrait_metric_summary = metric_summary,
    multitrait_prediction_counts = counts,
    multitrait_cv_plots = gp_multitrait_dl_cv_plot(pred_all),
    multitrait_uncertainty_summaries = gp_multitrait_dl_uncertainty_summary(pred_all),
    multitrait_rank_risk_summaries = rank_risk_summary,
    multitrait_rank_risk_plot = gp_multitrait_dl_risk_plot(rank_risk_summary)
  )
}

gp_multitrait_ml_gaussian_cv <- function(model_type = c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
                                         pheno_object,
                                         response,
                                         geno_omic_object,
                                         gen_name,
                                         cross_validation_meth = "K-Folds",
                                         nfolds = 5L,
                                         replication = 1L,
                                         eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                         ntree = 500,
                                         mtry = NULL,
                                         maxnodes = NULL,
                                         nodesize = NULL,
                                         lambda_rr = NULL,
                                         ncomp = 3L,
                                         svm_type = "eps-regression",
                                         svm_kernel = "Gaussian",
                                         sigma_value = NULL,
                                         C_value = 1,
                                         degree_value = 3,
                                         scale_value = 1,
                                         offset_value = 0,
                                         gamma_value = NULL,
                                         iteration = 100,
                                         learning_rate = 0.01,
                                         max_depth_xgb = 6,
                                         subsample_xgb = 0.7,
                                         xgb_booster = "gbtree",
                                         xgb_alpha = 0.001,
                                         xgb_lambda = 1.0,
                                         xgb_gamma = 0.01,
                                         min_child_weight = 1,
                                         xgb_nthread = 1L,
                                         colsample_bytree = 0.7,
                                         catboost_iterations = 500,
                                         catboost_depth = 6,
                                         catboost_learning_rate = 0.03,
                                         catboost_l2_leaf_reg = 3,
                                         catboost_thread_count = 1L,
                                         rf_n_jobs = 1L,
                                         random_seed = 123L,
                                         feature_score_metadata = NULL,
                                         feature_k = NULL,
                                         feature_scoring_cv = "fixed",
                                         feature_scoring_model = "Ridge_Regression",
                                         feature_source_map = NULL,
                                         feature_ridge_lambda = 1,
                                         feature_bayes_nIter = 1500L,
                                         feature_bayes_burnIn = 500L,
                                         feature_bayes_thin = 5L) {
  model_type <- match.arg(model_type)
  cv_token <- normalize_cv_token(cross_validation_meth)
  if (!cv_token %in% c("k_folds", "repeated_k_folds")) {
    stop("Multi-trait ML CV currently supports K-Folds only.", call. = FALSE)
  }

  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  ids <- intersect(as.character(ph[[gen_name]]), rownames(geno_omic_object))
  if (!length(ids)) {
    stop("No overlapping genotype IDs between phenotypes and processed features.", call. = FALSE)
  }
  ph <- ph[match(ids, as.character(ph[[gen_name]])), c(gen_name, response), drop = FALSE]
  x_base <- as.matrix(geno_omic_object[ids, , drop = FALSE])
  storage.mode(x_base) <- "double"
  feature_active <- !is.null(feature_score_metadata) && !is.null(feature_k)
  feature_scoring_cv <- gp_feature_normalize_cv_policy(feature_scoring_cv, allow_both = FALSE)

  cv_results <- vector("list", length = replication)
  for (rep_idx in seq_len(replication)) {
    pred_blocks <- list()
    target_blocks <- list()
    feature_records <- list()
    block_idx <- 1L
    target_block_idx <- 1L
    pending <- list()   # one fit job per trait x fold, run together below
    for (tr in response) {
      y <- as.numeric(ph[[tr]])
      obs_idx <- which(!is.na(y))
      if (!length(obs_idx)) next
      folds_rel <- gp_tuning_folds(
        y = seq_along(obs_idx),
        response_family = "gaussian",
        nfolds = min(nfolds, length(obs_idx)),
        random_state = as.integer((random_seed %||% 123L) + rep_idx - 1L)
      )
      for (fold_idx in seq_along(folds_rel)) {
        test_idx <- obs_idx[sort(unique(as.integer(folds_rel[[fold_idx]])))]
        train_idx <- setdiff(obs_idx, test_idx)
        x_feature <- x_base
        if (isTRUE(feature_active)) {
          feature_view <- gp_feature_multitrait_view(
            predictor_data = x_base,
            pheno_data = ph,
            response = tr,
            gen_name = gen_name,
            k = feature_k,
            selection_mode = feature_scoring_cv,
            feature_score_metadata = feature_score_metadata,
            test_rows = test_idx,
            scoring_model = feature_scoring_model,
            seed = as.integer((random_seed %||% 123L) + rep_idx * 10000L + fold_idx),
            source_block = feature_source_map,
            response_family = "gaussian",
            replication = rep_idx,
            fold = fold_idx,
            ridge_lambda = feature_ridge_lambda,
            bayes_nIter = feature_bayes_nIter,
            bayes_burnIn = feature_bayes_burnIn,
            bayes_thin = feature_bayes_thin,
            ntree = ntree,
            mtry = mtry,
            nodesize = nodesize,
            rf_n_jobs = rf_n_jobs
          )
          x_feature <- feature_view$predictor_data
          feature_records[[length(feature_records) + 1L]] <- feature_view$metadata
        }
        aux_traits <- setdiff(response, tr)
        aux_df <- ph[, aux_traits, drop = FALSE]
        aux_ready <- gp_multitrait_ml_prepare_aux_features(aux_df, train_idx = train_idx)
        x_all <- if (is.null(aux_ready)) x_feature else as.matrix(cbind(x_feature, aux_ready))
        storage.mode(x_all) <- "double"
        pending[[length(pending) + 1L]] <- c(gp_multitrait_ml_fit_trait(
          model_type = model_type,
          X_train = x_all[train_idx, , drop = FALSE],
          y_train = y[train_idx],
          X_test = x_all,
          return_job = TRUE,
          ntree = ntree,
          mtry = mtry,
          maxnodes = maxnodes,
          nodesize = nodesize,
          lambda_rr = lambda_rr,
          ncomp = ncomp,
          svm_type = svm_type,
          svm_kernel = svm_kernel,
          sigma_value = sigma_value,
          C_value = C_value,
          degree_value = degree_value,
          scale_value = scale_value,
          offset_value = offset_value,
          gamma_value = gamma_value,
          iteration = iteration,
          learning_rate = learning_rate,
          max_depth_xgb = max_depth_xgb,
          subsample_xgb = subsample_xgb,
          xgb_booster = xgb_booster,
          xgb_alpha = xgb_alpha,
          xgb_lambda = xgb_lambda,
          xgb_gamma = xgb_gamma,
          min_child_weight = min_child_weight,
          xgb_nthread = xgb_nthread,
          colsample_bytree = colsample_bytree,
          catboost_iterations = catboost_iterations,
          catboost_depth = catboost_depth,
          catboost_learning_rate = catboost_learning_rate,
          catboost_l2_leaf_reg = catboost_l2_leaf_reg,
          catboost_thread_count = catboost_thread_count,
          rf_n_jobs = rf_n_jobs
        ), list(tr = tr, fold_idx = fold_idx, test_idx = test_idx, y = y))
        # Keep the queued job light: rebuild its matrices from the shared
        # marker matrix when the batch runner reaches it.
        last <- length(pending)
        pending[[last]]$job$X_train <- NULL
        pending[[last]]$job$X_test <- NULL
        pending[[last]]$job$build <- local({
          xf <- x_feature; ar <- aux_ready; ti <- train_idx
          function() {
            xa <- if (is.null(ar)) xf else as.matrix(cbind(xf, ar))
            storage.mode(xa) <- "double"
            list(X_train = xa[ti, , drop = FALSE], X_test = xa)
          }
        })
      }
    }
    # All fold fits of this replication in one Python process (one process
    # per fit cost ~2 s each: 135 s for 2 traits x 5 folds x 3 replications).
    batch_results <- gp_py_ml_fit_predict_batch(lapply(pending, `[[`, "job"))
    for (k in seq_along(pending)) {
      tr <- pending[[k]]$tr
      fold_idx <- pending[[k]]$fold_idx
      test_idx <- pending[[k]]$test_idx
      y <- pending[[k]]$y
      {
        pred_all_targets <- gp_multitrait_ml_job_predictions(batch_results[[k]], pending[[k]]$scaler)
        pred_vals <- pred_all_targets[test_idx]
        pred_blocks[[block_idx]] <- data.frame(
          stringsAsFactors = FALSE,
          id_value = ids[test_idx],
          Trait = tr,
          model = model_type,
          Observed_value = y[test_idx],
          Predicted_value = pred_vals,
          Train_Test_Label = "Test",
          fold = fold_idx,
          rep = rep_idx,
          cv_role = "test"
        )
        names(pred_blocks[[block_idx]])[1] <- gen_name
        if (isTRUE(feature_active)) {
          pred_blocks[[block_idx]][["feature_k"]] <- as.integer(feature_k)
          pred_blocks[[block_idx]][["feature_scoring_model"]] <- feature_scoring_model
          pred_blocks[[block_idx]][["feature_scoring_cv"]] <- feature_scoring_cv
        }
        block_idx <- block_idx + 1L

        target_blocks[[target_block_idx]] <- data.frame(
          stringsAsFactors = FALSE,
          id_value = ids,
          Trait = tr,
          model = model_type,
          Observed_value = y,
          Predicted_value = pred_all_targets,
          Train_Test_Label = ifelse(is.na(y), "Test", "Train"),
          fold = fold_idx,
          rep = rep_idx,
          cv_role = "target_resampling"
        )
        names(target_blocks[[target_block_idx]])[1] <- gen_name
        target_block_idx <- target_block_idx + 1L
      }
    }
    pred_long <- if (length(pred_blocks)) do.call(rbind, pred_blocks) else data.frame()
    target_pred_long <- if (length(target_blocks)) do.call(rbind, target_blocks) else data.frame()
    metric_table <- gp_multitrait_ml_metric_table(
      pred_long = pred_long,
      eval_metrics = eval_metrics,
      model_type = model_type,
      rep = rep_idx
    )
    if (isTRUE(feature_active) && nrow(metric_table)) {
      metric_table[["feature_k"]] <- as.integer(feature_k)
      metric_table[["feature_scoring_model"]] <- feature_scoring_model
      metric_table[["feature_scoring_cv"]] <- feature_scoring_cv
    }
    cv_results[[rep_idx]] <- list(
      trait = paste(response, collapse = ","),
      model = model_type,
      rep = rep_idx,
      eval_metrics_reps = metric_table,
      ypred_cv_Reps_all = pred_long,
      target_prediction_Reps_all = target_pred_long,
      yprob_cv_Reps_all = NULL,
      feature_selection_metadata = gp_feature_bind_metadata(feature_records)
    )
  }

  processed <- gp_multitrait_ml_cv_process(cv_results)
  processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(cv_results, `[[`, "feature_selection_metadata")
  )
  list(cv_results = cv_results, cv_results_processed = processed)
}
