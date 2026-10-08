gp_hybrid_ml_supported_models <- function() {
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

gp_hybrid_bind_feature_blocks <- function(train_block, test_block = NULL) {
  train_mat <- if (is.null(train_block)) NULL else as.matrix(train_block)
  test_mat <- if (is.null(test_block)) NULL else as.matrix(test_block)

  if (is.null(train_mat) && is.null(test_mat)) {
    return(NULL)
  }
  if (is.null(train_mat)) {
    return(test_mat)
  }
  if (is.null(test_mat) || !nrow(test_mat)) {
    return(train_mat)
  }

  if (!identical(colnames(train_mat), colnames(test_mat))) {
    stop("Hybrid feature train/test blocks must have identical columns.", call. = FALSE)
  }

  out <- rbind(train_mat, test_mat)
  out <- out[!duplicated(rownames(out)), , drop = FALSE]
  storage.mode(out) <- "double"
  out
}

gp_hybrid_ml_summary_statistics <- function(pred_df,
                                            response,
                                            gen_name,
                                            female_parent,
                                            male_parent,
                                            model_type) {
  data.frame(
    stat = c(
      "mode",
      "response_family",
      "model_type",
      "scope",
      "response",
      "observed_hybrids",
      "predicted_hybrids",
      "n_female_parents",
      "n_male_parents",
      "female_parent_column",
      "male_parent_column",
      "feature_basis",
      "predictive_components"
    ),
    summary = c(
      "hybrid_ml",
      "gaussian",
      model_type,
      "hybrid gaussian true prediction from hybrid-level features with parent-aware CV scenarios",
      response,
      sum(pred_df$Train_Test_Label == "Train", na.rm = TRUE),
      sum(pred_df$Train_Test_Label == "Test", na.rm = TRUE),
      length(unique(as.character(pred_df[[female_parent]]))),
      length(unique(as.character(pred_df[[male_parent]]))),
      female_parent,
      male_parent,
      "hybrid_genotype_rows",
      "intercept;female_additive;male_additive;interaction"
    ),
    stringsAsFactors = FALSE
  )
}

gp_hybrid_ml_diagnostic_plot <- function(pred_df, model_label) {
  obs <- pred_df[!is.na(pred_df$Observed_value), , drop = FALSE]
  if (!nrow(obs)) {
    return(NULL)
  }
  ggplot2::ggplot(
    obs,
    ggplot2::aes(x = Observed_value, y = Predicted_value, color = Train_Test_Label)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Hybrid Gaussian", model_label, "observed vs predicted"),
      subtitle = "Hybrid ML uses hybrid-level genotype features; parent columns drive hybrid CV scenarios",
      x = "Observed value",
      y = "Predicted value"
    )
}

gp_hybrid_ml_component_profiles <- function(ph,
                                            response,
                                            female_parent,
                                            male_parent,
                                            shrinkage_k = 2) {
  y <- as.numeric(ph[[response]])
  keep <- !is.na(y)
  if (!any(keep)) {
    return(list(
      intercept = 0,
      female_delta = numeric(),
      male_delta = numeric(),
      basis = "predictive_parent_mean_decomposition"
    ))
  }

  train_df <- data.frame(
    y = y[keep],
    female = as.character(ph[[female_parent]][keep]),
    male = as.character(ph[[male_parent]][keep]),
    stringsAsFactors = FALSE
  )
  intercept <- mean(train_df$y, na.rm = TRUE)

  female_mean <- tapply(train_df$y, train_df$female, mean, na.rm = TRUE)
  female_n <- table(train_df$female)
  female_shrink <- pmax(0, as.numeric(female_n) / (as.numeric(female_n) + shrinkage_k))
  names(female_shrink) <- names(female_n)
  female_delta <- (female_mean - intercept) * female_shrink[names(female_mean)]
  female_delta[!is.finite(female_delta)] <- 0

  male_mean <- tapply(train_df$y, train_df$male, mean, na.rm = TRUE)
  male_n <- table(train_df$male)
  male_shrink <- pmax(0, as.numeric(male_n) / (as.numeric(male_n) + shrinkage_k))
  names(male_shrink) <- names(male_n)
  male_delta <- (male_mean - intercept) * male_shrink[names(male_mean)]
  male_delta[!is.finite(male_delta)] <- 0

  list(
    intercept = intercept,
    female_delta = female_delta,
    male_delta = male_delta,
    basis = "predictive_parent_mean_decomposition"
  )
}

gp_hybrid_ml_apply_component_decomposition <- function(pred_df,
                                                       female_parent,
                                                       male_parent,
                                                       profiles) {
  if (is.null(pred_df) || !nrow(pred_df)) {
    return(pred_df)
  }
  female_levels <- as.character(pred_df[[female_parent]])
  male_levels <- as.character(pred_df[[male_parent]])
  female_add <- as.numeric(unname(profiles$female_delta[female_levels]))
  male_add <- as.numeric(unname(profiles$male_delta[male_levels]))
  female_add[is.na(female_add)] <- 0
  male_add[is.na(male_add)] <- 0
  baseline <- as.numeric(rep(as.numeric(profiles$intercept %||% 0), nrow(pred_df)))
  interaction <- as.numeric(pred_df$Predicted_value) - baseline - female_add - male_add
  interaction[!is.finite(interaction)] <- NA_real_

  pred_df$Prediction_intercept <- baseline
  pred_df$Female_additive_contribution <- female_add
  pred_df$Male_additive_contribution <- male_add
  pred_df$Hybrid_interaction_contribution <- interaction
  pred_df$Female_GCA_like_contribution <- female_add
  pred_df$Male_GCA_like_contribution <- male_add
  pred_df$SCA_like_contribution <- interaction
  pred_df$Predictive_decomposition_basis <- profiles$basis %||% "predictive_parent_mean_decomposition"
  pred_df$Predictive_decomposition_remarks <- paste(
    "Predictive decomposition only:",
    "baseline + female additive + male additive + interaction = predicted value"
  )
  pred_df
}

gp_hybrid_ml_align_inputs <- function(pheno_object,
                                      response,
                                      gen_name,
                                      female_parent,
                                      male_parent,
                                      geno_omic_object,
                                      heter_groups = NULL) {
  ph <- as.data.frame(pheno_object, stringsAsFactors = FALSE)
  required_cols <- c(gen_name, female_parent, male_parent, response)
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    required_cols <- c(required_cols, heter_groups)
  }
  missing_cols <- setdiff(required_cols, names(ph))
  if (length(missing_cols)) {
    stop(
      paste("Hybrid ML phenotype data is missing required columns:", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  ids <- as.character(ph[[gen_name]])
  keep <- ids %in% rownames(geno_omic_object)
  ids <- ids[keep]
  ph <- ph[keep, required_cols, drop = FALSE]
  if (!length(ids)) {
    stop("No overlapping hybrid IDs between phenotypes and processed features.", call. = FALSE)
  }
  x <- as.matrix(geno_omic_object[match(ids, rownames(geno_omic_object)), , drop = FALSE])
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    env_block <- gp_one_hot_env_block(ph[[heter_groups]])
    x <- cbind(x, as.matrix(env_block))
  }
  storage.mode(x) <- "double"
  list(ph = ph, ids = ids, x = x)
}

gp_hybrid_ml_fit_response <- function(model_type,
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
                                      rf_n_jobs = 1L) {
  zero_test <- is.null(dim(X_test)) || nrow(X_test) == 0L
  if (isTRUE(zero_test)) {
    X_test <- X_train[1, , drop = FALSE]
  }
  out <- gp_multitrait_ml_fit_trait(
    model_type = model_type,
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
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
  if (isTRUE(zero_test)) {
    out$test_pred <- numeric(0)
  }
  out
}

gp_hybrid_ml_metric_table <- function(pred_df,
                                      eval_metrics,
                                      model_type) {
  if (is.null(pred_df) || !nrow(pred_df)) {
    return(data.frame())
  }
  split_keys <- unique(pred_df[, c("cv_scenario", "fold", "rep"), drop = FALSE])
  out <- lapply(seq_len(nrow(split_keys)), function(i) {
    key <- split_keys[i, , drop = FALSE]
    idx <- pred_df$cv_scenario == key$cv_scenario &
      pred_df$fold == key$fold &
      pred_df$rep == key$rep
    sub <- pred_df[idx, , drop = FALSE]
    row <- data.frame(
      cv_scenario = key$cv_scenario,
      fold = key$fold,
      rep = key$rep,
      model = model_type,
      stringsAsFactors = FALSE
    )
    for (m in eval_metrics) {
      row[[m]] <- safe_metric_value(
        y_true = sub$Observed_value,
        y_pred = sub$Predicted_value,
        metric = m,
        response_family = "gaussian"
      )
    }
    row
  })
  do.call(rbind, out)
}

gp_hybrid_ml_cv_process <- function(pred_df,
                                    eval_df) {
  metric_summary <- if (is.null(eval_df) || !nrow(eval_df)) {
    data.frame()
  } else {
    numeric_cols <- names(eval_df)[vapply(eval_df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, "rep")
    stats::aggregate(
      eval_df[, numeric_cols, drop = FALSE],
      by = list(cv_scenario = eval_df$cv_scenario, model = eval_df$model),
      FUN = mean,
      na.rm = TRUE
    )
  }

  counts <- if (is.null(pred_df) || !nrow(pred_df)) {
    data.frame()
  } else {
    stats::aggregate(
      Predicted_value ~ cv_scenario + Train_Test_Label,
      data = pred_df,
      FUN = length
    )
  }

  list(
    hybrid_metric_summary = metric_summary,
    hybrid_prediction_counts = counts,
    hybrid_cv_predictions = pred_df,
    hybrid_cv_plot = gp_hybrid_asreml_cv_plot(pred_df)
  )
}

gp_hybrid_ml_gaussian_model <- function(model_type = c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
                                        pheno_object,
                                        response,
                                        gen_name,
                                        female_parent,
                                        male_parent,
                                        heter_groups = NULL,
                                        geno_omic_object,
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
                                         n_bootstrap = 30L,
                                         random_seed = 123L,
                                         confidence_level = 0.95,
                                         internal_cv_nfolds = 5L,
                                         scaling = TRUE,
                                         centering = FALSE) {
  model_type <- match.arg(model_type)
  aligned <- gp_hybrid_ml_align_inputs(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    geno_omic_object = geno_omic_object,
    heter_groups = heter_groups
  )
  ph <- aligned$ph
  ids <- aligned$ids
  x <- aligned$x
  y <- as.numeric(ph[[response]])
  train_idx <- which(!is.na(y))
  if (!length(train_idx)) {
    stop("Hybrid ML requires at least one observed response value.", call. = FALSE)
  }
  test_idx <- setdiff(seq_len(nrow(ph)), train_idx)
  internal_id <- sprintf("hybrid_row_%06d", seq_len(nrow(ph)))
  x_internal <- x
  rownames(x_internal) <- internal_id
  model_key <- gp_multitrait_ml_python_model_key(model_type)
  model_params <- switch(
    model_type,
    "RandomForest" = list(
      ntree = ntree,
      mtry = mtry %||% max(1L, floor(sqrt(ncol(x)))),
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
    "PartialLeastSquare" = list(
      ncomp = max(1L, min(as.integer(ncomp %||% 3L), ncol(x), length(train_idx) - 1L))
    ),
    list()
  )
  model_pheno <- data.frame(
    .HybridRowID = internal_id[train_idx],
    .HybridResponse = y[train_idx],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  base_model <- gp_python_tabular_model(
    pheno_object = model_pheno,
    geno_omic_object = x_internal[train_idx, , drop = FALSE],
    geno_omic_test_object = if (length(test_idx)) {
      x_internal[test_idx, , drop = FALSE]
    } else {
      NULL
    },
    response = ".HybridResponse",
    gen_name = ".HybridRowID",
    response_family = "gaussian",
    model_type = model_key,
    model_label = model_type,
    model_params = model_params,
    scaling = scaling,
    centering = centering,
    tune_folds = internal_cv_nfolds,
    n_bootstrap = n_bootstrap,
    system_database = TRUE
  )
  base_pred <- as.data.frame(
    base_model[["predicted_values"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  base_match <- match(internal_id, as.character(base_pred[[".HybridRowID"]]))
  if (anyNA(base_match)) {
    stop("Hybrid ML predictions did not align with the hybrid rows.", call. = FALSE)
  }
  base_pred <- base_pred[base_match, , drop = FALSE]

  pred_df <- data.frame(
    stringsAsFactors = FALSE,
    id_value = ids,
    Female = as.character(ph[[female_parent]]),
    Male = as.character(ph[[male_parent]]),
    Observed_value = y,
    Predicted_value = suppressWarnings(as.numeric(base_pred[["Predicted_value"]])),
    Train_Test_Label = ifelse(is.na(y), "Test", "Train")
  )
  names(pred_df)[1] <- gen_name
  names(pred_df)[2] <- female_parent
  names(pred_df)[3] <- male_parent
  if (!is.null(heter_groups) && nzchar(heter_groups)) {
    pred_df[[heter_groups]] <- ph[[heter_groups]]
  }
  uncertainty_cols <- setdiff(
    names(base_pred),
    c(".HybridRowID", "Predicted_value", "Train_Test_Label", "Observed_value")
  )
  for (nm in uncertainty_cols) pred_df[[nm]] <- base_pred[[nm]]
  pred_df <- gp_hybrid_ml_apply_component_decomposition(
    pred_df = pred_df,
    female_parent = female_parent,
    male_parent = male_parent,
    profiles = gp_hybrid_ml_component_profiles(
      ph = ph,
      response = response,
      female_parent = female_parent,
      male_parent = male_parent
    )
  )

  predictive_matrices <- if (!is.null(heter_groups) && nzchar(heter_groups)) {
    gp_ml_predictive_matrix_bundle(
      pred = pred_df,
      id_col = gen_name,
      dimension_col = heter_groups,
      error_pred = pred_df
    )
  } else {
    list()
  }
  c(list(
    predicted_values = pred_df,
    diagnostic_plots = gp_hybrid_ml_diagnostic_plot(pred_df, model_label = model_type),
    prediction_uncertainty_source = unique(pred_df[["Prediction_uncertainty_source"]]),
    predictive_covariance_basis = if (length(predictive_matrices)) {
      "descriptive predictions and in-sample prediction errors; not genetic or residual covariance"
    } else {
      NULL
    }
  ), predictive_matrices)
}

gp_hybrid_ml_gaussian_cv <- function(model_type = c("RandomForest", "Ridge_Regression", "PartialLeastSquare", "Lasso", "SupportVectorMachine", "Xgboost", "CatBoost"),
                                     pheno_object,
                                     response,
                                     gen_name,
                                     female_parent,
                                     male_parent,
                                     heter_groups = NULL,
                                     geno_omic_object,
                                     cross_validation_meth = "Hybrid_Known_Parents",
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
                                     nfolds = 5L,
                                     random_state = 123L,
                                     replication = 1L,
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
  aligned <- gp_hybrid_ml_align_inputs(
    pheno_object = pheno_object,
    response = response,
    gen_name = gen_name,
    female_parent = female_parent,
    male_parent = male_parent,
    geno_omic_object = geno_omic_object,
    heter_groups = heter_groups
  )
  ph <- aligned$ph
  ids <- aligned$ids
  x <- aligned$x
  feature_active <- !is.null(feature_score_metadata) && !is.null(feature_k)
  feature_scoring_cv <- gp_feature_normalize_cv_policy(feature_scoring_cv, allow_both = FALSE)

  cv_results <- vector("list", length = replication)
  for (rep_idx in seq_len(replication)) {
    scenario_defs <- gp_hybrid_build_cv_scenarios(
      pheno_object = ph,
      response = response,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      random_state = random_state,
      replication = rep_idx
    )
    pred_parts <- list()
    feature_records <- list()
    part_idx <- 1L
    for (sc in scenario_defs) {
      scenario_rows <- gp_hybrid_cv_scenario_rows(sc, ph, gen_name)
      train_idx <- scenario_rows$train
      test_idx <- scenario_rows$test
      if (!length(train_idx) || !length(test_idx)) {
        next
      }
      x_feature <- x
      if (isTRUE(feature_active)) {
        feature_view <- gp_feature_multitrait_view(
          predictor_data = x,
          pheno_data = ph,
          response = response,
          gen_name = gen_name,
          k = feature_k,
          selection_mode = feature_scoring_cv,
          feature_score_metadata = feature_score_metadata,
          test_rows = test_idx,
          scoring_model = feature_scoring_model,
          seed = as.integer((random_state %||% 123L) + rep_idx * 10000L + part_idx),
          source_block = feature_source_map,
          response_family = "gaussian",
          replication = rep_idx,
          fold = part_idx,
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
      y <- as.numeric(ph[[response]])
      train_idx <- intersect(train_idx, which(!is.na(y)))
      test_idx <- intersect(test_idx, which(!is.na(y)))
      if (!length(train_idx) || !length(test_idx)) {
        next
      }
      train_cols <- c(gen_name, female_parent, male_parent, response)
      if (!is.null(heter_groups) && nzchar(heter_groups)) {
        train_cols <- c(train_cols, heter_groups)
      }
      train_ph <- ph[c(train_idx, test_idx), train_cols, drop = FALSE]
      train_ph[[response]][seq.int(length(train_idx) + 1L, nrow(train_ph))] <- NA_real_
      profiles <- gp_hybrid_ml_component_profiles(
        ph = train_ph,
        response = response,
        female_parent = female_parent,
        male_parent = male_parent
      )
      fit_res <- gp_hybrid_ml_fit_response(
        model_type = model_type,
        X_train = x_feature[train_idx, , drop = FALSE],
        y_train = y[train_idx],
        X_test = x_feature[test_idx, , drop = FALSE],
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
      pred_parts[[part_idx]] <- data.frame(
        stringsAsFactors = FALSE,
        id_value = ids[test_idx],
        Female = as.character(ph[[female_parent]][test_idx]),
        Male = as.character(ph[[male_parent]][test_idx]),
        Observed_value = y[test_idx],
        Predicted_value = fit_res$test_pred,
        Train_Test_Label = "Test",
        cv_scenario = as.character(sc$scenario),
        fold = as.character(sc$fold),
        rep = rep_idx,
        cv_role = "test"
      )
      names(pred_parts[[part_idx]])[1] <- gen_name
      names(pred_parts[[part_idx]])[2] <- female_parent
      names(pred_parts[[part_idx]])[3] <- male_parent
      if (!is.null(heter_groups) && nzchar(heter_groups)) {
        pred_parts[[part_idx]][[heter_groups]] <- ph[[heter_groups]][test_idx]
      }
      if (isTRUE(feature_active)) {
        pred_parts[[part_idx]][["feature_k"]] <- as.integer(feature_k)
        pred_parts[[part_idx]][["feature_scoring_model"]] <- feature_scoring_model
        pred_parts[[part_idx]][["feature_scoring_cv"]] <- feature_scoring_cv
      }
      pred_parts[[part_idx]] <- gp_hybrid_ml_apply_component_decomposition(
        pred_df = pred_parts[[part_idx]],
        female_parent = female_parent,
        male_parent = male_parent,
        profiles = profiles
      )
      part_idx <- part_idx + 1L
    }
    pred_df <- if (length(pred_parts)) do.call(rbind, pred_parts) else data.frame()
    metric_table <- gp_hybrid_ml_metric_table(
      pred_df = pred_df,
      eval_metrics = eval_metrics,
      model_type = model_type
    )
    if (isTRUE(feature_active) && nrow(metric_table)) {
      metric_table[["feature_k"]] <- as.integer(feature_k)
      metric_table[["feature_scoring_model"]] <- feature_scoring_model
      metric_table[["feature_scoring_cv"]] <- feature_scoring_cv
    }
    cv_results[[rep_idx]] <- list(
      trait = response,
      model = model_type,
      rep = rep_idx,
      eval_metrics_reps = metric_table,
      ypred_cv_Reps_all = pred_df,
      yprob_cv_Reps_all = NULL,
      feature_selection_metadata = gp_feature_bind_metadata(feature_records)
    )
  }

  pred_all <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, `[[`, "ypred_cv_Reps_all"))
  pred_all <- if (length(pred_all)) do.call(rbind, pred_all) else data.frame()
  eval_all <- Filter(function(x) !is.null(x) && nrow(x) > 0, lapply(cv_results, `[[`, "eval_metrics_reps"))
  eval_all <- if (length(eval_all)) do.call(rbind, eval_all) else data.frame()

  processed <- gp_hybrid_ml_cv_process(pred_all, eval_all)
  processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(cv_results, `[[`, "feature_selection_metadata")
  )
  list(cv_results = cv_results, cv_results_processed = processed)
}

# Hybrid ML/DL features are hybrid-level genotypes. When the inbred parents'
# genotypes are supplied, apply the heterozygosity filter to the parents
# (where it is a meaningful QC rule for DH/RIL lines) and drop the failing
# markers from the hybrid features. Without hybrid-level geno_data, the hybrid
# features are built as the expected F1 dosage (female + male) / 2.
gp_hybrid_ml_parent_geno_prepare <- function(pheno_data,
                                             geno_data = NULL,
                                             female_geno_data = NULL,
                                             male_geno_data = NULL,
                                             gen_name,
                                             female_parent,
                                             male_parent,
                                             het_threshold = 0.1,
                                             ploidy = "auto",
                                             message = TRUE) {
  if (is.null(female_geno_data) && is.null(male_geno_data)) {
    return(list(geno_data = geno_data, summary = NULL))
  }
  female_geno_data <- female_geno_data %||% male_geno_data
  male_geno_data <- male_geno_data %||% female_geno_data
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  id_cols <- c(gen_name, female_parent, male_parent)
  if (length(id_cols) != 3L || !all(id_cols %in% names(ph))) {
    stop(
      "Hybrid ML/DL with parent genotypes needs gen_name, female_parent and ",
      "male_parent columns in pheno_data.",
      call. = FALSE
    )
  }
  as_parent_matrix <- function(x, arg) {
    x <- as.matrix(x)
    if (is.null(rownames(x)) || is.null(colnames(x))) {
      stop(arg, " must have parent IDs as row names and marker names as column names.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    x
  }
  fg <- as_parent_matrix(female_geno_data, "female_geno_data")
  mg <- as_parent_matrix(male_geno_data, "male_geno_data")
  females <- unique(as.character(stats::na.omit(ph[[female_parent]])))
  males <- unique(as.character(stats::na.omit(ph[[male_parent]])))
  fg <- fg[rownames(fg) %in% females, , drop = FALSE]
  mg <- mg[rownames(mg) %in% males, , drop = FALSE]
  if (!nrow(fg) || !nrow(mg)) {
    stop(
      "female_geno_data/male_geno_data row names do not match the Female/Male ",
      "parent IDs in pheno_data.",
      call. = FALSE
    )
  }
  markers <- intersect(colnames(fg), colnames(mg))
  if (!length(markers)) {
    stop("female_geno_data and male_geno_data share no marker columns.", call. = FALSE)
  }
  fg <- fg[, markers, drop = FALSE]
  mg <- mg[, markers, drop = FALSE]

  het_fail <- rep(FALSE, length(markers))
  names(het_fail) <- markers
  if (!is.null(het_threshold)) {
    het_f <- as.numeric(geno_qc_metrics(fg, ploidy = ploidy)$heterozygosity)
    het_m <- as.numeric(geno_qc_metrics(mg, ploidy = ploidy)$heterozygosity)
    het_fail[] <- (!is.na(het_f) & het_f > het_threshold) |
      (!is.na(het_m) & het_m > het_threshold)
  }
  failing <- markers[het_fail]
  kept_parent_markers <- markers[!het_fail]

  if (!is.null(geno_data)) {
    hyb <- as.matrix(geno_data)
    not_in_parents <- setdiff(colnames(hyb), markers)
    keep <- setdiff(colnames(hyb), failing)
    if (!length(keep)) {
      stop(
        "Every hybrid marker failed the parent heterozygosity filter ",
        "(het_threshold = ", het_threshold, ").",
        call. = FALSE
      )
    }
    out <- hyb[, keep, drop = FALSE]
    source <- "hybrid_geno_data"
  } else {
    if (!length(kept_parent_markers)) {
      stop(
        "Every marker failed the parent heterozygosity filter ",
        "(het_threshold = ", het_threshold, ").",
        call. = FALSE
      )
    }
    trio <- unique(ph[!is.na(ph[[gen_name]]), id_cols, drop = FALSE])
    trio[] <- lapply(trio, as.character)
    dup <- unique(trio[[gen_name]][duplicated(trio[[gen_name]])])
    if (length(dup)) {
      stop(
        "Each hybrid must map to one female/male parent pair; conflicting ",
        "parents for: ", paste(utils::head(dup, 10), collapse = ", "),
        call. = FALSE
      )
    }
    missing_f <- setdiff(trio[[female_parent]], rownames(fg))
    missing_m <- setdiff(trio[[male_parent]], rownames(mg))
    if (length(missing_f) || length(missing_m)) {
      stop(
        "Parent genotypes are missing for ",
        paste(c(
          if (length(missing_f)) paste("female:", paste(utils::head(missing_f, 10), collapse = ", ")),
          if (length(missing_m)) paste("male:", paste(utils::head(missing_m, 10), collapse = ", "))
        ), collapse = "; "),
        call. = FALSE
      )
    }
    out <- (fg[trio[[female_parent]], kept_parent_markers, drop = FALSE] +
      mg[trio[[male_parent]], kept_parent_markers, drop = FALSE]) / 2
    rownames(out) <- trio[[gen_name]]
    not_in_parents <- character()
    source <- "expected_F1_from_parents"
  }

  summary <- data.frame(
    item = c(
      "hybrid_feature_source", "parent_het_threshold", "female_parents_checked",
      "male_parents_checked", "parent_markers_checked",
      "markers_removed_parent_heterozygosity",
      "hybrid_markers_not_in_parent_genotypes", "hybrid_markers_kept"
    ),
    value = c(
      source, if (is.null(het_threshold)) "off" else as.character(het_threshold),
      nrow(fg), nrow(mg), length(markers), length(failing),
      length(not_in_parents), ncol(out)
    ),
    stringsAsFactors = FALSE
  )
  if (isTRUE(message)) {
    base::message(
      "Hybrid ML/DL parent QC: ", length(failing), " of ", length(markers),
      " markers removed for parent heterozygosity > ",
      if (is.null(het_threshold)) "off" else het_threshold,
      "; ", ncol(out), " hybrid features kept (", source, ")."
    )
  }
  list(geno_data = out, summary = summary)
}