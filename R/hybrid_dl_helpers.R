gp_hybrid_dl_supported_models <- function() c("mlp", "ft_transformer")

gp_hybrid_dl_summary_statistics <- function(pred_df,
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
      "hybrid_dl",
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

gp_hybrid_dl_fit_response <- function(model_type,
                                      X_train,
                                      y_train,
                                      X_test,
                                      scaling = TRUE,
                                      centering = TRUE,
                                      compile_model = TRUE,
                                      deterministic = FALSE,
                                      random_seed = 123L,
                                      device = NULL,
                                      use_amp = TRUE,
                                      batch_norm = TRUE,
                                      validation_split = 0.2,
                                      epochs = 20L,
                                      batch_size = 64L,
                                      mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                      mlp_learning_rate = 1e-3,
                                      ft_d_model = 128L,
                                      ft_heads = 8L,
                                      ft_layers = 3L,
                                      ft_ff_mult = 4L,
                                      ft_dropout = 0.1,
                                      ft_token_dropout = 0,
                                       ft_use_cls = TRUE,
                                       dropout = 0.2,
                                       l2_weight_decay = 1e-4,
                                       optimizer_name = "adam",
                                       max_grad_norm = 1) {
  zero_test <- is.null(dim(X_test)) || nrow(X_test) == 0L
  if (isTRUE(zero_test)) {
    X_test <- X_train[1, , drop = FALSE]
  }

  x_prep <- NULL
  X_tr <- as.matrix(X_train)
  X_te <- as.matrix(X_test)
  storage.mode(X_tr) <- "double"
  storage.mode(X_te) <- "double"
  if (ncol(X_tr) > 0L) {
    keep_cols <- vapply(
      seq_len(ncol(X_tr)),
      function(j) {
        s <- stats::sd(X_tr[, j], na.rm = TRUE)
        is.finite(s) && s > 0
      },
      logical(1)
    )
    if (any(keep_cols)) {
      X_tr <- X_tr[, keep_cols, drop = FALSE]
      X_te <- X_te[, keep_cols, drop = FALSE]
    }
  }
  if (isTRUE(scaling) || isTRUE(centering)) {
    x_prep <- gp_fast_center_scale_fit(
      X_tr,
      center = isTRUE(centering),
      scale = isTRUE(scaling)
    )
    X_tr <- as.matrix(x_prep$data)
    X_te <- as.matrix(gp_fast_center_scale_apply(X_te, x_prep))
  }

  y_scaler <- gp_fast_y_scale_fit(y_train)
  y_tr <- y_scaler$scaled

  args_common <- compact(list(
    optimizer_name = optimizer_name,
    use_amp = use_amp,
    max_grad_norm = max_grad_norm,
    compile_model = compile_model,
    deterministic = deterministic,
    random_seed = random_seed,
    device = device,
    batch_norm = batch_norm,
    validation_split = validation_split
  ))

  args_model <- switch(
    tolower(model_type),
    "mlp" = compact(list(
      mlp_neurons_per_layer = mlp_neurons_per_layer,
      mlp_final_attention = FALSE,
      mlp_attention_across_layers = FALSE,
      mlp_learning_rate = mlp_learning_rate,
      mlp_epochs = epochs,
      mlp_batch_size = batch_size,
      mlp_dropout_rate = dropout,
      mlp_l2_weight_decay = l2_weight_decay
    )),
    "ft_transformer" = compact(list(
      ft_d_model = ft_d_model,
      ft_heads = ft_heads,
      ft_layers = ft_layers,
      ft_ff_mult = ft_ff_mult,
      ft_dropout = ft_dropout,
      ft_token_dropout = ft_token_dropout,
      ft_use_cls = ft_use_cls,
      ft_epochs = epochs,
      ft_batch_size = batch_size,
      ft_l2_weight_decay = l2_weight_decay
    )),
    stop("Unsupported hybrid DL model: ", model_type, call. = FALSE)
  )

  X_pred <- if (isTRUE(zero_test)) X_tr else rbind(X_tr, X_te)
  pred_scaled_all <- gp_dl_bridge_fit_predict(
    model_type = tolower(model_type),
    X_train = X_tr,
    y_train = y_tr,
    X_test = X_pred,
    response_family = "gaussian",
    dl_args = c(args_model, args_common)
  )
  pred_scaled_all <- as.numeric(pred_scaled_all)
  pred_train <- revert_scaling(pred_scaled_all[seq_len(nrow(X_tr))], y_scaler)
  pred_test <- if (isTRUE(zero_test)) {
    numeric(0)
  } else {
    revert_scaling(pred_scaled_all[nrow(X_tr) + seq_len(nrow(X_te))], y_scaler)
  }

  list(
    train_pred = as.numeric(pred_train),
    test_pred = as.numeric(pred_test),
    trained_model = NULL
  )
}

gp_hybrid_dl_metric_table <- function(pred_df,
                                      eval_metrics,
                                      model_type) {
  gp_hybrid_ml_metric_table(
    pred_df = pred_df,
    eval_metrics = eval_metrics,
    model_type = model_type
  )
}

gp_hybrid_dl_cv_process <- function(pred_df,
                                    eval_df) {
  gp_hybrid_ml_cv_process(pred_df = pred_df, eval_df = eval_df)
}

gp_hybrid_dl_gaussian_model <- function(model_type = c("mlp", "ft_transformer"),
                                        pheno_object,
                                        response,
                                        gen_name,
                                        female_parent,
                                        male_parent,
                                        heter_groups = NULL,
                                        geno_omic_object,
                                        scaling = TRUE,
                                        centering = TRUE,
                                        compile_model = TRUE,
                                        deterministic = FALSE,
                                        random_seed = 123L,
                                        device = NULL,
                                        use_amp = TRUE,
                                        batch_norm = TRUE,
                                        validation_split = 0.2,
                                        epochs = 20L,
                                        batch_size = 64L,
                                        mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                        mlp_learning_rate = 1e-3,
                                        ft_d_model = 128L,
                                        ft_heads = 8L,
                                        ft_layers = 3L,
                                        ft_ff_mult = 4L,
                                        ft_dropout = 0.1,
                                        ft_token_dropout = 0,
                                        ft_use_cls = TRUE,
                                         dropout = 0.2,
                                         l2_weight_decay = 1e-4,
                                         optimizer_name = "adam",
                                         max_grad_norm = 1,
                                          n_bootstrap = 30L,
                                          dl_internal_calibration = TRUE,
                                          confidence_level = 0.95,
                                         early_stop = TRUE) {
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
    stop("Hybrid DL requires at least one observed response value.", call. = FALSE)
  }
  test_idx <- setdiff(seq_len(nrow(ph)), train_idx)
  internal_id <- sprintf("hybrid_row_%06d", seq_len(nrow(ph)))
  x_internal <- x
  rownames(x_internal) <- internal_id
  model_pheno <- data.frame(
    .HybridRowID = internal_id[train_idx],
    .HybridResponse = y[train_idx],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  base_model <- deep_learning_model(
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
    model_type = model_type,
    message = FALSE,
    scaling = scaling,
    centering = centering,
    n_bootstrap = n_bootstrap,
    dl_internal_calibration = dl_internal_calibration,
    system_database = TRUE,
    optimizer_name = optimizer_name,
    use_amp = use_amp,
    max_grad_norm = max_grad_norm,
    early_stop = early_stop,
    compile_model = compile_model,
    deterministic = deterministic,
    random_seed = random_seed,
    device = device,
    batch_norm = batch_norm,
    validation_split = validation_split,
    epochs = epochs,
    batch_size = batch_size,
    mlp_neurons_per_layer = mlp_neurons_per_layer,
    mlp_learning_rate = mlp_learning_rate,
    ft_d_model = ft_d_model,
    ft_heads = ft_heads,
    ft_layers = ft_layers,
    ft_ff_mult = ft_ff_mult,
    ft_dropout = ft_dropout,
    ft_token_dropout = ft_token_dropout,
    ft_use_cls = ft_use_cls,
    dropout = dropout,
    l2_weight_decay = l2_weight_decay
  )
  base_pred <- as.data.frame(
    base_model[["predicted_values"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  base_match <- match(internal_id, as.character(base_pred[[".HybridRowID"]]))
  if (anyNA(base_match)) {
    stop("Hybrid DL predictions did not align with the hybrid rows.", call. = FALSE)
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
      error_pred = pred_df,
      random_seed = random_seed
    )
  } else {
    list()
  }
  c(list(
    predicted_values = pred_df,
    model_parameters = base_model[["model_parameters"]],
    trained_model = NULL,
    diagnostic_plots = gp_hybrid_ml_diagnostic_plot(pred_df, model_label = model_type),
    prediction_uncertainty_source = unique(pred_df[["Prediction_uncertainty_source"]]),
    predictive_covariance_basis = if (!dl_internal_calibration) {
      "descriptive predictions; heldout prediction errors unavailable; not genetic or residual covariance"
    } else if (length(predictive_matrices)) {
      "descriptive predictions and in-sample prediction errors; not genetic or residual covariance"
    } else {
      NULL
    }
  ), predictive_matrices)
}

gp_hybrid_dl_gaussian_cv <- function(model_type = c("mlp", "ft_transformer"),
                                     pheno_object,
                                     response,
                                     gen_name,
                                     female_parent,
                                     male_parent,
                                     heter_groups = NULL,
                                     geno_omic_object,
                                     cross_validation_meth = "Hybrid_Known_Parents",
                                     eval_metrics = c("root_mean_squared_error", "mean_absolute_error", "kendalls_tau"),
                                     nfolds = 5L,
                                     random_state = 123L,
                                     replication = 1L,
                                     scaling = TRUE,
                                     centering = TRUE,
                                     compile_model = TRUE,
                                     deterministic = FALSE,
                                     device = NULL,
                                     use_amp = TRUE,
                                     batch_norm = TRUE,
                                     validation_split = 0.2,
                                     epochs = 20L,
                                     batch_size = 64L,
                                     mlp_neurons_per_layer = as.integer(c(128L, 64L)),
                                     mlp_learning_rate = 1e-3,
                                     ft_d_model = 128L,
                                     ft_heads = 8L,
                                     ft_layers = 3L,
                                     ft_ff_mult = 4L,
                                     ft_dropout = 0.1,
                                     ft_token_dropout = 0,
                                     ft_use_cls = TRUE,
                                     dropout = 0.2,
                                     l2_weight_decay = 1e-4,
                                     optimizer_name = "adam",
                                     max_grad_norm = 1,
                                     feature_score_metadata = NULL,
                                     feature_k = NULL,
                                     feature_scoring_cv = "fixed",
                                     feature_scoring_model = "Ridge_Regression",
                                     feature_source_map = NULL,
                                     feature_ridge_lambda = 1,
                                     feature_bayes_nIter = 1500L,
                                     feature_bayes_burnIn = 500L,
                                     feature_bayes_thin = 5L,
                                     ntree = 500L,
                                     mtry = NULL,
                                     nodesize = NULL,
                                     rf_n_jobs = 1L) {
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
      if (!length(train_idx) || !length(test_idx)) next
      y <- as.numeric(ph[[response]])
      train_idx <- intersect(train_idx, which(!is.na(y)))
      test_idx <- intersect(test_idx, which(!is.na(y)))
      if (!length(train_idx) || !length(test_idx)) next
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
      fit_res <- gp_hybrid_dl_fit_response(
        model_type = model_type,
        X_train = x_feature[train_idx, , drop = FALSE],
        y_train = y[train_idx],
        X_test = x_feature[test_idx, , drop = FALSE],
        scaling = scaling,
        centering = centering,
        compile_model = compile_model,
        deterministic = deterministic,
        random_seed = as.integer(random_state + rep_idx - 1L),
        device = device,
        use_amp = use_amp,
        batch_norm = batch_norm,
        validation_split = validation_split,
        epochs = epochs,
        batch_size = batch_size,
        mlp_neurons_per_layer = mlp_neurons_per_layer,
        mlp_learning_rate = mlp_learning_rate,
        ft_d_model = ft_d_model,
        ft_heads = ft_heads,
        ft_layers = ft_layers,
        ft_ff_mult = ft_ff_mult,
        ft_dropout = ft_dropout,
        ft_token_dropout = ft_token_dropout,
        ft_use_cls = ft_use_cls,
        dropout = dropout,
        l2_weight_decay = l2_weight_decay,
        optimizer_name = optimizer_name,
        max_grad_norm = max_grad_norm
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
    metric_table <- gp_hybrid_dl_metric_table(
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

  processed <- gp_hybrid_dl_cv_process(pred_all, eval_all)
  processed[["feature_selection_metadata"]] <- gp_feature_bind_metadata(
    lapply(cv_results, `[[`, "feature_selection_metadata")
  )
  list(cv_results = cv_results, cv_results_processed = processed)
}
