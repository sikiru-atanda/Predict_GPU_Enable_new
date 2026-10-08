gp_route_cross_validation_specialized_models <- function(ctx) {
  gp_model_execute_import_context(ctx)
  if (isTRUE(multi_trait_gp)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Multi-trait GP CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) < 2L || !identical(response_family, "gaussian")) {
      stop("Multi-trait GP CV requires at least two gaussian response columns.", call. = FALSE)
    }
    model_canonical <- gp_canonicalize_supported_model_names(GS_model_cv)
    if (length(model_canonical) != 1L || !model_canonical %in% gp_multitrait_gp_supported_models()) {
      stop(
        "Multi-trait GP CV supports: ",
        paste(gp_display_supported_model_names(gp_multitrait_gp_supported_models()), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    mt_pheno <- pheno_clean[["pheno_clean_data"]]
    mt_gp_is_met <- gp_is_multi_environment_trait_panel(
      pheno_data = mt_pheno,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response = response,
      response_family = response_family
    )
    if (!is.null(feature_score_metadata) && !is.null(feature_k)) {
      stop(
        "Multi-trait GP CV with response-derived feature selection is not yet supported. Disable feature_scoring or supply an externally fixed kernel.",
        call. = FALSE
      )
    }
    mt_kernel_bank <- gmatrix_kernel_model_ready_list %||% list()
    mt_gmatrix <- mt_kernel_bank[["gmatrix_model_ready"]] %||% gmatrix
    if (is.null(mt_gmatrix)) {
      stop("Multi-trait GP CV requires one model-ready genomic relationship matrix.", call. = FALSE)
    }
    cv_pipeline <- gp_multitrait_gp_gaussian_cv(
      estimate_kernel_weights = isTRUE(gp_estimate_kernel_weights),
      force_prediction_se = isTRUE(gp_force_prediction_se),
      model_type = model_canonical,
      pheno_object = mt_pheno,
      response = response,
      gen_name = gen_name,
      gmatrix = mt_gmatrix,
      kernel_list = mt_kernel_bank,
      kernel_weights = lowrank_kernel_weights,
      heter_groups = if (isTRUE(mt_gp_is_met)) heter_groups else NULL,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      eval_metrics = eval_metrics,
      random_seed = random_seed %||% random_state %||% 123L,
      gp_factor_cache = gp_factor_cache,
      gp_fa_rank = gp_fa_rank %||% 1L,
      gp_return_se = gp_return_se,
      gp_full_vc = gp_full_vc,
      gp_return_trait_correlations = gp_return_trait_correlations,
      gp_iters = gp_iters,
      para_tunning = para_tunning,
      env_similarity = env_similarity,
      env_ids = env_ids,
      env_covariates = env_covariates,
      reaction_norm_feature_qc = reaction_norm_feature_qc,
      kenv_kernel = kenv_kernel,
      kenv_bandwidth = kenv_bandwidth,
      kenv_kernel_kwargs = kenv_kernel_kwargs
    )
    model_public <- gp_public_model_label(model_canonical)[1L]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = gp_detect_python(),
      preferred_python = gp_preferred_python(purpose = "gp"),
      context = if (isTRUE(mt_gp_is_met)) {
        "cross_validation_multi_trait_multi_environment_gp"
      } else {
        "cross_validation_multi_trait_gp"
      },
      extra_fields = list(
        task_unit = if (isTRUE(mt_gp_is_met)) "multi_trait_multi_environment_gp_cv" else "multi_trait_gp_cv",
        gp_backend_method = cv_pipeline$cv_results_processed$gp_route$backend_method,
        all_traits_masked_together = TRUE
      )
    )
    return(results_handling(
      GS_model = model_public,
      cv_results_raw = cv_pipeline[["cv_results"]],
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_pipeline[["cv_results_processed"]],
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_multi_trait_gp",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(multi_trait_bayes)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Multi-trait Bayesian CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) < 2L || !identical(response_family, "gaussian")) {
      stop("Multi-trait Bayesian CV requires at least two gaussian response columns.", call. = FALSE)
    }
    model_canonical <- gp_canonicalize_supported_model_names(GS_model_cv)
    if (length(model_canonical) != 1L || !model_canonical %in% gp_multitrait_bayes_supported_models()) {
      stop(
        "Multi-trait Bayesian CV supports: ",
        paste(gp_multitrait_bayes_supported_models(), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    mt_pheno <- pheno_clean[["pheno_clean_data"]]
    mt_bayes_is_met <- gp_is_multi_environment_trait_panel(
      pheno_data = mt_pheno,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response = response,
      response_family = response_family
    )
    if (!is.null(feature_score_metadata) && !is.null(feature_k)) {
      stop(
        "Multi-trait Bayesian CV with response-derived feature selection is not yet supported. Disable feature_scoring or supply an externally fixed kernel.",
        call. = FALSE
      )
    }
    mt_kernel_bank <- gmatrix_kernel_model_ready_list %||% list()
    kernels <- gp_collect_kernel_inputs(
      gmatrix = mt_kernel_bank[["gmatrix_model_ready"]] %||% gmatrix,
      omic1_kernel = mt_kernel_bank[["omic1_kernel_model_ready"]] %||% NULL,
      omic2_kernel = mt_kernel_bank[["omic2_kernel_model_ready"]] %||% NULL,
      omic3_kernel = mt_kernel_bank[["omic3_kernel_model_ready"]] %||% NULL,
      kernel_list = mt_kernel_bank
    )
    if (!length(kernels)) {
      stop("Multi-trait Bayesian CV requires at least one model-ready relationship/kernel matrix.", call. = FALSE)
    }
    bayes_para <- bayes_parameter_check(nIter = nIter, burnIn = burnIn, thin = thin)
    cv_pipeline <- gp_multitrait_bayes_gaussian_cv(
      model_type = model_canonical,
      pheno_object = mt_pheno,
      response = response,
      gen_name = gen_name,
      kernels = kernels,
      bayes_para = bayes_para,
      heter_groups = if (isTRUE(mt_bayes_is_met)) heter_groups else NULL,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      eval_metrics = eval_metrics,
      random_seed = random_seed %||% random_state %||% 123L,
      heter_resid = heter_resid,
      confidence_level = confidence_level %||% 0.95,
      CI_width_thresholds = CI_width_thresholds %||% c(0.33, 0.66),
      high_reliability_thres = high_reliability_thres %||% 0.9,
      low_reliability_thres = low_reliability_thres %||% 0.5
    )
    model_public <- gp_public_model_label(model_canonical)[1L]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = if (isTRUE(mt_bayes_is_met)) {
        "cross_validation_multi_trait_multi_environment_bayes"
      } else {
        "cross_validation_multi_trait_bayes"
      },
      extra_fields = list(
        task_unit = if (isTRUE(mt_bayes_is_met)) "multi_trait_multi_environment_bayes_cv" else "multi_trait_bayes_cv",
        bglr_function = "Multitrait",
        all_traits_masked_together = TRUE
      )
    )
    return(results_handling(
      GS_model = model_public,
      cv_results_raw = cv_pipeline[["cv_results"]],
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_pipeline[["cv_results_processed"]],
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_multi_trait_bayes",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(hybrid_dl)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Hybrid DL CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) != 1L) {
      stop("Hybrid DL currently supports one gaussian response at a time.", call. = FALSE)
    }
    if (!identical(response_family, "gaussian")) {
      stop("Hybrid DL currently supports gaussian traits only.", call. = FALSE)
    }
    if (length(GS_model_cv) != 1L || !GS_model_cv %in% gp_hybrid_dl_supported_models()) {
      stop(
        paste(
          "Hybrid DL CV currently supports:",
          paste(gp_hybrid_dl_supported_models(), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
      stop("Provide both female_parent and male_parent column names for hybrid DL CV.", call. = FALSE)
    }
    if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
      stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
    }
    if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
      stop("Hybrid DL CV requires processed geno/omic features.", call. = FALSE)
    }
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_k,
      feature_k_grid = feature_k_grid,
      feature_scoring_cv = feature_scoring_cv,
      context = "Hybrid DL CV"
    )
    cv_pipeline <- gp_hybrid_dl_gaussian_cv(
      model_type = GS_model_cv,
      pheno_object = pheno_data %||% pheno_clean[["pheno_clean_data"]],
      response = response[[1]],
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
      cross_validation_meth = cross_validation_meth,
      eval_metrics = eval_metrics,
      nfolds = nfolds,
      random_state = random_state,
      replication = replication,
      scaling = scaling,
      centering = centering,
      compile_model = compile_model,
      deterministic = deterministic,
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
      max_grad_norm = max_grad_norm,
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_cfg$k,
      feature_scoring_cv = feature_cfg$mode,
      feature_scoring_model = feature_scoring_model,
      feature_source_map = feature_source_map,
      feature_ridge_lambda = feature_ridge_lambda,
      feature_bayes_nIter = feature_bayes_nIter,
      feature_bayes_burnIn = feature_bayes_burnIn,
      feature_bayes_thin = feature_bayes_thin,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    )
    cv_results <- cv_pipeline[["cv_results"]]
    cv_results_processed <- cv_pipeline[["cv_results_processed"]]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = gp_detect_python(),
      preferred_python = gp_preferred_python(purpose = "dl"),
      context = "cross_validation_hybrid_dl"
    )
    return(results_handling(
      GS_model = GS_model_cv,
      cv_results_raw = cv_results,
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_results_processed,
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_hybrid_dl",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(hybrid_ml)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Hybrid ML CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) != 1L) {
      stop("Hybrid ML currently supports one gaussian response at a time.", call. = FALSE)
    }
    if (!identical(response_family, "gaussian")) {
      stop("Hybrid ML currently supports gaussian traits only.", call. = FALSE)
    }
    if (length(GS_model_cv) != 1L || !GS_model_cv %in% gp_hybrid_ml_supported_models()) {
      stop(
        paste(
          "Hybrid ML CV currently supports:",
          paste(gp_hybrid_ml_supported_models(), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
      stop("Provide both female_parent and male_parent column names for hybrid ML CV.", call. = FALSE)
    }
    if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
      stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
    }
    if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
      stop("Hybrid ML CV requires processed geno/omic features.", call. = FALSE)
    }
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_k,
      feature_k_grid = feature_k_grid,
      feature_scoring_cv = feature_scoring_cv,
      context = "Hybrid ML CV"
    )
    cv_pipeline <- gp_hybrid_ml_gaussian_cv(
      model_type = GS_model_cv,
      pheno_object = pheno_data %||% pheno_clean[["pheno_clean_data"]],
      response = response[[1]],
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      eval_metrics = eval_metrics,
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs,
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
      max_depth_xgb = max_depth,
      subsample_xgb = subsample,
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
      random_state = random_state,
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_cfg$k,
      feature_scoring_cv = feature_cfg$mode,
      feature_scoring_model = feature_scoring_model,
      feature_source_map = feature_source_map,
      feature_ridge_lambda = feature_ridge_lambda,
      feature_bayes_nIter = feature_bayes_nIter,
      feature_bayes_burnIn = feature_bayes_burnIn,
      feature_bayes_thin = feature_bayes_thin
    )
    cv_results <- cv_pipeline[["cv_results"]]
    cv_results_processed <- cv_pipeline[["cv_results_processed"]]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = "cross_validation_hybrid_ml"
    )
    return(results_handling(
      GS_model = GS_model_cv,
      cv_results_raw = cv_results,
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_results_processed,
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_hybrid_ml",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(multi_trait_ml)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Multi-trait ML CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) < 2L) {
      stop("Multi-trait ML requires at least two response columns.", call. = FALSE)
    }
    if (!identical(response_family, "gaussian")) {
      stop("Multi-trait ML currently supports gaussian traits only.", call. = FALSE)
    }
    if (length(GS_model_cv) != 1L || !GS_model_cv %in% gp_multitrait_ml_supported_models()) {
      stop(
        paste(
          "Multi-trait ML CV currently supports:",
          paste(gp_multitrait_ml_supported_models(), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
      stop("Multi-trait ML CV requires processed geno/omic features.", call. = FALSE)
    }
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_k,
      feature_k_grid = feature_k_grid,
      feature_scoring_cv = feature_scoring_cv,
      context = "Multi-trait ML CV"
    )
    cv_pipeline <- gp_multitrait_ml_gaussian_cv(
      model_type = GS_model_cv,
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = response,
      geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
      gen_name = gen_name,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      eval_metrics = eval_metrics,
      ntree = ntree,
      mtry = mtry,
      maxnodes = maxnodes,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs,
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
      max_depth_xgb = max_depth,
      subsample_xgb = subsample,
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
      random_seed = random_seed,
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_cfg$k,
      feature_scoring_cv = feature_cfg$mode,
      feature_scoring_model = feature_scoring_model,
      feature_source_map = feature_source_map,
      feature_ridge_lambda = feature_ridge_lambda,
      feature_bayes_nIter = feature_bayes_nIter,
      feature_bayes_burnIn = feature_bayes_burnIn,
      feature_bayes_thin = feature_bayes_thin
    )
    cv_results <- cv_pipeline[["cv_results"]]
    cv_results_processed <- cv_pipeline[["cv_results_processed"]]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = "cross_validation_multi_trait_ml"
    )
    return(results_handling(
      GS_model = GS_model_cv,
      cv_results_raw = cv_results,
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_results_processed,
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_multi_trait_ml",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(multi_trait_asreml)) {
    # Multi-trait ASReml CV uses the same predict-or-extract fallback chain
    # the true-prediction path uses (predict.asreml -> extract via
    # mod$coefficients$random parsing on failure). Mirrors the multi_trait_ml
    # CV scaffolding so the downstream cv_results_processed plumbing is reused.
    if (!isTRUE(cv_evaluation_only)) {
      stop("Multi-trait ASReml CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) < 2L) {
      stop("Multi-trait ASReml-R requires at least two response columns.", call. = FALSE)
    }
    if (!identical(response_family, "gaussian")) {
      stop("Multi-trait ASReml-R currently supports gaussian traits only.", call. = FALSE)
    }
    asreml_kernel_bank <- gmatrix_kernel_model_ready_list %||% list()
    asreml_kernels <- gp_collect_kernel_inputs(
      gmatrix = asreml_kernel_bank[["gmatrix_model_ready"]] %||% gmatrix %||% NULL,
      omic1_kernel = asreml_kernel_bank[["omic1_kernel_model_ready"]] %||% NULL,
      omic2_kernel = asreml_kernel_bank[["omic2_kernel_model_ready"]] %||% NULL,
      omic3_kernel = asreml_kernel_bank[["omic3_kernel_model_ready"]] %||% NULL,
      kernel_list = asreml_kernel_bank
    )
    if (!length(asreml_kernels)) {
      stop("Multi-trait ASReml-R requires at least one relationship/kernel matrix.", call. = FALSE)
    }
    mt_asreml_pheno <- pheno_clean[["pheno_clean_data"]]
    mt_asreml_is_met <- gp_is_multi_environment_trait_panel(
      pheno_data = mt_asreml_pheno,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response = response,
      response_family = response_family
    )
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_k,
      feature_k_grid = feature_k_grid,
      feature_scoring_cv = feature_scoring_cv,
      context = "Multi-trait ASReml CV"
    )
    if (isTRUE(mt_asreml_is_met) && !is.null(feature_cfg$k)) {
      stop(
        "Grouped MT-MET ASReml CV does not yet support response-derived feature selection; ",
        "supply a fixed gmatrix or disable feature_scoring.",
        call. = FALSE
      )
    }
    if (isTRUE(mt_asreml_is_met) && !is.null(ctx$weights)) {
      stop(
        "Grouped MT-MET ASReml CV does not yet accept observation weights; weights were not silently ignored.",
        call. = FALSE
      )
    }
    cv_pipeline <- if (isTRUE(mt_asreml_is_met)) {
      gp_multitrait_asreml_mtmet_gaussian_cv(
        pheno_object = mt_asreml_pheno,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups,
        gmatrix = asreml_kernels[[1L]],
        kernel_list = if (length(asreml_kernels) > 1L) asreml_kernels[-1L] else NULL,
        cross_validation_meth = cross_validation_meth,
        nfolds = nfolds,
        replication = replication,
        eval_metrics = eval_metrics,
        heter_resid = heter_resid,
        var_cov_str = var_cov_str,
        inverse = inverse,
        epsilon = epsilon,
        engine = engine,
        workspace = workspace %||% 1e08,
        pworkspace = ctx$pworkspace %||% 1e06,
        maxit = maxit %||% 50L,
        ai_sing = TRUE,
        random_seed = random_seed %||% random_state %||% 123L
      )
    } else {
      gp_multitrait_asreml_gaussian_cv(
        pheno_object = mt_asreml_pheno,
        response = response,
        gen_name = gen_name,
        gmatrix = asreml_kernels[[1L]],
        kernel_list = if (length(asreml_kernels) > 1L) asreml_kernels[-1L] else NULL,
        cross_validation_meth = cross_validation_meth,
        nfolds = nfolds,
        replication = replication,
        eval_metrics = eval_metrics,
        trait_covariance = gp_multitrait_asreml_trait_covariance(var_cov_str),
        workspace = workspace %||% 1e08,
        pworkspace = ctx$pworkspace %||% 1e06,
        maxit = maxit %||% 50L,
        ai_sing = TRUE,
        random_seed = random_seed,
        feature_predictor_data = ml_dat_res[["merged_data"]][["merge_data"]] %||% NULL,
        feature_score_metadata = feature_score_metadata,
        feature_k = feature_cfg$k,
        feature_scoring_cv = feature_cfg$mode,
        feature_scoring_model = feature_scoring_model,
        feature_source_map = feature_source_map,
        gmatrix_method = gmatrix_method,
        ploidy = ploidy %||% "auto",
        feature_ridge_lambda = feature_ridge_lambda,
        feature_bayes_nIter = feature_bayes_nIter,
        feature_bayes_burnIn = feature_bayes_burnIn,
        feature_bayes_thin = feature_bayes_thin,
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        rf_n_jobs = rf_n_jobs
      )
    }
    cv_results <- cv_pipeline[["cv_results"]]
    cv_results_processed <- cv_pipeline[["cv_results_processed"]]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = if (isTRUE(mt_asreml_is_met)) {
        "cross_validation_multi_trait_multi_environment_asreml"
      } else {
        "cross_validation_multi_trait_asreml"
      }
    )
    return(results_handling(
      GS_model = GS_model_cv,
      cv_results_raw = cv_results,
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_results_processed,
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = if (isTRUE(mt_asreml_is_met)) {
        "CV_results_multi_trait_multi_environment_asreml"
      } else {
        "CV_results_multi_trait_asreml"
      },
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }
  if (isTRUE(multi_trait_dl)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Multi-trait DL CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    if (length(response) < 2L) {
      stop("Multi-trait DL requires at least two response columns.", call. = FALSE)
    }
    if (!identical(response_family, "gaussian")) {
      stop("Multi-trait DL currently supports gaussian traits only.", call. = FALSE)
    }
    if (length(GS_model_cv) != 1L || !GS_model_cv %in% gp_multitrait_dl_supported_models()) {
      stop(
        paste(
          "Multi-trait DL CV currently supports:",
          paste(gp_multitrait_dl_supported_models(), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
      stop("Multi-trait DL CV requires processed geno/omic features.", call. = FALSE)
    }
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_k,
      feature_k_grid = feature_k_grid,
      feature_scoring_cv = feature_scoring_cv,
      context = "Multi-trait DL CV"
    )
    cv_pipeline <- gp_multitrait_dl_gaussian_cv(
      model_type = GS_model_cv,
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = response,
      geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
      gen_name = gen_name,
      cross_validation_meth = cross_validation_meth,
      nfolds = nfolds,
      replication = replication,
      eval_metrics = eval_metrics,
      scaling = scaling,
      centering = centering,
      compile_model = compile_model,
      deterministic = deterministic,
      random_seed = random_seed,
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
      final_attention = final_attention,
      attention_across_multiple_layers = attention_across_multiple_layers,
      optimizer_name = optimizer_name,
      max_grad_norm = max_grad_norm,
      feature_score_metadata = feature_score_metadata,
      feature_k = feature_cfg$k,
      feature_scoring_cv = feature_cfg$mode,
      feature_scoring_model = feature_scoring_model,
      feature_source_map = feature_source_map,
      feature_ridge_lambda = feature_ridge_lambda,
      feature_bayes_nIter = feature_bayes_nIter,
      feature_bayes_burnIn = feature_bayes_burnIn,
      feature_bayes_thin = feature_bayes_thin,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      rf_n_jobs = rf_n_jobs
    )
    cv_results <- cv_pipeline[["cv_results"]]
    cv_results_processed <- cv_pipeline[["cv_results_processed"]]
    cv_run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = gp_detect_python(),
      preferred_python = gp_preferred_python(purpose = "dl"),
      context = "cross_validation_multi_trait_dl"
    )
    return(results_handling(
      GS_model = GS_model_cv,
      cv_results_raw = cv_results,
      res_model_output = NULL,
      res_summary_stat = NULL,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = cv_results_processed,
      run_metadata = cv_run_metadata,
      system_database = system_database,
      plot_filename = "CV_results_multi_trait_dl",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected,
      feature_score_metadata = feature_score_metadata
    ))
  }

  NULL
}
