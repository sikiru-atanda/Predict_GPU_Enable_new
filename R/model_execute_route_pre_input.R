gp_route_pre_input_specialized_models <- function(ctx) {
  parent_hybrid_mode <- isTRUE(ctx$hybrid_asreml) || isTRUE(ctx$hybrid_bayes) || isTRUE(ctx$hybrid_gp)
  feature_request <- isTRUE(ctx$feature_scoring) ||
    !is.null(ctx$feature_score_metadata) ||
    !is.null(ctx$feature_k) ||
    !is.null(ctx$feature_k_grid)
  if (isTRUE(parent_hybrid_mode) && isTRUE(feature_request)) {
    if (is.null(ctx$feature_score_metadata)) {
      stop(
        "Automatic phenotype-based feature scoring is not defined for parent-kernel hybrid models because phenotype rows are hybrids while marker rows are parents. Supply `feature_selected`, or supply externally calculated `feature_score_metadata` from an explicitly documented hybrid feature encoding.",
        call. = FALSE
      )
    }
    feature_cfg <- gp_feature_specialized_cv_config(
      feature_score_metadata = ctx$feature_score_metadata,
      feature_k = ctx$feature_k,
      feature_k_grid = ctx$feature_k_grid,
      feature_scoring_cv = ctx$feature_scoring_cv,
      context = "Parent-kernel hybrid prediction"
    )
    if (!identical(feature_cfg$mode, "fixed")) {
      stop(
        "Parent-kernel hybrid prediction currently accepts fixed external feature rankings only; fold-internal selection requires an explicit hybrid feature encoding.",
        call. = FALSE
      )
    }
    if (any(vapply(
      list(ctx$gmatrix, ctx$gkernel, ctx$female_gmatrix, ctx$male_gmatrix),
      Negate(is.null),
      logical(1L)
    ))) {
      stop(
        "A fixed parent-hybrid feature ranking cannot subset a supplied parent relationship matrix. Provide raw parent genotype matrices only, or externally rebuild and pass relationship matrices matching `feature_selected`.",
        call. = FALSE
      )
    }
    selected <- gp_feature_selected_by_source(
      ctx$feature_score_metadata,
      trait = as.character(ctx$response[[1L]]),
      k = feature_cfg$k
    )
    if (length(selected) == 1L && identical(names(selected), "merged")) {
      parent_columns <- unique(c(
        colnames(ctx$geno_data %||% NULL),
        colnames(ctx$female_geno_data %||% NULL),
        colnames(ctx$male_geno_data %||% NULL)
      ))
      if (all(selected[[1L]] %in% parent_columns)) {
        names(selected) <- "geno_data"
      }
    }
    ctx$feature_selected <- selected
  }
  ctx <- gp_feature_apply_raw_selection_to_inputs(ctx)
  gp_model_execute_import_context(ctx)
if (isTRUE(hybrid_asreml)) {
  hybrid_model_name <- as.character((GS_model %||% GS_model_cv)[1] %||% "")
  if (sum(
    isTRUE(multi_trait_asreml),
    isTRUE(multi_trait_gp),
    isTRUE(multi_trait_bayes),
    isTRUE(multi_trait_ml),
    isTRUE(multi_trait_dl),
    isTRUE(met_ml_dl),
    isTRUE(hybrid_bayes),
    isTRUE(hybrid_gp),
    isTRUE(hybrid_ml),
    isTRUE(hybrid_dl)
  ) > 0L) {
    stop("hybrid_asreml is a separate path; do not combine it with hybrid_bayes, hybrid_gp, hybrid_ml, hybrid_dl, multi-trait, or MET ML/DL paths.", call. = FALSE)
  }
  if (length(response) != 1L) {
    stop("Hybrid ASReml-R currently supports one gaussian response at a time.", call. = FALSE)
  }
  if (!identical(response_family, "gaussian")) {
    stop("Hybrid ASReml-R currently supports gaussian traits only.", call. = FALSE)
  }
  if (!identical(hybrid_model_name, "GBLUP")) {
    stop("Hybrid ASReml-R currently supports GS_model = 'GBLUP' only.", call. = FALSE)
  }
  if (length(engine) != 1L || !identical(as.character(engine), "asreml")) {
    stop("Hybrid ASReml-R requires engine = 'asreml'.", call. = FALSE)
  }
  if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
    stop("Provide both female_parent and male_parent column names for hybrid ASReml-R.", call. = FALSE)
  }
  if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
    stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
  }
    if (gp_formula_is_intercept_only(fixed)) fixed <- NULL  # ~1 is the default fixed part
    if (!gp_formula_is_intercept_only(fixed) || !is.null(random) || !is.null(cova) || !is.null(var_cov_str) || isTRUE(heter_resid)) {
      stop("Hybrid ASReml-R currently supports environment-aware prediction via heter_groups, but not custom fixed/random/variance structures yet.", call. = FALSE)
  }
  if (is.null(gmatrix) && is.null(female_gmatrix) && is.null(male_gmatrix) &&
      is.null(geno_data) && (is.null(female_geno_data) || is.null(male_geno_data))) {
    stop("Provide gmatrix, female_gmatrix/male_gmatrix, shared parent geno_data, or both female_geno_data and male_geno_data for hybrid ASReml-R.", call. = FALSE)
  }

  if (isTRUE(cross_validation)) {
    hybrid_cv_method <- cross_validation_meth %||% "Hybrid_Known_Parents"
    res_model_output <- gp_hybrid_asreml_gaussian_cv(
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = response,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      cross_validation_meth = hybrid_cv_method,
      eval_metrics = eval_metrics,
      gmatrix = gmatrix,
      female_gmatrix = female_gmatrix,
      male_gmatrix = male_gmatrix,
      geno_data = geno_data,
      female_geno_data = female_geno_data,
      male_geno_data = male_geno_data,
      gmatrix_method = gmatrix_method,
      ploidy = ploidy %||% "auto",
      include_sca = hybrid_include_sca,
      inverse = inverse,
      epsilon = epsilon,
      engine = engine,
      workspace = workspace,
      maxit = maxit,
      nfolds = nfolds,
      random_state = random_state,
      replication = replication
    )
    res_summary_stat <- list(
      hybrid_metric_summary = res_model_output[["cv_results_processed"]][["hybrid_metric_summary"]],
      hybrid_prediction_counts = res_model_output[["cv_results_processed"]][["hybrid_prediction_counts"]]
    )
    run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = "model_execute_hybrid_asreml_cv",
      extra_fields = c(
        gp_runtime_profile_metadata_fields(run_profile),
        list(
          task_unit = "hybrid_asreml_cv",
          cv_method = hybrid_cv_method
        )
      )
    )
    return(results_handling(
      GS_model = hybrid_model_name,
      res_model_output = NULL,
      res_summary_stat = res_summary_stat,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = res_model_output[["cv_results_processed"]],
      cv_results_raw = res_model_output[["cv_results_raw"]],
      run_metadata = run_metadata,
      system_database = system_database,
      plot_filename = "hybrid_asreml_cv",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected
    ))
  }

  res_model_output <- gp_hybrid_asreml_gaussian_model(
    pheno_object = pheno_clean[["pheno_clean_data"]],
    response = response,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      gmatrix = gmatrix,
    female_gmatrix = female_gmatrix,
    male_gmatrix = male_gmatrix,
    geno_data = geno_data,
    female_geno_data = female_geno_data,
    male_geno_data = male_geno_data,
    gmatrix_method = gmatrix_method,
    ploidy = ploidy %||% "auto",
    include_sca = hybrid_include_sca,
    inverse = inverse,
    epsilon = epsilon,
    engine = engine,
    workspace = workspace,
    maxit = maxit
  )
  res_summary_stat <- summary_statistics_hybrid(
    predicted_object = res_model_output[["predicted_values"]],
    response = response,
    female_parent = female_parent,
    male_parent = male_parent,
    mode = "hybrid_asreml",
    model_type = hybrid_model_name
  )
  run_metadata <- gp_runtime_metadata(
    execution_policy = NULL,
    python_path = NULL,
    preferred_python = NULL,
    context = "model_execute_hybrid_asreml",
    extra_fields = c(
      gp_runtime_profile_metadata_fields(run_profile),
      list(task_unit = "hybrid_asreml")
    )
  )
  return(results_handling(
    GS_model = hybrid_model_name,
    res_model_output = res_model_output,
    res_summary_stat = res_summary_stat,
    res_plot = NULL,
    res_plot_mean = NULL,
    res_plot_result_diagnostic = NULL,
    test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
    res_mod_results_cv_per_trait_model = NULL,
    res_plot_result_diagnostic_cv_only = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = NULL,
    run_metadata = run_metadata,
    system_database = system_database,
    plot_filename = "hybrid_asreml",
    plot_extension = plot_extension,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_units = plot_units,
    plot_dpi = plot_dpi,
    feature_selected = feature_selected
  ))
}

if (isTRUE(hybrid_bayes)) {
  hybrid_model_name <- as.character((GS_model %||% GS_model_cv)[1] %||% "")
  if (sum(
    isTRUE(multi_trait_asreml),
    isTRUE(multi_trait_gp),
    isTRUE(multi_trait_bayes),
    isTRUE(multi_trait_ml),
    isTRUE(multi_trait_dl),
    isTRUE(met_ml_dl),
    isTRUE(hybrid_asreml),
    isTRUE(hybrid_gp),
    isTRUE(hybrid_ml),
    isTRUE(hybrid_dl)
  ) > 0L) {
    stop("hybrid_bayes is a separate path; do not combine it with hybrid_asreml, hybrid_gp, hybrid_ml, hybrid_dl, multi-trait, or MET ML/DL paths.", call. = FALSE)
  }
  if (length(response) != 1L) {
    stop("Hybrid Bayesian models currently support one gaussian response at a time.", call. = FALSE)
  }
  if (!identical(response_family, "gaussian")) {
    stop("Hybrid Bayesian models currently support gaussian traits only.", call. = FALSE)
  }
  if (!hybrid_model_name %in% gp_hybrid_bayes_supported_models()) {
    stop(
      paste(
        "Hybrid Bayesian models currently support GS_model or GS_model_cv in:",
        paste(gp_hybrid_bayes_supported_models(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
    stop("Provide both female_parent and male_parent column names for hybrid Bayesian models.", call. = FALSE)
  }
  if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
    stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
  }
    if (gp_formula_is_intercept_only(fixed)) fixed <- NULL  # ~1 is the default fixed part
    if (!gp_formula_is_intercept_only(fixed) || !is.null(random) || !is.null(cova) || !is.null(var_cov_str) || isTRUE(heter_resid)) {
      stop("Hybrid Bayesian models currently support environment-aware prediction via heter_groups, but not custom fixed/random/variance structures yet.", call. = FALSE)
  }
  if (is.null(gmatrix) && is.null(female_gmatrix) && is.null(male_gmatrix) &&
      is.null(geno_data) && (is.null(female_geno_data) || is.null(male_geno_data))) {
    stop("Provide gmatrix, female_gmatrix/male_gmatrix, shared parent geno_data, or both female_geno_data and male_geno_data for hybrid Bayesian models.", call. = FALSE)
  }

  if (isTRUE(cross_validation)) {
    if (!isTRUE(cv_evaluation_only)) {
      stop("Hybrid Bayesian CV currently supports cv_evaluation_only = TRUE only.", call. = FALSE)
    }
    hybrid_cv_method <- cross_validation_meth %||% "Hybrid_Known_Parents"
    res_model_output <- gp_hybrid_bayes_gaussian_cv(
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = response,
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        cross_validation_meth = hybrid_cv_method,
      eval_metrics = eval_metrics,
      model_type = hybrid_model_name,
      gmatrix = gmatrix,
      female_gmatrix = female_gmatrix,
      male_gmatrix = male_gmatrix,
      geno_data = geno_data,
      female_geno_data = female_geno_data,
      male_geno_data = male_geno_data,
      gmatrix_method = gmatrix_method,
      include_sca = hybrid_include_sca,
      inverse = inverse,
      epsilon = epsilon,
      nIter = nIter,
      burnIn = burnIn,
      thin = thin,
      confidence_level = confidence_level,
      CI_width_thresholds = CI_width_thresholds,
      interval_width_high_threshold = interval_width_high_threshold,
      interval_width_low_threshold = interval_width_low_threshold,
      high_reliability_thres = high_reliability_thres,
      low_reliability_thres = low_reliability_thres,
      nfolds = nfolds,
      random_state = random_state,
      replication = replication
    )
    res_summary_stat <- list(
      hybrid_metric_summary = res_model_output[["cv_results_processed"]][["hybrid_metric_summary"]],
      hybrid_prediction_counts = res_model_output[["cv_results_processed"]][["hybrid_prediction_counts"]]
    )
    run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = "model_execute_hybrid_bayes_cv",
      extra_fields = c(
        gp_runtime_profile_metadata_fields(run_profile),
        list(
          task_unit = "hybrid_bayes_cv",
          cv_method = hybrid_cv_method
        )
      )
    )
    return(results_handling(
      GS_model = hybrid_model_name,
      res_model_output = NULL,
      res_summary_stat = res_summary_stat,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = res_model_output[["cv_results_processed"]],
      cv_results_raw = res_model_output[["cv_results_raw"]],
      run_metadata = run_metadata,
      system_database = system_database,
      plot_filename = "hybrid_bayes_cv",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected
    ))
  }

  res_model_output <- gp_hybrid_bayes_gaussian_model(
    pheno_object = pheno_clean[["pheno_clean_data"]],
    response = response,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      model_type = hybrid_model_name,
    gmatrix = gmatrix,
    female_gmatrix = female_gmatrix,
    male_gmatrix = male_gmatrix,
    geno_data = geno_data,
    female_geno_data = female_geno_data,
    male_geno_data = male_geno_data,
    gmatrix_method = gmatrix_method,
    include_sca = hybrid_include_sca,
    inverse = inverse,
    epsilon = epsilon,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    confidence_level = confidence_level,
    CI_width_thresholds = CI_width_thresholds,
    interval_width_high_threshold = interval_width_high_threshold,
    interval_width_low_threshold = interval_width_low_threshold,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    random_state = if (exists("random_state", inherits = TRUE) && !is.null(random_state)) random_state else 123L
  )
  res_summary_stat <- summary_statistics_hybrid(
    predicted_object = res_model_output[["predicted_values"]],
    response = response,
    female_parent = female_parent,
    male_parent = male_parent,
    mode = "hybrid_bayes",
    model_type = hybrid_model_name
  )
  run_metadata <- gp_runtime_metadata(
    execution_policy = NULL,
    python_path = NULL,
    preferred_python = NULL,
    context = "model_execute_hybrid_bayes",
    extra_fields = c(
      gp_runtime_profile_metadata_fields(run_profile),
      list(task_unit = "hybrid_bayes")
    )
  )
  return(results_handling(
    GS_model = hybrid_model_name,
    res_model_output = res_model_output,
    res_summary_stat = res_summary_stat,
    res_plot = NULL,
    res_plot_mean = NULL,
    res_plot_result_diagnostic = NULL,
    test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
    res_mod_results_cv_per_trait_model = NULL,
    res_plot_result_diagnostic_cv_only = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = NULL,
    run_metadata = run_metadata,
    system_database = system_database,
    plot_filename = "hybrid_bayes",
    plot_extension = plot_extension,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_units = plot_units,
    plot_dpi = plot_dpi,
    feature_selected = feature_selected
  ))
}

if (isTRUE(hybrid_gp)) {
  hybrid_model_name <- as.character((GS_model %||% GS_model_cv)[1] %||% "GP")
  if (sum(
    isTRUE(multi_trait_asreml),
    isTRUE(multi_trait_gp),
    isTRUE(multi_trait_bayes),
    isTRUE(multi_trait_ml),
    isTRUE(multi_trait_dl),
    isTRUE(met_ml_dl),
    isTRUE(hybrid_asreml),
    isTRUE(hybrid_bayes),
    isTRUE(hybrid_ml),
    isTRUE(hybrid_dl)
  ) > 0L) {
    stop("hybrid_gp is a separate path; do not combine it with hybrid_asreml, hybrid_bayes, hybrid_ml, hybrid_dl, dedicated multi-trait, or MET ML/DL flags. For hybrid GP multi-trait use response = c(...); for hybrid GP MET use heter_groups.", call. = FALSE)
  }
  hybrid_responses <- unique(as.character(response))
  hybrid_responses <- hybrid_responses[nzchar(hybrid_responses)]
  if (!identical(response_family, "gaussian")) {
    stop("Hybrid GP currently supports gaussian traits only.", call. = FALSE)
  }
  if (!hybrid_model_name %in% gp_hybrid_gp_supported_models()) {
    stop(
      paste(
        "Hybrid GP currently supports GS_model or GS_model_cv in:",
        paste(gp_hybrid_gp_supported_models(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
    stop("Provide both female_parent and male_parent column names for hybrid GP.", call. = FALSE)
  }
  if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
    stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
  }
  if (gp_formula_is_intercept_only(fixed)) fixed <- NULL  # ~1 is the default fixed part
  if (!gp_formula_is_intercept_only(fixed) || !is.null(random) || !is.null(cova) || !is.null(var_cov_str) || isTRUE(heter_resid)) {
    stop("Hybrid GP currently supports environment-aware prediction via heter_groups, but not custom fixed/random/variance structures yet.", call. = FALSE)
  }
  if (is.null(gmatrix) && is.null(female_gmatrix) && is.null(male_gmatrix) &&
      is.null(geno_data) && (is.null(female_geno_data) || is.null(male_geno_data))) {
    stop("Provide gmatrix, female_gmatrix/male_gmatrix, shared parent geno_data, or both female_geno_data and male_geno_data for hybrid GP.", call. = FALSE)
  }

  if (isTRUE(cross_validation)) {
    hybrid_cv_method <- cross_validation_meth %||% "Hybrid_Known_Parents"
    res_model_output <- if (length(hybrid_responses) > 1L) {
      gp_hybrid_gp_multi_response_cv(
        pheno_object = pheno_clean[["pheno_clean_data"]],
        response = hybrid_responses,
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        cross_validation_meth = hybrid_cv_method,
        eval_metrics = eval_metrics,
        model_type = hybrid_model_name,
        gmatrix = gmatrix,
        female_gmatrix = female_gmatrix,
        male_gmatrix = male_gmatrix,
        geno_data = geno_data,
        female_geno_data = female_geno_data,
        male_geno_data = male_geno_data,
        gmatrix_method = gmatrix_method,
        include_sca = hybrid_include_sca,
        lambda = hybrid_gp_lambda,
        lambda_grid = hybrid_gp_lambda_grid,
        component_weights = hybrid_gp_component_weights,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        gp_backend = gp_backend,
        gp_auto_policy = TRUE,
        nfolds = nfolds,
        random_state = random_state,
        replication = replication
      )
    } else {
      gp_hybrid_gp_gaussian_cv(
        pheno_object = pheno_clean[["pheno_clean_data"]],
        response = hybrid_responses[[1L]],
        gen_name = gen_name,
        female_parent = female_parent,
        male_parent = male_parent,
        heter_groups = heter_groups,
        cross_validation_meth = hybrid_cv_method,
        eval_metrics = eval_metrics,
        model_type = hybrid_model_name,
        gmatrix = gmatrix,
        female_gmatrix = female_gmatrix,
        male_gmatrix = male_gmatrix,
        geno_data = geno_data,
        female_geno_data = female_geno_data,
        male_geno_data = male_geno_data,
        gmatrix_method = gmatrix_method,
        include_sca = hybrid_include_sca,
        lambda = hybrid_gp_lambda,
        lambda_grid = hybrid_gp_lambda_grid,
        component_weights = hybrid_gp_component_weights,
        env_similarity = env_similarity,
        env_ids = env_ids,
        env_covariates = env_covariates,
        reaction_norm_feature_qc = reaction_norm_feature_qc,
        kenv_kernel = kenv_kernel,
        kenv_bandwidth = kenv_bandwidth,
        kenv_kernel_kwargs = kenv_kernel_kwargs,
        gp_backend = gp_backend,
        gp_auto_policy = TRUE,
        nfolds = nfolds,
        random_state = random_state,
        replication = replication
      )
    }
    res_summary_stat <- list(
      hybrid_metric_summary = res_model_output[["cv_results_processed"]][["hybrid_metric_summary"]],
      hybrid_prediction_counts = res_model_output[["cv_results_processed"]][["hybrid_prediction_counts"]]
    )
    run_metadata <- gp_runtime_metadata(
      execution_policy = NULL,
      python_path = NULL,
      preferred_python = NULL,
      context = "model_execute_hybrid_gp_cv",
      extra_fields = c(
        gp_runtime_profile_metadata_fields(run_profile),
        list(
          task_unit = "hybrid_gp_cv",
          cv_method = hybrid_cv_method,
          gp_model = gp_hybrid_gp_model_label(hybrid_model_name),
          trait_mode = if (length(hybrid_responses) > 1L) "joint_cross_trait_hybrid_gp" else "single_trait_hybrid_gp",
          trait_count = length(hybrid_responses)
        )
      )
    )
    return(results_handling(
      GS_model = gp_hybrid_gp_model_label(hybrid_model_name),
      res_model_output = NULL,
      res_summary_stat = res_summary_stat,
      res_plot = NULL,
      res_plot_mean = NULL,
      res_plot_result_diagnostic = NULL,
      test_diagonistic_plots = NULL,
      res_mod_results_cv_per_trait_model = NULL,
      res_plot_result_diagnostic_cv_only = NULL,
      geno_qc_stat = NULL,
      cv_results_processed = res_model_output[["cv_results_processed"]],
      cv_results_raw = res_model_output[["cv_results_raw"]],
      run_metadata = run_metadata,
      system_database = system_database,
      plot_filename = "hybrid_gp_cv",
      plot_extension = plot_extension,
      plot_width = plot_width,
      plot_height = plot_height,
      plot_units = plot_units,
      plot_dpi = plot_dpi,
      feature_selected = feature_selected
    ))
  }

  res_model_output <- if (length(hybrid_responses) > 1L) {
    gp_hybrid_gp_multi_response_model(
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = hybrid_responses,
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      model_type = hybrid_model_name,
      gmatrix = gmatrix,
      female_gmatrix = female_gmatrix,
      male_gmatrix = male_gmatrix,
      geno_data = geno_data,
      female_geno_data = female_geno_data,
      male_geno_data = male_geno_data,
      gmatrix_method = gmatrix_method,
      include_sca = hybrid_include_sca,
      lambda = hybrid_gp_lambda,
      lambda_grid = hybrid_gp_lambda_grid,
      component_weights = hybrid_gp_component_weights,
      env_similarity = env_similarity,
      env_ids = env_ids,
      env_covariates = env_covariates,
      reaction_norm_feature_qc = reaction_norm_feature_qc,
      kenv_kernel = kenv_kernel,
      kenv_bandwidth = kenv_bandwidth,
      kenv_kernel_kwargs = kenv_kernel_kwargs,
      gp_backend = gp_backend,
      gp_auto_policy = TRUE,
      random_state = random_state
    )
  } else {
    gp_hybrid_gp_gaussian_model(
      pheno_object = pheno_clean[["pheno_clean_data"]],
      response = hybrid_responses[[1L]],
      gen_name = gen_name,
      female_parent = female_parent,
      male_parent = male_parent,
      heter_groups = heter_groups,
      model_type = hybrid_model_name,
      gmatrix = gmatrix,
      female_gmatrix = female_gmatrix,
      male_gmatrix = male_gmatrix,
      geno_data = geno_data,
      female_geno_data = female_geno_data,
      male_geno_data = male_geno_data,
      gmatrix_method = gmatrix_method,
      include_sca = hybrid_include_sca,
      lambda = hybrid_gp_lambda,
      lambda_grid = hybrid_gp_lambda_grid,
      component_weights = hybrid_gp_component_weights,
      env_similarity = env_similarity,
      env_ids = env_ids,
      env_covariates = env_covariates,
      reaction_norm_feature_qc = reaction_norm_feature_qc,
      kenv_kernel = kenv_kernel,
      kenv_bandwidth = kenv_bandwidth,
      kenv_kernel_kwargs = kenv_kernel_kwargs,
      gp_backend = gp_backend,
      gp_auto_policy = TRUE,
      random_state = random_state
    )
  }
  res_summary_stat <- if (length(hybrid_responses) > 1L) {
    gp_hybrid_gp_multi_response_summary(
      predicted_object = res_model_output[["predicted_values"]],
      response = hybrid_responses,
      female_parent = female_parent,
      male_parent = male_parent,
      model_type = gp_hybrid_gp_model_label(hybrid_model_name)
    )
  } else {
    summary_statistics_hybrid(
      predicted_object = res_model_output[["predicted_values"]],
      response = hybrid_responses[[1L]],
      female_parent = female_parent,
      male_parent = male_parent,
      mode = "hybrid_gp",
      model_type = gp_hybrid_gp_model_label(hybrid_model_name)
    )
  }
  run_metadata <- gp_runtime_metadata(
    execution_policy = NULL,
    python_path = NULL,
    preferred_python = NULL,
    context = "model_execute_hybrid_gp",
    extra_fields = c(
      gp_runtime_profile_metadata_fields(run_profile),
      list(
        task_unit = "hybrid_gp",
        gp_model = gp_hybrid_gp_model_label(hybrid_model_name),
        trait_mode = if (length(hybrid_responses) > 1L) "joint_cross_trait_hybrid_gp" else "single_trait_hybrid_gp",
        trait_count = length(hybrid_responses)
      )
    )
  )
  return(results_handling(
    GS_model = gp_hybrid_gp_model_label(hybrid_model_name),
    res_model_output = res_model_output,
    res_summary_stat = res_summary_stat,
    res_plot = NULL,
    res_plot_mean = NULL,
    res_plot_result_diagnostic = NULL,
    test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
    res_mod_results_cv_per_trait_model = NULL,
    res_plot_result_diagnostic_cv_only = NULL,
    geno_qc_stat = NULL,
    cv_results_processed = NULL,
    run_metadata = run_metadata,
    system_database = system_database,
    plot_filename = "hybrid_gp",
    plot_extension = plot_extension,
    plot_width = plot_width,
    plot_height = plot_height,
    plot_units = plot_units,
    plot_dpi = plot_dpi,
    feature_selected = feature_selected
  ))
}


  NULL
}
