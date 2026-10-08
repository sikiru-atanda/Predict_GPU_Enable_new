gp_route_true_prediction_specialized_models <- function(ctx) {
  gp_model_execute_import_context(ctx)
  selected_joint_kernel_bank <- function(context) {
    if (is.null(feature_score_metadata) || is.null(feature_k)) return(NULL)
    if (any(vapply(
      list(ctx$gmatrix, ctx$gkernel, ctx$omic1_kernel, ctx$omic2_kernel, ctx$omic3_kernel, ctx$kernel_list),
      Negate(is.null),
      logical(1L)
    ))) {
      stop(
        context,
        " cannot feature-select user-supplied precomputed kernels. Supply the corresponding raw marker/omics inputs without those kernels, or provide externally rebuilt matching kernels.",
        call. = FALSE
      )
    }
    joint <- gp_feature_joint_kernel_bank(
      predictor_data = ml_dat_res[["merged_data"]][["merge_data"]],
      pheno_data = pheno_clean[["pheno_clean_data"]],
      response = response,
      gen_name = gen_name,
      feature_score_metadata = feature_score_metadata,
      k = feature_k,
      source_map = feature_source_map,
      source_matrices = list(
        geno_data = geno_omic_model_ready_list[["geno_model_ready"]] %||% NULL,
        omic1_data = geno_omic_model_ready_list[["omic1_model_ready"]] %||% NULL,
        omic2_data = geno_omic_model_ready_list[["omic2_model_ready"]] %||% NULL,
        omic3_data = geno_omic_model_ready_list[["omic3_model_ready"]] %||% NULL
      ),
      gmatrix_method = gmatrix_method,
      kernel_method = kernel_method,
      ploidy = ploidy %||% "auto",
      scaling = scaling,
      centering = centering,
      context = context
    )
    joint$kernel_bank
  }
 if (isTRUE(multi_trait_bayes)) {
   if (!identical(response_family, "gaussian")) {
     stop("Multi-trait Bayesian BGLR currently supports gaussian traits only.", call. = FALSE)
   }
   if (isTRUE(cross_validation)) {
     stop("Multi-trait Bayesian BGLR currently supports true prediction only.", call. = FALSE)
   }
   if (length(GS_model) != 1L || !GS_model %in% gp_multitrait_bayes_supported_models()) {
     stop(
       paste(
         "Multi-trait Bayesian BGLR supports GS_model in:",
         paste(gp_multitrait_bayes_supported_models(), collapse = ", ")
       ),
       call. = FALSE
     )
   }
   mt_bayes_pheno <- pheno_clean[["pheno_clean_data"]]
   mt_bayes_kernel_bank <- selected_joint_kernel_bank("Multi-trait Bayesian true prediction") %||%
     gmatrix_kernel_model_ready_list
   kernels <- gp_collect_kernel_inputs(
     gmatrix = mt_bayes_kernel_bank[["gmatrix_model_ready"]] %||% gmatrix,
     omic1_kernel = mt_bayes_kernel_bank[["omic1_kernel_model_ready"]] %||% NULL,
     omic2_kernel = mt_bayes_kernel_bank[["omic2_kernel_model_ready"]] %||% NULL,
     omic3_kernel = mt_bayes_kernel_bank[["omic3_kernel_model_ready"]] %||% NULL,
     kernel_list = mt_bayes_kernel_bank
   )
   if (!length(kernels)) {
     stop("Multi-trait Bayesian BGLR requires at least one model-ready relationship/kernel matrix.", call. = FALSE)
   }
   mt_bayes_is_met <- gp_is_multi_environment_trait_panel(
     pheno_data = mt_bayes_pheno,
     gen_name = gen_name,
     heter_groups = heter_groups,
     response = response,
     response_family = response_family
   )
   bayes_para <- bayes_parameter_check(
     nIter = nIter,
     burnIn = burnIn,
     thin = thin
   )
   mt_bayes_fit <- bayes_multitrait_joint_fit(
     pheno_data = mt_bayes_pheno,
     response = response,
     gen_name = gen_name,
     heter_groups = if (isTRUE(mt_bayes_is_met)) heter_groups else NULL,
     kernels = kernels,
     GS_model = GS_model,
     bayes_para = bayes_para,
     heter_resid = heter_resid,
     confidence_level = confidence_level %||% 0.95,
     CI_width_thresholds = CI_width_thresholds %||% c(0.33, 0.66),
     high_reliability_thres = high_reliability_thres %||% 0.9,
     low_reliability_thres = low_reliability_thres %||% 0.5,
     verbose = FALSE,
     random_state = if (exists("random_state", inherits = TRUE) && !is.null(random_state)) random_state else 123L
   )
   res_model_output <- mt_bayes_fit[["bayes_result"]]
  res_model_output <- gp_reconcile_multitrait_model_output(
    res_model_output,
    gen_name = gen_name,
    heter_groups = if (isTRUE(mt_bayes_is_met)) heter_groups else NULL,
    model = GS_model,
    response_family = response_family
   )
   pred_long <- res_model_output[["predicted_values"]]
   model_public <- gp_public_model_label(GS_model)[1]
   mode <- if (isTRUE(mt_bayes_is_met)) "multi_trait_multi_environment_bayes" else "multi_trait_bayes"
   res_summary_stat <- gp_model_execute_multitrait_summary(
     pred_long = pred_long,
     response = response,
     model_type = model_public,
     mode = mode
   )
   diag_plot <- gp_multitrait_ml_diagnostic_plot(
     pred_long,
     model_label = model_public
   )
   res_model_output[["diagnostic_plots"]] <- diag_plot
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = NULL,
     preferred_python = NULL,
     context = paste0("model_execute_", mode),
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(
         task_unit = mode,
         task_traits = length(response),
         bglr_function = "Multitrait",
         eta_cov_argument = "Cov"
       )
     )
   )
   handled_result <- results_handling(
     GS_model = model_public,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = diag_plot,
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = mode,
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = model_public,
      trait_name = if (length(response) > 1L) "multi_trait" else response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(gp_public_multitrait_run)) {
   if (!identical(response_family, "gaussian")) {
     stop("Multi-trait GP currently supports gaussian traits only.", call. = FALSE)
   }
   mt_gp_kernel_bank <- selected_joint_kernel_bank("Multi-trait GP true prediction") %||%
     gmatrix_kernel_model_ready_list
   mt_gp_gmatrix <- mt_gp_kernel_bank[["gmatrix_model_ready"]] %||% gmatrix
   if (is.null(mt_gp_gmatrix)) {
     stop("Multi-trait GP requires gmatrix input.", call. = FALSE)
   }
   mt_gp_pheno <- pheno_clean[["pheno_clean_data"]]
   mt_gp_is_met <- gp_is_multi_environment_trait_panel(
     pheno_data = mt_gp_pheno,
     gen_name = gen_name,
     heter_groups = heter_groups,
     response = response,
     response_family = response_family
   )
    mt_gp_mode <- if (isTRUE(mt_gp_is_met)) {
      "multi_trait_multi_environment_gp"
    } else {
      "multi_trait_gp"
    }
    mt_gp_route <- gp_multitrait_gp_route_spec(
      model_name = GS_model,
      is_met = mt_gp_is_met,
      fa_rank = gp_fa_rank
    )
    mt_gp_fit <- if (isTRUE(mt_gp_is_met)) {
      gp_multi_trait_met_model(
       force_prediction_se = isTRUE(gp_force_prediction_se),
       pheno_data = mt_gp_pheno,
       gmatrix = mt_gp_gmatrix,
       response = response,
       gen_name = gen_name,
       heter_groups = heter_groups,
       test_set = test_set,
       kernel_list = mt_gp_kernel_bank,
       kernel_weights = lowrank_kernel_weights,
       env_similarity = env_similarity,
       env_ids = env_ids,
       env_covariates = env_covariates,
       reaction_norm_feature_qc = reaction_norm_feature_qc,
       kenv_kernel = kenv_kernel,
       kenv_bandwidth = kenv_bandwidth,
       kenv_kernel_kwargs = kenv_kernel_kwargs,
       varcomp_mode = mt_gp_route$varcomp_mode,
       gp_factor_cache = gp_factor_cache,
       return_se = isTRUE(gp_return_se) || isTRUE(gp_full_vc),
       prediction_output = "all",
       return_trait_correlations = isTRUE(gp_return_trait_correlations),
        trait_structure = mt_gp_route$trait_structure,
        trait_fa_rank = mt_gp_route$trait_fa_rank,
        gxe_trait_structure = mt_gp_route$gxe_trait_structure,
        gxe_trait_fa_rank = mt_gp_route$gxe_trait_fa_rank,
        seed = random_seed %||% random_state %||% 12345L
      )
   } else {
     gp_multi_trait_model(
       estimate_kernel_weights = isTRUE(gp_estimate_kernel_weights),
       force_prediction_se = isTRUE(gp_force_prediction_se),
       pheno_data = mt_gp_pheno,
       gmatrix = mt_gp_gmatrix,
       response = response,
       gen_name = gen_name,
       test_set = test_set,
       kernel_list = mt_gp_kernel_bank,
       kernel_weights = lowrank_kernel_weights,
        varcomp_mode = mt_gp_route$varcomp_mode,
        gp_factor_cache = gp_factor_cache,
        return_se = isTRUE(gp_return_se) || isTRUE(gp_full_vc),
        prediction_output = "all",
        return_trait_correlations = isTRUE(gp_return_trait_correlations),
        trait_structure = mt_gp_route$trait_structure,
        trait_fa_rank = mt_gp_route$trait_fa_rank,
        seed = random_seed %||% random_state %||% 12345L
      )
   }
   pred_long <- gp_model_execute_multitrait_predicted_values(
     predictions = mt_gp_fit[["predictions"]],
     pheno_data = mt_gp_pheno,
     response = response,
     gen_name = gen_name,
     heter_groups = if (isTRUE(mt_gp_is_met)) heter_groups else NULL
   )
   total_pred_long <- if (isTRUE(mt_gp_is_met) &&
                          exists("gp_multitrait_across_environment_prediction", mode = "function")) {
     gp_multitrait_across_environment_prediction(pred_long, gen_name = gen_name,
                                                  heter_groups = heter_groups)
   } else {
     NULL
   }
   GS_model_public <- gp_public_model_label(GS_model)[1]
   GS_model_canonical <- gp_canonicalize_supported_model_names(GS_model)[1]
   model_parameters <- data.frame(
     stat = c(
       "mode",
       "gp_model",
        "gp_model_canonical",
        "gp_backend_method",
        "gp_varcomp_mode_effective",
        "gp_trait_structure",
        "gp_trait_fa_rank",
        "gp_gxe_trait_structure",
        "gp_gxe_trait_fa_rank",
        "gp_execution",
       "bridge_command",
       "gp_return_se",
       "gp_return_trait_correlations"
     ),
     summary = c(
       mt_gp_mode,
        GS_model_public,
        GS_model_canonical,
        mt_gp_route$backend_method,
        mt_gp_route$varcomp_mode,
        mt_gp_route$trait_structure,
        mt_gp_route$trait_fa_rank %||% NA_integer_,
        mt_gp_route$gxe_trait_structure %||% NA_character_,
        mt_gp_route$gxe_trait_fa_rank %||% NA_integer_,
        mt_gp_fit[["info"]][["gp_execution"]] %||% mt_gp_fit[["info"]][["execution"]] %||% NA_character_,
       mt_gp_fit[["info"]][["bridge_command"]] %||% NA_character_,
       isTRUE(gp_return_se) || isTRUE(gp_full_vc),
       isTRUE(gp_return_trait_correlations)
     ),
     stringsAsFactors = FALSE
    )
    mt_kernel_names <- as.character(mt_gp_fit[["info"]][["kernel_names"]] %||% character())
    mt_kernel_weights <- mt_gp_fit[["info"]][["kernel_weights"]] %||% numeric()
    model_parameters <- rbind(
      model_parameters,
      gp_multi_kernel_parameter_rows(
        kernel_names = mt_kernel_names,
        kernel_count = length(mt_kernel_names),
        strategy = "GP_weighted_additive_kernel_bank_shared_trait_covariance"
      ),
      data.frame(
        stat = c(
          "gp_kernel_weights_effective",
          "gp_kernel_covariance_basis",
          "gp_independent_kernel_covariances_estimated"
        ),
        summary = c(
          if (length(mt_kernel_weights)) {
            paste0(
              names(mt_kernel_weights), "=",
              format(as.numeric(mt_kernel_weights), digits = 15L, trim = TRUE),
              collapse = ";"
            )
          } else {
            NA_character_
          },
          mt_gp_fit[["info"]][["kernel_covariance_basis"]] %||% NA_character_,
          as.character(isTRUE(mt_gp_fit[["info"]][["independent_kernel_covariances_estimated"]]))
        ),
        stringsAsFactors = FALSE
      )
    )
   res_model_output <- list(
     predicted_values = pred_long,
     multitrait_prediction_wide = gp_model_execute_multitrait_prediction_wide(
       pred_long,
       gen_name = gen_name,
       heter_groups = if (isTRUE(mt_gp_is_met)) heter_groups else NULL
     ),
     multitrait_trait_counts = gp_multitrait_ml_trait_counts(pred_long, response),
     model_parameters = model_parameters,
     gp_result = mt_gp_fit[["result"]],
     gp_info = mt_gp_fit[["info"]],
     # Private fitted covariance payload used by the public variance-component
     # formatter.  results_handling() removes it before native-file export.
     gp_fit = mt_gp_fit[["fit"]],
     diagnostic_plots = gp_multitrait_ml_diagnostic_plot(
       pred_long,
       model_label = GS_model_public
     )
   )
   if (is.data.frame(total_pred_long) && nrow(total_pred_long)) {
     res_model_output[["across_environment_predicted_values"]] <- total_pred_long
     res_model_output[["Total_Predicted_value"]] <- total_pred_long
   }
   if (!is.null(mt_gp_fit[["result"]][["genetic_correlation"]])) {
     res_model_output[["genetic_correlation"]] <- mt_gp_fit[["result"]][["genetic_correlation"]]
   }
   if (!is.null(mt_gp_fit[["result"]][["genetic_covariance"]])) {
     res_model_output[["genetic_covariance"]] <- mt_gp_fit[["result"]][["genetic_covariance"]]
   }
   if (!is.null(mt_gp_fit[["result"]][["gxe_correlation"]])) {
     res_model_output[["gxe_correlation"]] <- mt_gp_fit[["result"]][["gxe_correlation"]]
   }
   if (!is.null(mt_gp_fit[["result"]][["gxe_covariance"]])) {
     res_model_output[["gxe_covariance"]] <- mt_gp_fit[["result"]][["gxe_covariance"]]
   }
   if (!is.null(mt_gp_fit[["result"]][["residual_correlation"]])) {
     res_model_output[["residual_correlation"]] <- mt_gp_fit[["result"]][["residual_correlation"]]
   }
   if (!is.null(mt_gp_fit[["result"]][["residual_covariance"]])) {
     res_model_output[["residual_covariance"]] <- mt_gp_fit[["result"]][["residual_covariance"]]
   }
   res_model_output <- gp_standardize_public_model_result(
    res_model_output,
    gen_name = gen_name,
    heter_groups = if (isTRUE(mt_gp_is_met)) heter_groups else NULL,
    model = GS_model_canonical,
    response_family = response_family
   )
   res_summary_stat <- gp_model_execute_multitrait_summary(
     pred_long = pred_long,
     response = response,
     model_type = GS_model_public,
     mode = mt_gp_mode
   )
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = gp_detect_gp_python(),
     preferred_python = gp_detect_gp_python(),
     purpose = "gp",
     include_accelerator = FALSE,
     context = paste0("model_execute_", mt_gp_mode),
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(
          task_unit = mt_gp_mode,
          task_traits = length(response),
          gp_model = GS_model_public,
          gp_model_canonical = GS_model_canonical,
          gp_backend_method = mt_gp_route$backend_method,
          gp_varcomp_mode_effective = mt_gp_route$varcomp_mode,
          gp_trait_structure = mt_gp_route$trait_structure,
          gp_trait_fa_rank = mt_gp_route$trait_fa_rank %||% NA_integer_,
          gp_gxe_trait_structure = mt_gp_route$gxe_trait_structure %||% NA_character_,
          gp_gxe_trait_fa_rank = mt_gp_route$gxe_trait_fa_rank %||% NA_integer_,
          gp_execution = mt_gp_fit[["info"]][["gp_execution"]] %||% mt_gp_fit[["info"]][["execution"]] %||% NA_character_,
         bridge_command = mt_gp_fit[["info"]][["bridge_command"]] %||% NA_character_
       )
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model_public,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = mt_gp_mode,
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model_public,
      trait_name = if (length(response) > 1L) "multi_trait" else response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(hybrid_ml)) {
   hybrid_model_name <- as.character((GS_model %||% GS_model_cv)[1] %||% "")
   if (isTRUE(met_ml_dl)) {
     stop("Hybrid ML and MET ML/DL are separate paths; combine them in a later implementation.", call. = FALSE)
   }
   if (length(response) != 1L) {
     stop("Hybrid ML currently supports one gaussian response at a time.", call. = FALSE)
   }
   if (!identical(response_family, "gaussian")) {
     stop("Hybrid ML currently supports gaussian traits only.", call. = FALSE)
   }
   if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
     stop("Provide both female_parent and male_parent column names for hybrid ML.", call. = FALSE)
   }
   if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
     stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
   }
   if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
     stop("Hybrid ML requires processed geno/omic features keyed by hybrid ID.", call. = FALSE)
   }
   ml_dat_res <- gp_feature_apply_to_ml_dat_res(
     ml_dat_res = ml_dat_res,
     feature_score_metadata = feature_score_metadata,
     trait = response[[1]],
     k = feature_k
   )
   if (isTRUE(cross_validation)) {
     if (length(GS_model_cv) != 1L || !GS_model_cv %in% gp_hybrid_ml_supported_models()) {
       stop(
         paste(
           "Hybrid ML currently supports:",
           paste(gp_hybrid_ml_supported_models(), collapse = ", ")
         ),
         call. = FALSE
       )
     }
     hybrid_cv_method <- cross_validation_meth %||% "Hybrid_Known_Parents"
     res_model_output <- gp_hybrid_ml_gaussian_cv(
       model_type = GS_model_cv,
       pheno_object = pheno_data %||% pheno_clean[["pheno_clean_data"]],
       response = response[[1]],
       gen_name = gen_name,
       female_parent = female_parent,
       male_parent = male_parent,
       heter_groups = heter_groups,
       geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
       cross_validation_meth = hybrid_cv_method,
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
       context = "model_execute_hybrid_ml_cv",
       extra_fields = c(
         gp_runtime_profile_metadata_fields(run_profile),
         list(
           task_unit = "hybrid_ml_cv",
           cv_method = hybrid_cv_method
         )
       )
     )
     return(results_handling(
       GS_model = GS_model_cv,
       res_model_output = NULL,
       res_summary_stat = res_summary_stat,
       res_plot = NULL,
       res_plot_mean = NULL,
       res_plot_result_diagnostic = NULL,
       test_diagonistic_plots = NULL,
       res_mod_results_cv_per_trait_model = NULL,
       res_plot_result_diagnostic_cv_only = NULL,
       geno_qc_stat = geno_qc_stat,
       cv_results_processed = res_model_output[["cv_results_processed"]],
       cv_results_raw = res_model_output[["cv_results"]],
       run_metadata = run_metadata,
       system_database = system_database,
       plot_filename = "hybrid_ml_cv",
       plot_extension = plot_extension,
       plot_width = plot_width,
       plot_height = plot_height,
       plot_units = plot_units,
       plot_dpi = plot_dpi,
       feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
     ))
   }
   if (length(GS_model) != 1L || !GS_model %in% gp_hybrid_ml_supported_models()) {
     stop(
       paste(
         "Hybrid ML currently supports:",
         paste(gp_hybrid_ml_supported_models(), collapse = ", ")
       ),
       call. = FALSE
     )
   }
   res_model_output <- gp_hybrid_ml_gaussian_model(
     model_type = GS_model,
     pheno_object = pheno_data %||% pheno_clean[["pheno_clean_data"]],
     response = response[[1]],
       gen_name = gen_name,
       female_parent = female_parent,
       male_parent = male_parent,
       heter_groups = heter_groups,
       geno_omic_object = gp_hybrid_bind_feature_blocks(
         train_block = ml_dat_res[["merged_data"]][["merge_data"]],
         test_block = ml_dat_res[["merged_data_test"]]
       ),
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
      n_bootstrap = n_bootstrap,
      random_seed = random_state,
      confidence_level = confidence_level %||% 0.95,
      internal_cv_nfolds = internal_cv_nfolds,
      scaling = scaling,
      centering = centering
   )
   res_summary_stat <- summary_statistics_hybrid(
     predicted_object = res_model_output[["predicted_values"]],
     response = response[[1]],
     female_parent = female_parent,
     male_parent = male_parent,
     mode = "hybrid_ml",
     model_type = GS_model
   )
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = NULL,
     preferred_python = NULL,
     context = "model_execute_hybrid_ml",
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(task_unit = "hybrid_ml")
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = "hybrid_ml",
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model,
      trait_name = response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(hybrid_dl)) {
   hybrid_model_name <- as.character((GS_model %||% GS_model_cv)[1] %||% "")
   if (isTRUE(met_ml_dl)) {
     stop("Hybrid DL and MET ML/DL are separate paths; combine them in a later implementation.", call. = FALSE)
   }
   if (length(response) != 1L) {
     stop("Hybrid DL currently supports one gaussian response at a time.", call. = FALSE)
   }
   if (!identical(response_family, "gaussian")) {
     stop("Hybrid DL currently supports gaussian traits only.", call. = FALSE)
   }
   if (is.null(female_parent) || !nzchar(female_parent) || is.null(male_parent) || !nzchar(male_parent)) {
     stop("Provide both female_parent and male_parent column names for hybrid DL.", call. = FALSE)
   }
   if (!all(c(female_parent, male_parent) %in% names(pheno_clean[["pheno_clean_data"]]))) {
     stop("female_parent and male_parent columns were not found in the phenotype data.", call. = FALSE)
   }
   if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
     stop("Hybrid DL requires processed geno/omic features keyed by hybrid ID.", call. = FALSE)
   }
   ml_dat_res <- gp_feature_apply_to_ml_dat_res(
     ml_dat_res = ml_dat_res,
     feature_score_metadata = feature_score_metadata,
     trait = response[[1]],
     k = feature_k
   )
   if (length(GS_model) != 1L || !GS_model %in% gp_hybrid_dl_supported_models()) {
     stop(
       paste(
         "Hybrid DL currently supports:",
         paste(gp_hybrid_dl_supported_models(), collapse = ", ")
       ),
       call. = FALSE
     )
   }
   res_model_output <- gp_hybrid_dl_gaussian_model(
     model_type = GS_model,
     pheno_object = pheno_data %||% pheno_clean[["pheno_clean_data"]],
       response = response[[1]],
       gen_name = gen_name,
       female_parent = female_parent,
       male_parent = male_parent,
       heter_groups = heter_groups,
       geno_omic_object = gp_hybrid_bind_feature_blocks(
         train_block = ml_dat_res[["merged_data"]][["merge_data"]],
         test_block = ml_dat_res[["merged_data_test"]]
       ),
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
      optimizer_name = optimizer_name,
      max_grad_norm = max_grad_norm,
       n_bootstrap = n_bootstrap,
       dl_internal_calibration = dl_internal_calibration,
       confidence_level = confidence_level %||% 0.95,
      early_stop = early_stop
   )
   res_summary_stat <- summary_statistics_hybrid(
     predicted_object = res_model_output[["predicted_values"]],
     response = response[[1]],
     female_parent = female_parent,
     male_parent = male_parent,
     mode = "hybrid_dl",
     model_type = GS_model
   )
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = gp_detect_python(),
     preferred_python = gp_preferred_python(purpose = "dl"),
     context = "model_execute_hybrid_dl",
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(task_unit = "hybrid_dl")
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = "hybrid_dl",
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model,
      trait_name = response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(multi_trait_ml)) {
   if (isTRUE(cross_validation)) {
     stop("Multi-trait ML currently supports true prediction only.", call. = FALSE)
   }
   if (isTRUE(met_ml_dl)) {
     stop("Multi-trait ML and MET ML/DL are separate paths; combine them in a later implementation.", call. = FALSE)
   }
   if (length(response) < 2L) {
     stop("Multi-trait ML requires at least two response columns.", call. = FALSE)
   }
   if (!identical(response_family, "gaussian")) {
     stop("Multi-trait ML currently supports gaussian traits only.", call. = FALSE)
   }
   if (length(GS_model) != 1L || !GS_model %in% gp_multitrait_ml_supported_models()) {
     stop(
       paste(
         "Multi-trait ML currently supports:",
         paste(gp_multitrait_ml_supported_models(), collapse = ", ")
       ),
       call. = FALSE
     )
   }
   if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
     stop("Multi-trait ML requires processed geno/omic features.", call. = FALSE)
   }
   ml_dat_res <- gp_feature_apply_multitrait_to_ml_dat_res(
     ml_dat_res = ml_dat_res,
     pheno_data = pheno_clean[["pheno_clean_data"]],
     response = response,
     gen_name = gen_name,
     feature_score_metadata = feature_score_metadata,
     k = feature_k
   )

   res_model_output <- gp_multitrait_ml_gaussian_model(
     model_type = GS_model,
     pheno_object = pheno_clean[["pheno_clean_data"]],
     response = response,
     geno_omic_object = gp_ml_merged_feature_matrix(ml_dat_res),
     gen_name = gen_name,
     ntree = ntree,
     mtry = mtry,
     maxnodes = maxnodes,
     nodesize = nodesize,
     rf_n_jobs = rf_n_jobs,
     lambda_rr = lambda_rr,
     ncomp = ncomp,
     svm_kernel = svm_kernel,
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
     internal_cv_nfolds = internal_cv_nfolds,
     internal_cv_replication = internal_cv_replication,
     random_seed = random_seed
   )
   res_summary_stat <- gp_multitrait_ml_summary_statistics(
     pred_long = res_model_output[["predicted_values"]],
     response = response,
     model_type = GS_model
   )
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = NULL,
     preferred_python = NULL,
     context = "model_execute_multi_trait_ml",
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(task_unit = "multi_trait_ml")
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = "multi_trait_ml",
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model,
      trait_name = if (length(response) > 1L) "multi_trait" else response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(multi_trait_asreml)) {
   mt_asreml_source_bank <- selected_joint_kernel_bank(
     "Multi-trait ASReml true prediction"
   ) %||% gmatrix_kernel_model_ready_list
   mt_asreml_kernels <- gp_collect_kernel_inputs(
     gmatrix = mt_asreml_source_bank[["gmatrix_model_ready"]] %||% gmatrix %||% NULL,
     omic1_kernel = mt_asreml_source_bank[["omic1_kernel_model_ready"]] %||% NULL,
     omic2_kernel = mt_asreml_source_bank[["omic2_kernel_model_ready"]] %||% NULL,
     omic3_kernel = mt_asreml_source_bank[["omic3_kernel_model_ready"]] %||% NULL,
     kernel_list = mt_asreml_source_bank
   )
   # cross_validation = TRUE for multi_trait_asreml is now handled by
   # model_execute_route_cross_validation.R via
   # gp_multitrait_asreml_gaussian_cv() with the predict-or-extract
   # fallback chain. Defensive guard kept here in case routing changes:
   if (isTRUE(cross_validation)) {
     stop("Multi-trait ASReml-R cross-validation is routed through model_execute_route_cross_validation.R (set cv_evaluation_only = TRUE).", call. = FALSE)
   }
   if (length(response) < 2L) {
     stop("Multi-trait ASReml-R requires at least two response columns.", call. = FALSE)
   }
   if (!identical(response_family, "gaussian")) {
     stop("Multi-trait ASReml-R currently supports gaussian traits only.", call. = FALSE)
   }
   if (length(GS_model) != 1L || !identical(as.character(GS_model), "GBLUP")) {
     stop("Multi-trait ASReml-R currently supports GS_model = 'GBLUP' only.", call. = FALSE)
   }
   if (length(engine) != 1L || !identical(as.character(engine), "asreml")) {
     stop("Multi-trait ASReml-R requires engine = 'asreml'.", call. = FALSE)
   }
   if (!length(mt_asreml_kernels)) {
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
    if (gp_formula_is_intercept_only(fixed)) fixed <- NULL  # ~1 is the default fixed part
    if (!gp_formula_is_intercept_only(fixed) || !is.null(random) || !is.null(cova)) {
      stop("Dedicated multi-trait ASReml-R routes do not accept custom fixed, random, or cova formulas.", call. = FALSE)
    }
    if (isTRUE(mt_asreml_is_met)) {
      if (is.null(heter_groups) || !nzchar(heter_groups)) {
        stop("Grouped MT-MET ASReml requires heter_groups as the environment column.", call. = FALSE)
      }
      if (!is.null(ctx$weights)) {
        stop(
          "Grouped MT-MET ASReml does not yet accept observation weights; weights were not silently ignored.",
          call. = FALSE
        )
      }
      res_model_output <- gp_multitrait_asreml_mtmet_gaussian_model(
        pheno_object = mt_asreml_pheno,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups,
        gmatrix = mt_asreml_kernels[[1L]],
        kernel_list = if (length(mt_asreml_kernels) > 1L) mt_asreml_kernels[-1L] else NULL,
        heter_resid = heter_resid,
        var_cov_str = var_cov_str,
        inverse = inverse,
        epsilon = epsilon,
        engine = engine,
        workspace = workspace %||% 1e08,
        pworkspace = ctx$pworkspace %||% 1e06,
        maxit = maxit %||% 50L,
        ai_sing = TRUE
      )
    } else {
      trait_covariance <- gp_multitrait_asreml_trait_covariance(var_cov_str)
      res_model_output <- gp_multitrait_asreml_gaussian_model(
        pheno_object = mt_asreml_pheno,
        response = response,
        gen_name = gen_name,
        gmatrix = mt_asreml_kernels[[1L]],
        kernel_list = if (length(mt_asreml_kernels) > 1L) mt_asreml_kernels[-1L] else NULL,
        inverse = inverse,
        epsilon = epsilon,
        engine = engine,
        workspace = workspace %||% 1e08,
        pworkspace = ctx$pworkspace %||% 1e06,
        maxit = maxit %||% 50L,
        trait_covariance = trait_covariance
      )
    }
  res_model_output <- gp_reconcile_multitrait_model_output(
    res_model_output,
    gen_name = gen_name,
    heter_groups = if (isTRUE(mt_asreml_is_met)) heter_groups else NULL,
    model = GS_model,
    response_family = response_family
   )
   mt_asreml_mode <- if (isTRUE(mt_asreml_is_met)) {
     "multi_trait_multi_environment_asreml"
   } else {
     "multi_trait_asreml"
   }
   res_summary_stat <- if (isTRUE(mt_asreml_is_met)) {
     gp_model_execute_multitrait_summary(
       pred_long = res_model_output[["predicted_values"]],
       response = response,
       model_type = GS_model,
       mode = mt_asreml_mode
     )
   } else {
     gp_multitrait_asreml_summary_statistics(
       pred_long = res_model_output[["predicted_values"]],
       response = response,
       model_type = GS_model
     )
   }
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = NULL,
     preferred_python = NULL,
     context = paste0("model_execute_", mt_asreml_mode),
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(task_unit = mt_asreml_mode)
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = mt_asreml_mode,
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model,
      trait_name = if (length(response) > 1L) "multi_trait" else response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }

 if (isTRUE(multi_trait_dl)) {
   if (isTRUE(cross_validation)) {
     stop("Multi-trait DL is currently implemented for true prediction only; cross-validation will be added next.", call. = FALSE)
   }
   if (isTRUE(met_ml_dl)) {
     stop("Multi-trait DL and MET ML/DL are separate paths; combine them in a later implementation.", call. = FALSE)
   }
   if (length(response) < 2L) {
     stop("Multi-trait DL requires at least two response columns.", call. = FALSE)
   }
   if (!identical(response_family, "gaussian")) {
     stop("Multi-trait DL currently supports gaussian traits only.", call. = FALSE)
   }
   if (length(GS_model) != 1L || !GS_model %in% gp_multitrait_dl_supported_models()) {
     stop(
       paste(
         "Multi-trait DL currently supports:",
         paste(gp_multitrait_dl_supported_models(), collapse = ", ")
       ),
       call. = FALSE
     )
   }
   if (is.null(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
     stop("Multi-trait DL requires processed geno/omic features.", call. = FALSE)
   }
   ml_dat_res <- gp_feature_apply_multitrait_to_ml_dat_res(
     ml_dat_res = ml_dat_res,
     pheno_data = pheno_clean[["pheno_clean_data"]],
     response = response,
     gen_name = gen_name,
     feature_score_metadata = feature_score_metadata,
     k = feature_k
   )

   res_model_output <- gp_multitrait_dl_gaussian_model(
     model_type = GS_model,
     pheno_object = pheno_clean[["pheno_clean_data"]],
     response = response,
     geno_omic_object = gp_ml_merged_feature_matrix(ml_dat_res),
     gen_name = gen_name,
     message = message,
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
      internal_cv_nfolds = internal_cv_nfolds,
      internal_cv_replication = internal_cv_replication,
      dl_internal_calibration = dl_internal_calibration,
      system_database = system_database
   )
   res_summary_stat <- gp_multitrait_dl_summary_statistics(
     pred_long = res_model_output[["predicted_values"]],
     response = response,
     model_type = GS_model
   )
   run_metadata <- gp_runtime_metadata(
     execution_policy = NULL,
     python_path = gp_detect_python(),
     preferred_python = gp_preferred_python(purpose = "dl"),
     context = "model_execute_multi_trait_dl",
     extra_fields = c(
       gp_runtime_profile_metadata_fields(run_profile),
       list(task_unit = "multi_trait_dl")
     )
   )
   handled_result <- results_handling(
     GS_model = GS_model,
     res_model_output = res_model_output,
     res_summary_stat = res_summary_stat,
     res_plot = NULL,
     res_plot_mean = NULL,
     res_plot_result_diagnostic = NULL,
     test_diagonistic_plots = res_model_output[["diagnostic_plots"]],
     res_mod_results_cv_per_trait_model = NULL,
     res_plot_result_diagnostic_cv_only = NULL,
     geno_qc_stat = geno_qc_stat,
     cv_results_processed = NULL,
     run_metadata = run_metadata,
     system_database = system_database,
     plot_filename = paste(response, collapse = "_"),
     plot_extension = plot_extension,
     plot_width = plot_width,
     plot_height = plot_height,
     plot_units = plot_units,
     plot_dpi = plot_dpi,
     feature_selected = ml_dat_res[["feature_selected"]] %||% feature_selected
   )
   return(gp_index_specialized_model_execute_result(
     handled_result = handled_result,
      model_label = GS_model,
      trait_name = if (length(response) > 1L) "multi_trait" else response[[1]],
      gen_name = gen_name,
      heter_groups = heter_groups,
      run_metadata = run_metadata,
     feature_score_metadata = feature_score_metadata
   ))
 }


  NULL
}
