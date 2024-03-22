#
#
# geno_qc_stat <- if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL
#
# #is this correct?
# if(!is.null(best_models)){
#   n_trait <- length(best_models[["trait"]])
#   n_model <- length(best_models[["model"]])
#
# } else {
#   if (length(GS_model) > 1) {
#     if (length(GS_model) != length(response)) {
#       stop("When the number of models is more than one, the number of models should be the same as the number of traits.")
#     }
#   }
#
#   task <- data.frame(model = GS_model, trait = response, stringsAsFactors = FALSE)
#   n_trait <-  length(response)
#   n_model <- length(GS_model)
# }
# # Main logic
# sys_name <- Sys.info()["sysname"]
# if (!is.null(num_cores) && num_cores > 1) {
#   #sys_name <- Sys.info()["sysname"]
#   set_parallel_plan(n_trait = n_trait, n_model = n_model,
#                     sys_name = sys_name)
# } else {
#   # Automatically determine the number of cores and use half of them
#   detected_cores <- parallel::detectCores(logical = TRUE)
#   # For non-Windows systems, consider physical cores only
#   num_cores <- round(detected_cores * 0.5)
#
#   set_parallel_plan(n_trait, n_model,num_cores,sys_name)
# }
#
# #results_use = results
#
# results <- future.apply::future_lapply(seq_len(nrow(best_models)), function(i) {
#   task_row <- best_models[i, ]
#
#   response <- as.character(task_row$trait)
#   GS_model <- as.character(task_row$model)
#
#   #### These models only works with one environment/location
#   if(length(pheno_clean[["pheno_clean_data"]][,gen_name])==length(unique(pheno_clean[["pheno_clean_data"]][,gen_name]))){
#
#     if ((GS_model %in% bayes_valid_models && is.null(rand_term_model_bayesian)) ||
#         (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models)) ||
#         (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models))) {
#
#       res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
#                                                       random = random,
#                                                       GS_model = GS_model,
#                                                       response = response,
#                                                       weights = weights,
#                                                       fixed_term_model_bayesian = fixed_term_model_bayesian,
#                                                       rand_term_model_bayesian = rand_term_model_bayesian,
#                                                       pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                       geno_data = if("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
#                                                       omic1_data = if("omic1_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic1_model_ready"]] else NULL,
#                                                       omic2_data = if("omic2_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic2_model_ready"]] else NULL,
#                                                       omic3_data = if("omic3_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic3_model_ready"]] else NULL,
#                                                       gen_name = gen_name,
#                                                       nIter = nIter,
#                                                       burnIn = burnIn,
#                                                       thin = thin,
#                                                       omics_data_label = omics_data_label,
#                                                       scaling = scaling)
#
#       # Compute summary statistics and plot accuracy
#       res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
#       #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
#       res_model_output <- res_model_output[["bayes_result"]]
#
#       # output <- list(GS_model = GS_model,
#       #                res_model_output = res_model_output,
#       #                res_summary_stat = res_summary_stat,
#       #                geno_qc_stat =geno_qc_stat
#       # )
#
#     }
#   } ## End of  Bayes A, B, C, BRR, BL
#
#   ##########################################################################
#   #########################################################################
#   ## Start of Reproducing Kernel Hilbert Spaces Regression RKHS,         ##
#   ## (BRR- Bayesian GBLUP ) and GBLUP (asreml) Model                     ##
#   ## for Single Location and multiple loc                                ##
#   ##                                                                     ##
#   ##                                                                     ##
#   ##########################################################################
#   #######################################################################
#
#   ## NOTE
#   ## BRR is changed to G-BRR to make distinction between BRR for marker matrix and GBLUP
#   # Check conditions for GS_model and rand_term_model_bayesian
#   if ((GS_model %in% c("RKHS", "GBLUP_BRR", "GBLUP") && is.null(rand_term_model_bayesian)) ||
#       (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models)) ||
#       (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models))) {
#
#     # Rename GS_model for GBLUP_BRR case
#     if (GS_model == "GBLUP_BRR") {
#       GS_modeluse <- GS_model
#       GS_model <- "BRR"
#     }
#
#     if (GS_model %in% c("BRR", "RKHS")) {
#       # Run Bayesian model for BRR and RKHS
#       res_model_output <- bayes_finalize_RKHS_GBLUPBRR(fixed = fixed,
#                                                        random = random,
#                                                        GS_model = GS_model,
#                                                        response = response,
#                                                        weights = weights,
#                                                        fixed_term_model_bayesian = fixed_term_model_bayesian,
#                                                        rand_term_model_bayesian = rand_term_model_bayesian,
#                                                        pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                        gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#                                                        omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#                                                        omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#                                                        omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#                                                        gen_name = gen_name,
#                                                        nIter = nIter,
#                                                        burnIn = burnIn,
#                                                        thin = thin,
#                                                        heter_groups = heter_groups,
#                                                        omics_kernel_label = omics_kernel_label)
#
#       # Compute summary statistics and plot accuracy
#       res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
#       #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
#       res_model_output <- res_model_output[["bayes_result"]]
#       ### This part is for GBLUP_BRR
#       # if(exists("GS_modeluse")){
#       #   GS_model <-  GS_modeluse
#       # }
#
#       # output <- list(GS_model = GS_model,
#       #                res_model_output = res_model_output,
#       #                res_summary_stat = res_summary_stat,
#       #                geno_qc_stat =geno_qc_stat
#       # )
#
#     } else {
#       if (GS_model == "GBLUP" && engine == 'asreml') {
#
#         #  # Run GBLUP model with ASReml
#         mod <- asreml_utilis_new(fixed = fixed,
#                                  random = random,
#                                  cova = cova,
#                                  GS_model = GS_model,
#                                  response = response,
#                                  pheno_data = pheno_clean[["pheno_clean_data"]],
#                                  gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#                                  omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#                                  omic2_kernel = if("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#                                  omic3_kernel = if("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#                                  gen_name = gen_name,
#                                  heter_groups = heter_groups,
#                                  heter_resid = heter_resid,
#                                  var_cov_str = var_cov_str,
#                                  weights = weights,
#                                  pworkspace = pworkspace,
#                                  workspace = workspace,
#                                  maxit = maxit,
#                                  inverse = inverse,
#                                  epsilon = epsilon,
#                                  engine = engine)
#
#         # Extract model output for ASReml
#
#         res_model_output <- asreml_mod_output_new(
#           mod_asreml = mod,
#           pheno_data = pheno_clean[["pheno_clean_data"]],
#           gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#           omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#           omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#           omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#           heter_groups = heter_groups,
#           gen_name = gen_name,
#           response = response,
#           var_cov_str = var_cov_str,
#           heter_resid = heter_resid,
#           pworkspace = pworkspace,
#           workspace = workspace,
#           maxit = maxit
#         )
#
#
#         res_summary_stat <- summary_statistics_asreml(mod =  res_model_output[["Asreml_model"]],
#                                                       response = response,
#                                                       pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                       heter_groups = heter_groups,
#                                                       predicted_value =  res_model_output[["Predicted_value"]],
#                                                       pred_heter_groups = NULL,
#                                                       variance_components = res_model_output[["Variance_components"]],
#                                                       eval_metrics = eval_metrics)
#
#         # output <- list(GS_model = GS_model,
#         #                res_model_output = res_model_output,
#         #                res_summary_stat = res_summary_stat,
#         #                geno_qc_stat =geno_qc_stat
#         # )
#
#       }
#     }
#   }
#   #### END GBLUP_RKHS, GBLUP_BRR and GBLUP (asreml)
#
#   ######################################################
#   ######################################################
#   ##                                                 ###
#   ## Machine Learning Models                         ###
#   ##                                                 ###
#   ######################################################
#   ######################################################
#
#   if (GS_model %in% AI_valid_models) {
#     if (length(unique(pheno_clean[["pheno_clean_data"]][, gen_name])) > length(pheno_clean[["pheno_clean_data"]][, gen_name])) {
#       stop(paste(msg, GS_model, 'only works for single location/enviroment.'), call. = FALSE)
#     }
#
#     switch(GS_model,
#            "Xgboost" = {
#              res_model_output <- AI_Xgb(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                         response = response,
#                                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                         message = message,
#                                         gen_name = gen_name,
#                                         scaling = scaling,
#                                         centering = centering,
#                                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                         AI_cv_nfolds = AI_cv_nfolds,
#                                         para_tunning = para_tunning,
#                                         xgb_paras_tunning = xgb_paras_tunning
#              )
#            },
#            "RandomForest" = {
#              res_model_output <- AI_randomForest(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                  response = response,
#                                                  geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                  geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                  message = message,
#                                                  gen_name = gen_name,
#                                                  scaling = scaling,
#                                                  centering = centering,
#                                                  omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                  AI_cv_nfolds = AI_cv_nfolds,
#                                                  para_tunning = para_tunning,
#                                                  rf_paras_tunning = rf_paras_tunning
#              )
#            },
#            "PartialLeastSquare" = {
#              res_model_output <-  AI_pls(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                          response = response,
#                                          geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                          geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                          message = message,
#                                          gen_name = gen_name,
#                                          scaling = scaling,
#                                          centering = centering,
#                                          omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                          para_tunning = para_tunning,
#                                          pls_paras_tunning = pls_paras_tunning)
#            },
#            "SupportVectorMachine" = {
#              res_model_output <- AI_svm(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                         response = response,
#                                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                         message = message,
#                                         gen_name = gen_name,
#                                         scaling = scaling,
#                                         centering = centering,
#                                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                         AI_cv_nfolds = AI_cv_nfolds,
#                                         para_tunning = para_tunning,
#                                         svm_paras_tunning = svm_paras_tunning
#              )
#            },
#            "K-NearestNeighbors" = {
#              res_model_output <- AI_knn(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                         response = response,
#                                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                         message = message,
#                                         gen_name = gen_name,
#                                         scaling = scaling,
#                                         centering = centering,
#                                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                         AI_cv_nfolds = AI_cv_nfolds,
#                                         para_tunning = para_tunning,
#                                         knn_paras_tunning = knn_paras_tunning
#              )
#            },
#            "Lasso" = {
#              res_model_output <- AI_RidgeRegression_Lasso(
#                pheno_object = ml_dat_res[["pheno_clean_data"]],
#                response = response,
#                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                gen_name = gen_name,
#                para_tunning = para_tunning,
#                AI_cv_nfolds = AI_cv_nfolds,
#                lasso_paras_tunning = lasso_paras_tunning,
#                message = message,
#                scaling = scaling,
#                centering = centering,
#                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                GS_model = GS_model
#              )
#            },
#            "Ridge_Regression" = {
#              res_model_output <- AI_RidgeRegression_Lasso(
#                pheno_object = ml_dat_res[["pheno_clean_data"]],
#                response = response,
#                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                gen_name = gen_name,
#                para_tunning = para_tunning,
#                AI_cv_nfolds = AI_cv_nfolds,
#                lasso_paras_tunning = rr_paras_tunning,
#                message = message,
#                scaling = scaling,
#                centering = centering,
#                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                GS_model = GS_model,
#              )
#            },
#            "deep_learning_model" = {
#
#              res_model_output <- deep_learning_model(
#                pheno_object=ml_dat_res[["pheno_clean_data"]],
#                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                response=response,
#                gen_name=gen_name,
#                num_hidden_layers = num_hidden_layers,
#                neurons_per_layer = neurons_per_layer,
#                learning_rate = learning_rate,
#                epochs = epochs,
#                batch_size = batch_size,
#                validation_split = validation_split,
#                early_stop = early_stop,
#                message = message,
#                scaling = scaling,
#                centering = centering,
#                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                para_tunning = para_tunning,
#                param_grid = dpl_paras_tunning
#              )
#
#            },
#            {
#              stop(paste(msg, "Select method to calculate geno_cleanmic relationship matrix"), call. = FALSE)
#            }
#     )
#
#     ##browser()
#     #View(res_model_output[["predicted_values"]])
#     res_summary_stat <- summary_statistics_AI(predicted_object = res_model_output[["predicted_values"]],
#                                               pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                               response = response,
#                                               test_set = ml_dat_res[["test_set"]],
#                                               geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                               eval_metrics = eval_metrics,
#                                               model_parameters = res_model_output[["model_parameters"]],
#                                               GS_model = GS_model
#     )
#
#     # output <- list(GS_model = GS_model,
#     #                res_model_output = res_model_output,
#     #                res_summary_stat = res_summary_stat,
#     #                geno_qc_stat = geno_qc_stat
#     # )
#
#   }
#   ### End machine learning
# list(GS_model = GS_model,res_model_output = res_model_output,
#      res_summary_stat = res_summary_stat, geno_qc_stat = geno_qc_stat)
#   #list(output =  output)
# }, future.seed = TRUE)
#
#
#
# if(!is.null(best_models)){
# names(results) <- best_models[["trait"]]
#
# } else {
#   names(results) <- response
# }
# ########
# best_models <- cv_results_processed[["best_models_list"]][[metric_for_ranking]]
#
# best_models_ggplot_rep <- cv_results_processed[["plot_reps_list"]][[metric_for_ranking]][["ggplot_boxplot_reps"]]
#
# best_models_ggplot_mean <- cv_results_processed[["plot_mean_list"]][[metric_for_ranking]][["ggplot_lineplot_mean"]]
#
#
# for (res in 1:length(results)) {
#
# #names(results[[1]])
#
#
# results_handling(GS_model = if("GS_model" %in% names(results[[res]])) results[[res]][["GS_model"]] else NULL,
#                         res_model_output = if("res_model_output" %in% names(results[[res]])) results[[res]][["res_model_output"]] else NULL,
#                         res_summary_stat = if("res_summary_stat" %in% names(results[[res]])) results[[res]][["res_summary_stat"]] else NULL,
#                         res_plot = best_models_ggplot_rep,
#                         res_plot_mean = best_models_ggplot_mean,
#                         geno_qc_stat = geno_qc_stat,
#                         cv_results_processed = cv_results_processed,
#                         system_database = system_database,
#                         plot_filename = if(!is.null(names(results)[res])) names(results)[res] else paste("trait", res, sep = "_"),
#                         plot_extension = plot_extension,
#                         plot_width = plot_width,
#                         plot_height = plot_height,
#                         plot_units = plot_units,
#                         plot_dpi = plot_dpi)
#
# }
