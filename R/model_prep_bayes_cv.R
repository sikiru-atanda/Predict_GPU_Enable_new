

model_prep_bayes_cv<- function(fixed = NULL,
                                   random = NULL,
                                   GS_model_cv = NULL,
                                   response = NULL,
                                   gen_name = NULL,
                                   pheno_data = NULL,
                                   weights = NULL,
                                   fixed_term_model_bayesian = NULL,
                                   rand_term_model_bayesian = NULL,
                                   nIter = NULL,
                                   burnIn = NULL,
                                   thin = NULL,
                                   geno_data = NULL,
                                   omic1_data = NULL,
                                   omic2_data = NULL,
                                   omic3_data = NULL,
                                   omics_data_label = NULL,
                                   gmatrix = NULL,
                                   omic1_kernel = NULL,
                                   omic2_kernel = NULL,
                                   omic3_kernel = NULL,
                                   heter_groups = NULL,
                                   omics_kernel_label = NULL,
                                   cross_validation = TRUE
                                   ){
#browser()
  bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
  res_model_output_list <- list()

  for (model in GS_model_cv) {

  if (model %in% bayes_valid_models) {
  res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
                                                  random = random,
                                                  GS_model = model,
                                                  response = response,
                                                  weights = weights,
                                                  fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                                  pheno_data = pheno_data,
                                                  geno_data = geno_data,
                                                  omic1_data = omic1_data,
                                                  omic2_data = omic2_data,
                                                  omic3_data = omic3_data,
                                                  gen_name = gen_name,
                                                  nIter = nIter,
                                                  burnIn = burnIn,
                                                  thin = thin,
                                                  omics_data_label = omics_data_label,
                                                  cross_validation = cross_validation)
  res_model_output_list[[model]] <-  res_model_output
  }

  if (model %in% bayes_gblup_valid_models) {

    if (model == "GBLUP_BRR") {
      model <- "BRR"
    }
  res_model_output <- bayes_finalize_RKHS_GBLUPBRR(fixed = fixed,
                                                  random = random,
                                                  GS_model = model,
                                                  response = response,
                                                  weights = weights,
                                                  fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                                  pheno_data = pheno_data,
                                                  gmatrix = gmatrix,
                                                  omic1_kernel = omic1_kernel,
                                                  omic2_kernel = omic2_kernel,
                                                  omic3_kernel = omic3_kernel,
                                                  gen_name = gen_name,
                                                  heter_groups = heter_groups,
                                                  nIter = nIter,
                                                  burnIn = burnIn,
                                                  thin = thin,
                                                  omics_kernel_label = omics_kernel_label,
                                                  cross_validation = cross_validation)


  res_model_output_list[[model]] <-  res_model_output
  }

  }

 return(res_model_output_list)

}
