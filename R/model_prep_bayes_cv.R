

#' Title
#'
#' @param fixed
#' @param random
#' @param GS_model_cv
#' @param response
#' @param gen_name
#' @param pheno_data
#' @param weights
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param test_set
#' @param nIter
#' @param burnIn
#' @param thin
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param omics_data_label
#' @param gmatrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param heter_groups
#' @param omics_kernel_label
#' @param cross_validation
#'
#' @return
#' @export
#'
#' @examples
model_prep_bayes_cv<- function(fixed = NULL,
                               random = NULL,
                               GS_model_cv = NULL,
                               response = NULL,
                               gen_name = NULL,
                               pheno_data = NULL,
                               weights = NULL,
                               fixed_term_model_bayesian = NULL,
                               rand_term_model_bayesian = NULL,
                               test_set = NULL,
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
  # if(!is.null(test_set)){
  # if(is.data.frame(test_set) | is.matrix(test_set)){
  #   test_set <-  test_set[, 1]
  #  }
  #
  # }

  for (model in GS_model_cv) {

  if (model %in% bayes_valid_models) {
  res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
                                                  random = random,
                                                  GS_model = model,
                                                  response = response,
                                                  weights = weights,
                                                  fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                                  #pheno_data = if(!is.null(test_set)) pheno_data[pheno_data[[gen_name]] %in% test_set, ] else pheno_data,
                                                  pheno_data = pheno_data,
                                                  geno_data = if(!is.null(test_set)) geno_data[rownames(geno_data) %in% test_set, ] else geno_data,
                                                  omic1_data = if(!is.null(test_set)) omic1_data[rownames(omic1_data) %in% test_set, ] else omic1_data,
                                                  omic2_data = if(!is.null(test_set)) omic2_data[rownames(omic2_data) %in% test_set, ] else omic2_data,
                                                  omic3_data = if(!is.null(test_set)) omic3_data[rownames(omic3_data) %in% test_set, ] else omic3_data,
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
                                                  #pheno_data = if(!is.null(test_set)) pheno_data[pheno_data[[gen_name]] %in% test_set, ] else pheno_data,
                                                  pheno_data = pheno_data,
                                                  gmatrix = if(!is.null(test_set)) gmatrix[rownames(gmatrix) %in% test_set, colnames(gmatrix) %in% test_set] else gmatrix,
                                                  omic1_kernel = if(!is.null(test_set)) omic1_kernel[rownames(omic1_kernel) %in% test_set, colnames(omic1_kernel) %in% test_set] else omic1_kernel,
                                                  omic2_kernel = if(!is.null(test_set)) omic2_kernel[rownames(omic2_kernel) %in% test_set, colnames(omic2_kernel) %in% test_set] else omic2_kernel,
                                                  omic3_kernel = if(!is.null(test_set)) omic3_kernel[rownames(omic3_kernel) %in% test_set, colnames(omic3_kernel) %in% test_set] else omic3_kernel,
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
