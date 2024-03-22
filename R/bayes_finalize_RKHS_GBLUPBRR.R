

#' Compile Bayesian RKHS, GBLUP_BRR Models
#'
#' This function prepares inputs for, and executes, Bayesian genomic selection models specifically designed for RKHS (Reproducing Kernel Hilbert Spaces), GBLUP (Genomic Best Linear Unbiased Prediction), and BRR (Bayesian Ridge Regression) approaches. It incorporates fixed and random effects, kernel matrices for genomic and omics data, and handles specific Bayesian parameters.
#'
#' @param fixed Formula or character vector specifying the fixed effects to be included in the model.
#' @param random Formula or character vector specifying the random effects to be included in the model.
#' @param GS_model Character string specifying the type of genomic selection model to be used: "RKHS", "GBLUP", or "BRR".
#' @param response Character string specifying the response variable in the phenotypic data.
#' @param weights Vector of weights for the observations in the model (optional).
#' @param fixed_term_model_bayesian Character vector specifying the fixed terms in the Bayesian model.
#' @param rand_term_model_bayesian Character vector specifying the random terms in the Bayesian model.
#' @param pheno_data Data frame containing the phenotypic data.
#' @param gkernel Kernel matrix for genomic data (optional).
#' @param gmatrix Genomic relationship matrix (optional).
#' @param omic1_kernel Kernel matrix for the first omics data type (optional).
#' @param omic2_kernel Kernel matrix for the second omics data type (optional).
#' @param omic3_kernel Kernel matrix for the third omics data type (optional).
#' @param gen_name Character string specifying the name of the genotype factor in the genomic data.
#' @param nIter Integer specifying the number of iterations for the Bayesian model.
#' @param burnIn Integer specifying the number of burn-in iterations for the Bayesian model.
#' @param thin Integer specifying the thinning interval for the Bayesian model.
#' @param heter_groups Character string identifying the grouping variable for heterogeneity (optional).
#' @param omics_kernel_label Character vector specifying labels for the omic data sets (optional).
#' @param ... Additional arguments for future extensions.
#'
#' @return A list containing two elements: `bayes_result`, which holds the processed output from the Bayesian model, and `bayes_model`, which contains the raw model object from the Bayesian analysis.
#'
#' @details
#' The function integrates steps required for running specific Bayesian genomic selection models, compiling necessary data and model specifications to execute the model and format the output. Designed to work with RKHS, GBLUP, and BRR approaches, it can handle complex scenarios involving multiple sources of omic data and varying model parameters.
#'
#' @export
#'

bayes_finalize_RKHS_GBLUPBRR <- function(fixed = NULL,
                                         random = NULL,
                                         GS_model = NULL,
                                         response = NULL,
                                         weights = NULL,
                                         fixed_term_model_bayesian = NULL,
                                         rand_term_model_bayesian = NULL,
                                         pheno_data = NULL,
                                         gkernel =  NULL,
                                         gmatrix = NULL,
                                         omic1_kernel = NULL,
                                         omic2_kernel = NULL,
                                         omic3_kernel = NULL,
                                         gen_name = NULL,
                                         nIter = NULL,
                                         burnIn = NULL,
                                         thin = NULL,
                                         heter_groups  = NULL,
                                         omics_kernel_label = NULL,
                                         cross_validation = FALSE,
                                         ...){


  ETA  <-  ETA_compiler_bayes_GBLUP(
                                    fixed = fixed,
                                    random = random,
                                    GS_model = GS_model,
                                    fixed_term_model_bayesian = fixed_term_model_bayesian,
                                    rand_term_model_bayesian = rand_term_model_bayesian,
                                    pheno_data = pheno_data,
                                    gkernel =  gkernel,
                                    gmatrix = gmatrix,
                                    omic1_kernel = omic1_kernel,
                                    omic2_kernel = omic2_kernel,
                                    omic3_kernel = omic3_kernel,
                                    heter_groups = heter_groups,
                                    gen_name = gen_name)


bayes_para <-  bayes_parameter_check(nIter = nIter,
                                     burnIn = burnIn,
                                     thin = thin)

if(isTRUE(cross_validation)){
  return(list(bayes_para = bayes_para,
              bayes_ETA = ETA))
}

mod <-  bayes_mod_execute(pheno_data = ETA[["pheno_data"]],
                          response = response,
                          weights = weights,
                          ETA = ETA[["ETA"]],
                          bayes_para = bayes_para
)


  if(GS_model%in%c("RKHS", "BRR")){
    # mod_output_bayes_RKHS
    res_model_output <- mod_output_bayes_gbluBRR_RKHS(mod = mod,
                                              ETA = ETA,
                                              gkernel =  gkernel,
                                              GS_model = GS_model,
                                              gmatrix = gmatrix,
                                              omic1_kernel = omic1_kernel,
                                              omic2_kernel = omic2_kernel,
                                              omic3_kernel = omic3_kernel,
                                              gen_name = gen_name,
                                              pheno_data = ETA[["pheno_data"]],
                                              heter_groups  = heter_groups,
                                              omics_kernel_label = omics_kernel_label,
                                              bayes_para = bayes_para
    )

  }

  # if(GS_model=="BRR"){
  #   # mod_output_bayes_BRRGBLUP
  #   res_model_output <- mod_output_bayes_gbluBRR_RKHS(mod = mod,
  #                                                 ETA = ETA,
  #                                                 gkernel =  gkernel,
  #                                                 gmatrix = gmatrix,
  #                                                 GS_model = GS_model,
  #                                                 omic1_kernel = omic1_kernel,
  #                                                 omic2_kernel = omic2_kernel,
  #                                                 omic3_kernel = omic3_kernel,
  #                                                 gen_name = gen_name,
  #                                                 pheno_data = ETA[["pheno_data"]],
  #                                                 heter_groups  = heter_groups,
  #                                                 omics_kernel_label = omics_kernel_label,
  #                                                 bayes_para = bayes_para
  #   )
  #
  # }

 return(list(bayes_result = res_model_output,
              bayes_model = mod))


}
