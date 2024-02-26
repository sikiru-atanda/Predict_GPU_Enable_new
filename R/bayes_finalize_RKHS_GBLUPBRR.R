

#' Title
#'
#' @param fixed
#' @param random
#' @param GS_model
#' @param response
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param pheno_data
#' @param gkernel
#' @param gmatrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param gen_name
#' @param nIter
#' @param burnIn
#' @param thin
#' @param heter_groups
#' @param omics_kernel_label
#' @param ...
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
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
                                    gen_name = gen_name)



bayes_para <-  bayes_parameter_check(nIter = nIter,
                                     burnIn = burnIn,
                                     thin = thin)

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

output <-  list(bayes_result = res_model_output,
              bayes_model = mod)

return(output)


}
