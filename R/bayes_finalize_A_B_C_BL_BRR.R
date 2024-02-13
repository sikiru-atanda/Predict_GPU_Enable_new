

#' Title
#'
#' @param fixed
#' @param random
#' @param GS_model
#' @param response
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param pheno_data
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param gen_name
#' @param nIter
#' @param burnIn
#' @param thin
#' @param omics_data_label
#'
#' @return
#' @export
#'
#' @examples
 bayes_finalize_A_B_C_BL_BRR <- function(fixed = NULL,
                                        random = NULL,
                                        GS_model = NULL,
                                        response = NULL,
                                        weights = NULL,
                                        fixed_term_model_bayesian = NULL,
                                        rand_term_model_bayesian = NULL,
                                        pheno_data = NULL,
                                        geno_data = NULL,
                                        omic1_data = NULL,
                                        omic2_data = NULL,
                                        omic3_data = NULL,
                                        gen_name = NULL,
                                        nIter = NULL,
                                        burnIn = NULL,
                                        thin = NULL,
                                        omics_data_label = NULL
                                             )
 {

            ETA  <-  ETA_compiler_bayes(
                                      fixed = fixed,
                                      random = random,
                                      GS_model = GS_model,
                                      fixed_term_model_bayesian = fixed_term_model_bayesian,
                                      rand_term_model_bayesian = rand_term_model_bayesian,
                                      pheno_data = pheno_data,
                                      geno_data = geno_data,
                                      omic1_data = omic1_data,
                                      omic2_data = omic2_data,
                                      omic3_data = omic3_data,
                                      gen_name = gen_name,
                                      omics_data_label = omics_data_label)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                      #files_key = "files_key"
            )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_data,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_data,
                                                 omic2_data = omic2_data,
                                                 omic3_data = omic3_data,
                                                 omics_data_label = omics_data_label,
                                                 bayes_para = bayes_para,
                                                 GS_model = GS_model)


            output = list(res_model_output, mod)

            return(output)

}
