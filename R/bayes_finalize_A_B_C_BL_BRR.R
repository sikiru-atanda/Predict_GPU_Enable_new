#' Compile Bayesian Genomic Selection Models BRR, BayesA, BayesB, BayesC, BL
#'
#' This function compiles inputs for Bayesian genomic selection models, executes the model,
#' and processes the output. It supports a range of genomic selection (GS) models, incorporating
#' fixed and random effects, multiple omic data sources, and specific Bayesian parameters.
#'
#' @param fixed Formula or character vector specifying the fixed effects to be included in the model.
#' @param random Formula or character vector specifying the random effects to be included in the model.
#' @param GS_model Character string specifying the type of genomic selection model to be used.
#' @param response Character string specifying the response variable in the phenotypic data.
#' @param weights Optional positive Stage 2 observation precisions, supplied as
#'   a numeric vector, one-column table, or phenotype column name. Gaussian BGLR
#'   fits receive their square roots to match BGLR's inverse-squared convention.
#' @param fixed_term_model_bayesian Character vector specifying the fixed terms in the Bayesian model.
#' @param rand_term_model_bayesian Character vector specifying the random terms in the Bayesian model.
#' @param pheno_data Data frame containing the phenotypic data.
#' @param geno_data Matrix or data frame containing the genotypic data.
#' @param omic1_data Matrix or data frame containing the first set of omic data (optional).
#' @param omic2_data Matrix or data frame containing the second set of omic data (optional).
#' @param omic3_data Matrix or data frame containing the third set of omic data (optional).
#' @param gen_name Character string specifying the name of the genotype factor in the genotypic data.
#' @param nIter Integer specifying the number of iterations for the Bayesian model.
#' @param burnIn Integer specifying the number of burn-in iterations for the Bayesian model.
#' @param thin Integer specifying the thinning interval for the Bayesian model.
#' @param omics_data_label Character vector specifying labels for the omic data sets (optional).
#' @param cross_validation Logical; when `TRUE`, the call is part of a
#'   cross-validation fold (CV bookkeeping is adjusted accordingly).
#' @param scaling Logical; if `TRUE`, centre and scale the marker / omics
#'   matrices before fitting (passed through to the ETA compiler).
#' @param response_family Response family: `"gaussian"`, `"binary"` (or
#'   alias `"binomial"`), or `"ordinal"`. Unordered nominal/multiclass outcomes
#'   are not supported by BGLR.
#' @param ... Additional arguments for future extensions.
#'
#' @return A list containing two elements: `bayes_result`, which holds the processed output from the
#' Bayesian model, and `bayes_model`, which contains the raw model object from the Bayesian analysis.
#'
#' @details
#' The function integrates various steps required for running Bayesian genomic selection models,
#' from compiling necessary data and model specifications to executing the model and formatting the output.
#' It is designed to work with multiple types of genomic selection models and can handle complex scenarios
#' involving multiple sources of omic data and varying model parameters.
#'
#' This is a complex function designed for specific genomic selection workflows
#' and requires extensive setup for data and model parameters. Example usage:

#' @export
#'
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
                                        omics_data_label = NULL,
                                        cross_validation = FALSE,
                                        scaling =TRUE,
                                        response_family = "gaussian",
                                        ...
                                             ) {

            dots <- list(...)

            valid_marker_bayes_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
            if (length(GS_model) != 1L || !(GS_model %in% valid_marker_bayes_models)) {
              stop(
                paste(
                  "Bayesian marker finalization supports GS_model values",
                  paste(valid_marker_bayes_models, collapse = ", "),
                  "only."
                ),
                call. = FALSE
              )
            }

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
                                      omics_data_label = omics_data_label,
                                      scaling = scaling)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)
            if(isTRUE(cross_validation)){
            return(list(bayes_ETA = ETA,
                        bayes_para = bayes_para))

            }

            mod <-  bayes_mod_execute(pheno_data = ETA[["pheno_data"]],
                                      response = response,
                                      weights = weights,
                                      ETA = ETA[["ETA"]],
                                      bayes_para = bayes_para,
                                      verbose = FALSE,
                                      GS_model = GS_model,
                                      response_family = response_family
                                      #files_key = "files_key"
            )


            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                  pheno_data = ETA[["pheno_data"]],
                                                  geno_data = geno_data,
                                                  gen_name = gen_name,
                                                  response = response,
                                                 omic1_data = omic1_data,
                                                 omic2_data = omic2_data,
                                                 omic3_data = omic3_data,
                                                 omics_data_label = omics_data_label,
                                                 bayes_para = bayes_para,
                                                 GS_model = GS_model,
                                                 CI_width_thresholds = dots[["CI_width_thresholds"]] %||% c(0.33, 0.66),
                                                 confidence_level = dots[["confidence_level"]] %||% 0.95,
                                                 high_reliability_thres = dots[["high_reliability_thres"]] %||% 0.9,
                                                 low_reliability_thres = dots[["low_reliability_thres"]] %||% 0.5,
                                                 n_components = dots[["n_components"]] %||% 20,
                                                 threshold = dots[["threshold"]] %||% 100,
                                                 target = dots[["target"]] %||% "test_set",
                                                 interval_width_high_threshold = dots[["interval_width_high_threshold"]] %||% NULL,
                                                 interval_width_low_threshold = dots[["interval_width_low_threshold"]] %||% NULL,
                                                 interval_width_moderate_threshold = dots[["interval_width_moderate_threshold"]] %||% NULL,
                                                 response_family = response_family)


            return(list(bayes_result = res_model_output,
                        bayes_model = mod))



}
