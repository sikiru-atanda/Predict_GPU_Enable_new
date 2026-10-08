

#' Compile Bayesian RKHS, GBLUP_BRR Models
#'
#' This function prepares inputs for, and executes, Bayesian genomic selection models specifically designed for RKHS (Reproducing Kernel Hilbert Spaces), GBLUP_BRR, and BRR (Bayesian Ridge Regression) approaches. It incorporates fixed and random effects, kernel matrices for genomic and omics data, and handles specific Bayesian parameters.
#'
#' @param fixed Formula or character vector specifying the fixed effects to be included in the model.
#' @param random Formula or character vector specifying the random effects to be included in the model.
#' @param GS_model Character string specifying the type of genomic selection model to be used: "RKHS", "GBLUP_BRR", or "BRR". Use "GBLUP" through the ASReml/GBLUP route, not this Bayesian finalizer.
#' @param response Character string specifying the response variable in the phenotypic data.
#' @param weights Optional positive Stage 2 observation precisions, supplied as
#'   a numeric vector, one-column table, or phenotype column name. Gaussian BGLR
#'   fits receive their square roots to match BGLR's inverse-squared convention.
#' @param fixed_term_model_bayesian Character vector specifying the fixed terms in the Bayesian model.
#' @param rand_term_model_bayesian Character vector specifying the random terms in the Bayesian model.
#' @param pheno_data Data frame containing the phenotypic data.
#' @param gkernel Kernel matrix for genomic data (optional).
#' @param gmatrix Genomic relationship matrix (optional).
#' @param omic1_kernel Kernel matrix for the first omics data type (optional).
#' @param omic2_kernel Kernel matrix for the second omics data type (optional).
#' @param omic3_kernel Kernel matrix for the third omics data type (optional).
#' @param kernel_list Optional named list of additional relationship or kernel matrices.
#' @param gen_name Character string specifying the name of the genotype factor in the genomic data.
#' @param nIter Integer specifying the number of iterations for the Bayesian model.
#' @param burnIn Integer specifying the number of burn-in iterations for the Bayesian model.
#' @param thin Integer specifying the thinning interval for the Bayesian model.
#' @param heter_groups Character string identifying the grouping variable for heterogeneity (optional).
#' @param heter_resid Logical indicating whether Gaussian Bayesian fits should estimate residual variance by `heter_groups` when supported.
#' @param omics_kernel_label Character vector specifying labels for the omic data sets (optional).
#' @param cross_validation Logical indicating whether to return cross-validation preparation objects instead of fitting the final model.
#' @param response_family Response family for the Bayesian fit.
#' @param ... Additional arguments for future extensions.
#'
#' @return A list containing two elements: `bayes_result`, which holds the processed output from the Bayesian model, and `bayes_model`, which contains the raw model object from the Bayesian analysis.
#'
#' @details
#' The function integrates steps required for running specific Bayesian genomic
#' selection models, compiling necessary data and model specifications to
#' execute the model and format the output. Designed to work with RKHS, GBLUP_BRR,
#' and BRR approaches, it can handle complex scenarios involving multiple
#' sources of omic data and varying model parameters.
#'
#' For multi-environment Gaussian `RKHS` fits with random structure such as
#' `~ GID + GID:Env`, PredictProR now defaults the internal Bayesian ETA to:
#' \itemize{
#' \item main `GID` term as `BRR`
#' \item interaction `GID:Env` term as `RKHS`
#' }
#'
#' This avoids fitting two highly redundant RKHS kernels on the same GRM
#' backbone for both main and interaction effects, which previously depressed
#' the fitted marginal genetic variance and reliability. Single-environment
#' `RKHS` fits are unchanged.
#'
#' In multi-environment RKHS, `heter_resid = TRUE` selects the
#' environment-as-trait `BGLR::Multitrait` fit with environment-specific
#' genetic variances and diagonal environment-specific residual variances. This
#' mode requires `weights = NULL`. With `heter_resid = FALSE`, RKHS uses the
#' environment-aware genomic-main plus GxE kernel fit with common genetic and
#' residual variance parameters; that route accepts Stage 2 precision weights.
#' The effective mode and covariance structures are recorded in
#' `model_parameters`.
#'
#' \strong{GBLUP_BRR vs RKHS equivalence in MET heter_resid:} when
#' `heter_resid = TRUE` with multiple environments, both `RKHS` and
#' `GBLUP_BRR` are routed through `bayes_multitrait_env_heter_fit`, which
#' calls `BGLR::Multitrait` with a kernel-prior genetic term (`model = "RKHS"`
#' internally). With `X = eigen-sqrt(K)` the BRR-induced genetic covariance
#' \eqn{\sigma^2 X X^\top} equals the RKHS prior \eqn{\sigma^2 K}; BGLR's
#' multitrait sampler only implements the kernel-prior form, so MET
#' `GBLUP_BRR` and MET `RKHS` are \emph{numerically identical given the same
#' RNG seed}. This is by design; the `GS_model` label is preserved only in
#' the save-prefix and output metadata for traceability. In single-environment
#' fits the two diverge (separate `BRR`-design vs `RKHS`-kernel branches in
#' `BGLR::BGLR` with different prior scales), both anchored to the ASReml
#' GBLUP reference.
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
                                         kernel_list = NULL,
                                         gen_name = NULL,
                                         nIter = NULL,
                                         burnIn = NULL,
                                         thin = NULL,
                                         heter_groups  = NULL,
                                         heter_resid = FALSE,
                                         omics_kernel_label = NULL,
                                         cross_validation = FALSE,
                                         response_family = "gaussian",
                                         ...){

  dots <- list(...)
  if (identical(GS_model, "GBLUP_BRR")) {
    GS_model <- "BRR"
  }
  if (length(GS_model) != 1L || !(GS_model %in% c("RKHS", "BRR"))) {
    stop("Bayesian kernel finalization supports GS_model values 'RKHS', 'BRR', or alias 'GBLUP_BRR'. Use 'GBLUP' through the ASReml/GBLUP route.", call. = FALSE)
  }
  heter_control <- gp_normalize_single_environment_heter_controls(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    heter_resid = heter_resid,
    var_cov_str = NULL,
    response = response,
    response_family = response_family
  )
  heter_groups <- heter_control$heter_groups
  heter_resid <- heter_control$heter_resid
  fam <- gp_resolve_response_family(response_family, y = pheno_data[[response]])
  stage2_precision <- gp_resolve_stage2_precision_weights(
    weights = weights,
    pheno_data = pheno_data,
    response = response,
    gen_name = gen_name,
    context = "Bayesian kernel Stage 2 observation weights"
  )

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
                                    kernel_list = kernel_list,
                                    heter_groups = heter_groups,
                                    gen_name = gen_name)


bayes_para <-  bayes_parameter_check(nIter = nIter,
                                     burnIn = burnIn,
                                     thin = thin)

if(isTRUE(cross_validation)){
  return(list(bayes_para = bayes_para,
              bayes_ETA = ETA))
}

if (isTRUE(heter_resid) &&
    identical(fam, "gaussian") &&
    GS_model %in% c("RKHS", "BRR") &&
    !is.null(heter_groups) &&
    heter_groups %in% names(ETA[["pheno_data"]]) &&
    gp_bayes_has_multiple_residual_groups(ETA[["pheno_data"]], heter_groups) &&
    anyDuplicated(as.character(ETA[["pheno_data"]][[gen_name]])) > 0) {
  if (!is.null(weights)) {
    if (identical(GS_model, "RKHS")) {
      stop(
        paste(
          "RKHS cannot combine Stage 2 observation weights with heterogeneous-environment mode",
          "because BGLR::Multitrait has no observation-weight argument.",
          "Remove `weights` to estimate environment-specific genetic and residual variances,",
          "or set `bayes_kernel_heter_resid = FALSE` for weighted environment-aware RKHS",
          "with common genetic and residual variance parameters."
        ),
        call. = FALSE
      )
    }
    stop(
      "Weighted MET BGLR with heterogeneous residuals is unavailable because BGLR::Multitrait has no observation-weight argument. Set `bayes_kernel_heter_resid = FALSE` to use the weighted univariate BGLR route, or use weighted GP/ASReml.",
      call. = FALSE
    )
  }
  # Note: the prior `!gp_bayes_eta_supports_groups(ETA[["ETA"]])` guard was
  # removed here. That heuristic was meant for *marker* BRR (BayesA/B/C/BL/BRR)
  # where BGLR's native `groups` argument handles heterogeneous residuals well.
  # In this finalizer every BRR fit is kernel-based (GBLUP_BRR via eigen-sqrt
  # of K). Empirically the native `groups` path with two large BRR ETA terms
  # (main + GxE Hadamard) loses the genetic signal because BGLR's default
  # sigma^2_beta prior scale doesn't calibrate to the eigen-sqrt design the
  # same way RKHS's prior scales to the kernel. The Multitrait kernel-prior
  # route below produces the same N(0, sigma^2 * K) prior for both RKHS and
  # GBLUP_BRR, and matches the corgh per-env decomposition.
  kernels <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    gkernel = gkernel,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  return(bayes_multitrait_env_heter_fit(
    pheno_data = ETA[["pheno_data"]],
    response = response,
    gen_name = gen_name,
    heter_groups = heter_groups,
    kernels = kernels,
    GS_model = GS_model,
    bayes_para = bayes_para,
    omics_kernel_label = omics_kernel_label,
    confidence_level = dots[["confidence_level"]] %||% 0.95,
    CI_width_thresholds = dots[["CI_width_thresholds"]] %||% c(0.33, 0.66),
    high_reliability_thres = dots[["high_reliability_thres"]] %||% 0.9,
    low_reliability_thres = dots[["low_reliability_thres"]] %||% 0.5
  ))
}

bayes_groups <- gp_bayes_residual_groups(
  pheno_data = ETA[["pheno_data"]],
  heter_groups = heter_groups,
  heter_resid = heter_resid,
  ETA = ETA[["ETA"]],
  response_family = response_family
)

mod <-  bayes_mod_execute(pheno_data = ETA[["pheno_data"]],
                          response = response,
                          weights = weights,
                          ETA = ETA[["ETA"]],
                          bayes_para = bayes_para,
                          GS_model = GS_model,
                          response_family = response_family,
                          groups = bayes_groups

)


  # mod_output_bayes_RKHS
  res_model_output <- mod_output_bayes_gbluBRR_RKHS(mod = mod,
                                              ETA = ETA,
                                              gkernel =  gkernel,
                                              GS_model = GS_model,
                                              gmatrix = gmatrix,
                                              omic1_kernel = omic1_kernel,
                                              omic2_kernel = omic2_kernel,
                                              omic3_kernel = omic3_kernel,
                                              kernel_list = kernel_list,
                                               gen_name = gen_name,
                                               pheno_data = ETA[["pheno_data"]],
                                               response = response,
                                               heter_groups  = heter_groups,
                                              omics_kernel_label = omics_kernel_label,
                                              bayes_para = bayes_para,
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
                                              response_family = response_family
    )

  precision_values <- if (isTRUE(stage2_precision$supplied)) {
    as.double(stage2_precision$precision)
  } else {
    rep(1, nrow(pheno_data))
  }
  res_model_output[["model_parameters"]] <- data.frame(
    stat = c(
      "bayesian_kernel_model",
      "stage2_observation_weighting",
      "stage2_weight_supplied",
      "stage2_weight_source",
      "stage2_weight_count",
      "stage2_weight_nonunit_count",
      "stage2_precision_min",
      "stage2_precision_max",
      "bglr_native_weight_transform"
    ),
    summary = c(
      if (identical(GS_model, "BRR")) "GBLUP_BRR" else "RKHS",
      if (isTRUE(stage2_precision$supplied)) {
        "precision_via_BGLR_inverse_squared_weight"
      } else {
        "none_equal_precision"
      },
      as.character(isTRUE(stage2_precision$supplied)),
      as.character(stage2_precision$source %||% "none"),
      as.character(length(precision_values)),
      as.character(sum(abs(precision_values - 1) > sqrt(.Machine$double.eps))),
      format(min(precision_values), digits = 17, trim = TRUE),
      format(max(precision_values), digits = 17, trim = TRUE),
      if (isTRUE(stage2_precision$supplied)) "sqrt(user_precision)" else "none"
    ),
    stringsAsFactors = FALSE
  )
  if (!identical(fam, "gaussian")) {
    res_model_output[["model_parameters"]] <- rbind(
      res_model_output[["model_parameters"]],
      data.frame(
        stat = c(
          "classification_variance_scale",
          "classification_link_function",
          "liability_residual_variance",
          "liability_residual_variance_estimated",
          "classification_reliability_estimand",
          "classification_variance_posterior_source"
        ),
        summary = c(
          "latent_liability",
          "probit",
          "1",
          "FALSE",
          "maximum_predicted_class_probability",
          "BGLR_varU_or_varB_posterior_draws_with_kernel_diagonal_scaling"
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  bayes_kernel_bank <- gp_collect_kernel_inputs(
    gmatrix = gmatrix,
    gkernel = gkernel,
    omic1_kernel = omic1_kernel,
    omic2_kernel = omic2_kernel,
    omic3_kernel = omic3_kernel,
    kernel_list = kernel_list
  )
  res_model_output[["model_parameters"]] <- rbind(
    res_model_output[["model_parameters"]],
    gp_multi_kernel_parameter_rows(
      kernel_names = names(bayes_kernel_bank),
      strategy = "separate_BGLR_kernel_ETA_components"
    )
  )
  if (identical(GS_model, "RKHS")) {
    res_model_output[["model_parameters"]] <- rbind(
      res_model_output[["model_parameters"]],
      data.frame(
        stat = c(
          "rkhs_environment_mode",
          "rkhs_heterogeneous_environment_variances",
          "rkhs_genetic_environment_structure",
          "rkhs_residual_environment_structure",
          "bayes_kernel_heter_resid_effective",
          "rkhs_stage2_weight_compatibility"
        ),
        summary = c(
          "environment_aware_common_variance",
          "FALSE",
          "genomic_main_plus_gxe_compound_symmetry",
          if (isTRUE(stage2_precision$supplied)) {
            "single_residual_scale_with_observation_precision"
          } else {
            "single_residual_scale"
          },
          "FALSE",
          "weights_supported_in_common_variance_mode"
        ),
        stringsAsFactors = FALSE
      )
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
