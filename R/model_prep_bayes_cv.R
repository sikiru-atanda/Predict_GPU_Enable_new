

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
#' @param kernel_list Optional named list of additional relationship or kernel matrices.
#' @param heter_groups
#' @param heter_resid Logical indicating whether Bayesian Gaussian fits should estimate residual variance by `heter_groups` when supported.
#' @param bayes_kernel_heter_resid Optional logical override for heterogeneous
#' residual routing in kernel-Bayesian models. `NULL` inherits `heter_resid`.
#' @param omics_kernel_label
#' @param cross_validation
#' @param response_family Response family used by the Bayesian model workflow.
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
                               kernel_list = NULL,
                               heter_groups = NULL,
                               heter_resid = FALSE,
                               bayes_kernel_heter_resid = NULL,
                               omics_kernel_label = NULL,
                               cross_validation = TRUE,
                               response_family = "gaussian"
                                   ){
  # Phase 3.21: kernel-Bayes (RKHS / GBLUP_BRR) routing is independent of
  # marker-Bayes residual-groups handling. When the user passes
  # `bayes_kernel_heter_resid` explicitly it overrides `heter_resid` for the
  # kernel-Bayes routing decision only -- marker-Bayes (BayesA/B/C/BL/BRR)
  # still follows `heter_resid` as before. NULL means "inherit".
  kernel_heter_resid <- if (!is.null(bayes_kernel_heter_resid)) {
    isTRUE(bayes_kernel_heter_resid)
  } else {
    heter_resid
  }
#browser()
  bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")
  res_model_output_list <- list()
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
  if (is.null(bayes_kernel_heter_resid)) {
    kernel_heter_resid <- heter_resid
  }
  if (is.null(test_set)) {
    test_ids <- NULL
  } else {
    test_ids <- test_set
    if (is.data.frame(test_ids) || is.matrix(test_ids)) {
      test_ids <- test_ids[, 1]
    }
    test_ids <- as.character(test_ids)
  }
  drop_test_rows <- function(x, label) {
    if (is.null(x) || is.null(test_ids)) {
      return(x)
    }
    rn <- rownames(x)
    if (is.null(rn)) {
      stop(paste(label, "must have row names before Bayesian CV test-set filtering."), call. = FALSE)
    }
    x[!rn %in% test_ids, , drop = FALSE]
  }
  drop_test_kernel <- function(x, label) {
    if (is.null(x) || is.null(test_ids)) {
      return(x)
    }
    rn <- rownames(x)
    cn <- colnames(x)
    if (is.null(rn) || is.null(cn)) {
      stop(paste(label, "must have row and column names before Bayesian CV test-set filtering."), call. = FALSE)
    }
    x[!rn %in% test_ids, !cn %in% test_ids, drop = FALSE]
  }
  # if(!is.null(test_set)){
  # if(is.data.frame(test_set) | is.matrix(test_set)){
  #   test_set <-  test_set[, 1]
  #  }
  #
  # }

  for (model in GS_model_cv) {
  model_key <- model

  if (model %in% bayes_valid_models) {
  if (gp_is_multi_environment_trait_panel(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response
  )) {
    stop(
      paste(
        "Bayesian marker-regression models",
        paste(bayes_valid_models, collapse = ", "),
        "are currently single-environment only in this package.",
        "Use RKHS or GBLUP_BRR for multi-environment Bayesian cross-validation."
      ),
      call. = FALSE
    )
  }
  res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
                                                  random = random,
                                                  GS_model = model,
                                                  response = response,
                                                  weights = weights,
                                                  fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                                  #pheno_data = if(!is.null(test_set)) pheno_data[!pheno_data[[gen_name]] %in% test_set, ] else pheno_data,
                                                  pheno_data = pheno_data,
                                                  geno_data = drop_test_rows(geno_data, "geno_data"),
                                                  omic1_data = drop_test_rows(omic1_data, "omic1_data"),
                                                  omic2_data = drop_test_rows(omic2_data, "omic2_data"),
                                                  omic3_data = drop_test_rows(omic3_data, "omic3_data"),
                                                  gen_name = gen_name,
                                                  nIter = nIter,
                                                  burnIn = burnIn,
                                                  thin = thin,
                                                  omics_data_label = omics_data_label,
                                                  cross_validation = cross_validation,
                                                  response_family = response_family)
  res_model_output[["bayes_groups"]] <- gp_bayes_residual_groups(
    pheno_data = pheno_data,
    heter_groups = heter_groups,
    heter_resid = heter_resid,
    ETA = res_model_output[["bayes_ETA"]][["ETA"]],
    response_family = response_family
  )
  res_model_output_list[[model_key]] <-  res_model_output
  }

  if (model %in% bayes_gblup_valid_models) {
    model_fit <- if (identical(model, "GBLUP_BRR")) "BRR" else model
  res_model_output <- bayes_finalize_RKHS_GBLUPBRR(fixed = fixed,
                                                  random = random,
                                                  GS_model = model_fit,
                                                  response = response,
                                                  weights = weights,
                                                  fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                  rand_term_model_bayesian = rand_term_model_bayesian,
                                                  #pheno_data = if(!is.null(test_set)) pheno_data[!pheno_data[[gen_name]] %in% test_set, ] else pheno_data,
                                                  pheno_data = pheno_data,
                                                  gmatrix = drop_test_kernel(gmatrix, "gmatrix"),
                                                  omic1_kernel = drop_test_kernel(omic1_kernel, "omic1_kernel"),
                                                  omic2_kernel = drop_test_kernel(omic2_kernel, "omic2_kernel"),
                                                  omic3_kernel = drop_test_kernel(omic3_kernel, "omic3_kernel"),
                                                  kernel_list = gp_drop_test_kernel_bank(kernel_list, test_ids, "kernel_list"),
                                                  gen_name = gen_name,
                                                  heter_groups = heter_groups,
                                                  heter_resid = kernel_heter_resid,
                                                  nIter = nIter,
                                                  burnIn = burnIn,
                                                  thin = thin,
                                                  omics_kernel_label = omics_kernel_label,
                                                  cross_validation = cross_validation,
                                                  response_family = response_family)
  res_model_output[["bayes_groups"]] <- gp_bayes_residual_groups(
    pheno_data = pheno_data,
    heter_groups = heter_groups,
    heter_resid = kernel_heter_resid,
    ETA = res_model_output[["bayes_ETA"]][["ETA"]],
    response_family = response_family
  )

  # Attach the multitrait kernel-prior CV metadata when the fold-level fit
  # should route through bayes_mod_cv_met_kernel_predict() instead of the
  # univariate BGLR::BGLR(groups=...) path. This mirrors the true-prediction
  # routing fix (commit 35a3bcb) for the CV side.
  if (gp_bayes_should_route_met_kernel_cv(
        GS_model = model_fit,
        heter_resid = kernel_heter_resid,
        pheno_data = pheno_data,
        gen_name = gen_name,
        heter_groups = heter_groups,
        response_family = response_family)) {
    res_model_output[["met_kernel_cv_meta"]] <- list(
      pheno_data = pheno_data,
      response = response,
      gen_name = gen_name,
      heter_groups = heter_groups,
      kernels = gp_collect_kernel_inputs(
        gmatrix = gmatrix,
        gkernel = NULL,
        omic1_kernel = omic1_kernel,
        omic2_kernel = omic2_kernel,
        omic3_kernel = omic3_kernel,
        kernel_list = kernel_list
      ),
      GS_model = model_fit,
      bayes_para = res_model_output[["bayes_para"]],
      omics_kernel_label = omics_kernel_label
    )
  }


  res_model_output_list[[model_key]] <-  res_model_output
  }

  }

 return(res_model_output_list)

}
