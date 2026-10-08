#' Describe the multi-trait prediction input standard
#'
#' @return A list describing multi-trait phenotype, genomic, and current model
#'   scope requirements.
#' @export
multi_trait_data_standard <- function() {
  list(
    phenotype = list(
      required = c("gen_name", "at least two response columns"),
      optional = c("heter_groups for multi-environment layouts"),
      notes = c(
        "Use one row per genotype for single-environment analysis, or repeated genotype rows with heter_groups for MT-MET analysis.",
        "At least two response columns are required.",
        "Unbalanced trait panels are allowed; joint GP/Bayesian CV retains at least two observed training values per trait in every fold.",
        "Trait overlap across rows is still required for useful transfer."
      )
    ),
    genomics = list(
      multi_trait_ml = c(
        "geno_data and/or omics data keyed by gen_name",
        "gmatrix, omics kernels, or kernel_list converted to named eigenfeature blocks"
      ),
      multi_trait_dl = c(
        "geno_data and/or omics data keyed by gen_name",
        "gmatrix, omics kernels, or kernel_list converted to named eigenfeature blocks"
      ),
      multi_trait_gp = c(
        "gmatrix or gkernel keyed by gen_name",
        "geno_data and/or omics data keyed by gen_name for internal kernel construction",
        "optional omics kernels or kernel_list keyed by gen_name"
      ),
      multi_trait_bayes = c(
        "gmatrix or gkernel keyed by gen_name",
        "geno_data keyed by gen_name for internal relationship-matrix construction",
        "optional omics kernels or kernel_list keyed by gen_name"
      ),
      multi_trait_asreml = c(
        "one or more gmatrix, gkernel, omics kernels, or kernel_list entries keyed by gen_name",
        "geno_data and/or omics data keyed by gen_name for internal kernel construction"
      ),
      notes = c(
        "Row names must match genotype IDs in the phenotype table.",
        "multi_trait_asreml fits one independent trait-covariance random term per kernel and reports their summed genetic covariance.",
        "ML/DL models concatenate raw multi-omics blocks and user-supplied kernel eigenfeature blocks; their uncertainty remains predictive rather than a genetic variance component.",
        paste0(
          "Joint multi-trait GP models are purpose-specific: ",
          "Gaussian-Process-GBLUP uses unstructured REML, FA-GBLUP uses ",
          "factor-analytic covariance, and Scalable-GBLUP uses operator/MoM."
        ),
        "Kernel-GBLUP is single-trait in PredictProR and is rejected in joint multi-trait mode.",
        "Multi-trait multi-environment true prediction and CV0/CV1/CV2 cross-validation support ASReml GBLUP, GP, and Bayesian routes.",
        "Single-environment joint multi-trait CV supports ASReml, GP, Bayesian, ML, and DL routes. GP/Bayesian folds hold out genotypes and mask all response traits together.",
        "multi_trait_asreml supports genotype-blocked K-Folds and Repeated_K-Folds. Fold fits provide out-of-fold predictions; a stable direct/final fit provides the biological variance components."
      )
    ),
    current_scope = list(
      response_family = "gaussian",
      multi_trait_gp_models = gp_display_supported_model_names(gp_multitrait_gp_supported_models()),
      multi_trait_bayes_models = gp_multitrait_bayes_supported_models(),
      multi_trait_ml_models = gp_multitrait_ml_supported_models(),
      multi_trait_dl_models = gp_multitrait_dl_supported_models(),
      multi_trait_asreml_model = "GBLUP + engine = asreml",
      multi_trait_asreml_cv = "single-environment K-Folds/Repeated_K-Folds; MT-MET CV0/CV1/CV2 and repeated variants",
      multi_trait_gp_cv = "single-environment K-Folds/Repeated_K-Folds; MT-MET CV0/CV1/CV2 and repeated variants",
      multi_trait_bayes_cv = "single-environment K-Folds/Repeated_K-Folds; MT-MET CV0/CV1/CV2 and repeated variants",
      multi_trait_multi_environment_models = c(
        gp_multitrait_asreml_supported_models(),
        gp_display_supported_model_names(gp_multitrait_gp_supported_models()),
        gp_multitrait_bayes_supported_models()
      )
    )
  )
}

gp_met_classification_supported_models <- function(response_family = "binary") {
  fam <- gp_resolve_response_family(response_family)
  predictive <- gp_met_supported_models()
  # GBLUP_BRR is an eigen-square-root BRR parameterization of the same kernel
  # prior used by RKHS. In categorical MET fits the two are covariance-
  # equivalent, so exposing both would create a duplicate model comparison.
  kernel_bayes <- "RKHS"
  switch(
    fam,
    binary = unique(c(predictive, kernel_bayes)),
    multiclass = predictive,
    ordinal = kernel_bayes,
    character()
  )
}

gp_met_supported_models_for_response_family <- function(response_family = "gaussian") {
  fam <- gp_resolve_response_family(response_family)
  if (identical(fam, "gaussian")) {
    return(gp_multi_environment_supported_models())
  }
  gp_met_classification_supported_models(fam)
}

gp_met_classification_supported_display_names <- function(response_family = "binary") {
  gp_display_supported_model_names(
    gp_met_classification_supported_models(response_family)
  )
}

#' Describe the single-trait multi-environment input standard
#'
#' @return A list describing single-trait MET phenotype, genomic, CV scenario,
#'   supported-model, and model-purpose requirements.
#' @export
met_data_standard <- function() {
  structured_gp <- gp_display_supported_model_names(
    gp_single_trait_multi_environment_gp_supported_models()
  )
  kernel_bayes <- c("GBLUP_BRR", "RKHS")
  predictive_ml_dl <- gp_display_supported_model_names(gp_met_supported_models())
  current_models <- gp_multi_environment_supported_display_names()
  list(
    phenotype = list(
      required = c("gen_name", "response column", "heter_groups environment column"),
      notes = c(
        "Use repeated genotype rows across environments.",
        "The environment column must be passed through heter_groups.",
        "This contract is for one response measured in at least two environments, not a single-environment table."
      )
    ),
    genomics = list(
      accepted = c(
        "geno_data keyed by gen_name (a GRM is auto-built; gmatrix_method = 'Yang' by default, override to e.g. 'VanRaden')",
        "omics data keyed by gen_name",
        "gmatrix / gkernel / omic1_kernel / omic2_kernel / omic3_kernel (precomputed NxN PSD kernels keyed by gen_name)",
        "kernel_list (named list of additional NxN kernels keyed by gen_name; each is eigen-decomposed into a feature block alongside gmatrix and omic kernels)"
      ),
      notes = c(
        "Row names must match genotype IDs in the phenotype table.",
        "Single-trait MET analysis requires at least one genomic, omic, or relationship/kernel source.",
        "Internally each kernel is eigen-decomposed and the top eigenvectors (var_explained = 0.95 by default) are used as features; multi-kernel inputs become concatenated feature blocks before the ML/DL fit.",
        "When only geno_data is supplied (no kernel, no explicit gmatrix_method), the kernel-input guardrail auto-sets gmatrix_method = 'Yang' so the kernel-feature pipeline has a kernel to PCA over."
      )
    ),
    cv_scenarios = c("CV0", "CV1", "CV2", "Repeated_CV0", "Repeated_CV1", "Repeated_CV2"),
    current_models = current_models,
    current_scope = list(
      response_families = c("gaussian", "binary", "multiclass", "ordinal"),
      classification_models = list(
        binary = gp_met_classification_supported_display_names("binary"),
        multiclass = gp_met_classification_supported_display_names("multiclass"),
        ordinal = gp_met_classification_supported_display_names("ordinal")
      ),
      classification_cv_selection = list(
        binary = "pooled out-of-fold balanced accuracy; log loss tie-breaker",
        multiclass = "pooled out-of-fold macro F1; log loss tie-breaker",
        ordinal = paste(
          "pooled out-of-fold quadratic-weighted kappa;",
          "mean absolute class error then log loss tie-breakers"
        )
      ),
      kernel_bayesian_models = kernel_bayes,
      rkhs_environment_modes = c(
        heterogeneous_environment = paste(
          "Set bayes_kernel_heter_resid = TRUE with weights = NULL to fit",
          "BGLR::Multitrait unstructured genetic covariance and diagonal",
          "environment-specific residual variances."
        ),
        weighted_common_variance = paste(
          "Set bayes_kernel_heter_resid = FALSE when Stage 2 weights are supplied;",
          "this fits environment-aware genomic main and GxE kernels with common",
          "genetic and residual variance parameters."
        ),
        incompatible_combination = paste(
          "RKHS with both Stage 2 weights and heterogeneous-environment mode",
          "is rejected before fitting because BGLR::Multitrait cannot use weights."
        )
      ),
      gp_models = structured_gp,
      asreml_model = "GBLUP with engine = 'asreml' and an explicit MET variance-covariance structure",
      predictive_ml_dl_models = predictive_ml_dl,
      genetic_model_note = paste(
        "Kernel-GBLUP and Gaussian-Process-GBLUP use homogeneous/compound-symmetry GxE structure;",
        "FA-GBLUP uses factor-analytic GxE; GBLUP_BRR fits its existing kernel route;",
        "RKHS is either unweighted heterogeneous-environment or weighted environment-aware common-variance;",
        "ASReml GBLUP uses the user-selected MET covariance structure."
      ),
      predictive_model_note = paste(
        "The listed ML/DL models use explicit environment-aware MET features for prediction.",
        "Their uncertainty is predictive and is not a genetic variance, mixed-model PEV, or heritability."
      ),
      classification_note = paste(
        "Binary and multiclass MET ML/DL routes combine kernel eigenfeatures with explicit environment indicators.",
        "Binary and ordinal RKHS routes require the environment as a fixed effect and both genotype main",
        "and genotype-by-environment random terms. GBLUP_BRR is not duplicated in categorical MET because its",
        "kernel eigen-BRR covariance matches RKHS. Multiclass ML/DL models are not relabelled as ordinal models."
      ),
      excluded_models = c(
        "Scalable-GBLUP is excluded from single-trait MET because it resolves to the same KRR/krr_exact estimator as Kernel-GBLUP.",
        "GBLUP_BRR is excluded from categorical MET comparisons because its kernel eigen-BRR parameterization is covariance-equivalent to RKHS.",
        "Marker-regression Bayes models and all unlisted ML/DL models are single-environment only and are rejected before fitting."
      )
    )
  )
}

gp_multitrait_contract_lines <- function(mode = NULL) {
  spec <- multi_trait_data_standard()
  mode_label <- if (is.null(mode) || !nzchar(mode)) "multi-trait prediction" else paste0(mode, " multi-trait prediction")
  c(
    paste("PredictProR detected that the uploaded data is not in the required", mode_label, "format."),
    "",
    "Multi-trait phenotype standard:",
    "  required: gen_name plus at least two response columns",
    "  one row per genotype for single environment, or repeated rows plus heter_groups for MT-MET",
    "  unbalanced trait panels are allowed for multi_trait_ml and multi_trait_dl",
    "",
    "Current genomic standards:",
    "  multi_trait_ml and multi_trait_dl: geno_data/omics data or relationship kernels keyed by gen_name",
    "  multi_trait_gp: gmatrix/gkernel, geno_data/omics data, or user kernels keyed by gen_name",
    "  multi_trait_bayes: gmatrix/gkernel, user kernels, or geno_data keyed by gen_name",
    "  multi_trait_asreml: one or more relationship kernels, or raw genomic/omics data keyed by gen_name",
    "",
    "Current scope:",
    paste("  response_family:", spec$current_scope$response_family),
    paste("  multi_trait_gp models:", paste(spec$current_scope$multi_trait_gp_models, collapse = ", ")),
    paste("  multi_trait_bayes models:", paste(spec$current_scope$multi_trait_bayes_models, collapse = ", ")),
    paste("  multi_trait_ml models:", paste(spec$current_scope$multi_trait_ml_models, collapse = ", ")),
    paste("  multi_trait_dl models:", paste(spec$current_scope$multi_trait_dl_models, collapse = ", ")),
    paste("  multi_trait_asreml:", spec$current_scope$multi_trait_asreml_model),
    paste("  multi-trait multi-environment models:", paste(spec$current_scope$multi_trait_multi_environment_models, collapse = ", ")),
    "",
    "Please reformat the uploaded data to match this standard before rerunning model_execute()."
  )
}

gp_met_contract_lines <- function() {
  spec <- met_data_standard()
  c(
    "PredictProR detected that the request is not in the supported single-trait multi-environment format.",
    "",
    "MET phenotype standard:",
    "  required: gen_name, one response column, and an environment column",
    "  repeated genotype rows across environments are required",
    "  pass the environment column through heter_groups",
    "",
    "MET genomic standard:",
    "  provide geno_data and/or omics data keyed by gen_name",
    "  or provide gmatrix/omic kernels or kernel_list keyed by gen_name",
    "  row names must match genotype IDs in the phenotype table",
    "",
    "Supported MET CV scenarios:",
    paste(" ", paste(spec$cv_scenarios, collapse = ", ")),
    "",
    "Current single-trait MET models:",
    paste(" ", paste(spec$current_models, collapse = ", ")),
    "",
    "Model-purpose boundary:",
    paste(" ", spec$current_scope$genetic_model_note),
    paste(" ", spec$current_scope$predictive_model_note),
    paste(" ", paste(spec$current_scope$excluded_models, collapse = " ")),
    "",
    "Please reformat the uploaded data to match this standard before rerunning model_execute()."
  )
}

gp_multitrait_stop_contract <- function(details = NULL, mode = NULL) {
  lines <- gp_multitrait_contract_lines(mode = mode)
  if (!is.null(details) && nzchar(details)) {
    lines <- c(details, "", lines)
  }
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

gp_met_stop_contract <- function(details = NULL) {
  lines <- gp_met_contract_lines()
  if (!is.null(details) && nzchar(details)) {
    lines <- c(details, "", lines)
  }
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

gp_response_observed_masks <- function(pheno_data, response = NULL) {
  n <- if (is.null(pheno_data)) 0L else nrow(pheno_data)
  if (is.null(pheno_data) || !n) {
    return(list(.all = logical(0)))
  }

  response <- as.character(response)
  response <- response[!is.na(response) & nzchar(response) & response %in% names(pheno_data)]
  if (!length(response)) {
    return(list(.all = rep(TRUE, n)))
  }

  stats::setNames(lapply(response, function(trait) {
    y <- pheno_data[[trait]]
    if (is.numeric(y)) {
      is.finite(y)
    } else {
      !is.na(y)
    }
  }), response)
}

gp_is_multi_environment_trait_panel <- function(pheno_data,
                                                gen_name,
                                                heter_groups = NULL,
                                                response = NULL,
                                                response_family = NULL) {
  if (is.null(pheno_data) || is.null(gen_name) || !gen_name %in% names(pheno_data)) {
    return(FALSE)
  }

  masks <- gp_response_observed_masks(pheno_data, response = response)
  any(vapply(masks, function(mask) {
    mask[is.na(mask)] <- FALSE
    if (!any(mask)) {
      return(FALSE)
    }
    observed <- pheno_data[mask, , drop = FALSE]
    if (!is.null(heter_groups) && nzchar(heter_groups) && heter_groups %in% names(observed)) {
      env_count <- length(unique(stats::na.omit(as.character(observed[[heter_groups]]))))
      if (env_count > 1L) {
        return(TRUE)
      }
    }
    gids <- as.character(observed[[gen_name]])
    gids <- gids[!is.na(gids) & nzchar(gids)]
    length(gids) > length(unique(gids))
  }, logical(1)))
}

gp_normalize_heter_groups_for_trait_panel <- function(pheno_data,
                                                      gen_name,
                                                      heter_groups = NULL,
                                                      heter_resid = FALSE,
                                                      var_cov_str = NULL,
                                                      response = NULL,
                                                      response_family = NULL,
                                                      preserve_heter_groups = FALSE) {
  is_multi_env <- gp_is_multi_environment_trait_panel(
    pheno_data = pheno_data,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response,
    response_family = response_family
  )
  if (!isTRUE(is_multi_env) && !isTRUE(preserve_heter_groups)) {
    return(list(
      heter_groups = NULL,
      heter_resid = FALSE,
      var_cov_str = NULL,
      multi_environment = FALSE
    ))
  }
  list(
    heter_groups = heter_groups,
    heter_resid = isTRUE(heter_resid),
    var_cov_str = var_cov_str,
    multi_environment = isTRUE(is_multi_env)
  )
}

gp_selected_model_for_stage <- function(GS_model = NULL, GS_model_cv = NULL, cross_validation = FALSE) {
  model <- if (isTRUE(cross_validation)) GS_model_cv else GS_model
  model <- unique(stats::na.omit(as.character(model %||% character())))
  model[nzchar(model)]
}

gp_multitrait_allowed_models_for_mode <- function(mode, is_met = FALSE, cross_validation = FALSE) {
  if (isTRUE(is_met)) {
    return(switch(
      mode,
      multi_trait_asreml = gp_multitrait_asreml_supported_models(),
      multi_trait_gp = gp_multitrait_gp_supported_models(),
      multi_trait_bayes = gp_multitrait_bayes_supported_models(),
      character()
    ))
  }
  switch(
    mode,
    multi_trait_gp = gp_multitrait_gp_supported_models(),
    multi_trait_bayes = gp_multitrait_bayes_supported_models(),
    multi_trait_asreml = gp_multitrait_asreml_supported_models(),
    multi_trait_ml = gp_multitrait_ml_supported_models(),
    multi_trait_dl = gp_multitrait_dl_supported_models(),
    character()
  )
}

gp_multitrait_model_scope_label <- function(mode, is_met = FALSE) {
  if (isTRUE(is_met)) {
    return("multi-trait multi-environment prediction")
  }
  switch(
    mode,
    multi_trait_gp = "Multi-trait GP",
    multi_trait_bayes = "Multi-trait Bayesian",
    multi_trait_asreml = "Multi-trait ASReml-R",
    multi_trait_ml = "Multi-trait ML",
    multi_trait_dl = "Multi-trait DL",
    "Multi-trait prediction"
  )
}

gp_validate_multitrait_model_selection <- function(mode,
                                                   selected_model,
                                                   is_met = FALSE,
                                                   cross_validation = FALSE) {
  selected_model <- gp_canonicalize_supported_model_names(selected_model)
  if (!length(selected_model)) {
    return(invisible(TRUE))
  }
  if (length(selected_model) != 1L) {
    gp_multitrait_stop_contract(
      details = paste(
        gp_multitrait_model_scope_label(mode, is_met = is_met),
        "currently expects one joint model at a time. Provided:",
        paste(selected_model, collapse = ", ")
      ),
      mode = mode
    )
  }
  allowed <- gp_multitrait_allowed_models_for_mode(
    mode = mode,
    is_met = is_met,
    cross_validation = cross_validation
  )
  if (!length(allowed) || !selected_model %in% allowed) {
    allowed_display <- gp_display_supported_model_names(allowed)
    mode_hint <- if (isTRUE(is_met)) {
      "MT-MET CV uses CV0, CV1, CV2, or their repeated variants."
    } else {
      ""
    }
    gp_multitrait_stop_contract(
      details = paste(
        gp_multitrait_model_scope_label(mode, is_met = is_met),
        "supports only:",
        paste(allowed_display, collapse = ", "),
        mode_hint
      ),
      mode = mode
    )
  }
  invisible(TRUE)
}

#' Validate the multi-trait prediction input standard
#'
#' @param pheno_data Optional phenotype data frame.
#' @param pheno_data_train Optional pre-split training phenotype data.
#' @param pheno_data_test Optional pre-split testing phenotype data.
#' @param geno_data Optional genomic marker matrix.
#' @param omic1_data Optional first omic feature matrix.
#' @param omic2_data Optional second omic feature matrix.
#' @param omic3_data Optional third omic feature matrix.
#' @param gmatrix Optional genomic relationship matrix.
#' @param gkernel Optional genomic kernel matrix.
#' @param omic1_kernel Optional first omic kernel matrix.
#' @param omic2_kernel Optional second omic kernel matrix.
#' @param omic3_kernel Optional third omic kernel matrix.
#' @param kernel_list Optional named list of user-supplied kernels.
#' @param response Response column names.
#' @param gen_name Genotype or sample identifier column name.
#' @param response_family Response family, usually \code{"auto"} or
#'   \code{"gaussian"}.
#' @param multi_trait_asreml Logical indicating whether ASReml multi-trait mode
#'   is active.
#' @param multi_trait_ml Logical indicating whether ML multi-trait mode is
#'   active.
#' @param multi_trait_dl Logical indicating whether DL multi-trait mode is
#'   active.
#' @param multi_trait_gp Logical indicating whether GP multi-trait mode is
#'   active.
#' @param multi_trait_bayes Logical indicating whether Bayesian BGLR
#'   multi-trait mode is active.
#' @param heter_groups Optional environment/grouping column name.
#' @param GS_model Optional selected true-prediction model.
#' @param GS_model_cv Optional selected cross-validation model.
#' @param cross_validation Logical indicating whether cross-validation is active.
#' @param cross_validation_meth Optional cross-validation method name. MT-MET
#'   workflows require CV0, CV1, CV2, or a repeated variant.
#'
#' @return Invisibly returns a validation list when inputs satisfy the contract.
#' @export
validate_multi_trait_input_standard <- function(pheno_data = NULL,
                                                pheno_data_train = NULL,
                                                pheno_data_test = NULL,
                                                geno_data = NULL,
                                                omic1_data = NULL,
                                                omic2_data = NULL,
                                                omic3_data = NULL,
                                                gmatrix = NULL,
                                                gkernel = NULL,
                                                omic1_kernel = NULL,
                                                omic2_kernel = NULL,
                                                omic3_kernel = NULL,
                                                kernel_list = NULL,
                                                response = NULL,
                                                gen_name = NULL,
                                                response_family = "auto",
                                                multi_trait_asreml = FALSE,
                                                multi_trait_ml = FALSE,
                                                multi_trait_dl = FALSE,
                                                multi_trait_gp = FALSE,
                                                multi_trait_bayes = FALSE,
                                                heter_groups = NULL,
                                                GS_model = NULL,
                                                GS_model_cv = NULL,
                                                cross_validation = FALSE,
                                                cross_validation_meth = NULL) {
  active <- c(
    multi_trait_asreml = isTRUE(multi_trait_asreml),
    multi_trait_ml = isTRUE(multi_trait_ml),
    multi_trait_dl = isTRUE(multi_trait_dl),
    multi_trait_gp = isTRUE(multi_trait_gp),
    multi_trait_bayes = isTRUE(multi_trait_bayes)
  )
  if (!any(active)) {
    return(invisible(list(valid = TRUE, mode = NULL, standard = multi_trait_data_standard())))
  }
  if (sum(active) > 1L) {
    gp_multitrait_stop_contract(
      details = paste("Select only one multi-trait workflow at a time. Current flags:", paste(names(active)[active], collapse = ", ")),
      mode = "general"
    )
  }
  active_mode <- names(active)[active]

  if (!is.null(pheno_data_train) || !is.null(pheno_data_test)) {
    gp_multitrait_stop_contract(
      details = "Multi-trait workflows currently expect the phenotype matrix in pheno_data, not split only across pheno_data_train/pheno_data_test.",
      mode = active_mode
    )
  }
  if (is.null(pheno_data)) {
    gp_multitrait_stop_contract(
      details = "Provide pheno_data for multi-trait prediction.",
      mode = active_mode
    )
  }
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)

  if (is.null(gen_name) || !nzchar(gen_name) || !gen_name %in% names(ph)) {
    gp_multitrait_stop_contract(
      details = paste("The multi-trait phenotype table must include the gen_name column. Current gen_name:", if (is.null(gen_name)) "NULL" else as.character(gen_name)),
      mode = active_mode
    )
  }
  if (length(response) < 2L || !all(response %in% names(ph))) {
    gp_multitrait_stop_contract(
      details = paste("Provide at least two response columns present in pheno_data. Current response:", paste(response, collapse = ", ")),
      mode = active_mode
    )
  }
  fam <- gp_normalize_response_family(response_family)
  if (!identical(fam, "gaussian") && !identical(fam, "auto")) {
    gp_multitrait_stop_contract(
      details = "Current multi-trait workflows support gaussian traits only.",
      mode = active_mode
    )
  }

  is_met <- gp_is_multi_environment_trait_panel(
    pheno_data = ph,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response,
    response_family = response_family
  )
  if (isTRUE(is_met) && isTRUE(cross_validation)) {
    gp_multitrait_auto_cv_method(
      cross_validation_meth = cross_validation_meth,
      is_met = TRUE
    )
  }
  selected_model <- gp_selected_model_for_stage(
    GS_model = GS_model,
    GS_model_cv = GS_model_cv,
    cross_validation = cross_validation
  )
  gp_validate_multitrait_model_selection(
    mode = active_mode,
    selected_model = selected_model,
    is_met = is_met,
    cross_validation = cross_validation
  )
  if (isTRUE(is_met) && !(active_mode %in% c("multi_trait_asreml", "multi_trait_gp", "multi_trait_bayes"))) {
    gp_multitrait_stop_contract(
      details = paste(
        "Current multi-trait multi-environment prediction is implemented through multi_trait_asreml, multi_trait_gp, or multi_trait_bayes.",
        "Use one of:",
        paste(c(
          gp_multitrait_asreml_supported_models(),
          gp_display_supported_model_names(gp_multitrait_gp_supported_models()),
          gp_multitrait_bayes_supported_models()
        ), collapse = ", ")
      ),
      mode = active_mode
    )
  }

  ids <- unique(as.character(ph[[gen_name]]))
  matrix_overlaps_ids <- function(x, require_square = TRUE) {
    if (is.null(x)) {
      return(FALSE)
    }
    rn <- rownames(x)
    cn <- colnames(x)
    has_rows <- !is.null(rn) && length(rn) && any(ids %in% rn)
    if (!isTRUE(require_square)) {
      return(isTRUE(has_rows))
    }
    has_cols <- !is.null(cn) && length(cn) && any(ids %in% cn)
    isTRUE(has_rows && has_cols)
  }
  feature_overlaps_ids <- function(x) {
    if (is.null(x)) {
      return(FALSE)
    }
    rn <- rownames(x)
    !is.null(rn) && length(rn) && any(ids %in% rn)
  }
  if (active_mode %in% c("multi_trait_ml", "multi_trait_dl")) {
    feature_sources <- list(geno_data, omic1_data, omic2_data, omic3_data)
    kernel_sources <- list(
      gmatrix = gmatrix, gkernel = gkernel,
      omic1_kernel = omic1_kernel, omic2_kernel = omic2_kernel,
      omic3_kernel = omic3_kernel
    )
    if (is.list(kernel_list) && length(kernel_list)) kernel_sources <- c(kernel_sources, kernel_list)
    has_feature_source <- any(vapply(feature_sources, feature_overlaps_ids, logical(1L)))
    has_kernel_source <- any(vapply(kernel_sources, matrix_overlaps_ids, logical(1L)))
    if (!isTRUE(has_feature_source) && !isTRUE(has_kernel_source)) {
      gp_multitrait_stop_contract(
        details = paste(
          active_mode,
          "requires raw genomic/omics features or one or more relationship kernels keyed by gen_name."
        ),
        mode = active_mode
      )
    }
    if (!is.null(geno_data)) {
      rn <- rownames(geno_data)
      if (is.null(rn) || !length(rn) || !sum(ids %in% rn)) {
        gp_multitrait_stop_contract(
          details = "geno_data must be keyed by gen_name and overlap the genotype IDs in pheno_data.",
          mode = active_mode
        )
      }
    }
    bad_kernels <- names(kernel_sources)[
      !vapply(kernel_sources, is.null, logical(1L)) &
        !vapply(kernel_sources, matrix_overlaps_ids, logical(1L))
    ]
    if (length(bad_kernels)) {
      gp_multitrait_stop_contract(
        details = paste(
          "All supplied multi-trait ML/DL kernels must have row and column names keyed by gen_name.",
          "Kernels without sufficient overlap:", paste(bad_kernels, collapse = ", ")
        ),
        mode = active_mode
      )
    }
  }

  if (identical(active_mode, "multi_trait_gp")) {
    kernel_sources <- list(
      gmatrix = gmatrix,
      gkernel = gkernel,
      omic1_kernel = omic1_kernel,
      omic2_kernel = omic2_kernel,
      omic3_kernel = omic3_kernel
    )
    if (is.list(kernel_list) && length(kernel_list)) {
      kernel_sources <- c(kernel_sources, kernel_list)
    }
    has_kernel_source <- any(vapply(kernel_sources, matrix_overlaps_ids, logical(1)))
    feature_sources <- list(
      geno_data = geno_data,
      omic1_data = omic1_data,
      omic2_data = omic2_data,
      omic3_data = omic3_data
    )
    has_feature_source <- any(vapply(feature_sources, feature_overlaps_ids, logical(1)))
    if (!isTRUE(has_kernel_source) && !isTRUE(has_feature_source)) {
      gp_multitrait_stop_contract(
        details = paste(
          "multi_trait_gp requires a relationship/kernel source keyed by gen_name",
          "(gmatrix, gkernel, omics kernel, or kernel_list), or geno_data/omics data",
          "so PredictProR can construct the kernel internally."
        ),
        mode = active_mode
      )
    }
    bad_kernels <- names(kernel_sources)[
      !vapply(kernel_sources, is.null, logical(1)) &
        !vapply(kernel_sources, matrix_overlaps_ids, logical(1))
    ]
    if (length(bad_kernels)) {
      gp_multitrait_stop_contract(
        details = paste(
          "All supplied multi_trait_gp kernels must have row and column names keyed by gen_name.",
          "Kernels without sufficient overlap:",
          paste(bad_kernels, collapse = ", ")
        ),
        mode = active_mode
      )
    }
    bad_features <- names(feature_sources)[
      !vapply(feature_sources, is.null, logical(1)) &
        !vapply(feature_sources, feature_overlaps_ids, logical(1))
    ]
    if (length(bad_features)) {
      gp_multitrait_stop_contract(
        details = paste(
          "All supplied multi_trait_gp feature matrices must have row names keyed by gen_name.",
          "Feature matrices without sufficient overlap:",
          paste(bad_features, collapse = ", ")
        ),
        mode = active_mode
      )
    }
    if (isTRUE(is_met) &&
        (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph))) {
      gp_multitrait_stop_contract(
        details = "Multi-trait multi-environment GP requires heter_groups as the environment column.",
        mode = active_mode
      )
    }
  }

  if (identical(active_mode, "multi_trait_bayes")) {
    kernel_sources <- list(
      gmatrix = gmatrix,
      gkernel = gkernel,
      omic1_kernel = omic1_kernel,
      omic2_kernel = omic2_kernel,
      omic3_kernel = omic3_kernel
    )
    if (is.list(kernel_list) && length(kernel_list)) {
      kernel_sources <- c(kernel_sources, kernel_list)
    }
    has_kernel_source <- any(vapply(kernel_sources, matrix_overlaps_ids, logical(1)))
    has_geno_source <- feature_overlaps_ids(geno_data)
    if (!isTRUE(has_kernel_source) && !isTRUE(has_geno_source)) {
      gp_multitrait_stop_contract(
        details = paste(
          "multi_trait_bayes requires a relationship/kernel source keyed by gen_name,",
          "or geno_data keyed by gen_name so PredictProR can construct the relationship matrix internally."
        ),
        mode = active_mode
      )
    }
    if (!is.null(geno_data) && !isTRUE(has_geno_source)) {
      gp_multitrait_stop_contract(
        details = "geno_data must have row names keyed by gen_name for multi_trait_bayes.",
        mode = active_mode
      )
    }
    bad_kernels <- names(kernel_sources)[
      !vapply(kernel_sources, is.null, logical(1)) &
        !vapply(kernel_sources, matrix_overlaps_ids, logical(1))
    ]
    if (length(bad_kernels)) {
      gp_multitrait_stop_contract(
        details = paste(
          "All supplied multi_trait_bayes kernels must have row and column names keyed by gen_name.",
          "Kernels without sufficient overlap:",
          paste(bad_kernels, collapse = ", ")
        ),
        mode = active_mode
      )
    }
    if (isTRUE(is_met) &&
        (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph))) {
      gp_multitrait_stop_contract(
        details = "Multi-trait multi-environment Bayesian prediction requires heter_groups as the environment column.",
        mode = active_mode
      )
    }
  }

  if (identical(active_mode, "multi_trait_asreml")) {
    feature_sources <- list(
      geno_data = geno_data, omic1_data = omic1_data,
      omic2_data = omic2_data, omic3_data = omic3_data
    )
    kernel_sources <- list(
      gmatrix = gmatrix, gkernel = gkernel,
      omic1_kernel = omic1_kernel, omic2_kernel = omic2_kernel,
      omic3_kernel = omic3_kernel
    )
    if (is.list(kernel_list) && length(kernel_list)) kernel_sources <- c(kernel_sources, kernel_list)
    has_feature_source <- any(vapply(feature_sources, feature_overlaps_ids, logical(1L)))
    has_kernel_source <- any(vapply(kernel_sources, matrix_overlaps_ids, logical(1L)))
    if (!isTRUE(has_kernel_source) && !isTRUE(has_feature_source)) {
      gp_multitrait_stop_contract(
        details = paste(
          "multi_trait_asreml requires one or more relationship kernels keyed by gen_name,",
          "or raw genomic/omics data from which model-ready kernels can be constructed."
        ),
        mode = active_mode
      )
    }
    bad_kernels <- names(kernel_sources)[
      !vapply(kernel_sources, is.null, logical(1L)) &
        !vapply(kernel_sources, matrix_overlaps_ids, logical(1L))
    ]
    if (length(bad_kernels)) {
      gp_multitrait_stop_contract(
        details = paste(
          "Every multi_trait_asreml kernel must be a named square matrix keyed by gen_name.",
          "Kernels without sufficient overlap:", paste(bad_kernels, collapse = ", ")
        ),
        mode = active_mode
      )
    }
    bad_features <- names(feature_sources)[
      !vapply(feature_sources, is.null, logical(1L)) &
        !vapply(feature_sources, feature_overlaps_ids, logical(1L))
    ]
    if (length(bad_features)) {
      gp_multitrait_stop_contract(
        details = paste(
          "Every raw multi_trait_asreml genomic/omics matrix must be keyed by gen_name.",
          "Feature matrices without sufficient overlap:", paste(bad_features, collapse = ", ")
        ),
        mode = active_mode
      )
    }
    if (isTRUE(is_met) &&
        (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph))) {
      gp_multitrait_stop_contract(
        details = "Grouped MT-MET ASReml requires heter_groups as the environment column.",
        mode = active_mode
      )
    }
  }

  invisible(list(valid = TRUE, mode = active_mode, standard = multi_trait_data_standard()))
}

#' Validate the MET ML/DL input standard
#'
#' @param pheno_data Optional phenotype data frame.
#' @param pheno_data_train Optional pre-split training phenotype data.
#' @param pheno_data_test Optional pre-split testing phenotype data.
#' @param geno_data Optional genomic marker matrix.
#' @param omic1_data Optional first omic feature matrix.
#' @param omic2_data Optional second omic feature matrix.
#' @param omic3_data Optional third omic feature matrix.
#' @param gmatrix Optional genomic relationship matrix.
#' @param omic1_kernel Optional first omic kernel matrix.
#' @param omic2_kernel Optional second omic kernel matrix.
#' @param omic3_kernel Optional third omic kernel matrix.
#' @param kernel_list Optional named list of user-supplied relationship/kernel
#'   matrices keyed by gen_name.
#' @param response Response column name.
#' @param response_family Response family for the MET request.
#' @param gen_name Genotype or sample identifier column name.
#' @param heter_groups Environment/grouping column name.
#' @param met_ml_dl Logical indicating whether MET ML/DL mode is active.
#' @param GS_model Optional selected true-prediction model.
#' @param GS_model_cv Optional selected cross-validation model.
#' @param cross_validation Logical indicating whether cross-validation is active.
#' @param cross_validation_meth Optional cross-validation method name.
#'
#' @return Invisibly returns a validation list when inputs satisfy the contract.
#' @export
validate_met_input_standard <- function(pheno_data = NULL,
                                        pheno_data_train = NULL,
                                        pheno_data_test = NULL,
                                        geno_data = NULL,
                                        omic1_data = NULL,
                                        omic2_data = NULL,
                                        omic3_data = NULL,
                                        gmatrix = NULL,
                                        omic1_kernel = NULL,
                                        omic2_kernel = NULL,
                                        omic3_kernel = NULL,
                                        kernel_list = NULL,
                                        response = NULL,
                                        response_family = "gaussian",
                                        gen_name = NULL,
                                        heter_groups = NULL,
                                        met_ml_dl = FALSE,
                                        GS_model = NULL,
                                        GS_model_cv = NULL,
                                        cross_validation = FALSE,
                                        cross_validation_meth = NULL) {
  if (!isTRUE(met_ml_dl)) {
    return(invisible(list(valid = TRUE, standard = met_data_standard())))
  }
  if (!is.null(pheno_data_train) || !is.null(pheno_data_test)) {
    gp_met_stop_contract(
      details = "MET ML/DL currently expects the phenotype table in pheno_data."
    )
  }
  if (is.null(pheno_data)) {
    gp_met_stop_contract(
      details = "Provide pheno_data for MET ML/DL."
    )
  }
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  if (is.null(gen_name) || !nzchar(gen_name) || !gen_name %in% names(ph)) {
    gp_met_stop_contract(
      details = paste("The MET phenotype table must include the gen_name column. Current gen_name:", if (is.null(gen_name)) "NULL" else as.character(gen_name))
    )
  }
  # MET ML/DL accepts one or more response columns; the CV dispatcher in
  # main_crossvalidation_execution_logic.R iterates per-trait per-fold and
  # passes a single string into gp_met_python_cv_predict / gp_met_dl_cv_predict,
  # so multi-trait MET reduces to a sequence of single-trait MET fits.
  response_chr <- as.character(response %||% character())
  if (!length(response_chr) || !all(nzchar(response_chr)) || !all(response_chr %in% names(ph))) {
    gp_met_stop_contract(
      details = paste("MET ML/DL response column(s) must be present in pheno_data. Current response:",
                      paste(response_chr, collapse = ", "))
    )
  }
  if (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph)) {
    gp_met_stop_contract(
      details = "Provide heter_groups as the environment column for MET ML/DL."
    )
  }
  if (!gp_is_multi_environment_trait_panel(
    pheno_data = ph,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response
  )) {
    gp_met_stop_contract(
      details = "MET ML/DL requires repeated genotype rows across environments."
    )
  }
  feature_sources <- list(geno_data, omic1_data, omic2_data, omic3_data)
  kernel_sources <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
  kernel_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
  if (is.list(kernel_list) && length(kernel_list)) {
    kernel_sources <- c(kernel_sources, kernel_list)
    kernel_names <- c(kernel_names, names(kernel_list) %||% paste0("kernel_list_", seq_along(kernel_list)))
  }
  # The MET ML/DL kernel-feature pipeline needs at least one feature or kernel
  # source. geno_data on its own IS valid here: gp_apply_model_input_guardrails
  # -> gp_auto_promote_met_ml_dl_kernel_build will auto-set gmatrix_method =
  # "Yang" downstream so process_geno_data() builds a GRM the pipeline can
  # PCA over. Override by passing gmatrix_method explicitly or a precomputed
  # gmatrix/kernel_list upfront.
  if (all(vapply(feature_sources, is.null, logical(1)))) {
    if (all(vapply(kernel_sources, is.null, logical(1)))) {
      gp_met_stop_contract(
        details = "MET ML/DL requires geno_data/omics data, gmatrix/gkernel/omic kernels, or kernel_list keyed by gen_name."
      )
    }
  }
  ids <- unique(as.character(ph[[gen_name]]))
  if (!is.null(geno_data)) {
    rn <- rownames(geno_data)
    if (is.null(rn) || !length(rn) || !sum(ids %in% rn)) {
      gp_met_stop_contract(
        details = "geno_data must be keyed by gen_name and overlap the genotype IDs in pheno_data."
      )
    }
  }
  for (i in seq_along(kernel_sources)) {
    K <- kernel_sources[[i]]
    if (is.null(K)) {
      next
    }
    rn <- rownames(K)
    cn <- colnames(K)
    if (is.null(rn) || is.null(cn) || !length(rn) || !length(cn) ||
        !sum(ids %in% rn) || !sum(ids %in% cn)) {
      gp_met_stop_contract(
        details = paste0(kernel_names[[i]], " must have row and column names keyed by gen_name and overlapping pheno_data.")
      )
    }
  }
  # Every listed model must have a real single-trait MET route. A mixed request
  # is rejected as a unit so unsupported models are never silently dropped
  # while an applicable model continues fitting.
  active_model <- unique(stats::na.omit(c(GS_model, GS_model_cv)))
  active_canon <- tryCatch(
    gp_canonicalize_supported_model_names(active_model),
    error = function(e) active_model
  )
  fam <- gp_resolve_response_family(
    response_family,
    y = if (length(response_chr) == 1L) ph[[response_chr]] else NULL
  )
  supported <- gp_met_supported_models_for_response_family(fam)
  unsupported <- setdiff(active_canon, supported)
  if (length(unsupported)) {
    gp_met_stop_contract(
      details = paste(
        paste0("Unsupported single-trait multi-environment response_family='", fam, "' model(s):"),
        paste(gp_display_supported_model_names(unsupported), collapse = ", "),
        "The request was stopped before fitting; no model was silently dropped. Supported models are:",
        paste(gp_display_supported_model_names(supported), collapse = ", ")
      )
    )
  }
  if (isTRUE(cross_validation)) {
    cv_ok <- c("CV0", "CV1", "CV2", "Repeated_CV0", "Repeated_CV1", "Repeated_CV2")
    cv_token <- normalize_cv_token(cross_validation_meth %||% "")
    if (is.null(cross_validation_meth) || !cv_token %in% normalize_cv_token(cv_ok)) {
      gp_met_stop_contract(
        details = paste("MET ML/DL cross-validation requires one of:", paste(cv_ok, collapse = ", "))
      )
    }
  }
  invisible(list(valid = TRUE, standard = met_data_standard()))
}
