#' Describe the general prediction input standard
#'
#' @return A list describing phenotype, model-family, genomic, and
#'   multi-environment input contracts.
#' @export
general_prediction_data_standard <- function() {
  list(
    phenotype = list(
      required = c("gen_name", "one or more response columns"),
      notes = c(
        "Use one row per genotype for single-environment workflows.",
        "If genotype rows repeat across environments, provide heter_groups.",
        "pheno_data is the preferred top-level phenotype input."
      )
    ),
    model_families = list(
      single_environment_ml_bayes = c("AI_valid_models", "bayes_valid_models"),
      single_response_single_environment_kernel_or_relationship_models =
        gp_multi_environment_kernel_display_names(
          gp_valid_models =
            gp_single_response_single_environment_gp_supported_models()
        ),
      structured_multi_response_gp_models = "FA-GBLUP"
    ),
    genomic_contracts = list(
      ml_bayes = c(
        "geno_data and/or omics data keyed by gen_name",
        "ML/DL may also use named PSD kernels through concatenated eigenfeature blocks"
      ),
      kernel_models = c(
        "gmatrix/gkernel/omic kernels or kernel_list keyed by gen_name",
        "or geno_data/omics data plus relationship/kernel construction instructions"
      )
    ),
    multi_environment_notes = c(
      "Repeated genotype rows across environments require heter_groups.",
      "Single-environment CV methods do not apply to repeated-row multi-environment tables.",
      paste(
        "FA-GBLUP requires a joint fit with at least two traits or at least",
        "two distinct environments; it is rejected for a single response in",
        "one environment."
      )
    )
  )
}

gp_general_contract_lines <- function() {
  spec <- general_prediction_data_standard()
  c(
    "PredictProR detected that the uploaded data is not in the required general prediction format.",
    "",
    "Phenotype standard:",
    "  required: gen_name plus one or more response columns",
    "  use one row per genotype for single-environment workflows",
    "  if genotype rows repeat across environments, provide heter_groups",
    "",
    "Genomic standard for single-environment ML/DL and Bayesian marker models:",
    "  provide geno_data and/or omics data keyed by gen_name",
    "  ML/DL may instead consume named PSD kernels through eigenfeature blocks",
    "",
    paste(
      "Genomic standard for kernel or relationship models",
      paste0("(", paste(gp_multi_environment_kernel_display_names(
        gp_valid_models =
          gp_single_response_single_environment_gp_supported_models()
      ), collapse = ", "), "):")
    ),
    "  provide gmatrix/gkernel/omic kernels or kernel_list keyed by gen_name",
    "  or provide geno_data/omics data together with relationship/kernel construction instructions",
    paste(
      "  FA-GBLUP requires a joint fit with at least two traits or at least",
      "two heter_groups levels"
    ),
    "",
    "Multi-environment note:",
    paste(" ", paste(spec$multi_environment_notes, collapse = " ")),
    "",
    "Please reformat the uploaded data to match this standard before rerunning model_execute()."
  )
}

gp_general_stop_contract <- function(details = NULL) {
  lines <- gp_general_contract_lines()
  if (!is.null(details) && nzchar(details)) {
    lines <- c(details, "", lines)
  }
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

#' Validate the general prediction input standard
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
#' @param kernel_list Optional named list of user-supplied relationship/kernel
#'   matrices keyed by gen_name.
#' @param response Response column names.
#' @param gen_name Genotype or sample identifier column name.
#' @param heter_groups Optional environment/grouping column name.
#' @param GS_model Optional selected true-prediction model.
#' @param GS_model_cv Optional selected cross-validation model.
#' @param cross_validation Logical indicating whether cross-validation is active.
#' @param cross_validation_meth Optional cross-validation method name.
#' @param AI_valid_models Valid machine-learning model names.
#' @param bayes_valid_models Valid Bayesian marker model names.
#' @param bayes_gblup_valid_models Valid Bayesian relationship/kernel model names.
#' @param gp_valid_models Valid GP relationship/kernel model names.
#' @param asreml_model Optional ASReml model name.
#'
#' @return Invisibly returns a validation list when inputs satisfy the contract.
#' @export
validate_general_input_standard <- function(pheno_data = NULL,
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
                                            heter_groups = NULL,
                                            GS_model = NULL,
                                            GS_model_cv = NULL,
                                            cross_validation = FALSE,
                                            cross_validation_meth = NULL,
                                            AI_valid_models = NULL,
                                            bayes_valid_models = NULL,
                                            bayes_gblup_valid_models = NULL,
                                            gp_valid_models = NULL,
                                            asreml_model = "GBLUP") {
  if (!is.null(pheno_data_train) || !is.null(pheno_data_test)) {
    if (is.null(pheno_data)) {
      gp_general_stop_contract(
        details = "Provide pheno_data, or ensure the split phenotype inputs follow the same gen_name/response contract."
      )
    }
  }
  if (is.null(pheno_data)) {
    return(invisible(list(valid = TRUE, standard = general_prediction_data_standard())))
  }
  ph <- as.data.frame(pheno_data, stringsAsFactors = FALSE)
  if (is.null(gen_name) || !nzchar(gen_name) || !gen_name %in% names(ph)) {
    gp_general_stop_contract(
      details = paste("The phenotype table must include the gen_name column. Current gen_name:", if (is.null(gen_name)) "NULL" else as.character(gen_name))
    )
  }
  if (!length(response) || !all(response %in% names(ph))) {
    gp_general_stop_contract(
      details = paste("Provide response column(s) present in pheno_data. Current response:", paste(response, collapse = ", "))
    )
  }

  active_models <- unique(stats::na.omit(gp_canonicalize_supported_model_names(c(GS_model, GS_model_cv))))
  if (!length(active_models)) {
    return(invisible(list(valid = TRUE, standard = general_prediction_data_standard())))
  }
  AI_valid_models <- gp_canonicalize_supported_model_names(AI_valid_models)
  gp_valid_models <- gp_valid_models %||% gp_lowrank_supported_models()
  gp_valid_models <- gp_canonicalize_supported_model_names(gp_valid_models)

  is_multi_env <- gp_is_multi_environment_trait_panel(
    pheno_data = ph,
    gen_name = gen_name,
    heter_groups = heter_groups,
    response = response
  )
  if (is_multi_env && (is.null(heter_groups) || !nzchar(heter_groups) || !heter_groups %in% names(ph))) {
    gp_general_stop_contract(
      details = "Repeated genotype rows were detected. Provide heter_groups as the environment/location column for multi-environment analysis."
    )
  }

  feature_inputs_present <- !all(vapply(list(geno_data, omic1_data, omic2_data, omic3_data), is.null, logical(1)))
  kernel_inputs_present <- !all(vapply(list(gmatrix, gkernel, omic1_kernel, omic2_kernel, omic3_kernel, kernel_list), is.null, logical(1)))

  ai_models <- active_models[active_models %in% AI_valid_models]
  bayes_marker_models <- active_models[active_models %in% bayes_valid_models]
  kernel_models <- active_models[active_models %in% gp_multi_environment_kernel_models(
    bayes_gblup_valid_models = bayes_gblup_valid_models,
    asreml_model = asreml_model,
    gp_valid_models = gp_valid_models
  )]

  if (length(ai_models)) {
    if (!feature_inputs_present && !kernel_inputs_present) {
      gp_general_stop_contract(
        details = paste(
          "The selected single-environment ML/DL models require raw genomic/omics features or named PSD kernels.",
          "Affected models:",
          paste(ai_models, collapse = ", ")
        )
      )
    }
  }
  if (length(bayes_marker_models) && !feature_inputs_present) {
    gp_general_stop_contract(
      details = paste(
        "The selected Bayesian marker-effect models require geno_data and/or omics feature data.",
        "Use GBLUP_BRR or RKHS when the inputs are relationship kernels only. Affected models:",
        paste(bayes_marker_models, collapse = ", ")
      )
    )
  }

  if (length(kernel_models)) {
    if (!feature_inputs_present && !kernel_inputs_present) {
      gp_general_stop_contract(
        details = paste(
          "The selected kernel/relationship models require either relationship/kernel inputs or feature data for relationship construction.",
          "Affected models:",
          paste(kernel_models, collapse = ", ")
        )
      )
    }
  }

  if (is_multi_env && isTRUE(cross_validation)) {
    cv_ok <- c("CV0", "CV1", "CV2", "Repeated_CV0", "Repeated_CV1", "Repeated_CV2")
    cv_token <- normalize_cv_token(cross_validation_meth %||% "")
    if (is.null(cross_validation_meth) || !cv_token %in% normalize_cv_token(cv_ok)) {
      gp_general_stop_contract(
        details = paste(
          "Repeated-row multi-environment phenotypes require a multi-environment CV method.",
          "Choose from:",
          paste(cv_ok, collapse = ", ")
        )
      )
    }
  }

  invisible(list(valid = TRUE, standard = general_prediction_data_standard()))
}
