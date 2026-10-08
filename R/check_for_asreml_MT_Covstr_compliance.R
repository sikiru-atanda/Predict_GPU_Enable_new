# Function to check for ASReml software requirement
checkForASReml <- function(engine, GS_model, GS_model_cv, cross_validation, msg) {
  msg <- ""
  if ((!isFALSE(cross_validation) && "GBLUP" %in% GS_model_cv) ||
      (isFALSE(cross_validation) && "GBLUP" %in% GS_model)) {
    if (engine != "asreml") {
      stop(paste(msg, "ASReml software is required to fit GBLUP for single or multi-environment."), call. = FALSE)
    }
  }
}

# Function to check variance-covariance structure inputs
checkVarianceCovarianceInputs <- function(heter_groups, heter_resid, var_cov_str, var_cov_str_available, msg) {
  msg <- ""
  if (!is.null(heter_groups) && !is.null(heter_resid) && is.null(var_cov_str)) {
    stop(msg, "Your data suggest multi-environment but variance-covariance structure is missing. Choose from: ", paste(var_cov_str_available, collapse = ", "), call. = FALSE)
  } else if (!is.null(heter_groups) && is.null(heter_resid) && !is.null(var_cov_str)) {
    stop(msg, "Your data suggest multi-environment but heter_resid must be TRUE.", call. = FALSE)
  } else if (!is.null(var_cov_str) && !(var_cov_str %in% var_cov_str_available)) {
    stop(msg, "Invalid output variance-covariance structure. Choose from: ", paste(var_cov_str_available, collapse = ", "), call. = FALSE)
  }
}

# Function to validate multi-environment structure and inputs
validateMultiEnvironment <- function(pheno_data, gen_name, heter_groups,
                                     heter_resid, var_cov_str, GS_model, var_cov_str_available,
                                     cross_validation, GS_model_cv, msg,
                                     response = NULL, response_family = NULL) {
  msg <- ""
  is_multi_environment <- length(pheno_data[[gen_name]]) > length(unique(pheno_data[[gen_name]]))
  if (exists("gp_is_multi_environment_trait_panel", mode = "function") && !is.null(response)) {
    is_multi_environment <- gp_is_multi_environment_trait_panel(
      pheno_data = pheno_data,
      gen_name = gen_name,
      heter_groups = heter_groups,
      response = response,
      response_family = response_family
    )
  }
  if (isTRUE(is_multi_environment)) {
    if (is.null(heter_groups)) {
      stop(paste(msg, "Your phenotypic data has a multi-environment structure, but the column containing the environment/location is missing. Please provide heter_groups parameter. For example: heter_groups = 'locations'. If you have location as column name in your phenotypic data."), call. = FALSE)
    }

    if(isFALSE(cross_validation)){
    if ("GBLUP" %in%GS_model) {
      checkVarianceCovarianceInputs(heter_groups, heter_resid, var_cov_str, var_cov_str_available, msg)
    }
    }
    if(isTRUE(cross_validation)){
      if ("GBLUP" %in% GS_model_cv) {
        checkVarianceCovarianceInputs(heter_groups, heter_resid, var_cov_str, var_cov_str_available, msg)
      }
    }

  } else {
    # Reset parameters for single environment
    heter_groups <- NULL
    heter_resid <- NULL
    var_cov_str <- NULL
  }
}



# Function to validate model requirements based on the environment structure
# validateModelRequirements <- function(GS_model, bayes_gblup_valid_models, condition1, condition1_1, condition2, msg) {
#   if (GS_model %in% valid_kernel_models) {
#     if (condition1 != condition1_1 && isTRUE(condition2)) {
#       stop(paste(msg, "To fit a relationship/kernel genomic prediction model, you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel), or an omics kernel. Additionally, you can provide genomic or omics data. Ensure you provide instructions on the genomic relationship matrix method or kernel method to calculate the relationship matrix."), call. = FALSE)
#     }
#   }
# }

# Function to validate model requirements based on the environment structure and cross-validation status
validateModelRequirements <- function(GS_model, GS_model_cv, bayes_gblup_valid_models, condition1, condition1_1, condition2, cross_validation, msg, gp_valid_models = NULL, asreml_model = "GBLUP") {
  # Determine the relevant model variable to use based on cross-validation status
  relevant_model <- ifelse(isTRUE(cross_validation), GS_model_cv, GS_model)
  valid_kernel_models <- c(bayes_gblup_valid_models, asreml_model)
  if (exists("gp_multi_environment_kernel_models", mode = "function")) {
    valid_kernel_models <- gp_multi_environment_kernel_models(
      bayes_gblup_valid_models = bayes_gblup_valid_models,
      asreml_model = asreml_model,
      gp_valid_models = gp_valid_models
    )
  }

  # Check if the relevant model meets the specified requirements
  if (!is.null(relevant_model) && any(relevant_model %in% valid_kernel_models)) {
    if (condition1 != condition1_1 && isTRUE(condition2)) {
      model_type_desc <- ifelse(isTRUE(cross_validation), "cross-validation", "")
      stop(paste(msg, sprintf("To fit a relationship/kernel genomic prediction model %s, you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel), or an omics kernel. Additionally, you can provide genomic or omics data. Ensure you provide instructions on the genomic relationship matrix method or kernel method to calculate the relationship matrix.", model_type_desc)), call. = FALSE)
    }
  }
}

