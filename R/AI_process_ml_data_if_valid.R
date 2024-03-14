# Adjust the function to accept external variables as parameters
AI_process_ml_data_if_valid <- function(model_check,
                                     geno_omic_model_ready_list,
                                     pheno_clean,
                                     response,
                                     gen_name) {
  if (isTRUE(model_check)) {
    geno_clean <- if ("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL
    omic_clean <- if (!"geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL

    # Call ML_data_processing with the necessary parameters
    ML_data_processing(pheno_clean = pheno_clean,
                       response = response,
                       gen_name = gen_name,
                       geno_clean = geno_clean,
                       omic_clean = omic_clean)
  } else {
    list() # Return an empty list if conditions are not met
  }
}




