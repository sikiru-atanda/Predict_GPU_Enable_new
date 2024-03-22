# Adjust the function to accept external variables as parameters
AI_process_ml_data_if_valid <- function(model_check,
                                     geno_omic_model_ready_list,
                                     pheno_clean,
                                     response,
                                     gen_name) {
  if (isTRUE(model_check)) {
    geno_clean <- if ("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL
    # Use grep to find keys that match "omic_model_ready" pattern, including "omic1_model_ready", "omic2_model_ready", etc.
    omic_keys <- grep("omic[0-9]*_model_ready", names(geno_omic_model_ready_list), value = TRUE)
    omic_clean <- lapply(omic_keys, function(key) geno_omic_model_ready_list[[key]])


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




