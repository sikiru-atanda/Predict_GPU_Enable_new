# Function to clean and process genomic data
#' Title
#'
#' @param geno_data
#' @param train_geno_data
#' @param test_geno_data
#' @param kernel_method
#' @param gmatrix_method
#' @param ...
#' @param test_set
#' @param train_set
#' @param gen_name
#' @param scale
#' @param map_data
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param qc_filtering
#' @param message
#' @param pheno_clean_list
#'
#' @return
#' @export
#'
#' @examples
process_geno_data <- function(geno_data = NULL,
                              train_geno_data = NULL,
                              test_geno_data = NULL,
                              test_set = NULL,
                              train_set = NULL,
                              pheno_clean_list = NULL,
                              gen_name = NULL,
                              kernel_method = NULL,
                              gmatrix_method = NULL,
                              scale = NULL,
                              map_data = NULL,
                              maf_threshold = NULL,
                              het_threshold = NULL,
                              ind_call_rate_threshold = NULL,
                              snp_call_rate_threshold = NULL,
                              impute = NULL,
                              qc_filtering = NULL,
                              message = NULL,
                              ...) {

  msg <- sprintf("==================================================\n")

# genomic data check ------------------------------------------------------

  cleaned_data <- geno_to_model(geno_data = geno_data,
                                train_geno_data = train_geno_data,
                                test_geno_data = test_geno_data,
                                maf_threshold = maf_threshold,
                                het_threshold = het_threshold,
                                ind_call_rate_threshold = ind_call_rate_threshold,
                                snp_call_rate_threshold = snp_call_rate_threshold,
                                impute = impute,
                                map_data = map_data,
                                qc_filtering = qc_filtering,
                                message = message)

  if(attr(cleaned_data[["snps_matrix"]], "cleared")!="for_model_fit" && all(class(cleaned_data[["snps_matrix"]])!=c("matrix", "array"))) {
    stop(print(paste(msg,'Data is not fit for model')), call. = FALSE)

  }

  #if ((exists('cleaned_data') & exists("pheno_clean"))) {
  if (!is.null(pheno_clean_list)) {
    pheno_match <- pheno_geno_match(object_pheno = pheno_clean_list[["pheno_clean_data"]],
                                    object_geno = cleaned_data[["snps_matrix"]],
                                    gen_name = gen_name,
                                    test_set = test_set,
                                    train_set = train_set,
                                    message = message)

    if (length(pheno_match) > 1) {
      model_ready <- pheno_match[["geno_pheno_match_data"]]
      test_set <- pheno_match[["test_data"]]
    } else {
      model_ready <- pheno_match[["geno_pheno_match_data"]]
    }

    if(is.null(kernel_method) && is.null(gmatrix_method)){
      return(list(geno_model_ready = model_ready,
                  clean_geno_qcstat = cleaned_data))

    }
    #rm(pheno_match)
  }

  if (!is.null(kernel_method)) {
    kernel <- kernel_calculation(M_matrix_clean = cleaned_data[["snps_matrix"]],
                                 scale = scale,
                                 method = kernel_method,
                                 message = message)

    return(list(gmatrix= kernel,
                geno_model_ready = model_ready,
                clean_geno_qcstat = cleaned_data))
    #rm(cleaned_data)
  }

  if (!is.null(gmatrix_method)) {
    gmatrix <- grm_calculation(geno_clean = cleaned_data[["snps_matrix"]],
                               method = gmatrix_method)

    return(list(gmatrix= gmatrix,
                geno_model_ready = model_ready,
                clean_geno_qcstat = cleaned_data))
    #rm(cleaned_data)
  }

}


# Function to clean and process omic data
#' Title
#'
#' @param omic_data
#' @param train_omic_data
#' @param test_omic_data
#' @param kernel_method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
process_omic_data <- function(omic_data = NULL,
                              train_omic_data = NULL,
                              test_omic_data = NULL,
                              kernel_method = NULL,
                              pheno_clean_list = NULL,
                              gen_name = NULL,
                              test_set = NULL,
                              train_set = NULL,
                              ...) {
  msg <- sprintf("==================================================\n")
  if(isFALSE(((is.null(omic_data) & is.null(train_omic_data)) & is.null(test_omic_data)))){
    #if (!is.null(data) && !is.null(train_data) && !is.null(test_data)) {
    cleaned_dataa <- omic_to_model(omic_data = omic_data,
                                   train_omic_data = train_omic_data,
                                   test_omic_data = test_omic_data,
                                   message = message)

    if (attr(cleaned_dataa, "cleared") != "for_model_fit" && all(class(cleaned_dataa) != c("matrix", "array"))) {
      stop(print(paste(msg, 'Data is not fit for model')), call. = FALSE)
    }

    #if ((exists('cleaned_data') & exists("pheno_clean"))) {
    if (!is.null(pheno_clean_list)) {
      pheno_match <- pheno_geno_match(object_pheno = pheno_clean_list[["pheno_clean_data"]],
                                      object_geno = cleaned_dataa,
                                      gen_name = gen_name,
                                      test_set = test_set,
                                      train_set = train_set,
                                      message = message)

      if (length(pheno_match) > 1) {
        model_ready <- pheno_match[["geno_pheno_match_data"]]
        test_set <- pheno_match[["test_data"]]
      } else {
        model_ready <- pheno_match[["geno_pheno_match_data"]]
      }

      if(is.null(kernel_method)){
        return(list(clean_omic = cleaned_dataa,
                    omic_model_ready = model_ready))

      }
      #rm(pheno_match, cleaned_data)
    }

    if (!is.null(kernel_method)) {
      kernel <- kernel_calculation(M_matrix_clean = cleaned_dataa,
                                   scale = scale,
                                   method = kernel_method,
                                   message = message)
      #rm(cleaned_data)
      return(list(kernel= kernel,
                  clean_omic = cleaned_dataa,
                  omic_model_ready = model_ready))
    }


  } else {
    return(NULL)
  }
}



