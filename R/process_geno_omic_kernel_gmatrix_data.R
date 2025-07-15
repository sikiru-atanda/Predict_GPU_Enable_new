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
                              scale = FALSE,
                              center = TRUE,
                              map_data = NULL,
                              maf_threshold = NULL,
                              het_threshold = NULL,
                              ind_call_rate_threshold = NULL,
                              snp_call_rate_threshold = NULL,
                              impute = FALSE,
                              imputation_method = "knn",
                              impute_knn_k = 5,
                              qc_filtering = NULL,
                              message = TRUE,
                              heter_groups = heter_groups,
                              ld_prunning_qc = TRUE,
                              ...) {

  msg <- "\n==================================================\n"
  low_call_rate_inds_removed <- NULL
# genomic data check ------------------------------------------------------
### geno_data will be a list when user supplied vcf/hampmap and it is recorded in the engine
  cleaned_data <- geno_to_model(geno_data = if (inherits(geno_data, "list")) geno_data[["snps_matrix"]] else geno_data,
                                train_geno_data = train_geno_data,
                                test_geno_data = test_geno_data,
                                maf_threshold = maf_threshold,
                                het_threshold = het_threshold,
                                ind_call_rate_threshold = ind_call_rate_threshold,
                                snp_call_rate_threshold = snp_call_rate_threshold,
                                impute = impute,
                                imputation_method = imputation_method,
                                impute_knn_k = impute_knn_k,
                                map_data = map_data,
                                qc_filtering = if (inherits(geno_data, "list")) NULL else qc_filtering,
                                message = message,
                                ld_prunning_qc = ld_prunning_qc)

  low_call_rate_inds_removed <- cleaned_data[["low_call_rate_inds_removed"]]
  ### This import the geno_qc from QC and recoding and add it for the final
  ## qc_metrics_and_summary_stat when raw snp data is provided
  if (inherits(geno_data, "list")) {
    metric_removed <- c("markers_callrate_removed",
                        "ind_callrate_removed",
                        "het_markers_removed",
                        "maf_markers_removed")
    index_metric_removed <- which(!(rownames(cleaned_data[["qc_metrics_and_summary_stat"]]) %in% metric_removed))

    cleaned_data[["qc_metrics_and_summary_stat"]] <- rbind(cleaned_data[["qc_metrics_and_summary_stat"]][index_metric_removed, ],
                                                           geno_data[["qc_metrics_and_summary_stat"]])
  }


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
                                    heter_groups = heter_groups,
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
                                 center = center,
                                 method = kernel_method,
                                 message = message)

    pheno_match_kernel <- pheno_geno_match(object_pheno = pheno_clean_list[["pheno_clean_data"]],
                                    object_geno = kernel,
                                    gen_name = gen_name,
                                    test_set = test_set,
                                    train_set = train_set,
                                    heter_groups = heter_groups,
                                    message = message)

    kernel <- pheno_match_kernel[["geno_pheno_match_data"]]

    return(list(gmatrix= kernel,
                geno_model_ready = model_ready,
                clean_geno_qcstat = cleaned_data))
    #rm(cleaned_data)
  }

  if (!is.null(gmatrix_method)) {
    gmatrix <- grm_calculation(geno_clean = cleaned_data[["snps_matrix"]],
                               method = gmatrix_method)

    pheno_match_kernel <- pheno_geno_match(object_pheno = pheno_clean_list[["pheno_clean_data"]],
                                           object_geno = gmatrix,
                                           gen_name = gen_name,
                                           test_set = test_set,
                                           train_set = train_set,
                                           heter_groups = heter_groups,
                                           message = message)

    gmatrix <- pheno_match_kernel[["geno_pheno_match_data"]]

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
                              heter_groups = NULL,
                              center = TRUE,
                              scale = FALSE,
                              message = message,
                              impute_omic = FALSE,
                              imputation_method = "knn",
                              impute_knn_k = 5,
                              na_threshold = 0.9,
                              ...) {

  ##browser()
  msg <- "\n==================================================\n"
  if(isFALSE(((is.null(omic_data) & is.null(train_omic_data)) & is.null(test_omic_data)))){
    #if (!is.null(data) && !is.null(train_data) && !is.null(test_data)) {
    cleaned_dataa <- omic_to_model(omic_data = omic_data,
                                   train_omic_data = train_omic_data,
                                   test_omic_data = test_omic_data,
                                   message = message,
                                   impute_omic = impute_omic,
                                   imputation_method = imputation_method,
                                   impute_knn_k = impute_knn_k,
                                   na_threshold = na_threshold)

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
                                      heter_groups = heter_groups,
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
                                   center = center,
                                   method = kernel_method,
                                   message = message)
      #rm(cleaned_data)

      pheno_match_kernel <- pheno_geno_match(object_pheno = pheno_clean_list[["pheno_clean_data"]],
                                             object_geno = kernel,
                                             gen_name = gen_name,
                                             test_set = test_set,
                                             train_set = train_set,
                                             heter_groups = heter_groups,
                                             message = message)

      kernel <- pheno_match_kernel[["geno_pheno_match_data"]]
      return(list(kernel= kernel,
                  clean_omic = cleaned_dataa,
                  omic_model_ready = model_ready))
    }


  } else {
    return(NULL)
  }
}



