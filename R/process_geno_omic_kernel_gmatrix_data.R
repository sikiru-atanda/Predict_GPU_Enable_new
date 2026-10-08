# Function to clean and process genomic data
#' Clean and process genomic data into kernels / relationship matrices
#'
#' Wraps the geno QC + kernel / GRM construction pipeline: applies the SNP-call
#' / individual-call / MAF / heterozygosity filters, optional imputation and
#' LD pruning, then builds the requested marker kernels (`kernel_method`) and /
#' or genomic relationship matrices (`gmatrix_method`). Handles a single
#' `geno_data` input or pre-split `train_geno_data` / `test_geno_data`.
#'
#' @param geno_data Genomic-marker matrix or data frame for the full set
#'   (rows = individuals, columns = SNPs).
#' @param train_geno_data Optional training-set genomic matrix.
#' @param test_geno_data Optional test-set genomic matrix matching
#'   `train_geno_data`.
#' @param test_set Optional integer / character vector identifying the test-set
#'   rows (otherwise inferred from missing phenotypes).
#' @param train_set Optional integer / character vector identifying the
#'   training-set rows.
#' @param pheno_clean_list Pre-checked phenotype list (from
#'   [phenotype_precheck]) used to align rows.
#' @param gen_name Name of the genotype-ID column in the phenotype data.
#' @param kernel_method Character vector of marker-kernel methods to calculate.
#'   When more than one method is supplied, the returned `gmatrix` is a named
#'   list of kernel matrices.
#' @param gmatrix_method Character vector of marker relationship-matrix methods
#'   to calculate. When more than one method is supplied, the returned
#'   `gmatrix` is a named list of relationship matrices.
#' @param scale Logical; scale columns to unit variance before computing the
#'   relationship.
#' @param center Logical; mean-centre columns before computing the
#'   relationship.
#' @param map_data Optional marker map (`Chrom`, `Pos`) used by LD pruning.
#' @param maf_threshold Minimum minor-allele frequency (MAF) below which SNPs
#'   are dropped.
#' @param het_threshold Maximum per-SNP heterozygosity above which SNPs are
#'   dropped.
#' @param ind_call_rate_threshold Minimum per-individual call rate; rows below
#'   this are dropped.
#' @param snp_call_rate_threshold Minimum per-SNP call rate; SNPs below this
#'   are dropped.
#' @param impute Logical; impute remaining missing genotype values.
#' @param imputation_method Imputation method when `impute = TRUE`; one of
#'   `"knn"`, `"mean"`, `"median"`, or `"mode"`.
#' @param impute_knn_k Integer; number of nearest neighbours when
#'   `imputation_method = "knn"`.
#' @param ploidy One positive integer or `"auto"`; propagated to QC and GRM
#'   construction.
#' @param qc_filtering Logical; apply the QC filters (MAF / het / call-rate).
#' @param message Logical; if `TRUE`, print summary messages.
#' @param heter_groups Optional heterogeneous-groups column for MET runs.
#' @param ld_prunning_qc Logical; if `TRUE`, apply LD-based pruning after the
#'   QC filters.
#' @param test_set_overrides_inferred Logical. When `TRUE`, the user-supplied
#'   `test_set` overrides inferred testing IDs from missing phenotypes or
#'   geno/omic-only records.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the processed genomic matrix, kernel / GRM(s), and the
#'   resolved train / test indices.
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
                              ploidy = "auto",
                              qc_filtering = NULL,
                              message = TRUE,
                              heter_groups = NULL,
                              test_set_overrides_inferred = FALSE,
                              ld_prunning_qc = TRUE,
                              ...) {

  msg <- ""
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
                                ploidy = if (inherits(geno_data, "list") && !is.null(geno_data$ploidy)) geno_data$ploidy else ploidy,
                                map_data = map_data,
                                qc_filtering = if (inherits(geno_data, "list")) NULL else qc_filtering,
                                message = message,
                                ld_prunning_qc = ld_prunning_qc)

  low_call_rate_inds_removed <- cleaned_data[["low_call_rate_inds_removed"]]
  low_call_rate_inds_removed <- check_test_set(low_call_rate_inds_removed)
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
                                    low_call_rate_inds_removed = low_call_rate_inds_removed,
                                    test_set_overrides_inferred = test_set_overrides_inferred,
                                    message = message)

    if (length(pheno_match) > 1) {
      model_ready <- pheno_match[["geno_pheno_match_data"]]
      test_set <- pheno_match[["test_data"]]
    } else {
      model_ready <- pheno_match[["geno_pheno_match_data"]]
    }
    resolved_ploidy <- cleaned_data$ploidy %||% gp_resolve_matrix_ploidy(
      cleaned_data[["snps_matrix"]], ploidy
    )
    attr(model_ready, "ploidy") <- resolved_ploidy

    if(is.null(kernel_method) && is.null(gmatrix_method)){
      out <- list(geno_model_ready = model_ready,
                  clean_geno_qcstat = cleaned_data,
                  ploidy = resolved_ploidy)
      if (!is.null(test_set)) out[["test_set"]] <- test_set
      if (!is.null(low_call_rate_inds_removed)) out[["low_call_rate_inds_removed"]] <- low_call_rate_inds_removed
      return(out)

    }
    #rm(pheno_match)
  }

  if (!is.null(kernel_method)) {
    kernel <- kernel_calculation(M_matrix_clean = cleaned_data[["snps_matrix"]],
                                 scale = scale,
                                 centering = center,
                                 method = kernel_method,
                                 message = message)

    pheno_match_kernel <- gp_match_kernel_bank(
      kernels = kernel,
      object_pheno = pheno_clean_list[["pheno_clean_data"]],
      gen_name = gen_name,
      test_set = test_set,
      train_set = train_set,
      heter_groups = heter_groups,
      low_call_rate_inds_removed = low_call_rate_inds_removed,
      test_set_overrides_inferred = test_set_overrides_inferred,
      message = message
    )

    kernel <- pheno_match_kernel[["geno_pheno_match_data"]]
    if (!is.null(pheno_match_kernel[["test_data"]])) {
      test_set <- pheno_match_kernel[["test_data"]]
    }
    out <- list(gmatrix= kernel,
                geno_model_ready = model_ready,
                clean_geno_qcstat = cleaned_data,
                ploidy = resolved_ploidy)
    if (!is.null(test_set)) out[["test_set"]] <- test_set
    if (!is.null(low_call_rate_inds_removed)) out[["low_call_rate_inds_removed"]] <- low_call_rate_inds_removed
    return(out)
    #rm(cleaned_data)
  }

  if (!is.null(gmatrix_method)) {
    gmatrix <- grm_calculation(
      geno_clean = cleaned_data[["snps_matrix"]],
      method = gmatrix_method,
      ploidy = cleaned_data$ploidy %||% ploidy
    )

    pheno_match_kernel <- gp_match_kernel_bank(
      kernels = gmatrix,
      object_pheno = pheno_clean_list[["pheno_clean_data"]],
      gen_name = gen_name,
      test_set = test_set,
      train_set = train_set,
      heter_groups = heter_groups,
      low_call_rate_inds_removed = low_call_rate_inds_removed,
      test_set_overrides_inferred = test_set_overrides_inferred,
      message = message
    )

    gmatrix <- pheno_match_kernel[["geno_pheno_match_data"]]
    if (!is.null(pheno_match_kernel[["test_data"]])) {
      test_set <- pheno_match_kernel[["test_data"]]
    }
    out <- list(gmatrix= gmatrix,
                geno_model_ready = model_ready,
                clean_geno_qcstat = cleaned_data,
                ploidy = resolved_ploidy)
    if (!is.null(test_set)) out[["test_set"]] <- test_set
    if (!is.null(low_call_rate_inds_removed)) out[["low_call_rate_inds_removed"]] <- low_call_rate_inds_removed
    return(out)
    #rm(cleaned_data)
  }

}


# Function to clean and process omic data
#' Clean and process omics data into kernels for model fitting
#'
#' Wraps the omics QC + kernel construction pipeline: drops high-missingness
#' columns, optionally imputes the rest, and builds the requested
#' omics kernel(s) for either a single `omic_data` matrix or pre-split
#' `train_omic_data` / `test_omic_data`.
#'
#' @param omic_data Omics matrix for the full set (rows = individuals).
#' @param train_omic_data Optional training omics matrix.
#' @param test_omic_data Optional test-set omics matrix matching
#'   `train_omic_data`.
#' @param pheno_clean_list Pre-checked phenotype list (from
#'   [phenotype_precheck]) used to align rows.
#' @param gen_name Name of the genotype-ID column in the phenotype data.
#' @param test_set Optional integer / character vector identifying the test-set
#'   rows.
#' @param train_set Optional integer / character vector identifying the
#'   training-set rows.
#' @param kernel_method Character vector of omics-kernel methods to calculate.
#'   When more than one method is supplied, the returned `kernel` is a named
#'   list of kernel matrices.
#' @param heter_groups Optional heterogeneous-groups column for MET runs.
#' @param center Logical; mean-centre columns before computing the kernel.
#' @param scale Logical; scale columns to unit variance before computing the
#'   kernel.
#' @param message Logical; if `TRUE`, print summary messages.
#' @param impute_omic Logical; impute remaining missing values.
#' @param imputation_method Imputation method when `impute_omic = TRUE`;
#'   one of `"knn"`, `"mean"`, `"median"`.
#' @param impute_knn_k Integer; number of nearest neighbours when
#'   `imputation_method = "knn"`.
#' @param na_threshold Numeric in `(0, 1]`; columns whose missing-value
#'   proportion exceeds this are dropped.
#' @param test_set_overrides_inferred Logical. When `TRUE`, the user-supplied
#'   `test_set` overrides inferred testing IDs from missing phenotypes or
#'   omic-only records.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list with the processed omics matrix, kernel(s), and the resolved
#'   train / test indices.
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
                              message = TRUE,
                              impute_omic = FALSE,
                              imputation_method = "knn",
                              impute_knn_k = 5,
                              na_threshold = 0.9,
                              test_set_overrides_inferred = FALSE,
                              ...) {

  ##browser()
  msg <- ""
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
                                      test_set_overrides_inferred = test_set_overrides_inferred,
                                      message = message)

      if (length(pheno_match) > 1) {
        model_ready <- pheno_match[["geno_pheno_match_data"]]
        test_set <- pheno_match[["test_data"]]
      } else {
        model_ready <- pheno_match[["geno_pheno_match_data"]]
      }

      if(is.null(kernel_method)){
        out <- list(clean_omic = cleaned_dataa,
                    omic_model_ready = model_ready)
        if (!is.null(test_set)) out[["test_set"]] <- test_set
        return(out)

      }
      #rm(pheno_match, cleaned_data)
    }

    if (!is.null(kernel_method)) {
      kernel <- kernel_calculation(M_matrix_clean = cleaned_dataa,
                                   scale = scale,
                                   centering = center,
                                   method = kernel_method,
                                   message = message)
      #rm(cleaned_data)

      pheno_match_kernel <- gp_match_kernel_bank(
        kernels = kernel,
        object_pheno = pheno_clean_list[["pheno_clean_data"]],
        gen_name = gen_name,
        test_set = test_set,
        train_set = train_set,
        heter_groups = heter_groups,
        test_set_overrides_inferred = test_set_overrides_inferred,
        message = message
      )

      kernel <- pheno_match_kernel[["geno_pheno_match_data"]]
      if (!is.null(pheno_match_kernel[["test_data"]])) {
        test_set <- pheno_match_kernel[["test_data"]]
      }
      out <- list(kernel= kernel,
                  clean_omic = cleaned_dataa,
                  omic_model_ready = model_ready)
      if (!is.null(test_set)) out[["test_set"]] <- test_set
      return(out)
    }


  } else {
    return(NULL)
  }
}



