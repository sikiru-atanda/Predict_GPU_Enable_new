
handle_missing_values <- function(data, na_threshold = 0.9, impute = TRUE, imputation_method = "knn",
                                  impute_knn_k = 5, genotype_ploidy = NULL) {

  if (is.null(impute_knn_k)) impute_knn_k <- 5
  msg <- ""
  if (is.null(rownames(data))) {
    stop(paste(msg, "The feature matrix or data.frame must have row names."), call. = FALSE)
  }
  if (any(is.na(rownames(data)))) {
    stop(paste(msg, "The feature matrix or data.frame row names must not contain NA values."), call. = FALSE)
  }

  is_matrix <- is.matrix(data)
  data_matrix <- if (is_matrix) data else as.matrix(data)
  if (!is.numeric(data_matrix)) {
    stop(paste(msg, "The feature matrix or data.frame contains non-numeric values."), call. = FALSE)
  }

  if (is.null(na_threshold)) na_threshold <- 0.9
  if (!is.null(na_threshold) && na_threshold < 0.5) {
    na_threshold <- 1 - na_threshold
  }

  marker_missing_rate <- if (geno_qc_cpp_available()) {
    geno_qc_metrics_cpp(data_matrix)$marker_missing_rate
  } else {
    colMeans(is.na(data_matrix))
  }
  low_call_rate_snps <- which(marker_missing_rate > na_threshold)

  if (length(low_call_rate_snps) > 0) {
    message(insight::print_color(paste(msg, paste("Removing snps/features with high missing rate:", length(low_call_rate_snps))), "blue"))
    data_matrix <- data_matrix[, -low_call_rate_snps, drop = FALSE]
  }

  if (isTRUE(impute) && anyNA(data_matrix)) {
    imputation_method <- tolower(trimws(imputation_method))
    if (!is.null(genotype_ploidy)) {
      imputed <- gp_impute_alt_dosage(
        dosage = t(data_matrix), ploidy = genotype_ploidy,
        method = imputation_method, k = impute_knn_k
      )
      data_matrix <- t(imputed$dosage)
    } else if (imputation_method %in% c("median", "mean")) {
      if (geno_impute_summary_cpp_available()) {
        data_matrix <- geno_impute_summary_cpp(data_matrix, method = imputation_method)
      } else {
        for (j in seq_len(ncol(data_matrix))) {
          tmp <- data_matrix[, j, drop = TRUE]
          replacement <- if (imputation_method == "median") {
            round(stats::median(tmp, na.rm = TRUE))
          } else {
            round(mean(tmp, na.rm = TRUE))
          }
          data_matrix[is.na(tmp), j] <- replacement
        }
      }
    } else if (imputation_method == "mode") {
      for (j in seq_len(ncol(data_matrix))) {
        tmp <- data_matrix[, j, drop = TRUE]
        observed <- tmp[is.finite(tmp)]
        if (!length(observed)) {
          stop(paste(msg, "Cannot mode-impute a feature with no observed values."), call. = FALSE)
        }
        counts <- table(observed)
        data_matrix[is.na(tmp), j] <- as.numeric(names(counts)[which.max(counts)])
      }
    } else if (imputation_method == "knn") {
      data_matrix <- handle_large_scale_knn(
        data = data_matrix,
        k = impute_knn_k,
        chunk_size = 50,
        num_cores = NULL
      )
    } else {
      stop(paste(msg, "Unknown imputation_method. Use 'knn', 'mean', 'median', or 'mode'."), call. = FALSE)
    }
  }

  if (is_matrix) {
    data_matrix
  } else {
    as.data.frame(data_matrix, check.names = FALSE)
  }
}



#' Quality Control for Genomic Data
#'
#' This function applies quality control (QC) filters to genomic data. It checks and adjusts SNP coding, removes monomorphic markers, and filters markers based on minor allele frequency (MAF), heterozygosity, and call rates for both individuals and SNPs.
#'
#' @param object_geno A matrix or data.frame containing genomic data with individuals in rows and SNPs in columns.
#' @param qc_filtering Logical, if TRUE, quality control filtering is applied based on the thresholds provided.
#' @param maf_threshold Numeric, the threshold for minor allele frequency below which SNPs are removed.
#' @param het_threshold Numeric, the threshold for heterozygosity above which SNPs are removed.
#' @param ind_call_rate_threshold Numeric, the threshold for individual call rate below which individuals are removed.
#' @param snp_call_rate_threshold Numeric, the threshold for SNP call rate below which SNPs are removed.
#' @param impute Logical, if TRUE, missing values are imputed. Currently, this parameter is not used in the function but could be implemented for imputation.
#' @param imputation_method Imputation method for missing genotype values when
#'   `impute = TRUE`; one of `"knn"`, `"mean"`, `"median"`, or `"mode"`.
#' @param impute_knn_k Integer; number of nearest neighbours used when
#'   `imputation_method = "knn"`.
#' @param ploidy One positive integer or `"auto"`. Polyploid numeric matrices
#'   must provide explicit ploidy or carry a `ploidy` attribute.
#' @param ld_prunning_qc Logical; if `TRUE`, apply LD-based pruning of markers
#'   after the QC filters.
#' @param message Logical, if TRUE, messages about the QC process are displayed.
#' @param ... Additional arguments affecting the QC process.
#'
#' @return A list containing:
#'   - \code{snps_matrix}: The genomic data matrix after applying QC filters.
#'   - \code{qc_metrics_and_summary_stat}: A data frame summarizing the QC process, including the number of markers and individuals removed.
#'
#' @examples
#' \dontrun{
#' # Assuming `genomic_data` is a matrix with SNP data
#' qc_results <- geno_precheck(object_geno = genomic_data, qc_filtering = TRUE,
#'                             maf_threshold = 0.05, het_threshold = 0.2,
#'                             ind_call_rate_threshold = 0.9, snp_call_rate_threshold = 0.9)
#' qc_genomic_data <- qc_results$snps_matrix
#' qc_summary <- qc_results$qc_metrics_and_summary_stat
#' }
#'
#' @importFrom stats colMeans
#' @importFrom dplyr rename rownames_to_column
#' @import tibble
#' @export
#'
geno_precheck <- function(object_geno = NULL,
                          qc_filtering = TRUE,
                          maf_threshold = 0.05,
                          het_threshold = 0.2,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          impute=TRUE,
                          imputation_method = "knn", #median, mean
                          impute_knn_k = 5,
                          ploidy = "auto",
                          ld_prunning_qc = TRUE,
                          message = TRUE,
                          ...) {
  # ... (input validation, if necessary)
  msg <- ""
  if (!is.null(object_geno)) {
    if("data.table" %in% class(object_geno)){
      stop(paste(msg,'Genomic data must be data.frame or matrix not character.'), call. = FALSE)
    }
    if(inherits(object_geno, "character")) stop(paste(msg,'Genomic data should be data.frame or matrix not character.'), call. = FALSE)
    if (!is.matrix(object_geno)) {
      object_geno <- as.matrix(object_geno)
    }

    # Check row and column names in object_geno.
    if (is.null(rownames(object_geno)) || is.null(colnames(object_geno))) {
      stop(paste(msg,"Individual or marker names not assigned to rows or columns of 'object_geno'."), call. = FALSE)
      }

    if(any(colnames(object_geno) %in% c("NA", "Na", "na"))){
      stop(paste(msg, "Column names can't contain NA."))
    }

   requested_ploidy <- gp_validate_ploidy(ploidy)
   AA <- if (identical(requested_ploidy, "auto")) {
     detect_genomic_coding(object_geno = object_geno)
   } else {
     "Explicit ALT dosage"
   }
   if(AA=="SNP (0, 1, 2, -1)") {
     stop(paste(msg,"SNP recoding is decoded wrongly as (0,1,2-1).\n Snp recode should either be SNP: (-1, 0, 1) or (0, 1, 2)."), call. = FALSE)

   } else if(AA=="SNP (-1, 0, 1)"){
     object_geno <- object_geno + 1
     if(isTRUE(message)) {
       message(insight::print_color(paste(msg, "The allele dosages are not in 0, 1, 2.\n We Fix it to required format."), "blue"))

     }
     gc()
     }else if(AA=="SNP (0, 0.5, 1)"){
       # Convert 0.5 to 1, and 1 to 2
       object_geno[object_geno == 1] <- 2
       object_geno[object_geno == 0.5] <- 1

       if(isTRUE(message)) {
         message(insight::print_color(paste(msg, "The allele dosages are not in 0, 1, 2.\n We Fix it to required format."), "blue"))

       }
       gc()
     } else{
     if(AA=="Presence/Absence (0, 1)"){
       # Convert 1s to 2s
       object_geno[object_geno == 1] <- 2
       if(isTRUE(message)) {
         message(insight::print_color(paste(msg, "The allele dosages are not in 0, 2.\n We Fix it to required format."), "blue"))
       }

     }
       gc()
     }

   ploidy <- gp_resolve_matrix_ploidy(object_geno, requested_ploidy)
   gp_validate_alt_dosage(object_geno, ploidy, hard_calls = FALSE)
   object_geno <- gp_attach_ploidy(object_geno, ploidy)

   if (!is.numeric(object_geno)) {
     stop(paste(msg, 'The geno data contains non-numeric values'), call. = FALSE)
   }
    # Check if allele dosage are not in  0, 1, 2 format but -1, 0, 1 format
  #   check_geno <- which(object_geno == -1)
  #   if (length(check_geno) != 0) {
  #     if(isTRUE(message)) {
  #       message("The allele dosages are not in 0, 1, 2. Fixing it.")
  #     }
  #     object_geno <- object_geno + 1
  #   } else {
  #
  #     ### Then check if the matrix is presence and absence where presence is 1 and absence is 0
  #     ## change it to 0, 2
  #   ## Check if it not coded 0,1 instead 0, 2
  #   unique_values <- unique(object_geno)
  #   if(is.character(unique_values)=="character"){
  #     unique_values = as.double(unique_values)
  #   }
  #
  #   if (all(unique_values %in% c(0, 1))) {
  #     object_geno[object_geno == 1] <- 2
  #   }
  #
  # }
    #######
    ### Set initial value for qc parameters
    markers_callrate_removed  <-  0
    ind_callrate_removed  <-  0
    het_markers_removed  <-  0
    maf_markers_removed <-  0
    low_call_rate_inds_removed <- NULL
    ####
    ###
    # Remove monomorphic markers
    qc_metrics <- geno_qc_metrics(object_geno, ploidy = ploidy)
    monomorphic_markers <- qc_metrics$monomorphic
    if (length(monomorphic_markers) > 0) {
      if(isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing monomorphic markers:", length(monomorphic_markers))), "blue"))
      }

        object_geno <- object_geno[, -monomorphic_markers, drop = FALSE]
        qc_metrics <- NULL


      #map_data <- map_data[-monomorphic_markers, ]
      total_mono <-  length(monomorphic_markers)
      rm(monomorphic_markers); gc()
    } else {
      if(isTRUE(message)) {
        message(insight::print_color(paste(msg, "No monomorphic markers to remove."), "blue"))
      }
      total_mono <-  0
    }
    ###
    if(isTRUE(qc_filtering)){
    # Minor alele frequency
    if (!is.null(maf_threshold)) {
      qc_metrics <- geno_qc_metrics(object_geno, ploidy = ploidy)
      MAF <- qc_metrics$maf
      maf_markers <- which(MAF < maf_threshold)

      if (length(maf_markers) > 0) {
        object_geno <- object_geno[, -maf_markers, drop = FALSE]
        maf_markers_removed <- length(maf_markers)
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing markers with MAF below threshold:", maf_markers_removed)), "blue"))
        }
      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No markers removed based on MAF threshold."), "blue"))
        }
      }
      rm(MAF, maf_markers, qc_metrics); gc()
    }

    ######

    ### Remove markers with high missing value based on desired thresold
    if(!is.null(snp_call_rate_threshold)) {

      # Calculate individual call rate
      qc_metrics <- geno_qc_metrics(object_geno, ploidy = ploidy)
      snp_call_rate <- qc_metrics$marker_missing_rate

      # Filter based on individual call rate
      #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
      low_call_rate_snps <- which(snp_call_rate > (1 - snp_call_rate_threshold))
      if (length(low_call_rate_snps) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing SNPs with low call rate:", length(low_call_rate_snps))), "blue"))
        }
        object_geno <- object_geno[, -low_call_rate_snps, drop = FALSE]
        #hapmap <- hapmap[low_call_rate_snps, ]
        #map_data <- map_data[low_call_rate_snps, ]
        markers_callrate_removed  <-  length(low_call_rate_snps)
      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No SNPs removed based on SNP call rate threshold."), "blue"))
        }

        markers_callrate_removed  <-  0
      }
      rm(snp_call_rate, low_call_rate_snps, qc_metrics); gc()
    }

    ##
    ### Remove markers with high missing value based on desired thresold
    if(!is.null(ind_call_rate_threshold)) {
      # Calculate individual call rate
      qc_metrics <- geno_qc_metrics(object_geno, ploidy = ploidy)
      ind_call_rate <- qc_metrics$individual_missing_rate

      # Filter based on SNP call rate
      #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
      low_call_rate_inds <- which(ind_call_rate > (1 - ind_call_rate_threshold))

      if (length(low_call_rate_inds) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing Individuals with low call rate:", length(low_call_rate_inds))), "blue"))
        }

        low_call_rate_inds_removed <- rownames(object_geno)[low_call_rate_inds]
        object_geno <- object_geno[-low_call_rate_inds, , drop = FALSE]

        ind_callrate_removed <- length(low_call_rate_inds)

      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No individuals removed based on call rate threshold."), "blue"))
        }

        ind_callrate_removed <-  0
      }

      rm(ind_call_rate, qc_metrics); gc()
    }
    ##
    # Calculate heterozygosity
    if(!is.null(het_threshold)){
      qc_metrics <- geno_qc_metrics(object_geno, ploidy = ploidy)
      heteroz <- qc_metrics$heterozygosity

      # Filter based on heterozygosity
      #het_markers <- which(heteroz >= het_threshold)
      het_markers <- which(heteroz > het_threshold)
      if (length(het_markers) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing markers with high heterozygosity:", length(het_markers))), "blue"))
        }
        object_geno <-  object_geno[, -het_markers, drop = FALSE]
        het_markers_removed  <- length(het_markers)
        #snp_data <- snp_data[-het_markers, ]

      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No markers removed based on heterozygosity threshold."), "blue"))
        }

        het_markers_removed <- 0
      }
      rm(het_markers, heteroz, qc_metrics); gc()

    }

  }
   ####

  } else {
    stop(paste(msg,"Marker/snp data cannot be empty."), call. = FALSE)

  }

  if (!NCOL(object_geno)) {
    stop(paste(
      msg,
      "Genotype QC removed every marker (removed: call rate", markers_callrate_removed,
      "| heterozygosity", het_markers_removed, "| MAF", maf_markers_removed,
      "| monomorphic", total_mono, "). Check the thresholds: het_threshold =",
      format(het_threshold %||% "NULL"), "is an inbred-line filter (set het_threshold = NULL",
      "for heterozygous or outbred material), maf_threshold =", format(maf_threshold %||% "NULL"), "."
    ), call. = FALSE)
  }

  object_geno <- handle_missing_values(data = object_geno,
                                       impute=impute, imputation_method = imputation_method,
                                       impute_knn_k = impute_knn_k,
                                       genotype_ploidy = ploidy)
  gp_validate_alt_dosage(object_geno, ploidy, hard_calls = FALSE)
  object_geno <- gp_attach_ploidy(object_geno, ploidy)
  ## check if there is duplicated snps
  duplicated_columns <- colnames(object_geno)[duplicated(colnames(object_geno))]
  if(length(duplicated_columns)>0){
    stop(paste(msg,"Marker/snp data contain duplicate snps."), call. = FALSE)
  }

  ## check duplicated rownames:
  duplicated_rownames <- rownames(object_geno)[duplicated(rownames(object_geno))]
  if(length(duplicated_rownames)>0){
    stop(paste(msg,"Marker/snp data contain duplicate genotypes."), call. = FALSE)
  }
  if(isTRUE(ld_prunning_qc)){
  keep_prunned_snp <- ld_prune_graph(object_geno, ploidy = ploidy)

  object_geno <- object_geno[, keep_prunned_snp, drop=FALSE]
  object_geno <- gp_attach_ploidy(object_geno, ploidy)

  }

  #### Aggregate all the maker data
  #########################
  # A NULL threshold means that filter was switched off; record it as NA so
  # the one-row summary can still be built.
  summary_stat_snp= data.frame(
    snp_call_rate_threshold = snp_call_rate_threshold %||% NA_real_,
    total_snp_removed = markers_callrate_removed,
    ind_call_rate_threshold = ind_call_rate_threshold %||% NA_real_,
    total_genotypes_removed = ind_callrate_removed,
    het_threshold = het_threshold %||% NA_real_,
    total_het_snps_removed = het_markers_removed,
    maf_threshold = maf_threshold %||% NA_real_,
    total_maf_snps_remove = maf_markers_removed,
    total_monomorphic_snps_removed = total_mono,
    ploidy = ploidy,
    stringsAsFactors = FALSE
  )

  summary_stat_snp <- summary_stat_snp |>
    t() |>
    as.data.frame()

  summary_stat_snp <- summary_stat_snp |>
    dplyr::mutate(metric = rownames(summary_stat_snp)) |>
    dplyr::select(metric, V1) |>
    dplyr::rename(stat = V1)

  # summary_stat_snp <- summary_stat_snp |>
  #   t() |>
  #   as.data.frame() |>
  #   tibble::rownames_to_column(var = "metric") |>
  #   dplyr::rename(stat = V1)
rownames(summary_stat_snp) <- NULL
  ######

  attr(object_geno, "cleared") <- "pass"
  attr(object_geno, "ploidy") <- ploidy

  return(list(
    snps_matrix = object_geno,
    qc_metrics_and_summary_stat = summary_stat_snp,
    low_call_rate_inds_removed = low_call_rate_inds_removed,
    ploidy = ploidy
  )
)


  }
