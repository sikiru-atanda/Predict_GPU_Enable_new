
handle_missing_values <- function(data, na_threshold = 0.2) {

  # Check if the input has row names
  if (is.null(rownames(data))) {
    stop("The feature matrix or data.frame must have row names.")
  }

  # Check for NA values in row names
  if (any(is.na(rownames(data)))) {
    stop("The feature matrix or data.frame row names must not contain NA values.")
  }
  # Check if the input is a matrix or a data frame
  is_matrix <- is.matrix(data)

  # Convert matrix to data frame for easier handling
  if (is_matrix) {
    data <- as.data.frame(data)
  }

  # Calculate the proportion of NA values in each column
  na_proportion <- colMeans(is.na(data))

  # Identify columns to remove (proportion of NAs > threshold)
  cols_to_remove <- na_proportion > na_threshold

  if (any(cols_to_remove)) {
    sapply(names(na_proportion[cols_to_remove]), function(col) {
      cat(sprintf("Column '%s' has %.2f%% NAs, which is above the threshold. Removing column.\n",
                  col, na_proportion[col] * 100))
    })
  }

  # Remove the columns with high proportion of NAs
  data <- data[, !cols_to_remove, drop = FALSE]

  # Impute remaining NAs with column medians
  cols_to_impute <- names(na_proportion[!cols_to_remove & na_proportion > 0])
  if (length(cols_to_impute) > 0) {
    sapply(cols_to_impute, function(col) {
      cat(sprintf("Column '%s' has %.2f%% NAs, which is within the threshold. Imputing NAs with median.\n",
                  col, na_proportion[col] * 100))
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
    })
  }

  # Convert back to matrix if the original input was a matrix
  if (is_matrix) {
    data <- as.matrix(data)
  }

  return(data)
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
#' @param message Logical, if TRUE, messages about the QC process are displayed.
#' @param ... Additional arguments affecting the QC process.
#'
#' @return A list containing:
#'   - \code{snps_matrix}: The genomic data matrix after applying QC filters.
#'   - \code{qc_metrics_and_summary_stat}: A data frame summarizing the QC process, including the number of markers and individuals removed.
#'
#' @examples
#' # Assuming `genomic_data` is a matrix with SNP data
#' qc_results <- geno_precheck(object_geno = genomic_data, qc_filtering = TRUE,
#'                             maf_threshold = 0.05, het_threshold = 0.2,
#'                             ind_call_rate_threshold = 0.9, snp_call_rate_threshold = 0.9)
#' qc_genomic_data <- qc_results$snps_matrix
#' qc_summary <- qc_results$qc_metrics_and_summary_stat
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
                          impute = TRUE,
                          message = TRUE,
                          ...) {
  # ... (input validation, if necessary)
  msg <- "\n==================================================\n"
  if (!is.null(object_geno)) {
    if("data.table" %in% class(object_geno)){
      stop(print(paste(msg,'Genomic data must be data.frame or matrix not character.')), call. = FALSE)
    }
    if(inherits(object_geno, "character")) stop(print(paste(msg,'Genomic data should be data.frame or matrix not character.')), call. = FALSE)
    if (!is.matrix(object_geno)) {
      object_geno <- as.matrix(object_geno)
    }

    # Check row and column names in object_geno.
    if (is.null(rownames(object_geno)) || is.null(colnames(object_geno))) {
      stop(print(paste(msg,"Individual or marker names not assigned to rows or columns of 'object_geno'.")), call. = FALSE)
      }

   AA <-  detect_genomic_coding(object_geno = object_geno)
   if(AA=="SNP (0, 1, 2, -1)") {
     stop(print(paste(msg,"SNP recoding is decoded wrongly as (0,1,2-1).\n Snp recode should either be SNP: (-1, 0, 1) or (0, 1, 2).")), call. = FALSE)

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

   all_numeric <- all(apply(object_geno, c(1, 2), is.numeric))

   # Stop execution if any element is not numeric
   if (!all_numeric) {
     stop('The geno data contains non-numeric values', call. = FALSE)
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
    ####
    ###
    # Remove monomorphic markers
    monomorphic_markers <- which(apply(object_geno, 2, function(x) length(table(x)) <= 1))
    if (length(monomorphic_markers) > 0) {
      if(isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing monomorphic markers:", length(monomorphic_markers))), "blue"))
      }

        object_geno <- object_geno[, -monomorphic_markers]      
        

      #map_data <- map_data[-monomorphic_markers, ]
      total_mono <-  length(monomorphic_markers)
      rm(monomorphic_markers); gc()
    } else {
      if(isTRUE(message)) {
        message(insight::print_color(paste(msg, "No monomorphic markers to remove."), "blue"))
      }
      total_mono = 0
    }
    ###
    if(isTRUE(qc_filtering)){
    # Minor alele frequency
    if (!is.null(maf_threshold)) {
      phat <- colMeans(object_geno, na.rm = TRUE) / 2
      MAF <- ifelse(phat < 0.5, phat, 1 - phat)

      if (any(MAF < maf_threshold)) {
        object_geno <- object_geno[, -which(MAF < maf_threshold)]
        maf_markers_removed <- length(which(MAF < maf_threshold))
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing markers with MAF below threshold:", maf_markers_removed)), "blue"))
        }
      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No markers removed based on MAF threshold."), "blue"))
        }
      }
    }

    ######

    ### Remove markers with high missing value based on desired thresold
    if(!is.null(snp_call_rate_threshold)) {

      # Calculate individual call rate
      snp_call_rate <- colMeans(is.na(object_geno))

      # Filter based on individual call rate
      #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
      low_call_rate_snps <- which(snp_call_rate > snp_call_rate_threshold)
      if (length(low_call_rate_snps) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing SNPs with low call rate:", length(low_call_rate_snps))), "blue"))
        }
        object_geno <- object_geno[, -low_call_rate_snps]
        #hapmap <- hapmap[low_call_rate_snps, ]
        #map_data <- map_data[low_call_rate_snps, ]
        markers_callrate_removed  <-  length(low_call_rate_snps)
      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No SNPs removed based on SNP call rate threshold."), "blue"))
        }

        markers_callrate_removed  <-  0
      }
      rm(snp_call_rate, low_call_rate_snps); gc()
    }

    ##
    ### Remove markers with high missing value based on desired thresold
    if(!is.null(ind_call_rate_threshold)) {
      # Calculate individual call rate
      ind_call_rate <- rowMeans(is.na(object_geno))

      # Filter based on SNP call rate
      #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
      low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)

      if (length(low_call_rate_inds) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing Individuals with low call rate:", length(low_call_rate_inds))), "blue"))
        }
        object_geno <- object_geno[-low_call_rate_inds, ]
      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No individuals removed based on call rate threshold."), "blue"))
        }

        ind_callrate_removed <-  0
      }

      rm(ind_call_rate, low_call_rate_inds); gc()
    }
    ##
    # Calculate heterozygosity
    if(!is.null(het_threshold)){
      heteroz <- apply(object_geno, 2, function(x) sum( x== 1, na.rm=T) / length(x))

      # Filter based on heterozygosity
      #het_markers <- which(heteroz >= het_threshold)
      het_markers <- which(heteroz > het_threshold)
      if (length(het_markers) > 0) {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, paste("Removing markers with high heterozygosity:", length(het_markers))), "blue"))
        }
        object_geno <-  object_geno[, -het_markers]
        het_markers_removed  <- length(het_markers)
        rm(het_markers, heteroz); gc()
        #snp_data <- snp_data[-het_markers, ]

      } else {
        if(isTRUE(message)) {
          message(insight::print_color(paste(msg, "No markers removed based on heterozygosity threshold."), "blue"))
        }

        het_markers_removed <- 0
      }

    }

  }
   ####

  } else {
    stop(print(paste(msg,"Marker/snp data cannot be empty.")), call. = FALSE)

  }

  object_geno <- handle_missing_values(data = object_geno)
  #### Aggregate all the maker data
  #########################
  summary_stat_snp= data.frame(
    snp_call_rate_threshold = snp_call_rate_threshold,
    total_snp_removed = markers_callrate_removed,
    ind_call_rate_threshold = ind_call_rate_threshold,
    total_genotypes_removed = ind_callrate_removed,
    het_threshold = het_threshold,
    total_het_snps_removed = het_markers_removed,
    maf_threshold = maf_threshold,
    total_maf_snps_remove = maf_markers_removed,
    total_monomorphic_snps_removed = total_mono,
    stringsAsFactors = FALSE
  )

  summary_stat_snp <- summary_stat_snp |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "metric") |>
    dplyr::rename(stat = V1)

  ######

  attr(object_geno, "cleared") <- "pass"

  return(list(
    snps_matrix = object_geno,
    qc_metrics_and_summary_stat = summary_stat_snp
  )
)


  }
