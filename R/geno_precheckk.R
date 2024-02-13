


#' Title
#'
#' @param object_geno
#' @param qc_filtering
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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

  if (!is.null(object_geno)) {
    if (!is.matrix(object_geno)) {
      object_geno <- as.matrix(object_geno)
    }

    # Check row and column names in object_geno.
    if (is.null(rownames(object_geno)) || is.null(colnames(object_geno))) {
      stop("Individual or marker names not assigned to rows or columns of 'object_geno'.")
    }

    # Check if allele dosage are not in  0, 1, 2 format but -1, 0, 1 format
    check_geno <- which(object_geno == -1)
    if (length(check_geno) != 0) {
      if (message) {
        message("The allele dosages are not in 0, 1, 2. Fixing it.")
      }
      object_geno <- object_geno + 1
    } else {

      ### Then check if the matrix is presence and absence where presence is 1 and absence is 0
      ## change it to 0, 2
    ## Check if it not coded 0,1 instead 0, 2
    unique_values <- unique(object_geno)
    if(is.character(unique_values)=="character"){
      unique_values = as.double(unique_values)
    }

    if (all(unique_values %in% c(0, 1))) {
      object_geno[object_geno == 1] <- 2
    }

  }
    ###
    # Remove monomorphic markers
    monomorphic_markers <- which(apply(object_geno, 2, function(x) length(table(x)) <= 1))
    if (length(monomorphic_markers) > 0) {
      if (message) {
        print(paste("Removing monomorphic markers:", length(monomorphic_markers)))
      }
      object_geno <- object_geno[, -monomorphic_markers, ]
      #map_data <- map_data[-monomorphic_markers, ]
      total_mono = length(monomorphic_markers)
      rm(monomorphic_markers); gc()
    } else {
      if (message) {
        print("No monomorphic markers to remove.")
      }
      total_mono = 0
    }
    ###
    if(isTRUE(qc_filtering)){
    # Minor alele frequency
    if (!is.null(maf_threshold)) {
      phat <- colMeans(object_geno, na.rm = TRUE) / 2
      MAF <- ifelse(phat < 0.5, phat, 1 - phat)

      if (message) {
        message("Minor allele frequency check (MAF).")
      }

      if (any(MAF < maf_threshold)) {
        object_geno <- object_geno[, -which(MAF < maf_threshold)]
        maf_markers_removed <- length(which(MAF < maf_threshold))
        if (message) {
          print(paste("Removing markers with MAF below threshold:", length(which(MAF < maf_threshold))))
        }
      } else {
        if (message) {
          message("No loci with MAF below threshold.")
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
        if (message) {
          print(paste("Removing SNPs with low call rate:", length(low_call_rate_snps)))
        }
        object_geno <- object_geno[, -low_call_rate_snps]
        #hapmap <- hapmap[low_call_rate_snps, ]
        #map_data <- map_data[low_call_rate_snps, ]
        markers_callrate_removed  = length(low_call_rate_snps)
      } else {
        if (message) {
          print("No SNPs removed based on SNP call rate threshold")
        }

        markers_callrate_removed  = 0
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
        if (message) {
          print(paste("Removing Individuals with low call rate:", length(low_call_rate_inds)))
        }
        object_geno <- object_geno[-low_call_rate_inds, ]
      } else {
        if (message) {
          print("No individuals removed based on call rate threshold.")
        }

        ind_callrate_removed = 0
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
        if (message) {
          print(paste("Removing markers with high heterozygosity:", length(het_markers)))
        }
        object_geno <-  object_geno[, -het_markers]
        het_markers_removed  <- length(het_markers)
        rm(het_markers, heteroz); gc()
        #snp_data <- snp_data[-het_markers, ]

      } else {
        if (message) {
          print("No markers removed based on heterozygosity threshold.")
        }

        het_markers_removed <- 0
      }

    }

  }
   ####

  } else {
    stop("Marker/snp data cannot be empty.")
  }

  #### Aggregate all the maker data

  if(!exists("markers_callrate_removed")) {
    markers_callrate_removed  = 0
  }
  if(!exists("ind_callrate_removed ")) {
    ind_callrate_removed  = 0
  }
  if(!exists("het_markers_removed ")) {
    het_markers_removed  = 0
  }
  if(!exists("maf_markers_removed")) {
    maf_markers_removed = 0
  }
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
