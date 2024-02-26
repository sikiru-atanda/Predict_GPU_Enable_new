#vcf_file_name = "2021_NDSU_AYT.vcf",
#vcf_file_path = "D:/PredictProR"
#' Title
#'
#' @param vcf_file_name
#' @param vcf_file_path
#' @param vcf_file
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param recode_format
#' @param out_put_map
#' @param message
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom data.table := .SD .SDcols lapply
#'
vcf_qc_recodeOLD <- function(vcf_file_name = NULL,
                          vcf_file_path = NULL,
                          vcf_file = NULL,
                          maf_threshold = 0.01,
                          het_threshold = 0.2,
                          ind_call_rate_threshold = 0.2,
                          snp_call_rate_threshold = 0.3,
                          #hwe_threshold = 0.001,  # Adjust as needed
                          impute = FALSE,
                          recode_format = "0,1,2",  # Specify "0,1,2" for 0/0, 0/1, 1/1 or "-1,0,1" for -1, 0, 1
                          out_put_map = TRUE,
                          #beagle_path = "D:/PredictProR",
                          message = TRUE) {

  if(!is.null(vcf_file_name) && !is.null(vcf_file_path)){
    # Construct the full file path
    full_file_path <- file.path(vcf_file_path, vcf_file_name)
    # Read the file using fread
    vcf_file <- data.table::fread(full_file_path,
                                  header = TRUE,
                                  skip = "#",
                                  check.names = FALSE,
                                  na.strings=c('./.', '.|.', 'NA', '.'))


  } else if (!is.null(vcf_file)){
### if user provide vcf_file it possible missing data which is usually
    ## ./. or .|. might not be converted to NA
    ## This used lappy with functionality in data.table
    vcf_file[, (10:ncol(vcf_file)) := lapply(.SD, function(x) {
      x[which(x %in% c('./.','.|.', 'NA', '.'))] <- NA
      return(x)
    }), .SDcols = 10:ncol(vcf_file)]


  }else{
    if(is.null(vcf_file)){

      stop("vcf file is missing.")

    }

  }


  if(isFALSE(data.table::is.data.table(vcf_file))){

    vcf_file <- data.table::as.data.table(vcf_file)
  }

  # Extract SNP data
  # snp_data <- vcf_data[, 10:ncol(vcf_data), with = FALSE]
  # #rownames(snp_data) <- as.character(vcf_data[, 3][[1]])
  # map_data <- vcf_data[, 1:9, with = FALSE]
  # rm(vcf_data); gc(); cat('0\14')

  # Remove monomorphic markers
  monomorphic_markers <- which(apply(vcf_file[, 10:ncol(vcf_file)], 1, function(x) length(table(x)) <= 1))
  if (length(monomorphic_markers) > 0) {
    if (message) {
      print(paste("Removing monomorphic markers:", length(monomorphic_markers)))
    }

    vcf_file <- vcf_file[-monomorphic_markers, ]
    #snp_data <- snp_data[-monomorphic_markers, ]
    total_mono <- length(monomorphic_markers)
  } else {
    if (message) {
      print("No monomorphic markers to remove.")
    }
    total_mono <-  0
  }


  # Calculate minor allele frequency
  allele_freq <- apply(vcf_file[, 10:ncol(vcf_file)], 1, function(x) {
    allele_counts <- table(x)
    minor_allele_count <- min(allele_counts)
    major_allele_count <- max(allele_counts)
    maf <- minor_allele_count / sum(allele_counts)
    return(maf)
  })

  # Filter based on minor allele frequency
  #maf_markers <- which(allele_freq >= maf_threshold)
  maf_markers <-  which(allele_freq < maf_threshold)
  if (length(maf_markers) > 0) {
    if (message) {
      print(paste("Removing markers with MAF below threshold:", length(maf_markers)))
    }

    vcf_file <-  vcf_file[-maf_markers, ]
    #snp_data <- snp_data[maf_markers, ]
    maf_markers_removed <- length(maf_markers)

  } else {
    if (message) {
      print("No markers removed based on MAF threshold.")
    }
    maf_markers_removed <- 0
  }

  # Calculate heterozygosity
  heteroz <- apply(vcf_file[, 10:ncol(vcf_file)], 1, function(x) sum(x %in% c("0|1", "1|0", "0/1", "1/0")) / length(x))

  # Filter based on heterozygosity
  #het_markers <- which(heteroz >= het_threshold)
  het_markers <- which(heteroz > het_threshold)
  if (length(het_markers) > 0) {
    if (message) {
      print(paste("Removing markers with high heterozygosity:", length(het_markers)))
    }
    vcf_file <-  vcf_file[-het_markers, ]
    het_markers_removed  <- length(het_markers)
    #snp_data <- snp_data[-het_markers, ]

  } else {
    if (message) {
      print("No markers removed based on heterozygosity threshold.")
    }

    het_markers_removed <- 0
  }

  # Calculate individual call rate
  snp_call_rate <- rowMeans(is.na(vcf_file[, 10:ncol(vcf_file)]))

  # Filter based on individual call rate
  #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
  low_call_rate_snps <- which(snp_call_rate > snp_call_rate_threshold)
  if (length(low_call_rate_snps) > 0) {
    if (message) {
      print(paste("Removing SNPs with low call rate:", length(low_call_rate_snps)))
    }
    vcf_file <-  vcf_file[-low_call_rate_snps, ]
    markers_callrate_removed <- length(low_call_rate_snps)
    #snp_data <- snp_data[low_call_rate_snps, ]
  } else {
    if (message) {
      print("No SNPs removed based on SNP call rate threshold")
    }

    markers_callrate_removed <- 0
  }


  # Calculate individual call rate
   ind_call_rate <- colMeans(is.na(vcf_file[, 10:ncol(vcf_file)]))

  # Filter based on SNP call rate
  #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
   low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)
   low_call_rate_inds <- low_call_rate_inds+9
  if (length(low_call_rate_inds) > 0) {
    if (message) {
      print(paste("Removing Individuals with low call rate:", length(low_call_rate_inds)))
    }
    vcf_file <- vcf_file[, -low_call_rate_inds, with = FALSE]

    #snp_data <- snp_data[, -low_call_rate_inds, with = FALSE]
    #snp_data <- snp_data[, low_call_rate_inds, with = FALSE]
  } else {
    if (message) {
      print("No individuals removed based on call rate threshold.")
    }
  }

  #################################
   ## Record to numeric
   if (recode_format %in% c("0,1,2", "-1,0,1")) {
     vcf_file[, (10:ncol(vcf_file)) := lapply(.SD, function(x) {
       if (recode_format == "0,1,2") {
         x[x %in% c("0|0", "0/0")] <- 0
         x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- 1
         x[x %in% c("1|1", "1/1")] <- 2
         x[x %in% c("./.", "NA", "NA/NA")] <- NA
       } else if (recode_format == "-1,0,1") {
         x[x %in% c("0|0", "0/0")] <- -1
         x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- 0
         x[x %in% c("1|1", "1/1")] <- 1
         x[x %in% c("./.", "NA", "NA/NA")] <- NA

       }
       x = as.double(as.character(unlist(x)))
       return(x)
     }), .SDcols = 10:ncol(vcf_file)]


   } else {
     stop("Invalid recode format. Use '0,1,2' or '-1,0,1'.")
   }
  ####
  # Impute missing values

  if(impute){
    if(any(is.na(vcf_file[, 10:ncol(vcf_file)]))) {

      vcf_file[, (10:ncol(vcf_file)) := lapply(.SD, function(x) {
        x = as.double(as.character(unlist(x)))
        #x = ifelse(is.na(x), round(mean(x, na.rm = TRUE)))
        x = ifelse(is.na(x), round(mean(x, na.rm = TRUE)), x)
        return(x)
      }), .SDcols = 10:ncol(vcf_file)]

      #This function will impute the missing values
      # for(j in 1:nrow(snp_data_numeric)){
      #   if (message) {
      #     cat("Imputing missing values marker...", "\n")
      #   }
      #   tmp <- snp_data_numeric[j, ]
      #   snp_data_numeric[j, ] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
      # }

    }

  }
#######
   ##################
   ## Aggregate all the info
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
   if(isTRUE(out_put_map)){
     map <-   vcf_file[, 1:9]

     vcf_file <- t(vcf_file[, 10:ncol(vcf_file)])
     colnames(vcf_file) <- as.character(map[, 3][[1]]) ### it assumed the snps name are in col 3
     vcf_file <- as.matrix(vcf_file)
     if(!is.null(recode_format)){
       class(vcf_file) <- c("matrix", "array", "genotype")
     }

     return(list(snps_matrix= vcf_file,
                 snp_map = map,
                 snp_call_rate_threshold = snp_call_rate_threshold,
                 total_snp_removed = markers_callrate_removed,
                 ind_call_rate_threshold = ind_call_rate_threshold,
                 total_genotypes_removed = ind_callrate_removed,
                 het_threshold = het_threshold,
                 total_het_snps_removed = het_markers_removed,
                 maf_threshold = maf_threshold,
                 total_maf_snps_remove = maf_markers_removed,
                 total_monomorphic_snps_removed = total_mono
     )
     )


   } else if(isFALSE(out_put_map)){
     map <-   vcf_file[, 1:9]

     vcf_file <- t(vcf_file[, 10:ncol(vcf_file)])
     colnames(vcf_file) <- as.character(map[, 3][[1]]) ### it assumed the snps name are in col 3
     vcf_file <- as.matrix(vcf_file)
     if(!is.null(recode_format)){
       class(vcf_file) <- c("matrix", "array", "genotype")
     }
     rm(map); gc()

     return(list(snps_matrix= vcf_file,
                 snp_call_rate_threshold = snp_call_rate_threshold,
                 total_snp_removed = markers_callrate_removed,
                 ind_call_rate_threshold = ind_call_rate_threshold,
                 total_genotypes_removed = ind_callrate_removed,
                 het_threshold = het_threshold,
                 total_het_snps_removed = het_markers_removed,
                 maf_threshold = maf_threshold,
                 total_maf_snps_remove = maf_markers_removed,
                 total_monomorphic_snps_removed = total_mono
     )
     )



   }
}

##############################
#########################
# Calculate Hardy-Weinberg Equilibrium p-values
# hwe_p_values <- apply(snp_data, 2, function(x) {
#   obs_counts <- table(x)
#   if (length(obs_counts) == 2) {  # Only for biallelic SNPs
#     obs_frequencies <- obs_counts / sum(obs_counts)
#     obs_hwe_stat <- sum((obs_frequencies - c(1, 2, 1))^2 / c(1, 2, 1))
#     obs_hwe_p_value <- 1 - pchisq(obs_hwe_stat, df = 1)
#     return(obs_hwe_p_value)
#   } else {
#     return(NA)  # Skip non-biallelic SNPs
#   }
# })

# Calculate Hardy-Weinberg equilibrium
# hwe_p_values <- apply(snp_data, 2, function(x) {
#   geno_counts <- table(x)
#   hwe_test <- chisq.test(c(geno_counts[1], geno_counts[2], geno_counts[3]))
#   return(hwe_test$p.value)
# })
#
# # Filter based on Hardy-Weinberg Equilibrium p-value
# hwe_filtered_snps <- which(hwe_p_values <= hwe_threshold)
# if (length(hwe_filtered_snps) > 0) {
#   if (message) {
#     print(paste("Removing SNPs not in Hardy-Weinberg Equilibrium:", length(hwe_filtered_snps)))
#   }
#   snp_data <- snp_data[, hwe_filtered_snps, with = FALSE]
# } else {
#   if (message) {
#     print("No SNPs removed based on Hardy-Weinberg Equilibrium.")
#   }
# }
################

#########################
# # Write SNP data to a file in PLINK format (you may need to adapt this)
# write.table(snp_data, file = "snp_data.plink", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#
# # Prepare BEAGLE parameter file (you may need to adapt this)
# beagle_params <- "beagle.params"
# writeLines(c(
#   "marker = snp_data.plink",
#   "impute = true",
#   "output = imputed_data",
#   "nthreads = 4",    # Set the number of threads based on your system
#   "burnin = 10",
#   "phase-its = 10",
#   "gprobs = true"
#   # Add other BEAGLE parameters as needed
# ), beagle_params)
#
#
# # Call BEAGLE for imputation
# beagle_cmd <- paste("java -jar \"D:/PredictProR/beagle.22Jul22.46e.jar\"", beagle_params, sep = " ")
#
# system(beagle_cmd)
#
# # Read imputed data back into R (you may need to adapt this)
# imputed_data <- fread("imputed_data.bgl.phased")
#########################

# Recode to numeric format
# if (recode_format %in% c("0,1,2", "-1,0,1")) {
#   snp_data_numeric <- apply(vcf_file[, 10:ncol(vcf_file)], c(1, 2), function(x) {
#     if (recode_format == "0,1,2") {
#       x[x %in% c("0|0", "0/0")] <- 0
#       x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- 1
#       x[x %in% c("1|1", "1/1")] <- 2
#       x[x %in% c("./.", "NA", "NA/NA")] <- NA
#     } else if (recode_format == "-1,0,1") {
#       x[x %in% c("0|0", "0/0")] <- -1
#       x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- 0
#       x[x %in% c("1|1", "1/1")] <- 1
#       x[x %in% c("./.", "NA", "NA/NA")] <- NA
#     }
#     return(as.numeric(x))
#   })
# } else {
#   stop("Invalid recode format. Use '0,1,2' or '-1,0,1'.")
# }
#########
## imput missing
# if (impute) {
#   if (any(is.na(snp_data_numeric))) {
#     if (message) {
#       cat("Imputing missing values...\n")
#     }
#     snp_data_numeric <- apply(snp_data_numeric, 1, function(x) ifelse(is.na(x), round(mean(x, na.rm = TRUE)), x))
#   }
# }
####
