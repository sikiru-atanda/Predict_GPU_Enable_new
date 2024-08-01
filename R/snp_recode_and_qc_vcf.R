
# Function to compute call rate
#' Title
#'
#' @param vcf_data
#' @param ind_call_rate_threshold
#'
#' @return
#' @export
#'
#' @examples
compute_call_rate <- function(vcf_data,
                              ind_call_rate_threshold) {

  if(!is.null(ind_call_rate_threshold)){
    total_genotypes_removed <- 0
    ind_call_rate <- colMeans(is.na(vcf_data[, 10:ncol(vcf_data)]))

    # Filter based on SNP call rate
    #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
    low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)
    low_call_rate_inds <- low_call_rate_inds+9
    if (length(low_call_rate_inds) > 0) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing Individuals with low call rate:", length(low_call_rate_inds))), "blue"))

      }
      vcf_data <- vcf_data[, -low_call_rate_inds, with = FALSE]

      total_genotypes_removed <- length(low_call_rate_inds)

      #snp_data <- snp_data[, -low_call_rate_inds, with = FALSE]
      #snp_data <- snp_data[, low_call_rate_inds, with = FALSE]
    } else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No individuals removed based on call rate threshold."), "blue"))

      }


    }

  }
  return(list(vcf_data = vcf_data,
              total_genotypes_removed = total_genotypes_removed))
}

# Calculate individual call rate



# Function to filter heterozygous genotypes based on a threshold
#' Title
#'
#' @param vcf_data
#' @param het_threshold
#'
#' @return
#' @export
#'
#' @examples
filter_heterozygous <- function(vcf_data,
                                het_threshold) {

  het_markers_removed <- 0

  heteroz <- apply(vcf_data[, 10:ncol(vcf_data)], 1, function(x) sum(x %in% c("0|1", "1|0", "0/1", "1/0")) / length(x))
  het_markers <- which(heteroz > het_threshold)

  if (length(het_markers) > 0) {
  vcf_data <-  vcf_data[-het_markers, ]

  het_markers_removed  <- length(het_markers)
  rm(het_markers, heteroz); gc()

  }

  return(list(vcf_data= vcf_data,
              het_markers_removed = het_markers_removed))
}

# Function to recode genotypes
recode_genotypes <- function(vcf_data, recode_format) {
  if(!is.null(recode_format)){
    if (recode_format %in% c("0,1,2", "-1,0,1")) {
      for (j in 10:ncol(vcf_data)) {
        vcf_data[[j]] <- lapply(vcf_data[[j]], function(x) {
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
        })
      }

    }
  }
  return(vcf_data)
}

# Function to process a single batch
#' Title
#'
#' @param vcf_data
#' @param het_threshold
#' @param recode_format
#'
#' @return
#' @export
#'
#' @examples
process_batch <- function(vcf_data,
                          het_threshold,
                          recode_format) {


  # Compute individual call rates
  #vcf_data <- compute_call_rate(vcf_data, ind_call_rate_threshold)

  # Filter heterozygous genotypes
  if(!is.null(het_threshold)) {
    vcf_dataa <- filter_heterozygous(vcf_data, het_threshold)
    vcf_data <- vcf_dataa[["vcf_data"]]
    het_markers_removed <- vcf_dataa[["het_markers_removed"]]

    rm(vcf_dataa); gc()

  } else {
    het_markers_removed <- 0
  }

  # Recode the genotypes
  if(!is.null(recode_format)) vcf_data <- recode_genotypes(vcf_data, recode_format)

  return(list(vcf_data = vcf_data,
              het_markers_removed = het_markers_removed))
}
##############
# Function to process VCF in batches using parallel processing
process_vcf_in_batches_parallel <- function(vcf_file_name_qc_plink,
                                            batch_size = 2000,
                                            het_threshold = 0.2,
                                            ind_call_rate_threshold = 0.9,
                                            recode_format = "0,1,2",
                                            num_cores = NULL) {
  # Read the entire VCF file header
  header_lines <- readLines(vcf_file_name_qc_plink)
  #header <- header_lines[grepl("^#", header_lines)]
  #col_header <- header_lines[grepl("^#CHROM", header_lines)]
  # Extract column names from the header
  #col_names <- strsplit(col_header, "\t")[[1]]
  header <- header_lines[grepl("^##", header_lines)]

  # Dynamically find the line that starts with a single # and contains "chrom"
  col_header <- header_lines[grepl("^#[^#]", header_lines, ignore.case = TRUE) & grepl("chrom", header_lines, ignore.case = TRUE)]

  # Extract column names from the header
  col_names <- strsplit(col_header, "\t")[[1]]


  # Count the number of variant lines
  total_lines <- length(header_lines) - length(header) - 1


  if(total_lines>batch_size){
    # Generate batch start points
    batch_starts <- seq(1, total_lines, by = batch_size)

    # Ensure last batch is handled correctly
    if (total_lines %% batch_size != 0) {
      batch_starts <- c(batch_starts, total_lines + 1)
    }


    #batch_starts <- batch_starts[1:10]
    sys_name <- Sys.info()["sysname"]
    plan_type <- ifelse(sys_name == "Windows", "multisession", "multicore")

    if (total_lines>batch_size) {
      if(is.null(num_cores)){

        num_cores <-  parallel::detectCores()
        num_cores <- num_cores*0.7
      }
      future::plan(plan_type, workers = num_cores)
    }

    # Parallel processing
    results <- future.apply::future_lapply(seq_along(batch_starts), function(i) {

      start_line <- batch_starts[i]
      skip_lines <- start_line + length(header) - 1

      if (i < length(batch_starts)) {
        nrows <- batch_size
      } else {
        nrows <- total_lines - batch_starts[i - 1] + 1
      }

      if(nrows!=0){

        vcf_data <- data.table::fread(vcf_file_name_qc_plink,
                                      skip = skip_lines,
                                      nrows = nrows,
                                      #header = TRUE,
                                      fill=TRUE,
                                      check.names = FALSE,
                                      na.strings=c('./.', '.|.', 'NA', '.'))

        colnames(vcf_data) <- col_names

      }

      # Process the batch
      batch_process <- process_batch(vcf_data,het_threshold, recode_format)

      return(batch_process)
    }, future.seed = TRUE)

    future::plan("sequential")

    # Combine all processed batches
    #results <- Filter(function(df) nrow(df) > 0, results)

    # Combine all processed batches
    #combined_geno <- do.call(rbind, results)

    # Extract vcf_data and filter out empty batches
    vcf_data_list <- lapply(results, function(res) res$vcf_data)
    vcf_data_list <- Filter(function(df) nrow(df) > 0, vcf_data_list)

    # Combine all processed vcf_data batches
    combined_geno <- do.call(rbind, vcf_data_list)

    # Sum het_markers_removed, handling NULL values
    het_markers_removed_list <- lapply(results, function(res) res$het_markers_removed)
    total_het_markers_removed <- sum(unlist(het_markers_removed_list), na.rm = TRUE)
    if(is.null(total_het_markers_removed)) total_het_markers_removed <- 0


  } else{
    vcf_data <- data.table::fread(vcf_file_name_qc_plink,
                                  skip = "#",
                                  header = TRUE,
                                  check.names = FALSE,
                                  na.strings=c('./.', '.|.', 'NA', '.'))

    results <- process_batch(vcf_data,het_threshold, recode_format)

    combined_geno <- results[["vcf_data"]]
    total_het_markers_removed <- results[["total_het_markers_removed"]]

    rm(vcf_data); gc()
  }

  #combined_geno <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  posindex <- grep("pos", colnames(combined_geno), ignore.case = TRUE)

  # Check if the position column is found
  if(length(posindex) != 0){
    pos_column <- colnames(combined_geno)[posindex]

    # Remove duplicated rows based on the position column
    combined_geno <- combined_geno[!duplicated(combined_geno[[pos_column]]), ]
  }

  combined_geno <- compute_call_rate(combined_geno, ind_call_rate_threshold)

  vcf_processed_data <- combined_geno[["vcf_data"]]
  total_genotypes_removed <- combined_geno[["total_genotypes_removed"]]


  return(list(vcf_processed_data = vcf_processed_data,
              total_het_markers_removed = total_het_markers_removed,
              total_genotypes_removed = total_genotypes_removed))
}



#' Quality Control and Recoding for VCF Data
#'
#' Performs quality control filters on Variant Call Format (VCF) data, including filtering based on minor allele frequency (MAF),
#' heterozygosity rate, individual and SNP call rates, and optionally recodes the genotype data.
#'
#' @param vcf_file_name Name of the VCF file (optional).
#' @param vcf_file_path Path to the directory containing the VCF file (optional).
#' @param vcf_file An object containing the VCF data (optional).
#' @param maf_threshold Threshold for minor allele frequency below which SNPs will be removed.
#' @param het_threshold Threshold for heterozygosity above which SNPs will be removed.
#' @param ind_call_rate_threshold Threshold for individual call rate below which individuals will be removed.
#' @param snp_call_rate_threshold Threshold for SNP call rate below which SNPs will be removed.
#' @param impute Logical indicating whether missing genotypes should be imputed.
#' @param recode_format String specifying the format for recoding genotypes. Can be "0,1,2" for homozygous reference,
#' heterozygous, and homozygous alternate, respectively, or "-1,0,1" for an alternative coding scheme.
#' @param out_put_map Logical indicating whether to output the SNP map along with the recoded data.
#' @param message Logical indicating whether messages about the QC process should be displayed.
#'
#' @details The function supports reading VCF data from a file or directly from an R object. It applies several QC
#' filters based on user-defined thresholds and can recode the genotype data into numeric format for further analysis.
#' The function allows for the imputation of missing data and the inclusion of a SNP map in the output.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{snps_matrix}{Matrix of the recoded genotype data.}
#'   \item{snp_map}{Data frame of the SNP map, if \code{out_put_map} is TRUE.}
#'   \item{qc_metrics_and_summary_stat}{Data frame summarizing the QC metrics and the number of SNPs/individuals removed.}
#' }
#'
#' @examples
#' # Assuming `vcf_data` is your VCF data frame
#' result <- vcf_qc_recode(vcf_file = vcf_data, maf_threshold = 0.05, het_threshold = 0.15,
#'                         ind_call_rate_threshold = 0.9, snp_call_rate_threshold = 0.9,
#'                         recode_format = "0,1,2", out_put_map = TRUE, message = TRUE)
#' @export
#' @importFrom data.table := .SD .SDcols lapply

vcf_qc_recode <-   function(vcf_file_name = NULL,
                           vcf_file_path = NULL,
                           #vcf_file = NULL,
                           maf_threshold = 0.05,
                           het_threshold = 0.2,
                           ind_call_rate_threshold = 0.9,
                           snp_call_rate_threshold = 0.9,
                           #hwe_threshold = 0.001,  # Adjust as needed
                           impute = FALSE,
                           recode_format = "0,1,2",  # Specify "0,1,2" for 0/0, 0/1, 1/1 or "-1,0,1" for -1, 0, 1
                           out_put_map = TRUE,
                           batch_size = 2000,
                           num_cores = NULL,
                           #beagle_path = "D:/PredictProR",
                           message = TRUE) {

  msg <- "\n==================================================\n"

  #### Place holders
  markers_callrate_removed  <-  0
  ind_callrate_removed  <-  0
  het_markers_removed  <-  0
  maf_markers_removed <-  0

  main_dir <- getwd()
  if(is.null(vcf_file_path)) vcf_file_path <- getwd()

  if(!is.null(vcf_file_name) && !is.null(vcf_file_path)){
    # Construct the full file path
    setwd(vcf_file_path)
    full_file_path <- file.path(vcf_file_path, vcf_file_name)
    systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_name <- paste("vcf_process", systime, sep = "_")

    #########################################
    # # Read the entire VCF file header
    # header_lines <- readLines(vcf_file_name)
    # header <- header_lines[grepl("^##", header_lines)]
    #
    # # Count the number of variant lines
    # total_lines <- length(header_lines) - length(header) - 1

    #############################################
    # if(total_lines<=5000){
    #   vcf_file <- data.table::fread(vcf_file_name,
    #                                 header = TRUE,
    #                                 skip = "#",
    #                                 check.names = FALSE,
    #                                 na.strings=c('./.', '.|.', 'NA', '.'))
    #
    #
    # }

  res <-   run_plink_qc(input_file = vcf_file_name,
                       output_name = output_name,
                       output_format = "vcf",
                       vcf_file_path = vcf_file_path,
                       remove_monomorphic = TRUE,
                       maf_threshold = maf_threshold,
                       heterozygosity = NULL,
                       snp_call_rate = snp_call_rate_threshold,
                       allow_extra_chr = TRUE,
                       individual_call_rate = NULL)

  if(res == "Failed"){
    stop(paste(msg,paste(paste("Processing of the vcf file failed.",
                               "Ensure you correctly specify the file that contain the vcf file",
                               " to the 'vcf_file_name'",
                               "Also, ensure you provide the right name of the vcf file to the 'vcf_file_name'."))), call. = FALSE)

  }

    log_content <- readLines(paste(output_name, "log", sep = "."))
    # Extract the number of variants before QC
    num_variants_before_qc <- as.numeric(sub(".*: ([0-9]+) variants scanned.*", "\\1", grep("variants scanned", log_content, value = TRUE)))

    # Extract the number of variants removed due to missing genotype data
    markers_callrate_removed <- as.numeric(sub(".*: ([0-9]+) variants removed due to missing genotype data.*", "\\1", grep("variants removed due to missing genotype data", log_content, value = TRUE)))

    # Extract the number of variants removed due to MAF threshold
    maf_markers_removed <- as.numeric(gsub("[^0-9]", "", sub(".* ([0-9]+) variants removed due to allele frequency threshold.*", "\\1", grep("variants removed due to allele frequency threshold", log_content, value = TRUE))))

    # Extract the number of variants remaining after QC
    num_variants_after_qc <- as.numeric(gsub("[^0-9]", "", sub(".* ([0-9]+) variants remaining after main filters.*", "\\1", grep("variants remaining after main filters", log_content, value = TRUE))))


    # Read the file using fread
    # vcf_file <- data.table::fread(full_file_path,
    #                               header = TRUE,
    #                               skip = "#",
    #                               check.names = FALSE,
    #                               na.strings=c('./.', '.|.', 'NA', '.'))


  }else{

    stop(paste(msg,paste("vcf file is missing.", "Provide the vcf_file_name",
                          "and vcf_file_path if different from the current working directory.")), call. = FALSE)


  }

  if(!is.null(snp_call_rate_threshold)){
    if (snp_call_rate_threshold < 0.5) {
      snp_call_rate_threshold <- 1 - snp_call_rate_threshold
    }

  }

  if(!is.null(ind_call_rate_threshold)){
    if (ind_call_rate_threshold < 0.5) {
      ind_call_rate_threshold <- 1 - ind_call_rate_threshold
    }

  }

  vcf_processed <-  process_vcf_in_batches_parallel(vcf_file_name_qc_plink = paste(output_name, "vcf", sep = "."),
                                                    batch_size = batch_size,
                                                    het_threshold = het_threshold,
                                                    ind_call_rate_threshold = ind_call_rate_threshold,
                                                    recode_format = recode_format,
                                                    num_cores = num_cores)



  vcf_processed_data <- vcf_processed[["vcf_processed_data"]]
  total_het_snps_removed <- vcf_processed[["total_het_markers_removed"]]
   total_genotypes_removed <- vcf_processed[["total_genotypes_removed"]]
   if(is.null(total_het_snps_removed) || is.na(total_het_snps_removed)) total_het_snps_removed <- 0
   if(is.null(total_genotypes_removed) || is.na(total_genotypes_removed)) total_genotypes_removed <- 0
  ##################

   if(isTRUE(impute)){
     #if(any(is.na(snp_data[, 12:ncol(snp_data)]))==T) {
     if(any(is.na(vcf_processed_data[, 10:ncol(vcf_processed_data)]))==T) {
       #This function will impute the missing values
       #for(j in 12:ncol(snp_data)){
       for(j in 10:ncol(vcf_processed_data)){
         tmp <- vcf_processed_data[,j, with=FALSE]
         #tmp <- vcf_processed_data[,j]
         tmp = as.double(as.character(unlist(tmp)))
         vcf_processed_data[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
       }

     }

   }

  ## Aggregate all the info
  #### Aggregate all the maker data
  summary_stat_snp= data.frame(
    snp_call_rate_threshold = snp_call_rate_threshold,
    total_snp_removed = markers_callrate_removed,
    ind_call_rate_threshold = ind_call_rate_threshold,
    total_genotypes_removed = total_genotypes_removed,
    het_threshold = het_threshold,
    total_het_snps_removed = total_het_snps_removed,
    maf_threshold = maf_threshold,
    total_maf_snps_remove = maf_markers_removed,
    #total_monomorphic_snps_removed = total_mono,
    stringsAsFactors = FALSE
  )

  summary_stat_snp <- summary_stat_snp |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "metric") |>
    dplyr::rename(stat = V1)

  rm(vcf_processed); gc()
  #########################
  if(isTRUE(out_put_map)){
    map <-   vcf_processed_data[, 1:9]
    map <-  as.data.frame(map)
    #colNAMES <- as.character(map[, 3][[1]])
    vcf_processed_data <- t(vcf_processed_data[, 10:ncol(vcf_processed_data)])
    vcf_processed_data <- as.matrix(vcf_processed_data)
    colnames(vcf_processed_data) <- as.character(map[, 3]) ### it assumed the snps name are in col 3

    if (!is.numeric(vcf_processed_data)) {
      # Apply as.numeric to each element in the matrix
      vcf_processed_data <- apply(vcf_processed_data, c(1, 2), as.numeric)
    }

    if(!is.null(recode_format)){
      class(vcf_processed_data) <- c("matrix", "array", "genotype")
    }

    return(list(snps_matrix= vcf_processed_data,
                snp_map = map,
                qc_metrics_and_summary_stat = summary_stat_snp

    )
    )


  } else {


    if(isFALSE(out_put_map)){
    map <-   vcf_processed_data[, 1:9]
    map <- as.data.frame(map)
    #colNAMES <- as.character(map[, 3][[1]])
    vcf_processed_data <- t(vcf_processed_data[, 10:ncol(vcf_processed_data)])
    #colnames(vcf_file) <- as.character(map[, 3][[1]]) ### it assumed the snps name are in col 3
    vcf_processed_data <- as.matrix(vcf_processed_data)
    colnames(vcf_processed_data) <-  as.character(map[, 3]) ### it assumed the snps name are in col 3
    #####
    if (!is.numeric(vcf_processed_data)) {
      # Apply as.numeric to each element in the matrix
      vcf_processed_data <- apply(vcf_processed_data, c(1, 2), as.numeric)
    }


    if(!is.null(recode_format)){
      class(vcf_processed_data) <- c("matrix", "array", "genotype")
    }
    rm(map); gc()

    return(list(snps_matrix= vcf_processed_data,
                qc_metrics_and_summary_stat = summary_stat_snp

    )
    )



  }
}

}
