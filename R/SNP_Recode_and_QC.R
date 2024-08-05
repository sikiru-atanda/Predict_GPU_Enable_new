

#' Convert Double Nucleotide Code to IUPAC Single-Letter Code
#'
#' Converts a double nucleotide code to its corresponding IUPAC single-letter code.
#'
#' @param double_code A string representing a double nucleotide code (e.g., "AA", "AT").
#'
#' @return A single-letter IUPAC code corresponding to the double nucleotide code.
#'
#' @examples
#' double_code_to_IUPAC("AA") # Returns "A"
#' double_code_to_IUPAC("AT") # Returns "W"
#'

double_code_to_IUPAC <- function(double_code) {
  iupac_map <- list(
    "AA" = "A", "TT" = "T",
    "CC" = "C", "GG" = "G",
    "AT" = "W", "TA" = "W",
    "CG" = "S", "GC" = "S",
    "AC" = "M", "CA" = "M",
    "GT" = "K", "TG" = "K",
    "AG" = "R", "GA" = "R",
    "TC" = "Y", "CT" = "Y",
    "TG" = "B", "GT" = "B",
    "TC" = "B", "CT" = "B",
    "AG" = "D", "GA" = "D",
    "AT" = "D", "TA" = "D",
    "AT" = "H", "TA" = "H",
    "AC" = "H", "CA" = "H",
    "AC" = "V", "CA" = "V",
    "CG" = "V", "GC" = "V",
    "NN" = "N"
  )
  return(iupac_map[[double_code]])
}

#' Standardize SNP Format
#'
#' Standardizes SNP format by converting double nucleotide codes to IUPAC single-letter codes.
#'
#' @param snp A string representing a SNP in either single or double format.
#'
#' @return A standardized SNP in IUPAC single-letter code.
#'
#' @examples
#' standardize_snp_format("A/T") # Returns "W"
#' standardize_snp_format("AA")  # Returns "A"
#'

standardize_snp_format <- function(snp) {
  if (is.na(snp)) {
    return(NA)
  }
  # Remove potential separator characters and convert to uppercase
  #snp <- gsub("[:/]", "", toupper(snp))
  snp <- gsub("[:/_]", "", toupper(snp))
  # If the SNP is a double code, convert it to IUPAC code
  if (nchar(snp) == 2) {
    snp <- double_code_to_IUPAC(snp)
  }
  return(snp)
}

#' Standardize SNP Format in HapMap Data Frame
#'
#' Standardizes SNP format in a HapMap data frame to be IUPAC compatible.
#'
#' @param hapmap A data frame containing HapMap genotype data starting from the 12th column.
#'
#' @return A HapMap data frame with standardized SNP format.
#'
#' @examples
#' hapmap <- data.frame(rs = 1:3, alleles = c("A/T", "C/G", "A/A"), X1 = c("AA", "CC", "AT"), X2 = c("TA", "GC", "AA"))
#' IUPAC_hapmap_compatible(hapmap)
#'

IUPAC_hapmap_compatible <- function(hapmap) {

  # Identify the SNP columns (starting from the 12th column)
  snp_cols <- names(hapmap)[12:ncol(hapmap)]

  # Standardize SNP format for each SNP column
  hapmap <- as.data.frame(hapmap)
  snp_data <- apply(hapmap[, colnames(hapmap)%in%snp_cols], 2, function(col) sapply(col, standardize_snp_format))
  rownames(snp_data) <- NULL
  hapmap <- cbind(hapmap[, 1:11], snp_data)
  hapmap <- data.table::as.data.table(hapmap)
  rm(snp_data)
  rm(snp_cols)
  gc()

  return(hapmap)

}

#' Remove Multiallelic Markers
#'
#' Removes markers with more than three alleles from a HapMap data frame.
#'
#' @param hapmap A data frame containing HapMap genotype data.
#'
#' @return A HapMap data frame with multiallelic markers removed.
#'
#' @examples
#' hapmap <- data.frame(rs = 1:3, alleles = c("A/T/G", "C/G", "A/A"), X1 = c("A", "C", "A"), X2 = c("T", "G", "A"))
#' remove_multiallelic_markers(hapmap)
#'
remove_multiallelic_markers <- function(hapmap) {
  # Ensure the input is a data.table
  #setDT(hapmap)

  # Identify the alleles column, regardless of capitalization
  alleles_col <- grep("^alleles$", names(hapmap), ignore.case = TRUE, value = TRUE)

  if (length(alleles_col) == 0) {
    stop("Alleles column not found")
  }

  # Function to check if a marker has more than three alleles
  has_more_than_three_alleles <- function(alleles) {
    return(length(unique(unlist(strsplit(as.character(alleles), split = "/")))) > 2)
  }

  # Apply the function to filter out markers with more than three alleles
  hapmap <- hapmap[!sapply(hapmap[[alleles_col]], has_more_than_three_alleles), ]

  return(hapmap)
}

#' Read HapMap File
#'
#' Reads a HapMap file from various formats (plain text, gzipped, or zipped).
#'
#' @param filepath A string specifying the path to the HapMap file.
#'
#' @return A data.table containing the HapMap data.
#'
#' @examples
#' # Assuming you have a HapMap file at the specified path
#' hapmap_data <- read_hapmap_file("path/to/hapmap_file.txt.gz")
#'

read_hapmap_file <- function(filepath) {
  msg <- "\n==================================================\n"
  file_extension <- tools::file_ext(filepath)

  if (file_extension == "gz") {
    # Read gzipped file
    hapmap_data <- data.table::fread(filepath, sep='\t', header=TRUE, check.names=FALSE,
                                     skip = "#",
                                     na.strings=c(NA,"N","NN","B","V","H","D",".","-"))
  } else if (file_extension == "zip") {
    # Create a unique temporary directory
    mainDir <- getwd()
    systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
    systime <- gsub("[-: ]", "_", systime)
    subDir <- paste("hapmap_unzip", systime, sep = "_")

    temp_dir <- file.path(mainDir, subDir)
    dir.create(temp_dir)

    # Extract the filename inside the zip (assuming only one file)
    unzip(filepath, exdir = temp_dir)
    unzipped_files <- list.files(temp_dir, full.names = TRUE)

    if (length(unzipped_files) != 1) {
      unlink(temp_dir, recursive = TRUE)  # Remove the temporary directory
      stop(paste(msg,"Zip file should contain exactly one file."), call. = FALSE)
    }

    hapmap_data <- data.table::fread(unzipped_files[1], sep='\t', header=TRUE,
                                     skip = "#",
                                     check.names = FALSE,
                                     na.strings=c(NA,"N","NN","B","V","H","D",".","-"))

    # Remove the temporary directory
    unlink(temp_dir, recursive = TRUE)

  } else if(file_extension == "txt"){

    hapmap_data <- data.table::fread(unzipped_files[1], sep='\t', header=TRUE,
                                     skip = "#",
                                     check.names = FALSE,
                                     na.strings=c(NA,"N","NN","B","V","H","D",".","-"))

  } else {
    # Read plain text file
    stop(paste(msg, "Unknown file extension."), call. = FALSE)
  }

  return(hapmap_data)
}

#' Check HapMap Columns
#'
#' Checks and ensures that the first 11 columns of the HapMap data match the required format.
#'
#' @param data A data frame containing HapMap data.
#'
#' @return Throws an error if the columns do not match; otherwise, returns invisibly.
#'
#' @examples
#' hapmap <- data.frame(rs = 1:3, alleles = c("A/T", "C/G", "A/A"), chrom = 1:3, pos = 1:3, strand = "+", assembly = NA,
#'                      center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA, X1 = c("AA", "CC", "AT"), X2 = c("TA", "GC", "AA"))
#' check_hapmap_columns(hapmap)
#'

check_hapmap_columns <- function(data) {

  msg <- "\n==================================================\n"

  required_columns <- c("rs#", "alleles", "chrom", "pos", "strand",
                        "assembly", "center", "protLSID", "assayLSID",
                        "panelLSID", "QCcode")

  # Convert both the required columns and data column names to lowercase for comparison
  data_columns <- tolower(colnames(data)[1:11])
  required_columns_lower <- tolower(required_columns)

  # Normalize "assembly#" to "assembly"
  data_columns <- gsub("assembly#", "assembly", data_columns)

  # Check if the first 11 columns match the required columns
  if (all(data_columns[1:11] == required_columns_lower)) {
    message(paste(msg,paste("The first 11 columns match the required HapMap format.",
                            "Data succefully loaded for processing.")))
    return(invisible())
  } else {
    # Generate an error message with a guided example
    error_message <- paste(
      "The first 11 columns do not match the required HapMap format.\n",
      "Expected columns (case-insensitive, in order):\n",
      paste(required_columns, collapse = ", "), "\n",
      "Your columns:\n",
      paste(colnames(data)[1:11], collapse = ", "), "\n",
      "Please ensure your data matches the required format."
    )
    stop(paste(msg,error_message), call. = FALSE)
  }
}

#' Quality Control and Recoding for HapMap Data
#'
#' Performs quality control checks and recodes HapMap genotype data. It allows for filtering based on minor allele frequency (MAF), heterozygosity, individual and SNP call rates, and optionally imputes missing data. The function can also recode genotype data into numeric formats.
#'
#' @param hapmap_file_name A string specifying the name of the HapMap file to be processed. If NULL, `hapmap` must be provided.
#' @param hapmap_file_path A string specifying the path to the directory containing the HapMap file. If NULL, `hapmap` must be provided.
#' @param maf_threshold A numeric value specifying the threshold for minor allele frequency (MAF). Markers below this threshold will be removed.
#' @param het_threshold A numeric value specifying the threshold for heterozygosity. Markers above this threshold will be removed.
#' @param ind_call_rate_threshold A numeric value specifying the threshold for individual call rate. Individuals below this threshold will be removed.
#' @param snp_call_rate_threshold A numeric value specifying the threshold for SNP call rate. SNPs below this threshold will be removed.
#' @param impute A logical value indicating whether missing data should be imputed. If TRUE, missing values will be imputed.
#' @param recode_format A string specifying the format for recoding genotype data. Can be "0,1,2" or "-1,0,1".
#' @param out_put_map A logical value indicating whether to output a map of SNP information.
#' @param message A logical value indicating whether to display messages during processing.
#' @param ... Additional arguments to be passed to underlying functions.
#' @importFrom data.table := .SD .SDcols apply lapply
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item{snps_matrix}{A matrix or data.table of the processed and possibly recoded SNP data.}
#'   \item{snp_map}{An optional data.table containing SNP mapping information, depending on `out_put_map`.}
#'   \item{qc_metrics_and_summary_stat}{A data.frame summarizing the quality control metrics and actions taken.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Save the example data to a file
#'   example_hapmap <- data.table::data.table(
#'     `rs#` = c("rs1", "rs2", "rs3", "rs4"),
#'     alleles = c("A/T", "C/G", "A/A", "G/T"),
#'     chrom = c(1, 1, 2, 2),
#'     pos = c(100, 200, 300, 400),
#'     strand = c("+", "+", "-", "-"),
#'     assembly = c("v1", "v1", "v1", "v1"),
#'     center = c("C1", "C1", "C1", "C1"),
#'     protLSID = c("p1", "p1", "p1", "p1"),
#'     assayLSID = c("a1", "a1", "a1", "a1"),
#'     panelLSID = c("panel1", "panel1", "panel1", "panel1"),
#'     QCcode = c(NA, NA, NA, NA),
#'     Sample1 = c("AA", "CC", "AA", "GG"),
#'     Sample2 = c("TT", "GG", "AA", "TT"),
#'     Sample3 = c("AT", "CG", "AA", "GT")
#'   )
#'   data.table::fwrite(example_hapmap, "example_hapmap.txt", sep = "\t")
#'
#'   # Run the function with the example data
#'   result <- hmp_qc_recode(hapmap_file_name = "example_hapmap.txt",
#'                           hapmap_file_path = ".",
#'                           maf_threshold = 0.05,
#'                           het_threshold = 0.1,
#'                           ind_call_rate_threshold = 0.95,
#'                           snp_call_rate_threshold = 0.95,
#'                           impute = TRUE,
#'                           recode_format = "0,1,2",
#'                           out_put_map = TRUE,
#'                           message = TRUE)
#' }
#'
#' @export

hmp_qc_recode <- function(hapmap_file_name = NULL,
                          hapmap_file_path = NULL,
                          maf_threshold = 0.01,
                          het_threshold = 0.1,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          impute = TRUE,
                          recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
                          out_put_map = TRUE,
                          message = TRUE,
                          ...
)
{
  msg <- "\n==================================================\n"

  markers_callrate_removed  <-  0
  ind_callrate_removed  <-  0
  het_markers_removed  <-  0
  maf_markers_removed <-  0
  snp_data <- NULL
  total_mono <-  0

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

  # Define replacement vectors
  heterozygous <- c('R', 'Y', 'S', 'W', 'K', 'M')
  missing_values <- c(NA, "NA", "N", "NN", "B", "V", "H", "D", ".", "-")
  valid_IUPAC <- c('A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N')
  ###############################

  if(!is.null(hapmap_file_name) & is.null(hapmap_file_path)) hapmap_file_path <- getwd()
  if(!is.null(hapmap_file_name) && !is.null(hapmap_file_path)){

    # Construct the full file path
    full_file_path <- file.path(hapmap_file_path, hapmap_file_name)
    hapmap <-  read_hapmap_file(full_file_path)
    # Read the file using fread
    # hapmap <- data.table::fread(full_file_path,
    #                             header = TRUE,
    #                             skip = "#",
    #                             check.names = FALSE,
    #                             na.strings=c(NA,"N","NN","B","V","H","D",".","-"))


  }else{
    if(is.null(hapmap_file_path) & is.null(hapmap_file_name)){
      stop((paste(msg,"Hapmap data is missing.")), call. = FALSE)

    }

  }

  # Check the columns
  check_hapmap_columns(hapmap)

  # if(isFALSE(data.table::is.data.table(hapmap))){
  #
  #   hapmap <- data.table::as.data.table(hapmap)
  # }
  #### Check and to be sure the expected header for snps name is present.
  ## it is usually rs#
  snp_names_check <- grep("rs", colnames(hapmap)[1:11], ignore.case = TRUE)
  if(length(snp_names_check)==0 | length(snp_names_check)>1){
    stop(paste(msg,"rs# header/column name is missing or provided wrongly."), call. = FALSE)
  }
  #########

  first_geno <- as.data.frame(hapmap[, 12])
  # IUPAC single-letter code
  if (!all(first_geno[, 1] %in% valid_IUPAC)) {

    # Construct an informative error message
    error_message_stop <- paste(
      msg,  # Custom message passed into the function
      "Please provide valid IUPAC codes.",
      "Examples of valid codes include:",
      paste(valid_IUPAC, collapse = ", ")
    )

    error_message <- paste(
      msg,  # Custom message passed into the function
      "Valid IUPAC codes was not provided. We fix it for you.")
    hapmap <- IUPAC_hapmap_compatible(hapmap)
    # Stop execution and return the error message
    first_geno <- as.data.frame(hapmap[, 12])
    if (!all(first_geno[, 1] %in% valid_IUPAC)) {

      stop(paste(msg, error_message_stop), call. = FALSE)
    } else {
      message(paste(msg, error_message))
    }

  }

  # Remove monomorphic markers
  monomorphic_markers <- which(apply(hapmap[, 12:ncol(hapmap)], 1, function(x) length(table(x)) <= 1))
  #monomorphic <- hapmap[, apply(.SD, 1, function(x) length(unique(x)) <= 1), .SDcols = 12:ncol(hapmap)]
  #monomorphic_markers <- which(monomorphic)

  #monomorphic_markers <- which(apply(hapmap[, 12:ncol(hapmap)], 1, function(x) length(table(x)) <= 1))
  if (length(monomorphic_markers) > 0) {
    if (isTRUE(message)) {
      message(insight::print_color(paste(msg, paste("Removing monomorphic markers:", length(monomorphic_markers))), "blue"))
    }
    hapmap <- hapmap[-monomorphic_markers, ]
    #map_data <- map_data[-monomorphic_markers, ]
    total_mono <-  length(monomorphic_markers)
    # rm(monomorphic_markers); gc()
  }

  rm(monomorphic_markers)
  #gc()


  ### Remove markers with high missing value based on desired thresold
  if(!is.null(snp_call_rate_threshold)) {
    # Calculate individual call rate
    snp_call_rate <- rowMeans(is.na(hapmap[, 12:ncol(hapmap)]))
    #snp_call_rate <- hapmap[, rowMeans(is.na(.SD)), .SDcols = 12:ncol(hapmap)]

    # Filter based on individual call rate
    #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
    low_call_rate_snps <- which(snp_call_rate > snp_call_rate_threshold)
    if (length(low_call_rate_snps) > 0) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing SNPs with low call rate:", length(low_call_rate_snps))), "blue"))
      }
      hapmap <- hapmap[-low_call_rate_snps, ]
      #hapmap <- hapmap[low_call_rate_snps, ]
      #map_data <- map_data[low_call_rate_snps, ]
      markers_callrate_removed  <-  length(low_call_rate_snps)
    } else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No SNPs removed based on SNP call rate threshold."), "blue"))
      }

      markers_callrate_removed  <-  0
    }
    rm(snp_call_rate, low_call_rate_snps)
    #; gc()
  }

  ##### Remove individual with bad call rate

  ### Remove markers with high missing value based on desired thresold
  if(!is.null(ind_call_rate_threshold)) {
    # Calculate individual call rate
    ind_call_rate <- colMeans(is.na(hapmap[, 12:ncol(hapmap)]))
    #ind_call_rate <- hapmap[, colMeans(is.na(.SD)), .SDcols = 12:ncol(hapmap)]

    # Filter based on SNP call rate
    #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
    low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)

    low_call_rate_inds <-  low_call_rate_inds+11
    if (length(low_call_rate_inds) > 0) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing Individuals with low call rate:", length(low_call_rate_inds))), "blue"))
      }
      hapmap <- hapmap[, -low_call_rate_inds, with = FALSE]
    } else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No individuals removed based on call rate threshold."), "blue"))
      }

      ind_callrate_removed <-  0
    }

    rm(ind_call_rate, low_call_rate_inds)
    #; gc()
  }


  ## Begin process to remove heterozygo
  if(!is.null(het_threshold)) {
    heteroz <- apply(hapmap[, 12:ncol(hapmap)], 1, function(x){
      return(length(which(x%in%heterozygous))/length(x))
    })

    # heteroz <- hapmap[, {
    #   row_heteroz <- apply(.SD, 1, function(x) sum(x %in% c('R', 'Y', 'S', 'W', 'K', 'M')) / length(x))
    #   list(heteroz = row_heteroz)
    # }, .SDcols = 12:ncol(hapmap)]
    #
    # # Convert to vector
    # heteroz <- heteroz$heteroz

    #het_markers <- which(heteroz >= het_threshold)
    het_markers <- which(heteroz > het_threshold)
    if (length(het_markers) > 0) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing markers with high heterozygosity:", length(het_markers))), "blue"))
      }
      hapmap <- hapmap[-het_markers, ]

      het_markers_removed <-  length(het_markers)
      rm(het_markers, heteroz)
      #; gc()
      #map_data <- map_data[-het_markers, ]
    } else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No markers removed based on heterozygosity threshold."), "blue"))
      }

      het_markers_removed <-  0
    }
  }

  #####
  ### Beging process to remove maf
  if(!is.null(maf_threshold)) {

    minor_allele_freq <- apply(hapmap[, c(2,12:ncol(hapmap)), with=FALSE], 1, function(x){
      allele1 <- substr(x[1], 1, 1)
      allele2 <- substr(x[1], 3, 3)
      if(length(which(x==allele1)) > length(which(x==allele2))){
        minor_allele <- allele2
      }
      else{minor_allele <- allele1}
      minor_allele_freq <- ((2*length(which(x[-1]==minor_allele))) + length(which(!x[-1]%in%c(allele1,allele2))))/(2*(length(x[-1])))
      return(minor_allele_freq)
    })

    #hapmap <- hapmap[minor_allele_freq >= maf_threshold,]
    maf_markers <-  which(minor_allele_freq < maf_threshold)

    if (length(maf_markers) > 0) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Removing markers with MAF below threshold:", length(maf_markers))), "blue"))
      }
      hapmap <- hapmap[-maf_markers,]

      maf_markers_removed <- length(maf_markers)
      rm(maf_markers)
      #; gc()

    }else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No markers removed based on MAF threshold."), "blue"))
      }

      maf_markers_removed <-  0
    }

  }

  ### remove markers that is with reference allele greater than 2
  hapmap <- remove_multiallelic_markers(hapmap)


  if (!is.null(recode_format)) {

    hapmap2numeric_meth2_1_0 <- function(hapmap){
      hapmap_numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
        # Replace heterozygous markers with 1
        x[which(x %in% heterozygous)] <- 1

        # Replace specific missing value indicators with NA
        x[which(x %in% missing_values)] <- NA

        allele1 <- substr(x[2], 1, 1)
        allele2 <- substr(x[2], 3, 3)
        if(length(which(x==allele1)) >= length(which(x==allele2))){
          x[which(x==allele1)] <- 2 # major allele is 2
          x[which(x==allele2)] <- 0
        }
        else{
          x[which(x==allele1)] <- 0
          x[which(x==allele2)] <- 2
        }
        return(x)
      })
      return(t(hapmap_numeric[-(1:2), ]))
      #return(hapmap_numeric)
    }

    #sik = unlist(hapmap_numeric)

    hapmap2numeric_meth1_1 <- function(hapmap){

      hapmap_numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
        x[which(x%in% heterozygous)] <- 0
        x[which(x%in% missing_values)] <- NA
        allele1 <- substr(x[2], 1, 1)
        allele2 <- substr(x[2], 3, 3)
        if(length(which(x==allele1)) >= length(which(x==allele2))){
          x[which(x==allele1)] <- 1 # major allele is 1
          x[which(x==allele2)] <- -1
        }
        else{
          x[which(x==allele1)] <- -1
          x[which(x==allele2)] <- 1
        }
        return(x)
      })


      return(t(hapmap_numeric[-(1:2), ]))
      #return(t(hapmap_numeric))
    }

    if (recode_format == "0,1,2"){
      # Here major allele is 2, minor is 0 and hetero is 1

      # Convert Hapmap allele into numeric
      snp_data <- hapmap2numeric_meth2_1_0(hapmap)
      #snp_data <-  cbind(hapmap[, 1:11], snp_data)
    } else if (recode_format == "-1,0,1") {
      snp_data <- hapmap2numeric_meth1_1(hapmap[, 12:ncol(hapmap)])
      #snp_data <-  cbind(hapmap[, 1:11], snp_data)
    } else {
      stop(paste(msg, "Invalid recode format. Use '0,1,2' or '-1,0,1'."), call. = FALSE)

    }

    #rm(hapmap); gc()

    if(isTRUE(impute)){
      #if(any(is.na(snp_data[, 12:ncol(snp_data)]))==T) {
      if(any(is.na(snp_data[, 1:ncol(snp_data)]))==T) {
        #This function will impute the missing values
        #for(j in 12:ncol(snp_data)){
        for(j in 1:ncol(snp_data)){
          #tmp <- snp_data[,j, with=FALSE]
          tmp <- snp_data[,j]
          tmp = as.double(as.character(unlist(tmp)))
          snp_data[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
        }

      }

    }
    # End
  }
  #original_row_names <- colnames(snp_data)[12:ncol(snp_data)]
  original_row_names <- colnames(snp_data)
  ##################
  ## Aggregate all the info
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

  if(isTRUE(out_put_map) & !is.null(snp_data)){
    #map <-   snp_data[, 1:11]
    map <-  hapmap[, 1:11]
    snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
    #snp_data <- t(snp_data[, 12:ncol(snp_data)])
    snp_data <- t(snp_data[, 1:ncol(snp_data)])
    if (!is.numeric(snp_data)) {
      # Apply as.numeric to each element in the matrix
      snp_data <- apply(snp_data, c(1, 2), as.numeric)
    }
    colnames(snp_data) <- as.data.frame(map)[, snp_names_index]
    #snp_data <- apply(snp_data, 2, as.double)

    rownames(snp_data) <- original_row_names

    class(snp_data) <- c("matrix", "array", "genotype")

    return(list(snps_matrix= snp_data,
                snp_map = map,
                qc_metrics_and_summary_stat = summary_stat_snp
    )
    )


  } else if(isFALSE(out_put_map) & !is.null(snp_data)){
    #map <-   snp_data[, 1:11]
    map <-  hapmap[, 1:11]
    snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
    #snp_data <- t(snp_data[, 12:ncol(snp_data)])
    snp_data <- t(snp_data[, 1:ncol(snp_data)])
    colnames(snp_data) <- as.data.frame(map)[, snp_names_index]

    #snp_data <- apply(snp_data, 2, as.double)

    rownames(snp_data) <- original_row_names

    class(snp_data) <- c("matrix", "array", "genotype")
    rm(map); gc()

    return(list(snps_matrix= snp_data,
                qc_metrics_and_summary_stat = summary_stat_snp))

  } else if(isTRUE(out_put_map) & is.null(recode_format)){
    map <-   hapmap[, 1:11]
    snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
    hapmap <- t(hapmap[, 12:ncol(hapmap)])
    colnames(hapmap) <- as.data.frame(map)[, snp_names_index]

    return(list(snps_matrix= hapmap,
                snp_map = map,
                qc_metrics_and_summary_stat = summary_stat_snp))

  } else {
    if(isFALSE(out_put_map) & is.null(recode_format)){
      map <-   hapmap[, 1:11]
      snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
      hapmap <- t(hapmap[, 12:ncol(hapmap)])
      colnames(hapmap) <- as.data.frame(map)[, snp_names_index]
      rm(map); gc()

      return(list(snps_matrix= hapmap,
                  qc_metrics_and_summary_stat = summary_stat_snp))

    }

  }

}


# convert_to_double_letter <- function(df, cols) {
#   # Define the mapping for single nucleotides and IUPAC codes
#   double_letter_map <- c(
#     "A" = "A:A", "C" = "C:C", "G" = "G:G", "T" = "T:T",
#     "M" = "A:C", "R" = "A:G", "W" = "A:T",
#     "S" = "C:G", "Y" = "C:T", "K" = "G:T"
#   )
#
#   # Apply the mapping to the specified columns
#   df[, cols] <- lapply(df[, cols], function(column) {
#     sapply(column, function(value) double_letter_map[value])
#   })
#   return(df)
# }

# # Define the function to convert single nucleotide values to double letters
# convert_to_double_letter <- function(df, cols) {
#   # Define the mapping for single nucleotides and IUPAC codes
#   double_letter_map <- c(
#     "A" = "AA", "C" = "CC", "G" = "GG", "T" = "TT",
#     "M" = "AC", "R" = "AG", "W" = "AT",
#     "S" = "CG", "Y" = "CT", "K" = "GT"
#   )
#
#   # Apply the mapping to the specified columns
#   df[, cols] <- lapply(df[, cols], function(column) {
#     sapply(column, function(value) double_letter_map[value])
#   })
#   return(df)
# }
#
# hapmap2 = as.data.frame(hapmap)
# cols=12:ncol(hapmap)
# #class(hapmap)
# # Apply the function to the specified columns
# hapmap_converted <- convert_to_double_letter(df=hapmap2, cols=cols)
#
# hapmap = hapmap_converted
