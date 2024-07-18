
# Function to read data from various file types
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

# Function to check and ensure the first 11 columns
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
#' @param hapmap An optional data.table object containing HapMap data. If NULL, `hapmap_file_name` and `hapmap_file_path` must be provided.
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
#'   result <- hmp_qc_recode(hapmap_file_name = "sample_data.txt",
#'                           hapmap_file_path = "path/to/data",
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


#hapmap_file_path = "D:/PEA_BARI"

#hapmap_file_name = 'HapMap_300_Aug25_Select_Chr_No_Filtering.hmp.txt.zip'

hmp_qc_recode <- function(hapmap_file_name = NULL,
                          hapmap_file_path = NULL,
                          hapmap = NULL,
                          maf_threshold = 0.01,
                          het_threshold = 0.1,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          #hwe_threshold = 0.001,  # Adjust as needed
                          impute = TRUE,
                          recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
                          out_put_map = TRUE,
                          message = TRUE,
                          ...
)
{
  msg <- "\n==================================================\n"
  ### Remove All loci with All NAs and monomorphic markers
  ####### TO DO
  # Develop function to accomodate when the data is already recoded but
  ## Not filtered
  #if(!is.null(hapmap) && !is.null(geno))#{
## Place holders
  markers_callrate_removed  <-  0
  ind_callrate_removed  <-  0
  het_markers_removed  <-  0
  maf_markers_removed <-  0
  snp_data <- NULL

  # Define replacement vectors
  heterozygous <- c('R', 'Y', 'S', 'W', 'K', 'M')
  missing_values <- c(NA, "NA", "N", "NN", "B", "V", "H", "D", ".", "-")

  ###############################

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


  } else if (!is.null(hapmap)){

    if(!inherits(hapmap, "data.table")) {

      hapmap <- data.table::as.data.table(hapmap)

    }


  }else{
    if(is.null(hapmap)){
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
    rm(monomorphic_markers); gc()
  } else {
    if (isTRUE(message)) {
      message(insight::print_color(paste(msg, "No monomorphic markers to remove."), "blue"))
    }
    total_mono <-  0
  }



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
    rm(snp_call_rate, low_call_rate_snps); gc()
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

    rm(ind_call_rate, low_call_rate_inds); gc()
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
      rm(het_markers, heteroz); gc()
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
  rm(maf_markers); gc()

    }else {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, "No markers removed based on MAF threshold."), "blue"))
      }

      maf_markers_removed <-  0
    }

  }


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
        x[which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M'))] <- 0
        x[which(x%in%c(NA,"NA","N","NN","B","V","H","D",".","-"))] <- NA
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


      #return(t(hapmap_numeric[-(1:2), ]))
      #return(t(hapmap_numeric))
    }

    if (recode_format == "0,1,2"){
      # Here major allele is 2, minor is 0 and hetero is 1

      # Convert Hapmap allele into numeric
      snp_data <- hapmap2numeric_meth2_1_0(hapmap)
      snp_data <-  cbind(hapmap[, 1:11], snp_data)
      } else if (recode_format == "-1,0,1") {
        snp_data <- hapmap2numeric_meth1_1(hapmap[, 12:ncol(hapmap)])
        snp_data <-  cbind(hapmap[, 1:11], snp_data)
    } else {
      stop("Invalid recode format. Use '0,1,2' or '-1,0,1'.")

    }

    #rm(hapmap); gc()

    if(isTRUE(impute)){
      if(any(is.na(snp_data[, 12:ncol(snp_data)]))==T) {
        #This function will impute the missing values
        for(j in 12:ncol(snp_data)){
          tmp <- snp_data[,j, with=FALSE]
          tmp = as.double(as.character(unlist(tmp)))
          snp_data[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
        }

      }

    }
    # End
  }
  original_row_names <- colnames(snp_data)[12:ncol(snp_data)]
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
   map <-   snp_data[, 1:11]
   snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
   snp_data <- t(snp_data[, 12:ncol(snp_data)])
   colnames(snp_data) <- as.data.frame(map)[, snp_names_index]
   snp_data <- apply(snp_data, 2, as.double)
   rownames(snp_data) <- original_row_names

   class(snp_data) <- c("matrix", "array", "genotype")

   return(list(snps_matrix= snp_data,
              snp_map = map,
              qc_metrics_and_summary_stat = summary_stat_snp
   )
   )


  } else if(isFALSE(out_put_map) & !is.null(snp_data)){
    map <-   snp_data[, 1:11]
    snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
    snp_data <- t(snp_data[, 12:ncol(snp_data)])
    colnames(snp_data) <- as.data.frame(map)[, snp_names_index]

    snp_data <- apply(snp_data, 2, as.double)

    rownames(snp_data) <- original_row_names

    class(snp_data) <- c("matrix", "array", "genotype")
    rm(map); gc()

    return(list(snps_matrix= snp_data,
                qc_metrics_and_summary_stat = summary_stat_snp
    )
    )

  } else if(isTRUE(out_put_map) & is.null(recode_format)){
    map <-   hapmap[, 1:11]
    snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
    hapmap <- t(hapmap[, 12:ncol(hapmap)])
    colnames(hapmap) <- as.data.frame(map)[, snp_names_index]

    return(list(snps_matrix= hapmap,
                snp_map = map,
                qc_metrics_and_summary_stat = summary_stat_snp
    )
    )

  } else {
    if(isFALSE(out_put_map) & is.null(recode_format)){
      map <-   hapmap[, 1:11]
      snp_names_index <- grep("rs", colnames(map), ignore.case = TRUE)
      hapmap <- t(hapmap[, 12:ncol(hapmap)])
      colnames(hapmap) <- as.data.frame(map)[, snp_names_index]
      rm(map); gc()

      return(list(snps_matrix= hapmap,
                  qc_metrics_and_summary_stat = summary_stat_snp
      )
      )

    }

  }

}


## -*- mode: R-4.1.2 -*-
###########################################################################
##
## hmp to numeric (AA=1,Aa=0,aa=-1) | (AA=2,Aa=1,aa=0)
##
##              Copyright (C) 2023 Sikiru
##
## ** Filename: snp_recode_and_qc.R
##
## ** This function perform QC and recoding of SNPs
##
## * Function Inputs :
## * hapmap = hapmap file
## * het = percentage of heterozyoug SNPs allowed (numberic)
## * maf = minor allele frequency (numeric)
## * call_rate = percentage of missing value permitted (numeric)
## * recode = recoding format for the SNP which can be either Meth2_1_0 that is (AA=2,Aa=1,aa=0) or "Meth1_0_-1 that is (AA=1,Aa=0,aa=-1)
##
##
## ** Authors: Sikiru
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but
##   WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##   General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software Foundation,
##   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##
#############################################################################
#############################################################################
## This function required installation of data.table library

# library(data.table)
#
# # CLEAR 'WorkSpace' (R  environment)
#
# rm(list = ls()); ls()
# gc()
# cat("\014")
#
#
# ## Set working directory where you have the hapmap file
# setwd("D:/PEA_BARI")
#
#
# ## Load Hapmap file
# Selec_ch.300 <- fread('HapMap_300_Aug25_Select_Chr_No_Filtering.hmp.txt', sep='\t', header=T, check.names=FALSE)

# Selec_ch.300 <- read.delim("HapMap_300_Aug25_Select_Chr_No_Filtering.hmp.txt",
#                            colClasses = "character",
#                            #nrows=50,
#                            comment.char="", check.names=FALSE,
#                            header=TRUE,na.strings=c(NA,"N","NN","B","V","H","D",".","-"))
#
# ### Data must be in data.table.
#
# Selec_ch.300 <-  as.data.table(Selec_ch.300)
#

## Start of the functions for QC and Recording of SNP



# hmp_qc_recode <- function(hapmap_file_name = "HapMap_300_Aug25_Select_Chr_No_Filtering.hmp.txt",
#                           hapmap_file_path = "D:/PEA_BARI",
#                           hapmap = NULL,
#                           maf_threshold = 0.01,
#                           het_threshold = 0.1,
#                           ind_call_rate_threshold = 0.9,
#                           snp_call_rate_threshold = 0.9,
#                           hwe_threshold = 0.001,  # Adjust as needed
#                           impute = FALSE,
#                           recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
#                           out_put_map = FALSE,
#                           message = TRUE,
#                           ...
# )
# {
#   ### Remove All loci with All NAs and monomorphic markers
#   ####### TO DO
#   # Develop function to accomodate when the data is already recoded but
#   ## Not filtered
#   #if(!is.null(hapmap) && !is.null(geno))#{
#
#   ###############################
#
#   if(!is.null(hapmap_file_name) && !is.null(hapmap_file_path)){
#     # Construct the full file path
#     full_file_path <- file.path(hapmap_file_path, hapmap_file_name)
#     # Read the file using fread
#     hapmap <- data.table::fread(full_file_path,
#                                 header = TRUE,
#                                 check.names = FALSE,
#                                 na.strings=c(NA,"N","NN","B","V","H","D",".","-"))
#
#
#   } else if (!is.null(hapmap)){
#
#     hapmap[, (12:ncol(hapmap)) := lapply(data.datable::.SD, function(x) {
#       x[which(x %in% c(NA, "NA", "N", "NN", "B", "V", "H", "D", ".", "-"))] <- NA
#       return(x)
#     }), .SDcols = 12:ncol(hapmap)]
#
#
#
#   }else{
#     if(is.null(hapmap)){
#
#       stop("Hapmap data is missing.")
#
#     }
#
#   }
#
#
#   if(isFALSE(data.table::is.data.table(hapmap))){
#
#     hapmap <- data.table::as.data.table(hapmap)
#   }
#   # snp_data = hapmap[, 12:ncol(hapmap)]
#   # rownames(snp_data) <- hapmap$`rs#`
#   #
#   # map = hapmap[, 1:11]
#
#   #########
#
#   # Remove monomorphic markers
#   monomorphic_markers <- which(apply(hapmap[, 12:ncol(hapmap)], 1, function(x) length(table(x)) <= 1))
#   if (length(monomorphic_markers) > 0) {
#     if (message) {
#       print(paste("Removing monomorphic markers:", length(monomorphic_markers)))
#     }
#     hapmap <- hapmap[-monomorphic_markers, ]
#     #map_data <- map_data[-monomorphic_markers, ]
#   } else {
#     if (message) {
#       print("No monomorphic markers to remove.")
#     }
#   }
#
#   total_mono = length(monomorphic_markers)
#   rm(monomorphic_markers); gc()
#
#   ### Remove markers with high missing value based on desired thresold
#   if(!is.null(snp_call_rate_threshold)) {
#
#     # Calculate individual call rate
#     snp_call_rate <- rowMeans(is.na(hapmap[, 12:ncol(hapmap)]))
#
#     # Filter based on individual call rate
#     #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
#     low_call_rate_snps <- which(snp_call_rate > snp_call_rate_threshold)
#     if (length(low_call_rate_snps) > 0) {
#       if (message) {
#         print(paste("Removing SNPs with low call rate:", length(low_call_rate_snps)))
#       }
#       hapmap <- hapmap[-low_call_rate_snps, ]
#       #hapmap <- hapmap[low_call_rate_snps, ]
#       #map_data <- map_data[low_call_rate_snps, ]
#     } else {
#       if (message) {
#         print("No SNPs removed based on SNP call rate threshold")
#       }
#     }
#     rm(snp_call_rate, low_call_rate_snps)
#   }
#
#   ##### Remove individual with bad call rate
#
#   ### Remove markers with high missing value based on desired thresold
#   if(!is.null(ind_call_rate_threshold)) {
#     # Calculate individual call rate
#     ind_call_rate <- colMeans(is.na(hapmap[, 12:ncol(hapmap)]))
#
#     # Filter based on SNP call rate
#     #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
#     low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)
#
#     low_call_rate_inds = low_call_rate_inds+11
#     if (length(low_call_rate_inds) > 0) {
#       if (message) {
#         print(paste("Removing Individuals with low call rate:", length(low_call_rate_inds)))
#       }
#       hapmap <- hapmap[, -low_call_rate_inds, with = FALSE]
#     } else {
#       if (message) {
#         print("No individuals removed based on call rate threshold.")
#       }
#     }
#
#     rm(ind_call_rate, low_call_rate_inds)
#   }
#
#
#   ## Begin process to remove heterozygo
#   if(!is.null(het_threshold)) {
#     heteroz <- apply(hapmap[, 12:ncol(hapmap)], 1, function(x){
#       return(length(which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M')))/length(x))
#     })
#
#         #het_markers <- which(heteroz >= het_threshold)
#     het_markers <- which(heteroz > het_threshold)
#     if (length(het_markers) > 0) {
#       if (message) {
#         print(paste("Removing markers with high heterozygosity:", length(het_markers)))
#       }
#       hapmap <- hapmap[-het_markers, ]
#       #map_data <- map_data[-het_markers, ]
#     } else {
#       if (message) {
#         print("No markers removed based on heterozygosity threshold.")
#       }
#     }
#
#     #snp_data <- snp_data[heteroz < het_threshold, ]
#
#     ## Gather meta data
#     het_markers_removed = length(het_markers)
#     rm(het_markers, het_markers); gc()
#
#   }
#
#   #####
#   ### Beging process to remove maf
#   if(!is.null(maf_threshold)) {
#
#     minor_allele_freq <- apply(hapmap[, c(2,12:ncol(hapmap)), with=FALSE], 1, function(x){
#       allele1 <- substr(x[1], 1, 1)
#       allele2 <- substr(x[1], 3, 3)
#       if(length(which(x==allele1)) > length(which(x==allele2))){
#         minor_allele <- allele2
#       }
#       else{minor_allele <- allele1}
#       minor_allele_freq <- ((2*length(which(x[-1]==minor_allele))) + length(which(!x[-1]%in%c(allele1,allele2))))/(2*(length(x[-1])))
#       return(minor_allele_freq)
#     })
#
#
#     #hapmap <- hapmap[minor_allele_freq >= maf_threshold,]
#     maf_markers <-  which(minor_allele_freq < maf_threshold)
#
#     hapmap <- hapmap[-maf_markers,]
#
#     maf_markers_removed <- length(maf_markers)
#  rm(maf_markers); gc()
#     #maf_markers_removed <-  (N_snp - ncol(hapmap))
#
#   }
#
#   if (recode_format %in% c("0,1,2", "-1,0,1")) {
#
#     hapmap2numeric_meth2_1_0 <- function(hapmap){
#       hapmap_numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
#         x[which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M'))] <- 1
#         x[which(x%in%c(NA,"NA","N","NN","B","V","H","D",".","-"))] <- NA
#         allele1 <- substr(x[2], 1, 1)
#         allele2 <- substr(x[2], 3, 3)
#         if(length(which(x==allele1)) >= length(which(x==allele2))){
#           x[which(x==allele1)] <- 2 # major allele is 2
#           x[which(x==allele2)] <- 0
#         }
#         else{
#           x[which(x==allele1)] <- 0
#           x[which(x==allele2)] <- 2
#         }
#         return(x)
#       })
#
#       return(t(hapmap_numeric[-(1:2), ]))
#       #return(t(hapmap_numeric))
#     }
#
#     hapmap2numeric_meth1_1 <- function(hapmap){
#       hapmap_numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
#         x[which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M'))] <- 0
#         x[which(x%in%c(NA,"NA","N","NN","B","V","H","D",".","-"))] <- NA
#         allele1 <- substr(x[2], 1, 1)
#         allele2 <- substr(x[2], 3, 3)
#         if(length(which(x==allele1)) >= length(which(x==allele2))){
#           x[which(x==allele1)] <- 1 # major allele is 1
#           x[which(x==allele2)] <- -1
#         }
#         else{
#           x[which(x==allele1)] <- -1
#           x[which(x==allele2)] <- 1
#         }
#         return(x)
#       })
#
#
#       return(t(hapmap_numeric[-(1:2), ]))
#       #return(t(hapmap_numeric))
#     }
#
#     if (recode_format == "0,1,2"){
#       # Here major allele is 2, minor is 0 and hetero is 1
#
#       # Convert Hapmap allele into numeric
#       snp_data <- hapmap2numeric_meth2_1_0(hapmap)
#       snp_data <-  cbind(hapmap[, 1:11], snp_data)
#       } else if (recode_format == "-1,0,1") {
#         snp_data <- hapmap2numeric_meth1_1(hapmap[, 12:ncol(hapmap)])
#         snp_data <-  cbind(hapmap[, 1:11], snp_data)
#     } else {
#       stop("Invalid recode format. Use '0,1,2' or '-1,0,1'.")
#
#     }
#
#     rm(hapmap); gc()
#
#     if(isTRUE(impute)){
#       if(any(is.na(snp_data[, 12:ncol(snp_data)]))==T) {
#         #This function will impute the missing values
#         for(j in 12:ncol(snp_data)){
#           tmp <- snp_data[,j, with=FALSE]
#           tmp = as.double(as.character(unlist(tmp)))
#           snp_data[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
#         }
#
#       }
#
#     }
#     # End
#
#
#   }
#   ##################
#   ## Aggregate all the info
#   #### Aggregate all the maker data
#
#   if(!exists("markers_callrate_removed")) {
#     markers_callrate_removed  = 0
#   } else {
#     markers_callrate_removed  = markers_callrate_removed
#
#   }
#   if(!exists("ind_callrate_removed ")) {
#     ind_callrate_removed  = 0
#   } else {
#     ind_callrate_removed = ind_callrate_removed
#
#   }
#
#
#   if(!exists("het_markers_removed ")) {
#     het_markers_removed  = 0
#   } else {
#     het_markers_removed  = het_markers_removed
#
#   }
#
#   if(!exists("maf_markers_removed")) {
#     maf_markers_removed = 0
#   } else {
#     maf_markers_removed = maf_markers_removed
#   }
#   #########################
#
#
#   if(isTRUE(out_put_map) & exists("snp_data")){
#    map <-   snp_data[, 1:11]
#
#    snp_data <- t(snp_data[, 12:ncol(snp_data)])
#    colnames(snp_data) <- map$`rs#`
#
#    snp_data <- apply(snp_data, 2, as.double)
#
#    class(snp_data) <- c("matrix", "array", "genotype")
#
#    return(list(snps_matrix= snp_data,
#               snp_map = map,
#               snp_call_rate_threshold = snp_call_rate_threshold,
#               total_snp_removed = markers_callrate_removed,
#               ind_call_rate_threshold = ind_call_rate_threshold,
#               total_genotypes_removed = ind_callrate_removed,
#               het_threshold = het_threshold,
#               total_het_snps_removed = het_markers_removed,
#               maf_threshold = maf_threshold,
#               total_maf_snps_remove = maf_markers_removed,
#               total_monomorphic_snps_removed = total_mono
#    )
#    )
#
#
#   } else if(isFALSE(out_put_map) & exists("snp_data")){
#     map <-   snp_data[, 1:11]
#
#     snp_data <- t(snp_data[, 12:ncol(snp_data)])
#     colnames(snp_data) <- map$`rs#`
#
#     snp_data <- apply(snp_data, 2, as.double)
#
#     class(snp_data) <- c("matrix", "array", "genotype")
#     rm(map); gc()
#
#     return(list(snps_matrix= snp_data,
#                 snp_call_rate_threshold = snp_call_rate_threshold,
#                 total_snp_removed = markers_callrate_removed,
#                 ind_call_rate_threshold = ind_call_rate_threshold,
#                 total_genotypes_removed = ind_callrate_removed,
#                 het_threshold = het_threshold,
#                 total_het_snps_removed = het_markers_removed,
#                 maf_threshold = maf_threshold,
#                 total_maf_snps_remove = maf_markers_removed,
#                 total_monomorphic_snps_removed = total_mono
#     )
#     )
#
#   } else if(isTRUE(out_put_map) & is.null(recode_format)){
#     map <-   hapmap[, 1:11]
#     hapmap <- t(hapmap[, 12:ncol(hapmap)])
#     colnames(hapmap) <- map$`rs#`
#
#     return(list(snps_matrix= hapmap,
#                 snp_map = map,
#                 snp_call_rate_threshold = snp_call_rate_threshold,
#                 total_snp_removed = markers_callrate_removed,
#                 ind_call_rate_threshold = ind_call_rate_threshold,
#                 total_genotypes_removed = ind_callrate_removed,
#                 het_threshold = het_threshold,
#                 total_het_snps_removed = het_markers_removed,
#                 maf_threshold = maf_threshold,
#                 total_maf_snps_remove = maf_markers_removed,
#                 total_monomorphic_snps_removed = total_mono
#     )
#     )
#
#   } else {
#     if(isFALSE(out_put_map) & is.null(recode_format)){
#       map <-   hapmap[, 1:11]
#       hapmap <- t(hapmap[, 12:ncol(hapmap)])
#       colnames(hapmap) <- map$`rs#`
#       rm(map); gc()
#
#       return(list(snps_matrix= hapmap,
#                   snp_call_rate_threshold = snp_call_rate_threshold,
#                   total_snp_removed = markers_callrate_removed,
#                   ind_call_rate_threshold = ind_call_rate_threshold,
#                   total_genotypes_removed = ind_callrate_removed,
#                   het_threshold = het_threshold,
#                   total_het_snps_removed = het_markers_removed,
#                   maf_threshold = maf_threshold,
#                   total_maf_snps_remove = maf_markers_removed,
#                   total_monomorphic_snps_removed = total_mono
#       )
#       )
#
#     }
#
#   }
#
#
#
#
#
#
#
#
#
#   #return(filter)
#
#
# }
#
#
