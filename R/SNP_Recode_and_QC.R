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


#' Title
#'
#' @param hapmap_file_name
#' @param hapmap_file_path
#' @param hapmap
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param hwe_threshold
#' @param impute
#' @param recode_format
#' @param out_put_map
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom data.table := .SD .SDcols lapply
#'
hmp_qc_recode <- function(hapmap_file_name = NULL,
                          hapmap_file_path = NULL,
                          hapmap = NULL,
                          maf_threshold = 0.01,
                          het_threshold = 0.1,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          #hwe_threshold = 0.001,  # Adjust as needed
                          impute = FALSE,
                          recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
                          out_put_map = FALSE,
                          message = TRUE,
                          ...
)
{
  ### Remove All loci with All NAs and monomorphic markers
  ####### TO DO
  # Develop function to accomodate when the data is already recoded but
  ## Not filtered
  #if(!is.null(hapmap) && !is.null(geno))#{

  ###############################

  if(!is.null(hapmap_file_name) && !is.null(hapmap_file_path)){
    # Construct the full file path
    full_file_path <- file.path(hapmap_file_path, hapmap_file_name)
    # Read the file using fread
    hapmap <- data.table::fread(full_file_path,
                                header = TRUE,
                                skip = "#",
                                check.names = FALSE,
                                na.strings=c(NA,"N","NN","B","V","H","D",".","-"))


  } else if (!is.null(hapmap)){

    for (j in 12:ncol(hapmap)) {
      hapmap[[j]] <- lapply(hapmap[[j]], function(x) {
        x[which(x %in% c(NA, "NA", "N", "NN", "B", "V", "H", "D", ".", "-"))] <- NA
        return(x)
      })
    }

  }else{
    if(is.null(hapmap)){

      stop("Hapmap data is missing.")

    }

  }


  if(isFALSE(data.table::is.data.table(hapmap))){

    hapmap <- data.table::as.data.table(hapmap)
  }
  # snp_data = hapmap[, 12:ncol(hapmap)]
  # rownames(snp_data) <- hapmap$`rs#`
  #
  # map = hapmap[, 1:11]

  #########

  # Remove monomorphic markers
  monomorphic_markers <- which(apply(hapmap[, 12:ncol(hapmap)], 1, function(x) length(table(x)) <= 1))
  if (length(monomorphic_markers) > 0) {
    if (message) {
      print(paste("Removing monomorphic markers:", length(monomorphic_markers)))
    }
    hapmap <- hapmap[-monomorphic_markers, ]
    #map_data <- map_data[-monomorphic_markers, ]
    total_mono = length(monomorphic_markers)
    rm(monomorphic_markers); gc()
  } else {
    if (message) {
      print("No monomorphic markers to remove.")
    }
    total_mono = 0
  }



  ### Remove markers with high missing value based on desired thresold
  if(!is.null(snp_call_rate_threshold)) {

    # Calculate individual call rate
    snp_call_rate <- rowMeans(is.na(hapmap[, 12:ncol(hapmap)]))

    # Filter based on individual call rate
    #low_call_rate_snps <- which(snp_call_rate < snp_call_rate_threshold)
    low_call_rate_snps <- which(snp_call_rate > snp_call_rate_threshold)
    if (length(low_call_rate_snps) > 0) {
      if (message) {
        print(paste("Removing SNPs with low call rate:", length(low_call_rate_snps)))
      }
      hapmap <- hapmap[-low_call_rate_snps, ]
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

  ##### Remove individual with bad call rate

  ### Remove markers with high missing value based on desired thresold
  if(!is.null(ind_call_rate_threshold)) {
    # Calculate individual call rate
    ind_call_rate <- colMeans(is.na(hapmap[, 12:ncol(hapmap)]))

    # Filter based on SNP call rate
    #low_call_rate_inds <- which(ind_call_rate < ind_call_rate_threshold)
    low_call_rate_inds <- which(ind_call_rate > ind_call_rate_threshold)

    low_call_rate_inds <-  low_call_rate_inds+11
    if (length(low_call_rate_inds) > 0) {
      if (message) {
        print(paste("Removing Individuals with low call rate:", length(low_call_rate_inds)))
      }
      hapmap <- hapmap[, -low_call_rate_inds, with = FALSE]
    } else {
      if (message) {
        print("No individuals removed based on call rate threshold.")
      }

      ind_callrate_removed = 0
    }

    rm(ind_call_rate, low_call_rate_inds); gc()
  }


  ## Begin process to remove heterozygo
  if(!is.null(het_threshold)) {
    heteroz <- apply(hapmap[, 12:ncol(hapmap)], 1, function(x){
      return(length(which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M')))/length(x))
    })

        #het_markers <- which(heteroz >= het_threshold)
    het_markers <- which(heteroz > het_threshold)
    if (length(het_markers) > 0) {
      if (message) {
        print(paste("Removing markers with high heterozygosity:", length(het_markers)))
      }
      hapmap <- hapmap[-het_markers, ]

      het_markers_removed = length(het_markers)
      rm(het_markers, heteroz); gc()
      #map_data <- map_data[-het_markers, ]
    } else {
      if (message) {
        print("No markers removed based on heterozygosity threshold.")
      }

      het_markers_removed = 0
    }

    #snp_data <- snp_data[heteroz < het_threshold, ]

    ## Gather meta data


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
      if (message) {
        print(paste("Removing markers with MAF below threshold:", length(maf_markers)))
      }
    hapmap <- hapmap[-maf_markers,]

    maf_markers_removed <- length(maf_markers)
  rm(maf_markers); gc()

    }else {
      if (message) {
        print("No markers removed based on MAF threshold.")
      }

      maf_markers_removed = 0
    }
    #maf_markers_removed <-  (N_snp - ncol(hapmap))

  }

  if (!is.null(recode_format)) {
  #if (recode_format %in% c("0,1,2", "-1,0,1")) {

    hapmap2numeric_meth2_1_0 <- function(hapmap){
      hapmap_numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
        x[which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M'))] <- 1
        x[which(x%in%c(NA,"NA","N","NN","B","V","H","D",".","-"))] <- NA
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
      #return(t(hapmap_numeric))
    }

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


      return(t(hapmap_numeric[-(1:2), ]))
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

    rm(hapmap); gc()

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

  if(isTRUE(out_put_map) & exists("snp_data")){
   map <-   snp_data[, 1:11]

   snp_data <- t(snp_data[, 12:ncol(snp_data)])
   colnames(snp_data) <- map$`rs#`

   snp_data <- apply(snp_data, 2, as.double)

   class(snp_data) <- c("matrix", "array", "genotype")

   return(list(snps_matrix= snp_data,
              snp_map = map,
              qc_metrics_and_summary_stat = summary_stat_snp
   )
   )


  } else if(isFALSE(out_put_map) & exists("snp_data")){
    map <-   snp_data[, 1:11]

    snp_data <- t(snp_data[, 12:ncol(snp_data)])
    colnames(snp_data) <- map$`rs#`

    snp_data <- apply(snp_data, 2, as.double)

    class(snp_data) <- c("matrix", "array", "genotype")
    rm(map); gc()

    return(list(snps_matrix= snp_data,
                qc_metrics_and_summary_stat = summary_stat_snp
    )
    )

  } else if(isTRUE(out_put_map) & is.null(recode_format)){
    map <-   hapmap[, 1:11]
    hapmap <- t(hapmap[, 12:ncol(hapmap)])
    colnames(hapmap) <- map$`rs#`

    return(list(snps_matrix= hapmap,
                snp_map = map,
                qc_metrics_and_summary_stat = summary_stat_snp
    )
    )

  } else {
    if(isFALSE(out_put_map) & is.null(recode_format)){
      map <-   hapmap[, 1:11]
      hapmap <- t(hapmap[, 12:ncol(hapmap)])
      colnames(hapmap) <- map$`rs#`
      rm(map); gc()

      return(list(snps_matrix= hapmap,
                  qc_metrics_and_summary_stat = summary_stat_snp
      )
      )

    }

  }









  #return(filter)


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
