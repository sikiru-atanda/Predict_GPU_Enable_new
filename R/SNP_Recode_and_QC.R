## -*- mode: R-4.1.2 -*-
###########################################################################
##
## hmp to numeric (AA=1,Aa=0,aa=-1) | (AA=2,Aa=1,aa=0)
##
##              Copyright (C) 2023 Sikiru
##
## ** Filename: SNP_QC_recodePipelineFinalize_GUD.R
##
## ** This function perform QC and recoding of SNPs
##
## * Function Inputs :
## * hapmap = hapmap file
## * het = percentage of heterozyoug SNPs allowed (numberic)
## * maf = minor allele frequency (numeric)
## * call_rate = percentage of missing value permitted (numeric)
## * recode = recoding format for the SNP which can be either Meth2_0 that is (AA=2,Aa=1,aa=0) or "Meth-1_1 that is (AA=1,Aa=0,aa=-1)
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
#' @param hapmap
#' @param het
#' @param maf
#' @param call_rate
#' @param ind_rate
#' @param recode
#' @param impute
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
marker_qc_recode <- function(hapmap = NULL,
                             het=NULL,
                             maf=NULL,
                             call_rate=NULL,
                             ind_rate = NULL,
                             recode = c("Meth2_0",
                                        "Meth-1_1"),
                             impute = NULL,
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
  Xa = hapmap[, 12:ncol(hapmap)]

  Xa <- t(Xa)

  colnames(Xa) <- hapmap$`rs#`
  #########
  # SNP data
  loc.list <- array(NA,ncol(Xa))
  nL <- ncol(Xa)
  n_ind <- nrow(Xa)
  lN <- 1:ncol(Xa)
  #l.names <- colnames(Xa)
  for (i in 1:nL){

    if (length(table(Xa[, i]))==1 |length(table(Xa[, i]))==0){
      loc.list[i] <- lN[i]

    }
  }


  # Remove NAs from list of monomorphic loci and loci with all NAs
  loc.list <- loc.list[!is.na(loc.list)]
  length(loc.list)

  # remove monomorphic loc and loci with all NAs

  if(length(loc.list > 0)){
    if(message){
    print(paste("Removing monomorphic loci", length(loc.list)))
    }
    # Remove loci flagged for deletion
    Xa <- Xa[,-loc.list]
    #SNP_base <- SNP_base[,!colnames(SNP_base)%in%loc.list]
  } else {
    if(message){
    print("No monomorphic loci to remove")
    }
  }

  ## Convert it back to SNP in the row and genotypes in the column for ease in the
  # later stage analysis
  #Xa <- t(Xa)

  #rm(loc.list,  nL, lN)

  rm(loc.list, lN)

  ### Remove markers with high missing value based on desired thresold
  if(!is.null(call_rate)) {


    CR <- colMeans(is.na(Xa))
    #poscr <- CR <= call_rate #selected by CR

    poscr <-  which(CR<call_rate)

    Xa <- Xa[, poscr]
    ## Gather Meta data
    markers_callrate_removed <- (nL - ncol(Xa))


  } else {

    if(message){
    print("Alert: SNPS with high missing value not removed")
}

    Xa <- Xa
  }

  ##### Remove individual with bad call rate

  ### Remove markers with high missing value based on desired thresold
  if(!is.null(ind_rate)) {


    indR <- rowMeans(is.na(Xa))
    #poscr <- CR <= call_rate #selected by CR

    posindR <-  which(indR<ind_rate)

    Xa <- Xa[, posindR]

    ### Gather Meta data
    ind_callrate_removed <- (n_ind - nrow(Xa))


  } else {

    if(message){
    print("Alert: Individual with high missing value not removed")
}

    Xa <- Xa
  }


  filter <- hapmap[hapmap$`rs#`%in%colnames(Xa), ]

  rm(Xa)
  ## Begin process to remove heterozygo
  if(!is.null(het)) {
    heteroz <- apply(filter[, 12:ncol(filter)], 1, function(x){
      return(length(which(x%in%c('R', 'Y', 'S', 'W', 'K', 'M')))/length(x))
    })

    N_filter = ncol(filter)

    filter <- filter[heteroz < het, ]

    ## Gather meta data
    het_markers_removed = N_filter - ncol(filter)

  } else {

    print("Alert:High heterozygous SNPS not removed")

    filter = filter
  }

  #####
  ### Beging process to remove maf
  if(!is.null(maf)) {

    minor.allele.freq <- apply(filter[, c(2,12:ncol(filter)), with=FALSE], 1, function(x){
      allele1 <- substr(x[1], 1, 1)
      allele2 <- substr(x[1], 3, 3)
      if(length(which(x==allele1)) > length(which(x==allele2))){
        minor.allele <- allele2
      }
      else{minor.allele <- allele1}
      minor.allele.freq <- ((2*length(which(x[-1]==minor.allele))) + length(which(!x[-1]%in%c(allele1,allele2))))/(2*(length(x[-1])))
      return(minor.allele.freq)
    })

    N_filter = ncol(filter)

    filter <- filter[minor.allele.freq >= maf,]

    maf_markers_removed <-  (N_filter - ncol(filter))

  } else {

    if(message){
    print("Alert:SNPS with low maf not removed")
    }

    filter = filter

  }

  if(!is.null(recode)) {

    hapmap2numeric_Meth2_0 <- function(hapmap){
      hapmap.numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
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
      return(t(hapmap.numeric))
    }

    hapmap2numeric_Meth1_1 <- function(hapmap){
      hapmap.numeric <- apply(hapmap[, c(1, 2,12:ncol(hapmap)), with=FALSE], 1, function(x){
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
      return(t(hapmap.numeric))
    }

    if(recode== "Meth2_0"){
      # Here major allele is 2, minor is 0 and hetero is 1


      ### Convert to numeric format
      # Convert Hapmap allele into numeric

      #filter <- as.data.frame(filter)

      filter <- hapmap2numeric_Meth2_0(filter)
    } else {
      if(recode== "Meth-1_1"){

        filter <- hapmap2numeric_Meth1_1(filter)
      }

    }

    SNP_names <- filter[, 1]

    filter <- filter[, -c(1:2)]

    filter <- t(filter)

    colnames(filter) <- SNP_names

    #class(filter) <- "integer"

    if(isTRUE(impute)){
      if(any(is.na(filter))==T) {
        #This function will impute the missing values
        for(j in 1:ncol(filter)){
          if(j%%1000==0) cat("Marker=",j,"\n")
          tmp <- filter[,j]
          filter[,j] <- ifelse(is.na(tmp),round(mean(tmp,na.rm=T)),tmp)
        }

      }

    }
    # End


  } else {

    filter = filter

  }

  if(!is.null(recode)){

    class(filter) <- c("matrix", "array", "genotype")
  } else {

    class(filter) <- c("data.table", "data.frame", "genotype")
  }


  #### Aggregate all the maker data
  total_number_genotype = n_ind
  marker_callrate = call_rate
  if(is.null(marker_callrate)) {marker_callrate  = 0}
  ind_callrate = ind_rate
  if(is.null(ind_rate)) {ind_callrate  = 0}
  if(!exists("markers_callrate_removed")) {markers_callrate_removed  = 0}
  if(!exists("ind_callrate_removed ")) {ind_callrate_removed  = 0}
  heterozygosity = het
  if(is.null(het)) {heterozygosity  = 0}
  if(!exists("het_markers_removed ")) {het_markers_removed  = 0}
  maf = maf
  if(is.null(maf)) {maf  = 0}
  if(!exists("maf_markers_removed")) {maf_markers_removed = 0}


  output = list(filter,
                total_number_genotype,
                marker_callrate,
                markers_callrate_removed,
                ind_callrate,
                ind_callrate_removed,
                heterozygosity,
                het_markers_removed,
                maf,
                maf_markers_removed
  )

  if(!is.null(recode)){
    names(output) <- c("marker_matrix",
                       "total_number_genotype",
                       "marker_callrate",
                       "markers_callrate_removed",
                       "ind_callrate",
                       "ind_callrate_removed",
                       "heterozygosity",
                       "het_markers_removed",
                       "maf",
                       "maf_markers_removed")

  } else {

    names(output) <- c("marker",
                       "total_number_genotype",
                       "marker_callrate",
                       "markers_callrate_removed",
                       "ind_callrate",
                       "ind_callrate_removed",
                       "heterozygosity",
                       "het_markers_removed",
                       "maf",
                       "maf_markers_removed")


  }

  rm(filter)
  return(output)

  #return(filter)


}


