
#' Calculate Genomic Relationship Matrix (GRM)
#'
#' This function computes the genomic relationship matrix using various methods,
#' including VanRaden, Weighted VanRaden, Yang, and methods that apply weights to
#' standardized genotypes or allele frequencies. It is designed to work with
#' matrix data representing genotypes and can incorporate weights to adjust
#' the contribution of each SNP.
#'
#' @param geno_clean A matrix of genotype data where rows represent individuals
#'   and columns represent SNPs. Genotype data must be coded as 0, 1, 2.
#' @param weight An optional matrix of weights for SNPs, used in weighted
#'   calculation methods. If provided, must match the dimensions and order
#'   of `geno_clean`.
#' @param method A character vector specifying the method to use for GRM
#'   calculation. Options include "VanRaden", "Weighted_VanRaden", "Yang",
#'   "weighted_GRM", and "weightedGRM_AlleleFreq".
#'
#' @return A matrix representing the genomic relationship matrix calculated
#'   based on the specified method.
#'
#' @examples
#' # Example genotype data (geno_clean)
#' geno_clean <- matrix(c(0, 1, 2, 1, 0, 1, 2, 2, 1), nrow = 3, byrow = TRUE)
#' colnames(geno_clean) <- c("SNP1", "SNP2", "SNP3")
#' rownames(geno_clean) <- c("Ind1", "Ind2", "Ind3")
#'
#' # Calculate GRM using VanRaden method
#' grm_vanraden <- grm_calculation(geno_clean, method = "VanRaden")
#'
#' # Calculate GRM using a weighted method
#' weights <- diag(c(1, 1, 1))
#' grm_weighted <- grm_calculation(geno_clean, weight = weights, method = "Weighted_VanRaden")
#'
#' @details The function supports several methods for GRM calculation, with
#'   "VanRaden" being the default. Weighted methods adjust the influence of
#'   each SNP based on provided weights, which can be useful for incorporating
#'   prior knowledge about SNP effects.
#'
#'   The "weighted_GRM" method applies weights directly to standardized genotypes
#'   before calculating the GRM, assuming the weights reflect the relative
#'   importance of each SNP.
#'
#'   The "weightedGRM_AlleleFreq" method adjusts allele frequencies by weights
#'   before calculating the GRM, assuming SNPs with different frequencies
#'   contribute differently to the genetic variance.
#'
#' @note Ensure that genotype data (`geno_clean`) is correctly formatted and
#'   coded as 0, 1, 2. The function will stop with an error if genotype data
#'   includes NA values or incorrect coding.
#'

grm_calculation <- function(
    geno_clean = NULL,
    weight = NULL,
    method= NULL
){

  msg <- "\n==================================================\n"

  ### iT important to check the name in the weight data is the match and the same
  ## order in the geno_clean data

  # Validate input
  if(!is.matrix(geno_clean)) stop(paste(msg, "snp/marker data must be a matrix."), call. = FALSE)
  if(!is.null(weight) && !is.matrix(weight)) stop(paste(msg, "weight must be a matrix if provided."), call. = FALSE)
  if(!is.null(weight) && !all(rownames(weight) %in% colnames(geno_clean))) {
    stop(paste(msg,"Not all SNPs in weight are present in snp/marker data."), call. = FALSE)
  }

  gmatrix_method_available <- c("VanRaden",
                                "Weighted_VanRaden",
                                "Yang",
                                "Epistasis")

  if(!is.null(method)){
    if (!(method %in% gmatrix_method_available)) {
      stop(paste(msg,"Invalid genomic relationship method. Choose from: ",
           paste(gmatrix_method_available, collapse = ", ")), call. = FALSE)
    }
  }

  if(!is.null(weight)){

    weight <- weight[match(colnames(geno_clean), rownames(weight)), , drop = FALSE]
    if(!identical(colnames(geno_clean), rownames(weight))) stop(paste(msg, "SNP order in weight does not match geno_clean."), call. = FALSE)
    #weight <- diag(as.vector(weight))

    # Literature
    # Weighting Strategies for Single-Step Genomic BLUP: An Iterative Approach
    # for Accurate Calculation of GEBV and  GWAS

    # D (weight) is a diagonal matrix of weights, where dii is the weight
    # for SNP i. In regular GBLUP-based methods, D = I, which gives
    # a weight of 1 to all SNP.

   ### It is possible when geno_clean was QC some snp did not meet the standard QC parameters
   ## and were dropped.This will fit that deficit
    ## This is lacking in AGHmatrix.
   # if(isFALSE(all(rownames(weight)%in%colnames(geno_clean)))){
   #
   #   ## This ensure that the order of the snp in geno_clean match that in the weight.
   #   weight = weight[colnames(geno_clean)%in%rownames(weight),]
   #
   #   weight <- diag(x= weight[,1], nrow = nrow(weight))
   #
   # }
   # ## This check if the snp name and that in the weight did not match or not the same order.
   # if (!identical(colnames(geno_clean), rownames(weight))) stop(print(paste(msg, 'Snp marker should be equivalent in both weight geno_clean.')), call. = FALSE)
 }
  #if(class(geno_clean)[1]!= "matrix") stop(print(paste(msg, 'object geno_clean must be matrix.')), call. = FALSE)

  freq <- colMeans(geno_clean)/2
  if (any(is.na(geno_clean) | freq == 0 | freq == 1)) geno_clean <- Remove_NA_Mono_SNP(geno_clean= geno_clean)

  ## Check if SNP data is coded 0, 1, 2

  checkG <- c(length(which(geno_clean == -1)))


  if (checkG!=0) stop(print(paste(msg, "SNP data must be coded 0, 1, 2")), call. = FALSE)


  if(missing(method)) stop(paste(msg, "Select either VanRaden or Yang to compute geno_cleanmic relationship matrix"), call. = FALSE)


  N_Individuals <- nrow(geno_clean)  ## Number of informative SNP
  N_marker <- ncol(geno_clean)  ## Number of geno_cleantypes

  ### Vanraden method
  VanRaden <- function(geno_clean, freq){
    #freq <- colMeans(geno_clean) / 2

    geno_clean <- scale(geno_clean, center=T, scale=F)
    return(tcrossprod(geno_clean) / sum(2*freq*(1-freq)))
  }

  Epistasis <- function(geno_clean, freq){
    #freq <- colMeans(geno_clean) / 2

    geno_clean <- scale(geno_clean, center=T, scale=F)

    G <- tcrossprod(geno_clean) / sum(2*freq*(1-freq))

    GG <- G*G

    return(GG)
  }

  ### Weighted_VanRaden method
  Weighted_VanRaden <- function(geno_clean, freq, weight){
    #freq <- colMeans(geno_clean) / 2
    weightt <- diag(as.vector(weight))
    geno_clean <- scale(geno_clean, center=T, scale=F)
    return(((geno_clean %*% weightt)%*% t(geno_clean))/sum(2*freq*(1-freq)))
  }
###
  # In this method, weights are directly applied to the standardized genotypes
  # before calculating the GRM. This approach assumes that the
  # weights reflect the relative importance or effect size of each SNP.

  weighted_GRM <- function(geno_clean, weight) {
    # Standardize genotypes (centering)
    standardizedGeno <- scale(geno_clean, center = TRUE, scale = FALSE)
    # Apply weights
    weightedGeno <- standardizedGeno * sqrt(weights)
    # Calculate GRM
    grm <- tcrossprod(weightedGeno) / ncol(weightedGeno)
    return(grm)
  }

  ###
  # This method involves adjusting allele frequencies by
  # weights before calculating the GRM, reflecting the
  # assumption that SNPs with different frequencies
  # contribute differently to the genetic variance.
  #
  weightedGRM_AlleleFreq <- function(geno_clean, weights) {
    freq <- colMeans(geno_clean) / 2
    # Apply weights to allele frequencies
    weightedFreq <- freq * weights
    # Calculate GRM using weighted allele frequencies
    p <- 2 * weightedFreq * (1 - weightedFreq)
    standardizedGeno <- scale(geno_clean, center = TRUE, scale = FALSE)
    grm <- tcrossprod(standardizedGeno) / sum(p)
    return(grm)
  }

  ####
  # Convert the genotype matrix to a dominance matrix
  convert_to_dominance <- function(geno_clean) {
    # Convert genotypes to dominance deviations
    # Heterozygotes (1) have dominance effect, homozygotes (0, 2) do not
    D <- ifelse(geno_clean == 1, 1, 0)

    # Calculate the dominance relationship matrix
    DDM <- tcrossprod(D) / ncol(D)

    return(DDM)
  }
  ##

  ### Yang method
  Yang <- function(geno){

    freq <- colMeans(geno)/2
    N_marker = ncol(geno)
    N_Individuals = nrow(geno)

    locusMat <- scale(x = geno, center = T, scale = F)
    locusMat <- (1/N_marker)*(locusMat %*% (t(locusMat) * (1/(2*freq*(1-freq)))))
    locusMat[lower.tri(locusMat, diag = T)] <- 0
    locusMat <- locusMat + t(locusMat)

    multiplier <- geno^2 - t(t(geno) * (1+2*freq)) + matrix(rep(2*freq^2, each=N_Individuals), ncol=N_marker)
    diag(locusMat) <- 1+(1/N_marker)*colSums(t(multiplier) *  (1/(2*freq*(1-freq))))

    return(locusMat)


  }

  Epistasis

  switch(method,
         "VanRaden" = {
           Ga <- VanRaden(geno_clean, freq)
         },
         "Weighted_VanRaden" = {
           Ga <- Weighted_VanRaden(geno_clean, freq, weight)
         },
         "Yang" = {
           Ga <- Yang(geno_clean)
         },
         "weighted_GRM"={
           Ga <- weighted_GRM(geno_clean, weight)
         },
         "Epistasis"={
           Ga <- Epistasis(geno_clean, freq)
         },
         "weightedGRM_AlleleFreq"={
           Ga <- weightedGRM_AlleleFreq(geno_clean, weight)
         },
         {

           stop(paste(msg, "Select method to calculate geno_cleanmic relationship matrix"), call. = FALSE)
         })


  return(Ga)

}

# pValues <- c(...)  # Vector of p-values from GWAS
# zScores <- qnorm(pValues / 2, lower.tail = FALSE)  # Convert p-values to Z-scores
#
# weights <- zScores^2
#
#
# Converting P-values to Z-scores
#
# First, convert the GWAS p-values to Z-scores.
# The Z-score represents the deviation of the observed association from
# the null hypothesis (no association), in units of the standard error.
# Larger Z-scores (either positive or negative) indicate stronger evidence
# against the null hypothesis.
#
# Using Z-scores as Weights
#
# One approach is to use the square of the Z-scores as weights.
# The rationale behind this is that the square of a Z-score is
# proportional to the chi-square statistic,
# which in turn is related to the
# variance explained by the SNP under certain assumptions.
