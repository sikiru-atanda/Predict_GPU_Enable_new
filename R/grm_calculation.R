

#' Title
#'
#' @param geno_clean
#' @param method
#'
#' @return
#'
#' @examples
#'
grm_calculation <- function(
    geno_clean = NULL,
    weight = NULL,
    method=c("VanRaden",
             "Weighted_VanRaden",
             "Yang",
             "Epistasis")
){

  msg <- sprintf("==================================================\n")

  ### iT important to check the name in the weight data is the match and the same
  ## order in the geno_clean data

 if(class(geno_clean)[1]!= "matrix") stop(print(paste(msg, 'object geno_clean must be matrix.')), call. = FALSE)

  if(!is.null(weight)){

    # Literature
    # Weighting Strategies for Single-Step Genomic BLUP: An Iterative Approach
    # for Accurate Calculation of GEBV and  GWAS

    # D (weight) is a diagonal matrix of weights, where dii is the weight
    # for SNP i. In regular GBLUP-based methods, D = I, which gives
    # a weight of 1 to all SNP.


   if(class(weight)[1]!= "matrix") stop(print(paste(msg, 'weight must be matrix.')), call. = FALSE)
   ### It is possible when geno_clean was QC some snp did not meet the standard QC parameters
   ## and were dropped.This will fit that deficit
    ## This is lacking in AGHmatrix.
   if(isFALSE(all(rownames(weight)%in%colnames(geno_clean)))){

     ## This ensure that the order of the snp in geno_clean match that in the weight.
     weight = weight[colnames(geno_clean)%in%rownames(weight),]

     weight <- diag(x= weight[,1], nrow = nrow(weight))

   }


   ## This check if the snp name and that in the weight did not match or not the same order.
   if (!identical(colnames(geno_clean), rownames(weight))) stop(print(paste(msg, 'Snp marker should be equivalent in both weight geno_data.')), call. = FALSE)
 }
  #if(class(geno_clean)[1]!= "matrix") stop(print(paste(msg, 'object geno_clean must be matrix.')), call. = FALSE)

  freq <- colMeans(geno_clean)/2
  if (any(is.na(geno_clean) | freq == 0 | freq == 1)) geno_clean= Remove_NA_Mono_SNP(geno_clean= geno_clean)

  ## Check if SNP data is coded 0, 1, 2

  checkG <- c(length(which(geno_clean == -1)))


  if (checkG!=0) stop(print(paste(msg, "SNP data must be coded 0, 1, 2")), call. = FALSE)


  if(missing(method)) stop(print(paste(msg, "Select either VanRaden or Yang to compute geno_cleanmic relationship matrix")), call. = FALSE)


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

    geno_clean <- scale(geno_clean, center=T, scale=F)
    return(((geno_clean %*% weight)%*% t(geno_clean))/sum(2*freq*(1-freq)))
  }

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
         {

           stop(print(paste(msg, "Select method to calculate geno_cleanmic relationship matrix")), call. = FALSE)
         })


  return(Ga)

}

