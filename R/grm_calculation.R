

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
    method=c("VanRaden",
             "Yang")
){

  msg <- sprintf("==================================================\n")


  if(class(geno_clean)[1]!= "matrix") stop(print(paste(msg, 'object geno_clean must be matrix.')), call. = FALSE)

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

    locusMat <- scale(geno_clean, center=T, scale=F)
    return(tcrossprod(geno_clean) / sum(2*freq*(1-freq)))
  }

  ### Yang method
  Yang <- function(geno_clean){

    freq <- colMeans(geno_clean)/2

    locusMat <- scale(x = geno_clean, center = T, scale = F)
    locusMat <- (1/N_marker)*(locusMat %*% (t(locusMat) * (1/(2*freq*(1-freq)))))
    locusMat[lower.tri(locusMat, diag = T)] <- 0
    locusMat <- locusMat + t(locusMat)

    multiplier <- geno_clean^2 - t(t(geno_clean) * (1+2*freq)) + matrix(rep(2*freq^2, each=N_Individuals), ncol=N_marker)
    diag(locusMat) <- 1+(1/N_marker)*colSums(t(multiplier) *  (1/(2*freq*(1-freq))))

    return(locusMat)


  }



  switch(method,
         "VanRaden" = {
           Ga <- VanRaden(geno_clean, freq)
         },
         "Yang" = {
           Ga <- Yang(geno_clean)
         },
         {

           stop(print(paste(msg, "Select method to calculate geno_cleanmic relationship matrix")), call. = FALSE)
         })


  return(Ga)

}

