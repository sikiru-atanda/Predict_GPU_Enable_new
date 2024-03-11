#' Calculate Standard Error, PEV, and Reliability for GBLUP
#'
#' This function calculates the standard errors, prediction error variances (PEVs),
#' and reliabilities for genomic estimated breeding values (GEBVs) using genomic best
#' linear unbiased prediction (GBLUP). It adjusts the genomic relationship matrix with
#' a shrinkage factor and computes these metrics based on the adjusted matrix.
#'
#' @param geno_object Numeric matrix representing the genomic relationship matrix (GRM).
#' @param X Optional; design matrix for fixed effects.
#' @param Z Optional; design matrix for random effects, defaults to an identity matrix if NULL.
#' @param va Numeric; the additive genetic variance component.
#' @param ve Numeric; the residual (environmental) variance component.
#' @param h2 Optional; heritability on the observed scale. If not provided, it is calculated from `va` and `ve`.
#' @param lambda Optional; shrinkage parameter. If not provided, it is calculated from heritability (`h2`).
#'
#' @return A data frame containing three columns:
#'   - `rel`: Reliability of the GEBVs.
#'   - `pev`: Prediction error variance of the GEBVs.
#'   - `sep`: Standard error of the prediction for the GEBVs.
#'
#' @examples
#' # Assuming geno_object, va, and ve are predefined:
#' results <- sep_pev_rel_gblup(geno_object = GRM,
#'                              va = 0.5,
#'                              ve = 0.5)
#' @export
#' @importFrom stats solve

sep_pev_rel_gblup <- function(geno_object = NULL,
                              X = NULL,
                              Z = NULL,
                              va = NULL,
                              ve = NULL,
                              h2 = NULL,
                              lambda = NULL){


  if(is.null(Z)){
     Z <- diag(1,nrow(geno_object))
}

dimnames(Z) <- list(rownames(geno_object),
                    colnames(geno_object))
if((is.null(h2) & is.null(lambda)) & (!is.null(va) & !is.null(ve))){
h2 <- va/(va+ve)
lambda <- (1-h2)/h2
} else if ((!is.null(h2) & is.null(lambda)) & (is.null(va) & is.null(ve))){
  lambda <- (1-h2)/h2
} else if ((is.null(h2) & !is.null(lambda)) & (is.null(va) & is.null(ve))){

  lambda <- lambda

} else {
  if ((is.null(h2) & !is.null(lambda)) & (is.null(va) & is.null(ve))){

    stop(message(paste(msg,'Provide genetic and residual variance.')), call. = FALSE)

  }

}
# Adjusting Genomic Relationship Matrix (C):
# Matrix CC is constructed to adjust the geno_object using
# the shrinkage factor.
# The inverse of the adjusted matrix is calculated.
ZZG <- crossprod(Z,Z) + solve(geno_object)*lambda
ZZG <- solve(ZZG)
ZZG_diag <- diag(ZZG)
##
# Standard errors (sepsep) are calculated as the square root of
# the product of individual variances and the residual variance.
# Prediction error variances (pevpev) are calculated as
# the product of individual variances and the residual variance.
# Reliability (relrel) is calculated as 1- minus the product
# of individual shrinkage factors and the shrinkage factor.
pev <- (ZZG_diag*ve)
sep <- sqrt(pev) ## Standard error
rel <- 1-ZZG_diag*lambda
rel = ifelse(rel<0, NA, rel)

output <- data.frame(rel=rel,pev=pev,sep=sep)


return(output)


}
