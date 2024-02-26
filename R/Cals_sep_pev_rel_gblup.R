
#' Title
#'
#' @param geno_object
#' @param X
#' @param Z
#' @param va
#' @param ve
#' @param h2
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
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
