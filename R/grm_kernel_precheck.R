
#' Title
#' pheno_data is a genomic relationship matrix or relationship matrix calculated using
#' other omics data. The relationship matrix can be calculated using any method.
#' NA is not allowed in the relationship matrix
#' ####
#' Checks
#' #######
#' 1. It check if the matrix is square matrix/symmetry, if not we fix it for the user
#' 2. It check if the matrix is positive definite, if not we fix it.
#' 3. It check for NA. If present the engine will stop further analysis.
#'
#'
#' @param pheno_data
#' @param ...
#' @param message
#'
#' @return
#' @export
#'
#' @examples
grm_kernel_precheck <- function(pheno_data= NULL,
                                message= TRUE,
                                bending = TRUE,
                                ...){

  msg <- sprintf("==================================================\n")

if (!is.null(pheno_data)){
  if(isTRUE(anyNA(pheno_data))){
    stop(print(paste(msg,'NA is not allowed in the grm or kernel matrix')), call. = FALSE)
    }
  ### Check if the Grm matrix is in class matrix if not convert to class matrix
  if (!is.matrix(pheno_data)) pheno_data <- as.matrix(pheno_data)
  ### Check if colname and rownames in grm/kernel matrix is the same
  if (!identical(colnames(pheno_data), rownames(pheno_data))) {stop(print(paste(msg,'colnames did not match rownames')), call. = FALSE)}
  ### Check if the grm/kernel matrix is symmetric

  #### To agree with the team, if we want to put a stop if the matrix is not
  # symmetric or fix it for the user at the expensive of few memory.
  # Keep in mind our users are poor, better and best
  #if(!isSymmetric.matrix(pheno_data)) {stop(print(paste(msg,'pheno_data is not symmetric')), call. = FALSE)}
  if(!isSymmetric.matrix(pheno_data)) {

    message(paste(msg,"Relationsip Matrix is not symmetric'. We fix it"))
    pheno_data <- Matrix::forceSymmetric(pheno_data)
    pheno_data <- Matrix::as.matrix(pheno_data)
  }

  if(isTRUE(bending)){
  if(isFALSE(matrixcalc::is.positive.definite(pheno_data))){

    message(paste(msg,"Relationsip Matrix is not positive definite. We fix it"))

    pheno_data <- as.matrix(Matrix::nearPD(pheno_data, posd.tol=1e-02, trace=FALSE)$mat)

  }

  } else {

    if(isFALSE(matrixcalc::is.positive.definite(pheno_data))){

      message(paste(msg,"Relationsip Matrix is not positive definite. Set bending = TRUE to fix it"))

    }

  }

}

  #### Declare it also as an pheno_data for final usage
  class(pheno_data) <-c("matrix", "array", "krm_data")

  attr(pheno_data, "cleared") <- "for_model_fit"
  ##############################################
  return(pheno_data)

}






