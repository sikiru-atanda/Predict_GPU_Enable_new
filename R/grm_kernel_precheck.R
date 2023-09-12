
#' Title
#' grm_kernel_data is a genomic relationship matrix or relationship matrix calculated using
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
#' @param grm_kernel_data
#' @param ...
#' @param message
#'
#' @return
#' @export
#'
#' @examples
grm_kernel_precheck <- function(grm_kernel_data= NULL,
                                message= TRUE,
                                bending = TRUE,
                                ...){

  msg <- sprintf("==================================================\n")

if (!is.null(grm_kernel_data)){
  if(isTRUE(anyNA(grm_kernel_data))){
    stop(print(paste(msg,'NA is not allowed in the grm or kernel matrix')), call. = FALSE)
    }
  ### Check if the Grm matrix is in class matrix if not convert to class matrix
  if (!is.matrix(grm_kernel_data)) grm_kernel_data <- as.matrix(grm_kernel_data)
  ### Check if colname and rownames in grm/kernel matrix is the same
  if (!identical(colnames(grm_kernel_data), rownames(grm_kernel_data))) {stop(print(paste(msg,'colnames did not match rownames')), call. = FALSE)}
  ### Check if the grm/kernel matrix is symmetric

  #### To agree with the team, if we want to put a stop if the matrix is not
  # symmetric or fix it for the user at the expensive of few memory.
  # Keep in mind our users are poor, better and best
  #if(!isSymmetric.matrix(grm_kernel_data)) {stop(print(paste(msg,'grm_kernel_data is not symmetric')), call. = FALSE)}
  if(!isSymmetric.matrix(grm_kernel_data)) {

    message(paste(msg,"Relationsip Matrix is not symmetric'. We fix it"))
    grm_kernel_data <- Matrix::forceSymmetric(grm_kernel_data)
    grm_kernel_data <- Matrix::as.matrix(grm_kernel_data)
  }

  if(isTRUE(bending)){
  if(isFALSE(matrixcalc::is.positive.definite(grm_kernel_data))){

    message(paste(msg,"Relationsip Matrix is not positive definite. We fix it"))

    grm_kernel_data <- as.matrix(Matrix::nearPD(grm_kernel_data, posd.tol=1e-02, trace=FALSE)$mat)

  }

  } else {

    if(isFALSE(matrixcalc::is.positive.definite(grm_kernel_data))){

      message(paste(msg,"Relationsip Matrix is not positive definite. Set bending = TRUE to fix it"))

    }

  }

}

  #### Declare it also as an grm_kernel_data for final usage
  class(grm_kernel_data) <-c("matrix", "array", "krm_data")

  attr(grm_kernel_data, "cleared") <- "for_model_fit"
  ##############################################
  return(grm_kernel_data)

}






