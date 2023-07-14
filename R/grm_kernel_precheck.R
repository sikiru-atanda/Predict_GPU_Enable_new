
#' Title
#'
#' @param object
#' @param ...
#' @param message
#'
#' @return
#' @export
#'
#' @examples
grm_kernel_precheck <- function(object= NULL,
                                message= TRUE,
                                ...){

  msg <- sprintf("==================================================\n")

if (!is.null(object)){
  if(isTRUE(anyNA(object))){
    stop(print(paste(msg,'NA is not allowed in the grm or kernel matrix')), call. = FALSE)
    }
  ### Check if the Grm matrix is in class matrix if not convert to class matrix
  if (!is.matrix(object)) object <- as.matrix(object)
  ### Check if colname and rownames in grm/kernel matrix is the same
  if (!identical(colnames(object), rownames(object))) {stop(print(paste(msg,'colnames did not match rownames')), call. = FALSE)}
  ### Check if the grm/kernel matrix is symmetric

  #### To agree with the team, if we want to put a stop if the matrix is not
  # symmetric or fix it for the user at the expensive of few memory.
  # Keep in mind our users are poor, better and best
  #if(!isSymmetric.matrix(object)) {stop(print(paste(msg,'object is not symmetric')), call. = FALSE)}
  if(!isSymmetric.matrix(object)) {
    object <- Matrix::forceSymmetric(object)
    object <- Matrix::as.matrix(object)
  }
}

  #### Declare it also as an object for final usage
  class(object) <-c("matrix", "array", "krm_data")

  attr(object, "cleared") <- "for_model_fit"
  ##############################################
  return(object)

}






