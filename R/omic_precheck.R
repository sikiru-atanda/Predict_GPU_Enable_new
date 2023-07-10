
#' Title
#'
#' @param object
#' @param message
#'
#' @return
#' @export
#'
#' @examples
omic_precheck <- function(object,
                          message = TRUE){
  msg <- sprintf("==================================================\n")
  if(!is.null(object)){
    if(class(object)[[1]]!="matrix"){

      if (message){
      message(paste(msg,"The data is not class matrix: we fix it." ))
      }

      object <- as.matrix(object)
    }



    if (is.numeric(object)==FALSE) {stop(print(paste(msg,'The data contains non-numeric value')), call. = FALSE)}

    ## Check for Na and remove
    Na_col.omit <- which((colSums(is.na(object))==0)==FALSE)

    if (length(Na_col.omit)!=0){
    #object = object[ , colSums(is.na(object))==0]

    object = object[ , -Na_col.omit]

    if (message){
      message("A total of ", length(Na_col.omit),
              " variable(s) / sample(s) were removed from due to missing value")

    }

    }

    #if (any(is.na(object))) {stop(print(paste(msg,'object contains Missing value')), call. = FALSE)}

    class(object) <-c("matrix", "array", "omic_matrix")

    attr(object, "cleared") <- "pass"

  } else {

    object = NULL

  }


  return(object)
}

