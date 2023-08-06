
#' Title
#'
#' @param object_pheno
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ML_undefined_test_train <- function(object_pheno = NULL,
                                    response = NULL,
                                    ...)
  {
    if(length(object_pheno)==1){

      object_pheno = object_pheno$pheno_data

      if(anyNA(object_pheno[, response])){
        Na_testing <- which(is.na(object_pheno[, response]))

      }

      if(exists("Na_testing")){

        test_set <- object_pheno[Na_testing, ]

        object_pheno <- object_pheno[-Na_testing, ]

      }


    } else {

      object_pheno <- object_pheno$pheno_data

      test_set <- object_pheno$test_set

    }

  if(exists("test_set")){

  #if(!is.null(test_set) & exists("object_pheno")){

    #if(!is.null(test_set)){

  output <- list(object_pheno, test_set)

  names(output) <- c("pheno_data", "test_set")
    #}

  } else {

    #if(!exists("test_set") & exists("object_pheno")){

  output <- list(object_pheno)

  names(output) <- "pheno_data"

    #}

  }



  return(output)

}
