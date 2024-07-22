
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

      object_pheno = object_pheno[["pheno_clean_data"]]

      # if(anyNA(object_pheno[, response])){
      #   Na_testing <- which(is.na(object_pheno[, response]))
      #
      # }
      if(any(is.na(object_pheno[, response]))) {
        # Get row indices where any of the specified columns have NAs
        Na_testing <- which(rowSums(is.na(object_pheno[, response])) > 0)

      }

      if(exists("Na_testing")){

        test_set <- object_pheno[Na_testing, ]

        object_pheno <- object_pheno[-Na_testing, ]

      } else {

        test_set <-  NULL
      }


    } else {

      object_pheno <- object_pheno[["pheno_clean_data"]]

      test_set <- object_pheno[["test_set"]]

    }

  if(exists("test_set") & !is.null(test_set)){

  #if(!is.null(test_set) & exists("object_pheno")){

    #if(!is.null(test_set)){

  output <- list(object_pheno, test_set)

  names(output) <- c("pheno_clean_data", "test_set")
    #}

  } else {

    #if(!exists("test_set") & exists("object_pheno")){

  output <- list(object_pheno)

  names(output) <- "pheno_clean_data"

    #}

  }



  return(output)

}
