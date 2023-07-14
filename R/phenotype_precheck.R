
#' Title
#'
#' @param pheno
#' @param gen_name
#' @param response
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
phenotype_precheck <- function(pheno = NULL,
                               gen_name = NULL,
                               response = NULL,
                               ...)
  {

      msg <- sprintf("==================================================\n")
    if(nrow(pheno)==0) { stop(print(paste(msg, 'No phenotypic records provided.')), call. = FALSE)

    }

    if (!inherits(pheno, what = 'data.frame')) {
      stop(print(paste(msg,"'phenotype' must be of class 'data.frame'")), call. = FALSE)

    }

      # if (!data.table::is.data.table(pheno)){
      #   pheno <- data.table::as.data.table(pheno)
      # }

      if (!is.data.frame(pheno)){
         pheno <- as.data.frame(pheno)
       }

      ### Check the provided response name correspond to the name in the data file

      if(!response%in%colnames(pheno)){

        stop(print(paste(msg,paste(paste("The specified ",  response),
                                   " did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

      }

    if(!gen_name%in%colnames(pheno)){
      stop(print(paste(msg,paste(paste("The specified column",  gen_name),
                                 "in the phenotypic data did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

    }
    ### Check to ensure no NA in the column GID/name
    if (anyNA(pheno[, gen_name]) || any(pheno[, gen_name]==-999)){

      stop(print(paste(msg,paste(paste('column',  gen_name),
                                 'should not have NA/missing'))), call. = FALSE)
    }

      # Assign appropriate class.
      class(pheno) <- c("data.frame", "phenotype")

      attr(pheno, "cleared") <- "pass"


    return(pheno)

}
