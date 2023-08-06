
#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param response
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
phenotype_precheck<- function(pheno_data = NULL,
                               gen_name = NULL,
                               response = NULL,
                               ...)
  {

      msg <- sprintf("==================================================\n")
    if(nrow(pheno_data)==0) { stop(print(paste(msg, 'No pheno_data records provided.')), call. = FALSE)

    }

    if (!inherits(pheno_data, what = 'data.frame')) {
      stop(print(paste(msg,"'pheno_data' must be of class 'data.frame'")), call. = FALSE)

    }

      # if (!data.table::is.data.table(pheno_data)){
      #   pheno_data <- data.table::as.data.table(pheno_data)
      # }

      if (!is.data.frame(pheno_data)){
         pheno_data <- as.data.frame(pheno_data)
       }

      ### Check the provided response name correspond to the name in the data file

      #if(!response%in%colnames(pheno_data)){
      if(sum(colnames(pheno_data)%in%response)<length(response)) {
        stop(print(paste(msg,paste(paste("The specified ",  response),
                                   " did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

      }

    if(!gen_name%in%colnames(pheno_data)){
      stop(print(paste(msg,paste(paste("The specified column",  gen_name),
                                 "in the pheno_data did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

    }
    ### Check to ensure no NA in the column GID/name
    if (anyNA(pheno_data[, gen_name]) || any(pheno_data[, gen_name]==-999)){

      stop(print(paste(msg,paste(paste('column',  gen_name),
                                 'should not have NA/missing'))), call. = FALSE)
    }

      if(!all(sapply(response, function(x, pheno_data) is.numeric(pheno_data[,x]),  pheno_data))) {
        pheno_data[, response] <-
          lapply(pheno_data[, response, drop = FALSE],
                 function(x) as.double(as.character(x)))

      }

      # Assign appropriate class.
      class(pheno_data) <- c("data.frame", "phenotype")

      attr(pheno_data, "cleared") <- "pass"


    return(pheno_data)

}
