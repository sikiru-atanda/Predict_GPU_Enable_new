
#' Title
#'
#' This function address the following objectives:
#' 1. Check if the phenotypic data file contain records
#' 2. check if the phenotypic data is in data.frame and if not we fit it for the user
#' 3. Check if the gen_name and response names provided by the user match the column
#'    names in the phenotypic data
#' 4. Check that the gen_name did not has NA as part of the genotypes ID
#' 5. Check that the response variable is numeric if not convert it to numeric
#' 7. Check for NA is not test here for response variables because testing set might have NA
#'    thus, NA is expected.
#' 6. If all these checks were fulfilled assign class(pheno_data) <- c("data.frame", "phenotype")
#' and define attribute as attr(pheno_data, "cleared") <- "pass"
#'
#' 7. This output will be input for phenotype_to_model function
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
phenotype_precheckOLD<- function(pheno_data = NULL,
                              gen_name = NULL,
                              response = NULL,
                              heter_groups = NULL,
                               ...)
  {

      msg <- "\n==================================================\n"
    if(nrow(pheno_data)==0) { stop(print(paste(msg, 'No pheno_data records provided.')), call. = FALSE)

    }

    if (!inherits(pheno_data, what = 'data.frame')) {
      #stop(print(paste(msg,"'pheno_data' must be of class 'data.frame'")), call. = FALSE)

      message(print(paste(msg,"'pheno_data' is not 'data.frame' type. We fix it.")))
      pheno_data <- as.data.frame(pheno_data)
    }

      # if (!data.table::is.data.table(pheno_data)){
      #   pheno_data <- data.table::as.data.table(pheno_data)
      # }

      # if (!is.data.frame(pheno_data)){
      #    pheno_data <- as.data.frame(pheno_data)
      #  }

      ### Check the provided response name correspond to the name in the data file

      #if(!response%in%colnames(pheno_data)){
      if(sum(colnames(pheno_data)%in%response)<length(response)) {
        stop(print(paste(msg,paste(paste("The specified ",  response),
                                   " did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

      }

      #### This is very much imp especially for BGLR
      if(!is.null(heter_groups)) {
        if(heter_groups%in%colnames(pheno_data))
        stop(print(paste(msg,paste(paste("The specified ",  heter_groups),
                                   " did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

        pheno_data = pheno_data[order(pheno_data[[heter_groups]]), ]

      }

    if(!gen_name%in%colnames(pheno_data)){
      stop(print(paste(msg,paste(paste("The specified column",  gen_name),
                                 "in the pheno_data did not match with your data.\n\t\t Please check and use apppropriatly"))), call. = FALSE)

    }

      if(var(pheno_data[[response]], na.rm = TRUE)==0){
        stop(print(paste(msg,paste(response, "variable has zero variance.\n\t This trait cannot be used for prediction model.\n\t\t Check the raw data and model that generate the BLUEs."))), call. = FALSE)

      }
    ### Check to ensure no NA in the column GID/name
    if (anyNA(pheno_data[[gen_name]]) || any(pheno_data[[gen_name]]==-999)){

      stop(print(paste(msg,paste(paste('column',  gen_name),
                                 'should not have NA/missing'))), call. = FALSE)
    }

      if(!all(sapply(response, function(x, pheno_data) is.numeric(pheno_data[,x]),  pheno_data))) {
        pheno_data[[response]] <-
          lapply(pheno_data[, response, drop = FALSE],
                 function(x) as.double(as.character(x)))

      }

      # Assign appropriate class.
      #class(pheno_data) <- c("data.frame", "phenotype")

      attr(pheno_data, "cleared") <- "pass"


    return(pheno_data)

}
