
#' Title
#'
#' @param object
#' @param message
#'
#' @return
#' @export
#'
#' @examples
omic_precheck <- function(object = NULL,
                          message = TRUE,
                          impute = FALSE){
  msg <- "\n==================================================\n"
  if(!is.null(object)){
    if("data.table" %in% class(object)){
      stop(print(paste(msg,'Omic data must be data.frame or matrix not character.')), call. = FALSE)
    }
    if(inherits(object, "character")) stop(print(paste(msg,'Omic data should be data.frame or matrix not character.')), call. = FALSE)
    #if(class(object)[[1]]!="matrix"){
    if(!inherits(object, "matrix")){
      if(isTRUE(message)){
      message(paste(msg,"The data is not class matrix: we fix it." ))
      }

      object <- as.matrix(object)
    }



    #if (is.numeric(object)==FALSE) {stop(print(paste(msg,'The data contains non-numeric value')), call. = FALSE)}
    # Check if all elements in the matrix are numeric
    all_numeric <- all(apply(object, c(1, 2), is.numeric))

    # Stop execution if any element is not numeric
    if (!all_numeric) {
      stop('The omic data contains non-numeric values', call. = FALSE)
    }
    ## Check for Na and remove
    object <- handle_missing_values(data = object) ## this fxn is present in geno_precheckk
    ## check if there is duplicated snps
    duplicated_columns <- colnames(object)[duplicated(colnames(object))]
    if(length(duplicated_columns)>0){
      stop(print(paste(msg,"Omic data contain duplicate features/predictors.")), call. = FALSE)
    }

    ## check duplicated rownames:
    duplicated_rownames <- rownames(object)[duplicated(rownames(object))]
    if(length(duplicated_rownames)>0){
      stop(print(paste(msg,"Omic data contain duplicate genotypes.")), call. = FALSE)
    }
    # ## Check for Na and remove
    # Na_col.omit <- which((colSums(is.na(object))==0)==FALSE)
    #
    # if (length(Na_col.omit)!=0){
    # #object = object[ , colSums(is.na(object))==0]
    #   if(isTRUE(impute)){
    #     object <- apply(object, 2, function(coll) {
    #       col_mean <- mean(coll, na.rm = TRUE)  # Calculate mean excluding NA
    #       coll[is.na(coll)] <- round(col_mean)  # Replace NA with the rounded mean
    #       return(coll)
    #     })
    #
    #     if (isTRUE(message)){
    #       message("A total of ", length(Na_col.omit),
    #               " variable(s) / sample(s) were detected with missing value and imputed.\n If you wish to remove set impute to False")
    #
    #     }
    #   } else {
    #
    # object <- object[ , -Na_col.omit]
    #
    # if (isTRUE(message)){
    #   message("A total of ", length(Na_col.omit),
    #           " variable(s) / sample(s) were removed from due to missing value")
    #
    # }
    #
    #   }
    # }

    #if (any(is.na(object))) {stop(print(paste(msg,'object contains Missing value')), call. = FALSE)}

    #class(object) <-c("matrix", "array", "omic_matrix")

    attr(object, "cleared") <- "pass"

  } else {

    object <- NULL

  }


  return(object)
}

