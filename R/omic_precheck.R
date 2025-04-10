
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
                          impute=TRUE,
                          na_threshold =0.9,
                          imputation_method = "knn", #median, mean
                          impute_knn_k = 5){
  msg <- "\n==================================================\n"
  if(!is.null(object)){
    if("data.table" %in% class(object)){
      stop(paste(msg,'Omic data must be data.frame or matrix not character.'), call. = FALSE)
    }
    if(inherits(object, "character")) stop(paste(msg,'Omic data should be data.frame or matrix not character.'), call. = FALSE)
    #if(class(object)[[1]]!="matrix"){
    if(!inherits(object, "matrix")){
      if(isTRUE(message)){
      message(paste(msg,"The data is not class matrix: we fix it." ))
      }

      object <- as.matrix(object)
    }

    if(any(colnames(object) %in% c("NA", "Na", "na"))){
      stop(paste(msg, "Column names can't contain NA."), call. = FALSE)
    }

    #if (is.numeric(object)==FALSE) {stop(print(paste(msg,'The data contains non-numeric value')), call. = FALSE)}
    # Check if all elements in the matrix are numeric
    all_numeric <- all(apply(object, c(1, 2), is.numeric))

    # Stop execution if any element is not numeric
    if (!all_numeric) {
      stop(paste(msg,'The omic data contains non-numeric values', call. = FALSE), call. = FALSE)
    }
    ## Check for Na and remove
    object <- handle_missing_values(data = object,  na_threshold =  na_threshold,
                                    impute=impute, imputation_method = imputation_method,
                                    impute_knn_k = impute_knn_k) ## this fxn is present in geno_precheckk
    ## check if there is duplicated snps
    duplicated_columns <- colnames(object)[duplicated(colnames(object))]
    if(length(duplicated_columns)>0){
      stop(paste(msg,"Omic data contain duplicate features/predictors."), call. = FALSE)
    }

    ## check duplicated rownames:
    duplicated_rownames <- rownames(object)[duplicated(rownames(object))]
    if(length(duplicated_rownames)>0){
      stop(paste(msg,"Omic data contain duplicate genotypes."), call. = FALSE)
    }

    zero_var_check <- caret::nearZeroVar(as.matrix(object))
    if (length(zero_var_check) > 0) {
      object <- object[, -zero_var_check]
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

