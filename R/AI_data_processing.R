
#' Title
#'
#' @param geno_data
#' @param omic_data  list of omics
#'
#' @return
#' @export
#'
#' @examples
#'

merge_data <- function(geno_data, omic_data) {
  # Ensure omic_data is a list for consistency
  # if (!is.list(omic_data)) {
  #   stop("omic_data should be a list.")
  # }

  # Check if both geno_data and omic_data are NULL or empty
  if (is.null(geno_data) && length(omic_data) == 0) {
    stop("Both geno_data and omic_data are NULL or empty.")
  }

  # Initialize merged_data
  merged_data <- NULL
  omic_count <- length(omic_data)

  # Merge geno_data if not NULL
  if (!is.null(geno_data)) {
    if(is.list(geno_data)){
      geno_data <- as.matrix(unlist(geno_data))
    }
    merged_data <- geno_data
    omic_count <- omic_count + 1  # Include geno_data in the count
  }

  # Merge omic_data
  if (length(omic_data) > 0) {
    merged_data <- if (is.null(merged_data)) {
      do.call(cbind, omic_data)
    } else {
      do.call(cbind, c(list(merged_data), omic_data))
    }
  }

  # Return merged data and omic_count if more than one source
  if (omic_count > 1) {
    return(list(merge_data = merged_data, omic_count = omic_count))
  } else {
    return(list(merge_data = merged_data))
  }
}


# merge_data <- function(geno_data,
#                        omic_data) {
#   # Check if both geno_data and omic_data are NULL
#   if (is.null(geno_data) && length(omic_data) == 0) { ## omic_data is a list of omics
#     stop("Both geno_data and omic_data are NULL.")
#   }
#
#   # If geno_data is NULL, merge omic_data only
#   if (is.null(geno_data)) {
#     if (length(omic_data)> 1){
#       omic_count <- length(omic_data)
#       return(list(merge_data = do.call(cbind, omic_data),
#                   omic_count = omic_count))
#     }
#     return(list(merge_data = do.call(cbind, omic_data)))
#   }
#
#   # If omic_data is empty, return geno_data
#   if (length(omic_data) == 0) {
#     return(list(merge_data = geno_data))
#     }
#
#   # Merge geno_data and omic_data
#   merged_data <- geno_data
#
#   for (i in seq_along(omic_data)) {
#     omic_count <- (length(omic_data)+1)
#     merged_data <- cbind(merged_data, omic_data[[i]])
#   }
#
#   return(list(merge_data = do.call(cbind, omic_data),
#               omic_count = omic_count))
# }

#' Title
#'
#' @param pheno_clean
#' @param response
#' @param geno_clean
#' @param omic_clean
#'
#' @return
#' @export
#'
#' @examples
ML_data_processing <- function(pheno_clean = NULL,
                               response = NULL,
                               geno_clean = NULL,
                               gen_name = NULL,
                               omic_clean = list()

) {
  # Separate train and test sets
  pheno_clean <- ML_undefined_test_train(object_pheno = pheno_clean,
                                         response = response)

  # Unpack pheno_clean if necessary
  if (length(pheno_clean) == 1 && is.list(pheno_clean[["pheno_clean_data"]])) {
    pheno_clean <- pheno_clean[["pheno_clean_data"]]
    test_set <- NULL
  } else{
    if (length(pheno_clean) == 2 && is.list(pheno_clean[["pheno_clean_data"]])) {
      pheno_clean <- pheno_clean[["pheno_clean_data"]]
      test_set <-pheno_clean[["test_set"]]
    }
  }


  # Merge geno and omic data
  merged_data <- merge_data(geno_clean, omic_clean)

  # Separate merged data into train and test sets
  if (!is.null(test_set)) {
    merged_data_test <- merged_data[rownames(merged_data[["merge_data"]]) %in% test_set[, gen_name], ]
    merged_data <- merged_data[!rownames(merged_data[["merge_data"]]) %in% test_set[, gen_name], ]
  }

  # Construct result list
  result <- list(pheno_clean_data = pheno_clean, merged_data = merged_data)
  if (length(merged_data)>1 && !is.null(merged_data[["omic_count"]])){
    result[["omic_count"]] <- merged_data[["omic_count"]]
  }

  # Add test set to result if it exists
  if (!is.null(test_set)) {
    result[["test_set"]] <- test_set
    result[["merged_data_test"]] <- merged_data_test
  } else {
    result[["test_set"]] <- NULL
    result[["merged_data_test"]] <- NULL
  }

  return(result)
}
