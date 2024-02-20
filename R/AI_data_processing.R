
#' Title
#'
#' @param geno_data
#' @param omic_data
#'
#' @return
#' @export
#'
#' @examples
merge_data <- function(geno_data,
                       omic_data) {
  # Check if both geno_data and omic_data are NULL
  if (is.null(geno_data) && length(omic_data) == 0) {
    stop("Both geno_data and omic_data are NULL.")
  }

  # If geno_data is NULL, merge omic_data only
  if (is.null(geno_data)) {
    return(do.call(cbind, omic_data))
  }

  # If omic_data is empty, return geno_data
  if (length(omic_data) == 0) {
    return(geno_data)
  }

  # Merge geno_data and omic_data
  merged_data <- geno_data

  for (i in seq_along(omic_data)) {
    merged_data <- cbind(merged_data, omic_data[[i]])
  }

  return(merged_data)
}

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
    if (length(pheno_clean) == 2 && is.list(pheno_clean[["pheno_clean_data"]])) {
      pheno_clean <- pheno_clean[["pheno_clean_data"]]
      test_set_ <-pheno_clean[["test_set"]]
    }
  }

  # Merge geno and omic data
  merged_data <- merge_data(geno_clean, omic_clean)

  # Separate merged data into train and test sets
  if (exists("test_set_")) {
    merged_data_test <- merged_data[rownames(merged_data) %in% test_set_[, gen_name], ]
    merged_data <- merged_data[!rownames(merged_data) %in% test_set_[, gen_name], ]
  }

  # Construct result list
  result <- list(pheno_clean_data = pheno_clean, merged_data = merged_data)

  # Add test set to result if it exists
  if (exists("test_set_")) {
    result$test_set <- test_set_
    result$merged_data_test <- merged_data_test
  } else {
    result$test_set <- NULL
    result$merged_data_test <- NULL
  }

  return(result)
}
