
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

  msg <- "\n==================================================\n"
  # Check if both geno_data and omic_data are NULL
  if (is.null(geno_data) && length(omic_data) == 0) {
    stop(paste(msg, "Both geno_data and omic_data are NULL."),call. = FALSE)
  }

  # Function to check row names
  check_row_names <- function(data_list) {
    row_names <- lapply(data_list, rownames)
    all(sapply(row_names, function(x) all(x == row_names[[1]])))
  }

  # If geno_data is NULL, merge omic_data only after checking row names
  if (is.null(geno_data)) {
    if (length(omic_data) > 1) {
      if (!check_row_names(omic_data)) {
        stop(paste(msg, "Row names of omic_data do not match."), call. = FALSE)
      }
      omic_count <- length(omic_data)
      return(list(merge_data = do.call(cbind, omic_data), omic_count = omic_count))
    }
    return(list(merge_data = do.call(cbind, omic_data), omic_count = 1))
  }

  # If omic_data is empty, return geno_data
  if (length(omic_data) == 0) {
    return(list(merge_data = geno_data, omic_count = 1))
  }

  # Check if row names of geno_data match with each omic_data
  for (i in seq_along(omic_data)) {
    if (!all(rownames(geno_data) == rownames(omic_data[[i]]))) {
      stop(paste(msg,paste("Row names of geno_data do not match with omic_data at index", i)), call. = FALSE)
    }
  }

  # Merge geno_data and omic_data
  merged_data <- geno_data
  for (i in seq_along(omic_data)) {
    merged_data <- cbind(merged_data, omic_data[[i]])
  }
  omic_count <- length(omic_data) + 1

  return(list(merge_data = merged_data, omic_count = omic_count))
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
#     return(list(merge_data = do.call(cbind, omic_data),
#                 omic_count = 1))
#   }
#
#   # If omic_data is empty, return geno_data
#   if (length(omic_data) == 0) {
#     return(list(merge_data = geno_data,
#                 omic_count = 1))
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

  msg <- "\n==================================================\n"
  test_set <- NULL
  # Separate train and test sets
  # pheno_clean <- ML_undefined_test_train(object_pheno = pheno_clean,
  #                                        response = response)
  pheno_data <- pheno_clean[["pheno_clean_data"]]
  # Unpack pheno_clean if necessary
  if("test_set"%in%names(pheno_clean)) {

    test_set <- pheno_clean[["test_set"]]
    if(is.data.frame(test_set) | is.matrix(test_set)){

      test_set <-  as.character(test_set[, 1])


    } else if (!is.list(test_set)){

      test_set <-  as.character(test_set)

    } else {
      if(is.list(test_set)){
        stop(paste(msg,'The testing set cannot be a list. Should be either dataframe, matrix or a vector.'), call. = FALSE)

      }

    }

    pheno_data <-  pheno_data[!pheno_data[[gen_name]] %in% test_set, ]

  } else {
    # Find the rows with NA in each response column if the user has NA as testing set
    na_rows <- lapply(response, function(col) which(is.na(pheno_data[[col]])))

    # Check if all vectors of NA rows are identical across the response, sparse is not allowed
    if (!all(sapply(na_rows, function(x) identical(x, na_rows[[1]])))) {
      stop(paste(msg, sprintf("Rows containing NA did not match across the response columns: %s.", paste(response, collapse = ", "))), call. = FALSE)
    }

    if (length(na_rows[[1]]) > 0) {

      test_set <- unique(na_rows[[1]])

      test_set <- as.character(pheno_data[[gen_name]][test_set])

      pheno_data <-  pheno_data[!pheno_data[[gen_name]] %in% test_set, ]

    }

  }

  # Merge geno and omic data
  merged_data <- merge_data(geno_clean, omic_clean)

  # Separate merged data into train and test sets
  if (!is.null(test_set)) {
    merged_data_test <- merged_data[["merge_data"]][rownames(merged_data[["merge_data"]]) %in% test_set, ]
    merged_data[["merge_data"]] <- merged_data[["merge_data"]][!rownames(merged_data[["merge_data"]]) %in% test_set, ]
  }

  # Construct result list
  result <- list(pheno_clean_data = pheno_data,
                 merged_data = merged_data)

  if (length(merged_data)>1 && !is.null(merged_data[["omic_count"]])){
    result[["omic_count"]] <- merged_data[["omic_count"]]
  }

  # Add test set to result if it exists
  if (!is.null(test_set)) {
    result[["test_set"]] <- test_set
    result[["merged_data_test"]] <- merged_data_test
  }

  return(result)
}


# merge_data <- function(geno_data, omic_data) {
#   # Ensure omic_data is a list for consistency
#   # if (!is.list(omic_data)) {
#   #   stop("omic_data should be a list.")
#   # }
#
#   # Check if both geno_data and omic_data are NULL or empty
#   # if (is.null(geno_data) && length(omic_data) == 0) {
#   #   stop("Both geno_data and omic_data are NULL or empty.")
#   # }
#
#   # Initialize merged_data
#   merged_data <- NULL
#   omic_count <- length(omic_data)
#
#   # Merge geno_data if not NULL
#   if (!is.null(geno_data)) {
#     if(is.list(geno_data)){
#       geno_data <- as.matrix(unlist(geno_data))
#     }
#     merged_data <- geno_data
#     omic_count <- omic_count + 1  # Include geno_data in the count
#   }
#
#   # Merge omic_data
#   if (length(omic_data) > 0) {
#     merged_data <- if (is.null(merged_data)) {
#       do.call(cbind, omic_data)
#     } else {
#       do.call(cbind, c(list(merged_data), omic_data))
#     }
#   }
#
#   # Return merged data and omic_count if more than one source
#   if (omic_count > 1) {
#     return(list(merge_data = merged_data, omic_count = omic_count))
#   } else {
#     return(list(merge_data = merged_data))
#   }
# }
