# Function to check row and column names consistency across datasets
check_names_consistency <- function(datasets, dataset_names) {
  # Filter out NULL datasets
  datasets_index <- which(!sapply(datasets, is.null))

  # If all datasets are NULL, return FALSE and exit
  if (length(datasets_index) == 0) {
    #message("All datasets are NULL. Nothing to check.")
    return(TRUE)
  }

  # If only one dataset is non-NULL, return TRUE
  if (length(datasets_index) == 1) {
    #message("Only one dataset is available. Nothing to check.")
    return(TRUE)
  }

  # Subset the datasets and names to exclude NULL entries
  datasets <- datasets[datasets_index]
  dataset_names <- dataset_names[datasets_index]

  # Check if all datasets are square matrices (n x n)
  is_square_matrix <- all(sapply(datasets, function(x) is.matrix(x) && nrow(x) == ncol(x)))

  # Define reference row and column names
  reference_row_names <- rownames(datasets[[1]])
  reference_col_names <- if (is_square_matrix) colnames(datasets[[1]]) else NULL

  # Check for consistency using `intersect`
  inconsistent_datasets <- sapply(seq_along(datasets), function(i) {
    row_mismatch <- length(intersect(reference_row_names, rownames(datasets[[i]]))) != length(reference_row_names)
    col_mismatch <- if (!is.null(reference_col_names)) {
      length(intersect(reference_col_names, colnames(datasets[[i]]))) != length(reference_col_names)
    } else {
      FALSE
    }
    row_mismatch || col_mismatch
  })

  # If inconsistencies are found, stop and report
  if (any(inconsistent_datasets)) {
    mismatched <- dataset_names[inconsistent_datasets]
    issues <- sapply(seq_along(datasets), function(i) {
      paste0(
        dataset_names[i], ": ",
        if (length(intersect(reference_row_names, rownames(datasets[[i]]))) != length(reference_row_names)) "Row names mismatch. " else "",
        if (!is.null(reference_col_names) && length(intersect(reference_col_names, colnames(datasets[[i]]))) != length(reference_col_names)) "Column names mismatch." else ""
      )
    })
    return(FALSE)
  }

  #message("All datasets have consistent row and column names (order-agnostic).")
  return(TRUE)
}



#check_names_consistency(datasets, dataset_names)
