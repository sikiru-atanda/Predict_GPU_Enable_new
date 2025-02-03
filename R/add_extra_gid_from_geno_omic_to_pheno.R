# Function to add missing hybrids to phenotypic data dynamically
add_extra_gid_from_geno_omic_to_pheno <- function(geno_data, pheno_data, gen_name, response_var, heter_group = NULL) {
  # Identify genotypes in the genotypic data but not in the phenotypic data
  msg <- "\n==================================================\n"

  diff_gid <- setdiff(rownames(geno_data), unique(pheno_data[[gen_name]]))
#   # Ensure response_var exist in df
#   if (!all(response_var %in% colnames(pheno_data))) {
#     stop(paste(msg, "Some response variables do not exist in the dataframe."), call. = FALSE)
#   }
#
#   # Create a logical matrix indicating NA positions for response variables
#   na_matrix <- is.na(pheno_data[response_var])
#
#   # Check if all rows have the same NA pattern
#   if (!all(rowSums(na_matrix) %in% c(0, length(response_var)))) {
#     stop(paste(msg, "Not all response variable columns have NAs in the same positions."), call. = FALSE)
#   }
#   if (all(na_matrix == FALSE)){
#   diff_gid <- setdiff(rownames(geno_data), unique(pheno_data[[gen_name]]))
#   } else{
#     diff_gid <- NULL
#     return(TRUE)
# }
  # If no missing genotypes, return the original phenotypic data
  if (length(diff_gid) == 0) {
    return(pheno_data)
  }

  # Dynamically determine column names of the phenotypic data
  col_names <- colnames(pheno_data)

  # Create a base data frame for the missing genotypes
  if (!is.null(heter_group)) {

    if(heter_group %in% col_names){
    # For multi-environment data: Expand genotypes to all environments
    environments <- as.character(unique(pheno_data[[heter_group]]))
    diff_data <- expand.grid(
      gen_name = diff_gid,
      heter_group = environments,
      stringsAsFactors = FALSE
    )

    colnames(diff_data) <- c(gen_name, heter_group)
    # Add remaining columns dynamically and fill them with NA
    for (col in setdiff(col_names, names(diff_data))) {
      diff_data[[col]] <- NA
    }

    # Ensure column order matches the original phenotypic data
    diff_data <- diff_data[, col_names, drop = FALSE]

    }
  } else {
    # For single-environment data: Create a data frame with all columns as NA
    diff_data <- data.frame(matrix(NA, nrow = length(diff_gid), ncol = length(col_names)))
    colnames(diff_data) <- col_names

    # Populate the `gen_name` column with missing genotypes
    diff_data[[gen_name]] <- diff_gid

    # Populate the `response_var` column with NA for the missing genotypes
    if (response_var %in% col_names) {
      diff_data[[response_var]] <- NA
    }
  }

  # Combine the phenotypic data with the new rows
  pheno_data <- rbind(pheno_data, diff_data)

  return(pheno_data)
}
