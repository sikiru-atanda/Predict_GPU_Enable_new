# Function to add missing hybrids to phenotypic data dynamically
add_extra_gid_from_geno_omic_to_pheno <- function(geno_data, pheno_data, gen_name, response_var, heter_group = NULL) {
  # Identify genotypes in the genotypic data but not in the phenotypic data
  diff_gid <- setdiff(rownames(geno_data), unique(pheno_data[[gen_name]]))

  # If no missing genotypes, return the original phenotypic data
  if (length(diff_gid) == 0) {
    return(pheno_data)
  }

  # Dynamically determine column names of the phenotypic data
  col_names <- colnames(pheno_data)

  # Create a base data frame for the missing genotypes
  if (!is.null(heter_group) && heter_group %in% col_names) {
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
