
# Validate terms in fixed and random effects are present in pheno_data
validate_terms <- function(term, data, term_type, gen_name, pheno_data) {
  msg <- "\n==================================================\n"

  if (length(pheno_data[[gen_name]]) == length(unique(pheno_data[[gen_name]]))){
    # Convert the formula to a character string
    formula_str <- as.character(term)

    # Check for the presence of a colon (:) in the formula string when the pheno data is a single location data
    contains_colon <- grepl(":", formula_str)

    if (any(contains_colon)) {
      stop(paste(msg, paste("The", term_type, "effect contains interaction term but your phenotypic data is a single environment.")), call. = FALSE)

    }
  }

  term_vars <- all.vars(term)
  if (!all(term_vars %in% names(data))) {
    stop(paste(msg, paste("All variables indicated in argument ", term_type, " should be present in phenotypic data.")), call. = FALSE)
  }
  #if(length(data[[gen_name]]) > length(unique(data[[gen_name]]))){
    missing_factors <- term_vars[!sapply(data[term_vars], is.factor)]
    data[missing_factors] <- lapply(data[missing_factors], factor)

  #}
    return(data)
}

#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param response
#' @param heter_groups
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
phenotype_precheck <- function(pheno_data = NULL,
                               gen_name = NULL,
                               response = NULL,
                               heter_groups = NULL,
                               random = NULL,
                               fixed = NULL,
                               type_pheno = NULL,
                               ...) {

  msg <- "\n==================================================\n"

  if (dplyr::is_grouped_df(pheno_data)) {
    pheno_data <- dplyr::ungroup(pheno_data)
  }

  # Check for empty data
  if (nrow(pheno_data) == 0) {
    stop(paste(msg, 'No pheno_data records provided.'), call. = FALSE)
  }

  # Ensure pheno_data is a data frame
  if (!inherits(pheno_data, 'data.frame')) {
    if (isTRUE(message)) warning(paste(msg, paste("pheno_data is not of class data.frame."," Converting it to a data frame.")))
    pheno_data <- as.data.frame(pheno_data)
  }

  # Check for presence of response variables
  missing_responses <- response[!response %in% colnames(pheno_data)]
  if (length(missing_responses) > 0) {
    stop(paste(msg, paste("The specified response variable(s) '", paste(missing_responses, collapse = "', '"), "' did not match with your data. Please check and use appropriately.")), call. = FALSE)
  }

  if(any(colnames(pheno_data) %in% c("NA", "Na", "na"))){
     stop(paste(msg, "Column names can't contain NA."))
  }
  # Order by heter_groups if specified and present
  if (!is.null(heter_groups)) {
    missing_heter_grps <- heter_groups[!heter_groups %in% colnames(pheno_data)]
    if (length(missing_heter_grps) > 0) {
      stop(paste(msg, paste("The variable '", paste(missing_heter_grps, collapse = "', '"), "' did not match with your data. Please check and use appropriately.")))
    } else {
      pheno_data <- pheno_data[order(pheno_data[[heter_groups]]), ]

      # Check the number of genotypes in each environment
      genotype_count <- pheno_data |>
        dplyr::group_by(!!rlang::sym(heter_groups)) |>
        dplyr::summarize(count = dplyr::n_distinct(!!rlang::sym(gen_name)))
      if(length(unique(genotype_count$count)) == 1) {
        stop(paste(msg, "Genotype counts are NOT identical across all environments."), call. = FALSE)
      }

      # Check if the genotype names are identical across all environments
      genotype_names_by_env <- pheno_data |>
        dplyr::group_by(!!rlang::sym(heter_groups)) |>
        dplyr::summarize(gen_names_list = list(sort(trimws(as.character(unique(!!rlang::sym(gen_name)))))))

      # Compare genotype names across environments element by element
      # Using base R to achieve the same functionality
      all_identical <- all(sapply(seq_along(genotype_names_by_env$gen_names_list[-1]), function(i) {
        identical(genotype_names_by_env$gen_names_list[[i]], genotype_names_by_env$gen_names_list[[1]])
      }))


      if(!all_identical) {
        stop(paste(msg,"Genotype names are NOT identical across all environments."), call. = FALSE)
      }

    }
  }

  # Check for gen_name presence
  if (!gen_name %in% colnames(pheno_data)) {
    stop(paste(msg, paste("The specified column", gen_name,"in the pheno_data did not match with your data. Please check and use appropriately.")), call. = FALSE)
  }

  # Check for NA in gen_name column
  if (anyNA(pheno_data[[gen_name]]) || any(pheno_data[[gen_name]] == -999)) {
    stop(paste(msg, paste("Column '", gen_name, "' should not have NA/missing values.")), call. = FALSE)
  }

  # Ensure response variables are numeric
  non_numeric_responses <- response[!sapply(pheno_data[response], is.numeric)]
  pheno_data[non_numeric_responses] <- lapply(pheno_data[non_numeric_responses], function(x) as.numeric(as.character(x)))

  # Check for zero variance in response variables
if(!is.null(type_pheno)){
  if(type_pheno !="test_set"){
  zero_variance_responses <- response[sapply(pheno_data[response], function(x) var(x, na.rm = TRUE) == 0)]
  if (length(zero_variance_responses) > 0) {
    msg <- "The following variable(s) have zero variance and cannot be used for prediction model: "
    stop(paste(msg, paste(paste(zero_variance_responses, collapse=", "), ". Check the raw data and model that generate the estimates.")), call. = FALSE)
  }

  }

}
  ## check that the fixed and random term are specified correctly

  if (!is.null(fixed)) {

    pheno_data <- validate_terms(fixed, pheno_data, "fixed", gen_name, pheno_data)
  }

  if (!is.null(random)) {
    pheno_data <- validate_terms(random, pheno_data, "random", gen_name, pheno_data)
  }

  attr(pheno_data, "cleared") <- "pass"
  return(pheno_data)
}
