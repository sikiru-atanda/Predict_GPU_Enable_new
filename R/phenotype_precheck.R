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
                               ...) {

  msg <- "==================================================\n"

  # Check for empty data
  if (nrow(pheno_data) == 0) {
    stop(msg, 'No pheno_data records provided.')
  }

  # Ensure pheno_data is a data frame
  if (!inherits(pheno_data, 'data.frame')) {
    if (message) warning(msg, "'pheno_data' is not of class 'data.frame'. Converting it to a data frame.")
    pheno_data <- as.data.frame(pheno_data)
  }

  # Check for presence of response variables
  missing_responses <- response[!response %in% colnames(pheno_data)]
  if (length(missing_responses) > 0) {
    stop(msg, "The specified response variable(s) '", paste(missing_responses, collapse = "', '"), "' did not match with your data. Please check and use appropriately.")
  }

  # Order by heter_groups if specified and present
  if (!is.null(heter_groups) && heter_groups %in% colnames(pheno_data)) {
    pheno_data <- pheno_data[order(pheno_data[[heter_groups]]), ]
  } else if (!is.null(heter_groups)) {
    stop(msg, "The specified heterogeneity groups '", heter_groups, "' did not match with your data. Please check and use appropriately.")
  }

  # Check for gen_name presence
  if (!gen_name %in% colnames(pheno_data)) {
    stop(msg, "The specified column '", gen_name, "' in the pheno_data did not match with your data. Please check and use appropriately.")
  }

  # Check for NA in gen_name column
  if (anyNA(pheno_data[[gen_name]]) || any(pheno_data[[gen_name]] == -999)) {
    stop(msg, "Column '", gen_name, "' should not have NA/missing values.")
  }

  # Ensure response variables are numeric
  non_numeric_responses <- response[!sapply(pheno_data[response], is.numeric)]
  pheno_data[non_numeric_responses] <- lapply(pheno_data[non_numeric_responses], function(x) as.numeric(as.character(x)))

  # Check for zero variance in response variables
  zero_variance_responses <- response[sapply(pheno_data[response], function(x) var(x, na.rm = TRUE) == 0)]
  if (length(zero_variance_responses) > 0) {
    stop(msg, paste(zero_variance_responses, "variable(s) have zero variance. These traits cannot be used for prediction model. Check the raw data and model that generate the BLUEs."))
  }

  # Validate terms in fixed and random effects are present in pheno_data
  validate_terms <- function(term, data, term_type, gen_name) {
    term_vars <- all.vars(term)
    if (!all(term_vars %in% names(data))) {
      stop(msg, "All variables indicated in argument ", term_type, " should be present in phenotypic data.")
    }
    if(length(data[[gen_name]]) > length(unique(data[[gen_name]]))){
      missing_factors <- term_vars[!sapply(data[term_vars], is.factor)]
      data[missing_factors] <- lapply(data[missing_factors], factor)

    }
  }

  if (!is.null(fixed)) {
    validate_terms(fixed, pheno_data, "fixed", gen_name)
  }

  if (!is.null(random)) {
    validate_terms(random, pheno_data, "random", gen_name)
  }

  attr(pheno_data, "cleared") <- "pass"
  return(pheno_data)
}
