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
                               ...) {

  msg <- "==================================================\n"

  if (nrow(pheno_data) == 0) {
    stop(print(paste(msg, 'No pheno_data records provided.')), call. = FALSE)
  }

  if (!inherits(pheno_data, 'data.frame')) {
    message(print(paste(msg, "'pheno_data' is not of class 'data.frame'. Converting it to a data frame.")))
    pheno_data <- as.data.frame(pheno_data)
  }


# check if all the response variables are present in the phenotypi --------

  if (sum(response%in%colnames(pheno_data)) < length(response)) {
    stop(print(paste(msg, paste("The specified response variable(s) '", paste(response, collapse = "', '"), "' did not match with your data. Please check and use appropriately."))), call. = FALSE)
  }

  if (!is.null(heter_groups)) {
    if(heter_groups %in% colnames(pheno_data)){
    pheno_data <- pheno_data[order(pheno_data[, heter_groups]), ]

    } else {
      stop(print(paste(msg, paste("The specified heterogeneity groups '", heter_groups, "' did not match with your data. Please check and use appropriately."))), call. = FALSE)
    }
  }

# Check if user defined GID/Name match the column name --------------------

  if (!gen_name %in% colnames(pheno_data)) {
    stop(print(paste(msg, paste("The specified column '", gen_name, "' in the pheno_data did not match with your data. Please check and use appropriately."))), call. = FALSE)
  }

# Check for NA in the name/GID --------------------------------------------

  if (anyNA(pheno_data[, gen_name]) || any(pheno_data[, gen_name] == -999)) {
    stop(print(paste(msg, paste("Column '", gen_name, "' should not have NA/missing values."))), call. = FALSE)
  }


# Check if the response variable if double/numeric if not convert  --------
  if (!all(sapply(response, function(x, pheno_data) is.numeric(pheno_data[, x]), pheno_data))) {
    pheno_data[, response] <- lapply(pheno_data[, response, drop = FALSE], function(x) as.double(as.character(x)))
  }

# check for the variance of the user defined trait or traits --------------
  if (any(sapply(response, function(x) var(pheno_data[, x], na.rm = TRUE) == 0))) {
    zero_variance_vars <- response[sapply(response, function(x) var(pheno_data[, x], na.rm = TRUE) == 0)]
    stop(print(paste(msg, paste(zero_variance_vars, "variable(s) have zero variance. These traits cannot be used for prediction model. Check the raw data and model that generate the BLUEs."))), call. = FALSE)
  }


  attr(pheno_data, "cleared") <- "pass"

  return(pheno_data)
}
