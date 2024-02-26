#' Title
#'
#' @param object_coeff1
#' @param object_coeff2
#' @param object_coeff3
#' @param object_coeff4
#' @param mu
#' @param object_geno
#' @param object_omics1
#' @param object_omics2
#' @param object_omics3
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict_from_coeff <- function(
    object_coeff1 = NULL,
    object_coeff2 = NULL,
    object_coeff3 = NULL,
    object_coeff4 = NULL,
    mu = NULL,
    object_geno = NULL,
    object_omics1 = NULL,
    object_omics2 = NULL,
    object_omics3 = NULL,
    message = TRUE,
    ...
) {
  # Create lists of coefficient objects and omics data objects
  coef_list <- list(object_coeff1, object_coeff2, object_coeff3, object_coeff4)
  omics_data_list <- list(object_geno, object_omics1, object_omics2, object_omics3)

  # Check if any coefficient object and omics data object are provided
  if (any(sapply(coef_list, Negate(is.null))) && any(sapply(omics_data_list, Negate(is.null)))) {

    # Get the index of the non-NULL coefficient object
    # Check if the corresponding omics data object is also provided
    coef_index <- which(!sapply(coef_list, is.null))
    omic_index <- which(!sapply(omics_data_list, is.null))
    if(!identical(coef_index, omic_index)){
      stop("Coefficients provided but corresponding omics data is missing.")
    }

    # Find the index of the non-NULL coefficient object for first
    coef_index <- which(!is.null(coef_list))[1]

    # Compute prediction based on the provided coefficient and omics data
    pred <- omics_data_list[[coef_index]] %*% coef_list[[coef_index]]

    # Sum up predictions from multiple omics data objects if available
    if(length(omics_data_list)>1){
      for (i in seq_along(omics_data_list)) {
        ## this ensure it not starting from 1 again
        if (i != coef_index && omics_data_list[[i]]) {
          pred <- pred + omics_data_list[[i]] %*% coef_list[[i]]
        }
      }

    }
    # Add intercept term if provided
    if (!is.null(mu)) {
      pred <- pred + mu
    }

    return(pred)

  } else {
    stop("Provide both coefficients and omics data.")
  }
}
