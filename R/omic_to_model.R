#' Title
#'
#' @param omic_data
#' @param train_omic_data
#' @param test_omic_data
#' @param message
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
omic_to_model <- function(omic_data = NULL,
                          train_omic_data = NULL,
                          test_omic_data = NULL,
                          message = TRUE,
                          impute_omic = FALSE,
                          ...) {

  msg <- "\n==================================================\n"

  if (is.null(omic_data) & is.null(train_omic_data) & is.null(test_omic_data)) {
    stop(paste(msg, 'Omic data is missing.'), call. = FALSE)
  }
  # Check if omic_data is provided
  if (!is.null(omic_data) & is.null(train_omic_data) & is.null(test_omic_data)) {
    omic_object <- omic_precheck(object = omic_data, message = message,
                                 impute = impute_omic)

    # Check if omic_object passed the checks and is of the correct class
    if (attr(omic_object, "cleared") != "pass" || !inherits(omic_object, c("matrix", "array"))) {
      stop(print(paste(msg, 'Data is not of class Omic_matrix')), call. = FALSE)
    }

  } else {
    # Check if train_omic_data is provided
    if (!is.null(train_omic_data) & is.null(test_omic_data)) {
      omic_object <- omic_precheck(object = train_omic_data, message = message)

      if(is.null(omic_object)){
        stop(paste(msg, 'Train_omic_data did not pass the required test. Check the data.'), call. = FALSE)
      }
    }

    # Check if test_omic_data is provided
    if (!is.null(test_omic_data) & is.null(train_omic_data)) {
      stop(paste(msg, 'Training set is missing.'), call. = FALSE)
      # test_omic_data <- omic_precheck(object = test_omic_data, message = message)
      #
      # if(is.null(test_omic_data)){
      #   stop(paste(msg, 'Test_omic_data did not pass the required test. Check the data.'), call. = FALSE)
      # }

    }
    if (!is.null(train_omic_data) && !is.null(test_omic_data)) {

      if (!identical(colnames(train_geno_data), colnames(test_omic_data))) {

        stop(paste(msg, 'Omic data not match in training and testing set data.'), call. = FALSE)
      }

      omic_object <- rbind(train_omic_data, test_omic_data)
      omic_object <- omic_precheck(object = omic_object, message = message)

      if(is.null(omic_object)){
        stop(paste(msg, 'Omic data did not pass the required test. Check the data.'), call. = FALSE)
      }


    }

  } ### Ends

  # Set attribute to indicate readiness for model fit
  attr(omic_object, "cleared") <- "for_model_fit"

  return(omic_object)
}
