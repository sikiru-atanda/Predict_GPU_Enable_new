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

  # Check if omic_data is provided
  if (!is.null(omic_data)) {
    omic_object <- omic_precheck(object = omic_data, message = message,
                                 impute = impute_omic)

    # Check if omic_object passed the checks and is of the correct class
    if (attr(omic_object, "cleared") != "pass" || !inherits(omic_object, c("matrix", "array"))) {
      stop(print(paste(msg, 'Data is not of class Omic_matrix')), call. = FALSE)
    }

  } else {
    # Check if train_omic_data is provided
    if (!is.null(train_omic_data)) {
      train_omic_data <- omic_precheck(object = train_omic_data, message = message)
    }

    # Check if test_omic_data is provided
    if (!is.null(test_omic_data)) {
      test_omic_data <- omic_precheck(object = test_omic_data, message = message)
    }

    # Check if both train_omic_data and test_omic_data passed the checks and have matching columns
    if (!is.null(train_omic_data) && !is.null(test_omic_data)) {
        if(attr(train_omic_data, "cleared") == "pass" &&
        attr(test_omic_data, "cleared") == "pass"){
      if(identical(colnames(train_omic_data), colnames(test_omic_data))){
        omic_object <- rbind(train_omic_data, test_omic_data)
      } else {
        stop(print(paste(msg,'Colnames for train_omic_data and test_omic_data did not match.')), call. = FALSE)
      }

        } else {
          stop(print(paste(msg,'Data is not object Omic_matrix')), call. = FALSE)
        }
          # Remove train_omic_data and test_omic_data from memory
      rm(train_omic_data, test_omic_data)

    } else if (!is.null(train_omic_data) && is.null(test_omic_data)) {
      # Check if only train_omic_data is provided
      message(paste(msg, 'Only train_omic_data is provided'))

      # Check if train_omic_data passed the checks
      if (attr(train_omic_data, "cleared") == "pass" && inherits(train_omic_data, c("matrix", "array"))) {
        omic_object <- train_omic_data

        # Remove train_omic_data from memory
        rm(train_omic_data)
      } else {
        stop(print(paste(msg, 'Data is not of class Omic_matrix')), call. = FALSE)
      }
    } else if (is.null(train_omic_data) && !is.null(test_omic_data)) {
      # Check if only test_omic_data is provided
      if (isTRUE(message)) {
        message(paste(msg, 'Only test_omic_data is provided'))
      }

      # Check if test_omic_data passed the checks
      if (attr(test_omic_data, "cleared") == "pass" && inherits(test_omic_data, c("matrix", "array"))) {
        omic_object <- test_omic_data

        # Remove test_omic_data from memory
        rm(test_omic_data)
      } else {
        stop(print(paste(msg, 'Data is not of class Omic_matrix')), call. = FALSE)
      }
    } else {
      # None of the omic data is provided
      stop(print(paste(msg, 'No omic data provided')), call. = FALSE)
    }
  }

  # Set attribute to indicate readiness for model fit
  attr(omic_object, "cleared") <- "for_model_fit"

  return(omic_object)
}
