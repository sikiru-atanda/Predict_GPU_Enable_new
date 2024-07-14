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
omic_to_modelOLD <- function(omic_data = NULL,
                          train_omic_data = NULL,
                          test_omic_data = NULL,
                          message = TRUE,
                          ...) {

  msg <- "\n==================================================\n"

  # Placeholder for omic_precheck function
  # omic_precheck <- function(object, message) {
  #   # Implement omic precheck logic
  #   # ...
  #
  #   return(omic_object)
  # }

  if (!is.null(omic_data)) {
    omic_object <- omic_precheck(object = omic_data, message = message)
    stopifnot(attr(omic_object, "cleared") == "pass" &&
                all(class(omic_object) == c("matrix", "array", "omic_matrix")),
              msg, 'Data is not class Omic_matrix1')
  } else {
    if (!is.null(train_omic_data)) {
      train_omic_data <- omic_precheck(object = train_omic_data, message = message)
      stopifnot(attr(train_omic_data, "cleared") == "pass" &&
                  all(class(train_omic_data) == c("matrix", "array", "omic_matrix")),
                msg, 'Data is not class Omic_matrix4')
    }
    if (!is.null(test_omic_data)) {
      test_omic_data <- omic_precheck(object = test_omic_data, message = message)
      stopifnot(attr(test_omic_data, "cleared") == "pass" &&
                  all(class(test_omic_data) == c("matrix", "array", "omic_matrix")),
                msg, 'Data is not class Omic_matrix5')
    }

    if (!is.null(train_omic_data) && !is.null(test_omic_data)) {
      if (all(class(test_omic_data) == c("matrix", "array", "omic_matrix")) &&
          attr(train_omic_data, "cleared") == "pass" &&
          all(class(train_omic_data) == c("matrix", "array", "omic_matrix"))) {
        stopifnot(identical(colnames(train_omic_data), colnames(test_omic_data)),
                  msg, 'Data is not class Omic_matrix2.')
        omic_object <- rbind(train_omic_data, test_omic_data)
      } else {
        stopifnot(FALSE, msg, 'Data is not object Omic_matrix3')
      }
    } else if (!is.null(train_omic_data) && is.null(test_omic_data)) {
      message(paste(msg, 'Only train_omic_data is provided'))
      stopifnot(attr(train_omic_data, "cleared") == "pass" &&
                  all(class(train_omic_data) == c("matrix", "array", "omic_matrix")),
                msg, 'Data is not class Omic_matrix4')
      omic_object <- train_omic_data
    } else if (is.null(train_omic_data) && !is.null(test_omic_data)) {
      if (message) {
        message(paste(msg, 'Only test_omic_data is provided'))
      }
      stopifnot(attr(test_omic_data, "cleared") == "pass" &&
                  all(class(test_omic_data) == c("matrix", "array", "omic_matrix")),
                msg, 'Data is not class Omic_matrix5')
      omic_object <- test_omic_data
    }
  }

  class(omic_object) <- c("matrix", "array", "omic_matrix")
  attr(omic_object, "cleared") <- "for_model_fit"

  return(omic_object)
}
