#' Prepare omics data for model fitting (pre-check, impute, align train/test)
#'
#' Wraps [omic_precheck] for the model-fitting workflow: accepts either a single
#' omics matrix `omic_data` or pre-split `train_omic_data` / `test_omic_data`,
#' applies the missingness threshold, optionally imputes, and returns the
#' matrix / matrices ready to be passed to a model.
#'
#' @param omic_data Single omics matrix for the full set (used when no separate
#'   train / test split is supplied).
#' @param train_omic_data Optional training omics matrix.
#' @param test_omic_data Optional test-set omics matrix matching
#'   `train_omic_data`.
#' @param message Logical; if `TRUE`, print summary messages.
#' @param impute_omic Logical; if `TRUE`, impute missing values.
#' @param imputation_method One of `"knn"`, `"mean"`, `"median"`.
#' @param impute_knn_k Integer; number of nearest neighbours when
#'   `imputation_method = "knn"`.
#' @param na_threshold Numeric in `(0, 1]`; columns whose missing-value
#'   proportion exceeds this are dropped.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return The pre-checked omics matrix or list of (train, test) matrices.
#' @export
#'
#' @examples
omic_to_model <- function(omic_data = NULL,
                          train_omic_data = NULL,
                          test_omic_data = NULL,
                          message = TRUE,
                          impute_omic = FALSE,
                          imputation_method = "knn",
                          impute_knn_k = 5,
                          na_threshold = 0.9,
                          ...) {

  msg <- ""

  if (is.null(omic_data) & is.null(train_omic_data) & is.null(test_omic_data)) {
    stop(paste(msg, 'Omic data is missing.'), call. = FALSE)
  }
  # Check if omic_data is provided
  if (!is.null(omic_data) & is.null(train_omic_data) & is.null(test_omic_data)) {
    omic_object <- omic_precheck(object = omic_data,
                                 impute = impute_omic,
                                 imputation_method = imputation_method,
                                 impute_knn_k = impute_knn_k,
                                 na_threshold = na_threshold,
                                 message = message)

    # Check if omic_object passed the checks and is of the correct class
    if (attr(omic_object, "cleared") != "pass" || !inherits(omic_object, c("matrix", "array"))) {
      stop(print(paste(msg, 'Data is not of class Omic_matrix')), call. = FALSE)
    }

  } else {
    # Check if train_omic_data is provided
    if (!is.null(train_omic_data) & is.null(test_omic_data)) {
      omic_object <- omic_precheck(object = train_omic_data,
                                   impute = impute_omic,
                                   imputation_method = imputation_method,
                                   impute_knn_k = impute_knn_k,
                                   na_threshold = na_threshold,
                                   message = message)

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

      if (!identical(colnames(train_omic_data), colnames(test_omic_data))) {

        stop(paste(msg, 'Omic data not match in training and testing set data.'), call. = FALSE)
      }

      omic_object <- rbind(train_omic_data, test_omic_data)
      omic_object <- omic_precheck(object = omic_object,
                                   impute = impute_omic,
                                   imputation_method = imputation_method,
                                   impute_knn_k = impute_knn_k,
                                   na_threshold = na_threshold,
                                   message = message)

      if(is.null(omic_object)){
        stop(paste(msg, 'Omic data did not pass the required test. Check the data.'), call. = FALSE)
      }


    }

  } ### Ends

  # Set attribute to indicate readiness for model fit
  attr(omic_object, "cleared") <- "for_model_fit"

  return(omic_object)
}
