#' Title
#' #' Title
#'
#'Overall, this serve as gateway between snp/marker-precheck function and readiness of
#' the snp/marker data for model fitting#'
#' The objective of this function is to do the following:
#' 1. Check the output from geno-precheck function for geno_data before declaring it for model fit
#' 2. If geno_data_train and geno_data_test were present and pass through the pre-check process,
#'   These will processed be processed that is:
#'    1) It check that column name (maker/snp) for both data match/the same
#'    2) Combined the dataset for model fit and prediction.
#'    It is assumed here that geno_data is missing/not provided by the user.
#' 3. The output will be a matrix(geno_data) declared for model fit.
#'
#'
#' @param geno_data
#' @param train_geno_data
#' @param test_geno_data
#' @param qc_filtering
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param map_data
#' @param message
#' @param ...


geno_to_model <- function(geno_data = NULL,
                          train_geno_data = NULL,
                          test_geno_data = NULL,
                          qc_filtering = TRUE,
                          maf_threshold = 0.01,
                          het_threshold = 0.2,
                          ind_call_rate_threshold = 0.9,
                          snp_call_rate_threshold = 0.9,
                          impute = TRUE,
                          map_data = NULL,
                          message = TRUE,
                          ...) {
  msg <- sprintf("==================================================\n")

  check_and_prepare <- function(data,...) {
    geno_object <- geno_precheck(object_geno = data,
                                 qc_filtering = qc_filtering,
                                 maf_threshold = maf_threshold,
                                 het_threshold = het_threshold,
                                 ind_call_rate_threshold = ind_call_rate_threshold,
                                 snp_call_rate_threshold = snp_call_rate_threshold,
                                 impute = impute,
                                 map_data  = map_data,
                                 message = message)



    if (attr(geno_object[[1]], "cleared") == "pass" && all(class(geno_object[[1]]) == c("matrix", "array"))) {
      #class(geno_object[[1]]) <- c("matrix", "array", "geno_data")
      attr(geno_object[[1]], "cleared") <- "for_model_fit"
    } else {
      geno_object <- NULL
    }

    return(geno_object)
  }

  if (!is.null(geno_data) && is.null(map_data) && is.null(test_geno_data) && is.null(train_geno_data)) {
    geno_object <- check_and_prepare(data = geno_data)

  } else if (!is.null(geno_data) && !is.null(map_data) && is.null(test_geno_data) && is.null(train_geno_data)) {
    geno_object <- check_and_prepare(data = geno_data)


    if (attr(geno_object[[1]], "cleared") != "pass" || !identical(class(geno_object[[1]]), c("matrix", "array"))) {
      stop(print(paste(msg, 'Geno or the omic data did not pass the required test. Check the data.')), call. = FALSE)
    }

  } else {
    if (!is.null(train_geno_data)) {
      train_geno_data <- check_and_prepare(data = train_geno_data)

    }

    if (!is.null(test_geno_data)) {
      test_geno_data <- check_and_prepare(data = test_geno_data)

    }

    if (!is.null(train_geno_data) && !is.null(test_geno_data)) {
      if (((attr(test_geno_data[[1]], "cleared") == "pass" && all(class(test_geno_data[[1]]) == c("matrix", "array"))) &&
           (attr(train_geno_data[[1]], "cleared") == "pass" && all(class(train_geno_data[[1]]) == c("matrix", "array"))))) {

        if (!identical(colnames(train_geno_data[[1]]), colnames(test_geno_data[[1]]))) {
          if (dim(train_geno_data[[1]])[2] != 0 && dim(test_geno_data[[1]])[2] != 0) {
            snp_names <- intersect(colnames(train_geno_data[[1]]), colnames(test_geno_data[[1]]))
            train_geno_data[[1]] <- train_geno_data[[1]][, colnames(train_geno_data[[1]]) %in% snp_names]
            test_geno_data[[1]] <- test_geno_data[[1]][, colnames(test_geno_data[[1]]) %in% snp_names]

            geno_object <- list(snps_matrix = rbind(train_geno_data[[1]], test_geno_data[[1]]))
            #class(geno_object) <- c("matrix", "array", "geno_data")
            attr(geno_object, "cleared") <- "for_model_fit"
            geno_object$qc_metrics_and_summary_stat <- train_geno_data[[2]]
          } else {
            stop(print(paste(msg, 'SNP/markers did not match in training and testing set data.')), call. = FALSE)
          }
        } else {
          geno_object <- rbind(train_geno_data[[1]], test_geno_data[[1]])
          #class(geno_object) <- c("matrix", "array", "geno_data")
          attr(geno_object, "cleared") <- "for_model_fit"
          geno_object <- list(snps_matrix = geno_object,
                              qc_metrics_and_summary_stat = train_geno_data[[2]])
        }
      } else {
        stop(print(paste(msg, 'SNP/marker data did not pass the required test. Check the data.')), call. = FALSE)
      }

    } else if (!is.null(train_geno_data) && is.null(test_geno_data)) {
      if (isTRUE(message)) {
        message(insight::print_color(paste(msg, paste("Only train_geno_data is provided.")), "blue"))
      }
      if ((attr(train_geno_data[[1]], "cleared") == "pass" && all(class(train_geno_data[[1]]) == c("matrix", "array")))) {
        geno_object <- train_geno_data[[1]]
        #class(geno_object) <- c("matrix", "array", "geno_data")
        attr(geno_object, "cleared") <- "for_model_fit"
        geno_object <- list(snps_matrix = geno_object,
                            qc_metrics_and_summary_stat = train_geno_data[[2]])
        rm(train_geno_data)
      }

    } else if (is.null(train_geno_data) && !is.null(test_geno_data)) {
      if (isTRUE(message)) {
        message(paste(insight::print_color("WARNINGS\n", "blue"),
                      insight::print_color(paste(msg, paste("Only the test_geno_data is provided.\n\t Check if this is correct.")), "blue")))
      }
      if ((attr(test_geno_data[[1]], "cleared") == "pass" && all(class(test_geno_data[[1]]) == c("matrix", "array")))) {
        geno_object <- test_geno_data[[1]]
        #class(geno_object) <- c("matrix", "array", "geno_data")
        attr(geno_object, "cleared") <- "for_model_fit"
        geno_object <- list(snps_matrix = geno_object,
                            qc_metrics_and_summary_stat = test_geno_data[[2]])
        rm(test_geno_data)
      }
    }
  }

  return(geno_object)
}
