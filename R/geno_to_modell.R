#' Prepare Genomic Data for Model Fitting
#'
#' This function prepares genomic data for model fitting by performing quality control (QC) checks, adjusting SNP coding, and ensuring compatibility between training and testing datasets. It can handle scenarios with only training data, only testing data, or both, and applies QC filters based on parameters such as minor allele frequency (MAF) and heterozygosity thresholds.
#'
#' @param geno_data A matrix or data.frame containing genomic data for the full dataset. Used if no separate training or testing datasets are provided.
#' @param train_geno_data A matrix or data.frame containing genomic data for the training set.
#' @param test_geno_data A matrix or data.frame containing genomic data for the testing set.
#' @param qc_filtering Logical, if TRUE, quality control filtering is applied.
#' @param maf_threshold Numeric, threshold for minor allele frequency below which SNPs are removed.
#' @param het_threshold Numeric, threshold for heterozygosity above which SNPs are removed.
#' @param ind_call_rate_threshold Numeric, threshold for individual call rate below which individuals are removed.
#' @param snp_call_rate_threshold Numeric, threshold for SNP call rate below which SNPs are removed.
#' @param impute Logical, if TRUE, indicates that missing values should be imputed. This feature is planned but not yet implemented.
#' @param map_data A data.frame or matrix containing marker information. Used alongside `geno_data` if provided.
#' @param message Logical, if TRUE, messages about the QC and data preparation process are displayed.
#' @param ... Additional arguments affecting the QC process.
#'
#' @return A list containing:
#'   - \code{snps_matrix}: The genomic data matrix after applying QC filters and ensuring compatibility between datasets.
#'   - \code{qc_metrics_and_summary_stat}: A data frame summarizing the QC process, including the number of markers and individuals removed.
#'
#' @examples
#' # Assuming `genomic_data` is a matrix with SNP data for the full dataset
#' result <- geno_to_model(geno_data = genomic_data, qc_filtering = TRUE)
#' prepared_genomic_data <- result$snps_matrix
#' qc_summary <- result$qc_metrics_and_summary_stat
#'
#' @importFrom stats colMeans
#' @importFrom dplyr rename rownames_to_column
#' @import tibble
#' @export


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
                                 #map_data  = map_data,
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
