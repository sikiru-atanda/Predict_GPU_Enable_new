
#' Title
#'
#' @param test_set
#' @param gen_name
#'
#' @return
#' @export
#'
#' @examples
check_test_set <- function(test_set = NULL
                           ) {
  msg <- sprintf("==================================================\n")
  if (is.null(test_set)) return(NULL)

  if (!is.data.frame(test_set) && !is.matrix(test_set) && !is.vector(test_set)) {
    stop(message(paste(msg, 'The testing set should be a dataframe, matrix, or a vector.')), call. = FALSE)
  }

  if (is.data.frame(test_set) || is.matrix(test_set)) {
    if(ncol(test_set)>1){
      stop(message(paste(msg, 'The testing set should be a dataframe, matrix with single column.')), call. = FALSE)
    }
    test_set <- unique(test_set[, 1])
  } else {
    test_set <- unique(test_set)
  }

  # if (length(test_set) == 0) {
  #   rm(test_set)
  #   message(insight::print_color(paste(msg, 'No test_set available.'), "red"))
  #   return(NULL)
  # }

  return(test_set)
}

pheno_geno_match <- function(object_geno = NULL,
                             object_pheno = NULL,
                             gen_name = NULL,
                             train_set = NULL,
                             test_set = NULL,
                             message = TRUE) {
  msg <- sprintf("==================================================\n")
  ID_pheno <- as.character(unique(object_pheno[, gen_name]))

  if (!isTRUE(all(ID_pheno %in% rownames(object_geno)))) {
    stop(print(paste(msg, 'Not all individuals with phenotypic records have genotypic/omic records.')), call. = FALSE)
  }

  if (!isTRUE(all(rownames(object_geno) %in% ID_pheno))) {
    if (isTRUE(message)) {
      message(insight::print_color(paste(msg, 'Not all individuals with genotypic/omic records have phenotypic records.'), "blue"))
    }

    test_set <- setdiff(rownames(object_geno), ID_pheno)
    test_set <- check_test_set(test_set)

    if (!is.null(test_set)) {

# check whether geno_omic data or grm/kernel ------------------------------

      if(nrow(object_geno)==ncol(object_geno)){
      object_geno <- object_geno[c(ID_pheno, test_set), c(ID_pheno, test_set)]
      } else {
        object_geno <- object_geno[c(ID_pheno, test_set), ]
      }
      attr(object_geno, "cleared") <- "model_ready_use"
      return(list(geno_pheno_match_data = object_geno,
                  test_data = test_set))
    } else {
      if (!is.null(train_set) && exists(ID_pheno)) {
        train_set <- check_test_set(test_set) ## For training set
        test_set <- setdiff(ID_pheno, train_set)
        test_set <- check_test_set(test_set)

        if (!is.null(test_set)) {
          if(nrow(object_geno)==ncol(object_geno)){
            object_geno <- object_geno[c(ID_pheno, test_set), c(ID_pheno, test_set)]
          } else {
            object_geno <- object_geno[c(ID_pheno, test_set), ]
          }
          attr(object_geno, "cleared") <- "model_ready_use"
          return(list(geno_pheno_match_data = object_geno,
                      test_data = test_set))
        }
      }
    }
  } else {
    if(nrow(object_geno)==ncol(object_geno)){
      object_geno <- object_geno[ID_pheno, ID_pheno]
    } else {
      object_geno <- object_geno[ID_pheno,  ]
    }
    attr(object_geno, "cleared") <- "model_ready_use"
    return(list(geno_pheno_match_data = object_geno
                ))
  }


}

