#' Title
#'
#' @param pheno
#' @param pheno_train
#' @param pheno_test
#' @param train_set
#' @param test_set
#' @param response
#' @param gen_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

phenotype_to_model <- function(
    pheno = NULL,
    pheno_train = NULL,
    pheno_test = NULL,
    train_set = NULL,
    test_set = NULL,
    response=NULL,
    gen_name=NULL,
    ...
) {

  msg <- sprintf("==================================================\n")


  ## Check availability of phenotypic data (training and testing set). This
  ## accommodate missing value with the assumption that testing will have NA

  if (!is.null(pheno)){

    pheno <- phenotype_precheck(pheno= pheno,
                                 gen_name = gen_name,
                                 response = response)

    ### The result object has to pass the test attribute before it can be stored
    if(attr(pheno, "cleared")=="pass" && all(class(pheno)==c("data.frame", "phenotype"))) {

    # Assign appropriate class.
    class(pheno) <- c("data.frame", "phenotype")

    attr(pheno, "cleared") <- "for_model_fit"

    output =  list(pheno)

    names(output) <- c("pheno_data")

    rm(pheno)

    } else {

      stop(print(paste(msg, "Data is not class phenotype.")), call. = FALSE)
    }

  } else {

    ### If pheno is missing, pheno traning should be available for traning and it is expected
    ## pheno testing is also available
    ## However, pheno testing might be missing if the user objective for cross-validation
    ## In general, missing value is not expected in the pheno training.
    if (!is.null(pheno_train)){

      pheno_train <- phenotype_precheck(pheno= pheno_train,
                                   gen_name = gen_name,
                                   response = response)

      if(attr(pheno_train, "cleared")=="pass" && all(class(pheno_train)==c("data.frame", "phenotype"))) {

        stop('pheno_train is not object phenotype')
      }

    }

    ### pheno test can be present or absent. It is expected this will contain
    ## missing value for the response variable

    if (!is.null(pheno_test)){

      pheno_test <- phenotype_precheck(pheno= pheno_test,
                                   gen_name = gen_name,
                                   response = response)

      if(attr(pheno_test, "cleared")=="pass" && all(class(pheno_test)==c("data.frame", "phenotype"))) {

        stop('pheno_train is not object phenotype')
      }
    }

    #### if both pheno_train and pheno_test are provided

    if ((!is.null(pheno_train) && !is.null(pheno_test)) & is.null(pheno)){

      if (!identical(colnames(pheno_train), colnames(pheno_test))){
        stop(print(paste(msg,'Columns name in the pheno_train not the same as pheno_test')), call. = FALSE)

      } else{

        pheno <- rbind(pheno_train, pheno_test)

        test_set <- as.character(unique(pheno_test[, gen_name]))

        ### The result object has to pass the test attribute before it can be stored

          # Assign appropriate class.
          class(pheno) <- c("data.frame", "phenotype")

          attr(pheno, "cleared") <- "for_model_fit"



      }

    }



  }



  if(!is.null(pheno) & (!is.null(test_set) && !is.null(train_set))){

    output =  list(pheno, test_set, train_set)

    names(output) <- c("pheno_data", "test_set", "train_set")



  } else {


    if (!is.null(pheno_train) & !is.null(pheno_test)){

    output =  list(pheno, test_set)

    names(output) <- c("pheno_data", "test_set")



    }

  }

  return(output)
}
