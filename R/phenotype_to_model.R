#' Title
#'
#' @param pheno_data
#' @param pheno_data_train
#' @param pheno_data_test
#' @param response
#' @param gen_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

phenotype_to_model <- function(
    pheno_data = NULL,
    pheno_data_train = NULL,
    pheno_data_test = NULL,
    #train_set = NULL,
    #test_set = NULL,
    response=NULL,
    gen_name=NULL,
    ...
) {

  msg <- sprintf("==================================================\n")


  ## Check availability of pheno_datatypic data (training and testing set). This
  ## accommodate missing value with the assumption that testing will have NA

  if (!is.null(pheno_data)){

    pheno_data <- phenotype_precheck(pheno_data= pheno_data,
                                 gen_name = gen_name,
                                 response = response)

    ### The result object has to pass the test attribute before it can be stored
    if(attr(pheno_data, "cleared")=="pass" && all(class(pheno_data)==c("data.frame", "phenotype"))) {

    # Assign appropriate class.
    class(pheno_data) <- c("data.frame", "phenotype")

    attr(pheno_data, "cleared") <- "for_model_fit"

    output =  list(pheno_data)

    names(output) <- c("pheno_data")

    #rm(pheno_data)

    } else {

      stop(print(paste(msg, "Data is not class phenotype.")), call. = FALSE)
    }

  } else {

    ### If pheno_data is missing, pheno_data traning should be available for traning and it is expected
    ## pheno_data testing is also available
    ## However, pheno_data testing might be missing if the user objective for cross-validation
    ## In general, missing value is not expected in the pheno_data training.
    if (!is.null(pheno_data_train)){

      pheno_data_train <- phenotype_precheck(pheno_data= pheno_data_train,
                                   gen_name = gen_name,
                                   response = response)

      if(attr(pheno_data_train, "cleared")!="pass" && all(class(pheno_data_train)!=c("data.frame", "phenotype"))) {

        stop(paste(msg, 'pheno_data_train is not object phenotype'))


      }

      if(anyNA(pheno_data_train)){ stop(paste(msg,"Missing value in not accepted in training set"))}

    }

    ### pheno_data test can be present or absent. It is expected this will contain
    ## missing value for the response variable

    if (!is.null(pheno_data_test)){

      pheno_data_test <- phenotype_precheck(pheno_data= pheno_data_test,
                                   gen_name = gen_name,
                                   response = response)

      if(attr(pheno_data_test, "cleared")!="pass" && all(class(pheno_data_test)!=c("data.frame", "phenotype"))) {

        stop(paste(msg, 'pheno_data_train is not object phenotype'))
      }
    }

    #### if both pheno_data_train and pheno_data_test are provided

    if ((!is.null(pheno_data_train) & !is.null(pheno_data_test)) & is.null(pheno_data)){

      if (!identical(colnames(pheno_data_train), colnames(pheno_data_test))){
        stop(print(paste(msg,'Columns name in the pheno_data_train not the same as pheno_data_test')), call. = FALSE)

      } else{

        pheno_data <- rbind(pheno_data_train, pheno_data_test)

        test_set <- data.frame(name = as.character(unique(pheno_data_test[, gen_name])))
        names(test_set) = gen_name

        ### The result object has to pass the test attribute before it can be stored

          # Assign appropriate class.
          class(pheno_data) <- c("data.frame", "phenotype")

          attr(pheno_data, "cleared") <- "for_model_fit"



      }

    }



  }



  #if(!is.null(pheno_data) & (!is.null(test_set) && !is.null(train_set))){

  if(!is.null(pheno_data) & exists('test_set')){

    if(!is.null(test_set)){

    output =  list(pheno_data, test_set)


    names(output) <- c("pheno_data", "test_set")

    }

    #rm(pheno_data, test_set)

    }else if (!is.null(pheno_data) & (is.null(pheno_data_train) & is.null(pheno_data_test))){

      output =  list(pheno_data)

      names(output) <- c("pheno_data")

      #rm(pheno_data)

  } else {


    if (!is.null(pheno_data_train) & !is.null(pheno_data_test)){

    output =  list(pheno_data, test_set)

    names(output) <- c("pheno_data", "test_set")

    #rm(pheno_data, test_set)

    }

  }

  return(output)
}
