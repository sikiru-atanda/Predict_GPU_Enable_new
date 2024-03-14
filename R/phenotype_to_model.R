#' Title
#' Overall, this serve as gateway between phenotype-precheck function and readiness of
#' the phenotypic data for model fitting
#' The objective of this function is to do the following:
#' 1. Check the output from phenotype-precheck function for pheno_data before declaring it for model fit.
#'    NA is allowed for response variable with the assumption that those individuals are the testing set.
#' 2. Check output from phenotype-precheck function for pheno_data_train and
#'    check for NA. If NA present it will not be declare for model fit because
#'    NA is not expected in training set. If all good it will be declared for model fit
#' 3. Check the output from phenotype-precheck function for pheno_data_test before declaring it for model fit
#'    NA is allowed for response variable here.
#' 4. If pheno_data_train and pheno_data_test were present and pass through the pre-check process,
#'    it will be processed that is.
#'    1) it check that column name for both data match/the same
#'    2) Combined the dataset for model fit and prediction.
#'    It is assumed here that pheno_data is missing/not provided by the user.
#' 5. The output will be pheno_data declared for model fit and test_set if pheno_data_test was provided.
#'
#' @param pheno_data phenotypic object NA is allowed
#' @param pheno_data_train phenotypic object for the training set. NA is allowed.
#' This is expected if pheno_data was not provided
#' @param pheno_data_test phenotypic object for the testing set. NA is allowed
#' @param response y variables/lables
#' @param gen_name column name containing individuals/genotypes
#' @param ...
#' @param train_set
#' @param test_set
#'
#' @return
#' @export
#'
#' @examples

phenotype_to_model <- function(
                              pheno_data = NULL,
                              pheno_data_train = NULL,
                              pheno_data_test = NULL,
                              train_set = NULL,
                              test_set = NULL,
                              response=NULL,
                              gen_name=NULL,
                              heter_groups = NULL,
                              random = NULL,
                              fixed = NULL,
                              ...
                          ) {

msg <- sprintf("==================================================\n")


  ## Check availability of pheno_datatypic data (training and testing set). This
  ## accommodate missing value with the assumption that testing will have NA
  ## Check for NA is not test here for response variables because testing set might have NA
  # thus, NA is expected in pheno_data.

  if ((!is.null(pheno_data) & (is.null(pheno_data_train) & is.null(pheno_data_test)))){

    pheno_data <- phenotype_precheck(pheno_data= pheno_data,
                                     gen_name = gen_name,
                                     response = response,
                                     heter_groups = heter_groups,
                                     random = random,
                                     fixed, fixed)

    if(!is.null(test_set)){

      if(is.data.frame(test_set) | is.matrix(test_set)){

       test_set = test_set[, 1]


      } else if (!is.list(test_set)){

        test_set = test_set

      } else {
        if(is.list(test_set)){
        stop(message(paste(msg,'The testing set cannot be a list. Should be either dataframe, matrix or a vector.')), call. = FALSE)

        }

      }

      if(length(test_set)>length(unique(as.character(pheno_data[, gen_name])))){
        stop(message(paste(msg, "The testing set size should be less than the unique genotypes in the pheno_data.")), call. = FALSE)
      }

      pheno_data[, response] <- ifelse(pheno_data[, gen_name]%in%test_set, NA,
                                       pheno_data[, response])

      test_set_ = test_set

      rm(test_set)

    } else{

      if(!is.null(train_set)){
        if(is.data.frame(train_set) | is.matrix(train_set)){
          train_set = train_set[, 1]

        } else if (!is.list(train_set)){

          train_set = train_set
        } else {
          if(is.list(train_set)){

            stop(message(paste(msg,'The training set cannot be a list. Should be either dataframe, matrix or a vector.')), call. = FALSE)
          }

          }

        pheno_data[, response] <- ifelse(!pheno_data[, gen_name]%in%train_set, NA,
                                         pheno_data[, response])

        test_set_ <- data.frame(name = as.character(unique(pheno[!pheno_data[, gen_name]%in%train_set, gen_name])), stringsAsFactors = FALSE)
        names(test_set) = gen_name

        if(nrow(test_set_)==0){

          rm(test_set_)

          message(paste( insight::print_color("WARNINGS\n", "blue"),
                         insight::print_color(paste(msg,paste("The training set size is the same size as the unique genotypes in the pheno_data.")), "blue")))

        }

      }

      }

    ### The result object has to pass the test attribute before it can be stored/
    ## pass through for the next step of check and declared good for model fit
    if(attr(pheno_data, "cleared")=="pass") {

    # Assign appropriate class.
    #class(pheno_data) <- c("data.frame", "phenotype")

    attr(pheno_data, "cleared") <- "for_model_fit"

    output <- list(pheno_clean_data = pheno_data)

    #names(output) <- c("pheno_data")

    #rm(pheno_data)

    } else {

      stop(print(paste(msg, "Data is not class phenotype.")), call. = FALSE)
    }

  } else {

    ### If pheno_data is missing, pheno_data traning should be available for traning and it is expected
    ## pheno_data testing is also available
    ## However, pheno_data testing might be missing if the user objective is for cross-validation
    ## In general, missing value is not expected in the pheno_data training.
    if (!is.null(pheno_data_train)){

      pheno_data_train_ <- phenotype_precheck(pheno_data= pheno_data_train,
                                              gen_name = gen_name,
                                              response = response,
                                              heter_groups = heter_groups,
                                              random = random,
                                              fixed, fixed)

      if(attr(pheno_data_train_, "cleared")!="pass") {

        stop(message(paste(msg, 'pheno_data_train is not object phenotype.')), call. = FALSE)



      }

      ## Though pass the pre-check test but NA is not expected in the pheno_data training.
      if(anyNA(pheno_data_train_)){

        stop(message(paste(msg, "Missing value in not accepted in training set.")), call. = FALSE)

        }

      rm(pheno_data_train)
    }

    ### pheno_data test can be present or absent. It is expected this will contain
    ## missing value for the response variable

    if (!is.null(pheno_data_test)){

      pheno_data_test_ <- phenotype_precheck(pheno_data= pheno_data_test,
                                             gen_name = gen_name,
                                             response = response)

      #if(attr(pheno_data_test_, "cleared")!="pass" && all(class(pheno_data_test_)!=c("data.frame", "phenotype"))) {
      if(attr(pheno_data_test_, "cleared")!="pass") {
        stop(print(paste(msg, 'pheno_data_train is not object phenotype.')), call. = FALSE)

      }

      rm(pheno_data_test)
    }

    #### if both pheno_data_train and pheno_data_test are provided

    if ((exists("pheno_data_train_") & exists("pheno_data_test_")) & is.null(pheno_data)){

      if (!identical(colnames(pheno_data_train_), colnames(pheno_data_test_))){
        stop(message(paste(msg,'Columns name in the pheno_data_train not the same as pheno_data_test.')), call. = FALSE)

      } else{

        pheno_data <- rbind(pheno_data_train, pheno_data_test)

        test_set_ <- data.frame(name = as.character(unique(pheno_data_test[, gen_name])), stringsAsFactors = FALSE)
        names(test_set) = gen_name

        ### The result object has to pass the test attribute before it can be stored

          # Assign appropriate class.
          #class(pheno_data) <- c("data.frame", "phenotype")

          attr(pheno_data, "cleared") <- "for_model_fit"



      }

    }



  }



  #if(!is.null(pheno_data) & (!is.null(test_set) && !is.null(train_set))){

  if(!is.null(pheno_data) & exists('test_set_')){

    if(!is.null(test_set_)){

    output <- list(pheno_clean_data = pheno_data,
                   test_set = test_set_)


    #names(output) <- c("pheno_data", "test_set")

    }

    #rm(pheno_data, test_set)

    }else if (!is.null(pheno_data) & (is.null(pheno_data_train) & is.null(pheno_data_test))){

      output <- list(pheno_clean_data = pheno_data)

      #names(output) <- c("pheno_data")

      #rm(pheno_data)

  } else {


    if (exists("pheno_data_train_") & exists("pheno_data_test_")){

    output <-   list(pheno_clean_data = pheno_data,
                   test_set = test_set_)

    #names(output) <- c("pheno_data", "test_set")

    #rm(pheno_data, test_set)

    }

  }

  return(output)
}
