
#' Title
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
#' @param geno_data
#' @param train_geno_data
#' @param test_geno_data
#' @param message
#' @param ...
#' @param map_data
#'
#' @return
#' @export
#'
#' @examples
geno_to_modelOLD <- function(geno_data = NULL,
                            train_geno_data = NULL,
                            test_geno_data = NULL,
                            map_data = NULL,
                            message= TRUE,
                            ...){

  msg <- sprintf("==================================================\n")

  if(((!is.null(geno_data) & is.null(map_data)) & (is.null(test_geno_data) & is.null(train_geno_data)))){
    geno_object = geno_precheck(object_geno = geno_data,
                                message = message)

    geno_object = geno_object[[1]]

    ### It has to pass test before it can be declared geno_object
    if(attr(geno_object, "cleared")=="pass" & all(class(geno_object)==c("matrix", "array", "geno_data"))){

      ### This step might look redundant, but is not
      #### Declare it also as an object for final usage
      class(geno_object) <-c("matrix", "array", "geno_data")

      attr(geno_object, "cleared") <- "for_model_fit"
      ##############################################

    } else {

      geno_object = NULL
    }

  } else if(((!is.null(geno_data) & !is.null(map_data)) & (is.null(test_geno_data) & is.null(train_geno_data)))){

    geno_object = geno_precheck(object_geno = geno_data,
                                message = message)

    geno_object = geno_object[[1]]


    geno_QC_summary_stat = geno_object[[2]]

    ### It has to pass test before it can be declared geno_object
    if(attr(geno_object, "cleared")=="pass" & all(class(geno_object)==c("matrix", "array", "geno_data"))){

      ### This step might look redundant, but is not
      #### Declare it also as an object for final usage
      class(geno_object) <-c("matrix", "array", "geno_data")

      attr(geno_object, "cleared") <- "for_model_fit"
      ##############################################

    } else {

      stop(print(paste(msg,'Geno or the omic data did not pass the required test. Check the data.')), call. = FALSE)
      #geno_object = NULL
    }


  }else {

    if(!is.null(train_geno_data)){
      train_geno_data = geno_precheck(object_geno = train_geno_data,
                                      message = message)

      train_geno_data =  train_geno_data[[1]]

      geno_QC_summary_stat = train_geno_data[[2]]
    }

    if(!is.null(test_geno_data)){
      test_geno_data = geno_precheck(object_geno = test_geno_data,
                                     message = message)

      test_geno_data =  test_geno_data[[1]]

      geno_QC_summary_stat = test_geno_data[[2]]

    }


    ### Both test_geno_data and train_geno_data has to pass through the filter
    ## test before they will be merge
    # if(attr(test_geno_data, "cleared")!="pass" & class(test_geno_data)!=c("matrix", "array", "geno_data")) {
    #
    #   msg <- sprintf("==================================================\n")
    #   stop(print(paste(msg,'Test_geno_data is not class formula.')), call. = FALSE)
    # }

#### Begin process to merge train_geno_data and test_geno_data.
    if(!is.null(train_geno_data) & !is.null(test_geno_data)){

      ## Both test_geno_data and train_geno_data has to pass through the filter
      ## test before they will be merge
    if(((attr(test_geno_data, "cleared")=="pass" && all(class(test_geno_data)==c("matrix", "array", "geno_data"))) & (attr(train_geno_data, "cleared")=="pass" && all(class(train_geno_data)==c("matrix", "array", "geno_data"))))){

      ## In some scenerio after QC the number of SNP in the test_geno_data and
      ## train_geno_data might differ
      if (!identical(colnames(train_geno_data), colnames(test_geno_data))) {

        ### Check if the SNP /columns is not zero
        if(dim(train_geno_data)[2]!=0 & dim(test_geno_data)[2]!=0){

          ## Identify SNP common to both train_geno_data and test_geno_data after QC
          snp_names <- intersect(colnames(train_geno_data),
                                 colnames(test_geno_data))

          train_geno_data <- train_geno_data[, colnames(train_geno_data)%in%snp_names]
          test_geno_data <-  test_geno_data[, colnames(test_geno_data)%in%snp_names]

          geno_object <- rbind(train_geno_data, test_geno_data)

          ### This step might look redundant, but is not
          #### Declare it also as an object for final usage
          class(geno_object) <-c("matrix", "array", "geno_data")

          attr(geno_object, "cleared") <- "for_model_fit"

          #test_set_ <- data.frame(name = rownames(test_geno_data), stringsAsFactors = FALSE)


        } else {

          stop(print(paste(msg,'SNP/markers did not match in training and testing set data.')), call. = FALSE)

        }



      } else {

        geno_object <- rbind(train_geno_data, test_geno_data)

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"
        ##############################################

      }

    } else{

      stop(print(paste(msg,'SNP/marker data did not pass the required test. Check the data.')), call. = FALSE)
      #geno_object <- NULL

    }

      ### If only the train_geno_data is provided by the user
    } else if(!is.null(train_geno_data) & is.null(test_geno_data)){

      if(isTRUE(message)){
        message(insight::print_color(paste(msg,paste("Only train_geno_data is provided.")), "blue"))

      }
      ### It has to pass test before it can be declared geno_object
      if(attr(train_geno_data, "cleared")=="pass" & all(class(train_geno_data)==c("matrix", "array", "geno_data"))){

        geno_object <- train_geno_data

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"

        rm(train_geno_data)
        ##############################################

      }


    } else {

      ## if only the test_geno_data is provided by the test_geno_data
      if(is.null(train_geno_data) & !is.null(test_geno_data)){

        message(paste(insight::print_color("WARNINGS\n", "blue"),
                         insight::print_color(paste(msg,paste("Only the test_geno_data is provided.\n\t Check if this is correct.")), "blue")))

        #message(paste(msg, 'Only test_geno_data is provided'))


        ### It has to pass test before it can be declared geno_object
        if(attr(test_geno_data, "cleared")=="pass" && all(class(test_geno_data)!=c("matrix", "array", "geno_data"))){

        geno_object <- test_geno_data

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"

        rm(test_geno_data)
        ##############################################
        }

      }

  }

# NOTE: We can include metadata is not yet included as part of the output but was generated
    ## Here is an instance where user provide snp/marker data directly
  }

return(geno_object)


}

