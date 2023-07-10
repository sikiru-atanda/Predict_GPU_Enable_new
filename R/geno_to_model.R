
#' Title
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
geno_to_model <- function(geno_data = NULL,
                          train_geno_data = NULL,
                          test_geno_data = NULL,
                          map_data = NULL,
                          message= TRUE,
                          ...){

  msg <- sprintf("==================================================\n")

  if(!is.null(geno_data) & is.null(map_data)){
    geno_object = geno_precheck(object = geno_data,
                                message = message)

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

  } else if(!is.null(geno_data) & is.null(map_data)){




  }else {

    if(!is.null(train_geno_data)){
      train_geno_data = geno_precheck(object = train_geno_data,
                                      message = message)
    }

    if(!is.null(test_geno_data)){
      test_geno_data = geno_precheck(object = test_geno_data,
                                     message = message)

    }


    ### Both test_geno_data and train_geno_data has to pass through the filter
    ## test before they will be merge
    # if(attr(test_geno_data, "cleared")!="pass" & class(test_geno_data)!=c("matrix", "array", "geno_data")) {
    #
    #   msg <- sprintf("==================================================\n")
    #   stop(print(paste(msg,'Test_geno_data is not class formula.')), call. = FALSE)
    # }


    if(!is.null(train_geno_data) & !is.null(test_geno_data)){

      ## Both test_geno_data and train_geno_data has to pass through the filter
      ## test before they will be merge
    if((attr(test_geno_data, "cleared")=="pass" && all(class(test_geno_data)==c("matrix", "array", "geno_data"))) & (attr(train_geno_data, "cleared")=="pass" && all(class(train_geno_data)==c("matrix", "array", "geno_data")))){

      if (!identical(colnames(train_geno_data), colnames(test_geno_data))) {

        stop(print(paste(msg,'SNP/markers did not match in training and testing set data')), call. = FALSE)

      } else {

        geno_object <- rbind(train_geno_data, test_geno_data)

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"
        ##############################################

      }

    } else{
      geno_object <- NULL

    }

    } else if(!is.null(train_geno_data) & is.null(test_geno_data)){

      message (paste(msg, 'Only train_geno_data is provided'))
      ### It has to pass test before it can be declared geno_object
      if(attr(train_geno_data, "cleared")=="pass" & all(class(train_geno_data)==c("matrix", "array", "geno_data"))){

        geno_object <- train_geno_data

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"

        rm(train_geno_data)
        ##############################################

      } else {

        geno_object = NULL
      }

    } else {

      if(is.null(train_geno_data) & !is.null(test_geno_data)){

        message(paste(msg, 'Only test_geno_data is provided'))

        ### It has to pass test before it can be declared geno_object
        if(attr(test_geno_data, "cleared")=="pass" && all(class(test_geno_data)!=c("matrix", "array", "geno_data"))){

        geno_object <- test_geno_data

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        class(geno_object) <-c("matrix", "array", "geno_data")

        attr(geno_object, "cleared") <- "for_model_fit"

        rm(test_geno_data)
        ##############################################
        } else {

          geno_object = NULL
        }

      }

  }


  }

return(geno_object)


}

