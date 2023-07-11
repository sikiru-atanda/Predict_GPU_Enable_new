
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
omic_to_model <- function(omic_data = NULL,
                          train_omic_data = NULL,
                          test_omic_data = NULL,
                          message= TRUE,
                          ...){

  msg <- sprintf("==================================================\n")

  if(!is.null(omic_data)){
    omic_object = omic_precheck(object = omic_data,
                                message = message)

    if(attr(omic_object, "cleared")=="pass" & all(class(omic_object)==c("matrix", "array", "omic_matrix"))){

      omic_object = omic_object
      ##############################################

    } else {

      stop(print(paste(msg,'Data is not class Omic_matrix1')), call. = FALSE)
      #omic_object = NULL
    }

    ### It has to pass test before it can be declared omic_object
    # if(attr(omic_object, "cleared")!="pass" & all(class(omic_object)!=c("matrix", "array", "omic_matrix"))){
    #
    #   omic_object <- omic_object
    #
    #   ### This step might look redundant, but is not
    #   #### Declare it also as an object for final usage
    #   class(omic_object) <-c("matrix", "array", "omic_matrix")
    #
    #   attr(omic_object, "cleared") <- "for_model_fit"
    #   ##############################################
    #
    # } else {
    #
    #   stop(print(paste(msg,'Data is not class Omic_matrix')), call. = FALSE)
    #   #omic_object = NULL
    # }

  } else {

    if(!is.null(train_omic_data)){
      train_omic_data = omic_precheck(object = train_omic_data,
                                      message = message)
    }

    if(!is.null(test_omic_data)){
      test_omic_data = omic_precheck(object = test_omic_data,
                                     message = message)

    }


    ### Both test_omic_data and train_omic_data has to pass through the filter
    ## test before they will be merge
    # if(attr(test_omic_data, "cleared")!="pass" & class(test_omic_data)!=c("matrix", "array", "omic_data")) {
    #
    #   msg <- sprintf("==================================================\n")
    #   stop(print(paste(msg,'Test_omic_data is not class formula.')), call. = FALSE)
    # }


    if(!is.null(train_omic_data) && !is.null(test_omic_data)){

      ## Both test_omic_data and train_omic_data has to pass through the filter
      ## test before they will be merge
      if((attr(test_omic_data, "cleared")=="pass" & all(class(test_omic_data)==c("matrix", "array", "omic_matrix"))) & (attr(train_omic_data, "cleared")!="pass" & all(class(train_omic_data)!=c("matrix", "array", "omic_matrix")))){

        if (!identical(colnames(train_omic_data), colnames(test_omic_data))) {

          stop(print(paste(msg,'Data is not class Omic_matrix2.')), call. = FALSE)

        } else {

          omic_object <- rbind(train_omic_data, test_omic_data)

          ### This step might look redundant, but is not
          #### Declare it also as an object for final usage
          # class(omic_object) <-c("matrix", "array", "omic_matrix")
          #
          # attr(omic_object, "cleared") <- "for_model_fit"
          ##############################################

        }

      } else{

        stop(print(paste(msg,'Data is not object Omic_matrix3')), call. = FALSE)
        #omic_object <- NULL

      }

    } else if(!is.null(train_omic_data) & is.null(test_omic_data)){

      message (paste(msg, 'Only train_omic_data is provided'))
      ### It has to pass test before it can be declared omic_object
      if(attr(train_omic_data, "cleared")=="pass" & all(class(train_omic_data)==c("matrix", "array", "omic_matrix"))){

        omic_object <- train_omic_data

        ### This step might look redundant, but is not
        #### Declare it also as an object for final usage
        # class(omic_object) <-c("matrix", "array", "omic_matrix")
        #
        # attr(omic_object, "cleared") <- "for_model_fit"
        #
        # rm(train_omic_data)
        ##############################################

      } else {

        stop(print(paste(msg,'Data is not class Omic_matrix4')), call. = FALSE)

        #omic_object = NULL
      }

    } else {

      if(is.null(train_omic_data) & !is.null(test_omic_data)){

        message(paste(msg, 'Only test_omic_data is provided'))

        ### It has to pass test before it can be declared omic_object
        if(attr(test_omic_data, "cleared")=="pass" & all(class(test_omic_data)==c("matrix", "array", "omic_matrix"))){

          omic_object <- test_omic_data

          # ### This step might look redundant, but is not
          # #### Declare it also as an object for final usage
          # class(omic_object) <-c("matrix", "array", "omic_matrix")
          #
          # attr(omic_object, "cleared") <- "for_model_fit"

          rm(test_omic_data)
          ##############################################
        } else {
          stop(print(paste(msg,'Data is not class Omic_matrix5')), call. = FALSE)
          #omic_object = NULL
        }

      }

    }


  }

    ### This step might look redundant, but is not
    #### Declare it also as an object for final usage
    class(omic_object) <-c("matrix", "array", "omic_matrix")

    attr(omic_object, "cleared") <- "for_model_fit"
    ##############################################


  return(omic_object)


}

