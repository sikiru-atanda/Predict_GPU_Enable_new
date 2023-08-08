
#' Title
#'
#' @param gen_name
#' @param core
#' @param message
#' @param center
#' @param para_tunning
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param knn_paras_tunning
#' @param k
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_knn <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   knn_paras_tunning= c(k_tune=NULL),
                   k = 5,
                   ...

){

  #msg <- sprintf("==================================================\n")
  ##############################################################
  ###################################################################
  ###  Start ML Analysis
  ###
 ######################################################################
  #####################################################################


# #### Initializing parallel
if(length(response)>1){
  if (is.null(core)){
    cl = parallel::detectCores()

    if (cl> 4){
      # Try in parallel
      cl <- parallel::makeCluster(4)
    } else{
      cl <- parallel::makeCluster(2)
    }

  } else {
    if(!is.null(core)){

      cl <- parallel::makeCluster(core)
    }
  }

  doParallel::registerDoParallel(cl)



  Univariate <- foreach::foreach(trait = 1:length(response),
                                 .errorhandling='pass') %dopar% {

  if(isTRUE(para_tunning)){
    # create hyperparameter grid

    AI_grid <- expand.grid(
      #k = data.frame(k = seq(11,85,by = 2)))

      k = data.frame(k = k_tune)

    )

    #here we do one better then a validation set, we use cross validation to
    #expand the amount of info we have!

    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    AI_trcontrol = caret::trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)




    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


      AI_fit = caret::train(x = geno_omic_object,
                             y = pheno_object[, response[trait]],
                             trControl = AI_trcontrol,
                             tuneGrid = AI_grid,
                             method = "knn")



      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      bestTune <- c(AI_fit$bestTune, AI_fit$method)

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)


    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {


          AI_fit = caret::knnreg(x = geno_omic_object,
                                 y = pheno_object[, response[trait]],
                                 k = k
                                    )


        }


        AI_preds <- stats::predict(AI_fit,
                                   geno_omic_object,
                                   reshape = TRUE)


        AI_preds <- as.data.frame(AI_preds)

        names(AI_preds) = response[trait]


        # res_feature = caret::varImp(AI_fit)[[1]]
        # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
        # res_feature$feature = rownames(res_feature)


       }


    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
  } else {


    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {

      AI_fit = caret::knnreg(x = geno_omic_object,
                             y = pheno_object[, response[trait]],
                             k = k
      )

      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

        AI_fit = caret::knnreg(x = geno_omic_object,
                               y = pheno_object[, response[trait]],
                               k = k
        )


      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]



      }

    }


  } ## End of when no need for tunning.


     model_para <- c(k = AI_fit$k
                    )

output = list(model_para,
              AI_preds,
              AI_fit)


names(output) = c("model_parameters",
                  "predicted_values",
                  "trained_model")


Univariate = output


    }

  names(Univariate) <- response

  output <-  Univariate

  rm(Univariate)

  ### End when length of response variable is more than 1

   } else {

## when length of response variable is 1
#########################

     if(isTRUE(para_tunning)){
       # create hyperparameter grid

       AI_grid <- expand.grid(
         #k = data.frame(k = seq(11,85,by = 2)))

         k = data.frame(k = k_tune)

       )

       #here we do one better then a validation set, we use cross validation to
       #expand the amount of info we have!

       # if(core){
       # cl <- parallel::makeCluster(core)
       # doParallel::registerDoParallel(cl)
       # }

       AI_trcontrol = caret::trainControl(method = "cv",
                                          number = 5,
                                          verboseIter = TRUE,
                                          returnData = FALSE,
                                          returnResamp = "all",
                                          allowParallel = TRUE)




       if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


         AI_fit = caret::train(x = geno_omic_object,
                               y = pheno_object[, response],
                               trControl = AI_trcontrol,
                               tuneGrid = AI_grid,
                               method = "knn")



         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         bestTune <- c(AI_fit$bestTune, AI_fit$method)

         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)


       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {


           AI_fit = caret::knnreg(x = geno_omic_object,
                                  y = pheno_object[, response],
                                  k = k
           )


         }


         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_object,
                                    reshape = TRUE)


         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response


         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)


       }


       #######################################################
       ##########################################################
       ####### Start when their is no need for tunning
       ########################################################
       ########################################################
     } else {


       if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {

         AI_fit = caret::knnreg(x = geno_omic_object,
                                y = pheno_object[, response],
                                k = k
         )

         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)



       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

           AI_fit = caret::knnreg(x = geno_omic_object,
                                  y = pheno_object[, response],
                                  k = k
           )


           AI_preds <- stats::predict(AI_fit,
                                      geno_omic_object,
                                      reshape = TRUE)

           AI_preds <- as.data.frame(AI_preds)

           names(AI_preds) = response



         }

       }


     } ## End of when no need for tunning.





     model_para <- c(k = AI_fit$k)


     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")


}


return(output)



 }
