
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
#' @param ranForest_paras_tunning
#' @param ntree
#' @param mtry
#' @param maxnodes
#' @param importance
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_randomForest <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   ranForest_paras_tunning= c(mtry_tune=NULL),
                   ntree=500,
                   mtry = NULL,
                   maxnodes = NULL,
                   importance=TRUE,
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
    n_features = ncol(geno_omic_object)

    AI_grid <- expand.grid(
      mtry = floor(n_features * ranForest_paras_tunning$mtry_tune)

    )


    #c(.05, .15)

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
                             method = "rf")

      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      bestTune <- c(AI_fit$bestTune, AI_fit$method)

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)

      res_feature = caret::varImp(AI_fit)[[1]]
      res_feature = data.frame(feature = colnames(geno_omic_object),
                               weight = res_feature)
      res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

         if(is.null(mtry) & is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              importance = importance)

        } else if (!is.null(mtry) & is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              mtry = mtry,
                                              importance = importance)


        } else if (!is.null(mtry) & !is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              mtry = mtry,
                                              maxnodes = maxnodes,
                                              importance = importance)


        } else {

          if (is.null(mtry) & !is.null(maxnodes)){
            AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                y = pheno_object[, response[trait]],
                                                ntree = ntree,
                                                maxnodes = maxnodes,
                                                importance = importance)

          }

        }


        AI_preds <- stats::predict(AI_fit,
                                   geno_omic_object,
                                   reshape = TRUE)


        AI_preds <- as.data.frame(AI_preds)

        names(AI_preds) = response[trait]


        # res_feature = caret::varImp(AI_fit)[[1]]
        # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
        # res_feature$feature = rownames(res_feature)

        res_feature = caret::varImp(AI_fit)[[1]]
        res_feature = data.frame(feature = colnames(geno_omic_object),
                                 weight = res_feature)
        res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]



      }

    }


    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
  } else {


    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


      if(is.null(mtry) & is.null(maxnodes)){
        AI_fit = randomForest::randomForest(x = geno_omic_object,
                                            y = pheno_object[, response[trait]],
                                            ntree = ntree,
                                            importance = importance)

      } else if (!is.null(mtry) & is.null(maxnodes)){
        AI_fit = randomForest::randomForest(x = geno_omic_object,
                                            y = pheno_object[, response[trait]],
                                            ntree = ntree,
                                            mtry = mtry,
                                            importance = importance)


      } else if (!is.null(mtry) & !is.null(maxnodes)){
        AI_fit = randomForest::randomForest(x = geno_omic_object,
                                            y = pheno_object[, response[trait]],
                                            ntree = ntree,
                                            mtry = mtry,
                                            maxnodes = maxnodes,
                                            importance = importance)


      } else {

        if (is.null(mtry) & !is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              maxnodes = maxnodes,
                                              importance = importance)

        }

      }



      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)
      res_feature = caret::varImp(AI_fit)[[1]]
      res_feature = data.frame(feature = colnames(geno_omic_object),
                               weight = res_feature)
      res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]




    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        if(is.null(mtry) & is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              importance = importance)

        } else if (!is.null(mtry) & is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              mtry = mtry,
                                              importance = importance)


        } else if (!is.null(mtry) & !is.null(maxnodes)){
          AI_fit = randomForest::randomForest(x = geno_omic_object,
                                              y = pheno_object[, response[trait]],
                                              ntree = ntree,
                                              mtry = mtry,
                                              maxnodes = maxnodes,
                                              importance = importance)


        } else {

          if (is.null(mtry) & !is.null(maxnodes)){
            AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                y = pheno_object[, response[trait]],
                                                ntree = ntree,
                                                maxnodes = maxnodes,
                                                importance = importance)

          }

        }


      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      res_feature = caret::varImp(AI_fit)[[1]]
      res_feature = data.frame(feature = colnames(geno_omic_object),
                               weight = res_feature)
      res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]




      }

    }


  } ## End of when no need for tunning.


     model_para <- c(ntree = AI_fit$ntree,
                     mtry = AI_fit$mtry
                    )

output = list(model_para,
              AI_preds,
              res_feature,
              AI_fit)


names(output) = c("model_parameters", "predicted_values",
                  "feature_weight", "trained_model")


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
       n_features = ncol(geno_omic_object)

       AI_grid <- expand.grid(
         mtry = floor(n_features * ranForest_paras_tunning$mtry_tune)

       )


       #c(.05, .15)

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
                               method = "rf")

         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         bestTune <- c(AI_fit$bestTune, AI_fit$method)

         res_feature = caret::varImp(AI_fit)[[1]]
         res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         res_feature$feature = rownames(res_feature)



       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

           if(is.null(mtry) & is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 importance = importance)

           } else if (!is.null(mtry) & is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 mtry = mtry,
                                                 importance = importance)


           } else if (!is.null(mtry) & !is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 mtry = mtry,
                                                 maxnodes = maxnodes,
                                                 importance = importance)


           } else {

             if (is.null(mtry) & !is.null(maxnodes)){
               AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                   y = pheno_object[, response],
                                                   ntree = ntree,
                                                   maxnodes = maxnodes,
                                                   importance = importance)

             }

           }


           AI_preds <- stats::predict(AI_fit,
                                      geno_omic_object,
                                      reshape = TRUE)


           AI_preds <- as.data.frame(AI_preds)

           names(AI_preds) = response


           res_feature = caret::varImp(AI_fit)[[1]]
           res_feature = data.frame(feature = colnames(geno_omic_object),
                                    weight = res_feature)
           res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]




         }

       }


       #######################################################
       ##########################################################
       ####### Start when their is no need for tunning
       ########################################################
       ########################################################
     } else {


       if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


         if(is.null(mtry) & is.null(maxnodes)){
           AI_fit = randomForest::randomForest(x = geno_omic_object,
                                               y = pheno_object[, response],
                                               ntree = ntree,
                                               importance = importance)

         } else if (!is.null(mtry) & is.null(maxnodes)){
           AI_fit = randomForest::randomForest(x = geno_omic_object,
                                               y = pheno_object[, response],
                                               ntree = ntree,
                                               mtry = mtry,
                                               importance = importance)


         } else if (!is.null(mtry) & !is.null(maxnodes)){
           AI_fit = randomForest::randomForest(x = geno_omic_object,
                                               y = pheno_object[, response],
                                               ntree = ntree,
                                               mtry = mtry,
                                               maxnodes = maxnodes,
                                               importance = importance)


         } else {

           if (is.null(mtry) & !is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 maxnodes = maxnodes,
                                                 importance = importance)

           }

         }



         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         res_feature = caret::varImp(AI_fit)[[1]]
         res_feature = data.frame(feature = colnames(geno_omic_object),
                                   weight = res_feature)
         res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]




       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
           if(is.null(mtry) & is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 importance = importance)

           } else if (!is.null(mtry) & is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 mtry = mtry,
                                                 importance = importance)


           } else if (!is.null(mtry) & !is.null(maxnodes)){
             AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                 y = pheno_object[, response],
                                                 ntree = ntree,
                                                 mtry = mtry,
                                                 maxnodes = maxnodes,
                                                 importance = importance)


           } else {

             if (is.null(mtry) & !is.null(maxnodes)){
               AI_fit = randomForest::randomForest(x = geno_omic_object,
                                                   y = pheno_object[, response],
                                                   ntree = ntree,
                                                   maxnodes = maxnodes,
                                                   importance = importance)

             }

           }


           AI_preds <- stats::predict(AI_fit,
                                      geno_omic_object,
                                      reshape = TRUE)

           AI_preds <- as.data.frame(AI_preds)

           names(AI_preds) = response

           res_feature = caret::varImp(AI_fit)[[1]]
           res_feature = data.frame(feature = colnames(geno_omic_object),
                                    weight = res_feature)
           res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]




         }

       }


     } ## End of when no need for tunning.


     model_para <- c(ntree = AI_fit$ntree,
                     mtry = AI_fit$mtry
     )

     output = list(model_para,
                   AI_preds,
                   res_feature,
                   AI_fit)


     names(output) = c("model_parameters", "predicted_values",
                       "feature_weight", "trained_model")


}


return(output)



 }
