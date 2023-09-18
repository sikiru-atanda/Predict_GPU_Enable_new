
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
#' @param svm_paras_tunning
#' @param c
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_svm <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   svm_paras_tunning= c(c_tune=NULL),
                   c = 1,
                   ...

){

  #msg <- sprintf("==================================================\n")
  ##############################################################
  ###################################################################
  ###  Start ML Analysis
  ###
 ######################################################################
  #####################################################################
# http://dataworldblog.blogspot.com/2017/08/support-vector-machines-svm-in-r.html

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

      C = data.frame(C = c_tune)

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
                             scaled  = FALSE,
                             method = "svmLinear")



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


          AI_fit = kernlab::ksvm(x = geno_omic_object,
                                 y = pheno_object[, response[trait]],
                                 C = c,
                                 scaled = FALSE
                                    )


        }


        AI_preds <- kernlab::predict(AI_fit,
                                   geno_omic_object)


        AI_preds <- as.data.frame(AI_preds)

        names(AI_preds) = response[trait]


        # res_feature = caret::varImp(AI_fit)[[1]]
        # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
        # res_feature$feature = rownames(res_feature)


       }

    model_para <- bestTune
    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
  } else {


    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {

      AI_fit = kernlab::ksvm(x = geno_omic_object,
                             y = pheno_object[, response[trait]],
                             scaled  = FALSE,
                             C = c
      )

      AI_preds <- kernlab::predict(AI_fit,
                                  geno_omic_test_object)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

        AI_fit = kernlab::ksvm(x = geno_omic_object,
                               y = pheno_object[, response[trait]],
                               scaled  = FALSE,
                               C = c
        )


      AI_preds <- kernlab::predict(AI_fit,
                                  geno_omic_object)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]



      }

    }

    model_para <- c(epsilon = AI_fit@param$epsilon,
                    C = AI_fit@param$C
    )

  } ## End of when no need for tunning.




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

         C = data.frame(C = c_tune)

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
                               scaled  = FALSE,
                               method = "svmLinear")



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


           AI_fit = kernlab::ksvm(x = geno_omic_object,
                                  y = pheno_object[, response],
                                  scaled  = FALSE,
                                  C = c
           )


         }


         AI_preds <- kernlab::predict(AI_fit,
                                    geno_omic_object)


         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response


         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)


       }


       model_para <- bestTune


       #######################################################
       ##########################################################
       ####### Start when their is no need for tunning
       ########################################################
       ########################################################
     } else {


       if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {

         AI_fit = kernlab::ksvm(x = geno_omic_object,
                                y = pheno_object[, response],
                                scaled  = FALSE,
                                C = c
         )

         AI_preds <- kernlab::predict(AI_fit,
                                    geno_omic_test_object)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)



       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {

           AI_fit = kernlab::ksvm(x = geno_omic_object,
                                  y = pheno_object[, response],
                                  scaled  = FALSE,
                                  C = c
           )


           AI_preds <-  kernlab::predict(AI_fit,
                                      geno_omic_object
                                      )

           AI_preds <- as.data.frame(AI_preds)

           names(AI_preds) = response



         }

       }

       model_para <- c(epsilon = AI_fit@param$epsilon,
                       C = AI_fit@param$C)

     } ## End of when no need for tunning.








     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")


}


return(output)



 }



