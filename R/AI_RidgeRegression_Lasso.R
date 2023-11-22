
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
#' @param ...
#' @param Lasso_paras_tunning
#' @param lambda
#' @param GS_model
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_RidgeRegression_Lasso <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   Lasso_paras_tunning= c(lambda_tune=NULL),
                   lambda = NULL,
                   GS_model = c("Lasso",
                              "Ridge_Regression"),
                   ...

){

  #msg <- sprintf("==================================================\n")
  ##############################################################
  ###################################################################
  ###  Start ML Analysis
  ###
 ######################################################################
  #####################################################################
## http://www.science.smith.edu/~jcrouser/SDS293/labs/lab10-r.html
  if (is.null(lambda)){

    lambda = seq(0.000001,0.9,length.out=100)^4
  }

  if(GS_model=="Ridge_Regression"){

    alpha = 0
  }

  if(GS_model=="Lasso"){

    alpha = 1
  }

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

      lambda = data.frame(alpha = alpha, lambda = lambda_tune)

    )

    #data.frame(alpha = 1, lambda = seq(0.01,3,.01))
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
                             GS_model = "glmnet")



      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]

      bestTune <- c(AI_fit$bestTune, AI_fit$GS_model)

      # res_feature = caret::varImp(AI_fit)[[1]]
      # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
      # res_feature$feature = rownames(res_feature)


    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {


        AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                     y=pheno_object[, response[trait]],
                                     nfolds = 5,
                                     alpha = alpha,
                                     standardize = FALSE,
                                     lambda = lambda)

        AI_fit=  glmnet::glmnet(x=geno_omic_object,
                                y=pheno_object[, response[trait]],
                                alpha = alpha,
                                standardize = FALSE,
                                lambda =AI_fit_CV$lambda.min)



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

      AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                   y=pheno_object[, response[trait]],
                                   nfolds = 5,
                                   alpha = alpha,
                                   standardize = FALSE,
                                   lambda = lambda)

      AI_fit=  glmnet::glmnet(x=geno_omic_object,
                              y=pheno_object[, response[trait]],
                              alpha = alpha,
                              standardize = FALSE,
                              lambda =AI_fit_CV$lambda.min)


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

        AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                     y=pheno_object[, response[trait]],
                                     nfolds = 5,
                                     alpha = alpha,
                                     standardize = FALSE,
                                     lambda = lambda)

        AI_fit=  glmnet::glmnet(x=geno_omic_object,
                                y=pheno_object[, response[trait]],
                                alpha = alpha,
                                standardize = FALSE,
                                lambda =AI_fit_CV$lambda.min)



      AI_preds <- stats::predict(AI_fit,
                                  geno_omic_object,
                                  reshape = TRUE)

      AI_preds <- as.data.frame(AI_preds)

      names(AI_preds) = response[trait]



      }

    }


  } ## End of when no need for tunning.


     model_para <- c(lambda = AI_fit$lambda,
                     alpha = alpha
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

         lambda = data.frame(alpha = alpha, lambda = lambda_tune)

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
                               GS_model = "glment")



         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         AI_preds <- as.data.frame(AI_preds)

         names(AI_preds) = response

         bestTune <- c(AI_fit$bestTune, AI_fit$GS_model)

         # res_feature = caret::varImp(AI_fit)[[1]]
         # res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         # res_feature$feature = rownames(res_feature)


       } else {

         if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {


           AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                        y=pheno_object[, response],
                                        nfolds = 5,
                                        alpha = alpha,
                                        standardize = FALSE,
                                        lambda = lambda)

           AI_fit=  glmnet::glmnet(x=geno_omic_object,
                                   y=pheno_object[, response],
                                   alpha = alpha,
                                   standardize = FALSE,
                                   lambda =AI_fit_CV$lambda.min)



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

         AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                      y=pheno_object[, response],
                                      nfolds = 5,
                                      alpha = alpha,
                                      standardize = FALSE,
                                      lambda = lambda)

         AI_fit=  glmnet::glmnet(x=geno_omic_object,
                                 y=pheno_object[, response],
                                 alpha = alpha,
                                 standardize = FALSE,
                                 lambda =AI_fit_CV$lambda.min)


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

           AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                        y=pheno_object[, response],
                                        nfolds = 5,
                                        alpha = alpha,
                                        standardize = FALSE,
                                        lambda = lambda)

           AI_fit=  glmnet::glmnet(x=geno_omic_object,
                                   y=pheno_object[, response],
                                   alpha = alpha,
                                   standardize = FALSE,
                                   lambda =AI_fit_CV$lambda.min)



           AI_preds <- stats::predict(AI_fit,
                                      geno_omic_object,
                                      reshape = TRUE)

           AI_preds <- as.data.frame(AI_preds)

           names(AI_preds) = response



         }

       }


     } ## End of when no need for tunning.





     model_para <- c(lambda = AI_fit$lambda,
                     alpha = alpha
     )


     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")


}


return(output)



 }
