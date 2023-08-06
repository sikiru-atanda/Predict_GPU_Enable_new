
#' Title
#'
#' @param gen_name
#' @param core
#' @param message
#' @param center
#' @param xgb_paras_tunning
#' @param learning_rate
#' @param max_depth
#' @param subsample
#' @param booster
#' @param iteration
#' @param para_tunning
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_Xgb <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   xgb_paras_tunning= c(Iter_tune = NULL, # number of boosting iterations
                                        learning_rate_tune = NULL, # learning rate, low value means model is more robust to overfitting
                                        L2_tune = NULL, # L2 Regularization (Ridge Regression)
                                        L1_tune = NULL),
                   learning_rate = 0.001,
                   max_depth = 6,
                   subsample = 0.5,
                   booster = "gblinear",
                   iteration = 5000,
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
    xgb_grid = expand.grid(nrounds = para_tunning$Iter_tune , # number of boosting iterations
                           eta = para_tunning$learning_rate_tune, # learning rate, low value means model is more robust to overfitting
                           lambda = para_tunning$L2_tune, # L2 Regularization (Ridge Regression)
                           alpha = para_tunning$L1_tune # L1 Regularization (Lasso Regression)
    )



    #here we do one better then a validation set, we use cross validation to
    #expand the amount of info we have!

    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    xgb_trcontrol = caret::trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)




    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


      xgb_fit = caret::train(x = geno_omic_object,
                             y = pheno_object[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object)


    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                            label = pheno_object[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) = response[trait]

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object)

      }

    }


    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
  } else {

    xgb_params <- list(
      booster = booster,
      eta = learning_rate,
      max_depth = max_depth, #This indicates how deep the built tree can be.
      #The deeper the tree, the more splits it has and it captures more
      #information about how the data. We fit a decision tree with depths
      #ranging from 1 to 32 and plot the training and test errors
      gamma = 4,
      subsample = subsample,
      colsample_bytree = 1,
      objective = "reg:squarederror",
      eval_metric = c("rmse", "rmsle", "mape")
    )





    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


      geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                     label = pheno_object[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic_object,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic_object)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                          label = pheno_object[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic_object,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_object,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object)

      }

    }


  } ## End of when no need for tunning.


     model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

output = list(model_para,
              xgb_preds,
              res_feature,
              xgb_fit)


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
xgb_grid = expand.grid(nrounds = xgb_paras_tunning$Iter_tune , # number of boosting iterations
                         eta = xgb_paras_tunning$learning_rate_tune, # learning rate, low value means model is more robust to overfitting
                         lambda = xgb_paras_tunning$L2_tune, # L2 Regularization (Ridge Regression)
                         alpha = xgb_paras_tunning$L1_tune # L1 Regularization (Lasso Regression)
                         )



#here we do one better then a validation set, we use cross validation to
#expand the amount of info we have!

    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    xgb_trcontrol = caret::trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)




    if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {



   xgb_fit = caret::train(x = geno_omic_object,
                   y = pheno_object[, response],
                   trControl = xgb_trcontrol,
                   tuneGrid = xgb_grid,
                   method = "xgbLinear")



   xgb_preds <- stats::predict(xgb_fit,
                               geno_omic_test_object,
                        reshape = TRUE)

   xgb_preds <- as.data.frame(xgb_preds)

   names(xgb_preds) = response

   bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

   res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                   X_train = geno_omic_object)


    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                 label = pheno_object[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) = response

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object)

      }

    }


    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
    } else {

    xgb_params <- list(
      booster = booster,
      eta = learning_rate,
      max_depth = max_depth, #This indicates how deep the built tree can be.
      #The deeper the tree, the more splits it has and it captures more
      #information about how the data. We fit a decision tree with depths
      #ranging from 1 to 32 and plot the training and test errors
      gamma = 4,
      subsample = subsample,
      colsample_bytree = 1,
      objective = "reg:squarederror",
      eval_metric = c("rmse", "rmsle", "mape")
    )





    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {


      geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                     label = pheno[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic_object,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic_test_object,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                 label = pheno_object[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) <-  response

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object)

      }

    }




    } ## End of when no need for tunning.





      #}

  model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

output = list(model_para,
              xgb_preds,
              res_feature,
              xgb_fit)


names(output) = c("model_parameters", "predicted_values",
                  "feature_weight", "trained_model")


}


return(output)



 }
