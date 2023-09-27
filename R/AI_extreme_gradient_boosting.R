
#' Title Extreme Gradient Boosting Machine Learning Genomic Selection Pipeline
#' Hyper-parameter tunning  of the parameters is allowed if desired by user
#' The parameters are:
#' Iter_tune : number of boosting iterations
#' learnining_rate_tune:  simply means how fast the model learns.
#' Each tree added modifies the overall model.
#' The magnitude of the modification is controlled by learning rate.
#' The lower the learning rate, the slower the model learns.
#' The advantage of slower learning rate is that the model becomes more robust
#' and efficient.
#'
#'
#' @param gen_name column name containing individuals/genotypes
#' @param message
#' @param center standardization of the x-variables
#' @param xgb_paras_tunning parameter to for tunning
#' @param learning_rate how slow/fast the model learn
#' @param max_depth max_depth refers to the number of leaves of each tree
#' @param subsample  This help to reduce the correlation between results from individual learners.
#' @param booster  to determine if the model is for regression or classification problem
#' @param iteration Number of iteration
#' @param para_tunning  if user required parameter tunning
#' @param pheno_object phenotypic object NA is allowed
#' @param geno_omic_object multi-omic data, NA not allowed
#' @param geno_omic_test_object multi-omic data for testing set if not present in geno_omic_object
#' @param response  y variables/lables
#' @param core number of ram for paralllel job
#' @param resample_method_tune
#' @param number_of_fold_tune
#' @param N_feature_impo number of feature/ x variables to extract based on the importance/weight
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
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   xgb_paras_tunning= c(Iter_tune = NULL, # number of boosting iterations
                                        learning_rate_tune = NULL, # learning rate, low value means model is more robust to overfitting
                                        L2_tune = NULL, # L2 Regularization (Ridge Regression)
                                        L1_tune = NULL),
                   resample_method_tune = "cv", # c("cv","boot")
                   number_of_fold_tune = 5,
                   learning_rate = 0.001,
                   max_depth = 6,
                   subsample = 0.5,
                   booster = "gblinear",
                   iteration = 5000,
                   N_feature_impo = 10,
                   core = NULL,
                   ...

){

  #msg <- sprintf("==================================================\n")
  ##############################################################
  ###################################################################
  ###  Start ML Analysis
  ###
 ######################################################################
  #####################################################################


# #### Initializing parallel for multiple response
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

#### if user interested in tunning the parameters
  if(isTRUE(para_tunning)){
    xgb_grid = expand.grid(nrounds = para_tunning$Iter_tune , # number of boosting iterations
                           eta = para_tunning$learning_rate_tune, # learning rate, low value means model is more robust to overfitting
                           lambda = para_tunning$L2_tune, # L2 Regularization (Ridge Regression)
                           alpha = para_tunning$L1_tune # L1 Regularization (Lasso Regression)
    )




    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    ### This function is used to specify the parameters for training using caret
    xgb_trcontrol = caret::trainControl(method = resample_method_tune,
                                        number = number_of_fold_tune,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)



## Here the user provide geno_omic_object as training set and
## and geno_omic_test_object as testing set. Thus the hyperparameter tunning
## is done with the geno_omic_object (training set)

    if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)) {

##  By default, the train function chooses the model with the largest
## performance value /or smallest, for mean squared error in regression models.
      xgb_fit = caret::train(x = geno_omic_object,
                             y = pheno_object[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")


## The train model based on the best hyperparameters is used for prediction
      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)
## convert the predicted value to dataframe
      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

## Extract the best hyperparamter values for META data purpose and
## subsequent prediction exercise the will make use of the training data
      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)
## Extract feature/x variables based on the importance/weight using the
## feature_impo_xgb function
      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object,
                                      N_feature_impo = N_feature_impo)


    } else {

      ## when user provide only the training set
      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                            label = pheno_object[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 0
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) = response[trait]

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object,
                                        N_feature_impo = N_feature_impo)

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
        verbose = 0
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_test_object,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object,
                                      N_feature_impo = N_feature_impo)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                          label = pheno_object[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic_object,
        nrounds = iteration,
        verbose = 0
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic_object,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response[trait]

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object,
                                      N_feature_impo = N_feature_impo)

      }

    }


  } ## End of when no need for tunning.

### Model paramerts for META data
     model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

     ## Output
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
                                   X_train = geno_omic_object,
                                   N_feature_impo = N_feature_impo)


    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                 label = pheno_object[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 0
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) = response

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object,
                                        N_feature_impo = N_feature_impo)

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
        verbose = 0
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic_test_object,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      names(xgb_preds) = response

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object,
                                      N_feature_impo = N_feature_impo)



    } else {

      if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
        geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                                 label = pheno_object[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic_object,
          nrounds = iteration,
          verbose = 0
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic_object,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        names(xgb_preds) <-  response

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic_object,
                                        N_feature_impo = N_feature_impo)

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
