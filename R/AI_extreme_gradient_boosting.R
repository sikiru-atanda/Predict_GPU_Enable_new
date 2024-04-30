
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
#' @param scale
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
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   xgb_paras_tunning= list(Iter_tune = seq(500, 5000, 500), # number of boosting iterations
                                        learning_rate_tune = c(0.01, 0.05, 0.1), # learning rate, low value means model is more robust to overfitting
                                        max_depth = c(3, 6, 9),
                                        gamma = c(0, 0.01, 0.1),
                                        colsample_bytree = c(0.5, 0.75, 1),
                                        min_child_weight = c(1, 3, 5),
                                        subsample = c(0.5, 0.75, 1),
                                        L2_tune = c(0, 0.5, 1), #  for linear gbL2 Regularization (Ridge Regression)
                                        L1_tune = c(0, 0.5, 1)), # for linear gb
                   resample_method_tune = "cv", # c("cv","boot")
                   number_of_fold_tune = 5,
                   learning_rate = 0.001,
                   max_depth = 6,
                   subsample = 0.5,
                   xgb_booster =  "gbtree", # "gblinear",
                   iteration = 5000,
                   N_feature_impo = 10,
                   ...

){

  if(!is.null(geno_omic_object)){
    if(isTRUE(scaling)){
      geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
    } else{
      if(isTRUE(centering) && !is.null(omic_count)){
        geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = FALSE)
      }
    }
  }

  if(!is.null(geno_omic_test_object)){
    if(isTRUE(scaling)){
      geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
    } else{
      if(isTRUE(centering) && !is.null(omic_count)){
        geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = FALSE)
      }
    }
  }
  #msg <- sprintf("==================================================\n")
## when length of response variable is 1
#########################

if(isTRUE(para_tunning)){

  xgb_grid_tree <- NULL

  xgb_grid_linear <- NULL

  xgb_trcontrol = caret::trainControl(method = "cv",
                                      number = AI_cv_nfolds,
                                      verboseIter = TRUE,
                                      returnData = FALSE,
                                      returnResamp = "all",
                                      allowParallel = TRUE)

  if (xgb_booster == "gbtree") {
    # Define the parameter grid for a tree-based model
    xgb_grid_tree <-  expand.grid(
      nrounds = xgb_paras_tunning$Iter_tune,   # Number of boosting rounds
      eta = xgb_paras_tunning$learning_rate_tune,  # Learning rate
      max_depth = xgb_paras_tunning$max_depth,  # Varying tree depths
      gamma = xgb_paras_tunning$gamma,  # Minimum loss reduction required for further partition
      colsample_bytree = xgb_paras_tunning$colsample_bytree,  # Subsample ratio of columns when constructing each tree
      min_child_weight = xgb_paras_tunning$min_child_weight,  # Minimum sum of instance weight needed in a child
      subsample = xgb_paras_tunning$subsample
    )

  }

if(xgb_booster=="gblinear"){
  # Define the parameter grid for a linear model
  xgb_grid_linear <-  expand.grid(
    nrounds = xgb_paras_tunning$Iter_tune,  # Number of boosting rounds
    eta = xgb_paras_tunning$learning_rate_tune,  # Learning rate
    lambda = xgb_paras_tunning$L2_tune,  # L2 Regularization
    alpha = xgb_paras_tunning$L1_tune   # L1 Regularization
  )
}


#here we do one better then a validation set, we use cross validation to
#expand the amount of info we have!

    # xgb_trcontrol = caret::trainControl(method = "cv",
    #                                     number = AI_cv_nfolds,
    #                                     verboseIter = TRUE,
    #                                     returnData = FALSE,
    #                                     returnResamp = "all",
    #                                     allowParallel = TRUE)

    if(!is.null(geno_omic_object)){

      GID <- rownames(geno_omic_object)
      if(!is.null(xgb_grid_linear)){
   xgb_fit = caret::train(x = geno_omic_object,
                   y = pheno_object[, response],
                   trControl = xgb_trcontrol,
                   tuneGrid = xgb_grid_linear,
                   method = "xgbLinear")
      } else {
        if(!is.null(xgb_grid_tree)){
          xgb_fit = caret::train(x = geno_omic_object,
                                 y = pheno_object[, response],
                                 trControl = xgb_trcontrol,
                                 tuneGrid = xgb_grid_tree,
                                 method = "xgbTree")
        }

      }

   if(!is.null(geno_omic_test_object)) {
     GID <- rownames(geno_omic_test_object)
     AI_preds <- stats::predict(xgb_fit,
                               geno_omic_test_object,
                        reshape = TRUE)

   } else {

     if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
       AI_preds <- stats::predict(xgb_fit,
                                   geno_omic_object,
                                   reshape = TRUE)

     }
   }


   AI_preds <- data.frame(name = GID,
                          Predicted_value = as.data.frame(AI_preds),
                          Standard_error = NA,
                          PEV = NA,
                          Reliability = NA,
                          stringsAsFactors = FALSE)


   names(AI_preds)[1:2] <-  c(gen_name, "Predicted_value")


   bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

   res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                   X_train = geno_omic_object,
                                   N_feature_impo = N_feature_impo)

   res_feature <- res_feature$feature_weight
    #} else {

    } else {

      stop("Missing training set")
    }

    ####### Start when their is no need for tunning
    } else {

      if (booster == "gbtree") {
        xgb_params <- list(
          booster = booster,
          eta = learning_rate,
          max_depth = max_depth,
          gamma = 4,
          subsample = subsample,
          colsample_bytree = 1,
          objective = "reg:squarederror",
          eval_metric = c("rmse", "rmsle", "mape")
        )
      } else if (booster == "gblinear") {
        xgb_params <- list(
          booster = booster,
          eta = learning_rate,
          objective = "reg:squarederror",
          eval_metric = c("rmse", "rmsle", "mape")
        )
      }
    # xgb_params <- list(
    #   booster = booster,
    #   eta = learning_rate,
    #   max_depth = max_depth, #This indicates how deep the built tree can be.
    #   #The deeper the tree, the more splits it has and it captures more
    #   #information about how the data. We fit a decision tree with depths
    #   #ranging from 1 to 32 and plot the training and test errors
    #   gamma = 4,
    #   subsample = subsample,
    #   colsample_bytree = 1,
    #   objective = "reg:squarederror",
    #   eval_metric = c("rmse", "rmsle", "mape")
    # )


    if(!is.null(geno_omic_object) & !is.null(pheno_object)) {

      GID <- rownames(geno_omic_object)
      geno_omic_object <- xgboost::xgb.DMatrix(data = geno_omic_object,
                                              label = pheno_object[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic_object,
        nrounds = iteration,
        verbose = 0
      )
if(!is.null(geno_omic_test_object)){
  GID <- rownames(geno_omic_test_object)
  AI_preds <- stats::predict(xgb_fit,
                           geno_omic_test_object,
                           reshape = TRUE)

} else {
  if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){

    AI_preds <- stats::predict(xgb_fit,
                                geno_omic_object,
                                reshape = TRUE)


  }

}

      AI_preds <- data.frame(name = GID,
                             Predicted_value = as.data.frame(AI_preds),
                             Standard_error = NA,
                             PEV = NA,
                             Reliability = NA,
                             stringsAsFactors = FALSE)

#View(AI_preds)
      names(AI_preds)[1:2] <-  c(gen_name, "Predicted_value")

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic_object,
                                      N_feature_impo = N_feature_impo)


      res_feature <- res_feature$feature_weight
    } else {

      stop("training data missing")
    }

    } ## End of when no need for tunning.

      #}

  model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

  model_para <- as.data.frame(unlist(model_para))
  model_para <-  data.frame(stat = rownames(model_para),
                            summary = model_para,
                  stringsAsFactors = FALSE)
  colnames(model_para)[1:2] <- c("stat", "summary")

  rownames(model_para) <- NULL

output = list(model_para,
              AI_preds,
              res_feature,
              xgb_fit)


names(output) = c("model_parameters", "predicted_values",
                  "feature_weight", "trained_model")

return(output)


 }
