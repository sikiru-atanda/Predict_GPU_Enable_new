

bayes_mod_cv <- function(y,
                        ETA,
                        weights,
                        bayes_para,
                        tst){

  fit <- BGLR::BGLR(
    y=y,
    ETA=ETA,
    weights = weights,
    nIter = bayes_para[["nIter"]],
    burnIn =  bayes_para[["burnIn"]],
    thin =  bayes_para[["thin"]],
    verbose = FALSE
    #saveAt =systime
  )

 return(fit$yHat[tst])

}

AI_xgboost_cv <- function(y,
                          omics,
                          tst,
                          eta = 0.001,
                          nrounds = 5000,
                          max_depth = 6,
                          scaling = FALSE,
                          centering = TRUE,
                          omic_count,
                          gamma = 4,
                          subsample = 0.5,
                          colsample_bytree = 1){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || !is.null(omic_count)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }
  ### Set the paramters and hyper parameters for extreme graident boosting
  suppressMessages({
  xgb_params <- list(
    booster = "gblinear",
    eta = eta,
    max_depth = max_depth, #This indicates how deep the built tree can be.
    #The deeper the tree, the more splits it has and it captures more
    #information about how the data. We fit a decision tree with depths
    #ranging from 1 to 32 and plot the training and test errors
    gamma = gamma,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    objective = "reg:squarederror",
    eval_metric = c("rmse", "rmsle", "mape")
  )

  omics_Xgb <- xgboost::xgb.DMatrix(data = omics[-tst, ],
                                    label = y[-tst])

  fit <- xgboost::xgb.train(
    params = xgb_params,
    data =  omics_Xgb,
    nrounds = nrounds,
    verbose = 0
  )

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)

  })

  preds <- as.data.frame(preds)
  return(preds[, 1])


}
######
AI_pls_cv <- function(y,
                      omics,
                      tst,
                      ncomp = 3,
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }

  pls_model <- pls::plsr(y~ omics,
                         scale = FALSE,
                         center = FALSE,
                         ncomp = ncomp,
                         validation = "none")
  cumulative_explained_variance <- cumsum(pls::explvar(pls_model))
  # Find the number of components explaining at least 90% of the variance
  num_components <- which(cumulative_explained_variance >= 90)[1]
  if(ncomp< num_components){
    ncomp <- num_components
    message(paste(msg, "The number of component provided explain less than 90% of the variance. We make adjustment as this might affect final result."))
  }

  pls_model <- pls::plsr(y[-tst]~ omics[-tst, ],
                         scale = FALSE,
                         center = FALSE,
                         ncomp = ncomp,
                         validation = "none")

  preds <- stats::predict(pls_model,
                          newdata =omics[tst, ],
                          ncomp = ncomp)

  preds <- as.data.frame(preds)
  return(preds[, 1])


}
######
AI_randomforest_cv <- function(y,
                               omics,
                               tst,
                               ntree = 500,
                               scaling = FALSE,
                               centering = TRUE,
                               omic_count){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || !is.null(omic_count)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }
  fit <- randomForest::randomForest(x = omics[-tst, ],
                                   y = y[-tst],
                                   ntree = ntree,
                                   importance = TRUE)

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)
  preds <- as.data.frame(preds)

  return(preds[, 1])


}
########
AI_ridge_regression_cv <- function(y,
                                   omics,
                                   tst,
                                   scaling = FALSE,
                                   centering = TRUE,
                                   omic_count){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }
  fit_CV<-glmnet::cv.glmnet(x= omics[-tst, ],
                            y = y[-tst],
                            #nfolds = 5,
                            alpha = 0,
                            standardize = FALSE
                            )
  # Optimal lambda for Ridge

  #lambda_1se_ridge <- cv_ridge$lambda.1se
  fit <-  glmnet::glmnet(x=omics[-tst, ],
                         y = y[-tst],
                         alpha = 0,
                         standardize = FALSE,
                         lambda =fit_CV$lambda.min)

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  return(preds[, 1])


}
######
AI_lasso_cv <- function(y,
                        omics,
                        tst,
                        scaling = FALSE,
                        centering = TRUE,
                        omic_count){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }
  fit_CV<-glmnet::cv.glmnet(x= omics[-tst, ],
                            y= y[-tst],
                            #nfolds = 5,
                            alpha = 1,
                            standardize = FALSE)

  fit <-  glmnet::glmnet(x= omics[-tst, ],
                         y = y[-tst],
                         alpha = 1,
                         standardize = FALSE,
                         lambda =fit_CV$lambda.min)

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  return(preds[, 1])


}
#####
AI_knn_cv <- function(y,
                      omics,
                      tst,
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count,
                      k = 5){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }

     fit <-  caret::knnreg(x = omics[-tst, ],
                           y = y[-tst],
                           k = k)


  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  return(preds[, 1])


}

######################
#####
AI_svm_cv <- function(y,
                      omics,
                      tst,
                      scaling = FALSE,
                      centering = TRUE,
                      omic_count,
                      c = 1){

  if(!is.null(omics)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      omics <- scale(omics, center = TRUE, scale = TRUE)
    }
  }

  fit <-  kernlab::ksvm(x = omics[-tst, ],
                        y = y[-tst],
                        scaled  = FALSE,
                        type = "nu-svr",
                        C = c)

  preds <- kernlab::predict(fit,
                            omics[tst, ])

  preds <- as.data.frame(preds)

  return(preds[, 1])


}
