
#' Title
#'
#' @param gen_name
#' @param core
#' @param message
#' @param para_tunning
#' @param pheno_object
#' @param geno_omic_object
#' @param geno_omic_test_object
#' @param response
#' @param ...
#' @param lasso_paras_tunning
#' @param lambda
#' @param GS_model
#' @param scale
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
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   AI_cv_nfolds = 5,
                   para_tunning = FALSE,
                   lasso_paras_tunning= list(lambda_tune=seq(0.000001,0.9,length.out=100)^4),
                   lambda_rr = NULL,
                   GS_model = c("Lasso",
                              "Ridge_Regression"),
                   ...

){

  if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
    if(isTRUE(scaling) || isFALSE(scaling)){
      geno_omic_object <- scale(geno_omic_object, center = TRUE, scale = TRUE)
    }

  }
  if(!is.null(geno_omic_test_object)){
    if(isTRUE(scaling) || isFALSE(scaling)){
      geno_omic_test_object <- scale(geno_omic_test_object, center = TRUE, scale = TRUE)
    }

  }
  #msg <- sprintf("==================================================\n")
  #####################################################################
  if (is.null(lambda_rr) & !inherits(lambda_rr, 'numeric')){

    lambda <-  seq(0.000001,0.9,length.out=100)^4
  }

  if(GS_model=="Ridge_Regression"){

    alpha <-  0
  }

  if(GS_model=="Lasso"){

    alpha <-  1
  }


## when length of response variable is 1
#########################

     if(isTRUE(para_tunning)){
       # create hyperparameter grid
       AI_grid <- expand.grid(
         #k = data.frame(k = seq(11,85,by = 2)))
         lambda = data.frame(alpha = alpha, lambda = lambda_tune)

       )

       AI_trcontrol <- caret::trainControl(method = "cv",
                                          number = AI_cv_nfolds,
                                          verboseIter = TRUE,
                                          returnData = FALSE,
                                          returnResamp = "all",
                                          allowParallel = TRUE)

       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         AI_fit <-  caret::train(x = geno_omic_object,
                               y = pheno_object[, response],
                               trControl = AI_trcontrol,
                               tuneGrid = AI_grid,
                               GS_model = "glment")

         if(!is.null(geno_omic_test_object)){
           GID <- rownames(geno_omic_test_object)
         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         } else {

           if(!is.null(geno_omic_object) & !is.null(geno_omic_test_object)){
             AI_preds <- stats::predict(AI_fit,
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

         bestTune <- c(AI_fit$bestTune, AI_fit$GS_model)


       } else {
         stop("Training set missing")
       }
       ####### Start when their is no need for tunning
     } else {


       if(!is.null(geno_omic_object)) {
         GID <- rownames(geno_omic_object)
         AI_fit_CV<-glmnet::cv.glmnet(x=geno_omic_object,
                                      y=pheno_object[, response],
                                      nfolds = AI_cv_nfolds,
                                      alpha = alpha,
                                      standardize = FALSE,
                                      lambda = lambda)

         AI_fit <-   glmnet::glmnet(x=geno_omic_object,
                                   y=pheno_object[, response],
                                   alpha = alpha,
                                   standardize = FALSE,
                                   lambda =AI_fit_CV$lambda.min)

         if(!is.null(geno_omic_test_object)){
           GID <- rownames(geno_omic_test_object)

         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         } else {

           if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
             AI_preds <- stats::predict(AI_fit,
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



       }


     } ## End of when no need for tunning.





     model_para <- data.frame(stat = c("lambda", "alpha"),
                              summary = c(AI_fit$lambda,alpha),
                              stringsAsFactors = FALSE)

     colnames(model_para)[1:2] <- c("stat", "summary")

     output <-  list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) <-  c("model_parameters",
                       "predicted_values",
                        "trained_model")

return(output)

 }
