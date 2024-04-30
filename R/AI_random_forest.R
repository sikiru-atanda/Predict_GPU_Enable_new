
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
#' @param rf_paras_tunning
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
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   rf_paras_tunning= list(mtry = TRUE,
                                          ntree = c(500, 1000, 1500),
                                          nodesize = c(1, 5, 10),
                                          maxnodes = c(30, 50, NULL)),  # NULL means no limit),
                   ntree=500,
                   mtry = NULL,
                   maxnodes = NULL,
                   importance=TRUE,
                   ...

){

  #msg <- sprintf("==================================================\n")

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
#########################

     if(isTRUE(para_tunning)){
       # create hyperparameter grid

       n_features = ncol(geno_omic_object)

       if(isTRUE(mtry)){
         mtry <- c(sqrt(ncol(geno_omic_object)), sqrt(ncol(geno_omic_object))/2, ncol(geno_omic_object)/3)
       } else {
         stop("set mtry_true")
       }

       AI_grid <- expand.grid(
         mtry = mtry,
         ntree = para_tunning$ntree,
         nodesize = para_tunning$nodesize,
         maxnodes = para_tunning$maxnodes  # NULL means no limit
       )

       AI_trcontrol = caret::trainControl(method = "cv",
                                          number = AI_cv_nfolds,
                                          verboseIter = TRUE,
                                          returnData = FALSE,
                                          returnResamp = "all",
                                          allowParallel = TRUE)

       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         AI_fit = caret::train(x = geno_omic_object,
                               y = pheno_object[, response],
                               trControl = AI_trcontrol,
                               tuneGrid = AI_grid,
                               method = "rf")

         if(!is.null(geno_omic_test_object)){
           GID <- rownames(geno_omic_test_object)
         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)

         } else {

           if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){

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


         names(AI_preds)[1:2] = c(gen_name, "Predicted_value")


         bestTune <- c(AI_fit$bestTune, AI_fit$method)

         res_feature = caret::varImp(AI_fit)[[1]]
         res_feature = res_feature[order(res_feature[, 1], decreasing = T)[1:20], drop=F,]
         res_feature$feature = rownames(res_feature)




       } else {

         stop("Training set missing")
       }

       ####### Start when their is no need for tunning

     } else {


       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {

         GID <- rownames(geno_omic_object)
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

if(!is.null(geno_omic_test_object)){
  GID <- rownames(geno_omic_test_object)
         AI_preds <- stats::predict(AI_fit,
                                    geno_omic_test_object,
                                    reshape = TRUE)
} else {
  if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)){

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
         names(AI_preds)[1:2] = c(gen_name, "Predicted_value")

         res_feature = caret::varImp(AI_fit)[[1]]
         res_feature = data.frame(feature = colnames(geno_omic_object),
                                   weight = res_feature)
         res_feature = res_feature[order(res_feature$weight, decreasing = T)[1:20], drop=F,]


       } else {

         stop('Training set missing')
       }


     } ## End of when no need for tunning.


     model_para <- data.frame(stat = c("ntree","mtry"),
                              summary = c(AI_fit$ntree, AI_fit$mtry),
                              stringsAsFactors = FALSE)

     colnames(model_para)[1:2] <- c("stat", "summary")
     output = list(model_para,
                   AI_preds,
                   res_feature,
                   AI_fit)


     names(output) = c("model_parameters", "predicted_values",
                       "feature_weight", "trained_model")


return(output)

 }
