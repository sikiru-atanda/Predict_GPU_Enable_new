
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
                   AI_cv_nfolds = 5,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   knn_paras_tunning= c(k_tune=NULL),
                   k = 5,
                   ...

){

  ## KNN relies on distance metrics to find the nearest neighbors,
  ## so scaling the features is critical to ensure that
  ## all dimensions contribute equally to the distance calculations.
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

     if(isTRUE(para_tunning)){
       # create hyperparameter grid

       AI_grid <- expand.grid(
         #k = data.frame(k = seq(11,85,by = 2)))

         k = data.frame(k = k_tune)

       )

       #here we do one better then a validation set, we use cross validation to
       #expand the amount of info we have!

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
                               method = "knn")

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




       } else {
         stop("Training set missing")
       }

       ####### Start when their is no need for tunning

     } else {


       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         AI_fit = caret::knnreg(x = geno_omic_object,
                                y = pheno_object[, response],
                                k = k
         )


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

       } else {
         stop("Training set missing")
       }


     } ## End of when no need for tunning.

     model_para <- data.frame(stat ="k",
                              summary = AI_fit$k)
     colnames(model_para)[1:2] <- c("stat", "summary")
     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")

return(output)



 }
