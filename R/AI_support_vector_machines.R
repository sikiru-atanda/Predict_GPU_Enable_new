
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
                   scale = TRUE,
                   para_tunning = FALSE,
                   AI_cv_nfolds = 5,
                   svm_paras_tunning= c(c_tune=NULL),
                   c = 1,
                   ...

){

  #msg <- sprintf("==================================================\n")
     if(isTRUE(para_tunning)){
       # create hyperparameter grid

       AI_grid <- expand.grid(
         #k = data.frame(k = seq(11,85,by = 2)))

         C = data.frame(C = c_tune)

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
                               scaled  = FALSE,
                               method = "svmLinear")

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
         names(AI_preds)[1] = c(gen_name)

         bestTune <- c(AI_fit$bestTune, AI_fit$method)

       } else {

         stop("Training set missing")
       }

       ####### Start when their is no need for tunning

     } else {


       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         AI_fit = kernlab::ksvm(x = geno_omic_object,
                                y = pheno_object[, response],
                                scaled  = FALSE,
                                C = c
         )

         if(!is.null(geno_omic_test_object)) {
           GID <- rownames(geno_omic_test_object)
           AI_preds <- kernlab::predict(AI_fit,
                                        geno_omic_test_object)
         } else {
           if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
             AI_preds <- kernlab::predict(AI_fit,
                                          geno_omic_test_object)
           }

         }

         AI_preds <- data.frame(name = GID,
                                Predicted_value = as.data.frame(AI_preds),
                                Standard_error = NA,
                                PEV = NA,
                                Reliability = NA,
                                stringsAsFactors = FALSE)
         names(AI_preds)[1] = c(gen_name)

       }

       model_para <- data.frame(parameters = c("epsilon", "C"),
                                value = c(AI_fit@param$epsilon,AI_fit@param$C)
       )

     } ## End of when no need for tunning.


     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")

return(output)


 }



