#'
#' # Function to dynamically select the best model based on available metrics
#' #' Title
#' #'
#' #' @param model_list
#' #' @param preferred_metrics
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' select_best_model <- function(model_list,
#'                               preferred_metrics=c("RMSE", "Accuracy", "MAE")) {
#'
#'   scores <- lapply(model_list, function(model) {
#'     available_metrics <- intersect(names(model$results), preferred_metrics)
#'     if (length(available_metrics) > 0) {
#'       # Assuming lower values are better for selected metrics, flip if needed
#'       metric_scores <- sapply(available_metrics, function(metric) {
#'         if (metric == "Accuracy") {
#'           max(model$results[[metric]])  # Higher is better
#'         } else if(metric == "MAE"){
#'           min(model$results[[metric]])
#'         }else {
#'           if(metric == "RMSE"){
#'             min(model$results[[metric]])
#'           }
#'         }
#'         # Lower is better for RMSE
#'       })
#'       c(BestScore=min(metric_scores), BestMetric=names(metric_scores)[which.min(metric_scores)])
#'     } else {
#'       c(BestScore=NA, BestMetric=NA)
#'     }
#'   })
#'
#'   # Remove models without any of the preferred metrics
#'   scores <- scores[!sapply(scores, function(x) is.na(x["BestScore"]))]
#'
#'   # Identify the model with the best overall score
#'   best_model_index <- which.min(sapply(scores, function(x) x["BestScore"]))
#'   best_model_name <- names(scores)[best_model_index]
#'
#'   list(BestModel=model_list[[best_model_name]], BestModelName=best_model_name, Metrics=scores[[best_model_index]])
#' }


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
AI_svmOLD <- function(pheno_object=NULL,
                   geno_omic_object = NULL,
                   geno_omic_test_object = NULL,
                   response=NULL,
                   gen_name=NULL,
                   message = TRUE,
                   scaling = TRUE,
                   centering = FALSE,
                   omic_count = NULL,
                   para_tunning = FALSE,
                   AI_cv_nfolds = 5,
                   svm_paras_tunning= list(
                     kernel = c("radial", "linear", "polynomial"),
                     #cost = 10^seq(-2, 2, by = 1),
                     sigma = c(0.01, 0.05, 0.1),
                     C = c(1, 10, 100), ## for radial kernel
                     degree = c(3, 4),  # Default values, used only for polynomial
                     scale = c(0.1, 1) # used only for polynomial
                   ),
                   svm_kernel = "Gaussian", # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
                   sigma_value  = 0.1,       # Default sigma value for RBF kernel
                   C_value  = 1,             # Default cost parameter
                   degree_value = 3,        # Default degree for polynomial kernel
                   scale_value  = 1,         # Default scale for polynomial kernel
                   offset_value = 1,        # Default offset for polynomial kernel
                   ...

){

  # Translate user-friendly kernel names to `kernlab` kernel function names
  kernel_type <- switch(svm_kernel,
                        Gaussian = "rbfdot",         # Radial Basis Function kernel
                        Polynomial = "polydot",      # Polynomial kernel
                        Linear = "vanilladot",       # Linear kernel
                        Hyperbolic_tangent = "tanhdot"  # Sigmoid kernel
  )
  # Create a list to store kernel-specific parameters
  kernel_params <- list()

  # Set kernel parameters based on user input or defaults
  switch(kernel_type,
         rbfdot = {kernel_params <- list(sigma = sigma_value)},
         polydot = {kernel_params <- list(degree = degree_value, scale = scale_value, offset = offset_value)},
         #vanilladot = {kernel_params <- list(C = C_value)},  # Linear kernel
         tanhdot = {kernel_params <- list(scale = scale_value, offset = offset_value)}  # Sigmoid kernel
  )

  ## SVM performance is sensitive to the scale of the input features,
  ## so scaling the data is recommended to ensure that all features contribute
  ## equally to the distance metric in the kernel space
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
  #msg <- "\n==================================================\n"
     if(isTRUE(para_tunning)){
       # create hyperparameter grid

       # Include gamma for kernels that need it (not needed for linear unless specified)
       if ("radial" %in% svm_paras_tunning$kernel) {

         tuning_grid_radial <- expand.grid(sigma = svm_paras_tunning$sigma,
                                           C = svm_paras_tunning$C)
       }
       ##
       if ("polynomial" %in% svm_paras_tunning$kernel) {
         tuning_grid_poly <- expand.grid(degree = svm_paras_tunning$degree,
                                         scale = svm_paras_tunning$scale,
                                         C = svm_paras_tunning$C)
       }
       ##
       if ("linear" %in% svm_paras_tunning$kernel) {
         tuning_grid_linear <- expand.grid(C = data.frame(C =svm_paras_tunning$C))
       }
###
       # List of SVM models to evaluate
       svm_models <- list(
         svmRadial = list(method = "svmRadial", tuneGrid = tuning_grid_radial),
         svmPoly = list(method = "svmPoly", tuneGrid = tuning_grid_poly),
         svmLinear = list(method = "svmLinear", tuneGrid = tuning_grid_linear)
       )

       AI_trcontrol = caret::trainControl(method = "cv",
                                          number = AI_cv_nfolds,
                                          verboseIter = TRUE,
                                          returnData = FALSE,
                                          returnResamp = "all",
                                          allowParallel = TRUE)

       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         # Function to train and evaluate models
         train_and_evaluate <- function(model_info, data, labels, control) {
           caret::train(
                       x = data,
                       y = labels,
                       method = model_info$method,
                       trControl = control,
                       tuneGrid = model_info$tuneGrid
                       #preProcess = c("center", "scale")
                     )
         }

         # Train and evaluate all models
##lapply is used to apply the train_and_evaluate function to each element in the svm_models list
## data, label and control, they are passed explicitly to train_and_evaluate through lapply which
## are needed for training the model
         results <- lapply(svm_models,
                           train_and_evaluate,
                           data = geno_omic_object,
                           labels = pheno_object[, response],
                           control = AI_trcontrol)

         # Apply the function to the results
         best_model_info <- select_best_model(results)

         AI_fit <- best_model_info$BestModel
         # E
         # AI_fit = caret::train(x = geno_omic_object,
         #                       y = pheno_object[, response],
         #                       trControl = AI_trcontrol,
         #                       tuneGrid = AI_grid,
         #                       scaled  = FALSE,
         #                       #method = "svmLinear"
         #                       method = "svm"
         #                       )

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

         #bestTune <- c(AI_fit$bestTune, AI_fit$method)

       } else {

         stop("Training set missing")
       }

       ####### Start when their is no need for tunning

     } else {

       if(!is.null(geno_omic_object) & !is.null(pheno_object)) {
         GID <- rownames(geno_omic_object)
         if(kernel_type!="vanilladot"){
         AI_fit = kernlab::ksvm(x = geno_omic_object,
                                y = pheno_object[, response],
                                kernel = kernel_type,
                                scaled = FALSE,
                                type = "nu-svr",
                                C = C_value,
                                kpar = kernel_params)
         } else {
           if(kernel_type=="vanilladot"){
             AI_fit = kernlab::ksvm(x = geno_omic_object,
                                    y = pheno_object[, response],
                                    kernel = kernel_type,
                                    scaled = FALSE,
                                    type = "nu-svr",
                                    C = C_value
                                   )
           }

         }


         if(!is.null(geno_omic_test_object)) {
           GID <- rownames(geno_omic_test_object)
           AI_preds <- kernlab::predict(AI_fit,
                                        geno_omic_test_object)
         } else {
           if(!is.null(geno_omic_object) & is.null(geno_omic_test_object)) {
             AI_preds <- kernlab::predict(AI_fit,
                                          geno_omic_object)
           }

         }

         AI_preds <- data.frame(name = GID,
                                Predicted_value = as.data.frame(AI_preds),
                                Standard_error = NA,
                                PEV = NA,
                                Reliability = NA,
                                stringsAsFactors = FALSE)
         names(AI_preds)[1:2] = c(gen_name, "Predicted_value")

       }


     } ## End of when no need for tunning.


  # model_para <- data.frame(stat = c("epsilon", "C"),
  #                          summary = c(AI_fit@param$epsilon,AI_fit@param$C)
  # )
  # colnames(model_para)[1:2] <- c("stat", "summary")

  model_para <- data.frame(stat = c("epsilon", "C"),
                           #summary = c(AI_fit@param$epsilon,AI_fit@param$C)
                           summary = c(AI_fit@param$epsilon,0.1)
  )
  colnames(model_para)[1:2] <- c("stat", "summary")

     output = list(model_para,
                   AI_preds,
                   AI_fit)


     names(output) = c("model_parameters",
                       "predicted_values",
                        "trained_model")

return(output)


 }



