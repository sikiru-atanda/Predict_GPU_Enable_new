
# Function to extract matching variance-covariance structures in asreml model
#' Title
#'
#' @param formula
#' @param var_cov_structures
#'
#' @return
#' @export
#'
#' @examples
extract_var_cov_structures <- function(formula, var_cov_structures) {
  # Use a regular expression to find matches
  pattern <- paste0("\\b(", paste(var_cov_structures, collapse = "|"), ")\\b")
  matches <- regmatches(formula, gregexpr(pattern, formula))
  unique(unlist(matches))
}


#' Title
#'
#' @param y
#' @param ETA
#' @param weights
#' @param bayes_para
#' @param tst
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param pheno_data
#' @param gen_name
#' @param heter_groups
#' @param asreml_models_prep_cv
#' @param response
#' @param tst
#'
#' @return
#' @export
#'
#' @examples
asreml_mod_cv <- function(pheno_data,
                          gen_name,
                          heter_groups,
                          #var_cov_str,
                          asreml_models_prep_cv,
                          response,
                          tst){

  #browser()
asreml_tst_model_cv <- asreml_cv_model(pheno_dataa = pheno_data,
                                       response = response,
                                       #heter_groups = heter_groups,
                                       asreml_models_prep_cv = asreml_models_prep_cv,
                                       tst = tst)

var_cov_str_available <- c("us","corgh","corgv",
                           "corh","corv","fa","rr")
code_asr_fit <- asreml_models_prep_cv$code_asr_fit

# Find the index where 'random' appears in the vector
random_index <- grep("random", code_asr_fit, ignore.case = TRUE)

# Check if 'random' term exists and extract it
if (length(random_index) > 0) {
  # Extract the random part of the model
  random_formula <- asreml_models_prep_cv$code_asr_fit[random_index]
  # Apply function and extract matches for the variance-covariance structure
  var_cov_str <- extract_var_cov_structures(random_formula, var_cov_str_available)

  ## This is when var_cov_str return 0 rather than NULL when
  if(length(var_cov_str)==0) var_cov_str <- NULL

} else {
  var_cov_str <- NULL
}

gen_pos <-  asreml_models_prep_cv[["gen_pos"]]
inter_gen_pos <-  asreml_models_prep_cv[["inter_gen_pos"]]
names_in_inv_list <-  asreml_models_prep_cv[["names_in_inv_list"]]
rand_term <-  asreml_models_prep_cv[["rand_term"]]

BLUP <- summary(asreml_tst_model_cv[["model_cv"]], coef=TRUE)$coef.random

colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

#heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

if (!is.null(heter_groups)){
  ### It possible the user provide the heter_groups while it actually a single environment,
  ## This will check and turn it off
  if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
    heter_groups <-  NULL
    inter_gen_pos <-  NULL
  } else{
    if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      heter_grp <- as.character(unique(pheno_data[, heter_groups]))
      all_envs_for_met <- as.character(pheno_data[, heter_groups])
    }
  }

}

##############################################
#### Extract Breeding values/genetic effect estimate for all omics
##############################################
estimated_breeding_value_list <- list()

##
for (b in seq_along(names_in_inv_list)) {
  estimated_breeding_value_list[[names_in_inv_list[b]]] <- BLUP[grep(paste(names_in_inv_list[b],"\\)", sep = ""),rownames(BLUP)),]

}
###
### For variance structure extraction
if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

  if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
    ## Extract the number of factors
    #N_fa = substr(var_cov_str, 3, 100)

    for (bb in seq_along(names_in_inv_list)) {
      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][!rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])%in%rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]][grep('Comp',rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])),]), ]

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

      estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,3])

    }

  } else {

    if (var_cov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {

      for (bb in seq_along(names_in_inv_list)) {

        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

        estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,2])

      }


    }


  } ## End

  #################################
  if(!is.null(inter_gen_pos)){
    for (bb in seq_along(names_in_inv_list)) {

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(4, 1:2)]

      estimated_breeding_value_list[[names_in_inv_list[bb]]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name])))

      estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(1, 4, 2:3)]

      colnames(estimated_breeding_value_list[[names_in_inv_list[bb]]])[1:3] <- c(gen_name, heter_groups, "BLUP")

      #estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Standard_error"]^2
      rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]) = NULL


    }

  }

  # Bind all data frames into a single data frame
  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Group by the user-specified environment/location and sum the BLUP values
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))
  summarized_blup <- as.data.frame(summarized_blup)
  return(summarized_blup[tst, "Summed_BLUP"])

}

### For compound symmetric
if(is.null(var_cov_str) & !is.null(inter_gen_pos)){
  for (bb in seq_along(names_in_inv_list)) {
    estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
    estimated_breeding_value_list[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[bb]]), "\\)_", 3)[,2])
    rownames(estimated_breeding_value_list[[bb]]) <-  NULL
  }
  #################################
  if(!is.null(inter_gen_pos)){
    for (bb in seq_along(names_in_inv_list)) {
      estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(4, 1:2)]
      estimated_breeding_value_list[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list[[bb]][, gen_name])))
      estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(1, 4, 2:3)]
      colnames(estimated_breeding_value_list[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
      #estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[bb]][, "Standard_error"]^2
      rownames(estimated_breeding_value_list[[bb]]) <- NULL
    }
  }

  # Bind all data frames into a single data frame
  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Group by the user-specified environment/location and sum the BLUP values
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))

  summarized_blup <- as.data.frame(summarized_blup)
  return(summarized_blup[tst, "Summed_BLUP"])

}

if(is.null(var_cov_str) & is.null(inter_gen_pos) ){
  for (bb in seq_along(names_in_inv_list)){
    estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
    estimated_breeding_value_list[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[bb]]), "\\)_", 3)[,2])
    rownames(estimated_breeding_value_list[[bb]]) <-  NULL

    # if(is.null(estimated_breeding_value_list[[bb]])){
    #   stop("it null")
    # }
    #colnames( estimated_breeding_value_list[[bb]])
    #message(gen_name)
    #message(colnames(estimated_breeding_value_list[[bb]]))
    estimated_breeding_value_list[[bb]] <-   estimated_breeding_value_list[[bb]][, c(4, 1:2)]
    colnames(estimated_breeding_value_list[[bb]])[1:2] <- c(gen_name, "BLUP")
  }

  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Group by GID and sum the BLUP values for each GID
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))

  summarized_blup <- as.data.frame(summarized_blup)
  return(summarized_blup[tst, "Summed_BLUP"])

}





}

#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param eta
#' @param nrounds
#' @param max_depth
#' @param scaling
#' @param centering
#' @param omic_count
#' @param gamma
#' @param subsample
#' @param colsample_bytree
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param ncomp
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param ntree
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#' @param k
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param y
#' @param omics
#' @param tst
#' @param scaling
#' @param centering
#' @param omic_count
#' @param c
#'
#' @return
#' @export
#'
#' @examples
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
