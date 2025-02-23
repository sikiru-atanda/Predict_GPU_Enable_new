
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
                        tst,
                        bayes_model,
                        bayes_trait){

  systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  systime <- gsub("[-: ]", "_", systime)
  systime <- paste(bayes_model, bayes_trait, systime, sep= "_")

  fit <- BGLR::BGLR(
    y=y,
    ETA=ETA,
    weights = weights,
    nIter = bayes_para[["nIter"]],
    burnIn =  bayes_para[["burnIn"]],
    thin =  bayes_para[["thin"]],
    verbose = FALSE,
    saveAt =systime
  )

  ## Get the name of all files stored by BGLR using the current name and time the analysis was performed
  output_files_names = list.files(pattern=systime)

  unlink(output_files_names)

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
                                       gen_name = gen_name,
                                       heter_groups = heter_groups,
                                       asreml_models_prep_cv = asreml_models_prep_cv,
                                       tst = tst)

modm <-  asreml_tst_model_cv[["model_cv"]]

### Asreml required the pheno data in the environment to execute the predict func
pheno_dataa <-  asreml_tst_model_cv[["pheno_dataa"]]

if(is.null(modm)){
  return(NULL)
}
GIDs <- as.character(pheno_data[[gen_name]])

GID_tst <- GIDs[tst]

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
####
# Create a reference dataframe from pheno_dataa
if(!is.null(heter_groups)){
reference_order <- pheno_dataa |>
  dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups))

# predicted_value <- reference_order |>
#   dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
#   dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)


} else{
  reference_order <- pheno_dataa |>
    dplyr::select(!!dplyr::sym(gen_name))

  # predicted_value <- reference_order |>
  #   dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
  #   dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)

}


result_met <- tryCatch(
  {
    if (is.null(heter_groups)) {
      predicted_value <- asreml::predict.asreml(modm, classify = gen_name, sed = FALSE)$pvals
      predicted_value <- reference_order |>
        dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
        dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)

    } else {
      predicted_value <- asreml::predict.asreml(modm, classify = rand_term[[inter_gen_pos]], sed = FALSE)$pvals
      predicted_value <- reference_order |>
        dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
        dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)

    }
    predicted_value
  },
  error = function(e) {
    # Handle the error, you can print a message or take other actions
    cat(paste("Error in prediction:", conditionMessage(e), "\n"))
    return(NULL)  # Return NULL or an appropriate value to indicate the failure
  }
)


if(is.null(result_met)) {
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
      heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
      all_envs_for_met <- as.character(pheno_data[[heter_groups]])
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
      rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]) <-  NULL


    }

  }

  # Bind all data frames into a single data frame
  combined_df <- dplyr::bind_rows(estimated_breeding_value_list)

  # Group by the user-specified environment/location and sum the BLUP values
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))
  summarized_blup <- as.data.frame(summarized_blup)

  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
    dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), Summed_BLUP)

  #colnames(summarized_blup)[colnames(summarized_blup)%in%heter_groups] <- "Env"
  return(as.double(summarized_blup[tst, c("Summed_BLUP")]))


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

  if(!is.null(heter_groups)){
  #combined_df <- combined_df[order(combined_df[[heter_groups]]), ]
  #colnames(combined_df)[colnames(combined_df)%in%heter_groups] <- "Env"

  # Group by GID and sum the BLUP values for each GID
  summarized_blup <- combined_df |>
    dplyr::group_by(!!rlang::sym(gen_name)) |>
    dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))

  summarized_blup <- as.data.frame(summarized_blup)
  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
    dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), Summed_BLUP)
  #summarized_blup <- summarized_blup[as.character(summarized_blup[[gen_name]])%in%GID_tst, ]
  #summarized_blup <- summarized_blup[order(as.character(summarized_blup[[gen_name]])%in%GID_tst), ]
  #### sik
  #gen_namess <- as.character(summarized_blup[[gen_name]])
  # Get the match positions of gen_names in GID_tst
  #order_indices <- match(gen_namess, GID_tst)
  # Order the data frame based on these match positions
  #summarized_blup <- summarized_blup[order(order_indices), ]
  #colnames(summarized_blup)[colnames(summarized_blup)%in%heter_groups] <- "Env"
  return(as.double(summarized_blup[tst, c("Summed_BLUP")]))
  }


  # Group by the user-specified environment/location and sum the BLUP values
  # summarized_blup <- combined_df |>
  #   dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups)) |>
  #   dplyr::summarise(Summed_BLUP = sum(BLUP, na.rm = TRUE))
  #
  # summarized_blup <- as.data.frame(summarized_blup)
  # return(summarized_blup[tst, "Summed_BLUP"])

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

  #summarized_blup <- summarized_blup[as.character(summarized_blup[[gen_name]])%in%GID_tst, ]
  #summarized_blup <- summarized_blup[order(as.character(summarized_blup[[gen_name]])%in%GID_tst), ]
  #### sik
  #gen_namess <- as.character(summarized_blup[[gen_name]])
  # Get the match positions of gen_names in GID_tst
  #order_indices <- match(gen_namess, GID_tst)
  # Order the data frame based on these match positions
  #summarized_blup <- summarized_blup[order(order_indices), ]


  summarized_blup <- reference_order |>
    dplyr::left_join(summarized_blup, by = stats::setNames(c(gen_name), c(gen_name))) |>
    dplyr::select(!!dplyr::sym(gen_name), Summed_BLUP)

  return(as.double(summarized_blup[tst, "Summed_BLUP"]))

}

} else {

  predicted_value <- result_met

  if(!is.null(heter_groups)){
    #predicted_value <- predicted_value[order(predicted_value[[heter_groups]]), ]

    # predicted_value <- predicted_value[as.character(predicted_value[[gen_name]])%in%GID_tst, ]
    # #predicted_value <-predicted_value[order(as.character(predicted_value[[gen_name]])%in%GID_tst), ]
    # #### sik
    # gen_namess <- as.character(predicted_value[[gen_name]])
    # # Get the match positions of gen_names in GID_tst
    # order_indices <- match(gen_namess, GID_tst)
    # # Order the data frame based on these match positions
    # predicted_value <- predicted_value[order(order_indices), ]

    predicted_value <- reference_order |>
      dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name, heter_groups), c(gen_name, heter_groups))) |>
      dplyr::select(!!dplyr::sym(gen_name), !!dplyr::sym(heter_groups), predicted.value, std.error, status)

    #colnames(predicted_value)[colnames(predicted_value)%in%heter_groups] <- "Env"
    return(as.double(predicted_value[tst, c("predicted.value")]))
  }else {
    predicted_value <- reference_order |>
      dplyr::left_join(predicted_value, by = stats::setNames(c(gen_name), c(gen_name))) |>
      dplyr::select(!!dplyr::sym(gen_name), predicted.value, std.error, status)


    # predicted_value <- predicted_value[as.character(predicted_value[[gen_name]])%in%GID_tst, ]
    # #predicted_value <-predicted_value[order(as.character(predicted_value[[gen_name]])%in%GID_tst), ]
    # #### sik
    # gen_namess <- as.character(predicted_value[[gen_name]])
    # # Get the match positions of gen_names in GID_tst
    # order_indices <- match(gen_namess, GID_tst)
    # # Order the data frame based on these match positions
    # predicted_value <- predicted_value[order(order_indices), ]
    return(as.double(predicted_value[tst, "predicted.value"]))
  }


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
                          xgb_booster = "dart", #"gtree", #
                          xgb_rate_drop = 0.1,
                          xgb_skip_drop = 0.5,
                          xgb_objective = "reg:squarederror",
                          xgb_sample_type = "uniform",
                          xgb_normalize_type = "tree",
                          eta = 0.1,
                          nrounds = 100,
                          max_depth = 6,
                          scaling = FALSE,
                          centering = TRUE,
                          omic_count,
                          xgb_gamma = 0.01, ## 4
                          min_child_weight = 1,
                          subsample = 0.5,
                          colsample_bytree = 0.8,
                          xgb_alpha = 0.001, ## gblinear
                          xgb_lambda = 1.0, # gblinear,
                          early_stop_for_iteration_xgb = FALSE ## use when training set is large
                          ){

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  ### Set the paramters and hyper parameters for extreme graident boosting
  if(xgb_booster == "gblinear"){
  xgb_params <- list(
    booster = "gblinear",
    eta = eta,
    alpha = xgb_alpha,
    lambda = xgb_lambda,
    objective = "reg:squarederror",
    eval_metric = c("rmse", "rmsle", "mape")
  )
  } else {

    if(xgb_booster == "gbtree"){
      xgb_params <- list(
        booster = "gbtree",
        objective = "reg:squarederror",
        eta = eta,
        max_depth = max_depth,
        min_child_weight = min_child_weight,
        subsample = subsample,
        colsample_bytree = colsample_bytree
      )
    }

    if (xgb_booster == "dart") {
      xgb_params <- list(
        booster = xgb_booster,
        sample_type = xgb_sample_type,
        normalize_type = xgb_normalize_type,
        rate_drop = xgb_rate_drop,
        skip_drop = xgb_skip_drop,
        alpha = xgb_alpha,
        lambda = xgb_lambda,
        eta = eta,
        objective = "reg:squarederror",
        eval_metric = c("rmse", "rmsle", "mape")
      )
    }
}

  omics_Xgb <- xgboost::xgb.DMatrix(data = omics[-tst, ],
                                    label = y[-tst])

  if(isTRUE(early_stop_for_iteration_xgb)){

    yy <- y[-tst]
    omics_yy <- omics[-tst, ]
    indices <- dplyr::ntile(yy, 10)  # Creating 10 bins based on quantiles

    #set.seed(123)
    train_indices <- caret::createDataPartition(indices, p = 0.8, list = FALSE)
    training <- yy[train_indices]
    validation <- yy[-train_indices]

    train_geno <- omics_yy[train_indices, ]
    val_geno <- omics_yy[-train_indices, ]

    early_stopping_fraction <- 0.2  # 10% of iterations
    early_stopping_rounds <- max(50, floor(nrounds * early_stopping_fraction))  # At least 50 rounds

    # Creating DMatrix objects
    train_dmatrix <- xgboost::xgb.DMatrix(data = train_geno, label = training)
    eval_dmatrix <- xgboost::xgb.DMatrix(data = val_geno, label = validation)


    watchlist <- list(train = train_dmatrix, eval = eval_dmatrix)

    fit <- xgboost::xgb.train(
      params = xgb_params,
      data = omics_Xgb,
      nrounds = nrounds,
      early_stopping_rounds = early_stopping_rounds,
      watchlist = watchlist,
      maximize = FALSE, ##since these areeval_metric = c("rmse", "rmsle", "mape")
      verbose = 0
    )

  } else {

    fit <- xgboost::xgb.train(
      params = xgb_params,
      data = omics_Xgb,
      nrounds = nrounds,
      #early_stopping_rounds = early_stopping_rounds,
      #watchlist = watchlist,
      #maximize = FALSE, ##since these areeval_metric = c("rmse", "rmsle", "mape")
      verbose = 0
    )
  }


  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)



  preds <- as.data.frame(preds)

  preds <- revert_scaling(preds[, 1], y_scaler)
  return(preds)


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

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

if(is.null(ncomp) | !is.numeric(ncomp)) ncomp <- 3
if(is.null(ncomp)){
  pls_model <- pls::plsr(y ~ as.matrix(omics), validation = "CV", segments = 5)

  # Get the cross-validated RMSEP values
  rmsep_values <- pls::RMSEP(pls_model)

  # Extract the RMSEP values for cross-validation
  rmsep_cv <- rmsep_values$val["CV", , ]

  # Find the optimal number of components
  optimal_components <- which.min(rmsep_cv)

} else {
  optimal_components <- ncomp
}

  pls_model <- pls::plsr(y[-tst]~ omics[-tst, ],
                         scale = FALSE,
                         center = FALSE,
                         ncomp = optimal_components
                         #validation = "none"
                         )

  preds <- stats::predict(pls_model,
                          newdata =omics[tst, ],
                          ncomp = optimal_components)

  preds <- as.data.frame(preds)
  preds <- revert_scaling(preds[, 1], y_scaler)
  return(preds)


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

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  if(!is.null(ntree)){

  fit <- randomForest::randomForest(x = omics[-tst, ],
                                   y = y[-tst],
                                   ntree = ntree,
                                   importance = TRUE)
  } else {
    fit <- randomForest::randomForest(x = omics[-tst, ],
                                      y = y[-tst],
                                      importance = TRUE)
  }

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)
  preds <- as.data.frame(preds)
  preds <- revert_scaling(preds[, 1], y_scaler)

  return(preds)


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

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

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

  preds <- stats::predict(object = fit,
                          s= fit_CV$lambda.min,
                          newx = omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  preds <- revert_scaling(preds[, 1], y_scaler)
  return(preds)


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

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  fit_CV<-glmnet::cv.glmnet(x= omics[-tst, ],
                            y= y[-tst],
                            #nfolds = 5,
                            alpha = 1,
                            standardize = FALSE
                            )

  fit <-  glmnet::glmnet(x= omics[-tst, ],
                         y = y[-tst],
                         alpha = 1,
                         standardize = FALSE,
                         lambda =fit_CV$lambda.min)

  preds <- stats::predict(object = fit,
                          s= fit_CV$lambda.min,
                          newx = omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  preds <- revert_scaling(preds[, 1], y_scaler)

  return(preds)


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
                      omic_count = NULL,
                      k = 5){


  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  # }
  # if(!is.null(omic_count)) {
  #   if(isTRUE(scaling) || isFALSE(scaling)){
  #     omics <- scale(omics, center = TRUE, scale = TRUE)
  #   }
  # }

     fit <-  caret::knnreg(x = omics[-tst, ],
                           y = y[-tst],
                           k = k)

  preds <- stats::predict(fit,
                          omics[tst, ],
                          reshape = TRUE)

  preds <- as.data.frame(preds)

  preds <- revert_scaling(preds[, 1], y_scaler)

  return(preds)


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
                      svm_kernel = "Gaussian", # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
                      sigma_value  = 0.1,       # Default sigma value for RBF kernel
                      C_value  = 1,             # Default cost parameter
                      degree_value = 3,        # Default degree for polynomial kernel
                      scale_value  = 1,         # Default scale for polynomial kernel
                      offset_value = 0,
                      svm_type = "eps-regression") {       # Default offset for polynomial kernel
## Scale is not used because it inherently

  if(!is.null(omics)){
    scaler <- caret::preProcess(omics, method = c("center", "scale"))
    omics <- stats::predict(scaler, omics)
  }

  y_scaler <- caret::preProcess(as.data.frame(as.matrix(y)), method = c("center", "scale"))

  # Predict on the training data and get the scaled values
  y <- stats::predict(y_scaler, as.data.frame(as.matrix(y)))[, 1]

  # Translate user-friendly kernel names to `kernlab` kernel function names
  kernel_type <- switch(svm_kernel,
                        Gaussian="radial",
                        Polynomial="polynomial",
                        Linear="linear",
                        Hyperbolic_tangent="sigmoid")
  # Create a list to store kernel-specific parameters
  kernel_params <- list()

  # Set kernel parameters based on user input or defaults
  switch(kernel_type,
         radial = {kernel_params <- list(cost = C_value)},
         polynomial = {kernel_params <- list(cost = C_value, degree = degree_value)},
         Linear = {kernel_params <- list(cost = C_value)},  # Linear kernel
         sigmoid = {kernel_params <- list(cost = C_value, coef0 = offset_value)}  # Sigmoid kernel
  )

  para_index <- which(!sapply(kernel_params, is.null))
  if(length(para_index)!=0){
    kernel_params <-  kernel_params[para_index]

  } else {
    kernel_params <- list()
    kernel_params[c("degree", "coef0", "cost")] <- c(3, 0, 1)
  }

  para_names <- names(kernel_params)

  # Ensure default cost is set if not already specified
  if (!"cost" %in% para_names) {
    kernel_params$cost <- 1
  }

  if ("linear" %in% kernel_type) {
    svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], kernel = kernel_type, cost = kernel_params$cost, scale = FALSE)
  } else if ("radial" %in% kernel_type) {
    if (!"gamma" %in% para_names) {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], cost = kernel_params$cost, scale = FALSE)
    } else {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], cost = kernel_params$cost,
                              gamma = kernel_params$gamma, scale = FALSE)
    }
  } else if ("polynomial" %in% kernel_type) {
    if (all(c("cost", "degree") %in% para_names)) {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], cost = kernel_params$cost, degree = kernel_params$degree, scale = FALSE)
    } else {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], scale = FALSE)
    }
  } else if ("sigmoid" %in% kernel_type) {
    if ("coef0" %in% para_names) {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], cost = kernel_params$cost, coef0 = kernel_params$coef0, scale = FALSE)
    } else {
      svm_model <- e1071::svm(x = omics[-tst, ], y = y[-tst], scale = FALSE)
    }
  }

  preds <- stats::predict(svm_model, omics[tst, ])
  preds <- revert_scaling(preds, y_scaler)
  return(preds)
  # if(!is.null(omics)) {
  #   if(isTRUE(scaling) || isFALSE(scaling)){
  #     omics <- scale(omics, center = TRUE, scale = TRUE)
  #   }
  # }
  # if(kernel_type!="vanilladot"){
  # fit <-  kernlab::ksvm(x = omics[-tst, ],
  #                       y = y[-tst],
  #                       kernel = kernel_type,
  #                       scaled = FALSE,
  #                       type = "nu-svr",
  #                       C = C_value,
  #                       kpar = kernel_params)
  #
  # } else {
  #   if(kernel_type=="vanilladot"){
  #     AI_fit = kernlab::ksvm(x = omics[-tst, ],
  #                            y = y[-tst],
  #                            kernel = kernel_type,
  #                            scaled = FALSE,
  #                            type = "nu-svr",
  #                            C = C_value
  #     )
  #   }
  #
  # }
  #
  # preds <- kernlab::predict(fit,
  #                           omics[tst, ])
  #
  # preds <- as.data.frame(preds)
  #
  # return(preds[, 1])


}
