
#' Title
#'
#' @param n_trait
#' @param n_model
#' @param replication
#' @param num_cores
#' @param sys_name
#'
#' @return
#' @export
#'
#' @examples
set_parallel_plan <- function(n_trait,
                              n_model = 1,
                              replication = 1,
                              num_cores = NULL,
                              sys_name) {
  # Define the plan based on the system
  plan_type <- ifelse(sys_name == "Windows", "multisession", "multicore")

  # Check if parallel execution is beneficial
  if (n_trait > 1 || n_model > 1 || replication > 1) {
    if(is.null(num_cores)){

      num_cores <-  parallel::detectCores()
      num_cores <- num_cores*0.7
    }
    future::plan(plan_type, workers = num_cores)
  } else {
    future::plan("sequential")
  }
}


#' Title
#'
#' @param model
#' @param y
#' @param omics_data
#' @param tst
#' @param additional_params
#'
#' @return
#' @export
#'
#' @examples
predict_with_model <- function(model = NULL,
                               y = NULL,
                               omics_data = NULL,
                               tst = NULL,
                               additional_params = NULL) {
  # Generalized function to handle predictions for various models
  # 'additional_params' is a list of additional parameters required for each model

  switch(model,
         "Xgboost" = AI_xgboost_cv(y = y, omics = omics_data, tst = tst, eta = additional_params$eta,
                                   nrounds = additional_params$nrounds, max_depth = additional_params$max_depth,
                                   scaling = additional_params$scaling, omic_count = additional_params$omic_count,
                                   centering = additional_params$centering, xgb_gamma = additional_params$xgb_gamma,
                                   colsample_bytree = additional_params$colsample_bytree,
                                   subsample = additional_params$subsample, min_child_weight = additional_params$min_child_weight,
                                   xgb_alpha = additional_params$xgb_alpha, xgb_lambda = additional_params$xgb_lambda,
                                   xgb_booster = additional_params$xgb_booster,
                                   xgb_rate_drop = additional_params$xgb_rate_drop,xgb_skip_drop = additional_params$xgb_skip_drop,
                                   xgb_objective = additional_params$xgb_objective,xgb_sample_type = additional_params$xgb_sample_type,
                                   xgb_normalize_type = additional_params$xgb_normalize_type,
                                   early_stop_for_iteration_xgb = additional_params$early_stop_for_iteration_xgb),
         "RandomForest" = AI_randomforest_cv(y = y, omics = omics_data, tst = tst,
                                             scaling = additional_params$scaling,
                                             centering = additional_params$centering, ntree = additional_params$ntree,
                                             omic_count = additional_params$omic_count),
         "PartialLeastSquare" = AI_pls_cv(y = y, omics = omics_data, tst = tst,
                                          scaling = additional_params$scaling,
                                          centering = additional_params$centering, ncomp = additional_params$ncomp,
                                          omic_count = additional_params$omic_count),
         "Ridge_Regression" = AI_ridge_regression_cv(y = y, omics = omics_data, tst = tst,
                                                     scaling = additional_params$scaling,
                                                     centering = additional_params$centering,
                                                     omic_count = additional_params$omic_count),
         "Lasso" = AI_lasso_cv(y = y, omics = omics_data, tst = tst, scaling = additional_params$scaling,
                               centering = additional_params$centering,
                               omic_count = additional_params$omic_count),
         "SupportVectorMachine" = AI_svm_cv(y = y, omics = omics_data, tst = tst,
                                            scaling = additional_params$scaling,
                                            centering = additional_params$centering,
                                            C_value = additional_params$C_value,
                                            degree_value = additional_params$degree_value,
                                            scale_value = additional_params$scale_value,
                                            offset_value = additional_params$offset_value,
                                            omic_count = additional_params$omic_count),
         "K-NearestNeighbors" = AI_knn_cv(y = y, omics = omics_data, tst = tst,
                                          scaling = additional_params$scaling,
                                          centering = additional_params$centering, k = additional_params$k,
                                          omic_count = additional_params$omic_count),
         "GBLUP" = asreml_mod_cv(asreml_models_prep_cv = additional_params$asreml_models_prep_cv, pheno_data = additional_params$pheno_data,
                                    response = additional_params$response, heter_groups = additional_params$heter_groups,
                                 gen_name = additional_params$gen_name, tst = tst),
         "Bayes" = bayes_mod_cv(y = y, ETA = additional_params$ETA, weights = additional_params$weights,
                                bayes_para = additional_params$bayes_para, tst = tst,
                                bayes_model = additional_params$bayes_model, bayes_trait = additional_params$bayes_trait)
  )
}

###########################

#' Title
#'
#' @param pheno_data
#' @param test_set
#' @param response
#' @param gen_name
#' @param test_size
#' @param random_state
#' @param replication
#' @param weights
#' @param model_prep_all_bayes_cv
#' @param asreml_models_prep_cv
#' @param ml_dat_res
#' @param heter_groups
#' @param verbose
#' @param num_cores
#' @param nfolds
#' @param cross_validation_meth
#' @param sampling_method
#' @param eval_metrics
#' @param GS_model_cv
#' @param scaling
#' @param centering
#' @param eta
#' @param nrounds
#' @param max_depth
#' @param gamma
#' @param subsample
#' @param colsample_bytree
#' @param ncomp
#' @param ntree
#' @param k
#' @param c
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
models_execute_crossval <- function(pheno_data = NULL,
                                    test_set = NULL,
                                    response = NULL,
                                    gen_name = NULL,
                                    test_size = NULL,
                                    random_state = NULL,
                                    replication = NULL,
                                    weights = NULL,
                                    engine = NULL,
                                    model_prep_all_bayes_cv = NULL,
                                    asreml_models_prep_cv = NULL,
                                    ml_dat_res = NULL,
                                    heter_groups = NULL,
                                    verbose = FALSE,
                                    num_cores = NULL,
                                    nfolds = NULL,
                                    cross_validation_meth = NULL, ## this handle th ETA for bayes model
                                    sampling_method = NULL,
                                    eval_metrics = NULL,
                                    bayes_model = NULL,
                                    GS_model_cv = NULL,
                                    scaling = FALSE,
                                    centering = TRUE,
                                    eta = 0.1, ## xgboost
                                    nrounds = 100, ## xgboost
                                    max_depth = 6, ## xgboost
                                    xgb_gamma = 4, ## xgboost
                                    subsample = 0.5, ## xgboost
                                    colsample_bytree = 1, ## xgboost
                                    xgb_alpha = 0.001, ## xgboost linear
                                    xgb_lambda = 1, ## xgboost linear
                                    min_child_weight = 1, ## xgboost
                                    early_stop_for_iteration_xgb = TRUE,
                                    xgb_booster = "dart", #"gbtree",
                                    xgb_rate_drop = 0.1,
                                    xgb_skip_drop = 0.5,
                                    xgb_objective = "reg:squarederror",
                                    xgb_sample_type = "uniform",
                                    xgb_normalize_type = "tree",
                                    ncomp = 3, #### pls
                                    ntree = 500, ### random forest
                                    k = 5, ## for knn
                                    svm_kernel = "Gaussian", # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
                                    sigma_value  = 0.1,       # Default sigma value for RBF kernel
                                    C_value  = 1,             # Default cost parameter
                                    degree_value = 3,        # Default degree for polynomial kernel
                                    scale_value  = 1,         # Default scale for polynomial kernel
                                    offset_value = 1,
                                    ...){

 #browser()

  if(!is.null(cross_validation_meth) & length(cross_validation_meth)>1){
    stop('use only one cross_validation method at a time')
  }
  ##### if user choose repeated and stratified CV strategy but
  ## forget to choose sampling stratgy or replication is not defined.
  patterns <- c("stratified", "Repeated")

  # Use sapply to apply grep to each pattern and return a named logical vector indicating presence.
  if(is.null(sampling_method) | is.null(replication)){
    matche_strings <- sapply(patterns, function(pattern) {
      length(grep(pattern, cross_validation_meth, ignore.case = TRUE)) > 0
    }, simplify = FALSE)

    # Name the list elements with the patterns.
    names(matche_strings) <- patterns

    # Filter to get only the patterns that were found.
    present_patterns <- names(matche_strings)[unlist(matche_strings)]

    if("stratified"%in%present_patterns) sampling_method <- "stratified"

    if("Repeated"%in%present_patterns) {
      stop("You select repeated cross-validation provide number of replications.\n For example, replication = 2")
    }


  }

  # Prepare additional parameters for model prediction
  ETA <-  NULL
  bayes_model <- NULL
  bayes_trait <- NULL
  bayes_para <- NULL

  omic_count <- if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL

  additional_params <- list(ETA = ETA, weights = weights, bayes_para = bayes_para,
                            bayes_model = bayes_model, bayes_trait = bayes_trait,
                            scaling = scaling,centering = centering,
                            eta = eta, nrounds = nrounds,
                            max_depth = max_depth, xgb_gamma = xgb_gamma,
                            early_stop_for_iteration_xgb = early_stop_for_iteration_xgb,
                            colsample_bytree = colsample_bytree, subsample = subsample, ntree = ntree,
                            xgb_alpha = xgb_alpha, xgb_lambda = xgb_lambda,
                            min_child_weight = min_child_weight,
                            xgb_booster = xgb_booster, xgb_rate_drop = xgb_rate_drop,
                            xgb_skip_drop = xgb_skip_drop,
                            xgb_objective = xgb_objective,
                            xgb_sample_type = xgb_sample_type,
                            xgb_normalize_type = xgb_normalize_type,
                            ncomp = ncomp, C_value = C_value,degree_value = degree_value,
                            scale_value = scale_value, offset_value = offset_value,
                            k = k, omic_count = omic_count,
                            asreml_models_prep_cv = asreml_models_prep_cv, gen_name = gen_name,
                            pheno_data = pheno_data, response = response, heter_groups = heter_groups)


  AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare",
                       "SupportVectorMachine", "K-NearestNeighbors", "Lasso",
                       "Ridge_Regression")

  bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")



  asreml_model <- "GBLUP"
  #
  # GS_model_cv = c(AI_valid_models,
  #                 bayes_valid_models,
  #                 bayes_gblup_valid_models,
  #                 asreml_model)

  # Check if "GBLUP_BRR" is in the list and replace it with "BRR"
  # if("GBLUP_BRR" %in% GS_model_cv) {
  #   GS_model_cv[GS_model_cv == "GBLUP_BRR"] <- "BRR"
  # }

  n_trait <- length(response)
  n_model <- length(GS_model_cv)

  msg <- sprintf("==================================================\n")

  holds_out_methods_avail <- c("Hold_Out",
                               "Stratified_Hold_Out",
                               "Repeated_Hold_Out",
                               "Repeated_Stratified_Hold_Out"
                               #"Leave_one_Out"
                               )

  Kfolds_methods_avail <- c("K-Folds",
                            "Stratified_K-Folds",
                            "Repeated_K-Folds",
                            "Repeated_Stratified_K-Folds")


  # Create a named vector to map hold-out methods to K-folds methods
  method_mapping <- setNames(Kfolds_methods_avail, holds_out_methods_avail)

  # Function to convert hold-out method to K-folds method
  convert_to_kfolds <- function(method) {
    if (method %in% names(method_mapping)) {
      return(method_mapping[method])
    } else {
      stop("Provided method is not available in hold-out methods.")
    }
  }

  if(cross_validation_meth %in% holds_out_methods_avail){

    cross_validation_meth <- as.character(convert_to_kfolds(cross_validation_meth))

    if(is.null(nfolds)) nfolds <- 5

  }


  CVs_multi_envs_methods_avail <- c("CV1",
                                    "CV2",
                                    "Repeated_CV1",
                                    "Repeated_CV2")


  #ypred_cv_Reps_all <- list()
  ###############################################################
  # Main logic
  sys_name <- Sys.info()["sysname"]
  if (!is.null(num_cores) && num_cores > 1) {
    #sys_name <- Sys.info()["sysname"]
    set_parallel_plan(n_trait = n_trait, n_model = n_model,
                      replication = replication,num_cores = num_cores,
                      sys_name = sys_name)
  } else {
    # Automatically determine the number of cores and use half of them
    detected_cores <- parallel::detectCores(logical = TRUE)
    # For non-Windows systems, consider physical cores only
    num_cores <- round(detected_cores * 0.5)

    set_parallel_plan(n_trait = n_trait,
                      n_model= n_model,
                      replication = replication,
                      num_cores = num_cores,
                      sys_name = sys_name)
  }

  # Create a list of all combinations of response variables and replications
  tasks <- expand.grid(response = response,
                       replication = seq_len(replication),
                       modell = GS_model_cv,
                       stringsAsFactors = FALSE)

  # Execute each task in parallel
  #process_task <- function(task) {
  results <- future.apply::future_lapply(seq_len(nrow(tasks)), function(i) {
    task_row <- tasks[i, ]

    trait <- as.character(task_row$response)
    rep <- as.integer(task_row$replication)
    model <- as.character(task_row$modell)

    y_scaler <- caret::preProcess(as.data.frame(as.matrix(pheno_data[[trait]])), method = c("center", "scale"))

    # Predict on the training data and get the scaled values
    pheno_data[, trait] <- stats::predict(y_scaler, as.data.frame(as.matrix(pheno_data[[trait]])))[, 1]

    #print(c(trait, rep))  # For diagnostic purposes

    repp <- 1
    #######
    # Check if a seed was provided and calculate a new seed based on the replication number
    if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1) {
      # Ensure random_state is an integer
      base_seed <- as.integer(random_state)
      # Generate a new, valid integer seed for each replication
      # The modulo operation ensures the seed stays within the integer range
      new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
      #new_seed <- base_seed + rep
    } else {
      # Fallback seed if random_state is not set
      base_seed <- 123L
      new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
    }

    # Apply the new seed for replication this is important to prevent same seed is used across replication

    #######
    if(cross_validation_meth %in% c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail)) {
      if(cross_validation_meth %in% holds_out_methods_avail) {
        test_set_val <- hold_out_stratified_and_un(pheno_data = pheno_data,
                                                   gen_name = gen_name,
                                                   response = trait,
                                                   test_size = test_size,
                                                   random_state = new_seed,
                                                   replication = repp,  # Use rep here if your function supports per-replication processing
                                                   sampling_method = sampling_method)
      } else if(cross_validation_meth %in% Kfolds_methods_avail) {
        test_set_val <- kfolds_stratified_un(pheno_data = pheno_data,
                                             gen_name = gen_name,
                                             response = trait,
                                             test_size = test_size,
                                             nfolds = nfolds,
                                             random_state = new_seed ,
                                             replication = repp,  # Use rep here if your function supports per-replication processing
                                             sampling_method = sampling_method)
      } else if(cross_validation_meth %in% CVs_multi_envs_methods_avail) {
        CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
        test_set_val <- CV1_CV2_for_multi_environment(pheno_data = pheno_data,
                                                      gen_name = gen_name,
                                                      response = trait,
                                                      test_size = test_size,
                                                      CV = CV,
                                                      nfolds = nfolds,
                                                      heter_groups = heter_groups,
                                                      random_state = new_seed,
                                                      replication = repp,  # Use rep here if your function supports per-replication processing
                                                      sampling_method = sampling_method)
      }
    } else {
      stop("Unsupported cross-validation method specified. Choose from: ",
           paste(c(holds_out_methods_avail,
                   Kfolds_methods_avail,
                   CVs_multi_envs_methods_avail), collapse = ", "), call. = FALSE)
    }

    len_y <-  nrow(pheno_data)

    y <- pheno_data[, trait]

    if(cross_validation_meth %in%c(Kfolds_methods_avail,
                                   holds_out_methods_avail)){

      ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
      colnames(ypred_cv) <- c("y", "yhat")
      ypred_cv[, "y"] <- pheno_data[, trait]
      ##################
      ## Add rep to the column in addition to the metrics
      results_eval_metrics_reps <- matrix(NA, nrow = repp, ncol = length(eval_metrics)+1)
      results_eval_metrics_reps[, 1] <-  c(1:repp)
      rownames(results_eval_metrics_reps) <-  paste("REP",1:repp, sep = "_")
      colnames(results_eval_metrics_reps) <-  c("Rep", eval_metrics)

      results_eval_metrics_reps <-  as.data.frame(results_eval_metrics_reps)

    }

    ########################################

    if (cross_validation_meth %in% CVs_multi_envs_methods_avail){

      if(is.null(heter_groups )){
        stop(message(paste(msg,'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)

      }
      ENV <-  as.character(unique(pheno_data[, heter_groups]))

      ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
      colnames(ypred_cv) <- c("y", "yhat", heter_groups)
      ypred_cv[, "y"] <-  pheno_data[, trait]
      ypred_cv[, heter_groups] <-  as.character(pheno_data[, heter_groups])

      results_eval_metrics_reps <- matrix(NA, nrow = length(ENV), ncol = length(eval_metrics)+2) ### Add rep and Env to the cols
      results_eval_metrics_reps[, 1] <-  rep(1, length(ENV))
      results_eval_metrics_reps[, 2] <-  ENV

      #rownames(results_eval_metrics_all) = paste("REP",1:(NRep*length(ENV)), sep = "_")
      colnames(results_eval_metrics_reps) <-  c("Rep", "Env",
                                                eval_metrics)
      results_eval_metrics_reps_use <-   results_eval_metrics_reps
      ## Convert it to empyt dataframe for final storage
      results_eval_metrics_reps <-  data.frame()

    }

    ################
    #y <-  pheno_data[, response]

    #for (k in 1:NRep) {

    group <- test_set_val[[repp]]

    if(!cross_validation_meth %in%holds_out_methods_avail){

      #folds = test_set_val[[k]]

      for(j in 1:nfolds){

        yNA <- y
        ## set the g fold to NA as the testing set
        for (g in 1:len_y) {
          if(group[g] == j) { yNA[g]<- NA }
        }

        ## extract the position NA which is the testing set
        tst = which(is.na(yNA))

        if(model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          #model_use <- model
          #model <- "Bayes" # Use a general term for Bayesian models for the switch function
          if(model == "GBLUP_BRR") {
            model_GBLUP <- "BRR"
            additional_params$bayes_model <- model
            additional_params$bayes_trait <- trait
            additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
            additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]

          } else{
            additional_params$bayes_model <- model
            additional_params$bayes_trait <- trait
            additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
            additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]

          }
          # additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
          # additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]

          ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
                                                      tst = tst, additional_params = additional_params)
          #model <- model_use
        }

        if(!is.null(engine)){
        if(model == "GBLUP" && engine == "asreml") {

          ypred_cv[tst, "yhat"] <- predict_with_model(model = model, tst = tst, additional_params = additional_params)

        }

        }

        if(model %in% c(AI_valid_models)) {
          omics_data <-  ml_dat_res[["merged_data"]][["merge_data"]]

          if(!is.null(test_set)) omics_data[rownames(omics_data) %in% test_set, ] else omics_data

        ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = omics_data,
                                                    tst = tst, additional_params = additional_params)
        }

      } ### Ends folds

    } else {
      if(cross_validation_meth %in%holds_out_methods_avail){
      tst <-  test_set_val[[repp]]
      yNA <- y
      yNA[tst] <- NA

      if(model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
        #model <- "Bayes" # Use a general term for Bayesian models for the switch function
        if(model == "GBLUP_BRR") {
          model_GBLUP <- "BRR"
          additional_params$bayes_model <- model
          additional_params$bayes_trait <- trait
          additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
          additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]

        } else{
          additional_params$bayes_model <- model
          additional_params$bayes_trait <- trait
          additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
          additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]

        }
        ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
                                                    tst = tst, additional_params = additional_params)
      }

      if(!is.null(engine)){
      if(model == "GBLUP" && engine == "asreml") {

        ypred_cv[tst, "yhat"] <- predict_with_model(model = model, tst = tst, additional_params = additional_params)

      }

      }

      if(model %in% c(AI_valid_models)) {
        ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = ml_dat_res[["merged_data"]][["merge_data"]],
                                                    tst = tst, additional_params = additional_params)
      }

      }

}## End

    if (cross_validation_meth %in%CVs_multi_envs_methods_avail){
      ypred_cv = as.data.frame(ypred_cv)
      ypred_cv[, 'y'] =  as.double(ypred_cv[, 'y'])
      ypred_cv[, 'yhat'] =  as.double(ypred_cv[, 'yhat'])

      results_eval_metrics_reps_use = as.data.frame(results_eval_metrics_reps_use)
      for (eva in 1:length(eval_metrics)) {

        sik <- unlist(doBy::lapplyBy(~Env,data=ypred_cv,
                                     function(x){evaluation_metrics(x$yhat,x$y,
                                                                    eval_metrics = eval_metrics[eva])}))
        for (s in 1:length(sik)) {
          results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups]%in%names(sik)[s], c("Rep", eval_metrics[eva])] = c(repp, sik[s])

        }
      }
      results_eval_metrics_reps = rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)
    } else {
      if (!cross_validation_meth%in%CVs_multi_envs_methods_avail){
        for (eva in 1:length(eval_metrics)) {
          results_eval_metrics_reps[repp, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
                                                                                   y_predicted = ypred_cv [tst, "yhat"],
                                                                                   eval_metrics = eval_metrics[eva])
        }

      }

    }


    #} ### End Replication

    ### Check if some value are not double or numeric and convert it
    if(!all(sapply(eval_metrics, function(x, results_eval_metrics_reps) is.numeric(results_eval_metrics_reps[,x]),  results_eval_metrics_reps))) {
      results_eval_metrics_reps[, eval_metrics] <-
        lapply(results_eval_metrics_reps[, eval_metrics, drop = FALSE],
               function(x) as.double(as.character(x)))
    }

    list(trait = trait, rep = rep, model = model, eval_metrics_reps = results_eval_metrics_reps,
         ypred_cv_Reps_all = ypred_cv)
  }, future.seed = TRUE)



  future::plan("sequential")  # Reset to default plan
  return(results) # Adjust depending on how you want to return or further process the results


}
####
# # Initialize empty lists to hold the data frames
# trait_dfs <- list()
# mean_metrics_dfs <- list()
#
# # Extract unique traits from the results
# traits <- unique(sapply(TT, function(x) x$trait))
#
# # Aggregate results by trait and calculate mean metrics
# for (trait in traits) {
#   # Filter results for the current trait
#   trait_results <- Filter(function(x) x$trait == trait, TT)
#
#   # Combine eval_metrics_reps into a single data frame for the trait
#   trait_df <- do.call(rbind, lapply(trait_results, function(x) x$eval_metrics_reps))
#   trait_dfs[[trait]] <- trait_df
#
#   # Calculate mean of evaluation metrics across all replications for the trait
#   mean_metrics <- colMeans(trait_df[, -1], na.rm = TRUE)  # Assuming first column is 'Rep'
#   mean_metrics_df <- as.data.frame(t(mean_metrics))
#   colnames(mean_metrics_df) <- names(mean_metrics)
#   mean_metrics_dfs[[trait]] <- mean_metrics_df
# }
#
# # Assign names to the lists based on trait names for easier reference
# names(trait_dfs) <- traits
# names(mean_metrics_dfs) <- traits
#
# # Now, trait_dfs contains a data frame for each trait with all replications
# # mean_metrics_dfs contains a data frame for each trait with mean metrics
#
# # To access the data frame for a specific trait, e.g., "Yield", use:
# yield_df <- trait_dfs[["Yield"]]
# yield_mean_metrics <- mean_metrics_dfs[["Yield"]]
#
# # To see the structure or print the results, you can use str() or print():
# str(yield_df)
# print(yield_mean_metrics)
#
#
