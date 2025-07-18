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
# set_parallel_plan <- function(n_trait,
#                               n_model = 1,
#                               replication = 1,
#                               num_cores = NULL,
#                               globals_max_GB   = 4,
#                               docker_override  = FALSE,
#                               sys_name) {
#
#   #on.exit(future::plan("sequential"), add = TRUE)
#   sys_name <- if (docker_override) "Windows" else Sys.info()[["sysname"]]
#   plan_type <- ifelse(sys_name == "Windows", "multisession", "multicore")
#
#   # Check if parallel execution is beneficial
#    if (n_trait > 1 || n_model > 1 || replication > 1) {
#   #   if(is.null(num_cores)){
#   #
#   #     num_cores <-  parallel::detectCores(logical = TRUE)
#   #     num_cores <- num_cores*0.5
#   #   }
#   #   # Ensure the number of workers does not exceed a reasonable limit
#   #   max_cores <- parallel::detectCores(logical = FALSE)  # Physical cores only
#   #   num_cores <- min(num_cores, max_cores)
#
#   phys  <- parallel::detectCores(logical = FALSE)
#   avail <- if (is.null(num_cores)) floor(phys * 0.5) else floor(num_cores)
#
#   ## keep it sane
#   workers <- max(1L, min(avail, phys, n_trait * n_model * replication))
#
#     options(future.globals.maxSize = globals_max_GB * 1024^3L,
#             future.rng.onMisuse    = "ignore")
#
#     future::plan(plan_type, workers = workers, gc = TRUE)
#
#     #message(paste0("Using ", num_cores, " workers for parallel execution"))
#   } else {
#     future::plan("sequential")
#     workers <- 1L
#   }
#
#   invisible(workers)
# }

log_thread_env_vars <- function() {
  # now log what’s actually set
  vars <- c(
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "TF_INTRA_OP_PARALLELISM_THREADS",
    "TF_INTER_OP_PARALLELISM_THREADS",
    "TF_ENABLE_ONEDNN_OPTS"
  )
  for (v in vars) {
    message(sprintf("→ %s = %s", v, Sys.getenv(v, unset = "<unset>")))
  }

  invisible(NULL)
}


set_parallel_plan <- function(n_trait,
                              n_model = 1,
                              replication = 1,
                              num_cores = NULL,
                              globals_max_GB = 4,
                              docker_override = FALSE,
                              mode = "cross_validation",
                              sys_name) {
  task_type <- if (mode == "cross_validation") "Cross Validation" else "True Prediction"
  message(sprintf(
    "\nSetting up parallel plan for %s\n",
    task_type
  ))
  sys_name  <- if (docker_override) "Windows" else Sys.info()[["sysname"]]
  plan_type <- if (sys_name == "Windows") "multisession" else "multicore"
  message(sprintf("→ OS: %s (docker_override=%s) → using future plan '%s'",
                  sys_name, docker_override, plan_type))

  phys <- parallel::detectCores(logical = FALSE)
  avail <- if (is.null(num_cores)) floor(phys * 0.5) else floor(num_cores)
  message(sprintf("→ Core budget: avail = %d (num_cores=%s)", avail,
                  if (is.null(num_cores)) "auto" else num_cores))

  max_needed <- n_trait * n_model * replication
  workers <- max(1L, min(avail, phys, n_trait * n_model * replication))
  message(sprintf("→ Workload: traits × models × reps = %d; spawning %d worker(s)",
                  max_needed, workers))

  options(future.globals.maxSize = globals_max_GB * 1024^3L,
          future.rng.onMisuse    = "ignore")
  message(sprintf("→ future.globals.maxSize = %.1f GB", globals_max_GB))

  if (workers > 1L) {
    future::plan(plan_type, workers = workers, gc = TRUE)
    message(sprintf("→ future plan set to '%s' with %d workers", plan_type, workers))
  } else {
    future::plan("sequential")
    message("→ Using sequential plan (workers = 1)")
  }

  # now log all your thread‑control env vars in one shot
  log_thread_env_vars()

  invisible(workers)
}

#' Compute per-future thread budget and apply it globally
#'
#' @return integer number of threads to use per model
set_per_worker_threads <- function() {
  # how many futures will be running in parallel?
  # workers <- future::plan()$workers
  # if (is.null(workers) || workers < 1L) workers <- 1L
  #
  # cores   <- parallel::detectCores(logical = TRUE)
  workers <- future::nbrOfWorkers()
  # how many hardware threads on the box?
  cores   <- parallel::detectCores(logical = TRUE)
  # give each future an equal share
  intra   <- max(1L, floor(cores / workers))

  # enforce in OpenMP / MKL / BLAS
  ## These also good whene using Eigen, TBB, etc.
  Sys.setenv(OMP_NUM_THREADS = intra)
  Sys.setenv(MKL_NUM_THREADS = intra)


  invisible(intra)

}

# memo_readRDS <- local({
#   cache <- new.env(parent = emptyenv())
#   function(path) {
#     if (!exists(path, envir = cache, inherits = FALSE))
#       assign(path, readRDS(path), envir = cache)
#     get(path, envir = cache, inherits = FALSE)
#   }
# })


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

  if (model == "deep_learning_model") {
    ####compute safe thread budget for DPL
    # workers <- future::plan()$workers
    # if (is.null(workers) || workers < 1L) workers <- 1L
    #
    # cores   <- parallel::detectCores(logical = TRUE)
    workers <- future::nbrOfWorkers()
    cores   <- parallel::detectCores(logical =TRUE)
    intra   <-  as.integer(max(1L, floor(cores / workers)))
    inter   <- 1L

    Sys.setenv(OMP_NUM_THREADS = intra)
    Sys.setenv(MKL_NUM_THREADS = intra)

    ##### Tell what  TF will obey for resource usage
    if (reticulate::py_module_available("tensorflow")) {
      tf <- reticulate::import("tensorflow", delay_load = TRUE)
      tf$config$threading$set_intra_op_parallelism_threads(intra)
      tf$config$threading$set_inter_op_parallelism_threads(inter)
    }

  }

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
         "deep_learning_model" = deep_learning_model(y = y, omics = omics_data, tst = tst,
                                                     scaling = additional_params$scaling,
                                                     centering = additional_params$centering,
                                                     num_hidden_layers = additional_params$num_hidden_layers,
                                                     neurons_per_layer = additional_params$neurons_per_layer,
                                                     learning_rate_dp = additional_params$learning_rate_dp,
                                                     epochs = additional_params$epochs,
                                                     batch_size = additional_params$batch_size,
                                                     l2_regularizer_dp = additional_params$l2_regularizer_dp,
                                                     dropout_rate = additional_params$dropout_rate,
                                                     crossval = additional_params$crossval,
                                                     omic_count = additional_params$omic_count,
                                                     early_stop = additional_params$early_stop,
                                                     deep_learning_model = additional_params$deep_learning_model,
                                                     n_blocks = additional_params$n_blocks,
                                                     dense_layers_cnn = additional_params$dense_layers_cnn,
                                                     kernel_size = additional_params$kernel_size,
                                                     n_neurons_per_block = additional_params$n_neurons_per_block,
                                                     attention_on_final_layer = additional_params$attention_on_final_layer,
                                                     attention_across_multiple_layers = additional_params$attention_across_multiple_layers,
                                                     batch_normalization = additional_params$batch_normalization
                                                     ),
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
         "GBLUP" = asreml_mod_cv(asreml_models_prep_cv = additional_params$asreml_models_prep_cv,
                                 pheno_data = additional_params$pheno_data,
                                 response = additional_params$response,
                                 heter_groups = additional_params$heter_groups,
                                 gen_name = additional_params$gen_name, tst = tst),
         "Bayes" = bayes_mod_cv(y = y, ETA = additional_params$ETA, weights = additional_params$weights,
                                bayes_para = additional_params$bayes_para, tst = tst,
                                bayes_model = additional_params$bayes_model, bayes_trait = additional_params$bayes_trait)
         # "ND_modes" = ND_modes_cv(y = y, omics = omics_data,
         #                          gam_method = additional_params$gam_method,
         #                          selected = additional_params$selected,
         #                          max_features = additional_params$max_features,
         #                          #k_value = additional_params$k_value,
         #                          var_explained =  additional_params$var_explained,
         #                          tst = tst)

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
                                    selected_raw, # nd_mods
                                    gam_method = NULL, # nd_mods
                                    max_features = 50,
                                    k_value = 5,
                                    var_explained = 0.9,
                                    engine = NULL,
                                    model_prep_all_bayes_cv = NULL,
                                    asreml_models_prep_cv = NULL,
                                    ml_dat_res = NULL,
                                    heter_groups = NULL,
                                    verbose = FALSE,
                                    num_cores = NULL,
                                    nfolds = 5,
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
                                    num_hidden_layers = 1,
                                    neurons_per_layer = 64,
                                    learning_rate_dp = 0.001,
                                    epochs = 10,
                                    batch_size = 32 ,
                                    l2_regularizer_dp = 0.001,
                                    dropout_rate = 0.5,
                                    validation_split = 0.2,
                                    early_stop = TRUE,
                                    deep_learning_model = "mlp_with_attention",
                                    n_blocks = 2,
                                    dense_layers_cnn = c(128, 64),
                                    kernel_size = 3,
                                    n_neurons_per_block = NULL,
                                    attention_on_final_layer = TRUE,
                                    attention_across_multiple_layers = FALSE,
                                    batch_normalization = TRUE,
                                    crossval = TRUE,
                                    docker_nd_usage = FALSE,
                                    globals_max_GB = 4,
                                    ...){

 #browser()
  on.exit(future::plan("sequential"), add = TRUE)
  msg <- "\n==================================================\n"

  if(!is.null(cross_validation_meth) & length(cross_validation_meth)>1){
    stop(paste(msg, 'use only one cross_validation method at a time.'), call. = FALSE)
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

    # if("Repeated"%in%present_patterns && is.null(replication)) {
    #   stop(paste(msg, "You select repeated cross-validation provide number of replications.\n For example, replication = 2."), call. = FALSE)
    # }


  }

  # Prepare additional parameters for model prediction
  ETA <-  NULL
  bayes_model <- NULL
  bayes_trait <- NULL
  bayes_para <- NULL

  if(!is.null(ml_dat_res)){
  omics_data <-  ml_dat_res[["merged_data"]][["merge_data"]]

  if(!"merged_data_test"%in%names(ml_dat_res)){
  if(!is.null(test_set)) {

    omics_data <-  omics_data[!rownames(omics_data) %in% test_set, ]

  }

  }


  if("omic_count"%in%names(ml_dat_res)){
    omic_count <- ml_dat_res[["omic_count"]]

  }else{
    omic_count <- NULL
  }

  #print(omics_data)
  }

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
                            asreml_models_prep_cv = asreml_models_prep_cv, ## asreml
                            gen_name = gen_name,
                            pheno_data = pheno_data,
                            #response = trait,
                            #selected_raw = selected_raw,
                            #gam_method = gam_method,
                            max_features = max_features,
                            heter_groups = heter_groups,
                            num_hidden_layers = num_hidden_layers,
                            neurons_per_layer = neurons_per_layer,
                            learning_rate_dp = learning_rate_dp,
                            epochs = epochs,
                            batch_size = batch_size,
                            l2_regularizer_dp = l2_regularizer_dp,
                            dropout_rate = dropout_rate,
                            crossval = crossval,
                            early_stop = early_stop,
                            deep_learning_model = deep_learning_model,
                            n_blocks = n_blocks,
                            dense_layers_cnn = dense_layers_cnn,
                            kernel_size = kernel_size,
                            n_neurons_per_block =n_neurons_per_block,
                            attention_on_final_layer = attention_on_final_layer,
                            attention_across_multiple_layers = attention_across_multiple_layers,
                            batch_normalization = batch_normalization
                            )


  AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare",
                       "SupportVectorMachine", "K-NearestNeighbors", "Lasso",
                       "Ridge_Regression", "deep_learning_model")

  # ND_method_use = c("ND_mod1", "ND_mod2",
  #                    "ND_mod3", "ND_mod4",
  #                    "ND_mod5", "ND_mod6",
  #                    "ND_mod7", "ND_mod8")

  #ND_method_use = c("ND_mod1", "ND_mod2")

  bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
  bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")



  dp_models <- c("mlp_with_attention", "mlp", "ResNet", "cnn")

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
      stop(paste(msg, "Provided method is not available in hold-out methods."), call. = FALSE)
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
  if(docker_nd_usage) sys_name <- "Windows"
  workers <- set_parallel_plan(n_trait = n_trait, n_model = n_model,
                    replication = replication,num_cores = num_cores,
                    sys_name = sys_name,
                    globals_max_GB = globals_max_GB,
                    docker_override = docker_nd_usage)
  # if (!is.null(num_cores) && num_cores > 1) {
  #   #sys_name <- Sys.info()["sysname"]
  #   set_parallel_plan(n_trait = n_trait, n_model = n_model,
  #                     replication = replication,num_cores = num_cores,
  #                     sys_name = sys_name,
  #                     globals_max_GB = globals_max_GB,
  #                     docker_override = docker_nd_usage)
  # } else {
  #   # Automatically determine the number of cores and use half of them
  #   # detected_cores <- parallel::detectCores(logical = TRUE)
  #   # # For non-Windows systems, consider physical cores only
  #   # num_cores <- round(detected_cores * 0.5)
  #
  #   set_parallel_plan(n_trait = n_trait,
  #                     n_model= n_model,
  #                     replication = replication,
  #                     num_cores = num_cores,
  #                     sys_name = sys_name)
  # }

  # big_paths <- list(bayes = "model_prep_all_bayes_cv.rds",
  #                   addp  = "additional_params.rds",
  #                   phe = "pheno_data.rds",
  #                   omic_dat = "omics_data.rds")

  # 2. Register your cleanup BEFORE any return()
  # on.exit({
  #   existing <- vapply(big_paths, file.exists, logical(1))
  #   if (any(existing)) {
  #     file.remove(big_paths[existing])
  #   }
  # }, add = TRUE)

  # saveRDS(model_prep_all_bayes_cv, big_paths$bayes)
  # saveRDS(additional_params,       big_paths$addp)
  # saveRDS(pheno_data,       big_paths$phe)
  # saveRDS(omics_data,       big_paths$omic_dat)

  # 3. Also reset the plan on exit
  #on.exit(future::plan("sequential"), add = TRUE)


  # if (inherits(future::plan(), "multisession")) {
  #   cl <- future::plan()$workers   # the actual cluster
  #   parallel::clusterExport(cl, varlist = "memo_readRDS", envir = environment())
  # }


  # Create a list of all combinations of response variables and replications
  tasks <- expand.grid(response = response,
                       replication = seq_len(replication),
                       modell = GS_model_cv,
                       stringsAsFactors = FALSE)


  chunk_size <- if (workers > 1) ceiling(nrow(tasks) / workers) else NULL

  results <- future.apply::future_lapply(seq_len(nrow(tasks)),
                                         future.packages   = c("dplyr"),
                                         future.seed       = TRUE,
                                         future.chunk.size = chunk_size,
                                         function(i) {
    task_row <- tasks[i, ]

    trait <- as.character(task_row$response)
    rep <- as.integer(task_row$replication)
    model <- as.character(task_row$modell)

    # if(model%in%c(bayes_valid_models, bayes_gblup_valid_models)){
    # model_prep_all_bayes_cv <- memo_readRDS(big_paths$bayes)
    # }
    # additional_params <- memo_readRDS(big_paths$addp)
    #
    # pheno_data <- memo_readRDS(big_paths$phe)
    # omics_data <- memo_readRDS(big_paths$phe)


    if(any(model%in%dp_models)){
     additional_params$deep_learning_model <- model
     #model_use <- model
     model <- "deep_learning_model"
    }

    repp <- 1

    if (!is.null(random_state) && is.numeric(random_state) && length(random_state) == 1) {
      base_seed <- as.integer(random_state)
      new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
    } else {
      base_seed <- 123L
      new_seed <- (base_seed + rep * 10000L) %% .Machine$integer.max
    }

    if (cross_validation_meth %in% c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail)) {
      if (cross_validation_meth %in% holds_out_methods_avail) {
        test_set_val <- hold_out_stratified_and_un(pheno_data = pheno_data,
                                                   gen_name = gen_name,
                                                   response = trait,
                                                   test_size = test_size,
                                                   random_state = new_seed,
                                                   replication = repp,
                                                   sampling_method = sampling_method)
      } else if (cross_validation_meth %in% Kfolds_methods_avail) {
        test_set_val <- kfolds_stratified_un(pheno_data = pheno_data,
                                             gen_name = gen_name,
                                             response = trait,
                                             test_size = test_size,
                                             nfolds = nfolds,
                                             random_state = new_seed,
                                             replication = repp,
                                             sampling_method = sampling_method)
      } else if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
        CV <- as.integer(strsplit(cross_validation_meth, "CV")[[1]][2])
        test_set_val <- CV1_CV2_for_multi_environment(pheno_data = pheno_data,
                                                      gen_name = gen_name,
                                                      response = trait,
                                                      test_size = test_size,
                                                      CV = CV,
                                                      nfolds = nfolds,
                                                      heter_groups = heter_groups,
                                                      random_state = new_seed,
                                                      replication = repp,
                                                      sampling_method = sampling_method)
      }
    } else {
      stop(paste(msg, "Unsupported cross-validation method specified. Choose from: ",
                 paste(c(holds_out_methods_avail, Kfolds_methods_avail, CVs_multi_envs_methods_avail), collapse = ", ")), call. = FALSE)
    }

    len_y <- nrow(pheno_data)
    y <- as.double(pheno_data[[trait]])

    if (cross_validation_meth %in% c(Kfolds_methods_avail, holds_out_methods_avail)) {
      ypred_cv <- matrix(data=NA, nrow=len_y, ncol=2)
      colnames(ypred_cv) <- c("y", "yhat")
      ypred_cv <- as.data.frame(ypred_cv)
      ypred_cv[, "y"] <- as.double(pheno_data[[trait]])

      results_eval_metrics_reps <- matrix(NA, nrow = repp, ncol = length(eval_metrics)+1)
      results_eval_metrics_reps[, 1] <- 1
      rownames(results_eval_metrics_reps) <- paste("REP", 1, sep = "_")
      colnames(results_eval_metrics_reps) <- c("Rep", eval_metrics)
      results_eval_metrics_reps <- as.data.frame(results_eval_metrics_reps)
    }

    if (cross_validation_meth %in% CVs_multi_envs_methods_avail) {
      if (is.null(heter_groups)) {
        stop(message(paste(msg, 'For CV1 or CV2 column name for environment/location is required.')), call. = FALSE)
      }
      ENV <- as.character(unique(pheno_data[[heter_groups]]))
      ypred_cv <- matrix(data=NA, nrow=len_y, ncol=3)
      colnames(ypred_cv) <- c("y", "yhat", heter_groups)
      ypred_cv <- as.data.frame(ypred_cv)
      ypred_cv[["y"]] <- as.double(pheno_data[[trait]])
      ypred_cv[[heter_groups]] <- as.character(pheno_data[[heter_groups]])

      results_eval_metrics_reps <- matrix(NA, nrow = length(ENV), ncol = length(eval_metrics)+2)
      results_eval_metrics_reps[, 1] <- rep(1, length(ENV))
      results_eval_metrics_reps[, 2] <- ENV
      colnames(results_eval_metrics_reps) <- c("Rep", heter_groups, eval_metrics)
      results_eval_metrics_reps_use <- results_eval_metrics_reps
      results_eval_metrics_reps <- data.frame()
    }

    group <- test_set_val[[repp]]

    handle_error <- FALSE

    if (!cross_validation_meth %in% holds_out_methods_avail) {
      for (j in 1:nfolds) {
        yNA <- y
        for (g in 1:len_y) {
          if (group[g] == j) { yNA[g] <- NA }
        }
        tst <- which(is.na(yNA))

        if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            if (model == "GBLUP_BRR") {
              model_GBLUP <- "BRR"
              additional_params$bayes_model <-  model_GBLUP #model
              additional_params$bayes_trait <- trait
              additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
              #rm(model_GBLUP)
            } else {
              additional_params$bayes_model <- model
              additional_params$bayes_trait <- trait
              additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
            }
            ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
                                                        tst = tst, additional_params = additional_params)
          }, error = function(e) {
            message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }

        if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
          tryCatch({
            additional_params$response <- trait
            preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
            if (!is.null(preds)) {
              ypred_cv[tst, "yhat"] <- preds
            } else {
              stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
            }
          }, error = function(e) {
            message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }

        if (model %in% c(AI_valid_models)) {
          tryCatch({
            ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = omics_data,
                                                        tst = tst, additional_params = additional_params)

          }, error = function(e) {
            message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }
        ##
      }
    } else {
      if (cross_validation_meth %in% holds_out_methods_avail) {
        tst <- test_set_val[[repp]]
        yNA <- y
        yNA[tst] <- NA

        if (model %in% c(bayes_valid_models, bayes_gblup_valid_models)) {
          tryCatch({
            if (model == "GBLUP_BRR") {
              model_GBLUP <- "BRR"
              additional_params$bayes_model <- model_GBLUP#model
              additional_params$bayes_trait <- trait
              additional_params$ETA <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para <- model_prep_all_bayes_cv[[model_GBLUP]][["bayes_para"]]
              #rm(model_GBLUP)
            } else {
              additional_params$bayes_model <- model
              additional_params$bayes_trait <- trait
              additional_params$ETA <- model_prep_all_bayes_cv[[model]][["bayes_ETA"]][["ETA"]]
              additional_params$bayes_para <- model_prep_all_bayes_cv[[model]][["bayes_para"]]
            }
            ypred_cv[tst, "yhat"] <- predict_with_model(model = "Bayes", y = yNA,
                                                        tst = tst, additional_params = additional_params)
          }, error = function(e) {
            message(paste("Error in processing Bayes model", model, "for", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }

        if (!is.null(engine) && model == "GBLUP" && engine == "asreml") {
          tryCatch({
            additional_params$response <- trait
            preds <- predict_with_model(model = model, tst = tst, additional_params = additional_params)
            if (!is.null(preds)) {
              ypred_cv[tst, "yhat"] <- preds
            } else {
              stop(paste(msg, "Prediction with GBLUP model failed.\n"), call. = FALSE)
            }
          }, error = function(e) {
            message(paste("Error in processing GBLUP model for ", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }

        if (model %in% c(AI_valid_models)) {
          tryCatch({
            ypred_cv[tst, "yhat"] <- predict_with_model(model = model, y = yNA, omics_data = omics_data,
                                                        tst = tst, additional_params = additional_params)

            }, error = function(e) {
            message(paste("Error in processing AI model", model, "for", trait, ": ", e$message))
            handle_error <<- TRUE
          })
        }
        ##
      }
    }

    if (isTRUE(handle_error)) {
      #return(NULL)
      ypred_cv <- NULL
      results_eval_metrics_reps <- NULL
    }

    if (cross_validation_meth %in% CVs_multi_envs_methods_avail & isFALSE(handle_error)) {
      ypred_cv <- as.data.frame(ypred_cv)
      ypred_cv[, 'y'] <- as.double(ypred_cv[, 'y'])
      ypred_cv[, 'yhat'] <- as.double(ypred_cv[, 'yhat'])

      results_eval_metrics_reps_use <- as.data.frame(results_eval_metrics_reps_use)

      for (eva in 1:length(eval_metrics)) {
        sik <- ypred_cv |>
          dplyr::group_by(!!dplyr::sym(heter_groups)) |>
          dplyr::summarise(
            eval_metric = tryCatch(
              evaluation_metrics(yhat, y, eval_metrics = eval_metrics[eva]),
              error = function(e) NULL
            )
          )

        if (!is.null(sik)) {
          for (ii in 1:nrow(sik)) {
            if (!is.null(sik$eval_metric[ii])) {
              results_eval_metrics_reps_use[results_eval_metrics_reps_use[, heter_groups] == as.character(sik[[heter_groups]])[ii], c("Rep", eval_metrics[eva])] <- c(repp, sik$eval_metric[ii])
            }
          }
        }
      }

      results_eval_metrics_reps <- rbind(results_eval_metrics_reps_use, results_eval_metrics_reps)
    } else {
      if (!cross_validation_meth %in% CVs_multi_envs_methods_avail) {
        for (eva in 1:length(eval_metrics)) {
          tryCatch({
            results_eval_metrics_reps[repp, eval_metrics[eva]] <- evaluation_metrics(y_observed = ypred_cv[tst, "y"],
                                                                                     y_predicted = ypred_cv[tst, "yhat"],
                                                                                     eval_metrics = eval_metrics[eva])
          }, error = function(e) {
            results_eval_metrics_reps[repp, eval_metrics[eva]] <- NA
          })
        }
      }
    }

    if(model == "deep_learning_model") model <- as.character(task_row$modell)

    list(trait = trait,
         rep = rep,
         model = model,
         eval_metrics_reps = results_eval_metrics_reps,
         ypred_cv_Reps_all = ypred_cv)
  })

  # , future.seed = TRUE)

  # Process the results
  filtered_results <- lapply(results, function(res) {
    if (!is.null(res$eval_metrics_reps) && !is.null(res$ypred_cv_Reps_all)) {
      # Check if any element in eval_metrics_reps or ypred_cv_Reps_all is NA
      if (all(!is.na(unlist(res$eval_metrics_reps))) && all(!is.na(unlist(res$ypred_cv_Reps_all)))) {
        return(res)
      }
    }
    return(NULL)
  })
  # Remove NULL elements
  filtered_results <- Filter(Negate(is.null), filtered_results)

  future::plan("sequential")  # Reset to default plan
  return(filtered_results) # Adjust depending on how you want to return or further process the results


}
####


