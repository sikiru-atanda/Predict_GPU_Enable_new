#' Title
#'
#' @param pheno_data phenotypic data, which can be contain both training and testing set. NA is allowed. Dataframe or matrix is allowed
#' @param geno_data Genomic/SNP/Marker data NA is allowed but not expected.
#' numeric 0, 1, 2 (where 0 is minor allele, 1 is hetero and 2 is major allele)
#' and -1, 0, 1 is also allowed (where -1 is minor allele, 0 is hetero and 1 is major allele).
#'  Dataframe or matrix is allowed
#'  We allowed up to 4 different omics data for model fit
#' @param omic1_data Omic data (transcriptomic, metabolic, proteomic, environment etc) NA is allowed but not expected. Dataframe or matrix is allowed
#' @param omic2_data Similar to Omic1_data
#' @param omic3_data Similar to Omic1_data
#' @param gmatrix    Genomic relationship matrix, NA not allowed. Dataframe or matrix is allowed
#' @param train_geno_data Genomic data for training set if geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param train_omic1_data Omic data for training set
#' @param train_omic2_data Omic data for training set
#' @param train_omic3_data Omic data for training set
#' @param test_geno_data  Genomic data for testing set if geno_data is not provided or not included in the geno_data by the user.
#' In that scenario geno_data is considered training set and training_geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param test_omic1_data Omic data for testing set if omic1_data is not provided or not included in the omic1_data by the user.
#' In that scenario omic1_data is considered training set and training_geno_data is not provided by the user. Dataframe or matrix is allowed
#' @param test_omic2_data same as test_omic1_data
#' @param test_omic3_data same as test_omic1_data
#' @param train_set Dataframe with column name of the individual in the training set. This is useful maining
#' for purpose of cross-validation exercise.
#' @param omics_data_label  lable/name of the omics data
#' @param omics_kernel_label  lable/name for the omics_kernel if any
#' @param train_omics_label
#' @param test_omics_label
#' @param test_set Dataframe with column name of the individual in the testing set. Not required
#' if pheno_data contain individuals (testing set) with no phenotypic record as NA.
#' @param gmatrix_method two methods are currently available to calculate the genomic relationship matrix
#' Yang and Van-raden
#' @param response trait(s) of interest to the user
#' @param gen_name Column name containing individuals/genotypes
#' @param cova covariate if any.Its epected in formula i.e cova  = ~ Rain + Temp
#' @param fixed fixed terms. Its expected in formula i.e fixed = ~ name + Env
#' @param random random terms. Its expected in formula i.r random = ~ name + Env
#' @param heter_resid True or False if user want heterogeneous residual variance or not
#' @param heter_groups  Column name for Environment or location
#' @param weights weight for the response variable. Only dataframe
#' @param nIter  number of iteration for Bayesian models
#' @param burnIn number of burnin  for Bayesian models
#' @param thin   number of thinning for Bayesian models
#' @param GS_model GS-model for fit. The following are available
#' BayesA, BayesB, BayesC, Baysian Ridge Regression (BRR). These models only work with
#' M-matrix(genomic data and omic data) in single location.
#' Bayesian reproducing kernel Hilbert spaces regressions (RKHS),
#' Bayesian Genomic Best linear unbias estimate (BGBLUP). Both RKHS and BGBLUP can
#' fit both single and multiple location using reaction norm.
#' Genomic Best linear unbias estimate using asreml-R package with different
#' variance structure such as (FA, RR, US, CORGH, CORGV, CORH, CORV)
#' for multi-location/environment.
#' Machine learning models include:
#' Extreme Gradiant Boosting, Random Forest, KNN, Lasso, Ridge Regression,
#' Partial Least Square, Support Vector Machine. All the machine learning only work
#' in single location.
#' @param fixed_term_model_bayesian model for the fixed term which is always fixed
#' @param rand_term_model_bayesian model for the random terms which can be any of the above mentioned model
#' @param core number of ram for paralllel job
#' @param message if message/warning should be displayed
#' @param gkernel relationship matrix using different kernel methods
#' @param kernel_method kernel methods to calculate relationship matrix for the
#' different omics. Currently available are Gaussian kernel, exponential kernel
#' Polynomia kernel (order 2, 3, 4) and linear
#' @param omic1_kernel  relationship matrix using different kernel methods
#' @param omic2_kernel  relationship matrix using different kernel methods
#' @param omic3_kernel  relationship matrix using different kernel methods
#' @param pheno_data_train phenotypic data for training set. NA not allowed. Dataframe or matrix is allowed
#' @param pheno_data_test phenotypic data for the testing set. NA is allowed. Dataframe or matrix is allowed
#' @param cross_validation
#' @param coefficient_1 coefficient for the training set using either genomic or any omics data.
#' We allowed up to 4 omics data for model fit
#' @param coefficient_2
#' @param coefficient_3
#' @param coefficient_4
#' @param eval_metrics
#' @param para_tunning
#' @param var_cov_str user defined variance-covariance structure
#' @param engine if user has asreml
#' @param workspace allocate memory for asreml model fit
#' @param pworkspace allocate memory for predict function in asreml
#' @param bending this is important when the relationship matrix is not positive definitive. It fix it for the user. it has be TRUE
#' @param maxit number of iteration for asreml
#' @param pedigree_matrix
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param rcn_cutoff
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param pheno_data
#' @param pheno_data_train
#' @param pheno_data_test
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param omics_data_label
#' @param gmatrix
#' @param gkernel
#' @param pedigree_matrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param omics_kernel_label
#' @param train_geno_data
#' @param train_omic1_data
#' @param train_omic2_data
#' @param train_omic3_data
#' @param train_omics_label
#' @param test_geno_data
#' @param test_omic1_data
#' @param test_omic2_data
#' @param test_omic3_data
#' @param test_omics_label
#' @param coefficient_1
#' @param coefficient_2
#' @param coefficient_3
#' @param coefficient_4
#' @param train_set
#' @param test_set
#' @param gmatrix_method
#' @param kernel_method
#' @param response
#' @param gen_name
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param var_cov_str
#' @param weights
#' @param nIter
#' @param burnIn
#' @param thin
#' @param GS_model
#' @param eval_metrics
#' @param para_tunning
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param core
#' @param engine
#' @param workspace
#' @param pworkspace
#' @param maxit
#' @param bending
#' @param bend_value
#' @param blending
#' @param blending_value
#' @param high_diag_cut_off
#' @param low_diag_cut_off
#' @param duplicate_cut_off
#' @param rcn_cutoff
#' @param optimize_diagonal
#' @param optimize_duplicate
#' @param message
#' @param system_database this dictate if the output will be created in a folder or as list
#'                         the default is FALSE. Thus output will be folder.
#' @param ...
#' @param scale
#' @param inverse
#' @param epsilon
#' @param vcf_file_name
#' @param vcf_file_path
#' @param vcf_file
#' @param hapmap_file_name
#' @param hapmap_file_path
#' @param hapmap
#' @param maf_threshold
#' @param het_threshold
#' @param ind_call_rate_threshold
#' @param snp_call_rate_threshold
#' @param impute
#' @param recode_format
#' @param out_put_map
#' @param map_data
#' @param qc_filtering
#' @param xgb_paras_tunning
#' @param rf_paras_tunning
#' @param pls_paras_tunning
#' @param svm_paras_tunning
#' @param knn_paras_tunning
#' @param lasso_paras_tunning
#' @param rr_paras_tunning
#' @param dpl_paras_tunning
#'
#' @return
#' @export
#'
#' @examples
#'
model_execute <- function(
    pheno_data = NULL,
    # pheno_file_name = NULL,
    # pheno_file_path= NULL,
    pheno_data_train = NULL,
    pheno_data_test = NULL,
    geno_data = NULL,
    omic1_data = NULL,
    omic2_data = NULL,
    omic3_data = NULL,
    omics_data_label = list(omic1_data = NULL,
                            omic2_data = NULL,
                            omic3_data = NULL),
    gmatrix= NULL,
    gkernel = NULL,
    pedigree_matrix = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    omics_kernel_label = list(omic1_kernel = NULL,
                              omic2_kernel = NULL,
                              omic3_kernel = NULL),
    train_geno_data = NULL,
    train_omic1_data = NULL,
    train_omic2_data = NULL,
    train_omic3_data = NULL,
    train_omics_label = list(train_omic1_data = NULL,
                             train_omic2_data = NULL,
                             train_omic3_data = NULL),
    test_geno_data = NULL,
    test_omic1_data = NULL,
    test_omic2_data = NULL,
    test_omic3_data = NULL,
    test_omics_label = list(test_omic1_data = NULL,
                            test_omic2_data = NULL,
                            test_omic3_data = NULL),
    impute_omic = FALSE,
    coefficient_1 = NULL,
    coefficient_2 = NULL,
    coefficient_3 = NULL,
    coefficient_4 = NULL,
    train_set = NULL,
    test_set = NULL,
    gmatrix_method = NULL,
    kernel_method = NULL,
    response=NULL,
    gen_name=NULL,
    cova=NULL,
    fixed=NULL,
    random=NULL,
    heter_resid=FALSE,
    heter_groups=NULL,
    var_cov_str = NULL,
    weights =NULL,
    nIter=NULL,
    burnIn=NULL,
    thin=NULL,
    GS_model = NULL,
    eval_metrics = NULL,
    fixed_term_model_bayesian = 'FIXED',
    rand_term_model_bayesian = NULL,
    #core = NULL,
    engine = NULL,
    scaling = TRUE,
    centering = FALSE,
    workspace = 1e08,
    pworkspace= 1e06,
    maxit = 50,
    inverse = TRUE,
    epsilon = 1e-6,
    bending = TRUE,
    bend_value = 0.01,
    blending = FALSE,
    blending_value = 0.02,
    high_diag_cut_off = 1.2,
    low_diag_cut_off = 0.8,
    duplicate_cut_off = 0.95,
    rcn_cutoff = 1e-12,
    optimize_diagonal = FALSE,
    optimize_duplicate = FALSE,
    vcf_file_name = NULL,
    vcf_file_path = NULL,
    vcf_file = NULL,
    hapmap_file_name = NULL,
    hapmap_file_path = NULL,
    hapmap = NULL,
    maf_threshold = 0.01,
    het_threshold = 0.1,
    ind_call_rate_threshold = 0.9,
    snp_call_rate_threshold = 0.9,
    impute = FALSE,
    recode_format = "0,1,2",  # Specify "0,1,2" or -1, 0, 1
    out_put_map = FALSE,
    map_data = NULL,
    qc_filtering = TRUE,
    message= TRUE,
    system_database = FALSE,
    num_hidden_layers = 1,
    neurons_per_layer = NULL,
    epochs = 10,
    batch_size = 32 ,
    para_tunning = FALSE,
    param_grid = NULL,
    validation_split = 0.2,
    early_stop = TRUE,
    xgb_paras_tunning= list(Iter_tune = seq(500, 5000, 500), # number of boosting iterations
                         learning_rate_tune = c(0.01, 0.05, 0.1), # learning rate, low value means model is more robust to overfitting
                         max_depth = c(3, 6, 9),
                         gamma = c(0, 0.01, 0.1),
                         colsample_bytree = c(0.5, 0.75, 1),
                         min_child_weight = c(1, 3, 5),
                         subsample = c(0.5, 0.75, 1),
                         L2_tune = c(0, 0.5, 1), #  for linear gbL2 Regularization (Ridge Regression)
                         L1_tune = c(0, 0.5, 1)),
    rf_paras_tunning = list(mtry = TRUE,
                            ntree = c(500, 1000, 1500),
                            nodesize = c(1, 5, 10),
                            maxnodes = c(30, 50, NULL)),
    pls_paras_tunning= list(ncomp = 10),
    svm_paras_tunning= list(
      kernel = c("radial", "linear", "polynomial"),
      #cost = 10^seq(-2, 2, by = 1),
      sigma = c(0.01, 0.05, 0.1),
      C = c(1, 10, 100), ## for radial kernel
      degree = c(2, 3, 4),  # Default values, used only for polynomial
      scale = c(0.1, 1) # used only for polynomial
    ),
    knn_paras_tunning = list(k = seq(3, 21, by = 2),
                             weight = c("uniform", "distance"),
                             metric = c("euclidean", "manhattan")),
    k = 5,
    lasso_paras_tunning= list(lambda_tune=seq(0.000001,0.9,length.out=100)^4),
    rr_paras_tunning = NULL,
    dpl_paras_tunning = NULL,
    learning_rate = 0.01, #xgboost
    max_depth = 6, #xgboost
    subsample = 0.5, #xgboost
    xgb_booster =  "dart", #"gbtree", # #xgboost "gblinear",
    iteration = 1000, #xgboost
    N_feature_impo = 10, #xgboost
    resample_method_tune = "cv", # c("cv","boot") #xgboost
    number_of_fold_tune = 5, #xgboost
    min_child_weight = 0.8, # xgboost,
    #eta = 0.001, ## xgboost
    #nrounds = 5000, ## xgboost
    colsample_bytree = 1, ## xgboost
    xgb_alpha = 0.001, ## xgboost linear
    xgb_gamma = 0.01, ## xgboost it acts as a regularization parameter for controlling tree complexity
    lambda_rr = NULL,
    xgb_lambda = 1.0,  # xgboost linear
    xgb_rate_drop = 0.1,
    xgb_skip_drop = 0.5,
    xgb_objective = "reg:squarederror",
    xgb_sample_type = "uniform",
    xgb_normalize_type = "tree",
    ntree=500, ## RF
    mtry = NULL, ## RF
    maxnodes = NULL, ## RF
    importance=TRUE, ## RF
    ncomp = 3, # pls
    svm_kernel = "Gaussian", #svm "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
    svm_type = "eps-regression",
    sigma_value  = 0.1,       #svm Default sigma value for RBF kernel
    C_value  = 1,             #svm Default cost parameter
    degree_value = 3,        #svm Default degree for polynomial kernel
    scale_value  = 1,         #svm Default scale for polynomial kernel
    offset_value = 0,
    AI_cv_nfolds = 5,
    n_bootstrap = 100,
    early_stop_for_iteration_xgb = FALSE,
    CI_width_thresholds = c(0.33, 0.66),
    confidence_level = 0.95,
    high_reliability_thres = 0.9,
    low_reliability_thres = 0.5,
    abs_very_close_threshold = 0.01,
    abs_close_threshold = 0.05,
    n_components = 20,
    threshold = 100,
    iqr_multiplier = 1.5,
    interval_width_high_threshold = NULL,
    interval_width_moderate_threshold = NULL,
    interval_width_low_threshold = NULL,
    cross_validation = FALSE,
    cv_evaluation_only = TRUE,
    GS_model_cv = NULL,
    nfolds = 5,
    sampling_method = NULL,
    num_cores = NULL,
    replication = 1,
    test_size = 0.3,
    cross_validation_meth = "Stratified_Hold_Out",
    random_state = 123,
    metric_for_ranking = "accuracy",
    plot_extension = "pdf",
    plot_width = 17,
    plot_height = 12,
    plot_units = "in",
    plot_dpi = 300,
    plot_filename = "trait",
    Plot_name_result_diagnostic = NULL,

    ...
) {

#browser()
    msg <- "\n==================================================\n"


    eval_metrics_available <- c("accuracy", "mean_squared_error", "bias",
                                "root_mean_squared_error", "relative_squared_error",
                                "mean_absolute_error", "mean_absolute_percent_error", "kendalls_tau")
    if(!is.null(eval_metrics)){
      if (!all(eval_metrics %in% eval_metrics_available)) {
        stop(paste(msg, "Invalid evaluation metrics for the model. Choose from: ",
             paste(eval_metrics_available, collapse = ", ")), call. = FALSE)
      }
    }

    # Define available models and variance structures
    var_cov_str_available <- c("us","corgh","corgv",
                               "corh","corv","fa1","fa2", "fa3", "fa4",
                               "rr1","rr2", "rr3", "rr4")

    AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare",
                         "SupportVectorMachine", "K-NearestNeighbors", "Lasso",
                         "Ridge_Regression", "deep_learning_model")
    bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
    bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")

    asreml_model <- "GBLUP"


    all_models_avail <- c(AI_valid_models, bayes_valid_models,
                          bayes_gblup_valid_models, asreml_model)

#####
    holds_out_methods_avail <- c("Hold_Out",
                                 "Stratified_Hold_Out",
                                 "Repeated_Hold_Out",
                                 "Repeated_Stratified_Hold_Out",
                                 "Leave_one_Out")

    Kfolds_methods_avail <- c("K-Folds",
                              "Stratified_K-Folds",
                              "Repeated_K-Folds",
                              "Repeated_Stratified_K-Folds")

    CVs_multi_envs_methods_avail <- c("CV1",
                                      "CV2",
                                      "Repeated_CV1",
                                      "Repeated_CV2")
    all_cv_methods_avail <- c(holds_out_methods_avail,
                              Kfolds_methods_avail,
                              CVs_multi_envs_methods_avail)

    if(is.null(heter_resid)) heter_resid <- FALSE
    ### Check for executing cross_validation
    if(isTRUE(cross_validation)) {
      if(is.null(GS_model_cv) || is.null(cross_validation_meth)) {
        stop(paste(msg, "GS_model_cv and cross_validation_meth cannot be null when cross_validation is TRUE"), call. = FALSE)
      }

        if(!all(GS_model_cv%in%all_models_avail)){
          stop(paste(msg,"Invalid model. Choose from: ",
               paste(all_models_avail, collapse = ", ")), call. = FALSE)
        }

      if(length(cross_validation_meth)>1){
        stop(paste(msg,'use only one cross_validation method at a time'), call. = FALSE)
      }

      if(!cross_validation_meth%in%all_cv_methods_avail){
        stop(paste(msg,"Invalid cross validation method. Choose from: ",
             paste(all_cv_methods_avail, collapse = ", ")), call. = FALSE)
      }

      if(is.null(eval_metrics)){
        stop(paste(msg,"Provide evaluation metrics for models comparison. Choose from: ",
             paste(eval_metrics_available, collapse = ", ")), call. = FALSE)
      }

      ###
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
          stop(paste(msg,"You select repeated cross-validation provide number of replications.\n For example, replication = 2"), call. = FALSE)
        }


      }

    }



    # if(!is.null(var_cov_str)){
    # if (!(var_cov_str %in% var_cov_str_available)) {
    #   stop("Invalid output variance-covariance structure. Choose from: ",
    #        paste(var_cov_str_available, collapse = ", "), call. = FALSE)
    #   }
    # }
#######################################



    kernel_method_avaliable <- c("Gaussian_kernel",
                                 "Linear_kernel",
                                 "Composite_kernel",
                                 "Poly2_kernel",
                                 "Poly3_kernel",
                                 "Poly4_kernel")

    if(!is.null(kernel_method)){
      if (!(kernel_method %in% kernel_method_avaliable)) {
        stop(paste(msg,"Invalid kernel method. Choose from: ",
             paste(kernel_method_avaliable, collapse = ", ")), call. = FALSE)
      }
    }

    gmatrix_method_available <- c("VanRaden",
                                  "Weighted_VanRaden",
                                  "Yang",
                                  "Epistasis")

    if(!is.null(gmatrix_method)){
      if (!(gmatrix_method %in% gmatrix_method_available)) {
        stop(paste(msg,"Invalid genomic relationship method. Choose from: ",
             paste(gmatrix_method_available, collapse = ", ")), call. = FALSE)
      }
    }
    #############################
    if(is.null(GS_model) & isFALSE(cross_validation)){
      stop(paste(msg, "Genomic prediction model is missing."), call. = FALSE)

    }else{
      if(!isFALSE(cross_validation)){
      if (!any(GS_model_cv %in% c(bayes_valid_models,
                            bayes_gblup_valid_models,
                            asreml_model,
                            AI_valid_models))) {
        stop(paste(msg,"Invalid genomic prediction model. Choose from:\n", paste(c(bayes_valid_models,
                                                                      bayes_gblup_valid_models,
                                                                      asreml_model,
                                                                      AI_valid_models), collapse = ", ")),
           call. = FALSE)
      }
      } else {
      if(!is.null(GS_model) & isFALSE(cross_validation)){
        if (!any(GS_model %in% c(bayes_valid_models,
                                 bayes_gblup_valid_models,
                                 asreml_model,
                                 AI_valid_models))) {
          stop(paste(msg,"Invalid genomic prediction model. Choose from:\n", paste(c(bayes_valid_models,
                                                                           bayes_gblup_valid_models,
                                                                           asreml_model,
                                                                           AI_valid_models), collapse = ", ")),
               call. = FALSE)
        }
      }
    }

    }

    # Check for mandatory phenotypic data. The chain of loop is important to ease of checking
    ## and the pheno_data_train and pheno_data_test are converted to pheno_data to make life eaier for checking
    if (is.null(pheno_data)) {
      if (is.null(pheno_data_train) && is.null(pheno_data_test)) {
        stop(paste(msg, "Phenotypic data is missing."), call. = FALSE)
      } else if (!is.null(pheno_data_train) && is.null(pheno_data_test)) {
        if ((isTRUE(cross_validation) && isFALSE(cv_evaluation_only)) || (isFALSE(cross_validation) && isTRUE(cv_evaluation_only))) {
          stop(paste(msg, "Provide a dataframe of phenotypic data for the testing set with NA in the response variable(s) if the interest is to do actual prediction after cross-validation. Otherwise set cross_validation = TRUE and cv_evaluation_only = TRUE.\n"), call. = FALSE)
        } else {
          if (isTRUE(cross_validation) && isTRUE(cv_evaluation_only)) {
            if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
              pheno_data_train <- as.data.frame(pheno_data_train)
              #message(sprintf("%s The pheno data is grouped or a tibble. We fix convert to data.frame.\n", msg))
              stop(paste(msg, "Your pheno data is grouped or a tibble. Convert data.frame.\n"), call. = FALSE)
            }
            pheno_data <- pheno_data_train
          }
        }
      } else if (is.null(pheno_data_train) && !is.null(pheno_data_test)) {
        stop(paste(msg, "Provide a dataframe of phenotypic data for the training set."), call. = FALSE)
      } else {
        if (!is.null(pheno_data_train) && !is.null(pheno_data_test)) {
          if (!identical(colnames(pheno_data_train), colnames(pheno_data_test))) {
            stop(paste(msg, "Column names do not match in the pheno_data_train and pheno_data_test.\n"), call. = FALSE)
          } else {
            if (length(as.character(pheno_data_train[[gen_name]])) > length(as.character(unique(pheno_data_train[[gen_name]]))) &&
                length(as.character(pheno_data_test[[gen_name]])) > length(as.character(unique(pheno_data_test[[gen_name]])))) {
              if (is.null(heter_groups)) {
                stop(paste(msg, paste(
                  "Your phenotypic data has a multi-environment structure,",
                  "but the column containing the environment/location is missing.",
                  "Please provide heter_groups parameter.",
                  "For example: heter_groups = 'locations'.",
                  "If you have 'location' as a column name in your phenotypic data."
                )), call. = FALSE)
              }

              # Step 1: Check if the length of unique values in heter_groups matches
              train_unique <- unique(as.character(pheno_data_train[[heter_groups]]))
              test_unique <- unique(as.character(pheno_data_test[[heter_groups]]))

              if (length(train_unique) != length(test_unique)) {
                stop(paste(msg, "Number of unique values in heter_groups is not the same."), call. = FALSE)
              }

              # Step 2: Check if the names of unique values in heter_groups are identical
              if (!identical(train_unique, test_unique)) {
                stop(paste(msg, "Names of unique values in heter_groups do not match."), call. = FALSE)
              }

              if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
                pheno_data_train <- as.data.frame(pheno_data_train)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ###
              if (inherits(pheno_data_test, "tbl_df") || inherits(pheno_data_test, "grouped_df") || inherits(pheno_data_test, "tbl")) {
                pheno_data_test <- as.data.frame(pheno_data_test)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ##

              pheno_data <- rbind(pheno_data_train, pheno_data_test)

            } else if (length(as.character(pheno_data_train[[gen_name]])) == length(as.character(unique(pheno_data_train[[gen_name]]))) &&
                       length(as.character(pheno_data_test[[gen_name]])) == length(as.character(unique(pheno_data_test[[gen_name]])))) {

              if (inherits(pheno_data_train, "tbl_df") || inherits(pheno_data_train, "grouped_df") || inherits(pheno_data_train, "tbl")) {
                pheno_data_train <- as.data.frame(pheno_data_train)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }
              ###
              if (inherits(pheno_data_test, "tbl_df") || inherits(pheno_data_test, "grouped_df") || inherits(pheno_data_test, "tbl")) {
                pheno_data_test <- as.data.frame(pheno_data_test)
                message(paste(msg, 'The pheno data is grouped or a tibble. We fix convert to data.frame.\n'))
                #stop("Your pheno data is grouped or a tibble. Convert data.frame")
              }

              pheno_data <- rbind(pheno_data_train, pheno_data_test)

            } else {
              stop(paste(msg, "The pheno_train_data and pheno_test_data do not match.\n"), call. = FALSE)
            }
          }
        }
      }
    }

    ######
    if(is.null(gen_name)) stop(paste(msg, "Provide gen_name."), call. = FALSE)
    # Check for gen_name presence
    if(!is.null(pheno_data)){
    if (!gen_name %in% colnames(pheno_data)) {
      stop(sprintf("The specified column '%s' in the pheno_data did not match with your data. Please check and use appropriately.", gen_name))
    }

    }

    if (inherits(pheno_data, "tbl_df") || inherits(pheno_data, "grouped_df") || inherits(pheno_data, "tbl")) {
      pheno_data <- as.data.frame(pheno_data)
      message(paste(msg,"The pheno data is grouped or a tibble. We fix convert to data.frame"))
      #stop("Your pheno data is grouped or a tibble. Convert data.frame")
    }
    # else{
    #   if(!is.null(pheno_data_train) & !is.null(pheno_data_test)){
    #     if (!gen_name %in% colnames(pheno_data_train) & !gen_name %in% colnames(pheno_data_test)) {
    #       stop(sprintf("The specified column '%s' in the pheno_data did not match with your data. Please check and use appropriately.", gen_name))
    #     }
    #   }
    # }

    ## Check for scenrio where user provide pheno_train and pheno_test.
    ## Corresponding train and test geno or omic data must be provided

    # Define the error message
    error_message <- "When providing both pheno_data_train and pheno_data_test, at least one of the following pairs must be provided:
                  (train_geno_data and test_geno_data),
                  (train_omic1_data and test_omic1_data),
                  (train_omic2_data and test_omic2_data),
                  (train_omic3_data and test_omic3_data)."

    # Check if GS_model or GS_model_cv contains valid models
    valid_models_using_X_variables <- c(AI_valid_models, bayes_valid_models)
    if ((!is.null(GS_model) && any(GS_model %in% valid_models_using_X_variables)) ||
        (!is.null(GS_model_cv) && any(GS_model_cv %in% valid_models_using_X_variables))) {

      # Check if pheno_data_train and pheno_data_test are provided
      if (!is.null(pheno_data_train) && !is.null(pheno_data_test)) {

        # Check if at least one of the required pairs is provided
        conditions_met <- sum(
          !is.null(train_geno_data) && !is.null(test_geno_data),
          !is.null(train_omic1_data) && !is.null(test_omic1_data),
          !is.null(train_omic2_data) && !is.null(test_omic2_data),
          !is.null(train_omic3_data) && !is.null(test_omic3_data)
        )

        # If none of the pairs are provided, stop with an error message
        if (conditions_met < 1) {
          stop(paste(msg,error_message), call. = FALSE)
        }
      }
    } ### End

    ## Check for when both train and test set data are present in single file
    ##############
    # Define conditions
    condition1 <- is.null(geno_data) && is.null(omic1_data) && is.null(omic2_data) && is.null(omic3_data)
    condition1_1 <- is.null(gmatrix_method) && is.null(kernel_method)

    condition2 <- is.null(gmatrix) && is.null(gkernel) && is.null(omic1_kernel) && is.null(omic2_kernel) && is.null(omic3_kernel)

    ####
    # Check if all required data sets are null
    condition11 <- is.null(geno_data) && is.null(train_geno_data) && is.null(test_geno_data)
    condition12 <- is.null(omic1_data) && is.null(train_omic1_data) && is.null(test_omic1_data)
    condition13 <- is.null(omic2_data) && is.null(train_omic2_data) && is.null(test_omic2_data)
    condition14 <- is.null(omic3_data) && is.null(train_omic3_data) && is.null(test_omic3_data)

    error_message <- paste(msg, "Bayes Alphabets and machine learning models require omics or geno data.")
    # Check if GS_model is not null and belongs to valid models
    if (!is.null(GS_model) && any(GS_model %in% c(AI_valid_models, bayes_valid_models))) {
      # Check if all conditions are true
      if (condition11 && condition12 && condition13 && condition14) {
        stop(paste(msg,error_message), call. = FALSE)
      }
    }

    if (!is.null(GS_model_cv) && any(GS_model_cv %in% c(AI_valid_models, bayes_valid_models))) {
      # Check if all conditions are true
      if (condition11 && condition12 && condition13 && condition14) {
        stop(error_message, call. = FALSE)
      }
    }
    ###################

    # Check for ASReml requirement for GBLUP
    ## Define error message
    error_message <- "ASReml software is required to fit GBLUP for single or multi-environment.\n"

    ### Check that when asreml is used required format for omic data is provided
    if(!is.null(engine)){
      # Create a list of omics data and kernel
      omics <- list(omic1_data, omic2_data, omic3_data, geno_data)
      omics_kernel <- list(omic1_kernel, omic2_kernel, omic3_kernel, gmatrix, gkernel)

      # Remove NULL elements from the list
      omics <- Filter(Negate(is.null), omics)
      omics_kernel <-  Filter(Negate(is.null), omics_kernel)
      # Define a function to check if a matrix is square
      is_square_matrix <- function(mat) {
        if (!is.matrix(mat)) {
          mat <- as.matrix(mat)
        }
        return(nrow(mat) == ncol(mat))
      }

      # Apply the function to each element of the list to check if it's a square matrix
      if (length(omics) != 0) {
        square_matrices <- sapply(omics, is_square_matrix)
        if (any(square_matrices)) {
          stop(paste(msg, "Omic data or geno data is a square matrix. If this is a relationship matrix, provide it as: omic_kernel or gmatrix.\n"), call. = FALSE)
          #stop("Omic data or geno data is a square matrix. If this is a relationship matrix, provide it as: omic_kernel or gmatrix.\n", call. = FALSE)
        } else {
          # Check if kernel_method is NULL
          if (is.null(kernel_method) & any(square_matrices)) {
            stop(paste(msg,"Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n"), call. = FALSE)
            #stop("Provide Kernel method to calculate relationship matrix for omic data to fit GBLUP model.\n", call. = FALSE)
          }
        }

      }
      ###
      if (length(omics_kernel) != 0) {
        square_matrices <- sapply(omics_kernel, is_square_matrix)
        if (!all(square_matrices)) {
          stop(paste(msg, "Omic kernel or gmatrix or gkernel should be a relationship matrix.\n"), call. = FALSE)
          #stop("Omic kernel or gmatrix or gkernel should be a relationship matrix.\n", call. = FALSE)
        }

      }
    if(isFALSE(cross_validation) & !is.null(GS_model)){
    if (any(GS_model == "GBLUP") && engine != "asreml") {
      stop(paste(msg,error_message), call. = FALSE)
     }
    } else {
      if(!is.null(GS_model_cv)){
      if ("GBLUP" %in%GS_model_cv && engine != "asreml") {
        stop(paste(msg,error_message), call. = FALSE)
      }
      }
    }

    } else {
     if(isFALSE(cross_validation) & !is.null(GS_model)){
      if(is.null(engine) & any(GS_model == "GBLUP")){
        stop(paste(msg,error_message), call. = FALSE)
      }

     } else {
       if(!is.null(GS_model_cv)){
       if(is.null(engine) & "GBLUP" %in%GS_model_cv){
         stop(paste(msg,error_message), call. = FALSE)
       }

       }
    }

    }


    # Check for multi-environment structure and required inputs for GBLUP
    if (length(pheno_data[[gen_name]]) > length(unique(pheno_data[[gen_name]]))){
      if ((any(!is.null(GS_model_cv) & GS_model_cv %in% bayes_valid_models)) ||
          (any(!is.null(GS_model) & GS_model %in% bayes_valid_models))) {
        stop(paste(
          msg,
          "Your phenotypic data has a multi-environment structure,",
          "thus, use GBLUP, RKHS, or GBLUP_BRR model.\n"
        ), call. = FALSE)
      }

      if (is.null(heter_groups)) {
        stop(paste(msg,paste(
          "Your phenotypic data has a multi-environment structure,",
          "but the column containing the environment/location is missing.",
          "Please provide heter_groups parameter.",
          "For example: heter_groups = 'locations'.",
          "If you have location as column name in your phenotypic data.\n"
        )), call. = FALSE)
      }

      if (!any(cross_validation_meth %in% CVs_multi_envs_methods_avail)) {
        stop(paste(msg, paste(
          "Provide any of the following as cross-validation strategy for multi-environment GS:",
          paste(CVs_multi_envs_methods_avail, collapse = ", ")
        )), call. = FALSE)
      }

      ### This is important for asreml for multi-environment analysis
      # Define the error messages
      missing_var_cov_str_msg <- paste(msg, "Your data suggest multi-environment but variance-covariance structure is missing. Choose from:", paste(var_cov_str_available, collapse = ", "), call. = FALSE)
      missing_heter_resid_msg <- paste(msg, "Your data suggest multi-environment but variance-covariance structure. heter_resid must be TRUE.", call. = FALSE)
      invalid_var_cov_str_msg <- paste(msg, "Invalid output variance-covariance structure. Choose from:", paste(var_cov_str_available, collapse = ", "), call. = FALSE)

      # Check GS_model
      if (!is.null(GS_model) && any(GS_model %in% c("GBLUP"))||
          !is.null(GS_model_cv) && any(GS_model_cv %in% c("GBLUP"))) {
        if (!is.null(heter_groups) && isTRUE(heter_resid) && is.null(var_cov_str)) {
          stop(missing_var_cov_str_msg)
        } else if (!is.null(heter_groups) && isFALSE(heter_resid) && !is.null(var_cov_str)) {
          stop(missing_heter_resid_msg)
        } else if (!is.null(var_cov_str) && !(var_cov_str %in% var_cov_str_available)) {
          stop(invalid_var_cov_str_msg)
        }
      }

      # Check GS_model_cv for cross-validation
      if (isTRUE(cross_validation) && !is.null(GS_model_cv) && any(GS_model_cv %in% c("GBLUP"))) {
        if (!is.null(heter_groups) && isTRUE(heter_resid) && is.null(var_cov_str)) {
          stop(missing_var_cov_str_msg)
        } else if (!is.null(heter_groups) && isFALSE(heter_resid) && !is.null(var_cov_str)) {
          stop(missing_heter_resid_msg)
        } else if (!is.null(var_cov_str) && !(var_cov_str %in% var_cov_str_available)) {
          stop(invalid_var_cov_str_msg)
        }
      }

      # Check for required inputs for multi-environment GBLUP
      # Check if condition1 is true
      # Define the error message
      error_messagee <- paste(
        msg,
        paste("To fit a Bayesian or ASReml multi-environment GBLUP model,",
        "you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel),",
        "or an omics kernel. Additionally, you can provide genomic or omics data.",
        "Ensure you provide instructions on the genomic relationship matrix method",
        "or kernel method to calculate the relationship matrix.\n"
      ))

      # Check GS_model
      if (!is.null(GS_model) && any(GS_model %in% c(bayes_gblup_valid_models, "GBLUP"))) {
        if (condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }

      # Check GS_model_cv for cross-validation
      if (isTRUE(cross_validation) && !is.null(GS_model_cv) && any(GS_model_cv %in% c(bayes_gblup_valid_models, "GBLUP"))) {
        if (condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }

      # if (GS_model | GS_model_cv %in% c(bayes_gblup_valid_models, "GBLUP")){
      #   if (condition1!=condition1_1 & isTRUE(condition2)) {
      #     #print('ok')
      #     stop(paste(msg, "To fit a Bayesian or ASReml multi-environment GBLUP model, you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel), or an omics kernel. Additionally, you can provide genomic or omics data. Ensure you provide instructions on the genomic relationship matrix method or kernel method to calculate the relationship matrix."), call. = FALSE)
      #   }
      #
      # }

    } else {

    # Check for single environment GBLUP and Bayesian models
    if (length(pheno_data[[gen_name]]) == length(unique(pheno_data[[gen_name]]))){
 #### In case user erroneously provide this while it is a single location
      if(!is.null(heter_groups)) heter_groups <- NULL
      if(!is.null(heter_resid)) heter_resid <- NULL
      if(!is.null(var_cov_str)) var_cov_str <- NULL
      if(cross_validation_meth%in%CVs_multi_envs_methods_avail){
        stop(paste(msg, paste("You selected cross validation method for multi-environment",
                   "but your phenotypic data is a single environment.\n")), call. = FALSE)
      }
      # Check for required inputs for Bayesian or ASReml single environment GBLUP models
      # Define the error message
      error_messagee <- paste(
        msg,
        paste("To fit a Bayesian or ASReml single-environment GBLUP model,",
        "you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel),",
        "or an omics kernel. Additionally, you can provide genomic or omics data.",
        "Ensure you provide instructions on the genomic relationship matrix method",
        "or kernel method to calculate the relationship matrix."
      ))

      # Check GS_model
      if (!is.null(GS_model) && any(GS_model %in% c(bayes_gblup_valid_models, "GBLUP"))) {
        if (condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }

      # Check GS_model_cv for cross-validation
      if (isTRUE(cross_validation) && !is.null(GS_model_cv) && any(GS_model_cv %in% c(bayes_gblup_valid_models, "GBLUP"))) {
        if (condition1 != condition1_1 && isTRUE(condition2)) {
          stop(error_messagee, call. = FALSE)
        }
      }

      # Check for required inputs for Bayesian models
      if(isTRUE(cross_validation)){
        if(any(GS_model_cv%in%c(bayes_valid_models, AI_valid_models)))
          if (isTRUE(condition1)) {
            #print('ok')
            stop(paste(msg, paste("To fit a Bayesian or machine learning model,",
                       "provide genomic or omics data.")), call. = FALSE)
          }
      } else{
      if (any(GS_model %in% c(bayes_valid_models, AI_valid_models))){
        if (isTRUE(condition1)) {
          #print('ok')
          stop(paste(msg, paste("To fit a Bayesian or machine learning model,",
                     "provide genomic or omics data.")), call. = FALSE)
        }

      }

    }

    }

    }
    ##

    # Main Script Not use again remove main scripts
    # checkForASReml(engine, GS_model, GS_model_cv, cross_validation, msg)
    # validateMultiEnvironment(pheno_data, gen_name, heter_groups, heter_resid, var_cov_str, GS_model, var_cov_str_available,cross_validation, GS_model_cv, msg)
    # validateModelRequirements(GS_model, GS_model_cv, bayes_gblup_valid_models, condition1, condition1_1, condition2, cross_validation, msg)

### Check phenotype_to_model for details
 #    This serve as gateway between phenotype-precheck function and readiness of
 #    the phenotypic data for model fitting.

 pheno_clean <- phenotype_to_model(pheno_data = pheno_data,
                                   pheno_data_train = pheno_data_train,
                                   pheno_data_test = pheno_data_test,
                                   train_set = train_set,
                                   test_set = test_set,
                                   response = response,
                                   gen_name = gen_name,
                                   heter_groups = heter_groups,
                                   random = random,
                                   fixed = fixed)

 if("test_set"%in%names(pheno_clean)){
   test_set <- pheno_clean[["test_set"]]
 }
 #if(length(pheno_clean)==0) stop("pheno is null")
## pheno_clean is a list that can have one or two elements
 ## One element if only pheno_data is provided
 ## Two elements if pheno_data/pheno_training and pheno_data_testing was provided as input.
 ##

 ## Check if the pheno_data in the pheno_clean is declared model fit
 #if(attr(pheno_clean[[1]], "cleared")!="for_model_fit" && all(class(pheno_clean[[1]])!=c("data.frame"))) {
 if(attr(pheno_clean[["pheno_clean_data"]], "cleared")!="for_model_fit") {
     stop(print(paste(msg,'pheno_data is not phenotype data')), call. = FALSE)

 }

 ######
 if(!is.null(vcf_file) | (!is.null(vcf_file_name) & !is.null(vcf_file_path))){
   geno_data <- vcf_qc_recode(vcf_file_name = vcf_file_name,
                              vcf_file_path = vcf_file_path,
                              vcf_file = vcf_file,
                              maf_threshold = maf_threshold,
                              het_threshold =het_threshold,
                              ind_call_rate_threshold = ind_call_rate_threshold,
                              snp_call_rate_threshold = snp_call_rate_threshold,
                              impute = impute,
                              recode_format = recode_format,
                              out_put_map = out_put_map,
                              message = message)
 } else {
   if(!is.null(hapmap) | (!is.null(hapmap_file_name) & !is.null(hapmap_file_path))){
     geno_data <- hmp_qc_recode(hapmap_file_name = hapmap_file_name,
                                hapmap_file_path = hapmap_file_path,
                                hapmap = hapmap,
                                maf_threshold = maf_threshold,
                                het_threshold =het_threshold,
                                ind_call_rate_threshold = ind_call_rate_threshold,
                                snp_call_rate_threshold = snp_call_rate_threshold,
                                impute = impute,
                                recode_format = recode_format,
                                out_put_map = out_put_map,
                                message = message)
   }

 }

 #### Get the clean geno_data ready for model fit
 ## The geno_to_model function depend on geno_precheck function. The expected
 ## output is clean genomic data with no missing and all QC control is checked.
 ## For details check geno_to_model and geno_precheck function description
 ##
 #if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){


 # Process genomic data
 if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){
     geno_res <- process_geno_data(geno_data = geno_data,
                                   train_geno_data = train_geno_data,
                                   test_geno_data = test_geno_data,
                                   test_set = test_set,
                                   pheno_clean_list = pheno_clean,
                                   train_set = train_set,
                                   gen_name = gen_name,
                                   kernel_method = kernel_method,
                                   gmatrix_method = gmatrix_method,
                                   scale = scale,
                                   map_data = map_data,
                                   maf_threshold = maf_threshold,
                                   het_threshold = het_threshold,
                                   ind_call_rate_threshold = ind_call_rate_threshold,
                                   snp_call_rate_threshold = snp_call_rate_threshold,
                                   impute = impute,
                                   qc_filtering = qc_filtering,
                                   message = message,
                                   heter_groups = heter_groups)

 } else {
     geno_res <-  list()
 }

 #if(length(geno_res)==0) { stop("geno_res is empty")}
 # Process omic1 data
 omic1_res <- process_omic_data(omic_data = omic1_data,
                                train_omic_data = train_omic1_data,
                                test_omic_data = test_omic1_data,
                                kernel_method = kernel_method,
                                pheno_clean_list = pheno_clean,
                                gen_name = gen_name,
                                test_set = test_set,
                                train_set = train_set,
                                message = message,
                                heter_groups = heter_groups,
                                impute_omic = impute_omic)

 # Process omic2 data
 omic2_res <- process_omic_data(omic_data = omic2_data,
                                train_omic_data = train_omic2_data,
                                test_omic_data = test_omic2_data,
                                kernel_method = kernel_method,
                                pheno_clean_list = pheno_clean,
                                gen_name = gen_name,
                                test_set = test_set,
                                train_set = train_set,
                                message = message,
                                heter_groups = heter_groups,
                                impute_omic = impute_omic)

 # Process omic3 data
 omic3_res <- process_omic_data(omic_data = omic3_data,
                                train_omic_data = train_omic3_data,
                                test_omic_data = test_omic3_data,
                                kernel_method = kernel_method,
                                pheno_clean_list = pheno_clean,
                                gen_name = gen_name,
                                test_set = test_set,
                                train_set = train_set,
                                message = message,
                                heter_groups = heter_groups,
                                impute_omic = impute_omic)
 #################
 geno_omic_model_ready_list <- list()
 gmatrix_kernel_model_ready_list <- list()

 if (length(geno_res)!=0 && all(c("gmatrix", "geno_model_ready") %in% names(geno_res))) {
   gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] <- geno_res[["gmatrix"]]
   geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
 } else {
   if (length(geno_res)!=0 && "geno_model_ready" %in% names(geno_res)) {
     geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
   }
 }
 ##
 if(!is.null(gmatrix)){
   gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] <- gmatrix
 } else {
   if (!is.null(gkernel)) {
     gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] <- gkernel
   }
 }

 ####
 if (length(omic1_res)!=0 && all(c("kernel", "omic_model_ready")%in%names(omic1_res))) {
   gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] <- omic1_res[["kernel"]]
   geno_omic_model_ready_list[["omic1_model_ready"]] <- omic1_res[["omic_model_ready"]]
 } else if(length(omic1_res)!=0 && "omic_model_ready"%in%names(omic1_res)){
   geno_omic_model_ready_list[["omic1_model_ready"]] <- omic1_res[["omic_model_ready"]]
 } else if(length(omic1_res)!=0 && "kernel"%in%names(omic1_res)){
       gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] <- omic1_res[["kernel"]]

 } else {
   if (!is.null(omic1_kernel)) {
     gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] <- omic1_kernel
   }
 }

 if (length(omic2_res)!=0 && all(c("kernel", "omic_model_ready")%in%names(omic2_res))) {
   gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] <- omic2_res[["kernel"]]
   geno_omic_model_ready_list[["omic2_model_ready"]] <- omic2_res[["omic_model_ready"]]
 } else if(length(omic2_res)!=0 && "omic_model_ready"%in%names(omic2_res)){
   geno_omic_model_ready_list[["omic2_model_ready"]] <- omic2_res[["omic_model_ready"]]
 } else if(length(omic2_res)!=0 && "kernel"%in%names(omic2_res)){
       gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] <- omic2_res[["kernel"]]
 } else {
   if (!is.null(omic2_kernel)) {
     gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] <- omic2_kernel
   }
 }


 if (length(omic3_res)!=0 && all(c("kernel", "omic_model_ready")%in%names(omic3_res))) {
   gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] <- omic3_res[["kernel"]]
   geno_omic_model_ready_list[["omic3_model_ready"]] <- omic3_res[["omic_model_ready"]]
 } else if(length(omic3_res)!=0 && "omic_model_ready"%in%names(omic3_res)){
   geno_omic_model_ready_list[["omic3_model_ready"]] <- omic3_res[["omic_model_ready"]]
 } else if(length(omic3_res)!=0 && "kernel"%in%names(omic3_res)){
     gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] <- omic3_res[["kernel"]]
 } else {
   if (!is.null(omic3_kernel)) {
     gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] <- omic3_kernel
   }
 }

 # Define a list of kernel variables
 kernel_vars <- c("gmatrix_model_ready",
                  "omic1_kernel_model_ready",
                  "omic2_kernel_model_ready",
                  "omic3_kernel_model_ready")


 # Define a list to store results
 #results_kernel_list <- list()

 # Iterate over each kernel variable
 for (kernel_var in kernel_vars) {
     #checked_var <- paste0(kernel_var, "_checked")
     gmatrix_kernel_model_ready_list[[kernel_var]]
     # Check if the kernel variable exists
     if (!is.null(gmatrix_kernel_model_ready_list[[kernel_var]])) {
         # Pre-check the kernel data
         #results_kernel_list[[kernel_var]] <- grm_kernel_precheck(
       gmatrix_kernel_model_ready_list[[kernel_var]] <- grm_kernel_precheck(
             grm_kernel_data = gmatrix_kernel_model_ready_list[[kernel_var]],
             pedigree_matrix = pedigree_matrix,
             bending = bending,
             bend_value = bend_value,
             blending = blending,
             blending_value = blending_value,
             high_diag_cut_off = high_diag_cut_off,
             low_diag_cut_off = low_diag_cut_off,
             duplicate_cut_off = duplicate_cut_off,
             rcn_cutoff = rcn_cutoff,
             optimize_diagonal = optimize_diagonal,
             optimize_duplicate = optimize_duplicate,
             message = message
         )

         # Remove the intermediate kernel variable
         #rm(list = kernel_var)
     }
 }
 #gmatrix_kernel_model_ready_list <-  results_kernel_list
 #rm(results_kernel_list)
 # Iterate over each checked kernel variable for pheno-geno match
 if(length(gmatrix_kernel_model_ready_list)!=0){
 for (checked_kernel_var_name in names(gmatrix_kernel_model_ready_list)) {
     # Check if the checked kernel variable exists
     if (!is.null(gmatrix_kernel_model_ready_list[[checked_kernel_var_name]])) {
         # Perform pheno-geno match
         match_result <- pheno_geno_match(
             object_pheno = pheno_clean[["pheno_clean_data"]],
             object_geno = gmatrix_kernel_model_ready_list[[checked_kernel_var_name]],
             gen_name = gen_name,
             test_set = test_set,
             train_set = train_set,
             message = message
         )

         # Assign model-ready variable and update test_set if necessary
         #results_list[[paste0(sub("_checked$", "", checked_kernel_var_name), "_model_ready")]] <- match_result[[1]]
         #model_ready_name <- paste0(checked_kernel_var_name, "_model_ready")
         gmatrix_kernel_model_ready_list[[checked_kernel_var_name]] <- match_result[[1]]
         if (length(match_result) > 1) {
             test_set <- match_result[[2]]
             if(is.data.frame(test_set) | is.matrix(test_set)){
               test_set <-  test_set[, 1]
               test_set <-  unique(test_set) ## incase of MET pheno data
             }
         }
     }
 }


 }
### Concatenation of omics for ML
 # When calling the function, pass the external variables as arguments
 if (!is.null(GS_model) & is.null(GS_model_cv)) {
   ml_dat_res <- AI_process_ml_data_if_valid(model_check = any(GS_model %in% AI_valid_models), geno_omic_model_ready_list, pheno_clean, response, gen_name)

 } else if (!is.null(GS_model_cv) && (is.null(GS_model) || !is.null(GS_model))) {
   ml_dat_res <- AI_process_ml_data_if_valid(model_check = any(GS_model_cv %in% AI_valid_models), geno_omic_model_ready_list, pheno_clean, response, gen_name)
 } else {
   ml_dat_res <- list()
 }

 ### Ends
#################################################################
############# Cross-Validation Start

 best_models <-  NULL
 best_models_ggplot_rep <- NULL
 best_models_ggplot_mean <- NULL
 cv_results_processed <-  NULL
 model_prep_all_bayes_cv <-  NULL
 asreml_models_prep_cv <- NULL
 res_plot_result_diagnostic <-  NULL
 res_mod_results_cv_per_trait_model <-  NULL
 cv_results_predicted_vs_observed <- NULL
 res_plot_result_diagnostic_cv_only <- NULL
 diagnostic_plots <- NULL

 if(isTRUE(cross_validation)){
   #model_prep_all_bayes_cv <-  NULL
     if("test_set"%in%names(pheno_clean)) {

       test_set <- pheno_clean[["test_set"]]

     }

   if(!is.null(test_set)){
     if(is.data.frame(test_set) | is.matrix(test_set)){
       test_set <-  test_set[, 1]
       test_set <-  unique(test_set) ## incase of MET pheno data
     }

   }
   pheno_data <-  pheno_clean[["pheno_clean_data"]]
   if(!is.null(test_set)) {

     pheno_data <-  pheno_data[!pheno_data[[gen_name]] %in% test_set, ]

   }

   if(any(GS_model_cv%in% c(bayes_valid_models, bayes_gblup_valid_models))){
 model_prep_all_bayes_cv <- model_prep_bayes_cv(fixed = fixed,
                                               random = random,
                                               GS_model_cv = GS_model_cv,
                                               response = response,
                                               gen_name = gen_name,
                                               pheno_data = pheno_data,
                                               test_set = test_set,
                                               weights = weights,
                                               fixed_term_model_bayesian = fixed_term_model_bayesian,
                                               rand_term_model_bayesian = rand_term_model_bayesian,
                                               nIter = nIter,
                                               burnIn = burnIn,
                                               thin = thin,
                                               geno_data = if("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
                                               omic1_data = if("omic1_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic1_model_ready"]] else NULL,
                                               omic2_data = if("omic2_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic2_model_ready"]] else NULL,
                                               omic3_data = if("omic3_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic3_model_ready"]] else NULL,
                                               omics_data_label = omics_data_label,
                                               gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
                                               omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
                                               omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
                                               omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
                                               heter_groups = heter_groups,
                                               omics_kernel_label = omics_kernel_label,
                                               cross_validation = TRUE)

   }

   if(any(GS_model_cv%in% c("GBLUP"))){
     asreml_models_prep_cv <- asreml_utilis_new( fixed = fixed,
                                                 random = random,
                                                 engine = engine,
                                                 cova= cova,
                                                 GS_model = "GBLUP",
                                                 response = response,
                                                 pheno_data = pheno_data,
                                                 gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
                                                 omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
                                                 omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
                                                 omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
                                                 inverse = inverse,
                                                 epsilon = epsilon,
                                                 gen_name = gen_name,
                                                 heter_groups = heter_groups,
                                                 heter_resid = heter_resid,
                                                 var_cov_str = var_cov_str,
                                                 weights = weights,
                                                 workspace = workspace,
                                                 pworkspace= pworkspace,
                                                 maxit = maxit,
                                                 cross_validation = TRUE)
   }

   cv_results <- models_execute_crossval(pheno_data = pheno_data,
                                        test_set = test_set,
                                        response = response,
                                        gen_name = gen_name,
                                        test_size = test_size,
                                        random_state = random_state,
                                        replication = replication,
                                        heter_groups = heter_groups,
                                        cross_validation_meth = cross_validation_meth,
                                        nfolds = nfolds,
                                        sampling_method = sampling_method,
                                        model_prep_all_bayes_cv = model_prep_all_bayes_cv,
                                        asreml_models_prep_cv = asreml_models_prep_cv,
                                        engine = engine,
                                        ml_dat_res = ml_dat_res,
                                        GS_model_cv = GS_model_cv,
                                        num_cores = num_cores,
                                        eval_metrics = eval_metrics,
                                        scaling = scaling,
                                        centering =centering,
                                        eta = learning_rate, ## xgboost
                                        nrounds = iteration, ## xgboost
                                        max_depth = max_depth, ## xgboost
                                        xgb_gamma = xgb_gamma, ## xgboost
                                        subsample = subsample, ## xgboost
                                        colsample_bytree = colsample_bytree, ## xgboost
                                        xgb_alpha = xgb_alpha, ## xgboost linear
                                        xgb_lambda = xgb_lambda, ## xgboost linear
                                        min_child_weight = min_child_weight, ## xgboost
                                        early_stop_for_iteration_xgb = early_stop_for_iteration_xgb, ## xgboost
                                        xgb_booster = xgb_booster,
                                        ncomp = ncomp, #### pls
                                        ntree = ntree, ### random forest
                                        k = k, ## for knn
                                        svm_kernel = svm_kernel, # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
                                        sigma_value  = sigma_value,       # Default sigma value for RBF kernel
                                        C_value  = C_value,             # Default cost parameter
                                        degree_value = degree_value,        # Default degree for polynomial kernel
                                        scale_value  = scale_value,         # Default scale for polynomial kernel
                                        offset_value = offset_value)

if(cross_validation_meth%in%c("CV1",
                              "CV2",
                              "Repeated_CV1",
                              "Repeated_CV2")){
cv_results_processed <- cv1_cv2_and_across_env_result_plot_process(cv_results_data=cv_results,
                                                                   eval_metrics = eval_metrics,
                                                                   metric_for_ranking = metric_for_ranking)

best_models <- cv_results_processed[["best_models_list"]][[metric_for_ranking]]

best_models_ggplot_rep <- cv_results_processed[["plot_reps_list"]][[metric_for_ranking]][["ggplot_boxplot_reps"]]

best_models_ggplot_mean <- cv_results_processed[["plot_mean_list"]][[metric_for_ranking]][["ggplot_lineplot_mean"]]


} else {

  cv_results_predicted_vs_observed <- single_predicted_vs_observed_result_plots_process(results = cv_results,
                                                                                        pheno_data = pheno_data,
                                                                                        abs_very_close_threshold = abs_very_close_threshold,
                                                                                        abs_close_threshold = abs_close_threshold)

cv_results_processed <- cv_single_loc_result_plot_process(cv_results_data=cv_results,
                                                          eval_metrics = eval_metrics)

best_models <- cv_results_processed[["best_models_list"]][[metric_for_ranking]]

best_models_ggplot_rep <- cv_results_processed[["plot_reps_list"]][[metric_for_ranking]][["ggplot_boxplot_reps"]]

best_models_ggplot_mean <- cv_results_processed[["plot_mean_list"]][[metric_for_ranking]][["ggplot_lineplot_mean"]]


  }

 }

 geno_qc_stat <- if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL


 if(isTRUE(cv_evaluation_only) && isTRUE(cross_validation)){

   return(results_handling(GS_model =  NULL,
                           res_model_output =  NULL,
                           res_summary_stat =  NULL,
                           res_plot = best_models_ggplot_rep,
                           res_plot_mean = best_models_ggplot_mean,
                           res_plot_result_diagnostic = NULL,
                           test_diagonistic_plots = NULL,
                           res_plot_result_diagnostic_cv_only = cv_results_predicted_vs_observed$predicted_vs_observed_plots,
                           res_mod_results_cv_per_trait_model = cv_results_predicted_vs_observed$mod_res_per_trait_per_model,
                           geno_qc_stat = geno_qc_stat,
                           cv_results_processed = cv_results_processed,
                           system_database = system_database,
                           plot_filename = "CV_results",
                           #Plot_name_result_diagnostic = if(!is.null(names(res_mod_results_cv_per_trait_model)[res])) names(results)[res] else paste("trait_diganostic", res, sep = "_"),
                           plot_extension = plot_extension,
                           plot_width = plot_width,
                           plot_height = plot_height,
                           plot_units = plot_units,
                           plot_dpi = plot_dpi))

 }

 if(is.null(best_models)){
   cat("There is a problem with the cross-vlidation process. Check the data and the model\n.")
   return(NULL)
 }

###############################################################

 ##########################################################################
 #########################################################################
 ## Start of Bayes A, B, C, BL and BRR Models for Single Location       ##
 ##  This only accommodate n x p matrix  not nxn                        ##
 ##                                                                     ##
 ##########################################################################
 #######################################################################

 ## sik ############################################################################################ sik


 # geno_qc_stat <- if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL

 #######
 if(!is.null(best_models)){
   n_trait <- length(best_models[["trait"]])
   n_model <- length(best_models[["model"]])

 } else {
   if (length(GS_model) > 1) {
     if (length(GS_model) != length(response)) {
       stop(paste(msg, "When the number of models is more than one, the number of models should be the same as the number of traits."), call. = FALSE)
     }
   }

   #task <- data.frame(model = GS_model, trait = response, stringsAsFactors = FALSE)
   best_models <- data.frame(model = GS_model, trait = response, stringsAsFactors = FALSE)
   n_trait <-  length(response)
   n_model <- length(GS_model)
 }
 # Main logic
 sys_name <- Sys.info()["sysname"]
 if (!is.null(num_cores) && num_cores > 1) {
   #sys_name <- Sys.info()["sysname"]
   set_parallel_plan(n_trait = n_trait, n_model = n_model,
                     sys_name = sys_name)
 } else {
   # Automatically determine the number of cores and use half of them
   detected_cores <- parallel::detectCores(logical = TRUE)
   # For non-Windows systems, consider physical cores only
   num_cores <- round(detected_cores * 0.7)

   set_parallel_plan(n_trait= n_trait, n_model = n_model,num_cores = num_cores,
                     sys_name = sys_name)
 }

 #results_use = results

 results <- future.apply::future_lapply(seq_len(nrow(best_models)), function(i) {
   task_row <- best_models[i, ]

   response <- as.character(task_row$trait)
   GS_model <- as.character(task_row$model)

   if(!GS_model%in%AI_valid_models) {
     model_for_CI_cal <-"Bayes"
   } else{
     model_for_CI_cal <- "ML"
   }

   #### These models only works with one environment/location
   if(length(pheno_clean[["pheno_clean_data"]][,gen_name])==length(unique(pheno_clean[["pheno_clean_data"]][,gen_name]))){

     if ((GS_model %in% bayes_valid_models && is.null(rand_term_model_bayesian)) ||
         (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models)) ||
         (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models))) {

       bayes_A_B_C_BRR_mod_process <- tryCatch({
       res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
                                                       random = random,
                                                       GS_model = GS_model,
                                                       response = response,
                                                       weights = weights,
                                                       fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                       rand_term_model_bayesian = rand_term_model_bayesian,
                                                       pheno_data = pheno_clean[["pheno_clean_data"]],
                                                       geno_data = if("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
                                                       omic1_data = if("omic1_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic1_model_ready"]] else NULL,
                                                       omic2_data = if("omic2_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic2_model_ready"]] else NULL,
                                                       omic3_data = if("omic3_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic3_model_ready"]] else NULL,
                                                       gen_name = gen_name,
                                                       nIter = nIter,
                                                       burnIn = burnIn,
                                                       thin = thin,
                                                       omics_data_label = omics_data_label,
                                                       scaling = scaling,
                                                       CI_width_thresholds = CI_width_thresholds,
                                                       confidence_level = confidence_level,
                                                       high_reliability_thres = high_reliability_thres,
                                                       low_reliability_thres = low_reliability_thres,
                                                       n_components = n_components,
                                                       threshold = threshold,
                                                       target = "test_set",
                                                       confidence_level = confidence_level,
                                                       #iqr_multiplier = iqr_multiplier,
                                                       interval_width_high_threshold = interval_width_high_threshold,
                                                       interval_width_low_threshold = interval_width_low_threshold,
                                                       interval_width_moderate_threshold = interval_width_moderate_threshold)

       res_model_output  # Return the model object if everything is successful
       }, error = function(e) {
         #
         cat("Error:", conditionMessage(e), "\n")
         return(NULL)  # Return NULL
       })
       # Compute summary statistics and plot accuracy
       if(!is.null(bayes_A_B_C_BRR_mod_process)){

         res_model_output <- bayes_A_B_C_BRR_mod_process
         #res_model_output <- res_model_output[["bayes_result"]]

         bayes_summary_stat_process <- tryCatch({
       res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]],
                                                    eval_metrics = eval_metrics,
                                                    model_result = res_model_output[["bayes_result"]],
                                                    GS_model = GS_model,
                                                    gen_name = gen_name,
                                                    CI_width_thresholds = CI_width_thresholds,
                                                    confidence_level = confidence_level,
                                                    high_reliability_thres = high_reliability_thres,
                                                    low_reliability_thres = low_reliability_thres,
                                                    system_database = system_database)
       #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
       res_summary_stat  # Return the model object if everything is successful
         }, error = function(e) {
           # Handle the error, you can print a message or take other actions
           cat("Error:", conditionMessage(e), "\n")
           return(NULL)  # Return NULL or an appropriate value to indicate the failure
         })

         if(!is.null(bayes_summary_stat_process)){
           res_summary_stat <- bayes_summary_stat_process
       if("diagnostic_tst_plot"%in%names(res_summary_stat)){

         res_model_output[["diagnostic_plots"]] <- res_summary_stat[["diagnostic_tst_plot"]]

         res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
       }

         } else{
           cat(sprintf("Bayesian %s summary statistics failed.\n", GS_model))
           res_summary_stat <- NULL
         }

       } else {
         cat(sprintf("Bayesian %s model failed.\n", GS_model))
         res_model_output <- NULL
         res_summary_stat <- NULL
       }

       # output <- list(GS_model = GS_model,
       #                res_model_output = res_model_output,
       #                res_summary_stat = res_summary_stat,
       #                geno_qc_stat =geno_qc_stat
       # )

     }
   } ## End of  Bayes A, B, C, BRR, BL

   ##########################################################################
   #########################################################################
   ## Start of Reproducing Kernel Hilbert Spaces Regression RKHS,         ##
   ## (BRR- Bayesian GBLUP ) and GBLUP (asreml) Model                     ##
   ## for Single Location and multiple loc                                ##
   ##                                                                     ##
   ##                                                                     ##
   ##########################################################################
   #######################################################################

   ## NOTE
   ## BRR is changed to G-BRR to make distinction between BRR for marker matrix and GBLUP
   # Check conditions for GS_model and rand_term_model_bayesian
   if ((GS_model %in% c("RKHS", "GBLUP_BRR", "GBLUP") && is.null(rand_term_model_bayesian)) ||
       (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models)) ||
       (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models))) {

     # Rename GS_model for GBLUP_BRR case
     if (GS_model == "GBLUP_BRR") {
       GS_modeluse <- GS_model
       GS_model <- "BRR"
     }

     if (GS_model %in% c("BRR", "RKHS")) {
       # Run Bayesian model for BRR and RKHS
       bayes_RKHS_GBLUP_BRR_mod_process <- tryCatch({
       res_model_output <- bayes_finalize_RKHS_GBLUPBRR(fixed = fixed,
                                                        random = random,
                                                        GS_model = GS_model,
                                                        response = response,
                                                        weights = weights,
                                                        fixed_term_model_bayesian = fixed_term_model_bayesian,
                                                        rand_term_model_bayesian = rand_term_model_bayesian,
                                                        pheno_data = pheno_clean[["pheno_clean_data"]],
                                                        gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
                                                        omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
                                                        omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
                                                        omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
                                                        gen_name = gen_name,
                                                        nIter = nIter,
                                                        burnIn = burnIn,
                                                        thin = thin,
                                                        heter_groups = heter_groups,
                                                        omics_kernel_label = omics_kernel_label,
                                                        CI_width_thresholds = CI_width_thresholds,
                                                        confidence_level = confidence_level,
                                                        high_reliability_thres = high_reliability_thres,
                                                        low_reliability_thres = low_reliability_thres,
                                                        n_components = n_components,
                                                        threshold = threshold,
                                                        target = "test_set",
                                                        cross_validation = FALSE,
                                                        confidence_level = confidence_level,
                                                        #iqr_multiplier = iqr_multiplier,
                                                        interval_width_low_threshold = interval_width_low_threshold,
                                                        interval_width_high_threshold = interval_width_high_threshold,
                                                        interval_width_moderate_threshold = interval_width_moderate_threshold)

       res_model_output  # Return the model object if everything is successful
       }, error = function(e) {
         #
         cat("Error:", conditionMessage(e), "\n")
         return(NULL)  # Return NULL
       })

       if(!is.null(bayes_RKHS_GBLUP_BRR_mod_process)){
       # Compute summary statistics and plot accuracy
         res_model_output <- bayes_RKHS_GBLUP_BRR_mod_process
         #[["bayes_result"]]
         bayes_GBLUP_summary_stat_process <- tryCatch({
       res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]],
                                                    eval_metrics = eval_metrics,
                                                    model_result = res_model_output[["bayes_result"]],
                                                    GS_model = GS_model,
                                                    gen_name = gen_name,
                                                    CI_width_thresholds = CI_width_thresholds,
                                                    confidence_level = confidence_level,
                                                    high_reliability_thres = high_reliability_thres,
                                                    low_reliability_thres = low_reliability_thres,
                                                    system_database = system_database)
       #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)

       res_summary_stat  # Return the model object if everything is successful
         }, error = function(e) {
           #
           cat("Error:", conditionMessage(e), "\n")
           return(NULL)  # Return NULL
         })

         if(!is.null(bayes_GBLUP_summary_stat_process)){
           res_summary_stat <- bayes_GBLUP_summary_stat_process
           if("diagnostic_tst_plot"%in%names(res_summary_stat)){

             res_model_output[["diagnostic_plots"]] <- res_summary_stat[["diagnostic_tst_plot"]]

             res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
           }
         } else {
           res_summary_stat <- NULL
           cat(sprintf("Bayesian %s summary statistics failed.\n", GS_model))
         }

       } else {
         cat(sprintf("Bayesian %s summary statistics failed.\n", GS_model))

         res_model_output <- NULL
         res_summary_stat <- NULL
       }

       ### This part is for GBLUP_BRR
       # if(exists("GS_modeluse")){
       #   GS_model <-  GS_modeluse
       # }

       # output <- list(GS_model = GS_model,
       #                res_model_output = res_model_output,
       #                res_summary_stat = res_summary_stat,
       #                geno_qc_stat =geno_qc_stat
       # )

     } else {
       if (GS_model == "GBLUP" && engine == 'asreml') {

         result_asreml_mod <- tryCatch({
           # Run GBLUP model with ASReml
           mod <- asreml_utilis_new(
             fixed = fixed,
             random = random,
             cova = cova,
             GS_model = GS_model,
             response = response,
             pheno_data = pheno_clean[["pheno_clean_data"]],
             gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
             omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
             omic2_kernel = if("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
             omic3_kernel = if("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
             gen_name = gen_name,
             heter_groups = heter_groups,
             heter_resid = heter_resid,
             var_cov_str = var_cov_str,
             weights = weights,
             pworkspace = pworkspace,
             workspace = workspace,
             maxit = maxit,
             inverse = inverse,
             epsilon = epsilon,
             cross_validation = FALSE,
             engine = engine
           )

           mod  # Return the model object if everything is successful
         }, error = function(e) {
           # this Handle the error, a message will be printed a message. though other actions can be taking
           cat("Error:", conditionMessage(e), "\n")
           return(NULL)  # Return NULL or an appropriate value to indicate the failure
         })

         if (is.null(result_asreml_mod)) {
           cat("Model fitting failed due to convergence problem.\n")
           res_model_output <- NULL
           res_summary_stat <- NULL
           res_comp_checkk <- NULL
         } else {
           mod <- result_asreml_mod

           # Extract model output for ASReml
           vc <-  summary(mod$model)$varcomp

           res_comp_checkk <- tryCatch({
             VAR_check_Pos <- which(vc$bound == "?" | vc$bound == "S")
             if (length(VAR_check_Pos) > 0) {
               unstable_components <- paste("variance component for", paste(rownames(vc)[VAR_check_Pos], collapse = " and"), "is unstable, refix the model")
               stop(unstable_components, call. = FALSE)
             }
             vc  # Return the variance components if no error
           },
           error = function(e) {
             # this Handle the error, a message will be printed a message. though other actions can be taking
             cat("Error:", conditionMessage(e), "\n")
             return(NULL)  # Return NULL
           })

         }


         if(!is.null(res_comp_checkk)){

    asreml_mod_output_process <- tryCatch({
         res_model_output <- asreml_mod_output_new(
           mod_asreml = mod,
           pheno_data = pheno_clean[["pheno_clean_data"]],
           gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
           omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
           omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
           omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
           heter_groups = heter_groups,
           gen_name = gen_name,
           response = response,
           var_cov_str = var_cov_str,
           heter_resid = heter_resid,
           pworkspace = pworkspace,
           workspace = workspace,
           maxit = maxit
         )

         res_model_output  # Return the model object if everything is successful
    }, error = function(e) {
      # Handle the error, you can print a message or take other actions
      cat("Error:", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or an appropriate value to indicate the failure
    })

    if(!is.null(asreml_mod_output_process)){
      res_model_output <- asreml_mod_output_process
      asreml_summary_stat_process <- tryCatch({
         res_summary_stat <- summary_statistics_asreml(mod =  res_model_output[["Asreml_model"]],
                                                       response = response,
                                                       pheno_data = pheno_clean[["pheno_clean_data"]],
                                                       heter_groups = heter_groups,
                                                       GID_names = res_model_output[["Predicted_value"]][gen_name],
                                                       predicted_value =  if("Predicted_value"%in%colnames(res_model_output[["Predicted_value"]])) res_model_output[["Predicted_value"]]["Predicted_value"] else res_model_output[["Predicted_value"]]["BLUP"],
                                                       standard_errors = res_model_output[["Predicted_value"]]["Standard_error"],
                                                       prediction_error_var = res_model_output[["Predicted_value"]]["Prediction_error_variance"],
                                                       genetic_var = var(res_model_output[["Predicted_value"]]["Predicted_value"]),
                                                       pred_heter_groups = NULL,
                                                       variance_components = res_model_output[["Variance_components"]],
                                                       eval_metrics = eval_metrics,
                                                       gen_name = gen_name,
                                                       CI_width_thresholds = CI_width_thresholds,
                                                       confidence_level = confidence_level,
                                                       high_reliability_thres = high_reliability_thres,
                                                       low_reliability_thres = low_reliability_thres,
                                                       system_database = system_database)

         res_summary_stat  # Return the model object if everything is successful
      }, error = function(e) {
        # Handle the error, you can print a message or take other actions
        cat("Error:", conditionMessage(e), "\n")
        return(NULL)  # Return NULL or an appropriate value to indicate the failure
      })

      if(!is.null(asreml_summary_stat_process)){
        res_summary_stat <- asreml_summary_stat_process
        if("diagnostic_tst_plot"%in%names(res_summary_stat)){

          res_model_output[["diagnostic_plots"]] <- res_summary_stat[["diagnostic_tst_plot"]]

          res_summary_stat <- res_summary_stat[!names(res_summary_stat) %in% "diagnostic_tst_plot"]
        }
      } else{
        cat("Error processing summary statistics for asreml result.\n")
        res_summary_stat <- NULL
      }

    } else {
      cat("Error processing the output of asreml result.\n")
      res_model_output <- NULL
      res_summary_stat <- NULL

    }

         } else {
           res_model_output <- NULL
           res_summary_stat <- NULL

         }

         # output <- list(GS_model = GS_model,
         #                res_model_output = res_model_output,
         #                res_summary_stat = res_summary_stat,
         #                geno_qc_stat =geno_qc_stat
         # )

       }
     }
   }
   #### END GBLUP_RKHS, GBLUP_BRR and GBLUP (asreml)

   ######################################################
   ######################################################
   ##                                                 ###
   ## Machine Learning Models                         ###
   ##                                                 ###
   ######################################################
   ######################################################

   if (GS_model %in% AI_valid_models) {
     if (length(unique(pheno_clean[["pheno_clean_data"]][, gen_name])) > length(pheno_clean[["pheno_clean_data"]][, gen_name])) {
       stop(paste(msg, GS_model, 'only works for single location/enviroment.'), call. = FALSE)
     }

     # print(names(ml_dat_res))
     # if (!"merged_data" %in% names(ml_dat_res) || is.null(ml_dat_res[["merged_data"]][["merge_data"]])) {
     #   stop("merge_data is missing from ml_dat_res")
     # }
     #
     # # Check if the required data is present and correctly formatted
     # if (!"pheno_clean_data" %in% names(ml_dat_res) || is.null(ml_dat_res[["pheno_clean_data"]])) {
     #   stop("pheno_clean_data is missing from ml_dat_res")
     # }


     # switch(GS_model,
     #        "Xgboost" = {
     #          res_model_output <- AI_Xgb(pheno_object = ml_dat_res[["pheno_clean_data"]],
     #                                     response = response,
     #                                     geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #                                     #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #                                     geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #                                     message = message,
     #                                     gen_name = gen_name,
     #                                     scaling = scaling,
     #                                     centering = centering,
     #                                     omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #                                     AI_cv_nfolds = AI_cv_nfolds,
     #                                     para_tunning = para_tunning,
     #                                     xgb_paras_tunning = xgb_paras_tunning,
     #                                     resample_method_tune = resample_method_tune, # c("cv","boot")
     #                                     number_of_fold_tune = number_of_fold_tune,
     #                                     learning_rate = learning_rate,
     #                                     xgb_gamma = xgb_gamma,
     #                                     xgb_lambda = xgb_lambda,
     #                                     xgb_alpha = xgb_alpha,
     #                                     max_depth = max_depth,
     #                                     subsample = subsample,
     #                                     xgb_booster =  xgb_booster, # "gblinear",
     #                                     colsample_bytree = colsample_bytree, ## xgboost
     #                                     alpha = alpha, ## xgboost linear
     #                                     lambda = lambda, ## xgboost linear
     #                                     iteration = iteration,
     #                                     xgb_rate_drop = xgb_rate_drop,
     #                                     xgb_skip_drop = xgb_skip_drop,
     #                                     xgb_objective = xgb_objective,
     #                                     xgb_sample_type = xgb_sample_type,
     #                                     xgb_normalize_type = xgb_normalize_type,
     #                                     early_stop_for_iteration_xgb = early_stop_for_iteration_xgb,
     #                                     N_feature_impo = N_feature_impo,
     #                                     CI_width_thresholds = CI_width_thresholds,
     #                                     high_reliability_thres = high_reliability_thres,
     #                                     low_reliability_thres = low_reliability_thres,
     #                                     n_components = n_components,
     #                                     threshold = threshold,
     #                                     target = "test_set",
     #                                     #iqr_multiplier = iqr_multiplier,
     #                                     interval_width_low_threshold = interval_width_low_threshold,
     #                                     interval_width_high_threshold = interval_width_high_threshold,
     #                                     interval_width_moderate_threshold = interval_width_moderate_threshold,
     #                                     n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "RandomForest" = {
     #          res_model_output <- AI_randomForest(pheno_object = ml_dat_res[["pheno_clean_data"]],
     #                                              response = response,
     #                                              geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #                                              #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #                                              geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #                                              message = message,
     #                                              gen_name = gen_name,
     #                                              scaling = scaling,
     #                                              centering = centering,
     #                                              omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #                                              AI_cv_nfolds = AI_cv_nfolds,
     #                                              para_tunning = para_tunning,
     #                                              rf_paras_tunning = rf_paras_tunning,
     #                                              ntree=ntree,
     #                                              mtry = mtry,
     #                                              maxnodes = maxnodes,
     #                                              importance=importance,
     #                                              CI_width_thresholds = CI_width_thresholds,
     #                                              high_reliability_thres = high_reliability_thres,
     #                                              low_reliability_thres = low_reliability_thres,
     #                                              n_components = n_components,
     #                                              threshold = threshold,
     #                                              target = "test_set",
     #                                              #iqr_multiplier = iqr_multiplier,
     #                                              interval_width_low_threshold = interval_width_low_threshold,
     #                                              interval_width_high_threshold = interval_width_high_threshold,
     #                                              interval_width_moderate_threshold = interval_width_moderate_threshold,
     #                                              n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "PartialLeastSquare" = {
     #          res_model_output <-  AI_pls(pheno_object = ml_dat_res[["pheno_clean_data"]],
     #                                      response = response,
     #                                      geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #                                      #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #                                      geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #                                      message = message,
     #                                      gen_name = gen_name,
     #                                      scaling = scaling,
     #                                      centering = centering,
     #                                      omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #                                      para_tunning = para_tunning,
     #                                      ncomp = ncomp,
     #                                      pls_paras_tunning = pls_paras_tunning,
     #                                      resample_method_tune = resample_method_tune,
     #                                      N_feature_impo = N_feature_impo,
     #                                      CI_width_thresholds = CI_width_thresholds,
     #                                      high_reliability_thres = high_reliability_thres,
     #                                      low_reliability_thres = low_reliability_thres,
     #                                      n_components = n_components,
     #                                      threshold = threshold,
     #                                      target = "test_set",
     #                                      #iqr_multiplier = iqr_multiplier,
     #                                      interval_width_low_threshold = interval_width_low_threshold,
     #                                      interval_width_high_threshold = interval_width_high_threshold,
     #                                      interval_width_moderate_threshold = interval_width_moderate_threshold,
     #                                      n_bootstrap = n_bootstrap)
     #        },
     #        "SupportVectorMachine" = {
     #          res_model_output <- AI_svm(pheno_object = ml_dat_res[["pheno_clean_data"]],
     #                                     response = response,
     #                                     geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #                                     #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #                                     geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #                                     message = message,
     #                                     gen_name = gen_name,
     #                                     scaling = scaling,
     #                                     centering = centering,
     #                                     omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #                                     AI_cv_nfolds = AI_cv_nfolds,
     #                                     para_tunning = para_tunning,
     #                                     svm_paras_tunning = svm_paras_tunning,
     #                                     svm_type = svm_type,
     #                                     svm_kernel = svm_kernel, # "Gaussian", "Linear","Hyperbolic_tangent", "Polynomial"
     #                                     sigma_value  = sigma_value,       # Default sigma value for RBF kernel
     #                                     C_value  = C_value,             # Default cost parameter
     #                                     degree_value = degree_value,        # Default degree for polynomial kernel
     #                                     scale_value  = scale_value,         # Default scale for polynomial kernel
     #                                     offset_value = offset_value,
     #                                     gamma_value = gamma_value,
     #                                     CI_width_thresholds = CI_width_thresholds,
     #                                     high_reliability_thres = high_reliability_thres,
     #                                     low_reliability_thres = low_reliability_thres,
     #                                     n_components = n_components,
     #                                     threshold = threshold,
     #                                     target = "test_set",
     #                                     #iqr_multiplier = iqr_multiplier,
     #                                     interval_width_low_threshold = interval_width_low_threshold,
     #                                     interval_width_high_threshold = interval_width_high_threshold,
     #                                     interval_width_moderate_threshold = interval_width_moderate_threshold,
     #                                     n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "K-NearestNeighbors" = {
     #          res_model_output <- AI_knn(pheno_object = ml_dat_res[["pheno_clean_data"]],
     #                                     response = response,
     #                                     geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #                                     #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #                                     geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #                                     message = message,
     #                                     gen_name = gen_name,
     #                                     scaling = scaling,
     #                                     centering = centering,
     #                                     omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #                                     AI_cv_nfolds = AI_cv_nfolds,
     #                                     para_tunning = para_tunning,
     #                                     knn_paras_tunning = knn_paras_tunning,
     #                                     k = k,
     #                                     CI_width_thresholds = CI_width_thresholds,
     #                                     high_reliability_thres = high_reliability_thres,
     #                                     low_reliability_thres = low_reliability_thres,
     #                                     n_components = n_components,
     #                                     threshold = threshold,
     #                                     target = "test_set",
     #                                     #iqr_multiplier = iqr_multiplier,
     #                                     interval_width_low_threshold = interval_width_low_threshold,
     #                                     interval_width_high_threshold = interval_width_high_threshold,
     #                                     interval_width_moderate_threshold = interval_width_moderate_threshold,
     #                                     n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "Lasso" = {
     #          res_model_output <- AI_RidgeRegression_Lasso(
     #            pheno_object = ml_dat_res[["pheno_clean_data"]],
     #            response = response,
     #            geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #            #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #            geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #            gen_name = gen_name,
     #            para_tunning = para_tunning,
     #            AI_cv_nfolds = AI_cv_nfolds,
     #            lasso_paras_tunning = lasso_paras_tunning,
     #            message = message,
     #            scaling = scaling,
     #            centering = centering,
     #            omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #            GS_model = GS_model,
     #            lambda_rr = lambda_rr,
     #            CI_width_thresholds = CI_width_thresholds,
     #            high_reliability_thres = high_reliability_thres,
     #            low_reliability_thres = low_reliability_thres,
     #            n_components = n_components,
     #            threshold = threshold,
     #            target = "test_set",
     #            #iqr_multiplier = iqr_multiplier,
     #            interval_width_low_threshold = interval_width_low_threshold,
     #            interval_width_high_threshold = interval_width_high_threshold,
     #            interval_width_moderate_threshold = interval_width_moderate_threshold,
     #            n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "Ridge_Regression" = {
     #          res_model_output <- AI_RidgeRegression_Lasso(
     #            pheno_object = ml_dat_res[["pheno_clean_data"]],
     #            response = response,
     #            geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #            #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #            geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #            gen_name = gen_name,
     #            para_tunning = para_tunning,
     #            AI_cv_nfolds = AI_cv_nfolds,
     #            lasso_paras_tunning = rr_paras_tunning,
     #            message = message,
     #            scaling = scaling,
     #            centering = centering,
     #            omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #            GS_model = GS_model,
     #            lambda_rr = lambda_rr,
     #            CI_width_thresholds = CI_width_thresholds,
     #            high_reliability_thres = high_reliability_thres,
     #            low_reliability_thres = low_reliability_thres,
     #            n_components = n_components,
     #            threshold = threshold,
     #            target = "test_set",
     #            #iqr_multiplier = iqr_multiplier,
     #            interval_width_low_threshold = interval_width_low_threshold,
     #            interval_width_high_threshold = interval_width_high_threshold,
     #            interval_width_moderate_threshold = interval_width_moderate_threshold,
     #            n_bootstrap = n_bootstrap
     #          )
     #        },
     #        "deep_learning_model" = {
     #
     #          res_model_output <- deep_learning_model(
     #            pheno_object=ml_dat_res[["pheno_clean_data"]],
     #            geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
     #            #geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
     #            geno_omic_test_object = ml_dat_res[["merged_data_test"]],
     #            response=response,
     #            gen_name=gen_name,
     #            num_hidden_layers = num_hidden_layers,
     #            neurons_per_layer = neurons_per_layer,
     #            learning_rate = learning_rate,
     #            epochs = epochs,
     #            batch_size = batch_size,
     #            validation_split = validation_split,
     #            early_stop = early_stop,
     #            message = message,
     #            scaling = scaling,
     #            centering = centering,
     #            omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
     #            para_tunning = para_tunning,
     #            param_grid = dpl_paras_tunning
     #          )
     #
     #        },
     #        {
     #          stop(paste(msg, "Select method to calculate geno_cleanmic relationship matrix"), call. = FALSE)
     #        }
     # )

     # res_model_output <- tryCatch({

       switch(GS_model,
              "Xgboost" = {
                tryCatch({
                  res_model_output <- AI_Xgb(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    message = message,
                    gen_name = gen_name,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    AI_cv_nfolds = AI_cv_nfolds,
                    para_tunning = para_tunning,
                    xgb_paras_tunning = xgb_paras_tunning,
                    resample_method_tune = resample_method_tune,
                    number_of_fold_tune = number_of_fold_tune,
                    learning_rate = learning_rate,
                    xgb_gamma = xgb_gamma,
                    xgb_lambda = xgb_lambda,
                    xgb_alpha = xgb_alpha,
                    max_depth = max_depth,
                    subsample = subsample,
                    xgb_booster = xgb_booster,
                    colsample_bytree = colsample_bytree,
                    alpha = alpha,
                    lambda = lambda,
                    iteration = iteration,
                    xgb_rate_drop = xgb_rate_drop,
                    xgb_skip_drop = xgb_skip_drop,
                    xgb_objective = xgb_objective,
                    xgb_sample_type = xgb_sample_type,
                    xgb_normalize_type = xgb_normalize_type,
                    early_stop_for_iteration_xgb = early_stop_for_iteration_xgb,
                    N_feature_impo = N_feature_impo,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in Xgboost model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "RandomForest" = {
                tryCatch({
                  res_model_output <- AI_randomForest(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    message = message,
                    gen_name = gen_name,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    AI_cv_nfolds = AI_cv_nfolds,
                    para_tunning = para_tunning,
                    rf_paras_tunning = rf_paras_tunning,
                    ntree = ntree,
                    mtry = mtry,
                    maxnodes = maxnodes,
                    importance = importance,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in RandomForest model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "PartialLeastSquare" = {
                tryCatch({
                  res_model_output <- AI_pls(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    message = message,
                    gen_name = gen_name,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    para_tunning = para_tunning,
                    ncomp = ncomp,
                    pls_paras_tunning = pls_paras_tunning,
                    resample_method_tune = resample_method_tune,
                    N_feature_impo = N_feature_impo,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in PartialLeastSquare model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "SupportVectorMachine" = {
                tryCatch({
                  res_model_output <- AI_svm(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    message = message,
                    gen_name = gen_name,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    AI_cv_nfolds = AI_cv_nfolds,
                    para_tunning = para_tunning,
                    svm_paras_tunning = svm_paras_tunning,
                    svm_type = svm_type,
                    svm_kernel = svm_kernel,
                    sigma_value = sigma_value,
                    C_value = C_value,
                    degree_value = degree_value,
                    scale_value = scale_value,
                    offset_value = offset_value,
                    gamma_value = gamma_value,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in SupportVectorMachine model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "K-NearestNeighbors" = {
                tryCatch({
                  res_model_output <- AI_knn(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    message = message,
                    gen_name = gen_name,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    AI_cv_nfolds = AI_cv_nfolds,
                    para_tunning = para_tunning,
                    knn_paras_tunning = knn_paras_tunning,
                    k = k,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in K-NearestNeighbors model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "Lasso" = {
                tryCatch({
                  res_model_output <- AI_RidgeRegression_Lasso(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    gen_name = gen_name,
                    para_tunning = para_tunning,
                    AI_cv_nfolds = AI_cv_nfolds,
                    lasso_paras_tunning = lasso_paras_tunning,
                    message = message,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    GS_model = GS_model,
                    lambda_rr = lambda_rr,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in Lasso model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "Ridge_Regression" = {
                tryCatch({
                  res_model_output <- AI_RidgeRegression_Lasso(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    response = response,
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    gen_name = gen_name,
                    para_tunning = para_tunning,
                    AI_cv_nfolds = AI_cv_nfolds,
                    lasso_paras_tunning = rr_paras_tunning,
                    message = message,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    GS_model = GS_model,
                    lambda_rr = lambda_rr,
                    CI_width_thresholds = CI_width_thresholds,
                    high_reliability_thres = high_reliability_thres,
                    low_reliability_thres = low_reliability_thres,
                    n_components = n_components,
                    threshold = threshold,
                    target = "test_set",
                    interval_width_low_threshold = interval_width_low_threshold,
                    interval_width_high_threshold = interval_width_high_threshold,
                    interval_width_moderate_threshold = interval_width_moderate_threshold,
                    n_bootstrap = n_bootstrap
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in Ridge_Regression model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              "deep_learning_model" = {
                tryCatch({
                  res_model_output <- deep_learning_model(
                    pheno_object = ml_dat_res[["pheno_clean_data"]],
                    geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                    geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                    response = response,
                    gen_name = gen_name,
                    num_hidden_layers = num_hidden_layers,
                    neurons_per_layer = neurons_per_layer,
                    learning_rate = learning_rate,
                    epochs = epochs,
                    batch_size = batch_size,
                    validation_split = validation_split,
                    early_stop = early_stop,
                    message = message,
                    scaling = scaling,
                    centering = centering,
                    omic_count = ml_dat_res[["omic_count"]],
                    para_tunning = para_tunning,
                    param_grid = dpl_paras_tunning
                  )
                  res_model_output
                }, error = function(e) {
                  cat("Error in deep_learning_model:", conditionMessage(e), "\n")
                  return(NULL)
                })
              },
              {
                stop("Select method to calculate geno_omic relationship matrix", call. = FALSE)
              }
       )
     # }, error = function(e) {
     #   # Handle the error, you can print a message or take other actions
     #   cat("Error:", conditionMessage(e), "\n")
     #   return(NULL)  # Return NULL or an appropriate value to indicate the failure
     # })

    ######
     # Check the result of the model fitting process
     if (is.null(res_model_output)) {
       cat(paste(GS_model, "model fitting for failed due to an error.\n"))
       res_summary_stat <- NULL
     } else {
       AI_summary_stat_process <- tryCatch({
         res_summary_stat <- summary_statistics_AI(predicted_object = res_model_output[["predicted_values"]],
                                                   pheno_object = ml_dat_res[["pheno_clean_data"]],
                                                   response = response,
                                                   test_set = ml_dat_res[["test_set"]],
                                                   geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
                                                   eval_metrics = eval_metrics,
                                                   model_parameters = res_model_output[["model_parameters"]],
                                                   GS_model = GS_model)

         res_summary_stat  # Return the model object if everything is successful
       }, error = function(e) {
         # Handle the error, you can print a message or take other actions
         cat("Error:", conditionMessage(e), "\n")
         return(NULL)  # Return NULL or an appropriate value to indicate the failure
       })

       if(!is.null(AI_summary_stat_process)){
         res_summary_stat <- AI_summary_stat_process
       }else{
         cat("Error processing the output of machine learning model result.\n")
         res_summary_stat <- NULL
       }
     }

     #View(res_model_output[["predicted_values"]])

     # output <- list(GS_model = GS_model,
     #                res_model_output = res_model_output,
     #                res_summary_stat = res_summary_stat,
     #                geno_qc_stat = geno_qc_stat
     # )

   }
   ### End machine learning
   list(GS_model = GS_model,res_model_output = res_model_output,
        res_summary_stat = res_summary_stat, geno_qc_stat = geno_qc_stat)
   #list(output =  output)
 }, future.seed = TRUE)



 if(!is.null(best_models)){
   names(results) <- best_models[["trait"]]

 } else {
   names(results) <- response
 }


 # Initialize a list to hold the results if returning as a list when system_database is TRUE
 all_results <- list()
 mainDirt <- getwd()

 for (res in seq_along(results)) {
   tryCatch({
     processed_result <- results_handling(
       GS_model = if("GS_model" %in% names(results[[res]])) results[[res]][["GS_model"]] else NULL,
       res_model_output = if("res_model_output" %in% names(results[[res]])) results[[res]][["res_model_output"]] else NULL,
       res_summary_stat = if("res_summary_stat" %in% names(results[[res]])) results[[res]][["res_summary_stat"]] else NULL,
       res_plot = best_models_ggplot_rep,
       res_plot_mean = best_models_ggplot_mean,
       res_plot_result_diagnostic = cv_results_predicted_vs_observed$predicted_vs_observed_plots[[res]],
       test_diagonistic_plots = if(!is.null(results[[res]]$res_model_output)) results[[res]]$res_model_output$diagnostic_plots else NULL,
       res_mod_results_cv_per_trait_model = cv_results_predicted_vs_observed$mod_res_per_trait_per_model,
       geno_qc_stat = geno_qc_stat,
       res_plot_result_diagnostic_cv_only = NULL,
       cv_results_processed = cv_results_processed,
       system_database = system_database,
       plot_filename = if(!is.null(names(results)[res])) names(results)[res] else paste("trait", res, sep = "_"),
       plot_extension = plot_extension,
       plot_width = plot_width,
       plot_height = plot_height,
       plot_units = plot_units,
       plot_dpi = plot_dpi
     )

     # If returning as a list, append the processed result to the all_results list
     if (isTRUE(system_database)) {
       all_results[[length(all_results) + 1]] <- processed_result
       names(all_results)[length(all_results)] <- if(!is.null(names(results)[res])) names(results)[res] else paste("trait", res, sep = "_")
     }
   }, error = function(e) {
     setwd(mainDirt)
     message(paste("Error processing result", res, ":", e$message))
   })
 }

 # for (res in seq_along(results)) {
 #
 #   #names(results[[1]])
 #
 #
 #   processed_result <- results_handling(GS_model = if("GS_model" %in% names(results[[res]])) results[[res]][["GS_model"]] else NULL,
 #                           res_model_output = if("res_model_output" %in% names(results[[res]])) results[[res]][["res_model_output"]] else NULL,
 #                           res_summary_stat = if("res_summary_stat" %in% names(results[[res]])) results[[res]][["res_summary_stat"]] else NULL,
 #                           res_plot = best_models_ggplot_rep,
 #                           res_plot_mean = best_models_ggplot_mean,
 #                           res_plot_result_diagnostic = cv_results_predicted_vs_observed$predicted_vs_observed_plots[[res]],
 #                           test_diagonistic_plots = results[[res]]$res_model_output$diagnostic_plots,
 #                           res_mod_results_cv_per_trait_model = cv_results_predicted_vs_observed$mod_res_per_trait_per_model,
 #                           geno_qc_stat = geno_qc_stat,
 #                           res_plot_result_diagnostic_cv_only = NULL,
 #                           cv_results_processed = cv_results_processed,
 #                           system_database = system_database,
 #                           plot_filename = if(!is.null(names(results)[res])) names(results)[res] else paste("trait", res, sep = "_"),
 #                           #Plot_name_result_diagnostic = if(!is.null(names(res_mod_results_cv_per_trait_model)[res])) names(results)[res] else paste("trait_diganostic", res, sep = "_"),
 #                           plot_extension = plot_extension,
 #                           plot_width = plot_width,
 #                           plot_height = plot_height,
 #                           plot_units = plot_units,
 #                           plot_dpi = plot_dpi)
 #
 #   # If returning as a list, append the processed result to the all_results list
 #   if(isTRUE(system_database)) {
 #     all_results[[length(all_results) + 1]] <- processed_result
 #     names(all_results)[res] <- if(!is.null(names(results)[res])) names(results)[res] else paste("trait", res, sep = "_")
 #   }
 #   # Otherwise, the results_handling function handles file creation and saving
 # }

 # Return the list of all results if system_database is TRUE
 if(isTRUE(system_database)) {
   return(all_results)
 } else {
   # system_database is TRUE is false,
   ## implies If results are not being returned as a list,
   # just return a success message or NULL for job done
   return(invisible(TRUE))
 }
 # sik end###############################################################################################
#if(is.null(best_models)){
 #### These models only works with one environment/location
#  if(length(pheno_clean[["pheno_clean_data"]][,gen_name])==length(unique(pheno_clean[["pheno_clean_data"]][,gen_name]))){
#
#      if ((GS_model %in% bayes_valid_models && is.null(rand_term_model_bayesian)) ||
#          (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models)) ||
#          (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models))) {
#
#        res_model_output <- bayes_finalize_A_B_C_BL_BRR(fixed = fixed,
#                                                   random = random,
#                                                   GS_model = GS_model,
#                                                   response = response,
#                                                   weights = weights,
#                                                   fixed_term_model_bayesian = fixed_term_model_bayesian,
#                                                   rand_term_model_bayesian = rand_term_model_bayesian,
#                                                   pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                   geno_data = if("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
#                                                   omic1_data = if("omic1_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic1_model_ready"]] else NULL,
#                                                   omic2_data = if("omic2_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic2_model_ready"]] else NULL,
#                                                   omic3_data = if("omic3_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["omic3_model_ready"]] else NULL,
#                                                   gen_name = gen_name,
#                                                   nIter = nIter,
#                                                   burnIn = burnIn,
#                                                   thin = thin,
#                                                   omics_data_label = omics_data_label,
#                                                   scaling = scaling)
#
#      # Compute summary statistics and plot accuracy
#      res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
#      #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
#      res_model_output <- res_model_output[["bayes_result"]]
#
#      # return(results_handling(GS_model = GS_model,
#      #                         res_model_output = res_model_output,
#      #                         res_summary_stat = res_summary_stat,
#      #                         res_plot = NULL,
#      #                         geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
#      #                         system_database = system_database))
#
#  }
# } ## End of  Bayes A, B, C, BRR, BL
#
#  ##########################################################################
#  #########################################################################
#  ## Start of Reproducing Kernel Hilbert Spaces Regression RKHS,         ##
#  ## (BRR- Bayesian GBLUP ) and GBLUP (asreml) Model                     ##
#  ## for Single Location and multiple loc                                ##
#  ##                                                                     ##
#  ##                                                                     ##
#  ##########################################################################
#  #######################################################################
#
#  ## NOTE
#  ## BRR is chaneg to G-BRR
#  ## This is to make distinction between BRR for marker matrix and GBLUP
#  # Check conditions for GS_model and rand_term_model_bayesian
#  if ((GS_model %in% c("RKHS", "GBLUP_BRR", "GBLUP") && is.null(rand_term_model_bayesian)) ||
#      (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models)) ||
#      (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_gblup_valid_models))) {
#
#      # Rename GS_model for GBLUP_BRR case
#      if (GS_model == "GBLUP_BRR") {
#          GS_modeluse <- GS_model
#          GS_model <- "BRR"
#      }
#
#      if (GS_model %in% c("BRR", "RKHS")) {
#          # Run Bayesian model for BRR and RKHS
#        res_model_output <- bayes_finalize_RKHS_GBLUPBRR(fixed = fixed,
#                                                    random = random,
#                                                    GS_model = GS_model,
#                                                    response = response,
#                                                    weights = weights,
#                                                    fixed_term_model_bayesian = fixed_term_model_bayesian,
#                                                    rand_term_model_bayesian = rand_term_model_bayesian,
#                                                    pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                    gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#                                                    omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#                                                    omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#                                                    omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#                                                    gen_name = gen_name,
#                                                    nIter = nIter,
#                                                    burnIn = burnIn,
#                                                    thin = thin,
#                                                    heter_groups = heter_groups,
#                                                    omics_kernel_label = omics_kernel_label)
#
#          # Compute summary statistics and plot accuracy
#          res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
#          #res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
#          res_model_output <- res_model_output[["bayes_result"]]
#          ### This part is for GBLUP_BRR
#          if(exists("GS_modeluse")){
#            GS_model <-  GS_modeluse
#          }
#
#          # return(results_handling(GS_model = GS_model,
#          #                         res_model_output = res_model_output,
#          #                         res_summary_stat = res_summary_stat,
#          #                         res_plot = NULL,
#          #                         geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
#          #                         system_database = system_database))
#
#      } else {
#        if (GS_model == "GBLUP" && engine == 'asreml') {
#        # gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL
#        # if(!matrixcalc::is.positive.definite(gmatrix)) stop("GRM issue")
#
#        #  # Run GBLUP model with ASReml
#          mod <- asreml_utilis_new(fixed = fixed,
#                               random = random,
#                               cova = cova,
#                               GS_model = GS_model,
#                               response = response,
#                               pheno_data = pheno_clean[["pheno_clean_data"]],
#                               gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#                               omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#                               omic2_kernel = if("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#                               omic3_kernel = if("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#                               gen_name = gen_name,
#                               heter_groups = heter_groups,
#                               heter_resid = heter_resid,
#                               var_cov_str = var_cov_str,
#                               weights = weights,
#                               pworkspace = pworkspace,
#                               workspace = workspace,
#                               maxit = maxit,
#                               inverse = inverse,
#                               epsilon = epsilon,
#                               engine = engine)
#
#          # Extract model output for ASReml
#
#          res_model_output <- asreml_mod_output_new(
#              mod_asreml = mod,
#              pheno_data = pheno_clean[["pheno_clean_data"]],
#              gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
#              omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
#              omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
#              omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
#              heter_groups = heter_groups,
#              gen_name = gen_name,
#              response = response,
#              var_cov_str = var_cov_str,
#              heter_resid = heter_resid,
#              pworkspace = pworkspace,
#              workspace = workspace,
#              maxit = maxit
#          )
#
#
#          res_summary_stat <- summary_statistics_asreml(mod =  res_model_output[["Asreml_model"]],
#                                                        response = response,
#                                                        pheno_data = pheno_clean[["pheno_clean_data"]],
#                                                        heter_groups = heter_groups,
#                                                        predicted_value =  res_model_output[["Predicted_value"]],
#                                                        pred_heter_groups = NULL,
#                                                        variance_components = res_model_output[["Variance_components"]],
#                                                        eval_metrics = eval_metrics)
#
#
#          # remove_from_global <- function(var_names) {
#          #   for (var_name in var_names) {
#          #     if(exists(var_name, envir = .GlobalEnv)) {
#          #       rm(list = var_name, envir = .GlobalEnv)
#          #       #print(paste("Object", var_name, "removed from global environment."))
#          #     } else {
#          #       #print(paste("Object", var_name, "not found in global environment."))
#          #     }
#          #   }
#          # }
#          #
#          # #rm inv_object
#          # my_variable <- c("G_inv", "omic1_inv", "omic2_inv", "omic3_inv")
#          # remove_from_global(my_variable)
#
#
#          # return(results_handling(GS_model = GS_model,
#          #                         res_model_output = res_model_output,
#          #                         res_summary_stat = res_summary_stat,
#          #                         res_plot =  NULL,
#          #                         geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
#          #                         system_database = system_database))
#
#        }
#      }
#  }
#  #### END GBLUP_RKHS, GBLUP_BRR and GBLUP (asreml)
#
#  ######################################################
#  ######################################################
#  ##                                                 ###
#  ## Machine Learning Models                         ###
#  ##                                                 ###
#  ######################################################
#  ######################################################
#
#  ###  Start ML Analysis
#  ###
#
#
#      # AI_valid_models <- c("Xgboost",
#      #                      "RandomForest",
#      #                      "PartialLeastSquare",
#      #                      "SupportVectorMachine",
#      #                      "K-NearestNeighbors",
#      #                      "Lasso",
#      #                      "Ridge_Regression",
#      #                      "deep_learning_model")
#
#      if (GS_model %in% AI_valid_models) {
#          if (length(unique(pheno_clean[["pheno_clean_data"]][, gen_name])) > length(pheno_clean[["pheno_clean_data"]][, gen_name])) {
#              stop(paste(msg, GS_model, 'only works for single location/enviroment.'), call. = FALSE)
#          }
#
#          # ml_dat_res <- ML_data_processing(pheno_clean = pheno_clean,
#          #                                  response = response,
#          #                                  gen_name = gen_name,
#          #                                  geno_clean = if ("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
#          #                                  omic_clean = if (!"geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL
#          #                                  )
#
#          switch(GS_model,
#                 "Xgboost" = {
#                     res_model_output <- AI_Xgb(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                response = response,
#                                                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                message = message,
#                                                gen_name = gen_name,
#                                                scaling = scaling,
#                                                centering = centering,
#                                                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                AI_cv_nfolds = AI_cv_nfolds,
#                                                para_tunning = para_tunning,
#                                                xgb_paras_tunning = xgb_paras_tunning
#                     )
#                 },
#                 "RandomForest" = {
#                     res_model_output <- AI_randomForest(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                         response = response,
#                                                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                         message = message,
#                                                         gen_name = gen_name,
#                                                         scaling = scaling,
#                                                         centering = centering,
#                                                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                         AI_cv_nfolds = AI_cv_nfolds,
#                                                         para_tunning = para_tunning,
#                                                         rf_paras_tunning = rf_paras_tunning
#                     )
#                 },
#                 "PartialLeastSquare" = {
#                     res_model_output <-  AI_pls(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                 response = response,
#                                                 geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                 geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                 message = message,
#                                                 gen_name = gen_name,
#                                                 scaling = scaling,
#                                                 centering = centering,
#                                                 omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                 para_tunning = para_tunning,
#                                                 pls_paras_tunning = pls_paras_tunning)
#                 },
#                 "SupportVectorMachine" = {
#                     res_model_output <- AI_svm(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                response = response,
#                                                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                message = message,
#                                                gen_name = gen_name,
#                                                scaling = scaling,
#                                                centering = centering,
#                                                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                AI_cv_nfolds = AI_cv_nfolds,
#                                                para_tunning = para_tunning,
#                                                svm_paras_tunning = svm_paras_tunning
#                     )
#                 },
#                 "K-NearestNeighbors" = {
#                     res_model_output <- AI_knn(pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                response = response,
#                                                geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                                                geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                                                message = message,
#                                                gen_name = gen_name,
#                                                scaling = scaling,
#                                                centering = centering,
#                                                omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                                                AI_cv_nfolds = AI_cv_nfolds,
#                                                para_tunning = para_tunning,
#                                                knn_paras_tunning = knn_paras_tunning
#                     )
#                 },
#                 "Lasso" = {
#                     res_model_output <- AI_RidgeRegression_Lasso(
#                         pheno_object = ml_dat_res[["pheno_clean_data"]],
#                         response = response,
#                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                         gen_name = gen_name,
#                         para_tunning = para_tunning,
#                         AI_cv_nfolds = AI_cv_nfolds,
#                         lasso_paras_tunning = lasso_paras_tunning,
#                         message = message,
#                         scaling = scaling,
#                         centering = centering,
#                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                         GS_model = GS_model
#                     )
#                 },
#                 "Ridge_Regression" = {
#                     res_model_output <- AI_RidgeRegression_Lasso(
#                         pheno_object = ml_dat_res[["pheno_clean_data"]],
#                         response = response,
#                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                         gen_name = gen_name,
#                         para_tunning = para_tunning,
#                         AI_cv_nfolds = AI_cv_nfolds,
#                         lasso_paras_tunning = rr_paras_tunning,
#                         message = message,
#                         scaling = scaling,
#                         centering = centering,
#                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                         GS_model = GS_model,
#                     )
#                 },
#                 "deep_learning_model" = {
#
#                     res_model_output <- deep_learning_model(
#                         pheno_object=ml_dat_res[["pheno_clean_data"]],
#                         geno_omic_object = ml_dat_res[["merged_data"]][["merge_data"]],
#                         geno_omic_test_object = ml_dat_res[["merged_data_test"]][["merge_data"]],
#                         response=response,
#                         gen_name=gen_name,
#                         num_hidden_layers = num_hidden_layers,
#                         neurons_per_layer = neurons_per_layer,
#                         learning_rate = learning_rate,
#                         epochs = epochs,
#                         batch_size = batch_size,
#                         validation_split = validation_split,
#                         early_stop = early_stop,
#                         message = message,
#                         scaling = scaling,
#                         centering = centering,
#                         omic_count = if("omic_count"%in%names(ml_dat_res)) ml_dat_res[["omic_count"]] else NULL,
#                         para_tunning = para_tunning,
#                         param_grid = dpl_paras_tunning
#                     )
#
#                 },
#                 {
#                     stop(paste(msg, "Select method to calculate geno_cleanmic relationship matrix"), call. = FALSE)
#                 }
#          )
#
# ##browser()
# #View(res_model_output[["predicted_values"]])
#          res_summary_stat <- summary_statistics_AI(predicted_object = res_model_output[["predicted_values"]],
#                                                    pheno_object = ml_dat_res[["pheno_clean_data"]],
#                                                    response = response,
#                                                    test_set = ml_dat_res[["test_set"]],
#                                                    geno_omic_object = ml_dat_res[["merged_data"]],
#                                                    eval_metrics = eval_metrics,
#                                                    model_parameters = res_model_output[["model_parameters"]],
#                                                    GS_model = GS_model
#          )
#
#          # res_plot <- plot_acc_AI(mod = res_model_output,
#          #                         pheno_object = ml_dat_res[["pheno_clean_data"]],
#          #                         response = response,
#          #                         test_set = ml_dat_res[["test_set"]],
#          #                         GS_model = GS_model
#          # )
#
#          # return(results_handling(GS_model = GS_model,
#          #                         res_model_output = res_model_output,
#          #                         res_summary_stat = res_summary_stat,
#          #                         res_plot = NULL,
#          #                         geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
#          #                         system_database = system_database))
#
#      }
#
#   ### End machine learning

} ## end of function
