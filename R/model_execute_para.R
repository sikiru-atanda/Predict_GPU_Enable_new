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
    coefficient_1 = NULL,
    coefficient_2 = NULL,
    coefficient_3 = NULL,
    coefficient_4 = NULL,
    train_set = NULL,
    test_set = NULL,
    gmatrix_method = NULL,
    kernel_method = NULL,
    # kernel_method = c("Gaussian_kernel",
    #                   "Linear_kernel",
    #                   "Poly2_kernel",
    #                   "Poly3_kernel",
    #                   "Poly4_kernel",
    #                   "spectral_kernel",
    #                   "Matern_kernel",
    #                   "Composite_kernel",
    #                   "Normalized_laplacian_kernel",
    #                   "Anova_radial_basis_kernel"
    # ),
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
    GS_model = c("GBLUP",
                 "GBLUP_BRR",
                 "RKHS",
                 "BRR",
                 "BayesA",
                 "BayesB",
                 "BayesC",
                 "BL",
                 "Xgboost",
                 "RandomForest",
                 "PartialLeastSquare",
                 "SupportVectorMachine",
                 "K-NearestNeighbors",
                 "Lasso",
                 "Ridge_Regression"
                 ),
    eval_metrics = c("Accuracy",
                     "Mean_Squared_Error",
                     "Bias",
                     "Root_Mean_Squared_Error",
                     "Relative_Squared_Error",
                     "Mean_Absolute_Error",
                     "Mean_Absolute_Percent_Error"),
    para_tunning = FALSE,
    fixed_term_model_bayesian = 'FIXED',
    rand_term_model_bayesian = NULL,
    core = NULL,
    engine = NULL,
    scale = TRUE,
    #scaled = TRUE,
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
    xgb_paras_tunning = NULL,
    rf_paras_tunning = NULL,
    pls_paras_tunning = NULL,
    svm_paras_tunning = NULL,
    knn_paras_tunning = NULL,
    lasso_paras_tunning = NULL,
    rr_paras_tunning = NULL,
    dpl_paras_tunning = NULL,
    AI_cv_nfolds = 5,
    ...
) {

#browser()
    msg <- sprintf("==================================================\n")

    # Define available models and variance structures
    var_cov_str_available <- c("us",
                               "corgh",
                               "corgv",
                               "corh",
                               "corv",
                               "fa",
                               "rr")
    AI_valid_models <- c("Xgboost", "RandomForest", "PartialLeastSquare", "SupportVectorMachine", "K-NearestNeighbors", "Lasso", "Ridge_Regression", "deep_learning_model")
    bayes_valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL")
    bayes_gblup_valid_models <- c("GBLUP_BRR", "RKHS")

    # Check for mandatory phenotypic data
    if (is.null(pheno_data)){
      stop(paste(msg, "Phenotypic data is missing."), call. = FALSE)
    }

    # Check for ASReml requirement for GBLUP
    if (GS_model == "GBLUP" && engine != "asreml") {
      stop(paste(msg, "ASReml software is required to fit GBLUP for single or multi-environment."), call. = FALSE)
    }

    # Check for multi-environment structure and required inputs for GBLUP
    if (length(pheno_data[,gen_name]) > length(unique(pheno_data[,gen_name])) & is.null(heter_groups)){
      stop(paste(msg, "Your phenotypic data has a multi-environment structure, but the column containing the environment/location is missing. Provide it in heter_groups."), call. = FALSE)

      # Check for required inputs for multi-environment GBLUP
      # if ((GS_model %in% c("GBLUP_BRR", "RKHS")) &&
      #     ((is.null(gmatrix) || is.null(gkernel) || is.null(omic1_kernel) || is.null(omic2_kernel) || is.null(omic3_kernel)) ||
      #      ((is.null(geno_data) || is.null(omic1_data) || is.null(omic2_data) || is.null(omic3_data)) &&
      #       (is.null(gmatrix_method) || is.null(kernel_method))))) {
      #   stop(paste(msg, "To fit a Bayesian multi-environment GBLUP model, you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel), or an omics kernel. Additionally, you can provide genomic or omics data. Ensure you provide instructions on the genomic relationship matrix method or kernel method to calculate the relationship matrix."), call. = FALSE)
      # }

    }

    # Check for single environment GBLUP and Bayesian models
    # if (length(pheno_data[,gen_name]) == length(unique(pheno_data[,gen_name]))){
    #
    #   # Check for required inputs for Bayesian or ASReml single environment GBLUP models
    #   if ((GS_model %in% c(bayes_gblup_valid_models, "GBLUP")) &&
    #       ((is.null(gmatrix) || is.null(gkernel) || is.null(omic1_kernel) || is.null(omic2_kernel) || is.null(omic3_kernel)) ||
    #        (is.null(geno_data) || is.null(omic1_data) || is.null(omic2_data) || is.null(omic3_data)) &&
    #        (is.null(gmatrix_method) || is.null(kernel_method)))) {
    #     stop(paste(msg, "To fit a Bayesian or ASReml single environment GBLUP model, you need either a genomic matrix (gmatrix) or a genomic kernel (gkernel), or an omics kernel. Additionally, you can provide genomic or omics data. Ensure you provide instructions on the genomic relationship matrix method or kernel method to calculate the relationship matrix."), call. = FALSE)
    #   }
    #
    #
    #   # Check for required inputs for Bayesian models
    #   if ((GS_model %in% bayes_valid_models) &&
    #       (is.null(geno_data) || is.null(omic1_data) || is.null(omic2_data) || is.null(omic3_data))) {
    #     stop(paste(msg, "To fit a Bayesian model, provide genomic or omics data."), call. = FALSE)
    #   }
    #
    # }

    ##
### Check phenotype_to_model for details
 #    This serve as gateway between phenotype-precheck function and readiness of
 #    the phenotypic data for model fitting.

 pheno_clean <- phenotype_to_model(pheno_data = pheno_data,
                                   pheno_data_train = pheno_data_train,
                                   pheno_data_test = pheno_data_test,
                                   train_set = train_set,
                                   test_set = test_set,
                                   response = response,
                                   gen_name = gen_name)

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
                                   message = message)

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
                                message = message)

 # Process omic2 data
 omic2_res <- process_omic_data(omic_data = omic2_data,
                                train_omic_data = train_omic2_data,
                                test_omic_data = test_omic2_data,
                                kernel_method = kernel_method,
                                pheno_clean_list = pheno_clean,
                                gen_name = gen_name,
                                test_set = test_set,
                                train_set = train_set,
                                message = message)

 # Process omic3 data
 omic3_res <- process_omic_data(omic_data = omic3_data,
                                train_omic_data = train_omic3_data,
                                test_omic_data = test_omic3_data,
                                kernel_method = kernel_method,
                                pheno_clean_list = pheno_clean,
                                gen_name = gen_name,
                                test_set = test_set,
                                train_set = train_set,
                                message = message)
 #################
 geno_omic_model_ready_list <- list()
 gmatrix_kernel_model_ready_list <- list()

 if (length(geno_res)!=0 && all(c("gmatrix", "geno_model_ready") %in% names(geno_res))) {
   gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] <- geno_res[["gmatrix"]]
   geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
 } else if (length(geno_res)!=0 && "geno_model_ready" %in% names(geno_res)) {
   geno_omic_model_ready_list[["geno_model_ready"]] <- geno_res[["geno_model_ready"]]
 } else if(!is.null(gmatrix)){
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
         }
     }
 }


 }
### Concatenation of omics for ML

 if (GS_model %in% AI_valid_models) {

   ml_dat_res <- ML_data_processing(pheno_clean = pheno_clean,
                                    response = response,
                                    gen_name = gen_name,
                                    geno_clean = if ("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
                                    omic_clean = if (!"geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL
   )


 }

 ### Ends
 ##########################################################################
 #########################################################################
 ## Start of Bayes A, B, C, BL and BRR Models for Single Location       ##
 ##  This only accommodate n x p matrix  not nxn                        ##
 ##                                                                     ##
 ##########################################################################
 #######################################################################


    ### Model BRR for single location
 #if(((GS_model=="BRR") & is.null(rand_term_model_bayesian)) || ((is.null(GS_model) & (rand_term_model_bayesian=="BRR")))){

    # if((isTRUE(GS_model== "BRR" | GS_model== "BayesA"|  GS_model== "BayesB"| GS_model== "BayesC" | GS_model== "BL") & is.null(rand_term_model_bayesian)) |
    #    ((is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL"))))) |
    #    ((!is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL")))))){

 #### These models only works with one environment/location
 if(length(pheno_clean[["pheno_clean_data"]][,gen_name])==length(unique(pheno_clean[["pheno_clean_data"]][,gen_name]))){

     if ((GS_model %in% bayes_valid_models && is.null(rand_term_model_bayesian)) ||
         (is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models)) ||
         (!is.null(GS_model) && any(rand_term_model_bayesian %in% bayes_valid_models))) {


       # if(is.null(geno_omic_model_ready_list[["geno_model_ready"]])){
       #
       #   stop("geno is null")
       # }
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
                                                  omics_data_label = omics_data_label
         )

     # Compute summary statistics and plot accuracy
     res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
     res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
     res_model_output <- res_model_output[["bayes_result"]]

     return(results_handling(GS_model = GS_model,
                             res_model_output = res_model_output,
                             res_summary_stat = res_summary_stat,
                             res_plot = res_plot,
                             geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
                             system_database = system_database))

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
 ## BRR is chaneg to G-BRR
 ## This is to make distinction between BRR for marker matrix and GBLUP
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
                                                   omics_kernel_label = omics_kernel_label)

         # Compute summary statistics and plot accuracy
         res_summary_stat <- summary_statistics_bayes(mod = res_model_output[["bayes_model"]], eval_metrics = eval_metrics)
         res_plot <- plot_acc(mod = res_model_output[["bayes_model"]], response = response)
         res_model_output <- res_model_output[["bayes_result"]]
         ### This part is for GBLUP_BRR
         if(exists("GS_modeluse")){
           GS_model <-  GS_modeluse
         }

         return(results_handling(GS_model = GS_model,
                                 res_model_output = res_model_output,
                                 res_summary_stat = res_summary_stat,
                                 res_plot = res_plot,
                                 geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
                                 system_database = system_database))

     } else if (GS_model == "GBLUP" && engine == 'asreml') {

       # gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL
       # if(!matrixcalc::is.positive.definite(gmatrix)) stop("GRM issue")

       #  # Run GBLUP model with ASReml
         mod <- asreml_utilis(fixed = fixed,
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
                              core = core,
                              pworkspace = pworkspace,
                              workspace = workspace,
                              maxit = maxit,
                              inverse = inverse,
                              epsilon = epsilon,
                              engine = engine)

         # Extract model output for ASReml
         res_model_output <- asreml_mod_output(
             mod_asreml = mod,
             pheno_data = pheno_clean[["pheno_clean_data"]],
             gmatrix = if("gmatrix_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["gmatrix_model_ready"]] else NULL,
             omic1_kernel = if("omic1_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic1_kernel_model_ready"]] else NULL,
             omic2_kernel = if ("omic2_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic2_kernel_model_ready"]] else NULL,
             omic3_kernel = if ("omic3_kernel_model_ready" %in% names(gmatrix_kernel_model_ready_list)) gmatrix_kernel_model_ready_list[["omic3_kernel_model_ready"]] else NULL,
             heter_groups = heter_groups,
             gen_name = gen_name,
             var_cov_str = var_cov_str,
             heter_resid = heter_resid,
             pworkspace = pworkspace,
             workspace = workspace,
             maxit = maxit
         )

         res_summary_stat <- summary_statistics_asreml(mod = res_model_output[["Asreml_model"]],
                                                       response = response,
                                                       heter_groups = heter_groups,
                                                       predicted_value =  res_model_output[["Predicted_value"]],
                                                       pred_heter_groups = NULL,
                                                       variance_components = res_model_output[["Variance_components"]],
                                                       eval_metrics = eval_metrics)


         remove_from_global <- function(var_names) {
           for (var_name in var_names) {
             if(exists(var_name, envir = .GlobalEnv)) {
               rm(list = var_name, envir = .GlobalEnv)
               #print(paste("Object", var_name, "removed from global environment."))
             } else {
               #print(paste("Object", var_name, "not found in global environment."))
             }
           }
         }

         #rm inv_object
         my_variable <- c("G_inv", "omic1_inv", "omic2_inv", "omic3_inv")
         remove_from_global(my_variable)


         return(results_handling(GS_model = GS_model,
                                 res_model_output = res_model_output,
                                 res_summary_stat = res_summary_stat,
                                 res_plot =  NULL,
                                 geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
                                 system_database = system_database))

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

 ###  Start ML Analysis
 ###


     # AI_valid_models <- c("Xgboost",
     #                      "RandomForest",
     #                      "PartialLeastSquare",
     #                      "SupportVectorMachine",
     #                      "K-NearestNeighbors",
     #                      "Lasso",
     #                      "Ridge_Regression",
     #                      "deep_learning_model")

     if (GS_model %in% AI_valid_models) {
         if (length(unique(pheno_clean[["pheno_clean_data"]][, gen_name])) > length(pheno_clean[["pheno_clean_data"]][, gen_name])) {
             stop(paste(msg, GS_model, 'only works for single location/enviroment.'), call. = FALSE)
         }

         # ml_dat_res <- ML_data_processing(pheno_clean = pheno_clean,
         #                                  response = response,
         #                                  gen_name = gen_name,
         #                                  geno_clean = if ("geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL,
         #                                  omic_clean = if (!"geno_model_ready" %in% names(geno_omic_model_ready_list)) geno_omic_model_ready_list[["geno_model_ready"]] else NULL
         #                                  )

         switch(GS_model,
                "Xgboost" = {
                    res_model_output <- AI_Xgb(pheno_object = ml_dat_res[["pheno_clean_data"]],
                                               response = response,
                                               geno_omic_object = ml_dat_res[["merged_data"]],
                                               geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                                               message = message,
                                               gen_name = gen_name,
                                               scale = scale,
                                               AI_cv_nfolds = AI_cv_nfolds,
                                               para_tunning = para_tunning,
                                               xgb_paras_tunning = xgb_paras_tunning
                    )
                },
                "RandomForest" = {
                    res_model_output <- AI_randomForest(pheno_object = ml_dat_res[["pheno_clean_data"]],
                                                        response = response,
                                                        geno_omic_object = ml_dat_res[["merged_data"]],
                                                        geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                                                        message = message,
                                                        gen_name = gen_name,
                                                        scale = scale,
                                                        AI_cv_nfolds = AI_cv_nfolds,
                                                        para_tunning = para_tunning,
                                                        rf_paras_tunning = rf_paras_tunning
                    )
                },
                "PartialLeastSquare" = {
                    res_model_output <-  AI_pls(pheno_object = ml_dat_res[["pheno_clean_data"]],
                                                response = response,
                                                geno_omic_object = ml_dat_res[["merged_data"]],
                                                geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                                                message = message,
                                                gen_name = gen_name,
                                                scale = scale,
                                                para_tunning = para_tunning,
                                                pls_paras_tunning = pls_paras_tunning)
                },
                "SupportVectorMachine" = {
                    res_model_output <- AI_svm(pheno_object = ml_dat_res[["pheno_clean_data"]],
                                               response = response,
                                               geno_omic_object = ml_dat_res[["merged_data"]],
                                               geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                                               message = message,
                                               gen_name = gen_name,
                                               scale = scale,
                                               AI_cv_nfolds = AI_cv_nfolds,
                                               para_tunning = para_tunning,
                                               svm_paras_tunning = svm_paras_tunning
                    )
                },
                "K-NearestNeighbors" = {
                    res_model_output <- AI_knn(pheno_object = ml_dat_res[["pheno_clean_data"]],
                                               response = response,
                                               geno_omic_object = ml_dat_res[["merged_data"]],
                                               geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                                               message = message,
                                               gen_name = gen_name,
                                               scale = scale,
                                               AI_cv_nfolds = AI_cv_nfolds,
                                               para_tunning = para_tunning,
                                               knn_paras_tunning = knn_paras_tunning
                    )
                },
                "Lasso" = {
                    res_model_output <- AI_RidgeRegression_Lasso(
                        pheno_object = ml_dat_res[["pheno_clean_data"]],
                        response = response,
                        geno_omic_object = ml_dat_res[["merged_data"]],
                        geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                        gen_name = gen_name,
                        para_tunning = para_tunning,
                        AI_cv_nfolds = AI_cv_nfolds,
                        lasso_paras_tunning = lasso_paras_tunning,
                        message = message,
                        scale = scale,
                        GS_model = GS_model
                    )
                },
                "Ridge_Regression" = {
                    res_model_output <- AI_RidgeRegression_Lasso(
                        pheno_object = ml_dat_res[["pheno_clean_data"]],
                        response = response,
                        geno_omic_object = ml_dat_res[["merged_data"]],
                        geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                        gen_name = gen_name,
                        para_tunning = para_tunning,
                        AI_cv_nfolds = AI_cv_nfolds,
                        lasso_paras_tunning = rr_paras_tunning,
                        message = message,
                        scale = scale,
                        GS_model = GS_model,
                    )
                },
                "deep_learning_model" = {

                    res_model_output <- deep_learning_model(
                        pheno_object=ml_dat_res[["pheno_clean_data"]],
                        geno_omic_object = ml_dat_res[["merged_data"]],
                        geno_omic_test_object = ml_dat_res[["merged_data_test"]],
                        response=response,
                        gen_name=gen_name,
                        message = message,
                        scale = scale,
                        para_tunning = para_tunning,
                        param_grid = dpl_paras_tunning
                    )

                },
                {
                    stop(paste(msg, "Select method to calculate geno_cleanmic relationship matrix"), call. = FALSE)
                }
         )

#browser()
#View(res_model_output[["predicted_values"]])
         res_summary_stat <- summary_statistics_AI(predicted_object = res_model_output[["predicted_values"]],
                                                   pheno_object = ml_dat_res[["pheno_clean_data"]],
                                                   response = response,
                                                   test_set = ml_dat_res[["test_set"]],
                                                   geno_omic_object = ml_dat_res[["merged_data"]],
                                                   eval_metrics = eval_metrics,
                                                   model_parameters = res_model_output[["model_parameters"]],
                                                   GS_model = GS_model
         )

         res_plot <- plot_acc_AI(mod = res_model_output,
                                 pheno_object = ml_dat_res[["pheno_clean_data"]],
                                 response = response,
                                 test_set = ml_dat_res[["test_set"]],
                                 GS_model = GS_model
         )

         return(results_handling(GS_model = GS_model,
                                 res_model_output = res_model_output,
                                 res_summary_stat = res_summary_stat,
                                 res_plot = NULL,
                                 geno_qc_stat =if("clean_geno_qcstat" %in% names(geno_res)) geno_res[["clean_geno_qcstat"]][["qc_metrics_and_summary_stat"]] else NULL,
                                 system_database = system_database))

     }

  ### End machine learning

} ## end of function
