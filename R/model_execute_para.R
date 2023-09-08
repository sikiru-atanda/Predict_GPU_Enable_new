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
#' @param center if X-variables should be standardized
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
#' @param VarCov_str user defined variance-covariance structure
#' @param engine if user has asreml
#' @param workspace allocate memory for asreml model fit
#' @param pworkspace allocate memory for predict function in asreml
#' @param bending this is important when the relationship matrix is not positive definitive. It fix it for the user. it has be TRUE
#' @param maxit number of iteration for asreml
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
model_execute <- function(
    pheno_data=NULL,
    pheno_data_train = NULL,
    pheno_data_test = NULL,
    geno_data = NULL,
    omic1_data = NULL,
    omic2_data = NULL,
    omic3_data = NULL,
    gmatrix= NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    train_geno_data = NULL,
    train_omic1_data = NULL,
    train_omic2_data = NULL,
    train_omic3_data = NULL,
    test_geno_data = NULL,
    test_omic1_data = NULL,
    test_omic2_data = NULL,
    test_omic3_data = NULL,
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
    VarCov_str = NULL,
    weights =NULL,
    nIter=NULL,
    burnIn=NULL,
    thin=NULL,
    GS_model = c("GBLUP",
                 "RKHS",
                 "BRR",
                 "BayesA",
                 "BayesB",
                 "BayesC",
                 "BayesL",
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
    message = TRUE,
    center = TRUE,
    workspace = 1e08,
    pworkspace= 1e06,
    maxit = 50,
    bending = TRUE,
    ...
) {

    msg <- sprintf("==================================================\n")
    ### Get clean pheno data for model fit

### Check phenotype_to_model for details
 #    This serve as gateway between phenotype-precheck function and readiness of
 #    the phenotypic data for model fitting.

 pheno_clean <- phenotype_to_model(
                    pheno_data = pheno_data,
                    pheno_data_train = pheno_data_train,
                    pheno_data_test = pheno_data_test,
                    #train_set = train_set,
                    #test_set = test_set,
                    response = response,
                    gen_name = gen_name)

## pheno_clean is a list that can have one or two elements
 ## One element if only pheno_data is provided
 ## Two elements if pheno_data/pheno_training and pheno_data_testing was provided as input.
 ##

 ## Check if the pheno_data in the pheno_clean is declared model fit
 if(attr(pheno_clean[[1]], "cleared")!="for_model_fit" && all(class(pheno_clean[[1]])!=c("data.frame", "phenotype"))) {

    stop('pheno_data is not phenotype data')
 }


 #### Get the clean geno_data ready for model fit
 #if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){
 if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){
     geno_clean <-  geno_to_model(geno_data = geno_data,
                                 train_geno_data = train_geno_data,
                                 test_geno_data = test_geno_data,
                                 message = message)

     if(attr(geno_clean, "cleared")!="for_model_fit" && all(class(geno_clean)!=c("matrix", "array", "geno_data"))) {

         stop(print(paste(msg,'Data is not fit for model')), call. = FALSE)

     }

     ## The final check is matching the pheno_clean with geno_clean if both exist.
     ## It is assumed both has to be present to build model unless the user is only
     ## interested in calculating GRM or others
     ## geno_clean and pheno_clean check to ensure they are in the same order.
     ## This is important for Baysiand and Machine learning. Not necessary for asreml
     ## The output is list with geno_data in the same order as the pheno_clean
     ## test_set if provided in the geno_clean data.

     #### Aside getting the geno data ready for model fit.
     ## It also allow user to use it for calculation of  genomic relationship matrix
     ### To pass this pheno_data for genomic relationship matrix calculation the user has to provide/supply
     ## method to calculate the matrix available in grm_calculation or kernel_calculation function
     ##

     if(((exists("geno_clean") & exists("pheno_clean"))) & (is.null(gmatrix_method) & is.null(kernel_method))){
     geno_pheno_match = pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                            object_geno = geno_clean,
                            gen_name = gen_name,
                            message = message)


         if(length(geno_pheno_match)>1){
         geno_model_ready <- geno_pheno_match[[1]]
         test_set <- geno_pheno_match[[2]]
         } else {
        geno_model_ready <- geno_pheno_match[[1]]

         }

     rm(geno_pheno_match, geno_clean)
     }

     #### To calculate the gkernel only geno_clean is acceptable
     if ((exists("geno_clean") & !is.null(kernel_method))) {
         gkernel <- kernel_calculation(
             M_matrix_clean = geno_clean,
             center=center,
             method = kernel_method,
             message = message )

         rm(geno_clean)
     }

     #### To calculate the gmatrix only geno_clean is acceptable
     if (exists("geno_clean") & !is.null(gmatrix_method)) {
        gmatrix <- grm_calculation(
        geno_clean = geno_clean,
        method=gmatrix_method)

        rm(geno_clean)
     }


 }




 # if(isFALSE(((is.null(geno_data) & is.null(train_geno_data)) & is.null(test_geno_data)))){
 #    print('ok')
 # }


######

 #### Get the clean omic1_data ready for model fit or calculation of relationship matrix
 ### To pass this pheno_data for relationship matrix calculation the user has to provide/supply
 ## method to calculate the matrix available in kernel_calculation function
 ##
 #if(isFALSE(((!is.null(omic1_data) & is.null(train_omic1_data)) & is.null(test_omic1_data)))){

 if(isFALSE(((is.null(omic1_data) & is.null(train_omic1_data)) & is.null(test_omic1_data)))){

     #print('ok')
     omic1_clean <-  omic_to_model(omic_data = omic1_data,
                                   train_omic_data = train_omic1_data,
                                   test_omic_data = test_omic1_data,
                                   message = message)

     if(attr(omic1_clean, "cleared")!="for_model_fit" && all(class(omic1_clean)!=c("matrix", "array", "omic_matrix"))) {

         stop(print(paste(msg,'Data is not fit for model')), call. = FALSE)

     }

     if((exists('omic1_clean') & exists("pheno_clean")) & (is.null(gmatrix_method) & is.null(kernel_method))){
         omic1_pheno_match = pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                             object_geno = omic1_clean,
                                             gen_name = gen_name,
                                             message = message)


         if(length(omic1_pheno_match)>1){
             omic1_model_ready <- omic1_pheno_match[[1]]
             test_set <- omic1_pheno_match[[2]]
         } else {

        omic1_model_ready <- omic1_pheno_match[[1]]
         }

         rm(omic1_pheno_match, omic1_clean)
     }

     #### To calculate the omic1_kernel only omic1_clean is acceptable
     if (exists("omic1_clean") & !is.null(kernel_method)) {
         omic1_kernel <-  kernel_calculation(
             M_matrix_clean = omic1_clean,
             center=center,
             method = kernel_method,
             message = message )

         rm(omic1_clean)
     }

     }
#########

 #### Get the clean omic2_data ready for model fit
 ### Similar condition as omic1_data applies.
 if(isFALSE(((is.null(omic2_data) & is.null(train_omic2_data)) & is.null(test_omic2_data)))){

     omic2_clean <-  omic_to_model(omic_data = omic2_data,
                                       train_omic_data = train_omic2_data,
                                       test_omic_data = test_omic2_data,
                                       message = message)

     if(attr(omic2_clean, "cleared")!="for_model_fit" && all(class(omic2_clean)!=c("matrix", "array", "omic_matrix"))) {

         stop(print(paste(msg,'Data is not fit for model')), call. = FALSE)

     }


     if((exists('omic2_clean') & exists("pheno_clean")) & (is.null(gmatrix_method) | is.null(kernel_method))){
         omic2_pheno_match = pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                              object_geno = omic2_clean,
                                              gen_name = gen_name,
                                              message = message)


         if(length(omic2_pheno_match)>1){
             omic2_model_ready <- omic2_pheno_match[[1]]
             test_set <- omic2_pheno_match[[2]]
         }else {

             omic2_model_ready <- omic2_pheno_match[[1]]
         }

         rm(omic2_pheno_match, omic2_clean)
     }

     #### To calculate the omic2_kernel only omic2_clean is acceptable
     if (exists("omic2_clean") & !is.null(kernel_method)) {
         omic2_kernel <-  kernel_calculation(
             M_matrix_clean = omic2_clean,
             center=center,
             method = kernel_method,
             message = message )

         rm(omic2_clean)
     }

 }
 #########

 #### Get the clean omic3_data ready for model fit
 if(isFALSE(((is.null(omic3_data) & is.null(train_omic3_data)) & is.null(test_omic3_data)))){

     omic3_clean <-  omic_to_model(omic_data = omic3_data,
                                       train_omic_data = train_omic3_data,
                                       test_omic_data = test_omic3_data,
                                       message = message)

     if(attr(omic3_clean, "cleared")!="for_model_fit" && all(class(omic3_clean)!=c("matrix", "array", "omic_matrix"))) {

         stop(print(paste(msg,'Data is not fit for model')), call. = FALSE)

     }


     if((exists('omic3_clean') & exists("pheno_clean")) & (is.null(gmatrix_method) & is.null(kernel_method))){
         omic3_pheno_match = pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                              object_geno = omic3_clean,
                                              gen_name = gen_name,
                                              message = message)


         if(length(omic3_pheno_match)>1){
             omic3_model_ready <- omic3_pheno_match[[1]]
             test_set <- omic3_pheno_match[[2]]
         } else {

             omic3_model_ready <- omic3_pheno_match[[1]]
         }

         rm(omic3_pheno_match, omic3_clean)
     }

     #### To calculate the omic3_kernel only omic3_clean is acceptable
     if (exists("omic3_clean") & !is.null(kernel_method)) {
         omic3_kernel <-  kernel_calculation(
             M_matrix_clean = omic3_clean,
             center=center,
             method = kernel_method,
             message = message )

         rm(omic3_clean)
     }

 }
###################################################################################
 # Pre-Check for grm/kernel matrix if calculated from the marker/omic data
 # or provided by the user.  It has to pass through this pre-check before going
 # to conditioning effect such as bend. This is important for stability of the matrix
 # during matrix inverse.
 # While it is very important and we strongly suggest user to use
 # we allow the user the opportunity to decide to use it or not.
 ##
 ## Checks
 #' #######
 #' 1. It check if the matrix is square matrix/symmetry, if not we fix it for the user
 #' 2. It check if the matrix is positive definite, if not we fix it.
 #' 3. It check for NA. If present the engine will stop further analysis.
######################################################################################3

     if(!is.null(gmatrix)){

         gmatrix_checked <- grm_kernel_precheck(pheno_data= gmatrix,
                                                message= message,
                                                bending = bending)
     }


     if(!is.null(gkernel)){

         gkernel_checked <- grm_kernel_precheck(pheno_data= gkernel,
                                                message= message,
                                                bending = bending)
     }


     if(!is.null(omic1_kernel)){

         omic1_kernel_checked <- grm_kernel_precheck(pheno_data= omic1_kernel,
                                                     message= message,
                                                     bending = bending)
     }

     if(!is.null(omic2_kernel)){

         omic2_kernel_checked <- grm_kernel_precheck(pheno_data= omic2_kernel,
                                                     message= message,
                                                     bending = bending)
     }


     if(!is.null(omic3_kernel)){

         omic3_kernel_checked <- grm_kernel_precheck(pheno_data= omic3_kernel,
                                                     message= message,
                                                     bending = bending)
     }
################################################################
 ##### Pheno to geno match
 ################################################
### This is important as the user might provide the relationship matrix differently
 if(exists("gkernel_checked")){


    gkernel_pheno_match <- pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                     object_geno = gkernel_checked,
                     gen_name = gen_name,
                     message = message)

    if(length(gkernel_pheno_match)>1){
    gkernel_model_ready <- gkernel_pheno_match[[1]]
    test_set <- gkernel_pheno_match[[2]]
    } else{
        gkernel_model_ready <- gkernel_pheno_match[[1]]

    }

    rm(gkernel_pheno_match, gkernel, gkernel_checked)
 }

 if(exists("gmatrix_checked")){


     gmatrix_pheno_match <- pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                             object_geno = gmatrix_checked,
                                             gen_name = gen_name,
                                             message = message)

     if(length(gmatrix_pheno_match)>1){
         gmatrix_model_ready <- gmatrix_pheno_match[[1]]
         test_set <- gmatrix_pheno_match[[2]]
     } else {
         gmatrix_model_ready <- gmatrix_pheno_match[[1]]

     }

     rm(gmatrix_pheno_match, gmatrix, gmatrix_checked)
 }


 if(exists("omic1_kernel_checked")){


     omic1_kernel_pheno_match <- pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                             object_geno = omic1_kernel_checked,
                                             gen_name = gen_name,
                                             message = message)

     if(length(omic1_kernel_pheno_match)>1){
         omic1_kernel_model_ready <- omic1_kernel_pheno_match[[1]]
         test_set <- omic1_kernel_pheno_match[[2]]
     } else {
         omic1_kernel_model_ready <- omic1_kernel_pheno_match[[1]]
     }

 rm(omic1_kernel_pheno_match, omic1_kernel, omic1_kernel_checked)
 }


 if(exists("omic2_kernel_checked")){


     omic2_kernel_pheno_match <- pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                                  object_geno = omic2_kernel_checked,
                                                  gen_name = gen_name,
                                                  message = message)

     if(length(omic2_kernel_pheno_match)>1){
         omic2_kernel_model_ready <- omic2_kernel_pheno_match[[1]]
         test_set <- omic2_kernel_pheno_match[[2]]
     } else {
         omic2_kernel_model_ready <- omic2_kernel_pheno_match[[1]]
     }

     rm(omic2_kernel_pheno_match, omic2_kernel, omic2_kernel_checked)
 }

 if(exists("omic3_kernel_checked")){


     omic3_kernel_pheno_match <- pheno_geno_match(object_pheno = pheno_clean$pheno_data,
                                                  object_geno = omic3_kernel_checked,
                                                  gen_name = gen_name,
                                                  message = message)

     if(length(omic3_kernel_pheno_match)>1){
         omic3_kernel_model_ready <- omic3_kernel_pheno_match[[1]]
         test_set <- omic3_kernel_pheno_match[[2]]
     } else {
         omic3_kernel_model_ready <- omic3_kernel_pheno_match[[1]]
     }

     rm(omic3_kernel_pheno_match, omic3_kernel, omic3_kernel_checked)
 }


 ### End
#####################################################################


 ##########################################################################
 #########################################################################
 ## Start of Bayes A, B, C, BL and BRR Models for Single Location       ##                     ##                          ##
 ##  This only accommodate n x p matrix  not nxn                        ##                                      ##
 ##                                                                     ##
 ##########################################################################
 #######################################################################


    ### Model BRR for single location
 #if(((GS_model=="BRR") & is.null(rand_term_model_bayesian)) || ((is.null(GS_model) & (rand_term_model_bayesian=="BRR")))){

    # if((isTRUE(GS_model== "BRR" | GS_model== "BayesA"|  GS_model== "BayesB"| GS_model== "BayesC" | GS_model== "BL") & is.null(rand_term_model_bayesian)) |
    #    ((is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL"))))) |
    #    ((!is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL")))))){

 #### These models only works with one environment/location
 if(length(pheno_clean[[1]][,gen_name])==length(unique(pheno_clean[[1]][,gen_name]))){

 if((isTRUE(GS_model== "BRR" | GS_model== "BayesA"|  GS_model== "BayesB"| GS_model== "BayesC" | GS_model== "BL") & is.null(rand_term_model_bayesian)) |
    (is.null(GS_model) & length(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL"))!=0) |
    (!is.null(GS_model) & length(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL"))!=0)){


        if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

### The ETA_compiler_bayes compile the linear predictors and set parameters for the model.
## Check the function for details.
   ETA  <-  ETA_compiler_bayes(
        fixed = fixed,
        random = random,
        GS_model = GS_model,
        fixed_term_model_bayesian = fixed_term_model_bayesian,
        rand_term_model_bayesian = rand_term_model_bayesian,
        pheno_data = pheno_clean[[1]],
        geno_data = geno_model_ready,
        omic1_data = NULL,
        omic2_data = NULL,
        omic3_data = NULL,
        gen_name = gen_name)

   bayes_para <-  bayes_parameter_check(nIter = nIter,
                                        burnIn = burnIn,
                                        thin = thin)

    mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                          response = response,
                                          weights = weights,
                                          ETA = ETA$ETA,
                                          bayes_para = bayes_para,
                                          verbose = FALSE
                                          #files_key = "files_key"
                                          )

    res_model_output <- mod_output_bayes(mod = mod,
                                         ETA = ETA,
                                         geno_data = geno_model_ready,
                                         gen_name = gen_name,
                                         omic1_data = NULL,
                                         omic2_data = NULL,
                                         omic3_data = NULL,
                                         GS_model = GS_model)

    res_summary_stat <- summary_statistics_bayes(mod = mod)

   res_plot <-  plot_acc(mod = mod, response = response)


        }
        ###### omic1_clean

        if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = omic1_model_ready,
                omic2_data = NULL,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)
 ## M_matrix_bayes_mod_single_loc
            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = NULL,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

        ##### omic2_model_ready

        if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = NULL,
                omic2_data = omic2_model_ready,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

        ####

        ##### omic3_model_ready

        if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = NULL,
                omic2_data = NULL,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = NULL,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### geno_model_ready and omic1_model_ready

        if((exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = omic1_model_ready,
                omic2_data = NULL,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = NULL,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### geno_clean and omic2_clean

        if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = NULL,
                omic2_data = omic2_model_ready,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

        ##### geno_clean and omic3_clean

        if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = NULL,
                omic2_data = NULL,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = NULL,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### omic1_model_ready and omic2_model_ready

        if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = omic1_model_ready,
                omic2_data = omic2_model_ready,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### omic1_model_ready and omic3_model_ready

        if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = omic1_model_ready,
                omic2_data = NULL,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = NULL,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ####
        ##### omic2_model_ready and omic3_model_ready

        if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = NULL,
                omic2_data = omic2_model_ready,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

        ########

        ##### geno_model_ready, omic1_model_ready and omic2_model_ready

        if((exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = omic1_model_ready,
                omic2_data = omic2_model_ready,
                omic3_data = NULL,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = NULL,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### geno_model_ready, omic1_model_ready and omic3_model_ready

        if((exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = omic1_model_ready,
                omic2_data = NULL,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights =weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = NULL,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }


    #####

        ##### geno_model_ready, omic2_model_ready and omic3_model_ready

        if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = NULL,
                omic2_data = omic2_model_ready,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = NULL,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

        ##### omic1_model_ready, omic2_model_ready and omic3_model_ready

        if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = NULL,
                omic1_data = omic1_model_ready,
                omic2_data = omic2_model_ready,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                                  response = response,
                                                  weights = weights,
                                                  ETA = ETA$ETA,
                                                  bayes_para = bayes_para,
                                                  verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = NULL,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


            res_plot <-  plot_acc(mod = mod, response = response)

        }


        ##### geno,  omic1_clean, omic2clean and omic3_clean

        if((exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

            ETA  <-  ETA_compiler_bayes(
                fixed = fixed,
                random = random,
                GS_model = GS_model,
                fixed_term_model_bayesian = fixed_term_model_bayesian,
                rand_term_model_bayesian = rand_term_model_bayesian,
                pheno_data = pheno_clean[[1]],
                geno_data = geno_model_ready,
                omic1_data = omic1_model_ready,
                omic2_data = omic2_model_ready,
                omic3_data = omic3_model_ready,
                gen_name = gen_name)

            bayes_para <-  bayes_parameter_check(nIter = nIter,
                                                 burnIn = burnIn,
                                                 thin = thin)

            mod <-  bayes_mod_execute(pheno_data = pheno_data,
                                      response = response,
                                      weights = weights,
                                      ETA = ETA$ETA,
                                      bayes_para = bayes_para,
                                      verbose = FALSE
                                                  )

            res_model_output <- mod_output_bayes(mod = mod,
                                                 ETA = ETA,
                                                 geno_data = geno_model_ready,
                                                 gen_name = gen_name,
                                                 omic1_data = omic1_model_ready,
                                                 omic2_data = omic2_model_ready,
                                                 omic3_data = omic3_model_ready,
                                                 GS_model = GS_model
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)

            res_plot <-  plot_acc(mod = mod, response = response)

        }

 }
} ## End of  Bayes A, B, C, BRR, BL

 ##########################################################################
 #########################################################################
 ## Start of RKHS, (BRR- Bayesian GBLUP ) and GBLUP (asreml) Model
 ## for Single Location and multiple loc                                ##
 ##                                                                     ##
 ##                                                                     ##
 ##########################################################################
 #######################################################################

 # if((isTRUE(GS_model== "RKHS") |
 #     ((is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%"RKHS"))) |
 #     ((!is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%"RKHS"))))){

 # if((isTRUE(GS_model== "RKHS" | isTRUE(GS_model== "BRR")) & is.null(rand_term_model_bayesian)) |
 #    ((is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%c("RKHS", "BRR")))) |
 #    ((!is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%c("RKHS", "BRR"))))){
 if((isTRUE(GS_model== "RKHS" | isTRUE(GS_model== "BRR") | isTRUE(GS_model== "GBLUP")) & is.null(rand_term_model_bayesian)) |
    (is.null(GS_model) & length(rand_term_model_bayesian%in%c("RKHS", "BRR"))!=0) |
    (!is.null(GS_model) & length(rand_term_model_bayesian%in%c("RKHS", "BRR"))!=0)){


     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model=="BRR" | GS_model=="RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel =  gkernel_model_ready,
                 gen_name = gen_name)

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     gen_name = gen_name)

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin)

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(exists("gkernel_model_ready")){
             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gkernel =  gkernel_model_ready,
                                                       gen_name = gen_name,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gkernel =  gkernel_model_ready,
                                                           gen_name = gen_name,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }

         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gmatrix = gmatrix_model_ready,
                                                           gen_name = gen_name,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }

                 if(GS_model=="BRR"){
                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gmatrix = gmatrix_model_ready,
                                                               gen_name = gen_name,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )
                 }
             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine == 'asreml'){

                 if(exists('gmatrix_model_ready')){
                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     gmatrix = gmatrix_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

        res_model_output <- asreml_mod_output(
                                    mod_asreml = mod,
                                    pheno_data = pheno_clean[[1]],
                                    gmatrix = gmatrix_model_ready,
                                    heter_groups = heter_groups,
                                    gen_name = gen_name,
                                    VarCov_str = VarCov_str,
                                    heter_groups = heter_groups,
                                    heter_resid = heter_resid,
                                    pworkspace= pworkspace,
                                    workspace = workspace,
                                    maxit = maxit
                                             )

                 } else {

                     if(exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             pworkspace= pworkspace,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             heter_groups = heter_groups,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid,
                             pworkspace= pworkspace,
                             workspace = workspace,
                             maxit = maxit
                         )
                     }

                 }

        # res_plot <-  plot_acc(mod = res_summary_stat$mod,
        #                       response = response)

             }

         }

     }
     ###### omic1_clean

     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(GS_model=="BRR" | GS_model=="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic1_kernel =  omic1_kernel_model_ready,
             gen_name = gen_name
         )

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic1_kernel = omic1_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )
         }
         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
           if(GS_model=="GBLUP" & engine =="asreml")  {

               mod = asreml_utilis(fixed = fixed,
                                   random = random,
                                   cova=cova,
                                   GS_model = GS_model,
                                   response = response,
                                   pheno_data = pheno_clean[[1]],
                                   omic1_kernel = omic1_kernel_model_ready,
                                   gen_name = gen_name,
                                   heter_groups = heter_groups,
                                   heter_resid = heter_resid,
                                   VarCov_str = VarCov_str,
                                   weights = weights,
                                   core = core,
                                   pworkspace= pworkspace,
                                   workspace = workspace,
                                   maxit = maxit)

               res_model_output <- asreml_mod_output(
                   mod_asreml = mod,
                   pheno_data = pheno_clean[[1]],
                   omic1_kernel = omic1_kernel_model_ready,
                   gen_name = gen_name,
                   VarCov_str = VarCov_str,
                   heter_groups = heter_groups,
                   heter_resid = heter_resid,
                   pworkspace= pworkspace,
                   workspace = workspace,
                   maxit = maxit
               )

               # res_plot <-  plot_acc(mod = res_summary_stat$mod,
               #                       response = response)
           }

         }

     }

     ##### omic2_model_ready

     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(GS_model=="BRR" | GS_model=="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic2_kernel = omic2_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic2_kernel = omic2_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }
         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )
         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic2_kernel = omic2_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic2_kernel = omic2_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }

     ####

     ##### omic3_model_ready

     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(GS_model == "BRR" | GS_model == "RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic3_kernel = omic3_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic3_kernel = omic3_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){

             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )
         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic3_kernel = omic3_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


     ##### geno_model_ready and omic1_model_ready

     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model=="BRR" | GS_model=="RKHS"){

         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic1_kernel = omic1_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para)

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )
             }

         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }


                 if(GS_model=="BRR"){

                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic1_kernel = omic1_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     gmatrix = gmatrix_model_ready,
                                     omic1_kernel = omic1_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 } else {

                if( exists('gkernel_model_ready')){
                    mod = asreml_utilis(fixed = fixed,
                                        random = random,
                                        cova=cova,
                                        GS_model = GS_model,
                                        response = response,
                                        pheno_data = pheno_clean[[1]],
                                        gkernel = gkernel_model_ready,
                                        omic1_kernel = omic1_kernel_model_ready,
                                        gen_name = gen_name,
                                        heter_groups = heter_groups,
                                        heter_resid = heter_resid,
                                        VarCov_str = VarCov_str,
                                        weights = weights,
                                        core = core,
                                        workspace = workspace,
                                        maxit = maxit)

                    res_model_output <- asreml_mod_output(
                        mod_asreml = mod,
                        pheno_data = pheno_clean[[1]],
                        gkernel = gkernel_model_ready,
                        omic1_kernel = omic1_kernel_model_ready,
                        gen_name = gen_name,
                        VarCov_str = VarCov_str,
                        heter_groups = heter_groups,
                        heter_resid = heter_resid,
                        pworkspace= pworkspace,
                        workspace = workspace,
                        maxit = maxit
                    )

                }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


     ##### geno_clean and omic2_clean

     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model =="BRR" | GS_model =="RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic2_kernel = omic2_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }

         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }

                 if(GS_model=="BRR"){
                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic2_kernel = omic2_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gmatrix = gmatrix_model_ready,
                                         omic2_kernel = omic2_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gmatrix = gmatrix_model_ready,
                         omic2_kernel = omic2_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 } else {

                     if( exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             omic2_kernel = omic2_kernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             pworkspace= pworkspace,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             omic2_kernel = omic2_kernel_model_ready,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid,
                             pworkspace= pworkspace,
                             workspace = workspace,
                             maxit = maxit
                         )

                     }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }
     ##### geno_clean and omic3_clean

     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model=="BRR" | GS_model =="RKHS"){

         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic3_kernel = omic3_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }

         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){

                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }


                 if(GS_model=="BRR"){

                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic3_kernel = omic3_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }


             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)


         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gmatrix = gmatrix_model_ready,
                                         omic3_kernel = omic3_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gmatrix = gmatrix_model_ready,
                         omic3_kernel = omic3_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 } else {

                     if( exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             omic3_kernel = omic3_kernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             omic3_kernel = omic3_kernel_model_ready,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid
                         )

                     }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


     ##### omic1_model_ready and omic2_model_ready

     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(GS_model =="BRR" | GS_model =="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic1_kernel = omic1_kernel_model_ready,
             omic2_kernel = omic2_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic1_kernel = omic1_kernel_model_ready,
                                                   omic2_kernel = omic2_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)


         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic1_kernel = omic1_kernel_model_ready,
                                     omic2_kernel = omic2_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic1_kernel = omic1_kernel_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }



     ##### omic1_model_ready and omic3_model_ready


     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(GS_model =="BRR" | GS_model =="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic1_kernel = omic1_kernel_model_ready,
             omic3_kernel = omic3_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic1_kernel = omic1_kernel_model_ready,
                                                   omic3_kernel = omic3_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

         }


         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic1_kernel = omic1_kernel_model_ready,
                                     omic3_kernel = omic3_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic1_kernel = omic1_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


     ####
     ##### omic2_model_ready and omic3_model_ready


     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(GS_model=="BRR" | GS_model =="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic2_kernel = omic2_kernel_model_ready,
             omic3_kernel = omic3_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic2_kernel = omic2_kernel_model_ready,
                                                   omic3_kernel = omic3_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic2_kernel = omic2_kernel_model_ready,
                                     omic3_kernel = omic3_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic2_kernel = omic2_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }

     ########

     ##### geno_model_ready, omic1_model_ready and omic2_model_ready

     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model =="BRR" | GS_model=="RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic1_kernel = omic1_kernel_model_ready,
                 omic2_kernel = omic2_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }


         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }

                 if(GS_model=="BRR"){
                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic1_kernel = omic1_kernel_model_ready,
                                                               omic2_kernel = omic2_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {

             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gmatrix = gmatrix_model_ready,
                                         omic1_kernel = omic1_kernel_model_ready,
                                         omic2_kernel = omic2_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gmatrix = gmatrix_model_ready,
                         omic1_kernel = omic1_kernel_model_ready,
                         omic2_kernel = omic2_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 } else {

                     if( exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             omic1_kernel = omic1_kernel_model_ready,
                                             omic2_kernel = omic2_kernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             pworkspace= pworkspace,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             omic1_kernel = omic1_kernel_model_ready,
                             omic2_kernel = omic2_kernel_model_ready,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid,
                             pworkspace= pworkspace,
                             workspace = workspace,
                             maxit = maxit
                         )

                     }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }
     }

     }


     ##### geno_model_ready, omic1_model_ready and omic3_model_ready


     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model == "BRR" | GS_model == "RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic1_kernel = omic1_kernel_model_ready,
                 omic3_kernel = omic3_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para)

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }

         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }


                 if(GS_model=="BRR"){
                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic1_kernel = omic1_kernel_model_ready,
                                                               omic3_kernel = omic3_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

     } else {
         if(GS_model=="GBLUP" & engine =="asreml")  {
             if(exists('gmatrix_model_ready')){
                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     gmatrix = gmatrix_model_ready,
                                     omic1_kernel = omic1_kernel_model_ready,
                                     omic3_kernel = omic3_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

             } else {

                 if( exists('gkernel_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gkernel = gkernel_model_ready,
                                         omic1_kernel = omic1_kernel_model_ready,
                                         omic3_kernel = omic3_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gkernel = gkernel_model_ready,
                         omic1_kernel = omic1_kernel_model_ready,
                         omic3_kernel = omic3_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 }
             }

             # res_plot <-  plot_acc(mod = res_summary_stat$mod,
             #                       response = response)
         }

     }

   }

     #####

     ##### geno_model_ready, omic2_model_ready and omic3_model_ready


     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model =="BRR" | GS_model == "RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic2_kernel = omic2_kernel_model_ready,
                 omic3_kernel = omic3_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para)

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }


             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }


         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){

                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }

                 if(GS_model=="BRR"){

                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic2_kernel = omic2_kernel_model_ready,
                                                               omic3_kernel = omic3_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gmatrix = gmatrix_model_ready,
                                         omic2_kernel = omic2_kernel_model_ready,
                                         omic3_kernel = omic3_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gmatrix = gmatrix_model_ready,
                         omic2_kernel = omic2_kernel_model_ready,
                         omic3_kernel = omic3_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 } else {

                     if( exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             omic2_kernel = omic2_kernel_model_ready,
                                             omic3_kernel = omic3_kernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             pworkspace= pworkspace,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             omic2_kernel = omic2_kernel_model_ready,
                             omic3_kernel = omic3_kernel_model_ready,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid,
                             pworkspace= pworkspace,
                             workspace = workspace,
                             maxit = maxit
                         )

                     }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


     ##### omic1_model_ready, omic2_model_ready and omic3_model_ready

     if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(GS_model =="BRR" | GS_model =="RKHS"){
         ETA  <-  ETA_compiler_bayes_GBLUP(
             fixed = fixed,
             random = random,
             GS_model = GS_model,
             fixed_term_model_bayesian = fixed_term_model_bayesian,
             rand_term_model_bayesian = rand_term_model_bayesian,
             pheno_data = pheno_clean[[1]],
             omic1_kernel = omic1_kernel_model_ready,
             omic2_kernel = omic2_kernel_model_ready,
             omic3_kernel = omic3_kernel_model_ready,
             gen_name = gen_name)

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(GS_model=="RKHS"){
         res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                   ETA = ETA,
                                                   gen_name = gen_name,
                                                   omic1_kernel = omic1_kernel_model_ready,
                                                   omic2_kernel = omic2_kernel_model_ready,
                                                   omic3_kernel = omic3_kernel_model_ready,
                                                   pheno_data = ETA$pheno_data,
                                                   heter_groups  = heter_groups
         )

         }

         if(GS_model=="BRR"){
             res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {

                 mod = asreml_utilis(fixed = fixed,
                                     random = random,
                                     cova=cova,
                                     GS_model = GS_model,
                                     response = response,
                                     pheno_data = pheno_clean[[1]],
                                     omic1_kernel = omic1_kernel_model_ready,
                                     omic2_kernel = omic2_kernel_model_ready,
                                     omic3_kernel = omic3_kernel_model_ready,
                                     gen_name = gen_name,
                                     heter_groups = heter_groups,
                                     heter_resid = heter_resid,
                                     VarCov_str = VarCov_str,
                                     weights = weights,
                                     core = core,
                                     pworkspace= pworkspace,
                                     workspace = workspace,
                                     maxit = maxit)

                 res_model_output <- asreml_mod_output(
                     mod_asreml = mod,
                     pheno_data = pheno_clean[[1]],
                     omic1_kernel = omic1_kernel_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name,
                     VarCov_str = VarCov_str,
                     heter_groups = heter_groups,
                     heter_resid = heter_resid,
                     pworkspace= pworkspace,
                     workspace = workspace,
                     maxit = maxit
                 )

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }


     }



     ##### geno, omic1_clean, omic2clean and omic3_clean

     if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

         if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

         if(GS_model=="BRR" | GS_model =="RKHS"){
         if(exists('gkernel_model_ready')){
             ETA  <-  ETA_compiler_bayes_GBLUP(
                 fixed = fixed,
                 random = random,
                 GS_model = GS_model,
                 fixed_term_model_bayesian = fixed_term_model_bayesian,
                 rand_term_model_bayesian = rand_term_model_bayesian,
                 pheno_data = pheno_clean[[1]],
                 gkernel = gkernel_model_ready,
                 omic1_kernel = omic1_kernel_model_ready,
                 omic2_kernel = omic2_kernel_model_ready,
                 omic3_kernel = omic3_kernel_model_ready,
                 gen_name = gen_name
             )

         } else {

             if(exists('gmatrix_model_ready')){
                 ETA  <-  ETA_compiler_bayes_GBLUP(
                     fixed = fixed,
                     random = random,
                     GS_model = GS_model,
                     fixed_term_model_bayesian = fixed_term_model_bayesian,
                     rand_term_model_bayesian = rand_term_model_bayesian,
                     pheno_data = pheno_clean[[1]],
                     gmatrix = gmatrix_model_ready,
                     omic1_kernel = omic1_kernel_model_ready,
                     omic2_kernel = omic2_kernel_model_ready,
                     omic3_kernel = omic3_kernel_model_ready,
                     gen_name = gen_name
                 )

             }

         }

         bayes_para <-  bayes_parameter_check(nIter = nIter,
                                              burnIn = burnIn,
                                              thin = thin
         )

         ## bayes_mod_execute

         mod <-  bayes_mod_execute(pheno_data = ETA$pheno_data,
                                   response = response,
                                   weights = weights,
                                   ETA = ETA$ETA,
                                   bayes_para = bayes_para
         )

         if(exists('gkernel_model_ready')){

             if(GS_model=="RKHS"){
             res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                       ETA = ETA,
                                                       gen_name = gen_name,
                                                       gkernel = gkernel_model_ready,
                                                       omic1_kernel = omic1_kernel_model_ready,
                                                       omic2_kernel = omic2_kernel_model_ready,
                                                       omic3_kernel = omic3_kernel_model_ready,
                                                       pheno_data = ETA$pheno_data,
                                                       heter_groups  = heter_groups
             )

             }

             if(GS_model=="BRR"){
                 res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gkernel = gkernel_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

             }


         } else {

             if(exists('gmatrix_model_ready')){

                 if(GS_model=="RKHS"){
                 res_model_output <- mod_output_bayes_RKHS(mod = mod,
                                                           ETA = ETA,
                                                           gen_name = gen_name,
                                                           gmatrix = gmatrix_model_ready,
                                                           omic1_kernel = omic1_kernel_model_ready,
                                                           omic2_kernel = omic2_kernel_model_ready,
                                                           omic3_kernel = omic3_kernel_model_ready,
                                                           pheno_data = ETA$pheno_data,
                                                           heter_groups  = heter_groups
                 )

                 }

                 if(GS_model=="BRR"){
                     res_model_output <- mod_output_bayes_BRRGBLUP(mod = mod,
                                                               ETA = ETA,
                                                               gen_name = gen_name,
                                                               gmatrix = gmatrix_model_ready,
                                                               omic1_kernel = omic1_kernel_model_ready,
                                                               omic2_kernel = omic2_kernel_model_ready,
                                                               omic3_kernel = omic3_kernel_model_ready,
                                                               pheno_data = ETA$pheno_data,
                                                               heter_groups  = heter_groups
                     )

                 }

             }

         }

         res_summary_stat <- summary_statistics_bayes(mod = mod)

         res_plot <-  plot_acc(mod = mod, response = response)

         } else {
             if(GS_model=="GBLUP" & engine =="asreml")  {
                 if(exists('gmatrix_model_ready')){
                     mod = asreml_utilis(fixed = fixed,
                                         random = random,
                                         cova=cova,
                                         GS_model = GS_model,
                                         response = response,
                                         pheno_data = pheno_clean[[1]],
                                         gmatrix = gmatrix_model_ready,
                                         omic1_kernel = omic1_kernel_model_ready,
                                         omic2_kernel = omic2_kernel_model_ready,
                                         omic3_kernel = omic3_kernel_model_ready,
                                         gen_name = gen_name,
                                         heter_groups = heter_groups,
                                         heter_resid = heter_resid,
                                         VarCov_str = VarCov_str,
                                         weights = weights,
                                         core = core,
                                         pworkspace= pworkspace,
                                         workspace = workspace,
                                         maxit = maxit)

                     res_model_output <- asreml_mod_output(
                         mod_asreml = mod,
                         pheno_data = pheno_clean[[1]],
                         gmatrix = gmatrix_model_ready,
                         omic1_kernel = omic1_kernel_model_ready,
                         omic2_kernel = omic2_kernel_model_ready,
                         omic3_kernel = omic3_kernel_model_ready,
                         gen_name = gen_name,
                         VarCov_str = VarCov_str,
                         heter_groups = heter_groups,
                         heter_resid = heter_resid,
                         pworkspace= pworkspace,
                         workspace = workspace,
                         maxit = maxit
                     )

                 } else {

                     if( exists('gkernel_model_ready')){
                         mod = asreml_utilis(fixed = fixed,
                                             random = random,
                                             cova=cova,
                                             GS_model = GS_model,
                                             response = response,
                                             pheno_data = pheno_clean[[1]],
                                             gkernel = gkernel_model_ready,
                                             omic1_kernel = omic1_kernel_model_ready,
                                             omic2_kernel = omic2_kernel_model_ready,
                                             omic3_kernel = omic3_kernel_model_ready,
                                             gen_name = gen_name,
                                             heter_groups = heter_groups,
                                             heter_resid = heter_resid,
                                             VarCov_str = VarCov_str,
                                             weights = weights,
                                             core = core,
                                             pworkspace= pworkspace,
                                             workspace = workspace,
                                             maxit = maxit)

                         res_model_output <- asreml_mod_output(
                             mod_asreml = mod,
                             pheno_data = pheno_clean[[1]],
                             gkernel = gkernel_model_ready,
                             omic1_kernel = omic1_kernel_model_ready,
                             omic2_kernel = omic2_kernel_model_ready,
                             omic3_kernel = omic3_kernel_model_ready,
                             gen_name = gen_name,
                             VarCov_str = VarCov_str,
                             heter_groups = heter_groups,
                             heter_resid = heter_resid,
                             pworkspace= pworkspace,
                             workspace = workspace,
                             maxit = maxit
                         )

                     }
                 }

                 # res_plot <-  plot_acc(mod = res_summary_stat$mod,
                 #                       response = response)
             }

         }

     }


 } #### END GBLUP_RKHS, GBLUP_BRR and GBLUP (asreml)

 ######################################################
 ######################################################
 ##                                                 ###
 ## Machine Learning Models                         ###
 ##                                                 ###
 ######################################################
 ######################################################

 ##########################################################################
 ############################################################################
 ## Start data organization for ML model fitting
 ##
 ##############################################################################
 ##############################################################################

 #if(exists("pheno_clean")){

 # if(length(pheno_clean)==2){
 #
 #   object_pheno <- object_pheno$pheno_data
 #
 #      test_set <- object_pheno$test_set
 # } else {

 #### These models only works with one environment/location
 if(length(pheno_clean[[1]][,gen_name])==length(unique(pheno_clean[[1]][,gen_name]))){

 if(isTRUE(GS_model== "Xgboost" | GS_model== "RandomForest" | GS_model== "PartialLeastSquare" | GS_model== "SupportVectorMachine" | GS_model== "K-NearestNeighbors" | GS_model=="Lasso" | GS_model=="Ridge_Regression"))   {


 pheno_clean <- ML_undefined_test_train(object_pheno = pheno_clean,
                                        response = response)

 if(length(pheno_clean)==1){
     pheno_clean <- pheno_clean$pheno_data
     if(length(pheno_clean)==2){

         pheno_clean <- pheno_clean$pheno_data

         test_set_ <- pheno_clean$test_set
     }

 }

 #}
 ####################################

 if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         rm(geno_model_ready)
     } else {

         geno_model_ready_train <-  geno_model_ready

         rm(geno_model_ready)
     }


 }
 #####
 if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         rm(omic1_model_ready)

     } else {

         omic1_model_ready_train <-  omic1_model_ready

         rm(omic1_model_ready)
     }


 }
 ############
 if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         rm(omic2_model_ready)
     } else {

         omic2_model_ready_train <-  omic2_model_ready

         rm(omic2_model_ready)
     }



 }
 #############

 if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         rm(omic3_model_ready)

     }  else {

         omic3_model_ready_train <-  omic3_model_ready

         rm(omic3_model_ready)
     }


 }

 ########

 ##### geno_model_ready and omic1_model_ready

 if((exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         geno_omic1_test = cbind(geno_model_ready_test, omic1_model_ready_test)

         geno_omic1_train = cbind(geno_model_ready_train, omic1_model_ready_train)

         rm(geno_model_ready, omic1_model_ready,
            geno_model_ready_test, omic1_model_ready_test,
            geno_model_ready_train, omic1_model_ready_train)

     } else {

         geno_omic1_train = cbind(geno_model_ready, omic1_model_ready)

         rm(geno_model_ready, omic1_model_ready)
     }

 }
 ################

 ##### geno_model_ready and omic2_model_ready

 if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         geno_omic2_test = cbind(geno_model_ready_test, omic2_model_ready_test)

         geno_omic2_train = cbind(geno_model_ready_train, omic2_model_ready_train)

         rm(geno_model_ready, omic2_model_ready,
            geno_model_ready_test, omic2_model_ready_test,
            geno_model_ready_train, omic2_model_ready_train)

     } else {

         geno_omic2_train = cbind(geno_model_ready, omic2_model_ready)

         rm(geno_model_ready, omic2_model_ready)
     }

 }


 ################
 ##### geno_model_ready and omic3_model_ready

 if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         geno_omic3_test = cbind(geno_model_ready_test, omic3_model_ready_test)

         geno_omic3_train = cbind(geno_model_ready_train, omic3_model_ready_train)

         rm(geno_model_ready, omic3_model_ready,
            geno_model_ready_test, omic3_model_ready_test,
            geno_model_ready_train, omic3_model_ready_train)

     } else {

         geno_omic3_train = cbind(geno_model_ready, omic3_model_ready)

         rm(geno_model_ready, omic3_model_ready)
     }

 }

 ####
 ################
 #####  omic1_model_ready and omic2_model_ready

 if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic1_omic2_test = cbind(omic1_model_ready_test, omic2_model_ready_test)

         omic1_omic2_train = cbind(omic1_model_ready_train, omic2_model_ready_train)

         rm(omic1_model_ready, omic2_model_ready,
            omic1_model_ready_test, omic2_model_ready_test,
            omic1_model_ready_train, omic2_model_ready_train)

     } else {

         omic1_omic2_train = cbind(omic1_model_ready, omic2_model_ready)

         rm(omic1_model_ready, omic2_model_ready)
     }

 }

 ################
 #####  omic1_model_ready and omic3_model_ready

 if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic1_omic3_test = cbind(omic1_model_ready_test, omic3_model_ready_test)

         omic1_omic3_train = cbind(omic1_model_ready_train, omic3_model_ready_train)

         rm(omic1_model_ready, omic3_model_ready,
            omic1_model_ready_test, omic3_model_ready_test,
            omic1_model_ready_train, omic3_model_ready_train)

     } else {

         omic1_omic3_train = cbind(omic1_model_ready, omic3_model_ready)

         rm(omic1_model_ready, omic3_model_ready)
     }

 }


 ################
 #####  omic1_model_ready and omic3_model_ready

 if((!exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic2_omic3_test = cbind(omic2_model_ready_test, omic3_model_ready_test)

         omic2_omic3_train = cbind(omic2_model_ready_train, omic3_model_ready_train)

         rm(omic2_model_ready, omic3_model_ready,
            omic2_model_ready_test, omic3_model_ready_test,
            omic2_model_ready_train, omic3_model_ready_train)

     } else {

         omic2_omic3_train = cbind(omic2_model_ready, omic3_model_ready)

         rm(omic2_model_ready, omic3_model_ready)
     }

 }

 ################
 ##### geno_model_ready,  omic1_model_ready and omic2_model_ready

 if((exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]
         ####

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         geno_omic1_omic2_test = scale(cbind(cbind(geno_model_ready_test, omic1_model_ready_test), omic2_model_ready_test))

         geno_omic1_omic2_train = scale(cbind(cbind(geno_model_ready_train,omic1_model_ready_train), omic2_model_ready_train))

         rm(geno_model_ready, omic1_model_ready, omic2_model_ready,
            geno_model_ready_test, omic1_model_ready_test, omic2_model_ready_test,
            geno_model_ready_train, omic1_model_ready_train, omic2_model_ready_train)

     } else {

         geno_omic1_omic2_train = scale(cbind(cbind(geno_model_ready, omic1_model_ready), omic2_model_ready))

         rm(geno_model_ready, omic1_model_ready, omic2_model_ready)
     }

 }


 ################
 ##### geno_model_ready,  omic1_model_ready and omic3_model_ready

 if((exists('geno_model_ready') & (exists('omic1_model_ready') & (!exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]
         ####

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         geno_omic1_omic3_test = scale(cbind(cbind(geno_model_ready_test, omic1_model_ready_test), omic3_model_ready_test))

         geno_omic1_omic3_train = scale(cbind(cbind(geno_model_ready_train,omic1_model_ready_train), omic3_model_ready_train))

         rm(geno_model_ready, omic1_model_ready, omic3_model_ready,
            geno_model_ready_test, omic1_model_ready_test, omic3_model_ready_test,
            geno_model_ready_train, omic1_model_ready_train, omic3_model_ready_train)

     } else {

         geno_omic1_omic3_train = scale(cbind(cbind(geno_model_ready, omic1_model_ready), omic3_model_ready))

         rm(geno_model_ready, omic1_model_ready, omic3_model_ready)
     }

 }


 ################
 ##### geno_model_ready,  omic2_model_ready and omic3_model_ready

 if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]
         ####

         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         geno_omic2_omic3_test = scale(cbind(cbind(geno_model_ready_test, omic2_model_ready_test), omic3_model_ready_test))

         geno_omic2_omic3_train = scale(cbind(cbind(geno_model_ready_train,omic2_model_ready_train), omic3_model_ready_train))

         rm(geno_model_ready, omic2_model_ready, omic3_model_ready,
            geno_model_ready_test, omic2_model_ready_test, omic3_model_ready_test,
            geno_model_ready_train, omic2_model_ready_train, omic3_model_ready_train)

     } else {

         geno_omic2_omic3_train = scale(cbind(cbind(geno_model_ready, omic2_model_ready), omic3_model_ready))

         rm(geno_model_ready, omic2_model_ready, omic3_model_ready)
     }

 }
 #########
 #### omic1, omic2, omic 3

 if((!exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){


         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         ####
         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic1_omic2_omic3_test = scale(cbind(cbind(omic1_model_ready_test, omic2_model_ready_test), omic3_model_ready_test))

         omic1_omic2_omic3_train = scale(cbind(cbind(omic1_model_ready_train,omic2_model_ready_train), omic3_model_ready_train))

         rm(omic1_model_ready, omic2_model_ready, omic3_model_ready,
            omic1_model_ready_test, omic2_model_ready_test, omic3_model_ready_test,
            omic1_model_ready_train, omic2_model_ready_train, omic3_model_ready_train)

     } else {

         geno_omic1_omic3_train = scale(cbind(cbind(geno_model_ready, omic1_model_ready), omic3_model_ready))

         rm(geno_model_ready, omic1_model_ready, omic3_model_ready)
     }

 }


 ################
 ##### geno_model_ready, omic1_model_ready, omic2_model_ready and omic3_model_ready

 if((exists('geno_model_ready') & (exists('omic1_model_ready') & (exists('omic2_model_ready') & exists('omic3_model_ready'))))){

     if(exists("test_set_")){

         geno_model_ready_test <-  geno_model_ready[rownames(geno_model_ready)%in%test_set_[, gen_name], ]

         geno_model_ready_train <-  geno_model_ready[!rownames(geno_model_ready)%in%test_set_[, gen_name], ]
         ####

         omic1_model_ready_test <-  omic1_model_ready[rownames(omic1_model_ready)%in%test_set_[, gen_name], ]

         omic1_model_ready_train <-  omic1_model_ready[!rownames(omic1_model_ready)%in%test_set_[, gen_name], ]
         ####

         omic2_model_ready_test <-  omic2_model_ready[rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         omic2_model_ready_train <-  omic2_model_ready[!rownames(omic2_model_ready)%in%test_set_[, gen_name], ]

         ###
         omic3_model_ready_test <-  omic3_model_ready[rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         omic3_model_ready_train <-  omic3_model_ready[!rownames(omic3_model_ready)%in%test_set_[, gen_name], ]

         geno_omic1_omic2_omic3_test = scale(cbind(cbind(geno_model_ready_test, omic1_model_ready_test),
                                                   cbind(omic2_model_ready_test, omic3_model_ready_test)))

         geno_omic1_omic2_omic3_train = scale(cbind(cbind(geno_model_ready_train,omic1_model_ready_train),
                                                    cbind(omic2_model_ready_train, omic3_model_ready_train)))

         rm(geno_model_ready, omic1_model_ready, omic2_model_ready, omic3_model_ready,
            geno_model_ready_test, omic1_model_ready_test, omic2_model_ready_test, omic3_model_ready_test,
            geno_model_ready_train, omic1_model_ready_train, omic2_model_ready_train, omic3_model_ready_train)

     } else {

         geno_omic1_omic2_omic3_train = scale(cbind(cbind(geno_model_ready, omic1_model_ready),
                                                    cbind(omic2_model_ready, omic3_model_ready)))

         rm(geno_model_ready, omic1_model_ready, omic2_model_ready, omic3_model_ready)
     }

 }


 ##############################################################
 ###################################################################
 ###  Start ML Analysis
 ###
 ######################################################################
 #####################################################################

 if(exists('geno_model_ready_test') & exists('geno_model_ready_train')) {

     if(GS_model=="Xgboost"){

     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_model_ready_train,
                                geno_omic_test_object = geno_model_ready_test,
                                para_tunning = para_tunning
                                )

     }

     if(GS_model=="RandomForest"){

         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_model_ready_train,
                                    geno_omic_test_object = geno_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){

         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_model_ready_train,
                                             geno_omic_test_object = geno_model_ready_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){

         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_model_ready_train,
                                    geno_omic_test_object = geno_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     ##### RR and Lasso
     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){

         res_model_output <- AI_RidgeRegression_Lasso(
             pheno_object = pheno_clean,
             response = response,
             geno_omic_object = geno_model_ready_train,
             geno_omic_test_object = geno_model_ready_test,
             para_tunning = para_tunning,
             GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                       pheno_object= pheno_clean,
                                                       response = response,
                                                       test_set = test_set_,
                                                       geno_model_ready_train = geno_model_ready_train,
                                                       eval_metrics = eval_metrics,
                                               GS_model=GS_model
                                                       )

     res_plot <- plot_acc_AI(mod=res_model_output,
                         pheno_object= pheno_clean,
                         response = response,
                         test_set = test_set_,
                         GS_model=GS_model)



 } else {


     if(!exists('geno_model_ready_test') & exists('geno_model_ready_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_model_ready_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_model_ready_train,
                                                 para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_model_ready_train,
                                        para_tunning = para_tunning
             )

         }


         ####
         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_model_ready_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = geno_model_ready_train,
                                                  eval_metrics = eval_metrics,
                                                  GS_model=GS_model
         )
        #} ## end of Xgboost

     }
 }

 if(exists('omic1_model_ready_test') & exists('omic1_model_ready_train'))  {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic1_model_ready_train,
                                geno_omic_test_object = omic1_model_ready_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_model_ready_train,
                                    geno_omic_test_object = omic1_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic1_model_ready_train,
                                             geno_omic_test_object = omic1_model_ready_test,
                                             para_tunning = para_tunning
         )

     }


     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_model_ready_train,
                                    geno_omic_test_object = omic1_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     #### RR and Lasso
     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_model_ready_train,
                                    geno_omic_test_object = omic1_model_ready_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                       pheno_object= pheno_clean,
                                                       response = response,
                                                       test_set = test_set_,
                                               geno_omic_object = omic1_model_ready_train,
                                                       eval_metrics = eval_metrics,
                                               GS_model=GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=GS_model)

     #} ## End of Xgboost



 } else {

     if(!exists('omic1_model_ready_test') & exists('omic1_model_ready_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_model_ready_train,
                                    para_tunning = para_tunning
         )

         }

         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = omic1_model_ready_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         ### RR and Lasso

         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_model_ready_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                           pheno_object= pheno_clean,
                                                           response = response,
                                                   geno_omic_object = omic1_model_ready_train,
                                                           eval_metrics = eval_metrics,
                                                   GS_model=GS_model
         )

     #} ## End Xgboost


     }

 }

 if(exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic2_model_ready_train,
                                geno_omic_test_object = omic2_model_ready_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_model_ready_train,
                                    geno_omic_test_object = omic2_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic2_model_ready_train,
                                             geno_omic_test_object = omic2_model_ready_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_model_ready_train,
                                    geno_omic_test_object = omic2_model_ready_test,
                                    para_tunning = para_tunning
         )

     }


     ### RR and Lasso

     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_model_ready_train,
                                    geno_omic_test_object = omic2_model_ready_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }
     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = omic2_model_ready_train,
                                               eval_metrics = eval_metrics,
                                               GS_model=GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=GS_model)

     #}


 } else {

     if(!exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_model_ready_train,
                                    para_tunning = para_tunning
         )

}

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                         response = response,
                                         geno_omic_object = omic2_model_ready_train,
                                         para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         ### RR and Lasso
         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_model_ready_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

            res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                    pheno_object= pheno_clean,
                                                    response = response,
                                                   geno_omic_object = omic2_model_ready_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model=GS_model
         )

         #}

     }

 }


 if(exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {

     if (GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic3_model_ready_train,
                                geno_omic_test_object = omic3_model_ready_test,
                                para_tunning = para_tunning
     )

     }

     if (GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic3_model_ready_train,
                                    geno_omic_test_object = omic3_model_ready_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic3_model_ready_train,
                                             geno_omic_test_object = omic3_model_ready_test,
                                             para_tunning = para_tunning
         )

     }

     if (GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic3_model_ready_train,
                                    geno_omic_test_object = omic3_model_ready_test,
                                    para_tunning = para_tunning
         )

     }


     ### RR and Lasso

     if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic3_model_ready_train,
                                    geno_omic_test_object = omic3_model_ready_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = omic3_model_ready_train,
                                               eval_metrics = eval_metrics,
                                               GS_model=GS_model
     )


     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=GS_model)

     #}

     ## Strat of ranopdm foes

 } else {

     if(!exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic3_model_ready_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic3_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                         response = response,
                                         geno_omic_object = omic3_model_ready_train,
                                         para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic3_model_ready_train,
                                        para_tunning = para_tunning
             )

         }

         ## RR and Lasso

         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic3_model_ready_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }


         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = omic3_model_ready_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model=GS_model
         )


         }



 }


 if(exists('geno_omic1_test') & exists('geno_omic1_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic1_train,
                                geno_omic_test_object = geno_omic1_test,
                                para_tunning = para_tunning
     )

     }


     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_train,
                                    geno_omic_test_object = geno_omic1_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_omic1_train,
                                             geno_omic_test_object = geno_omic1_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_train,
                                    geno_omic_test_object = geno_omic1_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_train,
                                    geno_omic_test_object = geno_omic1_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic1_train,
                                               eval_metrics = eval_metrics,
                                               GS_model=GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=GS_model)

     #}



 } else {

     if(!exists('geno_omic1_test') & exists('geno_omic1_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_train,
                                        para_tunning = para_tunning
             )

         }

         ### RR and lasso

         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = geno_omic1_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model = GS_model
         )


         #}


     }


 }


 if(exists('geno_omic2_test') & exists('geno_omic2_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic2_train,
                                geno_omic_test_object = geno_omic2_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_train,
                                    geno_omic_test_object = geno_omic2_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_train,
                                    geno_omic_test_object = geno_omic2_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_train,
                                    geno_omic_test_object = geno_omic2_test,
                                    para_tunning = para_tunning
         )

     }


     ## RR and lasso

     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_train,
                                    geno_omic_test_object = geno_omic2_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }


     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic2_train,
                                               eval_metrics = eval_metrics,
                                               GS_model=GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=GS_model)

     #}


 } else {

     if(!exists('geno_omic2_test') & exists('geno_omic2_train')) {

         if (GS_model=="Xgboost"){

         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_train,
                                    para_tunning = para_tunning
         )

         }


         if (GS_model=="RandomForest"){

             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){

             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){

             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_train,
                                        para_tunning = para_tunning
             )

         }


         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){

             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                 pheno_object= pheno_clean,
                                                 response = response,
                                                  geno_omic_object = geno_omic2_train,
                                                  eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )


         #}


     }


 }


 if(exists('geno_omic3_test') & exists('geno_omic3_train')) {

     if(GS_model=="Xgboost"){

     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic3_train,
                                geno_omic_test_object = geno_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){

         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic3_train,
                                    geno_omic_test_object = geno_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){

         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic3_train,
                                    geno_omic_test_object = geno_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){

         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic3_train,
                                    geno_omic_test_object = geno_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     ### RR and lasso
     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){

         res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic3_train,
                                    geno_omic_test_object = geno_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic3_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}

     ## Strat rf


 } else {

     if(!exists('geno_omic3_test') & exists('geno_omic3_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic3_train,
                                    para_tunning = para_tunning
         )

         }


         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_omic3_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         ##
         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = geno_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model=GS_model
         )

         #}

     }


 }


 if(exists('omic1_omic2_test') & exists('omic1_omic2_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic1_omic2_train,
                                geno_omic_test_object = omic1_omic2_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_train,
                                    geno_omic_test_object = omic1_omic2_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic1_omic2_train,
                                             geno_omic_test_object = omic1_omic2_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_train,
                                    geno_omic_test_object = omic1_omic2_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_train,
                                    geno_omic_test_object = omic1_omic2_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = omic1_omic2_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}


 } else {

     if(!exists('omic1_omic2_test') & exists('omic1_omic2_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_train,
                                        para_tunning = para_tunning
             )

         }


         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = omic1_omic2_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )


         #}


     }


 }


 if(exists('omic1_omic3_test') & exists('omic1_omic3_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic1_omic3_train,
                                geno_omic_test_object = omic1_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic3_train,
                                    geno_omic_test_object = omic1_omic3_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic1_omic3_train,
                                             geno_omic_test_object = omic1_omic3_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic3_train,
                                    geno_omic_test_object = omic1_omic3_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic3_train,
                                    geno_omic_test_object = omic1_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                            pheno_object= pheno_clean,
                                            response = response,
                                            test_set = test_set_,
                                            geno_omic_object = omic1_omic3_train,
                                            eval_metrics = eval_metrics,
                                            GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}


 } else {

     if(!exists('omic1_omic3_test') & exists('omic1_omic3_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = omic1_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model = GS_model
         )


         #}

     }


 }


 if(exists('omic2_omic3_test') & exists('omic2_omic3_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic2_omic3_train,
                                geno_omic_test_object = omic2_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_omic3_train,
                                    geno_omic_test_object = omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_omic3_train,
                                    geno_omic_test_object = omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_omic3_train,
                                    geno_omic_test_object = omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                   pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_omic3_train,
                                    geno_omic_test_object = omic2_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }
     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                              pheno_object= pheno_clean,
                                              response = response,
                                              test_set = test_set_,
                                              geno_omic_object = omic2_omic3_train,
                                              eval_metrics = eval_metrics,
                                               GS_model = GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model=  GS_model)

     #}



 } else {

     if(!exists('omic2_omic3_test') & exists('omic2_omic3_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic2_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic2_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }
         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = omic2_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )


         #}


     }


 }


 if(exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic1_omic2_train,
                                geno_omic_test_object = geno_omic1_omic2_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_train,
                                    geno_omic_test_object = geno_omic1_omic2_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_omic1_omic2_train,
                                             geno_omic_test_object = geno_omic1_omic2_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_train,
                                    geno_omic_test_object = geno_omic1_omic2_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_train,
                                    geno_omic_test_object = geno_omic1_omic2_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic1_omic2_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}

     ## strat rf


 } else {

     if(!exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {

         if(GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_train,
                                    para_tunning = para_tunning
         )

         }

         if(GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_train,
                                        para_tunning = para_tunning
             )

         }

         if(GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_omic1_omic2_train,
                                                 para_tunning = para_tunning
             )

         }

         if(GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_train,
                                        para_tunning = para_tunning
             )

         }


         if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }
         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                  pheno_object= pheno_clean,
                                                  response = response,
                                                  geno_omic_object = geno_omic1_omic2_train,
                                                  eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )


         #}

     }


 }


 if(exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {

     if (GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic1_omic3_train,
                                geno_omic_test_object = geno_omic1_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if (GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_omic1_omic3_train,
                                             geno_omic_test_object = geno_omic1_omic3_test,
                                             para_tunning = para_tunning
         )

     }

     if (GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic1_omic3_train,
                                               eval_metrics = eval_metrics,
                                               GS_model = GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}


 } else {

     if(!exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_omic1_omic3_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }
         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                  pheno_object= pheno_clean,
                                                  response = response,
                                                  geno_omic_object = geno_omic1_omic3_train,
                                                  eval_metrics = eval_metrics,
                                                  GS_model= GS_model
         )


         #}

     }


 }

 ######
 if(exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {

     if(GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic2_omic3_train,
                                geno_omic_test_object = geno_omic2_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if(GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if(GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_omic2_omic3_train,
                                             geno_omic_test_object = geno_omic2_omic3_test,
                                             para_tunning = para_tunning
         )

     }

     if(GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }


     if(GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic2_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic2_omic3_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)


     #}

 } else {

     if(!exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic2_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_omic2_omic3_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic2_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = geno_omic2_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )

         #}


     }


 }

 #####

 if(exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {

     if (GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = omic1_omic2_omic3_train,
                                geno_omic_test_object = omic1_omic2_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if (GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_omic3_train,
                                    geno_omic_test_object = omic1_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = omic1_omic2_omic3_train,
                                             geno_omic_test_object = omic1_omic2_omic3_test,
                                             para_tunning = para_tunning
         )

     }


     if (GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_omic3_train,
                                    geno_omic_test_object = omic1_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_omic3_train,
                                    geno_omic_test_object = omic1_omic2_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = omic1_omic2_omic3_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}


 } else {

     if(!exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = omic1_omic2_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = omic1_omic2_omic3_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }

         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = omic1_omic2_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }

         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = omic1_omic2_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )


         #}

     }


 }

 #####
 if(exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {

     if (GS_model=="Xgboost"){
     res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                response = response,
                                geno_omic_object = geno_omic1_omic2_omic3_train,
                                geno_omic_test_object = geno_omic1_omic2_omic3_test,
                                para_tunning = para_tunning
     )

     }

     if (GS_model=="RandomForest"){
         res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }

     if (GS_model=="K-NearestNeighbors"){
         res_model_output <- AI_knn(pheno_object = pheno_clean,
                                             response = response,
                                             geno_omic_object = geno_omic1_omic2_omic3_train,
                                             geno_omic_test_object = geno_omic1_omic2_omic3_test,
                                             para_tunning = para_tunning
         )

     }

     if (GS_model=="SupportVectorMachine"){
         res_model_output <- AI_svm(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic2_omic3_test,
                                    para_tunning = para_tunning
         )

     }


     if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
         res_model_output <- AI_RidgeRegression_Lasso(
                                    pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_omic3_train,
                                    geno_omic_test_object = geno_omic1_omic2_omic3_test,
                                    para_tunning = para_tunning,
                                    GS_model = GS_model
         )

     }

     res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                               pheno_object= pheno_clean,
                                               response = response,
                                               test_set = test_set_,
                                               geno_omic_object = geno_omic1_omic2_omic3_train,
                                               eval_metrics = eval_metrics,
                                               GS_model= GS_model
     )

     res_plot <- plot_acc_AI(mod=res_model_output,
                          pheno_object= pheno_clean,
                          response = response,
                          test_set = test_set_,
                          GS_model= GS_model)

     #}


 } else {

     if(!exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {

         if (GS_model=="Xgboost"){
         res_model_output <- AI_Xgb(pheno_object = pheno_clean,
                                    response = response,
                                    geno_omic_object = geno_omic1_omic2_omic3_train,
                                    para_tunning = para_tunning
         )

         }

         if (GS_model=="RandomForest"){
             res_model_output <- AI_randomForest(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         if (GS_model=="K-NearestNeighbors"){
             res_model_output <- AI_knn(pheno_object = pheno_clean,
                                                 response = response,
                                                 geno_omic_object = geno_omic1_omic2_omic3_train,
                                                 para_tunning = para_tunning
             )

         }

         if (GS_model=="SupportVectorMachine"){
             res_model_output <- AI_svm(pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_omic3_train,
                                        para_tunning = para_tunning
             )

         }


         if (GS_model=="Lasso" | GS_model=="Ridge_Regression"){
             res_model_output <- AI_RidgeRegression_Lasso(
                                        pheno_object = pheno_clean,
                                        response = response,
                                        geno_omic_object = geno_omic1_omic2_omic3_train,
                                        para_tunning = para_tunning,
                                        GS_model = GS_model
             )

         }


         res_summary_stat <- summary_statistics_AI(mod=res_model_output,
                                                   pheno_object= pheno_clean,
                                                   response = response,
                                                   geno_omic_object = geno_omic1_omic2_omic3_train,
                                                   eval_metrics = eval_metrics,
                                                   GS_model= GS_model
         )

#}


     }


 }



}  ### End machine learning

 }

 ### if user provide only

 #if(GS_model)
 if (GS_model=="Xgboost"){
     if(exists("test_set_")){
 output <- list(res_model_output, res_summary_stat, res_plot)

 names(output) <- c('model_results', 'summary_statistic', 'res_plot')

     } else {

         output <- list(res_model_output, res_summary_stat)

         names(output) <- c('model_results', 'summary_statistic')
     }

 } else if (GS_model=="GBLUP"){

     output <- list(res_model_output)

     names(output) <- c('model_results')

} else {

     output <- list(res_model_output, res_summary_stat)

     names(output) <- c('model_results', 'summary_statistic')

 }

 return(output)

} ## end of function
