
#' Title
#'
#' @param train_set
#' @param test_set
#' @param response
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param weights
#' @param nIter
#' @param burnIn
#' @param thin
#' @param GS_model
#' @param cross_validation
#' @param test_size
#' @param random_state
#' @param nFolds
#' @param CV
#' @param core
#' @param ...
#' @param pheno_data
#' @param pheno_data_train
#' @param pheno_data_test
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param gmatrix
#' @param gkernel
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param train_geno_data
#' @param train_omic1_data
#' @param train_omic2_data
#' @param train_omic3_data
#' @param test_geno_data
#' @param test_omic1_data
#' @param test_omic2_data
#' @param test_omic3_data
#' @param gmatrix_method
#' @param kernel_method
#' @param gen_name
#' @param var_cov_str
#' @param eval_metrics
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param replication
#' @param message
#' @param center
#'
#' @return
#' @export
#'
#' @examples
#'
CrossVal <- function(
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
    eval_metrics = c("Accuracy",
                     "Mean_Squared_Error",
                     "Bias",
                     "Root_Mean_Squared_Error",
                     "Relative_Squared_Error",
                     "Mean_Absolute_Error",
                     "Mean_Absolute_Percent_Error"),
    fixed_term_model_bayesian = 'FIXED',
    rand_term_model_bayesian = NULL,
    cross_validation = c("Hold_Out",
                         "Stratified_Hold_Out",
                         "Repeated_Hold_Out",
                         "Repeated_Stratified_Hold_Out",
                         "K-Folds",
                         "Stratified_K-Folds",
                         "Repeated_K-Folds",
                         "Repeated_Stratified_K-Folds",
                         "Leave_one_Out",
                         "CV1",
                         "Repeated_CV1",
                         "CV2",
                         "Repeated_CV2"),
    test_size = 0.5,
    random_state = NULL,
    replication = 1,
    nFolds = NULL,
    CV = NULL,
    core = NULL,
    message = TRUE,
    center = TRUE,
    ...
) {

  msg <- sprintf("==================================================\n")


    pheno_clean <- phenotype_to_model(
      pheno_data = pheno_data,
      pheno_data_train = pheno_data_train,
      pheno_data_test = pheno_data_test,
      response = response,
      gen_name = gen_name)

    ### pheno_clean is a list with three elements.
    ## First element is pheno_data
    ## Second element is test_set if user provide pheno_data_test as an input
    ## Third element is train_set if user provide pheno_data_train as an input.

    if(attr(pheno_clean[[1]], "cleared")!="for_model_fit" && all(class(pheno_clean[[1]])!=c("data.frame", "phenotype"))) {

      stop('pheno_data is not phenotype data')
    }
#############################################
## This cater for condition where user provide the whole pheno data
## and provide ID of individuals in the training set for the purpose of
## cross-validation.
## NB this set must have phenotypic and genotypic records.
    if (!is.null(train_set) & is.null(test_set)){

      pheno_clean <- pheno_clean[(as.character(pheno_clean[, gen_name]))%in%train_set[,1], ]

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

    #### Get the clean omic1_data ready for model fit
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
    # or provided by the user.  # It has to pass through this pre-check before going
    # to conditioning effect such as bend or blending.
    # The conditioning of the grm/kernel matrix is important especially the bend
    # but we going to give user the opportunity to decide to do it or not.
    ######################################################################################3

    if(!is.null(gmatrix)){

      gmatrix_checked <- grm_kernel_precheck(object= gmatrix,
                                             message= message)
    }


    if(!is.null(gkernel)){

      gkernel_checked <- grm_kernel_precheck(object= gkernel,
                                             message= message)
    }


    if(!is.null(omic1_kernel)){

      omic1_kernel_checked <- grm_kernel_precheck(object= omic1_kernel,
                                                  message= message)
    }

    if(!is.null(omic2_kernel)){

      omic2_kernel_checked <- grm_kernel_precheck(object= omic2_kernel,
                                                  message= message)
    }


    if(!is.null(omic3_kernel)){

      omic3_kernel_checked <- grm_kernel_precheck(object= omic3_kernel,
                                                  message= message)
    }
    ################################################################
    ##### Pheno to geno match
    ################################################

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

 ################################


    if (cross_validation=="Stratified_Hold_Out") {
      test_set_val <- train_test_split(
        pheno_data = pheno_clean$pheno_data,
        gen_name = gen_name,
        response = response,
        test_size =  test_size,
        random_state = random_state,
        replication = 1,
        method = "stratified"
      )
    }
    if (cross_validation=="Repeated_Stratified_Hold_Out") {
      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}
      test_set_val <- train_test_split(
        pheno_data = pheno_clean$pheno_data,
        gen_name = gen_name,
        response = response,
        test_size =  test_size,
        random_state = random_state,
        replication = replication,
        method = "stratified"
      )
    }
    if (cross_validation=="Hold_Out") {
      test_set_val <- train_test_split(
        pheno_data = pheno_clean$pheno_data,
        gen_name = gen_name,
        response = response,
        test_size =  test_size,
        random_state = random_state,
        replication = 1,
        method = "unstratified"
      )

    }

    if (cross_validation=="Repeated_Hold_Out") {
      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}
      test_set_val <- train_test_split(
        pheno_data = pheno_clean$pheno_data,
        gen_name = gen_name,
        response = response,
        test_size =  test_size,
        random_state = random_state,
        replication = replication,
        method = "unstratified"
      )

    }




    if (cross_validation=="Stratified_K-Folds") {

      test_set_val <- CV_nfolds(pheno_data = pheno_clean$pheno_data,
                     response = response,
                     gen_name = gen_name,
                     nFolds = nFolds,
                     random_state = random_state,
                     replication = 1,
                     method = 'stratified')
    }

    if (cross_validation=="Repeated_Stratified_K-Folds") {
      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}

      test_set_val <- CV_nfolds(pheno_data = pheno_clean$pheno_data,
                            response = response,
                            gen_name = gen_name,
                            nFolds = nFolds,
                            random_state = random_state,
                            replication = replication,
                            method = 'stratified')
    }

    if (cross_validation=="K-Folds") {

      test_set_val <- CV_nfolds(pheno_data = pheno_clean$pheno_data,
                            response = response,
                            gen_name = gen_name,
                            nFolds = nFolds,
                            random_state = random_state,
                            replication = 1,
                            method = 'unstratified')
    }

    if (cross_validation=="Repeated_K-Folds") {
      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}

      test_set_val <- CV_nfolds(pheno_data = pheno_clean$pheno_data,
                            response = response,
                            gen_name = gen_name,
                            nFolds = nFolds,
                            random_state = random_state,
                            replication = replication,
                            method = 'unstratified')
    }

    if (cross_validation=="CV1" & !is.null(CV)) {
      if(CV!=1 | is.null(CV)){stop(print(paste(msg, "CV must be 1")), call. = FALSE)}

      test_set_val <-CV1_CV2(pheno_data = pheno_clean$pheno_data,
                     gen_name = gen_name,
                     heter_groups = heter_groups,
                     CV = CV,
                     nFolds = nFolds,
                     random_state = random_state,
                     replication = 1)
    }

    if (cross_validation=="Repeated_CV1" & !is.null(CV)) {
      if(CV!=1 | is.null(CV)){stop(print(paste(msg, "CV must be 1")), call. = FALSE)}

      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}

      test_set_val <-CV1_CV2(pheno_data = pheno_clean$pheno_data,
                     gen_name = gen_name,
                     heter_groups = heter_groups,
                     CV = CV,
                     nFolds = nFolds,
                     random_state = random_state,
                     replication = replication)
    }

    if (cross_validation=="CV2" & !is.null(CV)) {

      if(CV!=2 | is.null(CV)){stop(print(paste(msg, "CV must be 2")), call. = FALSE)}

      test_set_val <-CV1_CV2(pheno_data = pheno_clean$pheno_data,
                     gen_name = gen_name,
                     heter_groups = heter_groups,
                     CV = CV,
                     nFolds = nFolds,
                     random_state = random_state,
                     replication = 1)
    }

    if (cross_validation=="Repeated_CV2" & !is.null(CV)) {

      if(CV!=2 | is.null(CV)){stop(print(paste(msg, "CV must be 2")), call. = FALSE)}

      if(replication==1) { warning(paste(msg,'Replication must be greater than 1.'),
                                   call. = FALSE)}

      test_set_val <-CV1_CV2(pheno_data = pheno_clean$pheno_data,
                     gen_name = gen_name,
                     heter_groups = heter_groups,
                     CV = CV,
                     nFolds = nFolds,
                     random_state = random_state,
                     replication = replication)
    }

    if (cross_validation=="Leave_one_Out") {

      # Test <-CV1_CV2(pheno = pheno,
      #                 genotype = genotype,
      #                 heter_groups = heter_groups,
      #                 CV = CV,
      #                 nFolds = nFolds,
      #                 random_state = random_state,
      #                 Replication = 1)
    }
    ###
    if(is.null(GS_model)) {stop(print(paste(msg,"Provide genomic selection model from the list avaliable")), call. = FALSE)}

    if((isTRUE(GS_model== "BRR" | GS_model== "BayesA"|  GS_model== "BayesB"| GS_model== "BayesC" | GS_model== "BL") & is.null(rand_term_model_bayesian)) |
       ((is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL"))))) |
       ((!is.null(GS_model) & isTRUE(all(rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL")))))){

      if((exists('geno_model_ready') & (!exists('omic1_model_ready') & (!exists('omic2_model_ready') & !exists('omic3_model_ready'))))){

        ETA  <-  ETA_compiler_bayes(
          fixed = fixed,
          random = random,
          GS_model = GS_model,
          fixed_term_model_bayesian = fixed_term_model_bayesian,
          rand_term_model_bayesian = rand_term_model_bayesian,
          pheno_data = pheno_clean[[1]],
          geno_data = geno_model_ready,
          gen_name = gen_name)

        bayes_para <-  bayes_parameter_check(nIter = nIter,
                                             burnIn = burnIn,
                                             thin = thin)

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


        #res_summary_stat <- summary_statistics_bayes(mod = mod)


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
        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )




        #res_summary_stat <- summary_statistics_bayes(mod = mod)


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

      }


    } ## End of  Bayes A, B, C, BRR, BL

  ####
    ##########################################################################
    #########################################################################
    ## Start of RKHS  Model for Single Location                            ##                          ##
    ##                                                                     ##
    ##                                                                     ##
    ##########################################################################
    #######################################################################

    # if((isTRUE(GS_model== "RKHS") |
    #     ((is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%"RKHS"))) |
    #     ((!is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%"RKHS"))))){

    if((isTRUE(GS_model== "RKHS" | isTRUE(GS_model== "BRR")) & is.null(rand_term_model_bayesian)) |
       ((is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%c("RKHS", "BRR")))) |
       ((!is.null(GS_model) & isTRUE(rand_term_model_bayesian%in%c("RKHS", "BRR"))))){


      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )



      }
      ###### omic1_clean

      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }

      ##### omic2_model_ready

      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }

      ####

      ##### omic3_model_ready

      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }


      ##### geno_model_ready and omic1_model_ready

      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }


      ##### geno_clean and omic2_clean

      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }
      ##### geno_clean and omic3_clean

      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }


      ##### omic1_model_ready and omic2_model_ready

      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }



      ##### omic1_model_ready and omic3_model_ready


      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }


      ####
      ##### omic2_model_ready and omic3_model_ready


      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }

      ########

      ##### geno_model_ready, omic1_model_ready and omic2_model_ready

      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & !exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }



      ##### geno_model_ready, omic1_model_ready and omic3_model_ready


      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (!exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }

      #####

      ##### geno_model_ready, omic2_model_ready and omic3_model_ready


      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (!exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }


      ##### omic1_model_ready, omic2_model_ready and omic3_model_ready

      if(((!exists('gkernel_model_ready') & !exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )


      }



      ##### geno, omic1_clean, omic2clean and omic3_clean

      if(((exists('gkernel_model_ready') | exists('gmatrix_model_ready')) & (exists('omic1_kernel_model_ready') & (exists('omic2_kernel_model_ready') & exists('omic3_kernel_model_ready'))))){

        if(exists('gkernel_model_ready') & exists('gmatrix_model_ready')){ stop(paste(msg, 'Either gmatrix or gkernel is expected not both at the same time.'))}

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

        mod <-  bayes_mod_execute_crossval(object = ETA$pheno_data,
                                           response = response,
                                           weights = weights,
                                           ETA = ETA$ETA,
                                           bayes_para = bayes_para,
                                           test_set_val = test_set_val,
                                           heter_groups = heter_groups,
                                           core = core,
                                           cross_validation = cross_validation,
                                           eval_metrics =  eval_metrics
        )

      }


    } #### END GBLUP_RKHS


    ### if user provide only

    output <- list(res_model_output, res_summary_stat)

    names(output) <- c('model_results', 'summary_statistic')

    return(output)

} ## end of function



