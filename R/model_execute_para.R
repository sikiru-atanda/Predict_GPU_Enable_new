#' Title
#'
#' @param pheno_data
#' @param pheno_train
#' @param pheno_test
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param gmatrix
#' @param train_geno_data
#' @param train_omic1_data
#' @param train_omic2_data
#' @param train_omic3_data
#' @param test_geno_data
#' @param test_omic1_data
#' @param test_omic2_data
#' @param test_omic3_data
#' @param train_coefficient
#' @param train_set
#' @param test_set
#' @param gmatrix_method
#' @param response
#' @param gen_name
#' @param cova
#' @param fixed
#' @param random
#' @param heter_resid
#' @param heter_groups
#' @param varcov_str
#' @param weights
#' @param nIter
#' @param burnIn
#' @param thin
#' @param GS_model
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param Cross_validation
#' @param test_size
#' @param random_state
#' @param replication
#' @param nFolds
#' @param CV
#' @param core
#' @param message
#' @param ...
#' @param gkernel
#' @param kernel_method
#' @param center
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#'
#' @return
#' @export
#'
#' @examples
model_execute <- function(
    pheno_data=NULL,
    pheno_train = NULL,
    pheno_test = NULL,
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
    train_coefficient = NULL,
    train_set = NULL,
    test_set = NULL,
    gmatrix_method = NULL,
    # gmatrix_method = c("VanRaden",
    #                    "Yang"),
    # kernel_method = c("Gaussian",
    #                          "Linear",
    #                          "Poly2",
    #                          "Poly3",
    #                          "Poly4"),
    kernel_method = NULL,
    response=NULL,
    gen_name=NULL,
    cova=NULL,
    fixed=NULL,
    random=NULL,
    heter_resid=FALSE,
    heter_groups=NULL,
    varcov_str = NULL,
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
                 "BayesL"),
    fixed_term_model_bayesian = 'FIXED',
    rand_term_model_bayesian = NULL,
    Cross_validation = c("Hold_Out",
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
    ### Get clean pheno data for model fit


 pheno_clean <- phenotype_to_model(
                    pheno = pheno_data,
                    pheno_train = pheno_train,
                    pheno_test = pheno_test,
                    train_set = train_set,
                    test_set = test_set,
                    response = response,
                    gen_name = gen_name)

 ### pheno_clean is a list with three elements.
 ## First element is pheno_data
 ## Second element is test_set
 ## Third element is train_set

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

     if((exists("geno_clean") & exists("pheno_clean"))){
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

     rm(geno_pheno_match)
     }

     #### To calculate the gkernel only geno_clean is acceptable
     if (exists("geno_clean") & !is.null(kernel_method)) {
         gkernel <- kernel_calculation(
             M_matrix_clean = geno_clean,
             center=center,
             method = kernel_method,
             message = message )
     }

     #### To calculate the gmatrix only geno_clean is acceptable
     if (exists("geno_clean") & !is.null(kernel_method)) {
        gmatrix <- grm_calculation(
        geno_clean = geno_clean,
        method=gmatrix_method)
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

     if(exists('omic1_clean') & exists("pheno_clean")){
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

         rm(omic1_pheno_match)
     }

     #### To calculate the omic1_kernel only omic1_clean is acceptable
     if (exists("omic1_clean") & !is.null(kernel_method)) {
         omic1_kernel <-  kernel_calculation(
             M_matrix_clean = omic1_clean,
             center=center,
             method = kernel_method,
             message = message )
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


     if(exists('omic2_clean') & exists("pheno_clean")){
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

         rm(omic2_pheno_match)
     }

     #### To calculate the omic2_kernel only omic2_clean is acceptable
     if (exists("omic2_clean") & !is.null(kernel_method)) {
         omic2_kernel <-  kernel_calculation(
             M_matrix_clean = omic2_clean,
             center=center,
             method = kernel_method,
             message = message )
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


     if(exists('omic3_clean') & exists("pheno_clean")){
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

         rm(omic3_pheno_match)
     }

     #### To calculate the omic3_kernel only omic3_clean is acceptable
     if (exists("omic3_clean") & !is.null(kernel_method)) {
         omic3_kernel <-  kernel_calculation(
             M_matrix_clean = omic3_clean,
             center=center,
             method = kernel_method,
             message = message )
     }

 }
###################################################################################
 #  Pre-Check for grm/kernel matrix if calculated from the marker/omic data
 # or provided by the user.
 # It has to pass through this pre-check before going to
 # conditioning effect such as bend or blending. The conditioning of the grm/kernel
 # matrix is important especially the bend but we going to give user the opportunity
 # to decide to do it or not.
######################################################################################3

     if(!is.null(gmatrix)){

         gmatrix_checked <- grm_kernel_precheck(object= gmatrix,
                                                message= message)
     }


     if(!is.null(gkernel)){

         gkernel_checked <- grm_kernel_precheck(object= gkernel,
                                                message= message)
     }



     if(!is.null(gmatrix)){

         gmatrix_checked <- grm_kernel_precheck(object= gmatrix,
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


 }


 ### End
#####################################################################



 #### TO DO put a condition to check the user is providing only one location
 ### If the user provide only the marker matrix


    ### Model BRR for single location
 #if(((GS_model=="BRR") & is.null(rand_term_model_bayesian)) || ((is.null(GS_model) & (rand_term_model_bayesian=="BRR")))){

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
        omic1_data = NULL,
        omic2_data = NULL,
        omic3_data = NULL,
        gen_name = gen_name)

   bayes_para <-  bayes_parameter_check(nIter = nIter,
                                        burnIn = burnIn,
                                        thin = thin)

    mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
                                          response = response,
                                          weights = weights,
                                          ETA = ETA$ETA,
                                          bayes_para = bayes_para,
                                          verbose = FALSE,
                                          #files_key = "files_key"
                                          )

    res_model_output <- mod_output_bayes(mod = mod,
                                         ETA = ETA,
                                         geno_data = geno_model_ready,
                                         gen_name = gen_name,
                                         omic1_data = NULL,
                                         omic2_data = NULL,
                                         omic3_data = NULL)

    res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = NULL
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


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

            mod <-  M_matrix_bayes_mod_single_loc(object = ETA$pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


        }


        ##### omic1_clean, omic2clean and omic3_clean

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

            mod <-  M_matrix_bayes_mod_single_loc(object = pheno_data,
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
                                                 omic3_data = omic3_model_ready
            )

            res_summary_stat <- summary_statistics_bayes(mod = mod)


        }
}

 ### if user provide only

 output <- list(res_model_output, res_summary_stat)

 names(output) <- c('model_results', 'summary_statistic')

 return(output)

} ## end of function
