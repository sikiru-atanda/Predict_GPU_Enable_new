
#' Title
#'
#' @param pheno_data
#' @param pheno_data_train
#' @param pheno_data_test
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param train_geno_data
#' @param train_omic1_data
#' @param train_omic2_data
#' @param train_omic3_data
#' @param test_geno_data
#' @param test_omic1_data
#' @param test_omic2_data
#' @param test_omic3_data
#' @param train_set
#' @param test_set
#' @param response
#' @param gen_name
#' @param core
#' @param message
#' @param center
#' @param xgb_paras_tunning
#' @param ...
#' @param learning_rate
#' @param max_depth
#' @param subsample
#' @param booster
#' @param iteration
#' @param gmatrix
#' @param gkernel
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param gmatrix_method
#' @param kernel_method
#' @param para_tunning
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom foreach %dopar%
AI_xgbOLD <- function(pheno_data=NULL,
                   pheno_data_train = NULL,
                   pheno_data_test = NULL,
                   geno_data = NULL,
                   omic1_data = NULL,
                   omic2_data = NULL,
                   omic3_data = NULL,
                   train_geno_data = NULL,
                   train_omic1_data = NULL,
                   train_omic2_data = NULL,
                   train_omic3_data = NULL,
                   test_geno_data = NULL,
                   test_omic1_data = NULL,
                   test_omic2_data = NULL,
                   test_omic3_data = NULL,
                   gmatrix= NULL,
                   gkernel = NULL,
                   omic1_kernel = NULL,
                   omic2_kernel = NULL,
                   omic3_kernel = NULL,
                   gmatrix_method = NULL,
                   kernel_method = NULL,
                   train_set = NULL,
                   test_set = NULL,
                   response=NULL,
                   gen_name=NULL,
                   core = NULL,
                   message = TRUE,
                   center = TRUE,
                   para_tunning = FALSE,
                   xgb_paras_tunning= c(Iter_tune = NULL, # number of boosting iterations
                                        learning_rate_tune = NULL, # learning rate, low value means model is more robust to overfitting
                                        L2_tune = NULL, # L2 Regularization (Ridge Regression)
                                        L1_tune = NULL),
                   learning_rate = 0.001,
                   max_depth = 6,
                   subsample = 0.5,
                   booster = "gblinear",
                   iteration = 5000,
                   ...

){

  msg <- sprintf("==================================================\n")
  ### Get clean pheno data for model fit


  pheno_clean <- phenotype_to_model(
    pheno_data = pheno_data,
    pheno_data_train = pheno_data_train,
    pheno_data_test = pheno_data_test,
    train_set = train_set,
    test_set = test_set,
    response = response,
    gen_name = gen_name)

  ### pheno_clean is a list with three elements.
  ## First element is pheno_data
  ## Second element is test_set if user provide it as input
  ## Third element is train_set if user provide it as input.

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

#}
 ##############################################################
  ###################################################################
  ###  Start ML Analysis
  ###
 ######################################################################
  #####################################################################


# #### Initializing parallel
if(length(response)>1){
  if (is.null(core)){
    cl = parallel::detectCores()

    if (cl> 4){
      # Try in parallel
      cl <- parallel::makeCluster(4)
    } else{
      cl <- parallel::makeCluster(2)
    }

  } else {
    if(!is.null(core)){

      cl <- parallel::makeCluster(core)
    }
  }

  doParallel::registerDoParallel(cl)



  Univariate <- foreach::foreach(trait = 1:length(response),
                                 .errorhandling='pass') %dopar% {

  if(isTRUE(para_tunning)){
    xgb_grid = expand.grid(nrounds = para_tunning$Iter_tune , # number of boosting iterations
                           eta = para_tunning$learning_rate_tune, # learning rate, low value means model is more robust to overfitting
                           lambda = para_tunning$L2_tune, # L2 Regularization (Ridge Regression)
                           alpha = para_tunning$L1_tune # L1 Regularization (Lasso Regression)
    )



    #here we do one better then a validation set, we use cross validation to
    #expand the amount of info we have!

    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    xgb_trcontrol = caret::trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)




    if(exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


      xgb_fit = caret::train(x = geno_model_ready_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_model_ready_train)


    } else {

      if(!exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


        xgb_fit = caret::train(x = geno_model_ready_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_model_ready_train)
      }


    }

    #####
    if((exists('omic1_model_ready_test') & exists('omic1_model_ready_train'))) {

      xgb_fit = caret::train(x = omic1_model_ready_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic1_model_ready_train)


    } else {

      if(!exists('omic1_model_ready_test') & exists('omic1_model_ready_train')) {


        xgb_fit = caret::train(x = omic1_model_ready_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic1_model_ready_train)
      }


    }

    #####
    if((exists('omic2_model_ready_test') & exists('omic2_model_ready_train'))) {

      xgb_fit = caret::train(x = omic2_model_ready_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic2_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic2_model_ready_train)


    } else {

      if(!exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {


        xgb_fit = caret::train(x = omic2_model_ready_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic2_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic2_model_ready_train)
      }


    }

    ####
    #####
    if((exists('omic3_model_ready_test') & exists('omic3_model_ready_train'))) {

      xgb_fit = caret::train(x = omic3_model_ready_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic3_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic3_model_ready_train)


    } else {

      if(!exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {


        xgb_fit = caret::train(x = omic3_model_ready_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic3_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic3_model_ready_train)

      }


    }


    #####
    if((exists('geno_omic1_test') & exists('geno_omic1_train'))) {

      xgb_fit = caret::train(x = geno_omic1_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic1_train)


    } else {

      if(!exists('geno_omic1_test') & exists('geno_omic1_train')) {


        xgb_fit = caret::train(x = geno_omic1_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic1_train)
      }


    }

    #####
    if((exists('geno_omic2_test') & exists('geno_omic2_train'))) {

      xgb_fit = caret::train(x = geno_omic2_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic2_train)


    } else {

      if(!exists('geno_omic2_test') & exists('geno_omic2_train')) {


        xgb_fit = caret::train(x = geno_omic2_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic2_train)

      }


    }

    #####
    #####
    if((exists('geno_omic3_test') & exists('geno_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic2_train)


    } else {

      if(!exists('geno_omic3_test') & exists('geno_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic3_train)

      }


    }

    #####
    if((exists('omic1_omic2_test') & exists('omic1_omic2_train'))) {

      xgb_fit = caret::train(x = omic1_omic2_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic3_train)


    } else {

      if(!exists('omic1_omic2_test') & exists('omic1_omic2_train')) {


        xgb_fit = caret::train(x = omic1_omic2_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic1_omic2_train)
      }


    }

    #####
    if((exists('omic1_omic3_test') & exists('omic1_omic3_train'))) {

      xgb_fit = caret::train(x = omic1_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic1_omic3_train)


    } else {

      if(!exists('omic1_omic3_test') & exists('omic1_omic3_train')) {


        xgb_fit = caret::train(x = omic1_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic1_omic3_train)
      }


    }

    #####
    if((exists('omic2_omic3_test') & exists('omic2_omic3_train'))) {

      xgb_fit = caret::train(x = omic2_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  as.matrix(omic2_omic3_test),
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic2_omic3_train)


    } else {

      if(!exists('omic2_omic3_test') & exists('omic2_omic3_train')) {


        xgb_fit = caret::train(x = omic2_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic2_omic3_train)

      }


    }

    #####
    if((exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic2_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic1_omic2_train)


    } else {

      if(!exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic2_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic1_omic2_train)

      }


    }

    #####
    if((exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic1_omic3_train)


    } else {

      if(!exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic1_omic3_train)

      }


    }

    #####
    if((exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic2_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic2_omic3_train)


    } else {

      if(!exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic2_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic2_omic3_train)
      }


    }

    #####
    if((exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = omic1_omic2_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic1_omic2_omic3_train)


    } else {

      if(!exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {


        xgb_fit = caret::train(x = omic1_omic2_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic1_omic2_omic3_train)
      }


    }

    #####
    if((exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic2_omic3_train,
                             y = pheno_clean[, response[trait]],
                             trControl = xgb_trcontrol,
                             tuneGrid = xgb_grid,
                             method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_omic1_omic2_omic3_train)


    } else {

      if(!exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic2_omic3_train,
                               y = pheno_clean[, response[trait]],
                               trControl = xgb_trcontrol,
                               tuneGrid = xgb_grid,
                               method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_omic1_omic2_omic3_train)
      }


    }
    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
  } else {

    xgb_params <- list(
      booster = booster,
      eta = learning_rate,
      max_depth = max_depth, #This indicates how deep the built tree can be.
      #The deeper the tree, the more splits it has and it captures more
      #information about how the data. We fit a decision tree with depths
      #ranging from 1 to 32 and plot the training and test errors
      gamma = 4,
      subsample = subsample,
      colsample_bytree = 1,
      objective = "reg:squarederror",
      eval_metric = c("rmse", "rmsle", "mape")
    )





    if(exists('geno_model_ready_test') & exists('geno_model_ready_train')) {

      geno_model_ready_train <- xgboost::xgb.DMatrix(data = geno_model_ready_train,
                                                     label = pheno[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = geno_model_ready_train)



    } else {

      if(!exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


        geno_model_ready_train <- xgboost::xgb.DMatrix(data = geno_model_ready_train,
                                                       label = pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = geno_model_ready_train)


      }


    }

    #####
    if((exists('omic1_model_ready_test') & exists('omic1_model_ready_train'))) {

      omic1_model_ready_train <- xgboost::xgb.DMatrix(data = omic1_model_ready_train,
                                                      label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic1_model_ready_train)


    } else {

      if(!exists('omic1_model_ready_test') & exists('omic1_model_ready_train')) {


        omic1_model_ready_train <- xgboost::xgb.DMatrix(data = omic1_model_ready_train,
                                                        label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic1_model_ready_train)


      }


    }

    #####
    if((exists('omic2_model_ready_test') & exists('omic2_model_ready_train'))) {

      omic2_model_ready_train <- xgboost::xgb.DMatrix(data = omic2_model_ready_train,
                                                      label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic2_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  omic2_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic2_model_ready_train)



    } else {

      if(!exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {


        omic2_model_ready_train <- xgboost::xgb.DMatrix(data = omic2_model_ready_train,
                                                        label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic2_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )

        rm(omic2_model_ready_train)

        xgb_preds <- stats::predict(xgb_fit,
                                    omic2_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic2_model_ready_train)


      }


    }

    ####
    #####
    if((exists('omic3_model_ready_test') & exists('omic3_model_ready_train'))) {

      omic3_model_ready_train <- xgboost::xgb.DMatrix(data = omic3_model_ready_train,
                                                      label =  pheno_clean[, response[trait]])


      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic3_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )

      xgb_preds <- stats::predict(xgb_fit,
                                  omic3_model_ready_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train = omic3_model_ready_train)



    } else {

      if(!exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {


        omic3_model_ready_train <- xgboost::xgb.DMatrix(data = omic3_model_ready_train,
                                                        label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic3_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )




        xgb_preds <- stats::predict(xgb_fit,
                                    omic3_model_ready_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train = omic3_model_ready_train)


      }


    }


    #####
    if((exists('geno_omic1_test') & exists('geno_omic1_train'))) {

      geno_omic1_train <- xgboost::xgb.DMatrix(data = geno_omic1_train,
                                               label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic1_train)



    } else {

      if(!exists('geno_omic1_test') & exists('geno_omic1_train')) {


        geno_omic1_train <- xgboost::xgb.DMatrix(data = geno_omic1_train,
                                                 label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic1_train)


      }


    }

    #####
    if((exists('geno_omic2_test') & exists('geno_omic2_train'))) {

      geno_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic2_train,
                                               label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic2_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic2_train)



    } else {

      if(!exists('geno_omic2_test') & exists('geno_omic2_train')) {


        geno_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic2_train,
                                                 label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic2_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic2_train)


      }


    }

    #####
    #####
    if((exists('geno_omic3_test') & exists('geno_omic3_train'))) {

      geno_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic3_train,
                                               label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic3_train,
        nrounds = iteration,
        verbose = 1
      )

      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic3_train)

    } else {

      if(!exists('geno_omic3_test') & exists('geno_omic3_train')) {


        geno_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic3_train,
                                                 label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic3_train)


      }


    }

    #####
    if((exists('omic1_omic2_test') & exists('omic1_omic2_train'))) {

      omic1_omic2_train <- xgboost::xgb.DMatrix(data = omic1_omic2_train,
                                                label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic2_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  omic1_omic2_train)


    } else {

      if(!exists('omic1_omic2_test') & exists('omic1_omic2_train')) {


        omic1_omic2_train <- xgboost::xgb.DMatrix(data = omic1_omic2_train,
                                                  label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic2_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  omic1_omic2_train)


      }


    }

    #####
    if((exists('omic1_omic3_test') & exists('omic1_omic3_train'))) {

      omic1_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic3_train,
                                                label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  omic1_omic3_train)



    } else {

      if(!exists('omic1_omic3_test') & exists('omic1_omic3_train')) {


        omic1_omic3_train_ <- xgboost::xgb.DMatrix(data = omic1_omic3_train,
                                                   label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic3_train_,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    omic1_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  omic1_omic3_train)


      }


    }

    #####
    if((exists('omic2_omic3_test') & exists('omic2_omic3_train'))) {

      omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic2_omic3_train,
                                                label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                                  omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  omic2_omic3_train)



    } else {

      if(!exists('omic2_omic3_test') & exists('omic2_omic3_train')) {


        omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic2_omic3_train,
                                                  label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  omic2_omic3_train)


      }


    }

    #####
    if((exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train'))) {

      geno_omic1_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_train,
                                                     label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic2_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_omic2_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic1_omic2_train)

    } else {

      if(!exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {


        geno_omic1_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_train,
                                                       label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic2_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic1_omic2_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic1_omic2_train)

      }


    }

    #####
    if((exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train'))) {

      geno_omic1_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic3_train,
                                                     label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic1_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic1_omic3_train)

    } else {

      if(!exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {


        geno_omic1_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic3_train,
                                                       label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic3_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                                    as.matrix(geno_omic1_omic3_train),
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic1_omic3_train)

      }


    }

    #####
    if((exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train'))) {

      geno_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic2_omic3_train,
                                                     label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  geno_omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic2_omic3_train)

    } else {

      if(!exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {


        geno_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic2_omic3_train,
                                                       label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    geno_omic2_omic3_train,
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic2_omic3_train)


      }


    }

    #####
    if((exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train'))) {

      omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic2_omic3_train,
                                                      label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  omic1_omic2_omic3_test,
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  omic1_omic2_omic3_train)

    } else {

      if(!exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {


        omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic2_omic3_train,
                                                        label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    as.matrix(omic1_omic2_omic3_train),
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  omic1_omic2_omic3_train)


      }


    }

    #####
    if((exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train'))) {

      geno_omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_omic3_train,
                                                           label =  pheno_clean[, response[trait]])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                                  as.matrix(geno_omic1_omic2_omic3_test),
                                  reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                               X_train =  geno_omic1_omic2_omic3_train)

    } else {

      if(!exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {


        geno_omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_omic3_train,
                                                             label =  pheno_clean[, response[trait]])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                                    as.matrix(geno_omic1_omic2_omic3_train),
                                    reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                                 X_train =  geno_omic1_omic2_omic3_train)

      }


    }

    #summary(xgb_fit$params)


  } ## End of when no need for tunning.


     model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

output = list(model_para,
              xgb_preds,
              res_feature,
              xgb_fit)


names(output) = c("model_parameters", "predicted_values",
                  "feature_weight", "trained_model")


Univariate = output


    }

  names(Univariate) <- response

  output <-  Univariate

  rm(Univariate)

  ### End when length of response variable is more than 1

   } else {

## when length of response variable is 1
#########################

if(isTRUE(para_tunning)){
xgb_grid = expand.grid(nrounds = xgb_paras_tunning$Iter_tune , # number of boosting iterations
                         eta = xgb_paras_tunning$learning_rate_tune, # learning rate, low value means model is more robust to overfitting
                         lambda = xgb_paras_tunning$L2_tune, # L2 Regularization (Ridge Regression)
                         alpha = xgb_paras_tunning$L1_tune # L1 Regularization (Lasso Regression)
                         )



#here we do one better then a validation set, we use cross validation to
#expand the amount of info we have!

    # if(core){
    # cl <- parallel::makeCluster(core)
    # doParallel::registerDoParallel(cl)
    # }

    xgb_trcontrol = caret::trainControl(method = "cv",
                                        number = 5,
                                        verboseIter = TRUE,
                                        returnData = FALSE,
                                        returnResamp = "all",
                                        allowParallel = TRUE)




 if(exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


   xgb_fit = caret::train(x = geno_model_ready_train,
                   y = pheno_clean[, response],
                   trControl = xgb_trcontrol,
                   tuneGrid = xgb_grid,
                   method = "xgbLinear")



   xgb_preds <- stats::predict(xgb_fit,
                        geno_model_ready_test,
                        reshape = TRUE)

   xgb_preds <- as.data.frame(xgb_preds)

   bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

   res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                   X_train = geno_model_ready_train)


 } else {

   if(!exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


     xgb_fit = caret::train(x = geno_model_ready_train,
                     y = pheno_clean[, response],
                     trControl = xgb_trcontrol,
                     tuneGrid = xgb_grid,
                     method = "xgbLinear")



     xgb_preds <- stats::predict(xgb_fit,
                          geno_model_ready_train,
                          reshape = TRUE)

     xgb_preds <- as.data.frame(xgb_preds)

     bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

     res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                     X_train = geno_model_ready_train)
   }


     }

    #####
    if((exists('omic1_model_ready_test') & exists('omic1_model_ready_train'))) {

      xgb_fit = caret::train(x = omic1_model_ready_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic1_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic1_model_ready_train)


    } else {

      if(!exists('omic1_model_ready_test') & exists('omic1_model_ready_train')) {


        xgb_fit = caret::train(x = omic1_model_ready_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic1_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic1_model_ready_train)
      }


    }

    #####
    if((exists('omic2_model_ready_test') & exists('omic2_model_ready_train'))) {

      xgb_fit = caret::train(x = omic2_model_ready_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic2_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic2_model_ready_train)


    } else {

      if(!exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {


        xgb_fit = caret::train(x = omic2_model_ready_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic2_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic2_model_ready_train)
      }


    }

  ####
    #####
    if((exists('omic3_model_ready_test') & exists('omic3_model_ready_train'))) {

      xgb_fit = caret::train(x = omic3_model_ready_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic3_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic3_model_ready_train)


    } else {

      if(!exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {


        xgb_fit = caret::train(x = omic3_model_ready_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic3_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic3_model_ready_train)

      }


    }


    #####
    if((exists('geno_omic1_test') & exists('geno_omic1_train'))) {

      xgb_fit = caret::train(x = geno_omic1_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic1_train)


    } else {

      if(!exists('geno_omic1_test') & exists('geno_omic1_train')) {


        xgb_fit = caret::train(x = geno_omic1_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic1_train)
      }


    }

    #####
    if((exists('geno_omic2_test') & exists('geno_omic2_train'))) {

      xgb_fit = caret::train(x = geno_omic2_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic2_train)


    } else {

      if(!exists('geno_omic2_test') & exists('geno_omic2_train')) {


        xgb_fit = caret::train(x = geno_omic2_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic2_train)

      }


    }

    #####
    #####
    if((exists('geno_omic3_test') & exists('geno_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic2_train)


    } else {

      if(!exists('geno_omic3_test') & exists('geno_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic3_train)

      }


    }

    #####
    if((exists('omic1_omic2_test') & exists('omic1_omic2_train'))) {

      xgb_fit = caret::train(x = omic1_omic2_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic3_train)


    } else {

      if(!exists('omic1_omic2_test') & exists('omic1_omic2_train')) {


        xgb_fit = caret::train(x = omic1_omic2_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic1_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic1_omic2_train)
      }


    }

    #####
    if((exists('omic1_omic3_test') & exists('omic1_omic3_train'))) {

      xgb_fit = caret::train(x = omic1_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic1_omic3_train)


    } else {

      if(!exists('omic1_omic3_test') & exists('omic1_omic3_train')) {


        xgb_fit = caret::train(x = omic1_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic1_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic1_omic3_train)
      }


    }

    #####
    if((exists('omic2_omic3_test') & exists('omic2_omic3_train'))) {

      xgb_fit = caret::train(x = omic2_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           as.matrix(omic2_omic3_test),
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic2_omic3_train)


    } else {

      if(!exists('omic2_omic3_test') & exists('omic2_omic3_train')) {


        xgb_fit = caret::train(x = omic2_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic2_omic3_train)

      }


    }

    #####
    if((exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic2_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic1_omic2_train)


    } else {

      if(!exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic2_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic1_omic2_train)

      }


    }

    #####
    if((exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic1_omic3_train)


    } else {

      if(!exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic1_omic3_train)

      }


    }

    #####
    if((exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic2_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic2_omic3_train)


    } else {

      if(!exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic2_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic2_omic3_train)
      }


    }

    #####
    if((exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = omic1_omic2_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic1_omic2_omic3_train)


    } else {

      if(!exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {


        xgb_fit = caret::train(x = omic1_omic2_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             omic1_omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic1_omic2_omic3_train)
      }


    }

    #####
    if((exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train'))) {

      xgb_fit = caret::train(x = geno_omic1_omic2_omic3_train,
                      y = pheno_clean[, response],
                      trControl = xgb_trcontrol,
                      tuneGrid = xgb_grid,
                      method = "xgbLinear")



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_omic1_omic2_omic3_train)


    } else {

      if(!exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {


        xgb_fit = caret::train(x = geno_omic1_omic2_omic3_train,
                        y = pheno_clean[, response],
                        trControl = xgb_trcontrol,
                        tuneGrid = xgb_grid,
                        method = "xgbLinear")



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        bestTune <- c(xgb_fit$bestTune, xgb_fit$method)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_omic1_omic2_omic3_train)
      }


    }
    #######################################################
    ##########################################################
    ####### Start when their is no need for tunning
    ########################################################
    ########################################################
    } else {

    xgb_params <- list(
      booster = booster,
      eta = learning_rate,
      max_depth = max_depth, #This indicates how deep the built tree can be.
      #The deeper the tree, the more splits it has and it captures more
      #information about how the data. We fit a decision tree with depths
      #ranging from 1 to 32 and plot the training and test errors
      gamma = 4,
      subsample = subsample,
      colsample_bytree = 1,
      objective = "reg:squarederror",
      eval_metric = c("rmse", "rmsle", "mape")
    )





    if(exists('geno_model_ready_test') & exists('geno_model_ready_train')) {

      geno_model_ready_train <- xgboost::xgb.DMatrix(data = geno_model_ready_train,
                                                     label = pheno[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = geno_model_ready_train)



    } else {

      if(!exists('geno_model_ready_test') & exists('geno_model_ready_train')) {


        geno_model_ready_train <- xgboost::xgb.DMatrix(data = geno_model_ready_train,
                                                        label = pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             geno_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = geno_model_ready_train)


      }


    }

    #####
    if((exists('omic1_model_ready_test') & exists('omic1_model_ready_train'))) {

      omic1_model_ready_train <- xgboost::xgb.DMatrix(data = omic1_model_ready_train,
                                                      label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                           omic1_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic1_model_ready_train)


    } else {

      if(!exists('omic1_model_ready_test') & exists('omic1_model_ready_train')) {


        omic1_model_ready_train <- xgboost::xgb.DMatrix(data = omic1_model_ready_train,
                                                         label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                             omic1_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic1_model_ready_train)


      }


    }

    #####
    if((exists('omic2_model_ready_test') & exists('omic2_model_ready_train'))) {

      omic2_model_ready_train <- xgboost::xgb.DMatrix(data = omic2_model_ready_train,
                                                      label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic2_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           omic2_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic2_model_ready_train)



    } else {

      if(!exists('omic2_model_ready_test') & exists('omic2_model_ready_train')) {


        omic2_model_ready_train <- xgboost::xgb.DMatrix(data = omic2_model_ready_train,
                                                         label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic2_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )

        rm(omic2_model_ready_train)

        xgb_preds <- stats::predict(xgb_fit,
                             omic2_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic2_model_ready_train)


      }


    }

    ####
    #####
    if((exists('omic3_model_ready_test') & exists('omic3_model_ready_train'))) {

      omic3_model_ready_train <- xgboost::xgb.DMatrix(data = omic3_model_ready_train,
                                                      label =  pheno_clean[, response])


      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic3_model_ready_train,
        nrounds = iteration,
        verbose = 1
      )

      xgb_preds <- stats::predict(xgb_fit,
                           omic3_model_ready_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train = omic3_model_ready_train)



    } else {

      if(!exists('omic3_model_ready_test') & exists('omic3_model_ready_train')) {


        omic3_model_ready_train <- xgboost::xgb.DMatrix(data = omic3_model_ready_train,
                                                         label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic3_model_ready_train,
          nrounds = iteration,
          verbose = 1
        )




        xgb_preds <- stats::predict(xgb_fit,
                             omic3_model_ready_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train = omic3_model_ready_train)


      }


    }


    #####
    if((exists('geno_omic1_test') & exists('geno_omic1_train'))) {

      geno_omic1_train <- xgboost::xgb.DMatrix(data = geno_omic1_train,
                                               label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic1_train)



    } else {

      if(!exists('geno_omic1_test') & exists('geno_omic1_train')) {


        geno_omic1_train <- xgboost::xgb.DMatrix(data = geno_omic1_train,
                                                  label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic1_train)


      }


    }

    #####
    if((exists('geno_omic2_test') & exists('geno_omic2_train'))) {

      geno_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic2_train,
                                               label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic2_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic2_train)



    } else {

      if(!exists('geno_omic2_test') & exists('geno_omic2_train')) {


        geno_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic2_train,
                                                  label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic2_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic2_train)


      }


    }

    #####
    #####
    if((exists('geno_omic3_test') & exists('geno_omic3_train'))) {

      geno_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic3_train,
                                               label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic3_train,
        nrounds = iteration,
        verbose = 1
      )

      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic3_train)

    } else {

      if(!exists('geno_omic3_test') & exists('geno_omic3_train')) {


        geno_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic3_train,
                                                  label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic3_train)


      }


    }

    #####
    if((exists('omic1_omic2_test') & exists('omic1_omic2_train'))) {

      omic1_omic2_train <- xgboost::xgb.DMatrix(data = omic1_omic2_train,
                                                label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic2_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  omic1_omic2_train)


    } else {

      if(!exists('omic1_omic2_test') & exists('omic1_omic2_train')) {


        omic1_omic2_train <- xgboost::xgb.DMatrix(data = omic1_omic2_train,
                                                   label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic2_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             omic1_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  omic1_omic2_train)


      }


    }

    #####
    if((exists('omic1_omic3_test') & exists('omic1_omic3_train'))) {

      omic1_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic3_train,
                                                label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  omic1_omic3_train)



    } else {

      if(!exists('omic1_omic3_test') & exists('omic1_omic3_train')) {


        omic1_omic3_train_ <- xgboost::xgb.DMatrix(data = omic1_omic3_train,
                                                   label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic3_train_,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             omic1_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  omic1_omic3_train)


      }


    }

    #####
    if((exists('omic2_omic3_test') & exists('omic2_omic3_train'))) {

      omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic2_omic3_train,
                                                label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )



      xgb_preds <- stats::predict(xgb_fit,
                           omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)

      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  omic2_omic3_train)



    } else {

      if(!exists('omic2_omic3_test') & exists('omic2_omic3_train')) {


        omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic2_omic3_train,
                                                   label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  omic2_omic3_train)


      }


    }

    #####
    if((exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train'))) {

      geno_omic1_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_train,
                                                     label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic2_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_omic2_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic1_omic2_train)

    } else {

      if(!exists('geno_omic1_omic2_test') & exists('geno_omic1_omic2_train')) {


        geno_omic1_omic2_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_train,
                                                        label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic2_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic1_omic2_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic1_omic2_train)

      }


    }

    #####
    if((exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train'))) {

      geno_omic1_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic3_train,
                                                     label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic1_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic1_omic3_train)

    } else {

      if(!exists('geno_omic1_omic3_test') & exists('geno_omic1_omic3_train')) {


        geno_omic1_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic3_train,
                                                        label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic3_train,
          nrounds = iteration,
          verbose = 1
        )



        xgb_preds <- stats::predict(xgb_fit,
                             as.matrix(geno_omic1_omic3_train),
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)


        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic1_omic3_train)

      }


    }

    #####
    if((exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train'))) {

      geno_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic2_omic3_train,
                                                     label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           geno_omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic2_omic3_train)

    } else {

      if(!exists('geno_omic2_omic3_test') & exists('geno_omic2_omic3_train')) {


        geno_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic2_omic3_train,
                                                        label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             geno_omic2_omic3_train,
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic2_omic3_train)


      }


    }

    #####
    if((exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train'))) {

      omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic2_omic3_train,
                                                      label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = omic1_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           omic1_omic2_omic3_test,
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  omic1_omic2_omic3_train)

    } else {

      if(!exists('omic1_omic2_omic3_test') & exists('omic1_omic2_omic3_train')) {


        omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = omic1_omic2_omic3_train,
                                                         label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = omic1_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             as.matrix(omic1_omic2_omic3_train),
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  omic1_omic2_omic3_train)


      }


    }

    #####
    if((exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train'))) {

      geno_omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_omic3_train,
                                                           label =  pheno_clean[, response])

      xgb_fit <- xgboost::xgb.train(
        params = xgb_params,
        data = geno_omic1_omic2_omic3_train,
        nrounds = iteration,
        verbose = 1
      )


      xgb_preds <- stats::predict(xgb_fit,
                           as.matrix(geno_omic1_omic2_omic3_test),
                           reshape = TRUE)

      xgb_preds <- as.data.frame(xgb_preds)


      res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                      X_train =  geno_omic1_omic2_omic3_train)

    } else {

      if(!exists('geno_omic1_omic2_omic3_test') & exists('geno_omic1_omic2_omic3_train')) {


        geno_omic1_omic2_omic3_train <- xgboost::xgb.DMatrix(data = geno_omic1_omic2_omic3_train,
                                                              label =  pheno_clean[, response])

        xgb_fit <- xgboost::xgb.train(
          params = xgb_params,
          data = geno_omic1_omic2_omic3_train,
          nrounds = iteration,
          verbose = 1
        )


        xgb_preds <- stats::predict(xgb_fit,
                             as.matrix(geno_omic1_omic2_omic3_train),
                             reshape = TRUE)

        xgb_preds <- as.data.frame(xgb_preds)

        res_feature <- feature_impo_xgb(xgb_fit = xgb_fit,
                                        X_train =  geno_omic1_omic2_omic3_train)

      }


    }

    #summary(xgb_fit$params)


    } ## End of when no need for tunning.





      #}

  model_para <- c(nrounds = xgb_fit$niter, xgb_fit$params)

output = list(model_para,
              xgb_preds,
              res_feature,
              xgb_fit)


names(output) = c("model_parameters", "predicted_values",
                  "feature_weight", "trained_model")


}


return(output)



 }
