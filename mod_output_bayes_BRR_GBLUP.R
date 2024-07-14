
#' Title
#'
#' @param mod
#' @param ETA
#' @param gen_name
#' @param gkernel
#' @param gmatrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param pheno_data
#' @param heter_groups
#' @param omics_kernel_label
#' @param ...
#'
#' @return
#' @export
#' @importFrom magrittr |>
#' @examples
mod_output_bayes_BRRGBLUP <- function(
                                 mod=NULL,
                                 ETA=NULL,
                                 gen_name=NULL,
                                 gkernel=NULL,
                                 gmatrix = NULL,
                                 omic1_kernel=NULL,
                                 omic2_kernel=NULL,
                                 omic3_kernel=NULL,
                                 omics_kernel_label = list(omic1_kernel = NULL,
                                                           omic2_kernel = NULL,
                                                           omic3_kernel = NULL),
                                 pheno_data = NULL,
                                 heter_groups = NULL,
                                 bayes_para = NULL,
                                 ...){
  ##############
  msg <- "\n==================================================\n"

  g_use <- NULL
  if(inherits(omics_kernel_label,'list')){
    if(!all(sapply(omics_kernel_label, function(x){ is.null(x)}))!=FALSE){

      label <-  which(sapply(omics_kernel_label, function(x) !is.null(x)))

      print_lable <-  omics_kernel_label[label]

    } else {
      print_lable <-  NULL
    }

  }else {
    if(inherits(omics_kernel_label, "character")){

      print_lable <-  omics_kernel_label

    }

    if(is.null(omics_kernel_label)){
      print_lable <-  NULL

    }
  }
  ######## Ends


  if(!is.null(gmatrix)){
    gid_name <- rownames(gmatrix)
    g_retain <- gmatrix
  } else if (!is.null(gkernel)){
    gid_name <- rownames(gkernel)
    g_retain <- gkernel
  } else if (!is.null(omic1_kernel)){
    gid_name <- rownames(omic1_kernel)
  } else if (!is.null(omic2_kernel)){
    gid_name <- rownames(omic2_kernel)
  } else {
    if (!is.null(omic3_kernel)){
      gid_name <- rownames(omic3_kernel)
    }
  }
#####
  ######
  ### This is important because omic_kernel are compromised when there is more than  one environment
  ## So we need to original matrix to store back.
  if (!is.null(omic1_kernel)){
    omic1_retain <- omic1_kernel
  }

  if (!is.null(omic2_kernel)){
    omic2_retain <- omic2_kernel
  }

  if (!is.null(omic3_kernel)){
    omic3_retain <- omic3_kernel
  }

  ### Check bayes_parameter_check function in bayesians_preprocess for details
  ## The value here are the default values and assumed to be used when user did
  ## not provide the nIter and burnIn

  nIter <-  bayes_para[["nIter"]]
  burnIn <-   bayes_para[["burnIn"]]

  posindex <- (burnIn + 1):nIter

  #########
  #### When the genotype are present in more than one environment/location
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    if(!is.null(gmatrix)){
      gmatrix_copy <- gmatrix  ## For across Predicted value PEV, SE and REL
      gmatrix <- Zg%*%gmatrix%*%t(Zg)

      gmatrix <- grm_kernel_precheck(gmatrix)
    }

    if(!is.null(gkernel)){
      gmatrix_copy <- gkernel
      gkernel <- Zg%*%gkernel%*%t(Zg)

      gkernel <- grm_kernel_precheck(gkernel)
    }

    if(!is.null(omic1_kernel)){
      omic1_kernel_copy <- omic1_kernel
      omic1_kernel <- Zg%*%omic1_kernel%*%t(Zg)

      omic1_kernel <- grm_kernel_precheck(omic1_kernel)

    }

    if(!is.null(omic2_kernel)){
      omic2_kernel_copy <- omic2_kernel
      omic2_kernel <- Zg%*%omic2_kernel%*%t(Zg)

      omic2_kernel <- grm_kernel_precheck(omic2_kernel)
    }

    if(!is.null(omic3_kernel)){
      omic3_kernel_copy <- omic3_kernel
      omic3_kernel <- Zg%*%omic3_kernel%*%t(Zg)

      omic3_kernel <- grm_kernel_precheck(omic3_kernel)

    }

    # if(!is.null(heter_groups)){
    #   ZE <- model.matrix(~factor(pheno_data[,heter_groups])-1)
    #   ZEZE<-tcrossprod(ZE)
    #
    # }
    ENV <-  as.character(pheno_data[,heter_groups])

    Predicted_value <-  data.frame(name = pheno_data[,gen_name], Env = ENV, Predicted_value = mod$model$yHat)
  names(Predicted_value)[1] <-  gen_name
    } else {

    Predicted_value <- data.frame(name = pheno_data[,gen_name], Predicted_value = mod$model$yHat)
    names(Predicted_value)[1] <-  gen_name

    ENV = NULL

    if(!is.null(heter_groups)){
      heter_groups = NULL
    }
  }



  BIN <-  mod[["output_files_names"]][grepl("bin", mod[["output_files_names"]])]

  ### Extract Error variance
  Var_E <- scan(mod[["output_files_names"]][grepl("varE.dat", mod[["output_files_names"]])],
                     what = numeric(),
                     sep = "\n",
                    quiet = TRUE)

  Var_E <- Var_E[posindex]

  # calculate standard error
  Var_E_Se <- sd(Var_E)/sqrt(length(Var_E))


  ### If more than one M-matrix is provided. Effect will have it own coefficient
  ### These lines of code exttract the genomic variance and the error term
  if(length(BIN)>1){

    varB_files <- mod[["output_files_names"]][grepl("varB.dat", mod[["output_files_names"]])]


    if(length(BIN)==4){
      Var_U_1 <-  scan(varB_files[1],
                       what = numeric(),
                       sep = "\n",
                       quiet = TRUE)

      Var_U_1 <- Var_U_1[posindex]

      # calculate standard error
      Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))

      Var_U_2 <-  scan(varB_files[2],
                       what = numeric(),
                       sep = "\n",
                       quiet = TRUE)
      Var_U_2 <- Var_U_2[posindex]

      # calculate standard error
      Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))

      Var_U_3 <-  scan(varB_files[3],
                         what = numeric(),
                         sep = "\n",
                       quiet = TRUE)

      Var_U_3 <- Var_U_3[posindex]

      # calculate standard error
      Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))

      Var_U_4 <-  scan(varB_files[4],
                         what = numeric(),
                         sep = "\n",
                       quiet = TRUE)

      Var_U_4 <- Var_U_4[posindex]

      # calculate standard error
      Var_U_4_Se <- sd(Var_U_4)/sqrt(length(Var_U_4))
      #############################

      Var_U <-  Var_U_1 + Var_U_2 + Var_U_3 + Var_U_4

      genomic_h2 <- Var_U/(Var_U+Var_E)

      genomic_h2_Se <-  sd(genomic_h2)/sqrt(length(genomic_h2))
      genomic_h2 <-  mean(genomic_h2)

      Var_U_Se <-  sd(Var_U)/sqrt(length(Var_U))
      Var_U <- mean(Var_U)
      Var_U_1 <- mean(Var_U_1)
      Var_U_2 <- mean(Var_U_2)
      Var_U_3 <- mean(Var_U_3)
      Var_U_4 <- mean(Var_U_4)
      Var_E <- mean(Var_E)


    }

    if(length(BIN)==3){
      Var_U_1 <- scan(varB_files[1],
                         what = numeric(),
                         sep = "\n",
                      quiet = TRUE)

      Var_U_1 <- Var_U_1[posindex]
      # calculate standard error
      Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))


      Var_U_2 <- scan(varB_files[2],
                         what = numeric(),
                         sep = "\n",
                      quiet = TRUE)

      Var_U_2 <- Var_U_2[posindex]
      Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))

      Var_U_3 <- scan(varB_files[3],
                         what = numeric(),
                         sep = "\n",
                      quiet = TRUE)

      Var_U_3 <- Var_U_3[posindex]
      Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))
###########
      Var_U <- Var_U_1 + Var_U_2 + Var_U_3
      genomic_h2 <-  (Var_U)/(Var_U+Var_E)
      genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
      genomic_h2 <- mean(genomic_h2)

      Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
      Var_U <- mean(Var_U)
      Var_U_1 <- mean(Var_U_1)
      Var_U_2 <- mean(Var_U_2)
      Var_U_3 <- mean(Var_U_3)
      Var_E <-   mean(Var_E)


    }

    if(length(BIN)==2){
      Var_U_1 <- scan(varB_files[1],
                         what = numeric(),
                         sep = "\n",
                      quiet = TRUE)

      Var_U_1 <- Var_U_1[posindex]
      Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))

      Var_U_2 <- scan(varB_files[2],
                         what = numeric(),
                         sep = "\n",
                      quiet = TRUE)
      Var_U_2 <- Var_U_2[posindex]
      Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))


      Var_U <- Var_U_1 + Var_U_2

      genomic_h2 <- Var_U/(Var_U+Var_E)
      genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
      genomic_h2 <- mean(genomic_h2)

      Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
      Var_U <- mean(Var_U)
      Var_U_1 <- mean(Var_U_1)
      Var_U_2 <- mean(Var_U_2)
      Var_E <- mean(Var_E)

    }

  } else {

    Var_U <-  scan(mod[["output_files_names"]][grepl("varB.dat", mod[["output_files_names"]])],
                     what = numeric(),
                     sep = "\n")

    Var_U <- Var_U[posindex]
    Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))

    Var_E <- scan(mod[["output_files_names"]][grepl("varE.dat", mod[["output_files_names"]])],
                       what = numeric(),
                       sep = "\n")

    Var_E <- Var_E[posindex]
    Var_E_Se <-  sd(Var_E)/sqrt(length(Var_E))


    genomic_h2 <- Var_U/(Var_U+Var_E)
    genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
    genomic_h2 <- mean(genomic_h2)

    Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
    Var_U <- mean(Var_U)
    Var_E <-  mean(Var_E)

  }


  #### When there is more than one M matrix

  if(length(BIN)>1){
    for (aa in 1:length(BIN)) {

      if(ETA[["ETA_element_name"]][aa]=="gkernel" | ETA[["ETA_element_name"]][aa]=="gmatrix"){

        if(ETA[["ETA_element_name"]][aa]=="gkernel"){
          g_use = gkernel

          rm(gkernel)

        } else{

          if(ETA[["ETA_element_name"]][aa]=="gmatrix"){
            g_use = gmatrix

            rm(gmatrix)
          }


        }

      #if(ETA$ETA_element_name[aa]=="geno_data"){

        Bb=BGLR::readBinMat(BIN[aa])

         var_u <-  scan(varB_files[aa],
                          what = numeric(),
                          sep = "\n",
                        quiet = TRUE)

         var_u <- mean(var_u[posindex])

        res_Coeff_EBV_PEV_Rel_SE_gen <- Cal_Coeff_EBV_PEV_Rel_SE(
          beta = Bb,
          x_variable = g_use,
          gen_name = gen_name,
          var_u = var_u,
          heter_groups = heter_groups,
          hetero = ENV,
          gid_name = gid_name)


      }

      ##
      if(ETA$ETA_element_name[aa]=="omic1_kernel"){

        Bb=BGLR::readBinMat(BIN[aa])

        var_u <-  scan(varB_files[aa],
                       what = numeric(),
                       sep = "\n",
                       quiet = TRUE)

        var_u <- mean(var_u[posindex])

        res_Coeff_EBV_PEV_Rel_SE_omic1 <- Cal_Coeff_EBV_PEV_Rel_SE(
          beta = Bb,
          x_variable = omic1_kernel,
          gen_name = gen_name,
          var_u = var_u,
          heter_groups = heter_groups,
          hetero = ENV,
          gid_name = gid_name)


      }

      if(ETA[["ETA_element_name"]][aa]=="omic2_kernel"){

        Bb=BGLR::readBinMat(BIN[aa])

        var_u <-  scan(varB_files[aa],
                       what = numeric(),
                       sep = "\n",
                       quiet = TRUE)

        var_u <- mean(var_u[posindex])

        res_Coeff_EBV_PEV_Rel_SE_omic2 <- Cal_Coeff_EBV_PEV_Rel_SE(
          beta = Bb,
          x_variable = omic2_kernel,
          gen_name = gen_name,
          var_u = var_u,
          heter_groups = heter_groups,
          hetero = ENV,
          gid_name = gid_name)

      }

      if(ETA[["ETA_element_name"]][aa]=="omic3_kernel"){

        Bb=BGLR::readBinMat(BIN[aa])

        var_u <-  scan(varB_files[aa],
                       what = numeric(),
                       sep = "\n",
                       quiet = TRUE)

        var_u <- mean(var_u[posindex])

        res_Coeff_EBV_PEV_Rel_SE_omic3 <- Cal_Coeff_EBV_PEV_Rel_SE(
          beta = Bb,
          x_variable = omic3_kernel,
          gen_name = gen_name,
          var_u = var_u,
          heter_groups = heter_groups,
          hetero = ENV,
          gid_name = gid_name)

      }

    }
    #### start from here

    ### When you have just one M_matrix/X_matrix
  } else {


    if(ETA[["ETA_element_name"]][1]=="gkernel" | ETA[["ETA_element_name"]][1]=="gmatrix"){

      if(ETA[["ETA_element_name"]][1]=="gkernel"){
        g_use = gkernel

        rm(gkernel)

      } else{

        if(ETA[["ETA_element_name"]][1]=="gmatrix"){
          g_use = gmatrix

          rm(gmatrix)
        }



      }
    #if(ETA$ETA_element_name[1]=="geno_data"){

      Bb=BGLR::readBinMat(BIN)

      res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = g_use,
        gen_name = gen_name,
        var_u = Var_U,
        heter_groups = heter_groups,
        hetero = ENV,
        gid_name = gid_name)
    }

    ##
    if(ETA[["ETA_element_name"]][1]=="omic1_kernel"){


      Bb=BGLR::readBinMat(BIN)

      res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic1_kernel,
        gen_name = gen_name,
        var_u = Var_U,
        heter_groups = heter_groups,
        hetero = ENV,
        gid_name = gid_name)

    }

    if(ETA[["ETA_element_name"]][1]=="omic2_kernel"){

      Bb=BGLR::readBinMat(BIN)

      res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic2_kernel,
        gen_name = gen_name,
        var_u = Var_U,
        heter_groups = heter_groups,
        hetero = ENV,
        gid_name = gid_name)
    }

    if(ETA[["ETA_element_name"]][1]=="omic3_kernel"){

      Bb=BGLR::readBinMat(BIN)

      res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic3_kernel,
        gen_name = gen_name,
        var_u = Var_U,
        heter_groups = heter_groups,
        hetero = ENV,
        gid_name = gid_name)

    }

  }
  #####
  ## Create output for predicted value and residual value.
  ## The residual value dataframe also contain predicted value for two reasons
  #1) For ease of plotting
  #2) When testing set is present in the real world it is expected to be


  ### Residual value is only estimable for response value without NA
  tst <- which(is.na(mod$model$y))

  if(length(tst)!=0){
    Residual_value <- data.frame(
      name = NA,
      Predicted_value = mod$model$yHat[tst],
      Residual_value = (mod$model$y[tst] - mod$model$yHat[tst]),
      stringsAsFactors = FALSE)

  } else {

    Residual_value <- data.frame(
      name = NA,
      Predicted_value = mod$model$yHat,
      Residual_value = (mod$model$y - mod$model$yHat),
      stringsAsFactors = FALSE)

  }

####
  ####
  if(length(BIN)==1){


    if(!is.null(g_use)){
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        sep_pev_rel <- sep_pev_rel_gblup(geno_object = gmatrix_copy,
                                         va = Var_U,
                                         ve = Var_E)

        Across_env_Predicted_value <-  as.data.frame(Predicted_value |>
                                                       dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                       dplyr::summarise(Pred = mean(Predicted_value)) |>
                                                       dplyr::mutate(Std_error = sep_pev_rel$sep,
                                                                     PEV = sep_pev_rel$pev,
                                                                     Reliability = sep_pev_rel$rel)
        )

      } else {

      Predicted_value[, 1] <- rownames(g_use)
      #PEV <- apply(g_ebv, 1, var)
      Predicted_value <- Predicted_value |>
                         dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                                 PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                                 Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

      }


      if(length(tst)!=0){
        Residual_value[, 1] <- rownames(g_use)[tst]
      }else {

        Residual_value[, 1] <- rownames(g_use)
      }

    } else if(!is.null(omic1_kernel)){

      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        sep_pev_rel <- sep_pev_rel_gblup(geno_object = omic1_kernel_copy,
                                         va = Var_U,
                                         ve = Var_E)

        Across_env_Predicted_value <-  as.data.frame(Predicted_value |>
                                                       dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                       dplyr::summarise(Pred = mean(Predicted_value)) |>
                                                       dplyr::mutate(Std_error = sep_pev_rel$sep,
                                                                     PEV = sep_pev_rel$pev,
                                                                     Reliability = sep_pev_rel$rel)
        )

      } else {

      Predicted_value[, 1] <- rownames(omic1_kernel)

      Predicted_value <- Predicted_value |>
        dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                      PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                      Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

      }

      if(length(tst)!=0){
        Residual_value[, 1] <- rownames(omic1_kernel)[tst]
      }else {

        Residual_value[, 1] <- rownames(omic1_kernel)
      }


    } else if(!is.null(omic2_kernel)){


      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        sep_pev_rel <- sep_pev_rel_gblup(geno_object = omic2_kernel_copy,
                                         va = Var_U,
                                         ve = Var_E)

        Across_env_Predicted_value <-  as.data.frame(Predicted_value |>
                                                       dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                       dplyr::summarise(Pred = mean(Predicted_value)) |>
                                                       dplyr::mutate(Std_error = sep_pev_rel$sep,
                                                                     PEV = sep_pev_rel$pev,
                                                                     Reliability = sep_pev_rel$rel)
        )

      } else {

      Predicted_value[, 1] <- rownames(omic2_kernel)

      Predicted_value <- Predicted_value |>
        dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                      PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                      Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

      }

      if(length(tst)!=0){
        Residual_value[, 1] <- rownames(omic2_kernel)[tst]
      }else {

        Residual_value[, 1] <- rownames(omic2_kernel)
      }


    } else {

      if(!is.null(omic3_kernel)){

        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
          sep_pev_rel <- sep_pev_rel_gblup(geno_object = omic3_kernel_copy,
                                           va = Var_U,
                                           ve = Var_E)

          Across_env_Predicted_value <-  as.data.frame(Predicted_value |>
                                                         dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                         dplyr::summarise(Pred = mean(Predicted_value)) |>
                                                         dplyr::mutate(Std_error = sep_pev_rel$sep,
                                                                       PEV = sep_pev_rel$pev,
                                                                       Reliability = sep_pev_rel$rel)
          )

        } else {

        Predicted_value[, 1] <- rownames(omic3_kernel)

        Predicted_value <- Predicted_value |>
          dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                        PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                        Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

        }

        if(length(tst)!=0){
          Residual_value[, 1] <- rownames(omic3_kernel)[tst]
        } else {

          Residual_value[, 1] <- rownames(omic3_kernel)
        }


      }

    }

  } else {

    if(length(BIN)>1){

      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){


        Across_env_Predicted_value <-  as.data.frame(Predicted_value |>
                                                       dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                       dplyr::summarise(Pred = mean(Predicted_value)) |>
                                                       dplyr::mutate(Std_error = NA,
                                                                     PEV = NA,
                                                                     Reliability = NA)
        )

      }


    }
  }




  if(!is.null(g_use) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
    Variance_components <- bayes_variance_components(Var_U = Var_U,
                                                     Var_E = Var_E,
                                                     genomic_h2 = genomic_h2,
                                                     Var_U_Se = Var_U_Se,
                                                     Var_E_Se = Var_E_Se,
                                                     genomic_h2_Se = genomic_h2_Se)


    ### Combined all results into list
    Res <-  list(
      coefficients = res_Coeff_EBV_PEV_Rel_SE$Coefficient,
      EBV = res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
      Predicted_value =  Predicted_value,
      Variance_components = Variance_components,
      list(Geno_model_ready = g_retain),## This is to make it compatible in output format for when the omics is more than one
      mu = mod$model$mu)

    ## Add attribute/ name to the list
    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

  } else if(is.null(g_use) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Variance_components <- bayes_variance_components(Var_U = Var_U,
                                                     Var_E = Var_E,
                                                     genomic_h2 = genomic_h2,
                                                     Var_U_Se = Var_U_Se,
                                                     Var_E_Se = Var_E_Se,
                                                     genomic_h2_Se = genomic_h2_Se)


    ### Combined all results into list
    Res <-  list(
      coefficients = res_Coeff_EBV_PEV_Rel_SE$Coefficient,
      EBV = res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
      Predicted_value =  Predicted_value,
      Variance_components = Variance_components,
      list(Omic_model_ready = omic1_retain),## This is to make it compatible in output format for when the omics is more than one
      mu = mod$model$mu)

    ## Add attribute/ name to the list
    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")
    ##
    ## Add attribute/ name to the list

    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){

        rownames(Res$Variance_components)[1] <- paste(print_lable, "variance", sep="_")
        names(Res) <- c(
          "Coefficients",
          paste(print_lable, "estimated_breeding_value", sep="_"),
          "Predicted_value",
          "Variance_components",
          "M_matrix_model_ready",
          "Intercept")

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))

    }


  } else if(is.null(g_use) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Variance_components <- bayes_variance_components(Var_U = Var_U,
                                                     Var_E = Var_E,
                                                     genomic_h2 = genomic_h2,
                                                     Var_U_Se = Var_U_Se,
                                                     Var_E_Se = Var_E_Se,
                                                     genomic_h2_Se = genomic_h2_Se)


    ### Combined all results into list
    Res <-  list(
      coefficients = res_Coeff_EBV_PEV_Rel_SE$Coefficient,
      EBV = res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
      Predicted_value =  Predicted_value,
      Variance_components = Variance_components,
      list(Omic_model_ready = omic2_retain),## This is to make it compatible in output format for when the omics is more than one
      mu = mod$model$mu)

    ## Add attribute/ name to the list
    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")
    ##
    ####
    ## Add attribute/ name to the list

    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){

        rownames(Res$Variance_components)[1] <- paste(print_lable, "variance", sep="_")
        names(Res) <- c(
          "Coefficients",
          paste(print_lable, "estimated_breeding_value", sep="_"),
          "Predicted_value",
          "Variance_components",
          "M_matrix_model_ready",
          "Intercept")

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))

    }

  } else if(is.null(g_use) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Variance_components <- bayes_variance_components(Var_U = Var_U,
                                                     Var_E = Var_E,
                                                     genomic_h2 = genomic_h2,
                                                     Var_U_Se = Var_U_Se,
                                                     Var_E_Se = Var_E_Se,
                                                     genomic_h2_Se = genomic_h2_Se)


    ### Combined all results into list
    Res <-  list(
      coefficients = res_Coeff_EBV_PEV_Rel_SE$Coefficient,
      EBV = res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
      Predicted_value =  Predicted_value,
      Variance_components = Variance_components,
      list(Omic_model_ready = omic3_retain),## This is to make it compatible in output format for when the omics is more than one
      mu = mod$model$mu)

    ## Add attribute/ name to the list
    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")
    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1] <- paste(print_lable, "variance", sep="_")
        names(Res) <- c(
          "Coefficients",
          paste(print_lable, "estimated_breeding_value", sep="_"),
          "Predicted_value",
          "Variance_components",
          "M_matrix_model_ready",
          "Intercept")

      } else{
        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))

    }


  } else if(!is.null(g_use) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name

   #####
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    #####
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(g_use)[tst]
    }else {

      Residual_value[, 1] <- rownames(g_use)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient),
      EBV = list(Genomic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Omic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic_model_ready = omic1_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####


  } else if(!is.null(g_use) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name

    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)
    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(g_use)[tst]
    }else {

      Residual_value[, 1] <- rownames(g_use)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),
      EBV = list(Genomic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Omic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic_model_ready = omic2_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if((!is.null(g_use)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(g_use)[tst]
    }else {

      Residual_value[, 1] <- rownames(g_use)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Genomic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Omic_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if(is.null(g_use) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Omic1_coeff = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   Omic2_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),
      EBV = list(Omic1_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 Omic2_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_omics = list(Omic1_kernel_model_ready = omic1_retain,
                        Omic2_kernel_model_ready = omic2_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if(is.null(g_use) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name

    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic2_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic2_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Omic1_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                   Omic2_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Omic1_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                 Omic2_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_omics = list(Omic1_kernel_model_ready = omic2_retain,
                         Omic2_kernel_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if(is.null(g_use) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name

######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Omic1_coeff = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   Omic2_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Omic1_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 Omic2_estimated_breeding_value= res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_omics = list(omic1_kernel_model_ready = omic1_retain,
                        omic2_kernel_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if(!is.null(g_use) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic2_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic2_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U_3 = Var_U_3,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_3_Se = Var_U_3_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                   coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Estimated_breeding_value_1= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Estimated_breeding_value_2= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                 Estimated_breeding_value_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic1_kernel_model_ready = omic2_retain,
                        Omic2_kernel_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

    #####

  } else if(!is.null(g_use) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)
    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U_3 = Var_U_3,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_3_Se = Var_U_3_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Estimated_breeding_value_1= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Estimated_breeding_value_2= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 Estimated_breeding_value_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic1_kernel_model_ready = omic1_retain,
                        Omic2_kernel_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }

  } else if(!is.null(g_use) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)
    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U_3 = Var_U_3,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_3_Se = Var_U_3_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   Coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   Coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),
      EBV = list(Estimated_breeding_value_1= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 Estimated_breeding_value_2= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 Estimated_breeding_value_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_geno = list(Geno_model_ready = g_retain,
                        Omic1_kernel_model_ready = omic1_retain,
                        Omic2_kernel_model_ready = omic2_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    if(!is.null(print_lable)){

      print_lable = unlist(print_lable)


      if(length(print_lable)==1){
        rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
        ######
        names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
        ####

      } else {

        if(length(print_lable)>1){
          message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }


  } else if(is.null(g_use) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    Reliability <- 1 - (PEV / Var_U)
    Reliability = ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                          stringsAsFactors = FALSE)

    colnames(sum_EBV)[1] <- gen_name


    ######
    sum_EBV <- sum_EBV |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)

    ######

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(PEV),
                    PEV = PEV,
                    Reliability = Reliability)
    ###
    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_kernel)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_kernel)
    }

    ### Ends


    Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                      Var_U_2 = Var_U_2,
                                                      Var_U_3 = Var_U_3,
                                                      Var_U = Var_U,
                                                      Var_E = Var_E,
                                                      genomic_h2 = genomic_h2,
                                                      Var_U_1_Se = Var_U_1_Se,
                                                      Var_U_2_Se = Var_U_2_Se,
                                                      Var_U_3_Se = Var_U_3_Se,
                                                      Var_U_Se = Var_U_Se,
                                                      Var_E_Se = Var_E_Se,
                                                      genomic_h2_Se = genomic_h2_Se)

    Res <-  list(
      coeff = list(Coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   Coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                   Coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
      EBV = list(Estimated_breeding_value_1= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 Estimated_breeding_value_2= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                 Estimated_breeding_value_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
      sum_EBV,
      Predicted_value,
      Variance_components = Variance_components,
      omics_omics = list(Omic1_kernel_model_ready = omic1_retain,
                        Omic2_kernel_model_ready = omic2_retain,
                        Omic3_kernel_model_ready = omic3_retain
      ),
      mu = mod$model$mu)

    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "Intercept")

    ###
    ###
    if(!is.null(print_lable)){

      print_lable <- unlist(print_lable)


      if(length(print_lable)==3){
        rownames(Res$Variance_components)[1:3] = c(paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"), paste(print_lable[3], "variance", sep="_"))
        names(Res$Coefficients) = c(paste(print_lable[1], "coefficient", sep="_"), paste(print_lable[2], "coefficient", sep="_"), paste(print_lable[3], "coefficient", sep="_"))
        names(Res$Estimated_breeding_value) = c(paste(print_lable[1], "estimated_breeding_value", sep="_"),  paste(print_lable[2], "estimated_breeding_value", sep="_"), paste(print_lable[3], "estimated_breeding_value", sep="_"))
        #names(Res$Genetic_variance) = c(paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"), paste(print_lable[3], "variance", sep="_"))
        names(Res$M_matrix_model_ready) = c(paste(print_lable[1], "model_ready", sep="_"),  paste(print_lable[2], "model_ready", sep="_"), paste(print_lable[3], "model_ready", sep="_"))
      } else {

        if(length(print_lable)>3 | length(print_lable)<3){
          message(insight::print_color(paste(msg,paste("More than three Omic lables were provided. Default name was applied.")), "blue"))
        }

      }

    } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


    }


  } else{

    if(!is.null(g_use) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      g_omic_ebv = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
      PEV <- apply(g_omic_ebv, 1, var)
      Reliability <- 1 - (PEV / Var_U)
      Reliability = ifelse(Reliability<0, "Alias", Reliability)

      sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                            Estimated_breeding_value = res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
                            stringsAsFactors = FALSE)

      colnames(sum_EBV)[1] <- gen_name


      ######
      sum_EBV <- sum_EBV |>
        dplyr::mutate(Std_error = sqrt(PEV),
                      PEV = PEV,
                      Reliability = Reliability)

      ######

      Predicted_value <- Predicted_value |>
        dplyr::mutate(Std_error = sqrt(PEV),
                      PEV = PEV,
                      Reliability = Reliability)
      ###
      ###
      if(length(tst)!=0){
        Residual_value[, 1] <- rownames(omic1_kernel)[tst]
      }else {

        Residual_value[, 1] <- rownames(omic1_kernel)
      }

      ### Ends


      Variance_components <- bayes_variance_components (Var_U_1 = Var_U_1,
                                                        Var_U_2 = Var_U_2,
                                                        Var_U_3 = Var_U_3,
                                                        Var_U_4 = Var_U_4,
                                                        Var_U = Var_U,
                                                        Var_E = Var_E,
                                                        genomic_h2 = genomic_h2,
                                                        Var_U_1_Se = Var_U_1_Se,
                                                        Var_U_2_Se = Var_U_2_Se,
                                                        Var_U_3_Se = Var_U_3_Se,
                                                        Var_U_4_Se = Var_U_4_Se,
                                                        Var_U_Se = Var_U_Se,
                                                        Var_E_Se = Var_E_Se,
                                                        genomic_h2_Se = genomic_h2_Se)

      Res <-  list(
        coeff = list(Coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                     Coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                     Coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                     Coefficients_4 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
        EBV = list(Estimated_breeding_value_1= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                   Estimated_breeding_value_2= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                   Estimated_breeding_value_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                   Estimated_breeding_value_4 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
        sum_EBV,
        Predicted_value,
        Variance_components = Variance_components,
        omics_geno = list(Geno_model_ready = g_retain,
                          Omic1_kernel_model_ready = omic1_retain,
                          Omic2_kernel_model_ready = omic2_retain,
                          Omic3_kernel_model_ready = omic3_retain
        ),
        mu = mod$model$mu)

      names(Res) <- c(
        "Coefficients",
        "Estimated_breeding_value",
        "Total_estimated_breeding_value",
        "Predicted_value",
        "Variance_components",
        "M_matrix_model_ready",
        "Intercept")

      ###
      ###
      ###
      if(!is.null(print_lable)){

        print_lable <- unlist(print_lable)


        if(length(print_lable)==3){
          rownames(Res$Variance_components)[1:4] = c("Genomic_variance", paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"), paste(print_lable[3], "variance", sep="_"))
          names(Res$Coefficients) = c("Geno_coefficient", paste(print_lable[1], "coefficient", sep="_"), paste(print_lable[2], "coefficient", sep="_"), paste(print_lable[3], "coefficient", sep="_"))
          names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value", paste(print_lable[1], "estimated_breeding_value", sep="_"), paste(print_lable[2], "estimated_breeding_value", sep="_"), paste(print_lable[3], "estimated_breeding_value", sep="_"))
          #names(Res$Genetic_variance) = c("Genomic_variance", paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"), paste(print_lable[3], "variance", sep="_"))
          names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable[1], "model_ready", sep="_"), paste(print_lable[2], "model_ready", sep="_"), paste(print_lable[3], "model_ready", sep="_"))
        } else {

          if(length(print_lable)>3 | length(print_lable)<3){
            message(insight::print_color(paste(msg,paste("More than three Omic lables were provided. Default name was applied.")), "blue"))
          }

        }

      } else {

        message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


      }


    }

    ####

  }

####
  if(exists("Across_env_Predicted_value")){
 ### Add the across prediction after the Predicted_value for each environment
    aft <-  which(names(Res)=="Predicted_value")

    Res <- append(Res, list(Across_env_Predicted_value = Across_env_Predicted_value), after = aft)

    NN <-  paste("Across", paste(heter_groups, "Predicted_value", sep = "_"), sep = "_")

    names(Res)[(aft+1)] <- NN
  }

  #}
  ### remove the generated output files from the working directory
  unlink(mod[["output_files_names"]])
  return(Res)
}





