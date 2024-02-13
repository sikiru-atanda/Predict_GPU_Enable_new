
#' Title
#'
#' @param mod
#' @param ETA
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param gen_name
#' @param GS_model
#' @param omics_data_label
#' @param bayes_para
#' @param ...
#'@importFrom magrittr |>
#' @return
#' @export
#'
#' @examples
mod_output_bayes <- function(mod=NULL,
                             ETA=NULL,
                             gen_name=NULL,
                             geno_data=NULL,
                             omic1_data=NULL,
                             omic2_data=NULL,
                             omic3_data=NULL,
                             omics_data_label = list(omic1_data = NULL,
                                                     omic2_data = NULL,
                                                     omic3_data = NULL),
                             bayes_para = NULL,
                             GS_model = NULL,
                             ...){

  msg <- sprintf("==================================================\n")
  ### Check if the user provide lable/name for the omics data

  if(typeof(omics_data_label)=='list'){
  if(!all(sapply(omics_data_label, function(x){ is.null(x)}))!=FALSE){

    label <-  which(sapply(omics_data_label, function(x) !is.null(x)))

    print_lable <-  omics_data_label[label]

  } else {
    print_lable <-  NULL
  }

  }else {
    if(typeof(omics_data_label)=="character"){

      print_lable <-  omics_data_label

    }

    if(is.null(omics_data_label)){
      print_lable <-  NULL

    }
  }

  if(!is.null(geno_data)){
    gid_name <- rownames(geno_data)
  } else if (!is.null(omic1_data)){
    gid_name <- rownames(omic1_data)
  } else if (!is.null(omic2_data)){
    gid_name <- rownames(omic2_data)
  } else {
    if (!is.null(omic3_data)){
      gid_name <- rownames(omic3_data)
    }
  }
# sik$ETA_element_name
# TT = DT$output_files_names

  # if(!is.null(omic1_data)){
  #
  #   H = omic1_data%*% solve(t(omic1_data) %*% omic1_data) %*% t(omic1_data)
  #
  #   vc_mY = H*var(mod$model$y, na.rm = TRUE)*t(H)
  #
  #   se_omic1 <- sqrt(diag(vc_mY))
  #
  #   rm(H, vc_mY)
  # }
  #
  # if(!is.null(omic2_data)){
  #
  #   H = omic2_data%*% solve(t(omic2_data) %*% omic2_data) %*% t(omic2_data)
  #
  #   vc_mY = H*var(mod$model$y, na.rm = TRUE)*t(H)
  #
  #   se_omic2 <- sqrt(diag(vc_mY))
  #
  #   rm(H, vc_mY)
  # }
  #
  # if(!is.null(omic3_data)){
  #
  #   H = omic3_data%*% solve(t(omic3_data) %*% omic3_data) %*% t(omic3_data)
  #
  #   vc_mY = H*var(mod$model$y, na.rm = TRUE)*t(H)
  #
  #   se_omic3 <- sqrt(diag(vc_mY))
  #
  #   rm(H, vc_mY)
  # }
  #
  # if(!is.null(geno_data)){
  #
  #   H = geno_data%*% solve(t(geno_data) %*% geno_data) %*% t(geno_data)
  #
  #   vc_mY = H*var(mod$model$y, na.rm = TRUE)*t(H)
  #
  #   se_geno <- sqrt(diag(vc_mY))
  #
  #   rm(H, vc_mY)
  # }


    ### Check bayes_parameter_check function in bayesians_preprocess for details
    ## The value here are the default values and assumed to be used when user did
    ## not provide the nIter and burnIn

    nIter <- bayes_para$nIter
    burnIn <- bayes_para$burnIn

    posindex <- (burnIn + 1):nIter




  if(GS_model == "BL") {
    GS_model <- "lambda"
  }

BIN <- mod$output_files_names[grepl("bin", mod$output_files_names)]

### Extract Error variance

Var_E <- scan(mod$output_files_names[grepl("varE.dat", mod$output_files_names)],
                   what = numeric(),
                   sep = "\n", quiet =TRUE)

Var_E <- Var_E[posindex]

# calculate standard error
Var_E_Se <- sd(Var_E)/sqrt(length(Var_E))



### If more than one M-matrix is provided. Effect will have it own coefficient
### These lines of code extract the genomic variance and the error term
if(length(BIN)>1){

  if(GS_model=="BRR"){
  varB_files <- mod$output_files_names[grepl("varB.dat", mod$output_files_names)]
  } else {
  varB_files <- mod$output_files_names[grepl(paste(GS_model,"dat", sep = "."), mod$output_files_names)]
}

  if(length(BIN)==4){
    #### This was done because the output frm BGLR different from each model espcially BRR from others
    if(GS_model== "BRR"){

            Var_U_1 <- scan(varB_files[1],
                                 what = numeric(),
                                 sep = "\n", quiet =TRUE)

            Var_U_1 <- Var_U_1[posindex]

            # calculate standard error
            Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))

         ########
            Var_U_2 <- scan(varB_files[2],
                            what = numeric(),
                            sep = "\n", quiet =TRUE)

            Var_U_2 <- Var_U_2[posindex]

            # calculate standard error
            Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))

        ######

            Var_U_3 <- scan(varB_files[3],
                            what = numeric(),
                            sep = "\n", quiet =TRUE)

            Var_U_3 <- Var_U_3[posindex]

            # calculate standard error
            Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))
        #####

            Var_U_4 <- scan(varB_files[4],
                            what = numeric(),
                            sep = "\n", quiet =TRUE)

            Var_U_4 <- Var_U_4[posindex]

            # calculate standard error
            Var_U_4_Se <- sd(Var_U_4)/sqrt(length(Var_U_4))


    } else {

            if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
              Var_U_1 <- utils::read.table(varB_files[1])
            } else{
             Var_U_1 <- utils::read.table(varB_files[1],
                                          skip = 1)
            }


           Var_U_1 <- as.data.frame(tidyr::separate_rows(Var_U_1))[posindex, 1]

           Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))
          #Var_U_1 <- mean(as.data.frame(tidyr::separate_rows(Var_U_1))[, 1])
          ####

           if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
             Var_U_2 <- utils::read.table(varB_files[2])
           } else{
             Var_U_2 <- utils::read.table(varB_files[2],
                                          skip = 1)
           }

           Var_U_2 <- as.data.frame(tidyr::separate_rows(Var_U_2))[posindex, 1]
           Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))
          #Var_U_2 <- mean(as.data.frame(tidyr::separate_rows(Var_U_2))[, 1])

           if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
             Var_U_3 <- utils::read.table(varB_files[3])
           } else{
             Var_U_3 <- utils::read.table(varB_files[3],
                                          skip = 1)
           }

           Var_U_3 <- as.data.frame(tidyr::separate_rows(Var_U_3))[posindex, 1]
           Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))

          #Var_U_3 <- mean(as.data.frame(tidyr::separate_rows(Var_U_3))[, 1])

           if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
             Var_U_4 <- utils::read.table(varB_files[4])
           } else{
             Var_U_4 <- utils::read.table(varB_files[4],
                                          skip = 1)
           }

           Var_U_4 <- as.data.frame(tidyr::separate_rows(Var_U_4))[posindex, 1]
           Var_U_4_Se <- sd(Var_U_4)/sqrt(length(Var_U_4))
          #Var_U_4 <- mean(as.data.frame(tidyr::separate_rows(Var_U_4))[, 1])

    }

          Var_U <- (Var_U_1 + Var_U_2 + Var_U_3 + Var_U_4)

          genomic_h2 <- Var_U/(Var_U+Var_E)
          genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
          genomic_h2 <- mean(genomic_h2)

          Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
          Var_U <- mean(Var_U)
          Var_U_1 <- mean(Var_U_1)
          Var_U_2 <- mean(Var_U_2)
          Var_U_3 <- mean(Var_U_3)
          Var_U_4 <- mean(Var_U_4)
          Var_E <- mean(Var_E)

  }

  if(length(BIN)==3){

      if(GS_model== "BRR"){
            Var_U_1 <- scan(varB_files[1],
                               what = numeric(),
                               sep = "\n", quiet =TRUE)

            Var_U_1 <- Var_U_1[posindex]
            Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))

            Var_U_2 <- scan(varB_files[2],
                               what = numeric(),
                               sep = "\n", quiet =TRUE)
            Var_U_2 <- Var_U_2[posindex]
            Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))

            Var_U_3 <- scan(varB_files[3],
                               what = numeric(),
                               sep = "\n", quiet =TRUE)
            Var_U_3 <- Var_U_3[posindex]
            Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))

        } else {

          if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
              Var_U_1 <- utils::read.table(varB_files[1])
            } else{
              Var_U_1 <- utils::read.table(varB_files[1],
                                           skip = 1)
            }

          Var_U_1 <- as.data.frame(tidyr::separate_rows(Var_U_1))[posindex, 1]
          Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))
          #Var_U_1 <- mean(as.data.frame(tidyr::separate_rows(Var_U_1))[, 1])
          ####

        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
            Var_U_2 <- utils::read.table(varB_files[2])
          } else{
            Var_U_2 <- utils::read.table(varB_files[2],
                                         skip = 1)
          }

          Var_U_2 <- as.data.frame(tidyr::separate_rows(Var_U_2))[posindex, 1]
          Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))
          #Var_U_2 <- mean(as.data.frame(tidyr::separate_rows(Var_U_2))[, 1])
      ######
        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
            Var_U_3 <- utils::read.table(varB_files[3])
          } else{
            Var_U_3 <- utils::read.table(varB_files[3],
                                         skip = 1)
          }

          Var_U_3 <- as.data.frame(tidyr::separate_rows(Var_U_3))[posindex, 1]
          Var_U_3_Se <- sd(Var_U_3)/sqrt(length(Var_U_3))
          #Var_U_3 <-  mean(as.data.frame(tidyr::separate_rows(Var_U_3))[, 1])

    }

          Var_U <- Var_U_1 + Var_U_2 + Var_U_3
          genomic_h2 <-  (Var_U)/(Var_U+Var_E)
          genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
          genomic_h2 <- mean(genomic_h2)

          Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
          Var_U <- mean(Var_U)
          Var_U_1 <- mean(Var_U_1)
          Var_U_2 <- mean(Var_U_2)
          Var_U_3 <- mean(Var_U_3)
          Var_E <- mean(Var_E)



  }

  if(length(BIN)==2){

    if(GS_model== "BRR"){
      Var_U_1 <- scan(varB_files[1],
                         what = numeric(),
                         sep = "\n", quiet =TRUE)
      Var_U_1 <- Var_U_1[posindex]
      Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))

      Var_U_2 <- scan(varB_files[2],
                         what = numeric(),
                         sep = "\n", quiet =TRUE)
      Var_U_2 <- Var_U_2[posindex]
      Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))

    } else {

      if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
        Var_U_1 <- utils::read.table(varB_files[1])
      } else{
        Var_U_1 <- utils::read.table(varB_files[1],
                                     skip = 1)
      }

    Var_U_1 <- as.data.frame(tidyr::separate_rows(Var_U_1))[posindex, 1]
    Var_U_1_Se <- sd(Var_U_1)/sqrt(length(Var_U_1))
    #Var_U_1 <- mean(as.data.frame(tidyr::separate_rows(Var_U_1))[, 1])
    ####

    if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
      Var_U_2 <- utils::read.table(varB_files[2])
    } else{
      Var_U_2 <- utils::read.table(varB_files[2],
                                   skip = 1)
    }

    Var_U_2 <- as.data.frame(tidyr::separate_rows(Var_U_2))[posindex, 1]
    Var_U_2_Se <- sd(Var_U_2)/sqrt(length(Var_U_2))
    #Var_U_2 <- mean(as.data.frame(tidyr::separate_rows(Var_U_2))[, 1])

}

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

  if(GS_model== "BRR"){
    Var_U <- scan(mod$output_files_names[grepl("varB.dat", mod$output_files_names)],
                     what = numeric(),
                     sep = "\n", quiet =TRUE)
    Var_U <- Var_U[posindex]
    Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
    ### The output for the variance was recently change in BGLR, this was fixed
    ## accordingly.
    ## It thus expedient to always crosscheck for any new update with the
    ## dependency package
  } else {
    if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
    Var_U <- utils::read.table(mod$output_files_names[grepl(paste(GS_model,"dat", sep = "."), mod$output_files_names)])
    } else {

      Var_U <- utils::read.table(mod$output_files_names[grepl(paste(GS_model,"dat", sep = "."), mod$output_files_names)], skip = 1)
  }
    Var_U <- as.data.frame(tidyr::separate_rows(Var_U))[posindex, 1]
    Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
    #Var_U <- mean(as.data.frame(tidyr::separate_rows(Var_U))[, 1])

  }
    Var_E <- scan(mod$output_files_names[grepl("varE.dat", mod$output_files_names)],
                       what = numeric(),
                       sep = "\n", quiet =TRUE)

    Var_E <- Var_E[posindex]
    Var_E_Se <-  sd(Var_E)/sqrt(length(Var_E))

    genomic_h2 <- Var_U/(Var_U+Var_E)
    genomic_h2_Se <- sd(genomic_h2)/sqrt(length(genomic_h2))
    genomic_h2 <- mean(genomic_h2)

    Var_U_Se <- sd(Var_U)/sqrt(length(Var_U))
    Var_U <- mean(Var_U)
    Var_E <-  mean(Var_E)

}


#### When there is more than one M matrix get the coefficient and BV

if(length(BIN)>1){
  for (aa in 1:length(BIN)) {

    if(ETA$ETA_element_name[aa]=="geno_data"){

      ##Read the stored the posterior samples of the random effect (markers/omics)

      Bb <- BGLR::readBinMat(BIN[aa])

      # g_ebv <- as.matrix(geno_data)%*%t(Bb)
      #
      # # Posterior means of marker effects/coefficient
      # coeff <- data.frame(X_variables = colnames(geno_data),
      #                     coeff =colMeans(Bb),
      #                     stringsAsFactors = FALSE)
      # # Genomic estimated breeding values
      # GEBV <- data.frame(names = rownames(geno_data),
      #                    GEBV=geno_data%*%coeff,
      #                    #std.erro = se_geno,
      #                    stringsAsFactors = FALSE)
      # colnames(GEBV)[1] <- gen_name
      # #GEBV = data.frame(rowMeans(ebv))
      #
      # PEV <- apply(g_ebv, 1, var)
      #
      # GEBV$Std_error <- sqrt(PEV)
      # GEBV$PEV <- PEV

      # var_u= mean(scan(varB_files[aa],
      #                  what = numeric(),
      #                  sep = "\n"))

      if(GS_model== "BRR"){
        var_u <- scan(varB_files[aa],
                         what = numeric(),
                         sep = "\n", quiet =TRUE)

        var_u <- mean(var_u[posindex])
        #var_u_Se = sd(var_u)/sqrt(length(var_u))

        ### The output for the variance was recently change in BGLR, this was fixed
        ## accordingly.
        ## It thus expedient to always crosscheck for any new update with the
        ## dependency package
      } else {
        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
        var_u <- utils::read.table(varB_files[aa])

        } else {

          var_u <- utils::read.table(varB_files[aa], skip = 1)
        }
        var_u <- mean(as.data.frame(tidyr::separate_rows(var_u))[posindex, 1])
        #var_u_Se = sd(var_u)/sqrt(length(var_u))
        #var_u <- mean(as.data.frame(tidyr::separate_rows(var_u))[, 1])

      }
      ##Estimate of reliability
      #GEBV$Reliability <- 1 - (PEV /var_u)
      ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
      ## and SE
      res_Coeff_EBV_PEV_Rel_SE_gen <- Cal_Coeff_EBV_PEV_Rel_SE(
                                                           beta = Bb,
                                                           x_variable = geno_data,
                                                           gen_name = gen_name,
                                                           var_u = var_u,
                                                           gid_name =gid_name)
    }

    ##
    if(ETA$ETA_element_name[aa]=="omic1_data"){

      Bb <- BGLR::readBinMat(BIN[aa])

      if(GS_model== "BRR"){
        var_u <- scan(varB_files[aa],
                           what = numeric(),
                           sep = "\n", quiet =TRUE)
        var_u <- mean(var_u[posindex])

        ### The output for the variance was recently change in BGLR, this was fixed
        ## accordingly.
        ## It thus expedient to always crosscheck for any new update with the
        ## dependency package
      } else {
        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
        var_u <- utils::read.table(varB_files[aa])

        } else {
        var_u <- utils::read.table(varB_files[aa], skip = 1)
        }

        var_u <- mean(as.data.frame(tidyr::separate_rows(var_u))[posindex, 1])

      }

      ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
      ## and SE
      res_Coeff_EBV_PEV_Rel_SE_omic1 <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic1_data,
        gen_name = gen_name,
        var_u = var_u,
        gid_name =gid_name)

    }

    if(ETA$ETA_element_name[aa]=="omic2_data"){

      Bb <- BGLR::readBinMat(BIN[aa])

      if(GS_model== "BRR"){
        var_u <- scan(varB_files[aa],
                         what = numeric(),
                         sep = "\n", quiet =TRUE)
        var_u <- mean(var_u[posindex])
        #var_u_Se = sd(var_u)/sqrt(length(var_u))
        ### The output for the variance was recently change in BGLR, this was fixed
        ## accordingly.
        ## It thus expedient to always crosscheck for any new update with the
        ## dependency package
      } else {
        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
        var_u <- utils::read.table(varB_files[aa])

        } else {
          var_u <- utils::read.table(varB_files[aa], skip = 1)

        }
        var_u <- mean(as.data.frame(tidyr::separate_rows(var_u))[posindex, 1])
        #var_u_Se = sd(var_u)/sqrt(length(var_u))

      }

      ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
      ## and SE
      res_Coeff_EBV_PEV_Rel_SE_omic2 <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic2_data,
        gen_name = gen_name,
        var_u = var_u,
        gid_name =gid_name)
    }

    if(ETA$ETA_element_name[aa]=="omic3_data"){

      Bb <- BGLR::readBinMat(BIN[aa])

      omic3_ebv <- as.matrix(omic3_data)%*%t(Bb)

      if(GS_model== "BRR"){
        var_u <- scan(varB_files[aa],
                         what = numeric(),
                         sep = "\n", quiet =TRUE)
        var_u <- mean(var_u[posindex])

        ### The output for the variance was recently change in BGLR, this was fixed
        ## accordingly.
        ## It thus expedient to always crosscheck for any new update with the
        ## dependency package
      } else {
        if(GS_model=="BayesA" | GS_model=="BayesC" | GS_model=="lambda"){
        var_u <- utils::read.table(varB_files[aa])

        } else {

          var_u <- utils::read.table(varB_files[aa], skip = 1)
        }
        var_u <- mean(as.data.frame(tidyr::separate_rows(var_u))[posindex, 1])
        #var_u_Se = sd(var_u)/sqrt(length(var_u))
      }

      ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
      ## and SE
      res_Coeff_EBV_PEV_Rel_SE_omic3 <- Cal_Coeff_EBV_PEV_Rel_SE(
        beta = Bb,
        x_variable = omic3_data,
        gen_name = gen_name,
        var_u = var_u,
        gid_name =gid_name)
    }

  }
  #### start from here

  ### When you have just one M_matrix/X_matrix
} else {


  if(ETA$ETA_element_name[1]=="geno_data"){

    Bb <- BGLR::readBinMat(BIN)

    ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
    ## and SE
    res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
      beta = Bb,
      x_variable = geno_data,
      gen_name = gen_name,
      var_u = Var_U,
      gid_name =gid_name)
  }

  ##
  if(ETA$ETA_element_name[1]=="omic1_data"){

    Bb <- BGLR::readBinMat(BIN)

    ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
    ## and SE
    res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
      beta = Bb,
      x_variable = omic1_data,
      gen_name = gen_name,
      var_u = Var_U,
      gid_name =gid_name)

  }

  if(ETA$ETA_element_name[1]=="omic2_data"){

    Bb <- BGLR::readBinMat(BIN)

    ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
    ## and SE
    res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
      beta = Bb,
      x_variable = omic2_data,
      gen_name = gen_name,
      var_u = Var_U,
      gid_name =gid_name)
  }

  if(ETA$ETA_element_name[1]=="omic3_data"){

    Bb <- BGLR::readBinMat(BIN)

    ### Calculate the Coefficient, Estimated breeding value, PEV, Reliability
    ## and SE
    res_Coeff_EBV_PEV_Rel_SE <- Cal_Coeff_EBV_PEV_Rel_SE(
      beta = Bb,
      x_variable = omic3_data,
      gen_name = gen_name,
      var_u = Var_U,
      gid_name =gid_name)

  }

}

########

## Create output for predicted value and residual value.
## The residual value dataframe also contain predicted value for two reasons
#1) For ease of plotting
#2) When testing set is present in the real world it is expected to be
Predicted_value <- data.frame(name = NA,
                            Predicted_value = mod$model$yHat,
                            stringsAsFactors = FALSE)

colnames(Predicted_value)[1] <- gen_name

### Residual value is only estimable for response value without NA
tst <- which(is.na(mod$model$y))

if(length(tst)!=0){
  Residual_value <- data.frame(name = NA,
                               Predicted_value = mod$model$yHat[tst],
                               Residual_value = (mod$model$y[tst] - mod$model$yHat[tst]),
                               stringsAsFactors = FALSE)

} else {

  Residual_value <- data.frame(name = NA,
                              Predicted_value = mod$model$yHat,
                              Residual_value = (mod$model$y - mod$model$yHat),
                              stringsAsFactors = FALSE)

}


####
if(length(BIN)==1){

  if(!is.null(geno_data)){

    Predicted_value[, 1] <- rownames(geno_data)
    #PEV <- apply(g_ebv, 1, var)
    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                    PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                    Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(geno_data)[tst]
    }else {

      Residual_value[, 1] <- rownames(geno_data)
    }

  } else if(!is.null(omic1_data)){

    Predicted_value[, 1] <- rownames(omic1_data)

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                    PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                    Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_data)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_data)
    }


  } else if(!is.null(omic2_data)){

    Predicted_value[, 1] <- rownames(omic2_data)

    Predicted_value <- Predicted_value |>
      dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                    PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                    Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic2_data)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic2_data)
    }


  } else {

    if(!is.null(omic3_data)){

      Predicted_value[, 1] <- rownames(omic3_data)

      Predicted_value <- Predicted_value |>
        dplyr::mutate(Std_error = sqrt(res_Coeff_EBV_PEV_Rel_SE$PEV),
                      PEV = res_Coeff_EBV_PEV_Rel_SE$PEV,
                      Reliability = res_Coeff_EBV_PEV_Rel_SE$Reliability)

      if(length(tst)!=0){
        Residual_value[, 1] <- rownames(omic3_data)[tst]
      } else {

        Residual_value[, 1] <- rownames(omic3_data)
      }


    }

  }

} else {

  if(length(BIN)>1){

    if(!is.null(geno_data)){
    Predicted_value[, 1] <- rownames(geno_data)
    } else if (!is.null(omic1_data)){
      Predicted_value[, 1] <- rownames(omic1_data)
    } else if (!is.null(omic1_data)){
      Predicted_value[, 1] <- rownames(omic1_data)
    } else if (!is.null(omic2_data)){
      Predicted_value[, 1] <- rownames(omic2_data)
    } else {
      if (!is.null(omic3_data)){
        Predicted_value[, 1] <- rownames(omic3_data)
      }

    }


  }



}








#################

if(!is.null(geno_data) & ((is.null(omic1_data) &  is.null(omic2_data)) & is.null(omic3_data))){

  Variance_components <- Bayes_variance_components(Var_U = Var_U,
                                                   Var_E = Var_E,
                                                   genomic_h2 = genomic_h2,
                                                   Var_U_Se = Var_U_Se,
                                                   Var_E_Se = Var_E_Se,
                                                   genomic_h2_Se = genomic_h2_Se)


  ### Combined all results into list
Res <-  list(
           res_Coeff_EBV_PEV_Rel_SE$Coefficient,
           res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
           Predicted_value,
           Variance_components,
           list(geno_clean_ready = geno_data),## This is to make it compatible in output format for when the omics is more than one
           mod$model$mu
           )

## Add attribute/ name to the list
names(Res) <- c(
                "Coefficients",
                "Estimated_breeding_value",
                "Predicted_value",
                "Variance_components",
                "M_matrix_model_ready",
                "Intercept")

} else if(is.null(geno_data) & ((!is.null(omic1_data) &  is.null(omic2_data)) & is.null(omic3_data))){

  Variance_components <- Bayes_variance_components(Var_U = Var_U,
                                                   Var_E = Var_E,
                                                   genomic_h2 = genomic_h2,
                                                   Var_U_Se = Var_U_Se,
                                                   Var_E_Se = Var_E_Se,
                                                   genomic_h2_Se = genomic_h2_Se)


  ### Combined all results into list

  Res <-  list(
              res_Coeff_EBV_PEV_Rel_SE$Coefficient,
              res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
              Predicted_value,
              Variance_components,
              list(omic_clean_ready = omic1_data),## This is to make it compatible in output format for when the omics is more than one
              mod$model$mu)

  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "Intercept")


  ## Add attribute/ name to the list

  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


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

} else if(is.null(geno_data) & ((is.null(omic1_data) &  !is.null(omic2_data)) & is.null(omic3_data))){

  Variance_components <- Bayes_variance_components(Var_U = Var_U,
                                                   Var_E = Var_E,
                                                   genomic_h2 = genomic_h2,
                                                   Var_U_Se = Var_U_Se,
                                                   Var_E_Se = Var_E_Se,
                                                   genomic_h2_Se = genomic_h2_Se)


  Res <-  list(
             res_Coeff_EBV_PEV_Rel_SE$Coefficient,
             res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
             Predicted_value,
             Variance_components,
             list(omic_clean_ready = omic2_data),## This is to make it compatible in output format for when the omics is more than one
             mod$model$mu
             )

  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "Intercept")


  ## Add attribute/ name to the list

  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


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

} else if(is.null(geno_data) & ((is.null(omic1_data) &  is.null(omic2_data)) & !is.null(omic3_data))){

  Variance_components <- Bayes_variance_components(Var_U = Var_U,
                                                   Var_E = Var_E,
                                                   genomic_h2 = genomic_h2,
                                                   Var_U_Se = Var_U_Se,
                                                   Var_E_Se = Var_E_Se,
                                                   genomic_h2_Se = genomic_h2_Se)


  Res <-  list(
             res_Coeff_EBV_PEV_Rel_SE$Coefficient,
             res_Coeff_EBV_PEV_Rel_SE$Estimated_breeding_value,
             Predicted_value,
             Variance_components,
             list(omic_clean_ready = omic3_data), ## This is to make it compatible in output format for when the omics is more than one
             mod$model$mu)

  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "Intercept")

  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


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

} else if(!is.null(geno_data) & ((!is.null(omic1_data) &  is.null(omic2_data)) & is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV,
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


  ###
  if(length(tst)!=0){
  Residual_value[, 1] <- rownames(geno_data)[tst]
}else {

  Residual_value[, 1] <- rownames(geno_data)
}

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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
             list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                  Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient),
             list(Genomic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                  Omic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value),
             sum_EBV,
             Predicted_value,
             Variance_components,
             list(Geno_model_ready = geno_data,
                  Omic_model_ready = omic1_data
                                         ),
             mod$model$mu)

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

    print_lable <- unlist(print_lable)


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

} else if(!is.null(geno_data) & ((is.null(omic1_data) &  !is.null(omic2_data)) & is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
                        stringsAsFactors = FALSE)

  colnames(sum_EBV)[1] <- gen_name

  ##Estimate of reliabilities
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(geno_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(geno_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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
                list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                     Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),

               list(Genomic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                    Omic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
              sum_EBV,
              Predicted_value,
              Variance_components,
              list(Geno_model_ready = geno_data,
                   Omic_model_ready = omic2_data
                 ),
              mod$model$mu)

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

    print_lable <- unlist(print_lable)


    if(length(print_lable)==1){
      rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
      names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
    } else{

      if(length(print_lable)>1){
        message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }
  #####

} else if(!is.null(geno_data) & ((is.null(omic1_data) &  is.null(omic2_data)) & !is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(geno_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(geno_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components(Var_U_1 = Var_U_1,
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
               list(Geno_coeff = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                    Omic_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),

              list(Genomic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                    Omic_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),

               sum_EBV,
               Predicted_value,
               Variance_components,
               list(Geno_model_ready = geno_data,
                    Omic_model_ready = omic3_data
                                ),
               mod$model$mu)

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

    print_lable <- unlist(print_lable)


    if(length(print_lable)==1){
      rownames(Res$Variance_components)[1:2] = c('Genomic_variance', paste(print_lable, "variance", sep="_"))
      names(Res$Coefficients) = c("Geno_coefficient",  paste(print_lable, "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable, "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable, "model_ready", sep="_"))
    } else {
      if(length(print_lable)>1){
        message(insight::print_color(paste(msg,paste("More than one Omics lable were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }



} else if(is.null(geno_data) & ((!is.null(omic1_data) &  !is.null(omic2_data)) & is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(omic1_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(omic1_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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
               list(Omic1_coeff =  res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                    Omic2_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),
               list(Omic1_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                    Omic2_estimated_BV = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
               sum_EBV,
               Predicted_value,
               Variance_components,

               list(omic1_data_model_ready= omic1_data,
                    omic2_data_model_ready= omic2_data),
               mod$model$mu)

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

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:2] = c(paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c(paste(print_lable[1], "coefficient", sep="_"),  paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c(paste(print_lable[1], "estimated_breeding_value", sep="_"),  paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c(paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c(paste(print_lable[1], "model_ready", sep="_"),  paste(print_lable[2], "model_ready", sep="_"))
    }else {
      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than two Omics lable were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

} else if(is.null(geno_data) & ((is.null(omic1_data) &  !is.null(omic2_data)) & !is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(omic2_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(omic2_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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
               list(Omic2_coeff = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                    Omic3_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
               list(Omic2_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                    Omic3_estimated_BV = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
               sum_EBV,
               Predicted_value,
               Variance_components,
               list(omic2_data_model_ready= omic2_data,
                    omic3_data_model_ready= omic3_data),
               mod$model$mu)

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

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:2] = c(paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c(paste(print_lable[1], "coefficient", sep="_"),  paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c(paste(print_lable[1], "estimated_breeding_value", sep="_"),  paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c(paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c(paste(print_lable[1], "model_ready", sep="_"),  paste(print_lable[2], "model_ready", sep="_"))
    } else {
      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than two Omic lables were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

} else if(is.null(geno_data) & ((!is.null(omic1_data) &  is.null(omic2_data)) & !is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(omic1_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(omic1_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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
               list(Omic1_coeff = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                    Omic3_coeff = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
               list(Omic1_estimated_BV= res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                    Omic3_estimated_BV = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
               sum_EBV,
               Predicted_value,
               Variance_components,
               list(omic1_data_model_ready= omic1_data,
                    omic2_data_model_ready= omic3_data),
               mod$model$mu)

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

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:2] = c(paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c(paste(print_lable[1], "coefficient", sep="_"),  paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c(paste(print_lable[1], "estimated_breeding_value", sep="_"),  paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c(paste(print_lable[1], "variance", sep="_"),  paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c(paste(print_lable[1], "model_ready", sep="_"),  paste(print_lable[2], "model_ready", sep="_"))
    } else {

      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than two Omics lable were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

} else if(!is.null(geno_data) & ((is.null(omic1_data) &  !is.null(omic2_data)) & !is.null(omic3_data))){


  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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
  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(geno_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(geno_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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

Res <-list(
            list(
               coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
               coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
               coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
            list(
               EBV_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
               EBV_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
               EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
           sum_EBV,
           Predicted_value,
           Variance_components,
          list(geno_model_ready = geno_data,
               omic1_model_ready = omic2_data,
               omic2_model_ready = omic3_data
               ),
          mod$model$mu)


  names(Res) <- c(
             "Coefficients",
             "Estimated_breeding_value",
             "Total_estimated_breeding_value",
             "Predicted_value",
             "Variance_components",
             "M_matrix_model_ready",
             "intercept")

  ###
  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:3] = c('Genomic_variance', paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c("Geno_coefficient", paste(print_lable[1], "coefficient", sep="_"), paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value",  paste(print_lable[1], "estimated_breeding_value", sep="_"), paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c("Geno_model_ready", paste(print_lable[1], "model_ready", sep="_"), paste(print_lable[2], "model_ready", sep="_"))
    } else {

      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than three Omic lables were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

  #####

} else if(!is.null(geno_data) & ((!is.null(omic1_data) &  is.null(omic2_data)) & !is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(geno_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(geno_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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


  Res <-list(
             list(
                   coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                   coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                   coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
            list(EBV_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                 EBV_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                 EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
            sum_EBV,
            Predicted_value,
            Variance_components,
            list(geno_model_ready = geno_data,
                 omic1_model_ready = omic1_data,
                 omic2_model_ready = omic3_data),
          mod$model$mu)


  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Total_estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "intercept")

  ###
  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:3] = c('Genomic_variance', paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c("Geno_coefficient", paste(print_lable[1], "coefficient", sep="_"), paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value", paste(print_lable[1], "estimated_breeding_value", sep="_"), paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c("Genomic_variance", paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable[2], "model_ready", sep="_"), paste(print_lable[3], "model_ready", sep="_"))
    } else {

      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than three Omic lables were provided. Default name was applied.")), "blue"))
      }
    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

} else if(!is.null(geno_data) & ((!is.null(omic1_data) &  !is.null(omic2_data)) & is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)

  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(geno_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(geno_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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

  Res <-list(
              list(
                    coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                    coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                    coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient),
              list(
                   EBV_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                   EBV_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                   EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value),
              sum_EBV,
              Predicted_value,
              Variance_components,

             list(
                   geno_model_ready = geno_data,
                   omic1_model_ready = omic1_data,
                   omic2_model_ready = omic2_data
                   ),
            mod$model$mu)


  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Total_estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "intercept")

  ###
  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)


    if(length(print_lable)==2){
      rownames(Res$Variance_components)[1:3] = c('Genomic_variance', paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$Coefficients) = c("Geno_coefficient", paste(print_lable[1], "coefficient", sep="_"), paste(print_lable[2], "coefficient", sep="_"))
      names(Res$Estimated_breeding_value) = c("Genomic_estimated_breeding_value", paste(print_lable[1], "estimated_breeding_value", sep="_"), paste(print_lable[2], "estimated_breeding_value", sep="_"))
      #names(Res$Genetic_variance) = c("Genomic_variance", paste(print_lable[1], "variance", sep="_"), paste(print_lable[2], "variance", sep="_"))
      names(Res$M_matrix_model_ready) = c("Geno_model_ready",  paste(print_lable[2], "model_ready", sep="_"), paste(print_lable[3], "model_ready", sep="_"))
    } else {

      if(length(print_lable)>2 | length(print_lable)<2){
        message(insight::print_color(paste(msg,paste("More than three Omic lables were provided. Default name was applied.")), "blue"))
      }

    }

  } else {

      message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }

  ##

} else if(is.null(geno_data) & ((!is.null(omic1_data) &  !is.null(omic2_data)) & !is.null(omic3_data))){

  g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
  PEV <- apply(g_omic_ebv, 1, var)
  ##Estimate of reliabilities
  Reliability <- 1 - (PEV / Var_U)
  Reliability <- ifelse(Reliability<0, "Alias", Reliability)


  sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                        EBV = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

  ###
  if(length(tst)!=0){
    Residual_value[, 1] <- rownames(omic1_data)[tst]
  }else {

    Residual_value[, 1] <- rownames(omic1_data)
  }

  ### Ends

  Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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



  Res <-list(
              list(
                    coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                    coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                    coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
             list(
                   EBV_1 = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                   EBV_2 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                   EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
            sum_EBV,
            Predicted_value,
            Variance_components,

           list(omic1_model_ready = omic1_data,
               omic2_model_ready = omic2_data,
               omic3_model_ready = omic3_data),
           mod$model$mu)


  names(Res) <- c(
    "Coefficients",
    "Estimated_breeding_value",
    "Total_estimated_breeding_value",
    "Predicted_value",
    "Variance_components",
    "M_matrix_model_ready",
    "intercept")

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

  ##

} else{

  if(!is.null(geno_data) & ((!is.null(omic1_data) &  !is.null(omic2_data)) & !is.null(omic3_data))){

    g_omic_ebv <- res_Coeff_EBV_PEV_Rel_SE_gen$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic1$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic2$Posterior + res_Coeff_EBV_PEV_Rel_SE_omic3$Posterior
    PEV <- apply(g_omic_ebv, 1, var)
    ##Estimate of reliability
    Reliability <- 1 - (PEV / Var_U)
    Reliability <- ifelse(Reliability<0, "Alias", Reliability)

    sum_EBV <- data.frame(name = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value[, 1],
                          EBV = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value$EBV + res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value$EBV+ res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value$EBV,
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

    ###
    if(length(tst)!=0){
      Residual_value[, 1] <- rownames(omic1_data)[tst]
    }else {

      Residual_value[, 1] <- rownames(omic1_data)
    }

    ### Ends

    Variance_components <- Bayes_variance_components (Var_U_1 = Var_U_1,
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


    Res <-list(
                list(
                      coefficients_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Coefficient,
                      coefficients_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Coefficient,
                      coefficients_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Coefficient,
                      coefficients_4 = res_Coeff_EBV_PEV_Rel_SE_omic3$Coefficient),
                list(
                     EBV_1 = res_Coeff_EBV_PEV_Rel_SE_gen$Estimated_breeding_value,
                     EBV_2 = res_Coeff_EBV_PEV_Rel_SE_omic1$Estimated_breeding_value,
                     EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic2$Estimated_breeding_value,
                     EBV_3 = res_Coeff_EBV_PEV_Rel_SE_omic3$Estimated_breeding_value),
                sum_EBV,
                Predicted_value,
                Variance_components,

                list(
                     geno_model_ready = geno_data,
                     omic1_model_ready = omic1_data,
                     omic2_model_ready = omic2_data,
                     omic3_model_ready = omic3_data
                   ),
                mod$model$mu)


    names(Res) <- c(
      "Coefficients",
      "Estimated_breeding_value",
      "Total_estimated_breeding_value",
      "Predicted_value",
      "Variance_components",
      "M_matrix_model_ready",
      "intercept")

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

    ##

  }

  ####

  }



 #}
### remove the generated output files from the working directory
unlink(mod$output_files_names)
if(exists("res_Coeff_EBV_PEV_Rel_SE_omic1")){
  rm(res_Coeff_EBV_PEV_Rel_SE_omic1)
}

if(exists("res_Coeff_EBV_PEV_Rel_SE_omic2")){
  rm(res_Coeff_EBV_PEV_Rel_SE_omic2)
}

if(exists("res_Coeff_EBV_PEV_Rel_SE_omic3")){
  rm(res_Coeff_EBV_PEV_Rel_SE_omic3)
}

if(exists("res_Coeff_EBV_PEV_Rel_SE_gen")){
  rm(res_Coeff_EBV_PEV_Rel_SE_gen)
}

if(exists("res_Coeff_EBV_PEV_Rel_SE")){
  rm(res_Coeff_EBV_PEV_Rel_SE)
}

return(Res)
}





