
#' Title
#'
#' @param mod
#' @param ETA
#' @param ...
#' @param gen_name
#' @param gkernel
#' @param gmatrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param pheno_data
#' @param heter_groups
#'
#' @return
#' @export
#'
#' @examples
mod_output_bayes_RKHS <- function(mod=NULL,
                             ETA=NULL,
                             gen_name=NULL,
                             gkernel=NULL,
                             gmatrix = NULL,
                             omic1_kernel=NULL,
                             omic2_kernel=NULL,
                             omic3_kernel=NULL,
                             pheno_data = NULL,
                             heter_groups = NULL,
                             ...){

  # sik$ETA_element_name
  # TT = DT$output_files_names
  #### When the genotype are present in more than one environment/location
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    if(!is.null(gmatrix)){
      gmatrix <- Zg%*%gmatrix%*%t(Zg)

      gmatrix <- grm_kernel_precheck(gmatrix)
    }

    if(!is.null(gkernel)){
      gkernel <- Zg%*%gkernel%*%t(Zg)

      gkernel <- grm_kernel_precheck(gkernel)
    }

    if(!is.null(omic1_kernel)){
      omic1_kernel <- Zg%*%omic1_kernel%*%t(Zg)

      omic1_kernel <- grm_kernel_precheck(omic1_kernel)

    }

    if(!is.null(omic2_kernel)){
      omic2_kernel <- Zg%*%omic2_kernel%*%t(Zg)

      omic2_kernel <- grm_kernel_precheck(omic2_kernel)
    }

    if(!is.null(omic3_kernel)){
      omic3_kernel <- Zg%*%omic3_kernel%*%t(Zg)

      omic3_kernel <- grm_kernel_precheck(omic3_kernel)

    }

    # if(!is.null(heter_groups)){
    #   ZE <- model.matrix(~factor(pheno_data[,heter_groups])-1)
    #   ZEZE<-tcrossprod(ZE)
    #
    # }
    ENV = pheno_data[,heter_groups]

    predict_value = data.frame(name = pheno_data[,gen_name], Env = ENV, predict_value = mod$model$yHat)
    names(predict_value)[1] = gen_name
    } else {

    predict_value = data.frame(name = pheno_data[,gen_name], predict_value = mod$model$yHat)
    names(predict_value)[1] = gen_name
     }


  BIN = mod$output_files_names[grepl("_varU.dat", mod$output_files_names)]


  mu <- mean(scan(mod$output_files_names[grepl("mu.dat", mod$output_files_names)],
                     what = numeric(),
                     sep = "\n"))


  ### Extract Error variance
  Var_E <- mean(scan(mod$output_files_names[grepl("varE.dat", mod$output_files_names)],
                     what = numeric(),
                     sep = "\n"))

  ### If more than one M-matrix is provided. Effect will have it own coefficient
  ### These lines of code exttract the genomic variance and the error term
  if(length(BIN)>1){

    varB_files <- mod$output_files_names[grepl("_varU.dat", mod$output_files_names)]


    if(length(BIN)==4){
      Var_U_1= mean(scan(varB_files[1],
                         what = numeric(),
                         sep = "\n"))

      Var_U_2= mean(scan(varB_files[2],
                         what = numeric(),
                         sep = "\n"))

      Var_U_3= mean(scan(varB_files[3],
                         what = numeric(),
                         sep = "\n"))

      Var_U_4= mean(scan(varB_files[4],
                         what = numeric(),
                         sep = "\n"))

      genomic_h2 = (Var_U_1 + Var_U_2 + Var_U_3 + Var_U_4)/((Var_U_1 + Var_U_2 + Var_U_3 + Var_U_4)+Var_E)

      Var_U = Var_U_1 + Var_U_2 + Var_U_3 + Var_U_4
    }

    if(length(BIN)==3){
      Var_U_1= mean(scan(varB_files[1],
                         what = numeric(),
                         sep = "\n"))

      Var_U_2= mean(scan(varB_files[2],
                         what = numeric(),
                         sep = "\n"))

      Var_U_3= mean(scan(varB_files[3],
                         what = numeric(),
                         sep = "\n"))

      genomic_h2 = (Var_U_1 + Var_U_2 + Var_U_3)/((Var_U_1 + Var_U_2 + Var_U_3)+Var_E)

      Var_U = Var_U_1 + Var_U_2 + Var_U_3

    }

    if(length(BIN)==2){
      Var_U_1= mean(scan(varB_files[1],
                         what = numeric(),
                         sep = "\n"))

      Var_U_2= mean(scan(varB_files[2],
                         what = numeric(),
                         sep = "\n"))


      genomic_h2 = (Var_U_1 + Var_U_2)/((Var_U_1 + Var_U_2)+Var_E)

      Var_U = Var_U_1 + Var_U_2
    }

  } else {

    Var_U= mean(scan(mod$output_files_names[grepl("_varU.dat", mod$output_files_names)],
                     what = numeric(),
                     sep = "\n"))

    Var_E <- mean(scan(mod$output_files_names[grepl("varE.dat", mod$output_files_names)],
                       what = numeric(),
                       sep = "\n"))


    genomic_h2 = Var_U/(Var_U+Var_E)

    # PEV = apply(ebv, 1, var)
    #
    # ##Estimate of reliabilities
    # reliability = 1 - PEV / Var_U

  }


  #### When there is more than one M matrix

  if(length(BIN)>1){
    for (aa in 1:length(BIN)) {

      if(ETA$ETA_element_name[aa]=="gkernel" | ETA$ETA_element_name[aa]=="gmatrix"){

        if(ETA$ETA_element_name[aa]=="gkernel"){
          g_use = gkernel

          rm(gkernel)

        } else{

          if(ETA$ETA_element_name[aa]=="gmatrix"){
            g_use = gmatrix

            rm(gmatrix)
          }


        }

        g_ebv= mod$model$yHat - mod$model$mu

        # Posterior means of marker effects/coefficient
        coeffRaw = solve(t(g_use)*g_use)*t(g_use)*g_ebv
        coeff <- colMeans(coeffRaw)
        # Genomic estimated breeding values
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        GEBV <- data.frame(names = rownames(g_use), Env = ENV, GEBV=g_ebv)
        } else {

          GEBV <- data.frame(names = rownames(g_use), GEBV=g_ebv)
        }

        colnames(GEBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV = apply(coeffRaw, 1, var)

        GEBV$PEV <-PEV

        var_u= mean(scan(varB_files[aa],
                         what = numeric(),
                         sep = "\n"))
        ##Estimate of reliability
        GEBV$reliability = 1 - PEV /var_u

        rm(PEV, coeffRaw)
      }

      ##
      if(ETA$ETA_element_name[aa]=="omic1_kernel"){

        omic1_ebv= mod$model$yHat - mod$model$mu

        # Posterior means of marker effects/coefficient
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*omic1_ebv
        coeff_omic1 <- colMeans(coeffRaw_omic1)

        # Genomic estimated breeding values
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        omic1_EBV <- data.frame(names = rownames(omic1_kernel), Env = ENV, omic1_EBV=omic1_ebv)
        } else{

          omic1_EBV <- data.frame(names = rownames(omic1_kernel), omic1_EBV=omic1_ebv)

        }
        colnames(omic1_EBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV = apply(coeffRaw_omic1, 1, var)

        omic1_EBV$PEV <- PEV

        ##Estimate of reliabilities
        omic1_EBV$reliability = 1 - PEV / Var_U_1

        rm(PEV, coeffRaw_omic1)
      }

      if(ETA$ETA_element_name[aa]=="omic2_kernel"){

        omic2_ebv= mod$model$yHat - mod$model$mu

        # Posterior means of marker effects/coefficient
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*omic2_ebv
        coeff_omic2 <- colMeans(coeffRaw_omic2)

        # Genomic estimated breeding values
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        omic2_EBV <- data.frame(names = rownames(omic2_kernel), Env = ENV, omic2_EBV=omic2_ebv)
        } else {
          omic2_EBV <- data.frame(names = rownames(omic2_kernel), omic2_EBV=omic2_ebv)

        }
        colnames(omic2_EBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV = apply(coeffRaw_omic2, 1, var)

        omic2_EBV$PEV <- PEV

        ##Estimate of reliabilities
        omic2_EBV$reliability = 1 - PEV / Var_U_1

        rm(PEV, coeffRaw_omic2)
      }

      if(ETA$ETA_element_name[aa]=="omic3_kernel"){

        omic3_ebv= mod$model$yHat - mod$model$mu

        # Posterior means of marker effects/coefficient
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*omic3_ebv
        coeff_omic3 <- colMeans(coeffRaw_omic3)

        # Genomic estimated breeding values
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        omic3_EBV <- data.frame(names = rownames(omic3_kernel), Env = ENV,  omic3_EBV=omic3_ebv)

        } else {
          omic3_EBV <- data.frame(names = rownames(omic3_kernel), omic3_EBV=omic3_ebv)

        }

        colnames(omic3_EBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV = apply(coeffRaw_omic3, 1, var)

        omic3_EBV$PEV <- PEV

        ##Estimate of reliabilities
        omic3_EBV$reliability = 1 - PEV / Var_U_1

        rm(PEV, coeffRaw_omic3)
      }

    }
    #### start from here

    ### When you have just one M_matrix/X_matrix
  } else {


      if(ETA$ETA_element_name[1]=="gkernel" | ETA$ETA_element_name[1]=="gmatrix"){

        if(ETA$ETA_element_name[1]=="gkernel"){
          g_use = gkernel

          rm(gkernel)

        } else{

          if(ETA$ETA_element_name[1]=="gmatrix"){
            g_use = gmatrix

            rm(gmatrix)
          }



        }

        g_ebv= mod$model$yHat - mod$model$mu

        # Posterior means of marker effects/coefficient
        coeffRaw = solve(t(g_use)*g_use)*t(g_use)*g_ebv
        coeff <- colMeans(coeffRaw)
        # Genomic estimated breeding values
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        GEBV <- data.frame(names = rownames(g_use), Env = ENV, GEBV=g_ebv)
        } else {
          GEBV <- data.frame(names = rownames(g_use), GEBV=g_ebv)

        }
        colnames(GEBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV = apply(coeffRaw, 1, var)

        GEBV$PEV <-PEV

        ##Estimate of reliability
        GEBV$reliability = 1 - PEV /Var_U

        rm(PEV, coeffRaw)
      }



    if(ETA$ETA_element_name[1]=="omic1_kernel"){

      omic1_ebv= mod$model$yHat - mod$model$mu

      # Posterior means of marker effects/coefficient
      coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*omic1_ebv
      coeff_omic1 <- colMeans(coeffRaw_omic1)

      # Genomic estimated breeding values
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      omic1_EBV <- data.frame(names = rownames(omic1_kernel), Env = ENV, omic1_EBV=omic1_ebv)

      } else {

        omic1_EBV <- data.frame(names = rownames(omic1_kernel), omic1_EBV=omic1_ebv)
      }
      colnames(omic1_EBV)[1] = gen_name
      #GEBV = data.frame(rowMeans(ebv))

      PEV = apply(coeffRaw_omic1, 1, var)

      omic1_EBV$PEV <- PEV

      ##Estimate of reliability
      omic1_EBV$reliability = 1 - PEV / Var_U

      rm(PEV, coeffRaw_omic1)
    }

    if(ETA$ETA_element_name[1]=="omic2_kernel"){

      omic2_ebv= mod$model$yHat - mod$model$mu

      # Posterior means of marker effects/coefficient
      coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*omic2_ebv
      coeff_omic2 <- colMeans(coeffRaw_omic2)

      # Genomic estimated breeding values
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      omic2_EBV <- data.frame(names = rownames(omic2_kernel), Env = ENV, omic2_EBV=omic2_ebv)
      } else {

        omic2_EBV <- data.frame(names = rownames(omic2_kernel), omic2_EBV=omic2_ebv)
      }
      colnames(omic2_EBV)[1] = gen_name
      #GEBV = data.frame(rowMeans(ebv))

      PEV = apply(coeffRaw_omic2, 1, var)

      omic2_EBV$PEV <- PEV

      ##Estimate of reliabilities
      omic2_EBV$reliability = 1 - PEV / Var_U

      rm(PEV, coeffRaw_omic2)
    }

    if(ETA$ETA_element_name[1]=="omic3_kernel"){

      omic3_ebv= mod$model$yHat - mod$model$mu

      # Posterior means of marker effects/coefficient
      coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*omic3_ebv
      coeff_omic3 <- colMeans(coeffRaw_omic3)

      # Genomic estimated breeding values
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
      omic3_EBV <- data.frame(names = rownames(omic3_kernel), Env = ENV,  omic3_EBV=omic3_ebv)
      } else {
        omic3_EBV <- data.frame(names = rownames(omic3_kernel), omic3_EBV=omic3_ebv)

      }
      colnames(omic3_EBV)[1] = gen_name
      #GEBV = data.frame(rowMeans(ebv))

      PEV = apply(coeffRaw_omic3, 1, var)

      omic3_EBV$PEV <- PEV

      ##Estimate of reliabilities
      omic3_EBV$reliability = 1 - PEV / Var_U

      rm(PEV, coeffRaw_omic3)
    }

  }

  #if((!is.null(gkernel) | !is.null(gmatrix)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
  if(exists("g_use") & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
    Res <-  list(coefficients = coeff,
                 GEBV = GEBV,
                 predict_value = predict_value,
                 genomic_variance = Var_U,
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 g_model_ready = g_use,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients",
                    "GEBV",
                    "predict_value",
                    "genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "g_model_ready",
                    "intercept")

   } else if(!exists("g_use") & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#  } else if((is.null(gkernel) & is.null(gmatrix)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients = coeff_omic1,
                 EBV = omic1_EBV,
                 predict_value = predict_value,
                 genomic_variance = Var_U,
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic1_kernel,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients",
                    "EBV",
                    "predict_value",
                    "genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "intercept")

  } else if(!exists("g_use") & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
 # } else if((is.null(gkernel) & is.null(gmatrix)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients = coeff_omic2,
                 EBV = omic2_EBV,
                 predict_value = predict_value,
                 genomic_variance = Var_U,
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic2_kernel,
                 mu = mod$model$mu)


    names(Res) <- c("coefficients",
                    "EBV",
                    "predict_value",
                    "genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "intercept")

  } else if(!exists("g_use") & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((is.null(gkernel) & is.null(gmatrix)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients = coeff_omic3,
                 EBV = omic3_EBV,
                 predict_value = predict_value,
                 genomic_variance = Var_U,
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic3_kernel,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients",
                    "EBV",
                    "predict_value",
                    "genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "intercept")

  } else if(exists("g_use") & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
  #} else if((!is.null(gkernel) | !is.null(gmatrix)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic1,
                 EBV_1 = GEBV,
                 EBV_2 = omic1_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = GEBV$GEBV+omic1_EBV$omic1_EBV),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 total_genomic_variance = sum(Var_U_1, Var_U_2),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic1_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "g_model_ready",
                    "intercept")

  } else if(exists("g_use") & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#  } else if((!is.null(gkernel) | !is.null(gmatrix)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic2,
                 EBV_1 = GEBV,
                 EBV_2 = omic2_EBV,
                 sum_EBV = data.frame(name = omic2_EBV[, 1], EBV = GEBV$GEBV+omic2_EBV$omic2_EBV),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 total_genomic_variance = sum(Var_U_1, Var_U_2),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic2_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "g_model_ready",
                    "intercept")

  } else if((exists("g_use")) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((!is.null(gkernel) | !is.null(gmatrix)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic3,
                 EBV_1 = GEBV,
                 EBV_2 = omic3_EBV,
                 sum_EBV = data.frame(name = omic3_EBV[, 1], EBV = GEBV$GEBV+omic3_EBV$omic3_EBV),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 total_genomic_variance = sum(Var_U_1, Var_U_2),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic_model_ready = omic3_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M_matrix_model_ready",
                    "g_model_ready",
                    "intercept")

  } else if(!exists("g_use") & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
 # } else if((is.null(gkernel) & is.null(gmatrix)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff_omic1,
                 coefficients_2 = coeff_omic2,
                 EBV_1 = omic1_EBV,
                 EBV_2 = omic2_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = omic1_EBV$omic1_EBV + omic2_EBV$omic2_EBV),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 total_genomic_variance = sum(Var_U_1, Var_U_2),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic1_kernel,
                 omic2_model_ready = omic2_kernel,
                 mu = mod$model$mu)


    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M1_model_ready",
                    "M2_model_ready",
                    "intercept")

  } else if(!exists("g_use") & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((is.null(gkernel) & is.null(gmatrix)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <- list(coefficients_1 = coeff_omic2,
                coefficients_2 = coeff_omic3,
                EBV_1 = omic2_EBV,
                EBV_2 = omic3_EBV,
                genomic_variance_1 = Var_U_1,
                genomic_variance_2 = Var_U_2,
                sum_EBV = data.frame(name = omic2_EBV[, 1], EBV = omic2_EBV$omic2_EBV + omic3_EBV$omic3_EBV),
                predict_value = predict_value,
                total_genomic_variance = sum(Var_U_1, Var_U_2),
                residual_error = Var_E,
                genomic_heritability = genomic_h2,
                omic1_model_ready = omic2_kernel,
                omic2_model_ready = omic3_kernel,
                mu = mod$model$mu)

    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M1_model_ready",
                    "M2_model_ready",
                    "intercept")

  } else if(!exists("g_use") & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((is.null(gkernel) & is.null(gmatrix)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff_omic1,
                 coefficients_2 = coeff_omic3,
                 EBV_1 = omic1_EBV,
                 EBV_2 = omic3_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = omic1_EBV$omic1_EBV + omic3_EBV$omic3_EBV),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 total_genomic_variance = sum(Var_U_1, Var_U_2),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic1_kernel,
                 omic2_model_ready = omic3_kernel,
                 mu = mod$model$mu)

    names(Res) <- c("coefficients_1",
                    "coefficients_2",
                    "EBV_1",
                    "EBV_2",
                    "total_EBV",
                    "predict_value",
                    "genomic_variance_1",
                    "genomic_variance_2",
                    "total_genomic_variance",
                    "residual_error",
                    "genomic_heritability",
                    "M1_model_ready",
                    "M2_model_ready",
                    "intercept")

  } else if(exists("g_use") & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((!is.null(gkernel) | !is.null(gmatrix)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic2,
                 coefficients_3 = coeff_omic3,
                 EBV_1 = GEBV,
                 EBV_2 = omic2_EBV,
                 EBV_3 = omic3_EBV,
                 sum_EBV = data.frame(name = omic2_EBV[, 1], EBV = (GEBV$GEBV + omic2_EBV$omic2_EBV + omic3_EBV$omic3_EBV)),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 genomic_variance_3 = Var_U_3,
                 total_genomic_variance = sum(Var_U_1, Var_U_2, Var_U_3),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic2_kernel,
                 omic2_model_ready = omic3_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)


    names(Res) <-  c("coefficients_1",
                     "coefficients_2",
                     "coefficients_3",
                     "EBV_1",
                     "EBV_2",
                     "EBV_3",
                     "total_EBV",
                     "predict_value",
                     "genomic_variance_1",
                     "genomic_variance_2",
                     "genomic_variance_3",
                     "total_genomic_variance",
                     "residual_error",
                     "genomic_heritability",
                     "M1_model_ready",
                     "M2_model_ready",
                     "g_model_ready",
                     "intercept")

    #####

  } else if(exists("g_use") & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((!is.null(gkernel) | !is.null(gmatrix)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic1,
                 coefficients_3 = coeff_omic3,
                 EBV_1 = GEBV,
                 EBV_2 = omic1_EBV,
                 EBV_3 = omic3_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = (GEBV$GEBV + omic1_EBV$omic1_EBV + omic3_EBV$omic3_EBV)),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 genomic_variance_3 = Var_U_3,
                 total_genomic_variance = sum(Var_U_1, Var_U_2, Var_U_3),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic1_kernel,
                 omic2_model_ready = omic3_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)


    names(Res) <-  c("coefficients_1",
                     "coefficients_2",
                     "coefficients_3",
                     "EBV_1",
                     "EBV_2",
                     "EBV_3",
                     "total_EBV",
                     "predict_value",
                     "genomic_variance_1",
                     "genomic_variance_2",
                     "genomic_variance_3",
                     "total_genomic_variance",
                     "residual_error",
                     "genomic_heritability",
                     "M1_model_ready",
                     "M2_model_ready",
                     "g_model_ready",
                     "intercept")

  } else if(exists("g_use") & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#  } else if((!is.null(gkernel) | !is.null(gmatrix)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff,
                 coefficients_2 = coeff_omic1,
                 coefficients_3 = coeff_omic2,
                 EBV_1 = GEBV,
                 EBV_2 = omic1_EBV,
                 EBV_3 = omic2_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = (GEBV$GEBV + omic1_EBV$omic1_EBV + omic2_EBV$omic2_EBV)),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 genomic_variance_3 = Var_U_3,
                 total_genomic_variance = sum(Var_U_1, Var_U_2, Var_U_3),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic1_kernel,
                 omic2_model_ready = omic2_kernel,
                 g_model_ready = g_use,
                 mu = mod$model$mu)

    names(Res) <-  c("coefficients_1",
                     "coefficients_2",
                     "coefficients_3",
                     "EBV_1",
                     "EBV_2",
                     "EBV_3",
                     "total_EBV",
                     "predict_value",
                     "genomic_variance_1",
                     "genomic_variance_2",
                     "genomic_variance_3",
                     "total_genomic_variance",
                     "residual_error",
                     "genomic_heritability",
                     "M1_model_ready",
                     "M2_model_ready",
                     "g_model_ready",
                     "intercept")

  } else if(!exists("g_use") & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#  } else if((is.null(gkernel) & is.null(gmatrix)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Res <-  list(coefficients_1 = coeff_omic1,
                 coefficients_2 = coeff_omic2,
                 coefficients_3 = coeff_omic3,
                 EBV_1 = omic1_EBV,
                 EBV_2 = omic2_EBV,
                 EBV_3 = omic3_EBV,
                 sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = (omic1_EBV$omic1_EBV + omic2_EBV$omic2_EBV + omic3_EBV$omic3_EBV)),
                 predict_value = predict_value,
                 genomic_variance_1 = Var_U_1,
                 genomic_variance_2 = Var_U_2,
                 genomic_variance_3 = Var_U_3,
                 total_genomic_variance = sum(Var_U_1, Var_U_2, Var_U_3),
                 residual_error = Var_E,
                 genomic_heritability = genomic_h2,
                 omic1_model_ready = omic1_kernel,
                 omic2_model_ready= omic2_kernel,
                 omic3_model_ready = omic3_kernel,
                 mu = mod$model$mu)

    names(Res) <-  c("coefficients_1",
                     "coefficients_2",
                     "coefficients_3",
                     "EBV_1",
                     "EBV_2",
                     "EBV_3",
                     "total_EBV",
                     "predict_value",
                     "genomic_variance_1",
                     "genomic_variance_2",
                     "genomic_variance_3",
                     "total_genomic_variance",
                     "residual_error",
                     "genomic_heritability",
                     "M1_model_ready",
                     "M2_model_ready",
                     "M3_model_ready",
                     "intercept")

  } else{

    if(exists("g_use") & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
  #   if((!is.null(gkernel) | !is.null(gmatrix)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      Res <-  list(coefficients_1 = coeff,
                   coefficients_2 = coeff_omic1,
                   coefficients_3 = coeff_omic2,
                   coefficients_4 = coeff_omic3,
                   EBV_1 = GEBV,
                   EBV_2 = omic1_EBV,
                   EBV_3 = omic2_EBV,
                   EBV_4 = omic3_EBV,
                   sum_EBV = data.frame(name = omic1_EBV[, 1], EBV = (GEBV$GEBV + omic1_EBV$omic1_EBV + omic2_EBV$omic2_EBV + omic3_EBV$omic3_EBV)),
                   predict_value = predict_value,
                   genomic_variance_1 = Var_U_1,
                   genomic_variance_2 = Var_U_2,
                   genomic_variance_3 = Var_U_3,
                   genomic_variance_4 = Var_U_4,
                   total_genomic_variance = sum(Var_U_1, Var_U_2, Var_U_3, Var_U_4),
                   residual_error = Var_E,
                   genomic_heritability = genomic_h2,
                   omic1_model_ready = omic1_kernel,
                   omic2_model_ready = omic2_kernel,
                   omic3_model_ready = omic3_kernel,
                   g_model_ready = g_use,
                   mu = mod$model$mu)


      names(Res) <- c("coefficients_1",
                      "coefficients_2",
                      "coefficients_3",
                      "coefficients_4",
                      "EBV_1",
                      "EBV_2",
                      "EBV_3",
                      "EBV_4",
                      "total_EBV",
                      "predict_value",
                      "genomic_variance_1",
                      "genomic_variance_2",
                      "genomic_variance_3",
                      "genomic_variance_4",
                      "total_genomic_variance",
                      "residual_error",
                      "genomic_heritability",
                      "M1_model_ready",
                      "M2_model_ready",
                      "M3_model_ready",
                      "g_model_ready",
                      "intercept")

    }

    ####

  }



  #}
  ### remove the generated output files from the working directory
  unlink(mod$output_files_names)
  return(Res)

}





