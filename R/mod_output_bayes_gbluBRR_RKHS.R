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
#' @param omics_kernel_label
#' @param pheno_data
#' @param heter_groups
#' @param bayes_para
#' @param ...
#' @param GS_model
#'
#' @return
#' @export
#'
#' @examples
mod_output_bayes_gbluBRR_RKHS <- function(mod=NULL,
                                         ETA=NULL,
                                         gen_name=NULL,
                                         GS_model = NULL,
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
                                         CI_width_thresholds = c(0.33, 0.66),
                                         confidence_level = 0.95,
                                         high_reliability_thres = 0.9,
                                         low_reliability_thres = 0.5,
                                         n_components = 20,
                                         threshold = 100,
                                         target = "test_set",
                                         iqr_multiplier = 1.5,
                                         interval_width_high_threshold = NULL,
                                         interval_width_low_threshold = NULL,
                                         interval_width_moderate_threshold = NULL,
                                         system_database = FALSE,
                                         ...){

  diagnostic_plots <- NULL
  tst <- which(is.na(mod$model$y))

  if(length(tst)>1){
  Standard_error <- mod$model$SD.yHat[tst]
  PEV <- (mod$model$SD.yHat[tst])^2
  Reliability <- 1 - (PEV/var(mod$model$yHat[tst]))

  } else {
    Standard_error <- mod$model$SD.yHat
    PEV <- (mod$model$SD.yHat)^2
    Reliability <- 1 - (PEV/var(mod$model$yHat))

  }

  msg <- "\n==================================================\n"
  ### Check if the user provide lable/name for the omics data

  if(inherits(omics_kernel_label,'list')){
    if(!all(sapply(omics_kernel_label, function(x){ is.null(x)}))!=FALSE){
      label <-  which(sapply(omics_kernel_label, function(x) !is.null(x)))
      print_lable <-  omics_kernel_label[label]
    } else {
      print_lable <-  NULL
    }

  }else {
    if(inherits(omics_kernel_label,"character")){
      print_lable <-  omics_kernel_label
    }
    if(is.null(omics_kernel_label)){
      print_lable <-  NULL
    }
  }

  datasets <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
  dataset_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
  datasets_index <- which(!sapply(datasets, is.null))
  datasets <-  datasets[datasets_index]
  dataset_names <- dataset_names[datasets_index]

  tst <- which(is.na(mod$model$y))
  tst_GID <- unique(as.character(pheno_data[tst, gen_name]))

  # Subset datasets if tst has more than 1 element
  data_trn <- NULL
  data_tst <- NULL

  if (length(tst) > 1) {
    if(!is.null(heter_groups)){
      data_trn <- lapply(datasets, function(dataset) dataset[!rownames(dataset)%in%tst_GID, !colnames(dataset)%in%tst_GID])
      data_tst <- lapply(datasets, function(dataset) dataset[rownames(dataset)%in%tst_GID, colnames(dataset)%in%tst_GID])
    }else{
    data_trn <- lapply(datasets, function(dataset) dataset[-tst, -tst])
    data_tst <- lapply(datasets, function(dataset) dataset[tst, tst])

    }
  }

  if(length(datasets)>1 & length(tst)>1){
    data_trn <- do.call(cbind, data_trn)
    data_trn <- scale(data_trn)
    #
    data_tst <- do.call(cbind, data_tst)
    data_tst <- scale(data_tst)

    dataset <- do.call(cbind, datasets)
  } else{
    dataset <- do.call(cbind, datasets)

  }

  ### Check bayes_parameter_check function in bayesians_preprocess for details
  nIter <- bayes_para[["nIter"]]
  burnIn <- bayes_para[["burnIn"]]

  posindex <- (burnIn + 1):nIter
  ##### For MET analysis
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){

    ### Residual value is only estimable for response value without NA
    #tst <- which(is.na(mod$model$y))
    if(length(tst)>1){
      ### incidence matrix for main eff. of the genotypes

      Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

      ### Extract all environments in MET
      all_envs_for_met <-  as.character(pheno_data[,heter_groups])

      ##
      result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                             predictions = mod$model$yHat[tst],
                                                             standard_errors = mod$model$SD.yHat[tst],
                                                             confidence_level = confidence_level,
                                                             model_for_CI_cal = "Bayes",
                                                             boot_results = NULL)

      result_rel <-  reliability_thresholds(prediction_error_var = (mod$model$SD.yHat[tst])^2,
                                            genetic_var = var(mod$model$yHat[tst]),
                                            high_reliability_thres = high_reliability_thres,
                                            low_reliability_thres = low_reliability_thres)

      composite_reliability <- composite_reliability_tst(geno_trn = data_trn,
                                                         geno_tst = data_tst,
                                                         geno_tst_trn = NULL,
                                                         names_tst = NULL,
                                                         names_trn = NULL,
                                                         n_components = n_components,
                                                         threshold = threshold,
                                                         target = target,
                                                         interval_width = result_rel_MPIW$Uncertainty,
                                                         CI_width_thresholds = CI_width_thresholds,
                                                         interval_width_high_threshold = interval_width_high_threshold,
                                                         interval_width_low_threshold = interval_width_low_threshold,
                                                         apply_pca = TRUE)

      # predicted_value <-  data.frame(name = pheno_data[tst,gen_name],
      #                                Env = all_envs_for_met,
      #                                Predicted_value = mod$model$yHat[tst],
      #                                stringsAsFactors = FALSE)
      #
      # names(predicted_value)[1] <-  gen_name
      predicted_value <- data.frame(name = pheno_data[tst,gen_name],
                                    Predicted_value = mod$model$yHat[tst],
                                    Env = all_envs_for_met[tst],
                                    Standard_error = mod$model$SD.yHat[tst],
                                    PEV = (mod$model$SD.yHat[tst])^2,
                                    lower_bound = result_rel_MPIW$lower_bound,
                                    upper_bound = result_rel_MPIW$upper_bound,
                                    Uncertainty = result_rel_MPIW$Uncertainty,
                                    Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
                                    Reliability = result_rel$reliability,
                                    Reliability_remarks = result_rel$remarks,
                                    Reliability_percentage = result_rel$reliability_percentage,
                                    #Composite_reliability = NA,
                                    #Composite_reliability_percentage = NA,
                                    stringsAsFactors = FALSE)

      colnames(predicted_value)[c(1, 3)] <- c(gen_name, heter_groups)


      # diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = NULL,
      #                                                     GID_names = pheno_data[tst,gen_name],
      #                                                     CI_width_thresholds = CI_width_thresholds,
      #                                                     predictions = mod$model$yHat[tst],
      #                                                     standard_errors = mod$model$SD.yHat[tst],
      #                                                     prediction_error_var = (mod$model$SD.yHat[tst])^2,
      #                                                     genetic_var = var(mod$model$yHat[tst]),
      #                                                     confidence_level = 0.95,
      #                                                     model_for_CI_cal = "Bayes",
      #                                                     composite_reliability_score = composite_reliability$reliability_score,
      #                                                     composite_reliability = composite_reliability$trustworthiness,
      #                                                     composite_reliability_percentage = composite_reliability$reliability_percentage,
      #                                                     #threshold = NULL,
      #                                                     high_reliability_thres = high_reliability_thres,
      #                                                     low_reliability_thres = low_reliability_thres,
      #                                                     system_database = system_database)


      residual_value <- data.frame(name = pheno_data[tst,gen_name],
                                   Env = all_envs_for_met[tst],
                                   Predicted_value = mod$model$yHat[tst],
                                   Residual_value = (mod$model$y[tst] - mod$model$yHat[tst]),
                                   stringsAsFactors = FALSE)

      colnames(residual_value)[c(1, 2)] <- c(gen_name, heter_groups)

    } else {

      ### incidence matrix for main eff. of the genotypes

      Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

      ### Extract all environments in MET
      all_envs_for_met <-  as.character(pheno_data[,heter_groups])

      # predicted_value = data.frame(name = pheno_data[,gen_name],
      #                              Env = all_envs_for_met,
      #                              Predicted_value = mod$model$yHat)
      # names(predicted_value)[1] <-  gen_name

      result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                             predictions = mod$model$yHat,
                                                             standard_errors = mod$model$SD.yHat,
                                                             confidence_level = confidence_level,
                                                             model_for_CI_cal = "Bayes",
                                                             boot_results = NULL)

      result_rel <-  reliability_thresholds(prediction_error_var = (mod$model$SD.yHat)^2,
                                            genetic_var = var(mod$model$yHat),
                                            high_reliability_thres = high_reliability_thres,
                                            low_reliability_thres = low_reliability_thres)


      # composite_reliability <- composite_reliability_tst(geno_trn = data_trn,
      #                                                    geno_tst = data_tst,
      #                                                    geno_tst_trn = NULL,
      #                                                    names_tst = NULL,
      #                                                    names_trn = NULL,
      #                                                    n_components = n_components,
      #                                                    threshold = threshold,
      #                                                    target = target,
      #                                                    iqr_multiplier = iqr_multiplier,
      #                                                    interval_width = result_rel_MPIW$Uncertainty,
      #                                                    interval_width_high_threshold = interval_width_high_threshold,
      #                                                    interval_width_moderate_threshold = interval_width_moderate_threshold,
      #                                                    apply_pca = TRUE)

      predicted_value <- data.frame(name = pheno_data[,gen_name],
                                    Predicted_value = mod$model$yHat,
                                    Env = all_envs_for_met,
                                    Standard_error = mod$model$SD.yHat,
                                    PEV = (mod$model$SD.yHat)^2,
                                    lower_bound = result_rel_MPIW$lower_bound,
                                    upper_bound = result_rel_MPIW$upper_bound,
                                    Uncertainty = result_rel_MPIW$Uncertainty,
                                    Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
                                    Reliability = result_rel$reliability,
                                    Reliability_remarks = result_rel$remarks,
                                    Reliability_percentage = result_rel$reliability_percentage,
                                    #Composite_reliability = NA,
                                    #Composite_reliability_percentage = NA,
                                    stringsAsFactors = FALSE)

      colnames(predicted_value)[c(1, 3)] <- c(gen_name, heter_groups)

      residual_value <- data.frame(name = pheno_data[,gen_name],
                                   Env = all_envs_for_met,
                                   Predicted_value = mod$model$yHat,
                                   Residual_value = (mod$model$y - mod$model$yHat),
                                   stringsAsFactors = FALSE)

      colnames(residual_value)[c(1, 2)] <- c(gen_name, heter_groups)
    }

    #### End MET
  } else{ ##  start single enviornment results
    ##################################
    ## Create output for predicted value and residual value.
    ## The residual value dataframe also contain predicted value for two reasons
    #1) For ease of plotting
    #2) When testing set is present in the real world it is expected to be

    ### Residual value is only estimable for response value without NA
    #tst <- which(is.na(mod$model$y))
    if(length(tst)>1){

      result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                             predictions = mod$model$yHat[tst],
                                                             standard_errors = mod$model$SD.yHat[tst],
                                                             confidence_level = confidence_level,
                                                             model_for_CI_cal = "Bayes",
                                                             boot_results = NULL)

      result_rel <-  reliability_thresholds(prediction_error_var = (mod$model$SD.yHat[tst])^2,
                                            genetic_var = var(mod$model$yHat[tst]),
                                            high_reliability_thres = high_reliability_thres,
                                            low_reliability_thres = low_reliability_thres)

      composite_reliability <- composite_reliability_tst(geno_trn = data_trn,
                                                         geno_tst = data_tst,
                                                         geno_tst_trn = NULL,
                                                         names_tst = NULL,
                                                         names_trn = NULL,
                                                         n_components = n_components,
                                                         threshold = threshold,
                                                         target = target,
                                                         interval_width = result_rel_MPIW$Uncertainty,
                                                         CI_width_thresholds = CI_width_thresholds,
                                                         interval_width_high_threshold = interval_width_high_threshold,
                                                         interval_width_low_threshold = interval_width_low_threshold,
                                                         apply_pca = TRUE)

      predicted_value <- data.frame(name = pheno_data[tst,gen_name],
                                    Predicted_value = mod$model$yHat[tst],
                                    #Env = all_envs_for_met[tst],
                                    Standard_error = mod$model$SD.yHat[tst],
                                    PEV = (mod$model$SD.yHat[tst])^2,
                                    lower_bound = result_rel_MPIW$lower_bound,
                                    upper_bound = result_rel_MPIW$upper_bound,
                                    Uncertainty = result_rel_MPIW$Uncertainty,
                                    Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
                                    Reliability = result_rel$reliability,
                                    Reliability_remarks = result_rel$remarks,
                                    Reliability_percentage = result_rel$reliability_percentage,
                                    #Composite_reliability = composite_reliability$trustworthiness,
                                    #Composite_reliability_percentage = composite_reliability$reliability_percentage,
                                    stringsAsFactors = FALSE)

      colnames(predicted_value)[1] <- c(gen_name)

      diagnostic_plots <- diagnostic_plot_true_prediction(boot_results = NULL,
                                                          GID_names = pheno_data[tst,gen_name],
                                                          CI_width_thresholds = CI_width_thresholds,
                                                          predictions = mod$model$yHat[tst],
                                                          standard_errors = mod$model$SD.yHat[tst],
                                                          prediction_error_var = (mod$model$SD.yHat[tst])^2,
                                                          genetic_var = var(mod$model$yHat[tst]),
                                                          confidence_level = 0.95,
                                                          model_for_CI_cal = "Bayes",
                                                          composite_reliability_score = composite_reliability$reliability_score,
                                                          composite_reliability = composite_reliability$trustworthiness,
                                                          composite_reliability_percentage = composite_reliability$reliability_percentage,
                                                          #threshold = NULL,
                                                          high_reliability_thres = high_reliability_thres,
                                                          low_reliability_thres = low_reliability_thres,
                                                          system_database = system_database)

      residual_value <- data.frame(name = pheno_data[tst,gen_name],
                                   #Env = all_envs_for_met[tst],
                                   Predicted_value = mod$model$yHat[tst],
                                   Residual_value = (mod$model$y[tst] - mod$model$yHat[tst]),
                                   stringsAsFactors = FALSE)

      colnames(residual_value)[1] <- c(gen_name)
    } else {
      result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(CI_width_thresholds = CI_width_thresholds,
                                                             predictions = mod$model$yHat,
                                                             standard_errors = mod$model$SD.yHat,
                                                             confidence_level = confidence_level,
                                                             model_for_CI_cal = "Bayes",
                                                             boot_results = NULL)

      result_rel <-  reliability_thresholds(prediction_error_var = (mod$model$SD.yHat)^2,
                                            genetic_var = var(mod$model$yHat),
                                            high_reliability_thres = high_reliability_thres,
                                            low_reliability_thres = low_reliability_thres)

      composite_reliability <- composite_reliability_tst(geno_trn = data_trn,
                                                         geno_tst = data_tst,
                                                         geno_tst_trn = NULL,
                                                         names_tst = NULL,
                                                         names_trn = NULL,
                                                         n_components = n_components,
                                                         threshold = threshold,
                                                         target = target,
                                                         interval_width = result_rel_MPIW$Uncertainty,
                                                         CI_width_thresholds = CI_width_thresholds,
                                                         interval_width_high_threshold = interval_width_high_threshold,
                                                         interval_width_low_threshold = interval_width_low_threshold,
                                                         apply_pca = TRUE)

      predicted_value <- data.frame(name = pheno_data[,gen_name],
                                    Predicted_value = mod$model$yHat,
                                    #Env = all_envs_for_met,
                                    Standard_error = mod$model$SD.yHat,
                                    PEV = (mod$model$SD.yHat)^2,
                                    lower_bound = result_rel_MPIW$lower_bound,
                                    upper_bound = result_rel_MPIW$upper_bound,
                                    Uncertainty = result_rel_MPIW$Uncertainty,
                                    Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
                                    Reliability = result_rel$reliability,
                                    Reliability_remarks = result_rel$remarks,
                                    Reliability_percentage = result_rel$reliability_percentage,
                                    #Composite_reliability = composite_reliability$trustworthiness,
                                    #Composite_reliability_percentage = composite_reliability$reliability_percentage,
                                    stringsAsFactors = FALSE)

      colnames(predicted_value)[1] <- c(gen_name)

      residual_value <- data.frame(name = pheno_data[,gen_name],
                                   Predicted_value = mod$model$yHat,
                                   Residual_value = (mod$model$y - mod$model$yHat),
                                   stringsAsFactors = FALSE)

      colnames(residual_value)[1] <- c(gen_name)

    }

    Zg <- NULL
    all_envs_for_met <- NULL
  } ##  End single enviornment results



  #######################################
  if(GS_model=="RKHS"){
  BIN <- mod[["output_files_names"]][grepl("_varU.dat", mod[["output_files_names"]])]
  } else {
    if(GS_model=="BRR"){
      BIN <- mod[["output_files_names"]][grepl("bin", mod[["output_files_names"]])]
    }
  }
  ### Extract Error variance
  var_residual <- scan(mod[["output_files_names"]][grepl("varE.dat", mod[["output_files_names"]])],
                       what = numeric(),
                       sep = "\n", quiet =TRUE)
  #var_residual <- var_residual[posindex]
  var_residual <- var_residual

  se_var_residual <- standard_deviation(var_residual)
  #########
  if (GS_model == "BRR") {
    varB_files <- mod[["output_files_names"]][grepl("varB.dat", mod[["output_files_names"]])]
  } else {
    if (GS_model == "RKHS") {
    varB_files <- mod$output_files_names[grepl("_varU.dat", mod[["output_files_names"]])]

    }
  }

  #########

  #Var_U_omics_list <- vector(mode = "list", length = length(BIN))
  var_u_mean_omics_list <- list()
  se_var_u_omics_list <-  list()
  var_u_total <- 0
  genomic_h2 <- 0
  se_genomic_h2 <- 0
  datasets <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
  dataset_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
  datasets_index <- which(!sapply(datasets, is.null))
  datasets <-  datasets[datasets_index]
  dataset_names <- dataset_names[datasets_index]
  posterior_list <-   list()
  res_coeff_ebv_pev_rel_se_list <-  list()
  coefficients_list <-   list()
  estimated_breeding_value_list <-  list()
  m_matrix_model_ready_list <-   list()
  sum_posterior <-  0
  sum_estimated_breeding_value <-  0
  extracted_names_from_ETA_list <-  ETA[["ETA_element_name"]]
  if(length(dataset_names)==length(extracted_names_from_ETA_list)) {
    # Reorder dataset_names based on the order of name in extracted_names_from_ETA_list
    dataset_names <- dataset_names[match(extracted_names_from_ETA_list, dataset_names)]

  } else {
    stop("names must be the same length")
  }

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    if(!is.null(Zg)){
    ZgZg <- Zg%*%dataset%*%t(Zg)

    suppressMessages({
    ZgZg <- grm_kernel_precheck(ZgZg)
    })

    }

    if (!is.null(dataset)) {
      ## the first and second will cannot togther because I am looking through the datasets can might contain NULL
      #if (length(ETA[["ETA_element_name"]]) <= i){
      if(ETA[["ETA_element_name"]][i]==dataset_names[i]){
        #Var_U_omics_list[[dataset_names[i]]] <- process_varU(varB_files[i], posindex, GS_model)
        #gid_name <- rownames(geno_data)
        var_u_omics <- process_var_u(varB_files[i], posindex, GS_model)
        #posterior_list[[dataset_names[i]]] <- BGLR::readBinMat(BIN[aa])
        if(GS_model=="BRR"){
        res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]] <- cal_coeff_ebv_pev_rel_se_new(
          beta = BGLR::readBinMat(BIN[i]),
          x_variable = if(!is.null(Zg)) ZgZg else dataset,
          gen_name = gen_name,
          hetero = all_envs_for_met,
          heter_groups = heter_groups,
          var_u = mean(var_u_omics),
          mod =  mod,
          gid_name =rownames(dataset))

        } else {
          if(GS_model=="RKHS"){
        res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]] <-   cal_coeff_ebv_pev_rel_RHKS_glub(
                                                        mod =  mod,
                                                        gmatrix = if(!is.null(Zg)) ZgZg else dataset,
                                                        gen_name = gen_name,
                                                        var_u = mean(var_u_omics),
                                                        var_E = mean(var_residual),
                                                        hetero = all_envs_for_met,
                                                        heter_groups = heter_groups,
                                                        gid_name = rownames(dataset),
                                                        )


          }
        }
        var_u_mean_omics_list[[dataset_names[i]]] <- mean(var_u_omics)
        se_var_u_omics_list[[dataset_names[i]]] <- standard_deviation(var_u_omics)
        coefficients_list[[paste("coefficient",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Coefficient"]]
        estimated_breeding_value_list[[paste("estimated_breeding_value",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]]
        if(GS_model=="BRR"){
        sum_posterior <- sum_posterior + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Posterior"]]
        }
        sum_estimated_breeding_value <- sum_estimated_breeding_value + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]][,"Estimated_breeding_value"]
        if(dataset_names[i]=="gmatrix"){
        m_matrix_model_ready_list[[paste(gsub("gmatrix", "geno", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
        } else{
          m_matrix_model_ready_list[[paste(gsub("_kernel", "", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
        }
        var_u_total <- var_u_total + var_u_omics

        if(length(datasets)==1){
          # predicted_value[, 1] <- rownames(dataset)
          # #PEV <- apply(g_ebv, 1, var)
          # predicted_value <- predicted_value |>
          #   dplyr::mutate(
          #                 #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
          #                 Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
          #                 Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
          #                 Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])

          if(is.null(Zg)){
          sum_ebv <- data.frame(name = rownames(dataset),
                                Estimated_breeding_value = sum_estimated_breeding_value,
                                stringsAsFactors = FALSE)

          colnames(sum_ebv)[1] <- gen_name

          } else{
            sum_ebv <- data.frame(name =rownames(dataset),
                                  Env = all_envs_for_met,
                                  Estimated_breeding_value = sum_estimated_breeding_value,
                                  stringsAsFactors = FALSE)

            colnames(sum_ebv)[1:2] <- c(gen_name, heter_groups)
          }
          sum_ebv <- sum_ebv |>
            dplyr::mutate(
                          #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
                          Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
                          Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
                          Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])


          if(length(tst)>1){
            #residual_value[, 1] <- rownames(dataset)[tst]
            sum_ebv <- sum_ebv[tst, ]


          }
          # else {
          #
          #   residual_value[, 1] <- rownames(dataset)
          # }
          #### MET
          if(!is.null(Zg)){
            sep_pev_rel <- sep_pev_rel_gblup(geno_object = dataset,
                                             va = mean(var_u_omics),
                                             ve = mean(var_residual))

            across_env_predicted_value <-  as.data.frame(predicted_value[tst, ] |>
                                                           dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                           dplyr::summarise(
                                                             Predicted_value = mean(Predicted_value, na.rm = TRUE),
                                                             Standard_error = mean(Standard_error, na.rm = TRUE),
                                                             PEV = mean(PEV, na.rm = TRUE),
                                                             Reliability = mean(Reliability, na.rm = TRUE)
                                                           )
                                                           # dplyr::summarise(Predicted_value = mean(Predicted_value)) |>
                                                           # dplyr::mutate(
                                                           #               Standard_error = sep_pev_rel[,"sep"],
                                                           #               #Standard_error = Standard_error,
                                                           #               Prediction_error_variance = sep_pev_rel[,"pev"],
                                                           #               Reliability = sep_pev_rel[,"rel"])
            )

          }
          ### MET End
        } else{

          if(length(datasets)>1){
            if(i==1) gid_name <- rownames(dataset)
            #predicted_value[, 1] <- gid_name
            #residual_value[, 1] <- gid_name
            ##### Treat sum_EBV
            if(i==length(datasets)){
            if(GS_model=="BRR"){
            pev <- apply(sum_posterior, 1, var)
            rel <- 1 - (pev / mean(var_u_total))
            #rel <- ifelse(rel<0, NA, rel)

            } else {
              pev <- NA
              rel <- NA


            }

            if(is.null(Zg)){
              sum_ebv <- data.frame(name = gid_name,
                                    Estimated_breeding_value = sum_estimated_breeding_value,
                                    stringsAsFactors = FALSE)

              colnames(sum_ebv)[1] <- gen_name

            } else{
              sum_ebv <- data.frame(name = gid_name,
                                    Env = all_envs_for_met,
                                    Estimated_breeding_value = sum_estimated_breeding_value,
                                    stringsAsFactors = FALSE)

              colnames(sum_ebv)[1:2] <- c(gen_name, heter_groups)
            }

              if(length(tst)>1) sum_ebv <- sum_ebv[tst, ]
            sum_ebv <- sum_ebv |>
              dplyr::mutate(
                           #Standard_error = ifelse(!is.na(pev), sqrt(pev), NA),
                            Standard_error = Standard_error,
                            Prediction_error_variance = PEV,
                            Reliability = Reliability)

            # predicted_value <- predicted_value |>
            #   dplyr::mutate(
            #                 #Standard_error = ifelse(!is.na(pev), sqrt(pev), NA),
            #                 Standard_error = Standard_error,
            #                 Prediction_error_variance = PEV,
            #                 Reliability = Reliability)

            ## MET
            if(!is.null(Zg) & length(datasets)>1){


                across_env_predicted_value <-  as.data.frame(predicted_value[tst, ] |>
                                                               dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>  # Use across() to refer to the column specified by gen_name
                                                               dplyr::summarise(
                                                                 Predicted_value = mean(Predicted_value, na.rm = TRUE),
                                                                 Standard_error = mean(Standard_error, na.rm = TRUE),
                                                                 PEV = mean(PEV, na.rm = TRUE),
                                                                 Reliability = mean(Reliability, na.rm = TRUE)
                                                               )

                                                               # dplyr::mutate(
                                                               #
                                                               #   Standard_error = sep_pev_rel[,"sep"],
                                                               #   #Standard_error = Standard_error,
                                                               #   Prediction_error_variance = sep_pev_rel[,"pev"],
                                                               #   Reliability = sep_pev_rel[,"rel"])
                )


            }        ## End MET

            }
          } ##
        }

      }

      #}
    } ## End

  }

  variance_components <- bayes_variance_componentsnew(var_u_mean_omics_list,
                                                      se_var_u_omics_list,
                                                      var_u_total,
                                                      var_residual,
                                                      se_var_residual)

  if(!"geno_model_ready" %in%names(m_matrix_model_ready_list)){
    if(length(m_matrix_model_ready_list)>1){
      names(m_matrix_model_ready_list) <- paste(paste("Omics", seq_along(m_matrix_model_ready_list), sep = ""), "model_ready", sep = "_")
      names(coefficients_list) <- paste(paste("Omics", seq_along(coefficients_list), sep = ""), "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- paste(paste("Omics", seq_along(estimated_breeding_value_list), sep = ""), "estimated_breeding_value", sep = "_")
    } else {
      names(m_matrix_model_ready_list) <- paste("Omics", "model_ready", sep = "_")
      names(coefficients_list) <- "Coefficient" # paste("Omics", "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- "Estimated_breeding_value"  # paste("Omics", "estimated_breeding_value", sep = "_")

    }
  } else if(length(grep("omic", names(m_matrix_model_ready_list)))==0 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
    geno_index <- grep("geno", names(m_matrix_model_ready_list))
    names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
    names(coefficients_list)[geno_index] <- "Coefficient" # paste("Geno", "coefficient", sep = "_")
    names(estimated_breeding_value_list)[geno_index] <- "Estimated_breeding_value" # paste("Geno", "estimated_breeding_value", sep = "_")

  } else{
    if(length(grep("omic", names(m_matrix_model_ready_list)))>=1 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
      omic_index <- grep("omic", names(m_matrix_model_ready_list))
      geno_index <- grep("geno", names(m_matrix_model_ready_list))
      if(length(omic_index)>1){
        names(m_matrix_model_ready_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "estimated_breeding_value", sep = "_")
      } else {
        names(m_matrix_model_ready_list)[omic_index] <- paste("Omics", "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste("Omics", "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste("Omics", "estimated_breeding_value", sep = "_")

      }
      ####
      names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
      names(coefficients_list)[geno_index] <- paste("Geno", "coefficient", sep = "_")
      names(estimated_breeding_value_list)[geno_index] <- paste("Geno", "estimated_breeding_value", sep = "_")

    }

  }


  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)

    if(!is.null(geno_data)){
      print_lable <- c("Genomic", print_lable)
    }
    if(length(print_lable)>length(m_matrix_model_ready_list) | length(print_lable)<length(m_matrix_model_ready_list)){
      message(insight::print_color(paste(msg,paste("More than two Omics lable were provided. Default name was applied.")), "blue"))



    } else {
      if(length(print_lable)==length(m_matrix_model_ready_list)){
        rownames(variance_components)[1:length(print_lable)] <- paste(print_lable, "variance", sep="_")
        ######
        names(coefficients_list) <-  paste(print_lable, "coefficient", sep="_")
        names(estimated_breeding_value_list) <- paste(print_lable, "estimated_breeding_value", sep="_")
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(m_matrix_model_ready_list) <- paste(print_lable, "model_ready", sep="_")
        ####
      }
    }

  } else {

    message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.")), "blue"))


  }


  if(is.null(Zg)){

    #if(!is.null(diagnostic_plots)){
    res <- list(Coefficients = coefficients_list,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Residual_value = residual_value,
                diagnostic_plots = diagnostic_plots,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list)
    # } else {
    #   res <- list(Coefficients = coefficients_list,
    #               Estimated_breeding_value = estimated_breeding_value_list,
    #               Total_estimated_breeding_value = sum_ebv,
    #               Predicted_value =  predicted_value,
    #               Residual_value = residual_value,
    #               #diagnostic_plots = diagnostic_plots,
    #               Variance_components = variance_components,
    #               M_matrix_model_ready =  m_matrix_model_ready_list)
    # }

  } else {
    res <- list(Coefficients = coefficients_list,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Total_Predicted_value = across_env_predicted_value,
                Residual_value = residual_value,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list,
                diagnostic_plots = diagnostic_plots
    )

  }
  ### remove the generated output files from the working directory
  unlink(mod[["output_files_names"]])

  return(res)

}
