
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
standard_deviation <- function(x ){
  stats::sd(x, na.rm = TRUE)
}

bayes_marker_compute_target_summary <- function(mod,
                                                 ETA,
                                                 output_files_names,
                                                 component_variance_means = NULL,
                                                 draw_indices = NULL,
                                                 high_reliability_thres = 0.9,
                                                 low_reliability_thres = 0.5) {
  eta_input <- ETA[["ETA"]]
  eta_types <- vapply(eta_input, function(x) x$model, character(1))
  random_idx <- which(eta_types != "FIXED")

  total_draws <- NULL
  genetic_mean <- rep(0, length(mod$model$yHat))
  prior_genetic_variance <- rep(0, length(mod$model$yHat))
  train_idx <- which(!is.na(mod$model$y))
  for (eta_idx in random_idx) {
    b_file <- output_files_names[grepl(paste0("ETA_", eta_idx, "_b\\.bin$"), basename(output_files_names))]
    if (length(b_file) == 0) {
      stop("Unable to locate Bayesian posterior draw file for ETA term.", call. = FALSE)
    }
    beta_draws <- as.matrix(BGLR::readBinMat(b_file[1]))
    if (length(draw_indices) && max(draw_indices) <= nrow(beta_draws)) {
      beta_draws <- beta_draws[draw_indices, , drop = FALSE]
    }
    X_eta <- as.matrix(eta_input[[eta_idx]]$X)
    eta_draws <- as.matrix(X_eta %*% t(beta_draws))
    total_draws <- if (is.null(total_draws)) eta_draws else total_draws + eta_draws
    genetic_mean <- genetic_mean + drop(X_eta %*% mod$model$ETA[[eta_idx]]$b)
  }

  if (!is.null(component_variance_means)) {
    comp_means <- as.numeric(component_variance_means)
    if (length(comp_means) >= length(random_idx)) {
      for (k in seq_along(random_idx)) {
        X_eta <- as.matrix(eta_input[[random_idx[k]]]$X)
        row_energy <- rowSums(X_eta ^ 2)
        reference_energy <- mean(row_energy[train_idx], na.rm = TRUE)
        if (!is.finite(reference_energy) || reference_energy <= 0) {
          reference_energy <- mean(row_energy, na.rm = TRUE)
        }
        k_diag <- if (is.finite(reference_energy) && reference_energy > 0) {
          row_energy / reference_energy
        } else {
          rep(1, nrow(X_eta))
        }
        prior_genetic_variance <- prior_genetic_variance + (comp_means[k] * k_diag)
      }
    }
  }

  mu_draws <- gp_bayes_read_mu_draws(
    output_files_names,
    draw_indices = draw_indices,
    n_draws = ncol(total_draws)
  )
  if (is.null(mu_draws)) {
    mu_draws <- rep(as.double(mod$model$mu), ncol(total_draws))
  }
  pred_draws <- as.matrix(sweep(total_draws, 2, mu_draws, FUN = "+"))
  pred_mean <- as.double(mod$model$yHat)
  standard_error <- apply(pred_draws, 1, stats::sd, na.rm = TRUE)
  prediction_error_var <- apply(pred_draws, 1, stats::var, na.rm = TRUE)
  fallback_genetic_variance <- mean(
    apply(total_draws[train_idx, , drop = FALSE], 2, stats::var, na.rm = TRUE),
    na.rm = TRUE
  )
  if (!is.finite(fallback_genetic_variance)) {
    fallback_genetic_variance <- 0
  }
  if (all(!is.finite(prior_genetic_variance)) || all(prior_genetic_variance <= 0, na.rm = TRUE)) {
    target_genetic_variance <- rep(fallback_genetic_variance, length(prediction_error_var))
  } else {
    target_genetic_variance <- prior_genetic_variance
    bad_idx <- !is.finite(target_genetic_variance) | target_genetic_variance <= 0
    target_genetic_variance[bad_idx] <- fallback_genetic_variance
  }

  reliability_res <- reliability_thresholds(
    prediction_error_var = prediction_error_var,
    genetic_var = target_genetic_variance,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres
  )

  list(
    predicted_mean = pred_mean,
    genetic_mean = genetic_mean,
    standard_error = standard_error,
    prediction_error_var = prediction_error_var,
    target_genetic_variance = target_genetic_variance,
    reliability = reliability_res$reliability,
    reliability_remarks = reliability_res$remarks,
    reliability_percentage = reliability_res$reliability_percentage,
    pred_draws = pred_draws
  )
}

#' Title
#'
#' @param file
#' @param posindex
#' @param GS_model
#'
#' @return
#' @export
#'
#' @examples
process_var_u <- function(file, posindex, GS_model) {
  if(GS_model%in%c("BRR", "RKHS")){
    #var_u <- scan(file, what = numeric(), sep = "\n", quiet = TRUE)[posindex]
    var_u <- scan(file, what = numeric(), sep = "\n", quiet = TRUE)

  }else {
    if(GS_model%in%c("BayesA", "BayesC", "BL")){ # lambda
      var_u <- utils::read.table(file)
    } else {
      if(GS_model=="BayesB"){
        var_u <- utils::read.table(file, skip = 1)  # Bayes B

      }

    }
    #var_u <- as.data.frame(tidyr::separate_rows(var_u))[posindex, 1]
    #var_u <- as.data.frame(tidyr::separate_rows(var_u))[, 1]
    var_u <- var_u[, 1]

  }
  if (length(posindex) && all(is.finite(posindex)) &&
      min(posindex) >= 1L && max(posindex) <= length(var_u)) {
    var_u <- var_u[posindex]
  }
  #Var_U_Se_omics <- standard_deviation(Var_U)

  # return(list(Var_U_omics = Var_U,
  #             Var_U_omics_mean  = mean(Var_U),
  #             Var_U_Se_omics = Var_U_Se_omics))

  return(var_u)
}

mod_output_bayes <- function(mod=NULL,
                             ETA=NULL,
                              pheno_data = NULL,
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
                              response_family = "gaussian",
                              response = NULL,
                              ...){

fam <- gp_resolve_response_family(response_family, y = mod$model$y)
if (!identical(fam, "gaussian")) {
  res <- gp_bayes_classification_output(
    mod = mod,
    pheno_data = pheno_data %||% ETA[["pheno_data"]],
    gen_name = gen_name,
    response = response,
    response_family = fam,
    high_reliability_thres = high_reliability_thres,
    low_reliability_thres = low_reliability_thres,
    ETA = ETA,
    bayes_para = bayes_para,
    GS_model = GS_model,
    confidence_level = confidence_level
  )
  unlink(mod[["output_files_names"]])
  return(res)
}

  #browser()
  # Standard_error = mod$model$SD.yHat
  # PEV <- (mod$model$SD.yHat)^2
  # Reliability <- 1 - (PEV/var(mod$model$yHat))


msg <- "\n ==================================================\n"
### Check if the user provide lable/name for the omics data

diagnostic_plots <- NULL
tst <- NULL

if(inherits(omics_data_label, 'list')){
  if(!all(sapply(omics_data_label, function(x){ is.null(x)}))!=FALSE){

    label <-  which(sapply(omics_data_label, function(x) !is.null(x)))

    print_lable <- omics_data_label[label]

  } else {
    print_lable <-  NULL
  }

}else {
  if(inherits(omics_data_label, "character")){

    print_lable <-  omics_data_label

  }

  if(is.null(omics_data_label)){
    print_lable <-  NULL

  }
}

# if(GS_model == "BL") {
#   GS_model <- "lambda"
# }

datasets <- list(geno_data, omic1_data, omic2_data, omic3_data)
dataset_names <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")
datasets_index <- which(!sapply(datasets, is.null))
datasets <-  datasets[datasets_index]
dataset_names <- dataset_names[datasets_index]

if(length(datasets)>1){
  datasets <- do.call(cbind, datasets)
  datasets <- scale(datasets)
} else{

  datasets <- do.call(cbind, datasets)
}

train_test_label <-  ifelse(is.na(mod$model$y), "Test", "Train")
tst <- which(is.na(mod$model$y))
if(length(tst)==0){
  tst <- NULL
}
### Check bayes_parameter_check function in bayesians_preprocess for details
nIter <- bayes_para[["nIter"]]
burnIn <- bayes_para[["burnIn"]]

posindex <- gp_bayes_saved_draw_indices(nIter = nIter, burnIn = burnIn, thin = bayes_para[["thin"]])
##################################
## Create output for predicted value and residual value.
## The residual value dataframe also contain predicted value for two reasons
#1) For ease of plotting
#2) When testing set is present in the real world it is expected to be


### Residual value is only estimable for response value without NA

if(length(tst)!=0){
  residual_value <- data.frame(name = rownames(datasets)[tst],
                               Predicted_value = mod$model$yHat[tst],
                               Residual_value = (mod$model$y[tst] - mod$model$yHat[tst]),
                               stringsAsFactors = FALSE)
} else {
  residual_value <- data.frame(name =rownames(datasets),
                               Predicted_value = mod$model$yHat,
                               Residual_value = (mod$model$y - mod$model$yHat),
                               stringsAsFactors = FALSE)
}

colnames(residual_value)[1] <- gen_name
#######################################
BIN <- mod[["output_files_names"]][grepl("bin", mod[["output_files_names"]])]
### Extract Error variance
var_residual <- gp_bayes_read_varE_draws(mod[["output_files_names"]], draw_indices = posindex)
#var_residual <- var_residual[posindex]
var_residual <- var_residual
se_var_residual <- standard_deviation(var_residual)
#########
  if (GS_model == "BRR") {
    varB_files <- mod[["output_files_names"]][grepl("varB.dat", mod[["output_files_names"]])]
  } else {
    varB_files <- mod[["output_files_names"]][grepl(paste(GS_model, "dat", sep = "."), mod[["output_files_names"]])]
  }

#########

#Var_U_omics_list <- vector(mode = "list", length = length(BIN))
var_u_mean_omics_list <- list()
se_var_u_omics_list <-  list()
var_u_total <- 0
genomic_h2 <- 0
se_genomic_h2 <- 0
Predicted_value_for_CI_total <- 0
datasets <- list(geno_data, omic1_data, omic2_data, omic3_data)
dataset_names <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")
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

datasets_cbind <- as.matrix(do.call(cbind, datasets))
eta_models <- vapply(ETA[["ETA"]], function(x) x$model, character(1))
random_eta_idx <- which(eta_models != "FIXED")
if (length(random_eta_idx) < length(datasets)) {
  stop(paste(msg, "Unable to align Bayesian ETA terms with the supplied data matrices."), call. = FALSE)
}
dataset_eta_idx <- random_eta_idx[seq_along(datasets)]

extracted_names_from_ETA_list <-  ETA[["ETA_element_name"]]
if(length(dataset_names)==length(extracted_names_from_ETA_list)) {
  # Reorder dataset_names based on the order of name in extracted_names_from_ETA_list
  dataset_names <- dataset_names[match(extracted_names_from_ETA_list, dataset_names)]

} else {
  stop(paste(msg, "names must be the same length"), call. = FALSE)
}


for (i in seq_along(datasets)) {
  dataset <- as.matrix(datasets[[i]])
  eta_idx <- dataset_eta_idx[i]
  model_matrix <- as.matrix(ETA[["ETA"]][[eta_idx]]$X)
  if (!is.null(dataset)) {
    ## the first and second will cannot togther because I am looking through the datasets can might contain NULL
    #if (length(ETA[["ETA_element_name"]]) <= i){
    if(ETA[["ETA_element_name"]][i]==dataset_names[i]){
    #Var_U_omics_list[[dataset_names[i]]] <- process_varU(varB_files[i], posindex, GS_model)
    #gid_name <- rownames(geno_data)
    #var_u_omics <- process_var_u(varB_files[i], posindex, GS_model)
      if(GS_model%in%c("BayesA", "BayesC", "BL", "BRR", "BayesB")){
    beta_draws <- as.matrix(BGLR::readBinMat(BIN[i]))
    if (length(posindex) && max(posindex) <= nrow(beta_draws)) {
      beta_draws <- beta_draws[posindex, , drop = FALSE]
    }
    var_u_and_others_omics <- process_var_u_new(GS_model = GS_model,
                                     geno_data = if (GS_model == "BRR") model_matrix else dataset,
                                     y = mod$model$y,
                                     B = beta_draws,
                                     tst = tst)
    var_u_omics <-  var_u_and_others_omics[["var_u_omics"]]
    Predicted_value_for_CI <- compute_predicted_value(geno_data = if (GS_model == "BRR") model_matrix else dataset,
                                                      bMat = beta_draws,
                                                      mu_values = mod$model$mu)

      } else{
        if(GS_model%in%c("RKHS")){
        var_u_omics <- process_var_u(varB_files[i], posindex, GS_model)
        }
      }
    #posterior_list[[dataset_names[i]]] <- BGLR::readBinMat(BIN[aa])
    res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]] <- cal_coeff_ebv_pev_rel_se_new(
                                                              beta = beta_draws,
                                                              x_variable = if (GS_model == "BRR") model_matrix else dataset,
                                                              gen_name = gen_name,
                                                              var_u = mean(var_u_omics),
                                                              gid_name =rownames(dataset),
                                                              mod =  mod)

    var_u_mean_omics_list[[dataset_names[i]]] <- mean(var_u_omics)
    se_var_u_omics_list[[dataset_names[i]]] <- standard_deviation(var_u_omics)
    coefficients_list[[paste("coefficient",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Coefficient"]]
    estimated_breeding_value_list[[paste("estimated_breeding_value",dataset_names[i], sep = "_")]] <- res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]]
    sum_posterior <- sum_posterior + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Posterior"]]
    sum_estimated_breeding_value <- sum_estimated_breeding_value + res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Estimated_breeding_value"]][,"Estimated_breeding_value"]
    m_matrix_model_ready_list[[paste(gsub("_data", "", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
    var_u_total <- var_u_total + var_u_omics
    Predicted_value_for_CI_total <- Predicted_value_for_CI_total  + Predicted_value_for_CI

    if(length(datasets)==1){
      #predicted_value[, 1] <- rownames(dataset)
      #PEV <- apply(g_ebv, 1, var)
      # predicted_value <- predicted_value |>
      #   dplyr::mutate(
      #                 #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
      #                 Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
      #                 Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
      #                 Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])

      sum_ebv <- data.frame(name = rownames(dataset),
                            Estimated_breeding_value = sum_estimated_breeding_value,
                            stringsAsFactors = FALSE)

      colnames(sum_ebv)[1] <- gen_name
      sum_ebv <- sum_ebv |>
        dplyr::mutate(
                      #Standard_error = sqrt(res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]]),
                      Standard_error = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Standard_error"]],
                      Prediction_error_variance = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["PEV"]],
                      Reliability = res_coeff_ebv_pev_rel_se_list[[dataset_names[i]]][["Reliability"]])

      if(length(tst)>0) sum_ebv <- sum_ebv[tst, ]
      if(length(tst)!=0){
        residual_value[, 1] <- rownames(dataset)[tst]

      }else {

        residual_value[, 1] <- rownames(dataset)
      }

    } else{

      if(length(datasets)>1){
        if(i==1) gid_name <- rownames(dataset)
        #predicted_value[, 1] <- gid_name
        if(length(tst)!=0){
        residual_value[, 1] <- gid_name[tst]
        } else {
          residual_value[, 1] <- gid_name
        }
        ##### Treat sum_EBV
        if(i==length(datasets)){
        sum_ebv <- data.frame(name = gid_name,
                              Estimated_breeding_value = sum_estimated_breeding_value,
                              stringsAsFactors = FALSE)
        colnames(sum_ebv)[1] <- gen_name
        if(length(tst)>0) sum_ebv <- sum_ebv[tst, ]
        }
      }
    }

    }

#
    #}
  } ## End

}

if (is.matrix(sum_posterior) && ncol(sum_posterior) > 1L) {
  variance_rows <- if (length(tst)) setdiff(seq_len(nrow(sum_posterior)), tst) else seq_len(nrow(sum_posterior))
  var_u_total <- apply(sum_posterior[variance_rows, , drop = FALSE], 2, stats::var, na.rm = TRUE)
}

variance_components <- bayes_variance_componentsnew(var_u_mean_omics_list,
                                                   se_var_u_omics_list,
                                                   var_u_total,
                                                   var_residual,
                                                   se_var_residual)

target_summary <- bayes_marker_compute_target_summary(
  mod = mod,
  ETA = ETA,
  output_files_names = mod[["output_files_names"]],
  component_variance_means = unlist(var_u_mean_omics_list, use.names = FALSE),
  draw_indices = posindex,
  high_reliability_thres = high_reliability_thres,
  low_reliability_thres = low_reliability_thres
)

if(!"geno_model_ready" %in%names(m_matrix_model_ready_list)){
  if(length(m_matrix_model_ready_list)>1){
  names(m_matrix_model_ready_list) <- paste(paste("Omics", seq_along(m_matrix_model_ready_list), sep = ""), "model_ready", sep = "_")
  names(coefficients_list) <- paste(paste("Omics", seq_along(coefficients_list), sep = ""), "coefficient", sep = "_")
  names(estimated_breeding_value_list) <- paste(paste("Omics", seq_along(estimated_breeding_value_list), sep = ""), "estimated_breeding_value", sep = "_")
  } else {
    names(m_matrix_model_ready_list) <- paste("Omics", "model_ready", sep = "_")
    names(coefficients_list) <- "Coefficient"  # paste("Omics", "coefficient", sep = "_")
    names(estimated_breeding_value_list) <- "Estimated_breeding_value"  # paste("Omics", "estimated_breeding_value", sep = "_")

}
  } else if(length(grep("omic", names(m_matrix_model_ready_list)))==0 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
  geno_index <- grep("geno", names(m_matrix_model_ready_list))
  names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
  names(coefficients_list)[geno_index] <- "Coefficient" # paste("Geno", "coefficient", sep = "_")
  names(estimated_breeding_value_list)[geno_index] <- "Estimated_breeding_value" #  paste("Geno", "estimated_breeding_value", sep = "_")

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
###############################
### New improved
result_rel_MPIW <- reliability_thresholds_MPIW_from_CI(
  CI_width_thresholds = CI_width_thresholds,
  predictions = target_summary$predicted_mean,
  Predicted_value_for_CI = target_summary$pred_draws - mod$model$mu,
  mod = mod,
  standard_errors = target_summary$standard_error,
  confidence_level = confidence_level,
  model_for_CI_cal = "Bayes",
  boot_results = NULL
)

if (length(tst) > 1) {
  composite_reliability <- composite_reliability_tst(
    geno_trn = datasets_cbind[-tst, , drop = FALSE],
    geno_tst = datasets_cbind[tst, , drop = FALSE],
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
    apply_pca = TRUE
  )
} else {
  composite_reliability <- composite_reliability_tst(
    geno_trn = datasets_cbind,
    geno_tst = NULL,
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
    apply_pca = TRUE
  )
}

predicted_value <- data.frame(
  name = rownames(dataset),
  Predicted_value = target_summary$predicted_mean,
  Train_Test_Label = train_test_label,
  Standard_error = target_summary$standard_error,
  PEV = target_summary$prediction_error_var,
  lower_bound = result_rel_MPIW$lower_bound,
  upper_bound = result_rel_MPIW$upper_bound,
  Uncertainty = result_rel_MPIW$Uncertainty,
  Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
  Reliability = target_summary$reliability,
  Reliability_remarks = target_summary$reliability_remarks,
  Reliability_percentage = target_summary$reliability_percentage,
  Genetic_variance = target_summary$target_genetic_variance,
  stringsAsFactors = FALSE
)

# AI_preds <- data.frame(name = GID,
#                        Predicted_value = AI_preds,
#                        Standard_error = pred_SE,
#                        PEV = pred_variances,
#                        lower_bound = result_rel_MPIW$lower_bound,
#                        upper_bound = result_rel_MPIW$upper_bound,
#                        Uncertainty = result_rel_MPIW$Uncertainty,
#                        Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#                        Reliability = result_rel$reliability,
#                        Reliability_remarks = result_rel$remarks,
#                        stringsAsFactors = FALSE)


# predicted_value <- data.frame(name = NA,
#                               Predicted_value = mod$model$yHat,
#                               Standard_error = mod$model$SD.yHat,
#                               PEV = (mod$model$SD.yHat)^2,
#                               lower_bound = result_rel_MPIW$lower_bound,
#                               upper_bound = result_rel_MPIW$upper_bound,
#                               Uncertainty = result_rel_MPIW$Uncertainty,
#                               Uncertainty_remarks = result_rel_MPIW$reliability_remarks,
#                               Reliability = result_rel$reliability,
#                               Reliability_remarks = result_rel$remarks,
#                               Reliability_percentage = result_rel$reliability_percentage,
#                               Composite_reliability = composite_reliability$trustworthiness,
#                               Composite_reliability_percentage = composite_reliability$reliability_percentage,
#                               stringsAsFactors = FALSE)

colnames(predicted_value)[colnames(predicted_value)%in%c("name")] <- c(gen_name)

predicted_value$Predicted_value <- target_summary$predicted_mean
predicted_value$Standard_error <- target_summary$standard_error
predicted_value$PEV <- target_summary$prediction_error_var
predicted_value$Reliability <- target_summary$reliability
predicted_value$Reliability_remarks <- target_summary$reliability_remarks
predicted_value$Reliability_percentage <- target_summary$reliability_percentage
predicted_value$Genetic_variance <- target_summary$target_genetic_variance
predicted_value <- gp_bayes_attach_interval_provenance(
  prediction_table = predicted_value,
  model_type = GS_model,
  confidence_level = confidence_level,
  posterior_draw_count = ncol(target_summary$pred_draws)
)

ebv_idx <- if (nrow(sum_ebv) == length(target_summary$genetic_mean)) {
  seq_along(target_summary$genetic_mean)
} else if (!is.null(tst) && nrow(sum_ebv) == length(tst)) {
  tst
} else {
  seq_len(min(nrow(sum_ebv), length(target_summary$genetic_mean)))
}

sum_ebv$Estimated_breeding_value <- target_summary$genetic_mean[ebv_idx]
sum_ebv$Standard_error <- target_summary$standard_error[ebv_idx]
sum_ebv$Prediction_error_variance <- target_summary$prediction_error_var[ebv_idx]
sum_ebv$Reliability <- target_summary$reliability[ebv_idx]
sum_ebv$Genetic_variance <- target_summary$target_genetic_variance[ebv_idx]

diagnostic_plots <- diagnostic_plot_true_prediction(
  boot_results = NULL,
  GID_names = predicted_value[[gen_name]],
  CI_width_thresholds = CI_width_thresholds,
  predictions = predicted_value$Predicted_value,
  Predicted_value_for_CI = target_summary$pred_draws - mod$model$mu,
  mod = mod,
  standard_errors = predicted_value$Standard_error,
  prediction_error_var = predicted_value$PEV,
  genetic_var = target_summary$target_genetic_variance,
  confidence_level = confidence_level,
  model_for_CI_cal = "Bayes",
  high_reliability_thres = high_reliability_thres,
  low_reliability_thres = low_reliability_thres,
  system_database = system_database
)


##############################

if(!is.null(print_lable)){

  print_lable <- unlist(print_lable)

  if(!is.null(geno_data)){
    print_lable <- c("Genomic", print_lable)
  }
  if(length(print_lable)>length(m_matrix_model_ready_list) | length(print_lable)<length(m_matrix_model_ready_list)){
    message("The number of omics labels does not match the number of data layers; default layer names were used.")



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

  if (length(m_matrix_model_ready_list) > 1L) message("No omics labels were provided; default layer names were used.")


}


# if(is.null(tst) | length(tst)==0){
#
#   res <- list(Coefficients = coefficients_list,
#               Estimated_breeding_value = estimated_breeding_value_list,
#               Total_estimated_breeding_value = sum_ebv,
#               Predicted_value =  predicted_value,
#               Residual_value = residual_value,
#               Variance_components = variance_components,
#               M_matrix_model_ready =  m_matrix_model_ready_list,
#               diagnostic_plots  = NULL
#
#   )
#
# } else {
#   if(!is.null(diagnostic_plots)){
  marker_source_labels <- gp_bayes_eta_source_labels(ETA, random_eta_idx)
  marker_model_parameters <- data.frame(
    stat = c(
      "bayesian_marker_model",
      "bayesian_feature_block_count",
      "bayesian_feature_block_names",
      "bayesian_feature_block_strategy",
      "bayesian_total_genetic_variance_estimand"
    ),
    summary = c(
      as.character(GS_model),
      as.character(length(marker_source_labels)),
      paste(marker_source_labels, collapse = ";"),
      "separate_BGLR_marker_ETA_components_summed_for_prediction",
      "posterior_variance_of_summed_genetic_effects"
    ),
    stringsAsFactors = FALSE
  )
  res <- list(Coefficients = coefficients_list,
              Estimated_breeding_value = estimated_breeding_value_list,
              Total_estimated_breeding_value = sum_ebv,
              Predicted_value =  predicted_value,
              Residual_value = residual_value,
              Variance_components = variance_components,
              model_parameters = marker_model_parameters,
              M_matrix_model_ready =  m_matrix_model_ready_list,
              diagnostic_tst_plot  = diagnostic_plots)
#   } else {
#     res <- list(Coefficients = coefficients_list,
#                 Estimated_breeding_value = estimated_breeding_value_list,
#                 Total_estimated_breeding_value = sum_ebv,
#                 Predicted_value =  predicted_value,
#                 Residual_value = residual_value,
#                 Variance_components = variance_components,
#                 M_matrix_model_ready =  m_matrix_model_ready_list,
#                 diagnostic_plots  = NULL
#                 )
#   }
# } # sik
 ### remove the generated output files from the working directory
 unlink(mod[["output_files_names"]])
return(res)

}
#####
########






