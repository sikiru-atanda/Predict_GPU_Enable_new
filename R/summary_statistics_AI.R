
#' Summarize Classical ML Prediction Outputs
#'
#' For classical machine-learning and deep-learning Gaussian prediction tables
#' in this package, inferential `Standard_error` and `PEV` use held-out
#' calibration when available. `Prediction_stability` is the bounded
#' marker-adjustment score `max(0, min(1, 1 - SE_i^2 / var(yhat_all)))`, where
#' `yhat_all` contains all final prediction rows labeled `Train` or `Test`;
#' the ML/DL `Reliability` field is its plotting-compatible, explicitly
#' non-genetic alias. `Prediction_confidence` and
#' `Prediction_risk_score` are calibrated decision-support summaries derived
#' from the calibrated prediction variance scale, not from a Bayesian or
#' mixed-model genetic variance decomposition. `Rank_instability_risk` is the
#' breeder-facing training-calibrated ranking-risk score. `Prediction_unfamiliarity`
#' is a predictor-space novelty score that reflects how far a genotype lies
#' from the training feature distribution. It is used as an additional
#' unfamiliarity channel for the final Gaussian ML/DL prediction-risk score,
#' but it does not overwrite the separate ranking-risk channel.
#' `Prediction_risk_basis` indicates whether a row's displayed risk came from
#' held-out training OOB behavior (`Train`) or from the deployed calibrated
#' final-prediction scorer (`Test`).
#'
#' @param predicted_object Predicted-values frame from a classical-ML / DL
#'   wrapper (rows = individuals, with `Predicted_value`, `Standard_error`,
#'   `PEV`, etc.).
#' @param pheno_object Phenotype object carrying the cleaned phenotypes used
#'   to score the predictions.
#' @param response Name of the response (trait) column in `pheno_object`.
#' @param test_set Integer / character vector identifying the test-set rows.
#' @param geno_omic_object Genomic / omics feature matrix used for the fit
#'   (passed through for unfamiliarity / OOD scoring).
#' @param model_parameters Named list of model hyper-parameters from the
#'   underlying ML / DL fit, included in the summary report.
#' @param eval_metrics Character vector of metrics to compute (e.g.
#'   `"accuracy"`).
#' @param response_family Response distribution family
#'   (`"gaussian"`, `"binomial"`, `"multinomial"`, `"ordinal"`).
#' @param GS_model Model name string (used for labels / routing).
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A list of accuracy / uncertainty / reliability / risk summary
#'   tables and diagnostic plots.
#' @export
#'
#' @examples
summary_statistics_AI <- function(predicted_object= NULL,
                                  pheno_object= NULL,
                                  response = NULL,
                                  test_set = NULL,
                                  geno_omic_object=NULL,
                                  model_parameters = NULL,
                                  eval_metrics = NULL,
                                  response_family = "gaussian",
                                  GS_model = NULL,
                                   ...){

#browser()
  if(!is.null(eval_metrics)){
    Eval_met <- matrix(NA, nrow = length(eval_metrics), ncol = 1)

    rownames(Eval_met) <- eval_metrics

  } else {
    Eval_met <- NULL
  }




  n_pheno <- nrow(pheno_object)


  nfeatures <- if (is.null(geno_omic_object)) {
    NA_integer_
  } else {
    ncol(geno_omic_object)
  }
  #feature_names = colnames(geno_omic_object)
  yhat <-  predicted_object
  prob_cols <- grep("^Prob_", colnames(yhat), value = TRUE)
  #Res <-  cat(tmp,'\n')

  prediction_labels <- if ("Train_Test_Label" %in% names(yhat)) {
    trimws(as.character(yhat[["Train_Test_Label"]]))
  } else {
    character()
  }
  if (length(prediction_labels) == nrow(yhat) &&
      any(prediction_labels %in% c("Train", "Test"))) {
    n_trn <- sum(prediction_labels == "Train", na.rm = TRUE)
    n_tst <- sum(prediction_labels == "Test", na.rm = TRUE)
  } else {
    n_trn <- n_pheno
    n_tst <- if (is.null(test_set)) {
      0L
    } else if (inherits(test_set, "data.frame") || inherits(test_set, "matrix")) {
      nrow(test_set)
    } else if (is.atomic(test_set) && is.null(dim(test_set))) {
      length(test_set)
    } else {
      0L
    }
  }

  fam <- gp_resolve_response_family(response_family, y = pheno_object[[response]])
  if (identical(fam, "gaussian")) {
    trn_min <- round(min(pheno_object[[response]],na.rm=TRUE), 3)
    trn_max <- round(max(pheno_object[[response]],na.rm=TRUE), 3)
    var_trn <- round(var(pheno_object[[response]],na.rm=TRUE),3)
  } else {
    trn_min <- NA
    trn_max <- NA
    var_trn <- NA
  }

  if(!is.null(test_set)){
  Res_trn <- NA

  }else {
  Res_trn <- if (identical(fam, "gaussian")) {
    round(var(pheno_object[[response]] - yhat[,"Predicted_value"], na.rm = TRUE),3)
  } else {
    NA
  }
  }
  #n<-length(mod$model$y)

  if(!is.null(test_set) & !is.null(eval_metrics)){

    Eval_met <- NULL
    #if(!is.null(eval_metrics)){
    # for (i in 1:length(eval_metrics)){
    #
    #   Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[, response],
    #                                       y_predicted = yhat[, "Predicted_value"],
    #                                       eval_metrics = eval_metrics[i])
    #
    # }

    #}


  }else{

    #pred_acc <- paste('Prediction Accu of Training =',round(cor(pheno_object[, response],yhat[, 1]),3))

    if(!is.null(eval_metrics)){
    for (i in 1:length(eval_metrics)){
      y_pred_eval <- if (fam == "multiclass" && tolower(eval_metrics[i]) == "log_loss" && length(prob_cols) > 1L) {
        as.matrix(yhat[, prob_cols, drop = FALSE])
      } else if (fam == "binary" && tolower(eval_metrics[i]) %in% c("log_loss", "brier_score") && length(prob_cols) == 1L) {
        yhat[, prob_cols, drop = TRUE]
      } else {
        yhat[, "Predicted_value"]
      }
      Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[[response]],
                                          y_predicted = y_pred_eval,
                                          eval_metrics = eval_metrics[i],
                                          response_family = response_family)

    }

    }

  }


  Stat_Res <- as.data.frame(t(data.frame(Min = trn_min,
                                        Max = trn_max,
                                        Phenotype_Variance = var_trn,
                                        Residual_Variance = Res_trn,
                                        Number_TrainingSet = n_trn,
                                        Number_TestingSet = n_tst,
                                        Number_Predictors = nfeatures,
                                        Response_Family = fam,
                                        GS_model = GS_model
                                        )))
  Stat_Res$stat <- rownames( Stat_Res)
  Stat_Res <- Stat_Res[, c(2,1)]
  names(Stat_Res)[2] <- "summary"
  rownames(Stat_Res) <- NULL
  #####
  if(!is.null(eval_metrics) & !is.null(Eval_met)){
    Eval_met <- data.frame(Eval_met)
    Eval_met$stat <- rownames(Eval_met)
    Eval_met <- Eval_met[, c(2,1)]
    names(Eval_met)[2] <- "summary"
    rownames(Eval_met) <- NULL

    ### This is specific where the pheno_data has missing value but no test_value for the
    ## validation exercise. proper validation exercise should be done with validation fxn
    Eval_met <- Eval_met[complete.cases(Eval_met), ]
    ####
    Stat_Res <- rbind(Stat_Res, Eval_met)

  }

  if(!is.null(model_parameters)){
    # mm <-  data.frame(stat = names(model_parameters),
    #                 summary = unlist(model_parameters))
    # Check if any entry in the summary column matches the specific type when number of neuron is greater than 1
    if(is.data.frame(model_parameters)){
      if (all(c("stat", "summary") %in% names(model_parameters))) {
        mm_rows <- lapply(seq_len(nrow(model_parameters)), function(i) {
          values <- unlist(model_parameters[["summary"]][[i]], use.names = FALSE)
          if (!length(values)) values <- NA_character_
          data.frame(
            stat = rep(as.character(model_parameters[["stat"]][[i]]), length(values)),
            summary = as.character(values),
            stringsAsFactors = FALSE
          )
        })
        mm <- do.call(rbind, mm_rows)
      } else {
        mm <- data.frame(
          stat = names(model_parameters),
          summary = vapply(model_parameters, function(x) paste(x, collapse = ","), character(1)),
          stringsAsFactors = FALSE
        )
      }
    } else {
      mm <-  data.frame(stat = names(model_parameters),
                        summary = unlist(model_parameters))
  }


    rownames(mm) <- NULL
    Stat_Res <- rbind(Stat_Res, mm)
  }
  #####

  if (identical(fam, "gaussian") &&
      all(c("Standard_error", "PEV", "Prediction_stability", "Prediction_confidence", "Prediction_risk_score") %in% names(yhat))) {
    cal_attr <- attr(yhat, "uncertainty_calibration")
    cal_rows <- NULL
    risk_rows <- NULL
    if (is.list(cal_attr)) {
      cal_rows <- data.frame(
        stat = c(
          "Uncertainty_Calibration_Factor",
          "Uncertainty_Calibration_RMSE",
          "Uncertainty_Calibration_Mean_PEV",
          "Uncertainty_Calibration_Empirical_Coverage"
        ),
        summary = as.character(c(
          cal_attr$calibration_factor,
          cal_attr$calibration_rmse,
          cal_attr$calibration_mean_pev,
          cal_attr$empirical_coverage
        )),
        stringsAsFactors = FALSE
      )
    }
    if (all(c("Train_Test_Label", "Prediction_risk_remarks", "Prediction_risk_basis") %in% names(yhat))) {
      train_basis <- unique(yhat$Prediction_risk_basis[yhat$Train_Test_Label != "Test"])
      test_basis <- unique(yhat$Prediction_risk_basis[yhat$Train_Test_Label == "Test"])
      train_basis <- train_basis[!is.na(train_basis) & nzchar(train_basis)]
      test_basis <- test_basis[!is.na(test_basis) & nzchar(test_basis)]

      train_risk <- yhat$Prediction_risk_remarks[yhat$Train_Test_Label != "Test"]
      test_risk <- yhat$Prediction_risk_remarks[yhat$Train_Test_Label == "Test"]
      risk_levels <- c("Low Risk", "Moderate Risk", "High Risk")
      risk_percentage <- function(x, level) {
        usable <- !is.na(x) & nzchar(as.character(x))
        if (!any(usable)) return(NA_character_)
        sprintf("%.2f", 100 * mean(x[usable] == level))
      }

      risk_rows <- data.frame(
        stat = c(
          "Train_Risk_Basis",
          "Test_Risk_Basis",
          paste0("Train_", gsub(" ", "_", risk_levels), "_Percentage"),
          paste0("Test_", gsub(" ", "_", risk_levels), "_Percentage")
        ),
        summary = c(
          if (length(train_basis)) paste(train_basis, collapse = "; ") else NA_character_,
          if (length(test_basis)) paste(test_basis, collapse = "; ") else NA_character_,
          risk_percentage(train_risk, risk_levels[1]),
          risk_percentage(train_risk, risk_levels[2]),
          risk_percentage(train_risk, risk_levels[3]),
          risk_percentage(test_risk, risk_levels[1]),
          risk_percentage(test_risk, risk_levels[2]),
          risk_percentage(test_risk, risk_levels[3])
        ),
        stringsAsFactors = FALSE
      )
    }
    unfamiliarity_rows <- NULL
    if (all(c("Train_Test_Label", "Prediction_unfamiliarity", "Prediction_unfamiliarity_remarks", "Prediction_unfamiliarity_basis") %in% names(yhat))) {
      train_unf_basis <- unique(yhat$Prediction_unfamiliarity_basis[yhat$Train_Test_Label != "Test"])
      test_unf_basis <- unique(yhat$Prediction_unfamiliarity_basis[yhat$Train_Test_Label == "Test"])
      train_unf_basis <- train_unf_basis[!is.na(train_unf_basis) & nzchar(train_unf_basis)]
      test_unf_basis <- test_unf_basis[!is.na(test_unf_basis) & nzchar(test_unf_basis)]

      train_unf <- suppressWarnings(as.numeric(yhat$Prediction_unfamiliarity[yhat$Train_Test_Label != "Test"]))
      test_unf <- suppressWarnings(as.numeric(yhat$Prediction_unfamiliarity[yhat$Train_Test_Label == "Test"]))
      train_unf_remark <- yhat$Prediction_unfamiliarity_remarks[yhat$Train_Test_Label != "Test"]
      test_unf_remark <- yhat$Prediction_unfamiliarity_remarks[yhat$Train_Test_Label == "Test"]
      unfamiliarity_levels <- c("Familiar", "Moderately Unfamiliar", "Highly Unfamiliar")

      unfamiliarity_rows <- data.frame(
        stat = c(
          "Train_Prediction_Unfamiliarity_Basis",
          "Test_Prediction_Unfamiliarity_Basis",
          "Train_Mean_Prediction_Unfamiliarity",
          "Test_Mean_Prediction_Unfamiliarity",
          paste0("Train_", gsub(" ", "_", unfamiliarity_levels), "_Percentage"),
          paste0("Test_", gsub(" ", "_", unfamiliarity_levels), "_Percentage")
        ),
        summary = c(
          if (length(train_unf_basis)) paste(train_unf_basis, collapse = "; ") else NA_character_,
          if (length(test_unf_basis)) paste(test_unf_basis, collapse = "; ") else NA_character_,
          if (length(train_unf)) sprintf("%.4f", mean(train_unf, na.rm = TRUE)) else NA_character_,
          if (length(test_unf)) sprintf("%.4f", mean(test_unf, na.rm = TRUE)) else NA_character_,
          if (length(train_unf_remark)) sprintf("%.2f", 100 * mean(train_unf_remark == unfamiliarity_levels[1], na.rm = TRUE)) else NA_character_,
          if (length(train_unf_remark)) sprintf("%.2f", 100 * mean(train_unf_remark == unfamiliarity_levels[2], na.rm = TRUE)) else NA_character_,
          if (length(train_unf_remark)) sprintf("%.2f", 100 * mean(train_unf_remark == unfamiliarity_levels[3], na.rm = TRUE)) else NA_character_,
          if (length(test_unf_remark)) sprintf("%.2f", 100 * mean(test_unf_remark == unfamiliarity_levels[1], na.rm = TRUE)) else NA_character_,
          if (length(test_unf_remark)) sprintf("%.2f", 100 * mean(test_unf_remark == unfamiliarity_levels[2], na.rm = TRUE)) else NA_character_,
          if (length(test_unf_remark)) sprintf("%.2f", 100 * mean(test_unf_remark == unfamiliarity_levels[3], na.rm = TRUE)) else NA_character_
        ),
        stringsAsFactors = FALSE
      )
    }
    Stat_Res <- rbind(
      Stat_Res,
      data.frame(
        stat = c(
          "Standard_Error_Definition",
          "PEV_Definition",
          "Prediction_Stability_Definition",
          "Prediction_Confidence_Definition",
          "Prediction_Risk_Score_Definition",
          "Prediction_Unfamiliarity_Definition",
          "Rank_Instability_Risk_Definition",
          "Prediction_Risk_Basis_Definition"
        ),
        summary = c(
          "Target-specific predictive standard error when substantively identifiable; otherwise NA while marginal interval calibration is retained",
          "Target-specific predictive error variance when substantively identifiable; otherwise NA rather than repeating marginal cross-fitted MSE",
          "Bootstrap consistency score computed from the raw bootstrap prediction variance scale",
          "Predictive stability computed as Var(y_train) / (Var(y_train) + PEV) using the calibrated prediction variance scale",
          "Combined predictive-risk score that blends held-out calibrated ranking risk with predictor-space unfamiliarity",
          "Predictor-space novelty score computed from local distance to the training feature distribution",
          "Breeder-facing ranking-risk score learned from training-only held-out behavior and reported on the same scale as Prediction_risk_score",
          "Prediction_risk_basis records Heldout_Training_OOB, Calibrated_Final_Prediction, or explicit unavailability; availability depends on the fitted backend"
        ),
        stringsAsFactors = FALSE
      ),
      cal_rows,
      risk_rows,
      unfamiliarity_rows
    )
  }



  output <-  list(Stat_Res)

  names(output) <- "summary_statistics"

  return(output)

}




# plot_acc_AI <-function(mod=NULL,
#                     pheno_object= NULL,
#                     response = NULL,
#                     test_set = NULL,
#                     GS_model = NULL,
#                     ...){
#
#   yhat <-  mod$predicted_values
#   # DT_ <- data.frame(y = c(test_set[, response],yhat),yhat = c(yhat, test_set[, response]))
#   #   # Scatter plot by group
#   # ggplot2::ggplot(DT_, aes(x = y, y = yhat)) +
#   #   ggplot2::geom_point()+
#   #   ggplot2::geom_smooth(method="lm") +
#   #   ggpmisc::stat_poly_line() +
#   #   ggpmisc::stat_poly_eq(use_label(c("R2")))
#   #
#
#   #DT = data.frame(y = mod$y, yhat = mod$yHat)
#   if(!is.null(test_set))
#   {
#
#   # grDevices::tiff(file="saving_plot3.tiff", units="in",
#   #                   width=8, height=5, res=300)
#
#
#   graphics::plot(test_set[, response]~I(yhat),ylab="Fitted Value",
#                          xlab="Predicted Value" ,cex=1,bty="L")
#     graphics::points(y=test_set[, response],x=yhat,col=c("red", 'blue'),cex=1,pch=21)
#     #points(y=mod$y,x=mod$yHat,col=c("red", 'blue'),cex=1,pch=21)
#     graphics::legend("topleft", legend=c("testing", "training"),bty="n",
#                      pch=c(1,19), col=c("red","blue"))
#     #abline(lm(I(mod$y[-tst])~I(mod$yHat[-tst]))$coef,col=1,lwd=2)
#     graphics::abline(stats::lm(I(test_set[, response])~I(yhat))$coef,col=2,lwd=2)
#
#     grDevices::dev.off()
#
#   }
#
#
#   return(grDevices::jpeg(filename=paste(paste(paste(GS_model, response, sep="_"), "predAccuracy", sep="_"), "jpg", sep = "."), units="in",
#                                width=8, height=5, res=300))
#
#   grDevices::dev.off()
#
# }
#
#
#
#
#
