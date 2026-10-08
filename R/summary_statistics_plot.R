
predictpror_format_bayes_residual_variance <- function(varE, digits = 3) {
  if (is.null(varE) || length(varE) < 1L) {
    return(NA_character_)
  }
  values <- suppressWarnings(as.numeric(varE))
  keep <- is.finite(values)
  if (!any(keep)) {
    return(NA_character_)
  }
  values <- round(values[keep], digits)
  labels <- names(varE)
  if (!is.null(labels) && length(labels) == length(varE)) {
    labels <- labels[keep]
  } else {
    labels <- rep(NA_character_, length(values))
  }
  values_chr <- format(values, trim = TRUE, scientific = FALSE)
  if (length(values_chr) == 1L) {
    return(values_chr)
  }
  has_labels <- !is.na(labels) & nzchar(labels)
  if (any(has_labels)) {
    labels[!has_labels] <- paste0("group_", seq_along(labels))[!has_labels]
    return(paste0(labels, "=", values_chr, collapse = "; "))
  }
  paste(values_chr, collapse = "; ")
}

#' Summarize Bayesian Prediction Outputs
#'
#' Returns a compact summary table for Bayesian model fits, including phenotype
#' range, residual variance, train/test counts, and user-facing definitions for
#' target-specific uncertainty metrics when Bayesian prediction outputs are
#' available.
#'
#' For Bayesian prediction tables in this package, `Genetic_variance`, `PEV`,
#' and `Reliability` are reported for the exact prediction target being shown.
#' In single-environment analyses this is the individual genotype target. In
#' multi-environment analyses this can be either the environment-specific
#' genotype target or the across-environment genotype-average target.
#'
#' A near-zero reliability, including for `BayesB`, does not by itself imply the
#' predictions are wrong. It means that the posterior target-specific prediction
#' error variance is close to or larger than the corresponding target-specific
#' genetic variance under the fitted model and data.
#'
#' @param mod A fitted Bayesian model object returned by the package.
#' @param eval_metrics Optional evaluation metrics.
#' @param GS_model Model name used for the fit.
#' @param model_result Extracted package prediction output for the model.
#' @param gen_name Genotype identifier column name.
#' @param CI_width_thresholds Confidence-interval width thresholds.
#' @param confidence_level Confidence level used for uncertainty summaries.
#' @param high_reliability_thres High reliability cutoff.
#' @param low_reliability_thres Low reliability cutoff.
#' @param system_database Logical; controls export behavior in higher-level calls.
#' @param heter_groups Optional environment/location grouping column.
#' @param ... Additional unused arguments.
#'
#' @return A list containing a `summary_statistics` data frame.
#' @export
#'
#' @examples
summary_statistics_bayes <- function(mod=NULL,
                                     eval_metrics = NULL,
                                     GS_model = NULL,
                                     model_result = NULL,
                                     gen_name = NULL,
                                     response_family = "gaussian",
                                     CI_width_thresholds = c(0.33, 0.66),
                                     confidence_level = 0.95,
                                     high_reliability_thres = 0.7,
                                     low_reliability_thres = 0.4,
                                     system_database = FALSE,
                                     heter_groups = NULL,
                                     ...){
  fam <- gp_resolve_response_family(response_family, y = mod$model$y)

  tst <- NULL
  if(!is.null(eval_metrics) & is.null(heter_groups)){
  Eval_met <- matrix(NA, nrow = length(eval_metrics), ncol = 1)

  rownames(Eval_met) <- eval_metrics

  } else {

    Eval_met <- NULL
  }

  n_pheno <- sum(!is.na(mod$model$y))
  trn_min <- if (identical(fam, "gaussian")) round(min(mod$model$y, na.rm = TRUE), 3) else NA
  trn_max <- if (identical(fam, "gaussian")) round(max(mod$model$y, na.rm = TRUE), 3) else NA
  var_trn <- if (identical(fam, "gaussian")) round(var(mod$model$y, na.rm = TRUE), 3) else NA
  Res_trn <- if (identical(fam, "gaussian")) {
    predictpror_format_bayes_residual_variance(mod$model$varE, 3)
  } else {
    NA_character_
  }
  n_trn <- NA
  n_tst <- NA

  n<-length(mod$model$y)

 if(is.null(heter_groups)) {
  if(any(is.na(mod$model$y))){
    tst <- which(is.na(mod$model$y))

    n_trn <- n-length(tst)

    n_tst <- length(tst)
    Eval_met <- NULL
    # if(!is.null(eval_metrics)){
    # for (i in 1:length(eval_metrics)){
    #
    #   Eval_met[i, ] <- evaluation_metrics(y_observed = mod$model$y[tst],
    #                                       y_predicted = mod$model$yHat[tst],
    #                                       eval_metrics = eval_metrics[i])
    #
    # }
    #
    # }
 #    if(isFALSE(anyNA(model_result[["Predicted_value"]]["Standard_error"]))){
 # diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = model_result[["Predicted_value"]][gen_name][, 1],
 #                                                        CI_width_thresholds = CI_width_thresholds,
 #                                                        predictions = as.double(model_result[["Predicted_value"]]["Predicted_value"][, 1]),
 #                                                        standard_errors = as.double(model_result[["Predicted_value"]]["Standard_error"][, 1]),
 #                                                        prediction_error_var = as.double(model_result[["Predicted_value"]]["PEV"][tst, 1]),
 #                                                        genetic_var = as.double(var(model_result[["Predicted_value"]]["Predicted_value"][, 1])),
 #                                                        confidence_level = confidence_level,
 #                                                        model_for_CI_cal = "Bayes",
 #                                                        #threshold = NULL,
 #                                                        high_reliability_thres = high_reliability_thres,
 #                                                        low_reliability_thres = low_reliability_thres)
 #
 #    }

  }else{

    n_trn <- n

    n_tst <- 0

    if(!is.null(eval_metrics)){
    for (i in 1:length(eval_metrics)){

      pred_input <- mod$model$yHat
      if (!identical(fam, "gaussian") && !is.null(model_result) &&
          "Predicted_value" %in% names(model_result) &&
          is.data.frame(model_result[["Predicted_value"]])) {
        pred_df <- model_result[["Predicted_value"]]
        prob_cols <- grep("^Prob_", names(pred_df), value = TRUE)
        if (identical(fam, "binary") && length(prob_cols) >= 1L) {
          pred_input <- pred_df[[prob_cols[length(prob_cols)]]]
        } else if (identical(fam, "ordinal") && length(prob_cols) >= 2L) {
          pred_input <- as.matrix(pred_df[, prob_cols, drop = FALSE])
          colnames(pred_input) <- gsub("^Prob_", "", prob_cols)
        } else if ("Predicted_value" %in% names(pred_df)) {
          pred_input <- pred_df[["Predicted_value"]]
        }
      }

      Eval_met[i, ] <- evaluation_metrics(y_observed = mod$model$y,
                                          y_predicted = pred_input,
                                          eval_metrics = eval_metrics[i],
                                          response_family = fam)

     }

    }
    # if(isFALSE(anyNA(model_result[["Predicted_value"]]["Standard_error"]))){
    # diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = model_result[["Predicted_value"]][gen_name][,1],
    #                                                        CI_width_thresholds = CI_width_thresholds,
    #                                                        predictions = as.double(model_result[["Predicted_value"]]["Predicted_value"][, 1]),
    #                                                        standard_errors = as.double(model_result[["Predicted_value"]]["Standard_error"][, 1]),
    #                                                        prediction_error_var = as.double(model_result[["Predicted_value"]]["PEV"][, 1]),
    #                                                        genetic_var = as.double(var(model_result[["Predicted_value"]]["Predicted_value"])),
    #                                                        confidence_level = confidence_level,
    #                                                        model_for_CI_cal = "Bayes",
    #                                                        #threshold = NULL,
    #                                                        high_reliability_thres = high_reliability_thres,
    #                                                        low_reliability_thres = low_reliability_thres)
    #
    #
    # }

  }

}

  ## put here
  #model = data.frame()

  model <- c()

  if (!is.null(mod$model$ETA) && length(mod$model$ETA) > 0) {
    for(k in seq_along(mod$model$ETA))
    {

          if(!is.null(mod$model$ETA[[k]]$model)){
            #cat(" Coefficientes in ETA[",k,"] (",names(mod$ETA)[k],") modeled as in ", mod$ETA[[k]]$model,"\n")

             #model <- rbind(model, mod$model$ETA[[k]]$model)
             model <- cbind(model, mod$model$ETA[[k]]$model)

          }


    }
  }

  Stat_Res = as.data.frame(t(data.frame(Min = trn_min,
                                        Max = trn_max,
                                        Phenotype_Variance = var_trn,
                                        Residual_Variance = Res_trn,
                                        Number_TrainingSet = n_trn,
                                        Number_TestingSet = n_tst,
                                        Response_Family = fam,
                                        GS_model = GS_model
  )))
  Stat_Res$stat <- rownames( Stat_Res)
  Stat_Res <- Stat_Res[, c(2,1)]
  names(Stat_Res)[2] <- "summary"
  rownames(Stat_Res) <- NULL
#####
  if(is.null(heter_groups)){
  if(!is.null(eval_metrics) && !is.null(Eval_met)){
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

  }

  if (!is.null(model_result) && "Predicted_value" %in% names(model_result) &&
      is.data.frame(model_result[["Predicted_value"]])) {
    pred_df <- model_result[["Predicted_value"]]
    target_variance_col <- if ("Genetic_variance" %in% names(pred_df)) {
      "Genetic_variance"
    } else if ("Reliability_reference_variance" %in% names(pred_df)) {
      "Reliability_reference_variance"
    } else {
      NA_character_
    }
    has_target_metrics <- !is.na(target_variance_col) &&
      all(c("PEV", "Reliability") %in% names(pred_df))
    has_class_metrics <- all(c("Prediction_confidence", "Classification_uncertainty") %in% names(pred_df))

    if (has_target_metrics || has_class_metrics) {
      if (has_target_metrics) {
        target_desc <- if (is.null(heter_groups)) {
          "Individual genotype prediction target"
        } else {
          "Environment-specific genotype prediction target"
        }
        variance_stat <- if (identical(target_variance_col, "Genetic_variance")) {
          "Genetic_Variance_Definition"
        } else {
          "Reliability_Reference_Variance_Definition"
        }
        variance_desc <- if (identical(target_variance_col, "Genetic_variance")) {
          "Posterior target-specific genetic variance used for each reported prediction target"
        } else {
          "Training-response variance used as the GP reference variance for reliability"
        }

        definition_rows <- data.frame(
          stat = c(
            "Prediction_Target",
            variance_stat,
            "PEV_Definition",
            "Reliability_Definition"
          ),
          summary = c(
            target_desc,
            variance_desc,
            "Posterior target-specific prediction error variance for each reported prediction target",
            if (identical(target_variance_col, "Genetic_variance")) {
              "Computed as 1 - PEV / target-specific genetic variance"
            } else {
              "Computed as 1 - PEV / GP reliability reference variance"
            }
          ),
          stringsAsFactors = FALSE
        )
      } else {
        definition_rows <- data.frame(
          stat = c(
            "Prediction_Target",
            "Class_Probability_Definition",
            "Prediction_Confidence_Definition",
            "Classification_Uncertainty_Definition"
          ),
          summary = c(
            if (is.null(heter_groups)) "Individual genotype class target" else "Environment-specific class target",
            "Posterior class probabilities from the fitted Bayesian ordinal model",
            "Maximum posterior class probability for the reported predicted class",
            "Computed as 1 - Prediction_confidence"
          ),
          stringsAsFactors = FALSE
        )
      }

      if ("Total_Predicted_value" %in% names(model_result) &&
          is.data.frame(model_result[["Total_Predicted_value"]])) {
        definition_rows <- rbind(
          definition_rows,
          data.frame(
            stat = "Across_Environment_Target",
            summary = "Across-environment genotype average target reported in Total_Predicted_value",
            stringsAsFactors = FALSE
          )
        )
      }

      Stat_Res <- rbind(Stat_Res, definition_rows)
    }
  }

if((is.null(tst) || length(tst)<=1) && !is.null(heter_groups)){
  output <-  list(summary_statistics = Stat_Res)
} else if((!is.null(tst) || length(tst)<1) && !is.null(heter_groups)){
  output <-  list(summary_statistics = Stat_Res)
} else {
  output <-  list(summary_statistics = Stat_Res
                  #diagnostic_tst_plot = diagnostic_tst_plot
  )
}


  return(output)

}




# plot_acc <- function(mod,response, ...){
#
#   # DT_ <- data.frame(y = c(mod$y,mod$yHat),yhat = c(mod$yhat, mod$y))
#   #   # Scatter plot by group
#   # ggplot2::ggplot(DT_, aes(x = y, y = yhat)) +
#   #   ggplot2::geom_point()+
#   #   ggplot2::geom_smooth(method="lm") +
#   #   ggpmisc::stat_poly_line() +
#   #   ggpmisc::stat_poly_eq(use_label(c("R2")))
#   #
#
#   #DT = data.frame(y = mod$y, yhat = mod$yHat)
#   if(any(is.na(mod$y)))
#   {
#     tst <- which(is.na(mod$y))
#
#     # grDevices::tiff(file="saving_plot3.tiff", units="in",
#     #      width=8, height=5, res=300)
#
#
#  graphics::plot(mod$y[tst]~I(mod$yHat[tst]),ylab="Fitted Value",
#        xlab="Predicted Value" ,cex=1,bty="L")
#     graphics::points(y=mod$y[tst],x=mod$yHat[tst],col=c("red", 'blue'),cex=1,pch=21)
#   #points(y=mod$y,x=mod$yHat,col=c("red", 'blue'),cex=1,pch=21)
#     graphics::legend("topleft", legend=c("testing", "training"),bty="n",
#          pch=c(1,19), col=c("red","blue"))
#   #abline(lm(I(mod$y[-tst])~I(mod$yHat[-tst]))$coef,col=1,lwd=2)
#     graphics::abline(stats::lm(I(mod$y[tst])~I(mod$yHat[tst]))$coef,col=2,lwd=2)
#
#     grDevices::dev.off()
#
#   }
#
#   return(grDevices::jpeg(filename=paste(paste(response, "predAccuracy", sep="_"), "jpg", sep = "."), units="in",
#                          width=8, height=5, res=300))
#
#   grDevices::dev.off()
#
# }





