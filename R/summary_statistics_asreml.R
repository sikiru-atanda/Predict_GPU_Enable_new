
#' Title
#'
#' @param mod
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary_statistics_asreml <- function(mod=NULL,
                                      response = NULL,
                                      GID_names = NULL,
                                      pheno_data = NULL,
                                      eval_metrics = NULL,
                                      heter_groups = NULL,
                                      predicted_value = NULL,
                                      pred_heter_groups = NULL,
                                      variance_components = NULL,
                                      standard_errors = NULL,
                                      prediction_error_var = NULL,
                                      genetic_var = NULL,
                                      gen_name = NULL,
                                      CI_width_thresholds = c(0.33, 0.66),
                                      confidence_level = 0.95,
                                      high_reliability_thres = 0.7,
                                      low_reliability_thres = 0.4,
                                      system_database = FALSE,
                                      ...){

##browser()
  if(!is.null(eval_metrics)){
    Eval_met <- matrix(NA, nrow = length(eval_metrics), ncol = 1)

    rownames(Eval_met) <- eval_metrics

  } else {
    Eval_met <- NULL
  }

  #pheno <-  as.data.frame(mod$mf)
  n_pheno <- sum(!is.na(pheno_data[, response]))
  trn_min <-  round(min(pheno_data[, response],na.rm=TRUE), 3)
  trn_max <- round(max(pheno_data[, response],na.rm=TRUE), 3)
  var_trn <- round(var(pheno_data[, response],na.rm=TRUE),3)
  n_trn <- n_pheno
  n_tst <- 0

  if(is.null(heter_groups)){
  Res_trn <- round(variance_components["residual_variance", 1],3)

  } else {
    Res_trn <-  NA
  }

  if(nrow(pheno_data)==length(unique(pheno_data[, gen_name]))){


  n<-length(pheno_data[, response])

  if(any(is.na(pheno_data[, response]))){
    tst <- which(is.na(pheno_data[, response]))

    n_trn <- n-length(tst)

    n_tst <- length(tst)

    Eval_met <- NULL
    ###
    # if(!is.null(eval_metrics)){
    #   for (i in 1:length(eval_metrics)){
    #
    #     Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_data[, response][tst],
    #                                         y_predicted = predicted_value$Predicted_value[tst],
    #                                         eval_metrics = eval_metrics[i])
    #
    #   }
    #
    # }

    if(isFALSE(anyNA(standard_errors))){
    diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = GID_names[tst],
                                                           CI_width_thresholds = CI_width_thresholds,
                                                           predictions = predicted_value[tst],
                                                           standard_errors = standard_errors[tst],
                                                           prediction_error_var = prediction_error_var[tst],
                                                           genetic_var = genetic_var,
                                                           confidence_level = confidence_level,
                                                           model_for_CI_cal = "GBLUP",
                                                           #threshold = NULL,
                                                           high_reliability_thres = high_reliability_thres,
                                                           low_reliability_thres = low_reliability_thres)

    }

  }else{

    n_trn <- n

    n_tst <- 0
    if(isFALSE(anyNA(standard_errors))){
      diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = GID_names,
                                                             CI_width_thresholds = CI_width_thresholds,
                                                             predictions = predicted_value,
                                                             standard_errors = standard_errors,
                                                             prediction_error_var = prediction_error_var,
                                                             genetic_var = genetic_var,
                                                             confidence_level = confidence_level,
                                                             model_for_CI_cal = "GBLUP",
                                                             #threshold = NULL,
                                                             high_reliability_thres = high_reliability_thres,
                                                             low_reliability_thres = low_reliability_thres)

    }
####
    if(!is.null(eval_metrics)){
      for (i in 1:length(eval_metrics)){

        Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_data[, response],
                                            y_predicted = predicted_value$Predicted_value,
                                            eval_metrics = eval_metrics[i])

      }

    }

  }

  } else {

    Eval_met <- NULL
}
####
  Stat_Res <-  as.data.frame(t(data.frame(Min = trn_min,
                                        Max = trn_max,
                                        Phenotype_Variance = var_trn,
                                        Residual_Variance = Res_trn,
                                        Number_TrainingSet = n_trn,
                                        Number_TestingSet = n_tst
                                        #GS_model = model
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


  output <-  list(summary_statistics = Stat_Res,
                  diagnostic_tst_plot = diagnostic_tst_plot
  )


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
#     graphics::plot(mod$y[tst]~I(mod$yHat[tst]),ylab="Fitted Value",
#                    xlab="Predicted Value" ,cex=1,bty="L")
#     graphics::points(y=mod$y[tst],x=mod$yHat[tst],col=c("red", 'blue'),cex=1,pch=21)
#     #points(y=mod$y,x=mod$yHat,col=c("red", 'blue'),cex=1,pch=21)
#     graphics::legend("topleft", legend=c("testing", "training"),bty="n",
#                      pch=c(1,19), col=c("red","blue"))
#     #abline(lm(I(mod$y[-tst])~I(mod$yHat[-tst]))$coef,col=1,lwd=2)
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
#
#
#
#
#
