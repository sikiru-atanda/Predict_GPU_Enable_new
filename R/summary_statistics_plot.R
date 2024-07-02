
#' Title
#'
#' @param mod
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary_statistics_bayes <- function(mod=NULL,
                                     eval_metrics = NULL,
                                     GS_model = NULL,
                                     model_result = NULL,
                                     gen_name = NULL,
                                     CI_width_thresholds = c(0.33, 0.66),
                                     confidence_level = 0.95,
                                     high_reliability_thres = 0.7,
                                     low_reliability_thres = 0.4,
                                     system_database = FALSE,
                                     ...){

  tst <- NULL
  if(!is.null(eval_metrics)){
  Eval_met <- matrix(NA, nrow = length(eval_metrics), ncol = 1)

  rownames(Eval_met) <- eval_metrics

  } else {

    Eval_met <- NULL
  }

  n_pheno <- sum(!is.na(mod$model$y))
  trn_min <-  round(min(mod$model$y,na.rm=TRUE), 3)
  trn_max <- round(max(mod$model$y,na.rm=TRUE), 3)
  var_trn <- round(var(mod$model$y,na.rm=TRUE),3)
  Res_trn <- round(mod$model$varE,3)

  n<-length(mod$model$y)

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
    if(isFALSE(anyNA(model_result[["Predicted_value"]]["Standard_error"]))){
 diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = model_result[["Predicted_value"]][gen_name][, 1],
                                                        CI_width_thresholds = CI_width_thresholds,
                                                        predictions = as.double(model_result[["Predicted_value"]]["Predicted_value"][, 1]),
                                                        standard_errors = as.double(model_result[["Predicted_value"]]["Standard_error"][, 1]),
                                                        prediction_error_var = as.double(model_result[["Predicted_value"]]["PEV"][tst, 1]),
                                                        genetic_var = as.double(var(model_result[["Predicted_value"]]["Predicted_value"][, 1])),
                                                        confidence_level = confidence_level,
                                                        model_for_CI_cal = "Bayes",
                                                        #threshold = NULL,
                                                        high_reliability_thres = high_reliability_thres,
                                                        low_reliability_thres = low_reliability_thres)

    }

  }else{

    n_trn <- n

    n_tst <- 0

    if(!is.null(eval_metrics)){
    for (i in 1:length(eval_metrics)){

      Eval_met[i, ] <- evaluation_metrics(y_observed = mod$model$y,
                                          y_predicted = mod$model$yHat,
                                          eval_metrics = eval_metrics[i])

     }

    }
    if(isFALSE(anyNA(model_result[["Predicted_value"]]["Standard_error"]))){
    diagnostic_tst_plot <- diagnostic_plot_true_prediction(GID_names = model_result[["Predicted_value"]][gen_name][,1],
                                                           CI_width_thresholds = CI_width_thresholds,
                                                           predictions = as.double(model_result[["Predicted_value"]]["Predicted_value"][, 1]),
                                                           standard_errors = as.double(model_result[["Predicted_value"]]["Standard_error"][, 1]),
                                                           prediction_error_var = as.double(model_result[["Predicted_value"]]["PEV"][, 1]),
                                                           genetic_var = as.double(var(model_result[["Predicted_value"]]["Predicted_value"])),
                                                           confidence_level = confidence_level,
                                                           model_for_CI_cal = "Bayes",
                                                           #threshold = NULL,
                                                           high_reliability_thres = high_reliability_thres,
                                                           low_reliability_thres = low_reliability_thres)


    }

  }

  #model = data.frame()

  model <- c()

  for(k in 1:length(mod$model$ETA))
  {

        if(!is.null(mod$model$ETA[[k]]$model)){
          #cat(" Coefficientes in ETA[",k,"] (",names(mod$ETA)[k],") modeled as in ", mod$ETA[[k]]$model,"\n")

           #model <- rbind(model, mod$model$ETA[[k]]$model)
           model <- cbind(model, mod$model$ETA[[k]]$model)

        }


  }

  Stat_Res = as.data.frame(t(data.frame(Min = trn_min,
                                        Max = trn_max,
                                        Phenotype_Variance = var_trn,
                                        Residual_Variance = Res_trn,
                                        Number_TrainingSet = n_trn,
                                        Number_TestingSet = n_tst,
                                        GS_model = GS_model
  )))
  Stat_Res$stat <- rownames( Stat_Res)
  Stat_Res <- Stat_Res[, c(2,1)]
  names(Stat_Res)[2] <- "summary"
  rownames(Stat_Res) <- NULL
#####
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

if(is.null(tst) || length(tst)<=1){
  output <-  list(summary_statistics = Stat_Res

                  )
} else {
  output <-  list(summary_statistics = Stat_Res,
                  diagnostic_tst_plot = diagnostic_tst_plot
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





