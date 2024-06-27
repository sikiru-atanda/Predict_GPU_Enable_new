
#' Title
#'
#' @param mod
#' @param ...
#' @param pheno_object
#' @param response
#' @param test_set
#' @param geno_model_ready_train
#' @param eval_metrics
#'
#' @return
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

  nfeatures = ncol(geno_omic_object)
  #feature_names = colnames(geno_omic_object)
  yhat <-  predicted_object
  #Res <-  cat(tmp,'\n')

  trn_min <- round(min(pheno_object[, response],na.rm=TRUE), 3)
  trn_max <- round(max(pheno_object[, response],na.rm=TRUE), 3)
  var_trn <- round(var(pheno_object[, response],na.rm=TRUE),3)
  Res_trn <- round(var(pheno_object[, response] - yhat[,"Predicted_value"]),3)

  #n<-length(mod$model$y)

  if(!is.null(test_set) & !is.null(eval_metrics)){

  ## Number of Traning
    n_trn <- n_pheno
    ## Number of Testing
    n_tst <- nrow(test_set)

    #if(!is.null(eval_metrics)){
    for (i in 1:length(eval_metrics)){

      Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[, response],
                                          y_predicted = yhat[, "Predicted_value"],
                                          eval_metrics = eval_metrics[i])

    }

    #}


  }else{

    n_trn <-  n_pheno

    n_tst <- 0

    #pred_acc <- paste('Prediction Accu of Training =',round(cor(pheno_object[, response],yhat[, 1]),3))

    if(!is.null(eval_metrics)){
    for (i in 1:length(eval_metrics)){

      Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[, response],
                                          y_predicted = yhat[, "Predicted_value"],
                                          eval_metrics = eval_metrics[i])

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
                                        GS_model = GS_model
                                        )))
  Stat_Res$stat <- rownames( Stat_Res)
  Stat_Res <- Stat_Res[, c(2,1)]
  names(Stat_Res)[2] <- "summary"
  rownames(Stat_Res) <- NULL
  #####
  if(!is.null(eval_metrics)){
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
    mm <-  data.frame(stat = names(model_parameters),
                    summary = unlist(model_parameters))
    rownames(mm) <- NULL
    Stat_Res <- rbind(Stat_Res, mm)
  }
  #####



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
