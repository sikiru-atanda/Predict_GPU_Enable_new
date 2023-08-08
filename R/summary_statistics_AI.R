
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
summary_statistics_AI <- function(mod=NULL,
                                  pheno_object= NULL,
                                  response = NULL,
                                  test_set = NULL,
                                  geno_omic_object=NULL,
                                  eval_metrics = c("Accuracy",
                                                   "Mean_Squared_Error",
                                                   "Bias",
                                                   "Root_Mean_Squared_Error",
                                                   "Relative_Squared_Error",
                                                   "Mean_Absolute_Error",
                                                   "Mean_Absolute_Percent_Error"),
                                  GS_model = NULL,
                                                   ...){

  Eval_met <- matrix(NA, nrow = length(eval_metrics), ncol = 1)

  rownames(Eval_met) <- eval_metrics

  #n_pheno <- paste('Number of phenotypes=', (nrow(pheno_object)))

  n_pheno <- nrow(pheno_object)

  # if(GS_model=="Xgboost"){
  # nfeatures <-  mod$trained_model$nfeatures
  #
  # feature_names <- mod$trained_model$feature_names
  # }

  nfeatures = ncol(geno_omic_object)
  feature_names = colnames(geno_omic_object)

rm(geno_omic_object)
  yhat <-  mod$predicted_values
  #Res <-  cat(tmp,'\n')


  #cat(' Min (Traning set)', response, '= ', min(mod$y,na.rm=TRUE),'\n')
  #cat(' Min', paste0(paste0("(",response),')'), '= ', min(mod$y,na.rm=TRUE),'\n')
  trn_min <- paste(paste('Min', '= '), round(min(pheno_object[, response],na.rm=TRUE), 3), sep = "")
  #cat(' Max (Traning set)',response, '= ', max(mod$y,na.rm=TRUE),'\n')
  #cat(' Max', paste0(paste0("(",response),')'), '= ', max(mod$y,na.rm=TRUE),'\n')
  trn_max <- paste(paste('Max', '= '), round(max(pheno_object[, response],na.rm=TRUE), 3), sep = "")
  #cat(' Variance of phenotypes (TRN)=', round(var(mod$y,na.rm=TRUE),4),'\n')
  var_trn <- paste('Variance of phenotypes =', round(var(pheno_object[, response],na.rm=TRUE),3))
  Res_trn <- paste('Residual variance=',round(var(pheno_object[, response] - yhat[,1]),3))

  #n<-length(mod$model$y)

  if(!is.null(test_set))
  {


    n_trn <- paste('Number of Traning =',n_pheno)

    n_tst <- paste('Number of Testing =',nrow(test_set))

    #pred_acc <-  paste('Prediction Accuarcy =',round(cor(pheno_object[, response],yhat[,1]),3))

    for (i in 1:length(eval_metrics)){

      Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[, response],
                                          y_predicted = yhat,
                                          eval_metrics = eval_metrics[i])

    }

  }else{

    n_trn <- paste('Number of Traning =', n_pheno)

    n_tst <- paste('Number of Testing =',0)

    #pred_acc <- paste('Prediction Accu of Training =',round(cor(pheno_object[, response],yhat[, 1]),3))

    for (i in 1:length(eval_metrics)){

      Eval_met[i, ] <- evaluation_metrics(y_observed = pheno_object[, response],
                                          y_predicted = yhat[, 1],
                                          eval_metrics = eval_metrics[i])

    }

  }


  Stat_Res = as.data.frame(t(data.frame(Min = trn_min,
                                        Max = trn_max,
                                        Variance = var_trn,
                                        Residual = Res_trn,
                                        n_features = nfeatures,
                                        n_trn = n_trn,
                                        n_tst = n_tst
                                        )))
  Stat_Res$stat = rownames( Stat_Res)
  Stat_Res =  Stat_Res[, c(2,1)]
  names(Stat_Res)[2] <- "summary"

  Eval_met = data.frame(Eval_met)
  Eval_met$stat = rownames(Eval_met)
  Eval_met =  Eval_met[, c(2,1)]
  names(Eval_met)[2] <- "summary"

  Stat_Res = rbind(Stat_Res, Eval_met)
  #####
if(GS_model=="Xgboost"){
  model_para = do.call(rbind, mod$trained_model$params[-c(8, 9)])
  model_para = data.frame(rbind(model_para, eval_metric = mod$model_parameters[[9]][1]))
  model_para$Parameters = rownames(model_para)
  model_para = model_para[, c(2,1)]
  names(model_para)[2] <- "Value"

}

  if(GS_model=="RandomForest" | GS_model== "K-NearestNeighbors" | GS_model== "SupportVectorMachine" | GS_model=="Lasso" | GS_model=="Ridge_Regression"){

    model_para = data.frame(Value= mod$model_parameters)
    model_para$Parameters = rownames(model_para)
    model_para = model_para[, c(2,1)]


  }

 #  output <- list(Min = trn_min, Max = trn_max, Variance = var_trn,
 #                 Residual = Res_trn, n_trn = n_trn, n_tst = n_tst,
 #                 pred_acc = pred_acc,
 #                 model_para = mod$model_results$model_parameters[-length(mod$model_results$model_parameters)]
 # )
 #
 #
 #
 #  names(output) <-  c("trn_min", "trn_max", "variance_trn",
 #                      "Residual", "n_trn", "n_tst",
 #                      "pred_acc", "model_para")


  output <-  list(Stat_Res, model_para)

  names(output) <- c("Statics_summary", "model_parameters")

  return(output)

}



#' Title
#'
#' @param mod
#' @param pheno_object
#' @param response
#' @param test_set
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_acc_AI <-function(mod=NULL,
                    pheno_object= NULL,
                    response = NULL,
                    test_set = NULL,
                    GS_model = NULL,
                    ...){

  yhat <-  mod$predicted_values
  # DT_ <- data.frame(y = c(test_set[, response],yhat),yhat = c(yhat, test_set[, response]))
  #   # Scatter plot by group
  # ggplot2::ggplot(DT_, aes(x = y, y = yhat)) +
  #   ggplot2::geom_point()+
  #   ggplot2::geom_smooth(method="lm") +
  #   ggpmisc::stat_poly_line() +
  #   ggpmisc::stat_poly_eq(use_label(c("R2")))
  #

  #DT = data.frame(y = mod$y, yhat = mod$yHat)
  if(!is.null(test_set))
  {

  # grDevices::tiff(file="saving_plot3.tiff", units="in",
  #                   width=8, height=5, res=300)


  graphics::plot(test_set[, response]~I(yhat),ylab="Fitted Value",
                         xlab="Predicted Value" ,cex=1,bty="L")
    graphics::points(y=test_set[, response],x=yhat,col=c("red", 'blue'),cex=1,pch=21)
    #points(y=mod$y,x=mod$yHat,col=c("red", 'blue'),cex=1,pch=21)
    graphics::legend("topleft", legend=c("testing", "training"),bty="n",
                     pch=c(1,19), col=c("red","blue"))
    #abline(lm(I(mod$y[-tst])~I(mod$yHat[-tst]))$coef,col=1,lwd=2)
    graphics::abline(stats::lm(I(test_set[, response])~I(yhat))$coef,col=2,lwd=2)

    grDevices::dev.off()

  }


  return(grDevices::jpeg(filename=paste(paste(paste(GS_model, response, sep="_"), "predAccuracy", sep="_"), "jpg", sep = "."), units="in",
                               width=8, height=5, res=300))

  grDevices::dev.off()

}





