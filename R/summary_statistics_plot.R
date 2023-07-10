
#' Title
#'
#' @param mod
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary_statistics <- function(mod,...){

  n_pheno <-paste('Number of phenotypes=', (sum(!is.na(mod$model$y))))
  #Res <-  cat(tmp,'\n')


  #cat(' Min (Traning set)', response, '= ', min(mod$y,na.rm=TRUE),'\n')
  #cat(' Min', paste0(paste0("(",response),')'), '= ', min(mod$y,na.rm=TRUE),'\n')
  trn_min <- paste(paste('Min', '= '), round(min(mod$y,na.rm=TRUE), 3), sep = "")
  #cat(' Max (Traning set)',response, '= ', max(mod$y,na.rm=TRUE),'\n')
  #cat(' Max', paste0(paste0("(",response),')'), '= ', max(mod$y,na.rm=TRUE),'\n')
  trn_max <- paste(paste('Max', '= '), round(max(mod$y,na.rm=TRUE), 3), sep = "")
  #cat(' Variance of phenotypes (TRN)=', round(var(mod$y,na.rm=TRUE),4),'\n')
  var_trn <- paste('Variance of phenotypes =', round(var(mod$y,na.rm=TRUE),3))
  Res_trn <- paste('Residual variance=',round(mod$varE,3))

  n<-length(mod$y)

  if(any(is.na(mod$y)))
  {
    tst <- which(is.na(mod$y))

    n_trn <- paste('Number of Traning =',n-length(tst))

    n_tst <- paste('Number of Testing =',length(tst))

    pred_acc <-  paste('Prediction Accuarcy =',round(cor(mod$y[-tst],mod$yHat[-tst]),3))

  }else{

    n_trn <- paste('Number of Traning =',n)

    n_tst <- paste('Number of Testing =',0)

    pred_acc <- paste('Prediction Accu of Training =',round(cor(mod$y,mod$yHat),3))

  }


  for(k in 1:length(mod$ETA))
  {

        if(!is.null(names(mod$ETA)[k])){
          #cat(" Coefficientes in ETA[",k,"] (",names(mod$ETA)[k],") modeled as in ", mod$ETA[[k]]$model,"\n")

           model <- names(mod$ETA)[k]

        }


  }

  output <- list(Min = trn_min, Max = trn_max, Variance = var_trn,
                 Residual = Res_trn, n_trn = n_trn, n_tst = n_tst,
                 pred_acc = pred_acc, model_type = model)

  names(output) <-  c("trn_min", "trn_max", "variance_trn",
                      "Residual", "n_trn", "n_tst",
                      "pred_acc", "model_type")

  return(output)

}





# plot_acc <- function(mod,...){
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
#  PL <- plot(mod$y[tst]~I(mod$yHat[tst]),ylab="Fitted Value",
#        xlab="Predicted Value" ,cex=1,bty="L")
#   points(y=mod$y[tst],x=mod$yHat[tst],col=c("red", 'blue'),cex=1,pch=21)
#   #points(y=mod$y,x=mod$yHat,col=c("red", 'blue'),cex=1,pch=21)
#   legend("topleft", legend=c("testing", "training"),bty="n",
#          pch=c(1,19), col=c("red","blue"))
#   #abline(lm(I(mod$y[-tst])~I(mod$yHat[-tst]))$coef,col=1,lwd=2)
#   abline(lm(I(mod$y[tst])~I(mod$yHat[tst]))$coef,col=2,lwd=2)
#
#   }
#
#   return()
#
# }





