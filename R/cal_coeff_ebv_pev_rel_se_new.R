
#' Title
#'
#' ### This function calculate the Coefficient, EBV, PEV, SE, Reliability
#'
#' @param beta  the posterior samples of the random effects (marker/omic)
#' @param x_variable marker/omic data
#' @param gen_name
#' @param var_u genetic variance
#' @param hetero
#' @param heter_groups
#' @param gid_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cal_coeff_ebv_pev_rel_se_new <- function(beta,
                                     x_variable,
                                     gen_name,
                                     var_u,
                                     hetero = NULL,
                                     heter_groups = NULL,
                                     gid_name = NULL,
                                     GS_model = NULL,
                                     ...){

  # the breeding values (posterior)
  X_ebv <- as.matrix(x_variable)%*%t(beta)


  # Posterior means of marker effects/coefficient
  # coeff_omic3 <- colMeans(Bb)

  if(!is.null(heter_groups)){
    Coeff <- data.frame(X_variables = rep(gid_name, length(unique(hetero))),
                        Env = rep(unique(hetero), each= length(gid_name)),
                        Coeff =colMeans(beta),
                        stringsAsFactors = FALSE)
    names(Coeff)[2] <- heter_groups
  } else {

    if(!isSymmetric(x_variable)){
      Coeff <- data.frame(X_variables = colnames(x_variable),
                          Coeff =colMeans(beta),
                          stringsAsFactors = FALSE)

    } else {

      Coeff <- data.frame(X_variables = gid_name,
                          Coeff =colMeans(beta),
                          stringsAsFactors = FALSE)


    }

  }



  if(is.null(hetero)){
    # Genomic estimated breeding values
    EBV <- data.frame(names = gid_name,
                      Estimated_breeding_value=x_variable%*%Coeff[, "Coeff"],
                      stringsAsFactors = FALSE)

  } else {

    EBV <- data.frame(names = rep(gid_name, length(unique(hetero))),
                      Env = rep(unique(hetero), each= length(gid_name)),
                      Estimated_breeding_value=x_variable%*%Coeff[, "Coeff"],
                      stringsAsFactors = FALSE)

    names(EBV)[2] <- heter_groups



  }

  colnames(EBV)[1] = gen_name
  #GEBV = data.frame(rowMeans(ebv))

  PEV <- apply(X_ebv, 1, var)

  # EBV$Std_error <- sqrt(PEV)
  # EBV$PEV <-  PEV

  ##Estimate of reliabilities

  Reliability <- 1 - (PEV/var_u)
  Reliability <- ifelse(Reliability<0, NA, Reliability)

  EBV <- EBV |>
    dplyr::mutate(Standard_error = sqrt(PEV),
                  Prediction_error_variance = PEV,
                  Reliability =  Reliability)

  #EBV$Reliability <- Reliability

  output <- list(X_ebv, Coeff, EBV, PEV, Reliability)

  names(output) <- c('Posterior',
                     "Coefficient",
                     "Estimated_breeding_value",
                     "PEV",
                     "Reliability")

  return(output)

}
