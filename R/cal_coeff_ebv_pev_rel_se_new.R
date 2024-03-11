#' Calculate Coefficients, EBVs, PEVs, and Reliability with Standard Errors
#'
#' This function computes coefficients, estimated breeding values (EBVs), prediction
#' error variances (PEVs), and reliability for genomic selection models, incorporating
#' standard errors. It is designed for use with both single-environment and multi-environment
#' datasets, accommodating models like RKHS and GBLUP.
#'
#' @param beta Numeric vector or matrix of regression coefficients or marker effects.
#' @param x_variable Numeric matrix of predictor variables or the genomic relationship matrix.
#' @param gen_name Character string specifying the column name in `x_variable` that
#'   represents the genetic identifiers for the individuals.
#' @param var_u Numeric, the genetic variance component estimated from the model.
#' @param hetero Optional; a vector indicating the heterogeneity groups of individuals,
#'   if applicable.
#' @param heter_groups Optional; a character string specifying the name of the column
#'   representing heterogeneity groups in the output, if `hetero` is not NULL.
#' @param gid_name Character string specifying the names of genetic identifiers in the
#'   output data frame.
#' @param GS_model Optional; character string specifying the genomic selection model used.
#'   This parameter can influence certain calculations or outputs.
#' @param ... Additional arguments passed to the function.
#'
#' @return A list containing:
#'   - `Posterior`: Matrix of posterior estimates derived from `beta` and `x_variable`.
#'   - `Coefficient`: A data frame of coefficients for each variable.
#'   - `Estimated_breeding_value`: A data frame of estimated breeding values for individuals.
#'   - `PEV`: Prediction error variance associated with the EBVs.
#'   - `Reliability`: Reliability of the EBVs, calculated as 1 - (PEV/var_u).
#'
#' @examples
#' # Assume beta, x_variable, and var_u are already defined:
#' results <- cal_coeff_ebv_pev_rel_se_new(beta = beta_coefficients,
#'                                         x_variable = genomic_matrix,
#'                                         gen_name = "GenID",
#'                                         var_u = estimated_genetic_variance,
#'                                         gid_name = "GenID")
#' @export
#' @importFrom dplyr mutate
#' @importFrom stats var
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
