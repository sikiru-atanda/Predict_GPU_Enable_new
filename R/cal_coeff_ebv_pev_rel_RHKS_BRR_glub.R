#' Calculate Coefficients, EBV, PEV, and Reliability for RKHS and GBLUP Models
#'
#' This function computes the posterior means of marker effects (coefficients),
#' estimated breeding values (EBV), prediction error variance (PEV), and reliability
#' for genomic selection models based on RKHS and GBLUP approaches. It uses a genomic
#' relationship matrix, model outputs, and variance components to compute these values.
#'
#' @param mod A model object returned from a genomic selection analysis, which
#'   contains the fitted model components such as effects and predictions.
#' @param gmatrix A numeric matrix representing the genomic relationship matrix (GRM) or other kernels from omics
#'   used in the model. Row and column names should correspond to individual IDs.
#' @param gen_name A character string specifying the column name in `gmatrix` that
#'   represents the genetic identifiers for the individuals.
#' @param var_u The genetic variance component estimated from the model.
#' @param var_E The residual variance component estimated from the model.
#' @param hetero A vector indicating the heterogeneity groups of individuals, if applicable.
#' @param heter_groups A character string specifying the name of the column representing
#'   heterogeneity groups in the output, if `hetero` is not NULL.
#' @param gid_name A character string specifying the names of genetic identifiers in the
#'   output data frame.
#' @param ... Additional arguments passed to the function.
#'
#' @return A list containing:
#'   - `Posterior`: Raw posterior means of marker effects.
#'   - `Coefficient`: A data frame of coefficients for each variable.
#'   - `Estimated_breeding_value`: A data frame of estimated breeding values for individuals.
#'   - `PEV`: Prediction error variance associated with the EBVs.
#'   - `Reliability`: Reliability of the EBVs.
#'
#' @examples
#' # Assuming 'mod' is a model object from a genomic selection analysis,
#' # 'gmatrix' is a genomic relationship matrix, and 'var_u' and 'var_E'
#' # are estimated variance components:
#' result <- cal_coeff_ebv_pev_rel_RHKS_glub(mod = model_obj,
#'                                           gmatrix = genomic_matrix,
#'                                           gen_name = "IndividualID",
#'                                           var_u = 0.5,
#'                                           var_E = 0.3,
#'                                           gid_name = "GenID")
#' @export
#' @importFrom dplyr mutate
#' @importFrom stats solve t
cal_coeff_ebv_pev_rel_RHKS_glub <- function(mod = NULL,
                                            gmatrix = NULL,
                                            gen_name = NULL,
                                            var_u = NULL,
                                            var_E = NULL,
                                            hetero = NULL,
                                            heter_groups = NULL,
                                            gid_name = NULL,
                                            ...){

  g_ebv <- mod$model$yHat - mod$model$mu
  # the breeding values (posterior)
  # Posterior means of marker effects/coefficient approximation
  coeffRaw = solve(t(gmatrix)*gmatrix)*t(gmatrix)*gmatrix

  # Posterior means of marker effects/coefficient

  if(!is.null(heter_groups)){
    Coeff <- data.frame(x_variables = rep(gid_name , length(unique(hetero))),
                        Env = rep(unique(hetero), each= length(gid_name)),
                        coeff = colMeans(coeffRaw),
                        stringsAsFactors = FALSE)
    names(Coeff)[2] <- heter_groups
  } else {
    Coeff <- data.frame(x_variables = gid_name,
                        coeff = colMeans(coeffRaw),
                        stringsAsFactors = FALSE)

  }


  if(is.null(hetero)){
    # Genomic estimated breeding values
    EBV <- data.frame(names = rownames(gmatrix),
                      Estimated_breeding_value= g_ebv,
                      stringsAsFactors = FALSE)

  } else {
    EBV <- data.frame(names = rep(gid_name , length(unique(hetero))),
                      Env = rep(unique(hetero), each= length(gid_name)),
                      Estimated_breeding_value= g_ebv,
                      stringsAsFactors = FALSE)

    names(EBV)[2] <- heter_groups

  }

  colnames(EBV)[1] <- gen_name

  ### Calculate the SEP, PEV and Reliability
  # sep_pev_rel <- sep_pev_rel_gblup(geno_object = gmatrix,
  #                                  va = var_u,
  #                                  ve = var_E)

  # EBV <- EBV |>
  #   dplyr::mutate(Standard_error = sep_pev_rel$sep,
  #                 Prediction_error_variance = sep_pev_rel$pev,
  #                 Reliability = sep_pev_rel$rel)
  Standard_error = mod$model$SD.yHat
  PEV <- (mod$model$SD.yHat)^2

  Reliability <- 1 - (PEV/var(mod$model$yHat))
  EBV <- EBV |>
    dplyr::mutate(Standard_error = Standard_error,
                  Prediction_error_variance = PEV,
                  Reliability = Reliability)


  output <- list(Posterior = coeffRaw, Coefficient = Coeff,
                 Estimated_breeding_value = EBV, PEV = PEV,
                 Reliability = Reliability,
                 Standard_error = Standard_error)

  # names(output) <- c('Posterior',
  #                    "Coefficient",
  #                    "Estimated_breeding_value",
  #                    "PEV",
  #                    "Reliability",
  #                    "Standard_error")

  return(output)

}
