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
                                            eta_index = 1L,
                                            ...){

  eta_obj <- mod$model$ETA[[eta_index]]
  u_mean <- as.double(eta_obj$u)
  u_sd <- as.double(eta_obj$SD.u)
  pev <- u_sd^2

  if (is.null(hetero)) {
    Coeff <- data.frame(
      x_variables = gid_name,
      coeff = u_mean,
      stringsAsFactors = FALSE
    )
    EBV <- data.frame(
      names = gid_name,
      Estimated_breeding_value = u_mean,
      stringsAsFactors = FALSE
    )
  } else {
    if (length(gid_name) == length(u_mean) && length(hetero) == length(u_mean)) {
      # Record-level kernel (Z K Z'): one effect per phenotype record, in
      # record order. Label each record with its own genotype and group; a
      # genotype x group grid only matches this order for balanced data
      # sorted group-first in kernel order.
      coeff_ids <- as.character(gid_name)
      coeff_envs <- as.character(hetero)
    } else {
      coeff_ids <- rep(gid_name, length(unique(hetero)))
      coeff_envs <- rep(unique(hetero), each = length(gid_name))
    }
    Coeff <- data.frame(
      x_variables = coeff_ids,
      Env = coeff_envs,
      coeff = u_mean,
      stringsAsFactors = FALSE
    )
    names(Coeff)[2] <- heter_groups

    EBV <- data.frame(
      names = coeff_ids,
      Env = coeff_envs,
      Estimated_breeding_value = u_mean,
      stringsAsFactors = FALSE
    )
    names(EBV)[2] <- heter_groups
  }

  colnames(EBV)[1] <- gen_name

  Standard_error <- u_sd
  PEV <- pev
  Reliability <- 1 - (PEV / var_u)
  Reliability <- pmax(0, pmin(1, Reliability))
  EBV <- EBV |>
    dplyr::mutate(Standard_error = Standard_error,
                  Prediction_error_variance = PEV,
                  Reliability = Reliability)


  output <- list(Posterior = u_mean, Coefficient = Coeff,
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
