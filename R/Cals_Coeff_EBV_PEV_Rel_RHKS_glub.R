
#' Title
#'
#' ### This function calculate the Coefficient, EBV, PEV, SE, Reliability
#'
#' @param gen_name
#' @param var_u genetic variance
#' @param ...
#' @param mod
#' @param gmatrix
#' @param var_E
#' @param hetero
#' @param heter_groups
#' @param gid_name
#' @importFrom magrittr |>
#' @return
#' @export
#'
#' @examples
Cal_Coeff_EBV_PEV_Rel_RHKS_glub <- function(mod = NULL,
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
                    EBV= g_ebv,
                    stringsAsFactors = FALSE)

  } else {
    EBV <- data.frame(names = rep(gid_name , length(unique(hetero))),
                      Env = rep(unique(hetero), each= length(gid_name)),
                      EBV= g_ebv,
                      stringsAsFactors = FALSE)

    names(EBV)[2] <- heter_groups

  }

  colnames(EBV)[1] <- gen_name

  ### Calculate the SEP, PEV and Reliability
  sep_pev_rel <- sep_pev_rel_gblup(geno_object = gmatrix,
                                    va = var_u,
                                    ve = var_E)

  EBV <- EBV |>
    dplyr::mutate(Std_error = sep_pev_rel$sep,
                  PEV = sep_pev_rel$pev,
                  Reliability = sep_pev_rel$rel)


  output <- list(coeffRaw, Coeff, EBV, sep_pev_rel$pev, sep_pev_rel$rel)

  names(output) <- c('Posterior',
                     "Coefficient",
                     "Estimated_breeding_value",
                     "PEV",
                     "Reliability")

  return(output)

}
