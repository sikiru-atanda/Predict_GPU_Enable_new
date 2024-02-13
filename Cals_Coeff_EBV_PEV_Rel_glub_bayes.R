
#' Title
#'
#' ### This function calculate the Coefficient, EBV, PEV, SE, Reliability
#'
#' @param beta  the posterior samples of the random effects (marker/omic)
#' @param x_matrix marker/omic data
#' @param gen_name
#' @param var_u genetic variance
#'
#' @param mod
#' @param gmatrix
#' @param var_E
#' @param hetero
#' @param heter_groups
#' @param GS_model
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
Cal_Coeff_EBV_PEV_Rel_glub_bayes <- function(mod = NULL,
                                       beta = NULL,
                                       x_matrix = NULL,
                                       gmatrix = NULL,
                                       gen_name = NULL,
                                       var_u = NULL,
                                       var_E = NULL,
                                       hetero = NULL,
                                       heter_groups = NULL,
                                       GS_model = NULL,
                                       ...){
  if(GS_model == "RKHS"){
  g_ebv= mod$model$yHat - mod$model$mu
  # the breeding values (posterior)
  # Posterior means of marker effects/coefficient approximation
  coeffRaw = solve(t(gmatrix)*gmatrix)*t(gmatrix)*gmatrix

  # Posterior means of marker effects/coefficient

  Coeff <- data.frame(x_matrixs = colnames(gmatrix),
                      coeff = colMeans(coeffRaw),
                      stringsAsFactors = FALSE)


  if(is.null(hetero)){
    # Genomic estimated breeding values
    EBV <- data.frame(names = rownames(gmatrix),
                      EBV= g_ebv,
                      stringsAsFactors = FALSE)

  } else {
    EBV <- data.frame(names = rep(rownames(gmatrix), length(hetero)),
                      Env = rep(hetero, each= nrow(gmatrix)),
                      EBV= g_ebv,
                      stringsAsFactors = FALSE)

    names(EBV)[2] <- heter_groups

  }

  colnames(EBV)[1] = gen_name

  ### Calculate the SEP, PEV and Reliability
  sep_pev_rel <- sep_pev_rel_gblup(geno_object = gmatrix,
                                   va = var_u,
                                   ve = var_E)

  EBV <- EBV |>
    dplyr::mutate(Std_error = sep_pev_rel$sep,
                  PEV = sep_pev_rel$pev,
                  Reliability = sep_pev_rel$rel)


  output <- list(coeffRaw, Coeff, EBV, sep_pev_rel$pev, sep_pev_rel$rel)

  } else {
    if(GS_model=="BRR"){

      # the breeding values (posterior)
      X_ebv <- as.matrix(x_matrix)%*%t(beta)

      # Posterior means of marker effects/coefficient
      #coeff_omic3 <- colMeans(Bb)

      Coeff <- data.frame(x_matrixs = colnames(x_matrix),
                          coeff =colMeans(beta),
                          stringsAsFactors = FALSE)

      if(is.null(hetero)){
        # Genomic estimated breeding values
        EBV <- data.frame(names = rownames(x_matrix),
                          EBV=x_matrix%*%Coeff$coeff,
                          stringsAsFactors = FALSE)

      } else {

        EBV <- data.frame(names = rep(rownames(x_matrix), length(hetero)),
                          Env = rep(hetero, each= nrow(x_matrix)),
                          EBV=x_matrix%*%Coeff$coeff,
                          stringsAsFactors = FALSE)

        names(EBV)[2] <- heter_groups



      }

      colnames(EBV)[1] = gen_name
      #GEBV = data.frame(rowMeans(ebv))

      PEV <- apply(X_ebv, 1, var)

      Reliability <- 1 - (PEV/var_u)
      Reliability <- ifelse(Reliability<0, NA, Reliability)
      #####
      EBV <- EBV |>
        dplyr::mutate(Std_error = sqrt(PEV),
                      PEV = PEV,
                      Reliability = Reliability)


      output <- list(X_ebv, Coeff, EBV, PEV, Reliability)

    }

  }

  names(output) <- c('Posterior',
                     "Coefficient",
                     "Estimated_breeding_value",
                     "PEV",
                     "Reliability")

  return(output)

}
