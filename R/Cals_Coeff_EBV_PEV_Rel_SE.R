
#' Title
#'
#' ### This function calculate the Coefficient, EBV, PEV, SE, Reliability
#' @param beta  the posterior samples of the random effects (marker/omic)
#' @param x_variable marker/omic data
#' @param gen_name
#' @param var_u genetic variance
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
Cal_Coeff_EBV_PEV_Rel_SE <- function(beta,
                                     x_variable,
                                     gen_name,
                                     var_u,
                                     ...){

        # the breeding values (posterior)
        X_ebv <- as.matrix(x_variable)%*%t(beta)

        # Posterior means of marker effects/coefficient
        #coeff_omic3 <- colMeans(Bb)

        Coeff <- data.frame(x_variables = colnames(x_variable),
                            coeff =colMeans(beta),
                            stringsAsFactors = FALSE)

        # Genomic estimated breeding values
        EBV <- data.frame(names = rownames(x_variable),
                          EBV=x_variable%*%Coeff$coeff,
                          stringsAsFactors = FALSE)

        colnames(EBV)[1] = gen_name
        #GEBV = data.frame(rowMeans(ebv))

        PEV <- apply(X_ebv, 1, var)

        EBV$Std_error <- sqrt(PEV)
        EBV$PEV <-  PEV

        ##Estimate of reliabilities

        Reliability <- 1 - (PEV/var_u)
        EBV$Reliability <- Reliability

        output <- list(X_ebv, Coeff, EBV, PEV, Reliability)

        names(output) <- c('Posterior',
                           "Coefficient",
                           "Estimated_breeding_value",
                           "PEV",
                           "Reliability")

        return(output)

      }
