#' Calculate and Format Variance Components for Bayesian Models
#'
#' This function computes variance components, including genomic heritability, from Bayesian genomic selection model outputs. It aggregates mean variance across omics data, calculates total genetic variance, residual variance, genomic heritability, and their associated standard errors.
#'
#' @param var_u_mean_omics_list Numeric vector or list containing mean genetic variances for each omics dataset.
#' @param se_var_u_omics_list Numeric vector or list containing standard errors of genetic variances for each omics dataset.
#' @param var_u_total Numeric, mean total genetic variance across all omics datasets.
#' @param var_residual Numeric, mean residual variance.
#' @param se_var_residual Numeric, standard error of the residual variance.
#'
#' @return A data frame with rows for each variance component (genetic variance for each omics dataset, total genetic variance, residual variance, and heritability) and columns for the component values and their standard errors.
#'
#' @details
#' The function is designed to work with outputs from Bayesian genomic selection models. It simplifies the process of extracting key variance components and calculating genomic heritability from such models. The function assumes that input variances and standard errors are pre-computed and provided as inputs. It uses these inputs to calculate mean genomic heritability and its standard error.
#'
#' @examples
#' var_u_mean_omics_list <- c(0.2, 0.3)
#' se_var_u_omics_list <- c(0.05, 0.06)
#' var_u_total <- 0.5
#' var_residual <- 0.4
#' se_var_residual <- 0.07
#'
#' variance_components <- bayes_variance_componentsnew(
#'   var_u_mean_omics_list = var_u_mean_omics_list,
#'   se_var_u_omics_list = se_var_u_omics_list,
#'   var_u_total = var_u_total,
#'   var_residual = var_residual,
#'   se_var_residual = se_var_residual
#' )
#' print(variance_components)
#' @export

bayes_variance_componentsnew <-  function(var_u_mean_omics_list,
                                       se_var_u_omics_list,
                                       var_u_total,
                                       var_residual,
                                       se_var_residual){

  genomic_h2 <- (var_u_total)/(var_u_total+var_residual)
  se_genomic_h2 <- standard_deviation(genomic_h2)
  genomic_h2 <- mean(genomic_h2)
  var_residual <- mean(var_residual)
  se_var_u_total = standard_deviation(var_u_total)
  if(length(var_u_mean_omics_list)>1){
variance_components <- data.frame(Components = c(unlist(var_u_mean_omics_list), mean(var_u_total), var_residual, genomic_h2),
                                  Standard_error = c(unlist(se_var_u_omics_list), se_var_u_total, se_var_residual, se_genomic_h2),
                                  row.names = c(paste("genetic_variance", seq_along(var_u_mean_omics_list), sep = "_"), "total_genetic_variance",
                                                "residual_variance", "heritability"),
                                  stringsAsFactors = FALSE
)
  }

  if(length(var_u_mean_omics_list)==1){
    variance_components <- data.frame(Components = c(unlist(var_u_mean_omics_list), var_residual, genomic_h2),
                                      Standard_error = c(unlist(se_var_u_omics_list), se_var_residual, se_genomic_h2),
                                      row.names = c("genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )
  }

return(variance_components)

}
