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

  genetic_draws <- suppressWarnings(as.numeric(var_u_total))
  residual_draws <- suppressWarnings(as.numeric(var_residual))
  n_draws <- min(length(genetic_draws), length(residual_draws))
  if (n_draws > 0L) {
    genetic_draws <- utils::tail(genetic_draws, n_draws)
    residual_draws <- utils::tail(residual_draws, n_draws)
  }
  h2_denominator <- genetic_draws + residual_draws
  genomic_h2 <- ifelse(is.finite(h2_denominator) & h2_denominator > 0,
                       genetic_draws / h2_denominator,
                       NA_real_)
  se_genomic_h2 <- standard_deviation(genomic_h2)
  genomic_h2 <- mean(genomic_h2, na.rm = TRUE)
  var_residual <- mean(residual_draws, na.rm = TRUE)
  se_var_u_total = standard_deviation(genetic_draws)
  if(length(var_u_mean_omics_list)==0){
    variance_components <- data.frame(Components = c(mean(genetic_draws, na.rm = TRUE), var_residual, genomic_h2),
                                      Standard_error = c(se_var_u_total, se_var_residual, se_genomic_h2),
                                      row.names = c("total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )
  }
  if(length(var_u_mean_omics_list)>1){
component_names <- names(var_u_mean_omics_list)
if (is.null(component_names) || any(!nzchar(component_names))) {
  component_names <- as.character(seq_along(var_u_mean_omics_list))
}
component_names <- paste0("genetic_variance_", make.names(component_names))
variance_components <- data.frame(Components = c(unlist(var_u_mean_omics_list), mean(genetic_draws, na.rm = TRUE), var_residual, genomic_h2),
                                  Standard_error = c(unlist(se_var_u_omics_list), se_var_u_total, se_var_residual, se_genomic_h2),
                                  row.names = c(component_names, "total_genetic_variance",
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

#' Append per-environment Vg/Ve/h2 rows to a Bayesian MET variance-component table
#'
#' Internal. Brings the BGLR kernel-MET (`RKHS`, `GBLUP_BRR`) variance_components
#' into the same per-environment layout the ASReml and GP engines emit --
#' `genetic_variance_<env>`, `residual_variance_<env>`, `heritability_<env>` --
#' so MET variance output is consistent across model families. The legacy table
#' only carried per-ETA-term `genetic_variance_<i>` rows plus a single pooled
#' `heritability`; the index `i` counts random ETA terms (e.g. 2 kernels x
#' {GID, GID:env} = 4), NOT environments, which is opaque for a 5-environment
#' MET run.
#'
#' Per-env genetic variance is the mean marginal genetic variance from the
#' fitted covariance model within each environment. For kernel models this is
#' the corresponding diagonal of the total genetic covariance matrix; for BRR
#' it is the diagonal of the implied marker-effect covariance. This avoids the
#' downward shrinkage bias from taking the empirical variance of posterior mean
#' breeding values.
#' Per-env residual variance is the BGLR residual expanded to the observations:
#' a single `varE` -> constant across environments; a heterogeneous `varE` ->
#' the per-group residual. Heritability is `Vg_env / (Vg_env + Ve_env)`,
#' computed by [gp_variance_components_add_per_env].
#'
#' @param variance_components A Bayesian variance-component data frame (with
#'   `Components` and `Standard_error` columns).
#' @param env_labels Character/factor vector of per-observation environment
#'   labels (length = number of observations).
#' @param marginal_genetic_variance_per_obs Numeric per-observation marginal
#'   genetic variances from the fitted model covariance.
#' @param varE BGLR residual variance (scalar, or named per environment).
#' @param weights Optional per-observation weights (defaults to 1).
#' @return `variance_components` with per-environment rows appended.
#' @keywords internal
#' @noRd
bayes_variance_components_add_per_env <- function(variance_components,
                                                  env_labels,
                                                  marginal_genetic_variance_per_obs,
                                                  varE,
                                                  weights = NULL) {
  if (!is.data.frame(variance_components)) {
    return(variance_components)
  }
  env_labels <- as.character(env_labels)
  n <- length(env_labels)
  gvar <- suppressWarnings(as.numeric(marginal_genetic_variance_per_obs))
  if (length(gvar) != n || !any(is.finite(gvar))) {
    return(variance_components)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  ve_obs <- tryCatch(
    gp_bayes_observation_residual_variance(
      varE = varE,
      groups = env_labels,
      weights = weights,
      n = n
    ),
    error = function(e) {
      rep(mean(suppressWarnings(as.numeric(varE)), na.rm = TRUE), n)
    }
  )
  ve_obs <- suppressWarnings(as.numeric(ve_obs))

  # Average the observation-level marginal variances within each environment.
  # This remains defined for a one-observation environment and retains kernel
  # diagonal scaling when kernels have not been normalized to unit diagonal.
  env_genetic <- tapply(gvar, env_labels, function(z) {
    z <- z[is.finite(z)]
    if (length(z)) mean(z) else NA_real_
  })
  env_residual <- tapply(ve_obs, env_labels, function(z) mean(z, na.rm = TRUE))
  envs <- names(env_genetic)
  env_residual <- env_residual[envs]

  gp_variance_components_add_per_env(
    vc = variance_components,
    env_genetic  = stats::setNames(as.numeric(env_genetic),  envs),
    env_residual = stats::setNames(as.numeric(env_residual), envs)
  )
}
