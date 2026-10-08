#' Extract an ASReml variance component table
#'
#' @param model_or_summary ASReml model object or \code{summary.asreml} object.
#'
#' @return A data frame of ASReml variance components.
#' @keywords internal
asreml_varcomp_table <- function(model_or_summary) {
  summary_obj <- if (inherits(model_or_summary, "summary.asreml")) {
    model_or_summary
  } else {
    asreml::summary.asreml(model_or_summary)
  }

  vc <- summary_obj[["varcomp"]]
  if (is.null(vc)) {
    stop("ASReml summary does not contain a varcomp table.", call. = FALSE)
  }

  vc <- as.data.frame(vc, stringsAsFactors = FALSE)

  if (is.null(rownames(vc)) || !length(rownames(vc)) || all(!nzchar(rownames(vc)))) {
    model_obj <- if (inherits(model_or_summary, "summary.asreml")) NULL else model_or_summary
    vp_names <- NULL
    if (!is.null(model_obj) && !is.null(model_obj[["vparameters.type"]])) {
      vp_names <- names(model_obj[["vparameters.type"]])
    }
    if (!is.null(vp_names) && length(vp_names) == nrow(vc)) {
      rownames(vc) <- vp_names
    }
  }

  vc
}

#' Extract and format variance components from ASReml results
#'
#' Processes ASReml variance summaries for single- or multi-environment models.
#'
#' @param res_var_cov_h_ve A list containing genetic variances per omic source,
#'   total genetic variance, residual variance, heritability, and optional
#'   covariance/correlation structures.
#'
#' @return A variance component data frame, or a list of environment-specific
#'   component tables and covariance/correlation matrices when available.
#' @export
asreml_variance_components <-  function(res_var_cov_h_ve
                                          ){


  varG_per_omics <- res_var_cov_h_ve[["varG_per_omics"]] ## matrix of rows omics and columns env
  total_varG <- res_var_cov_h_ve[["Total_genetic_var"]] ## matrix of rows omics and columns env
  var_residual <- as.matrix(res_var_cov_h_ve[["Residual_Var"]]) # matrix of rows env just single column
  hertiability <- res_var_cov_h_ve[["Heritability"]]  # # matrix of 1 row and columns of envs
  varG_per_omics_se <- res_var_cov_h_ve[["varG_per_omics_SE"]]
  if (is.null(varG_per_omics_se) ||
      !identical(dim(as.matrix(varG_per_omics_se)), dim(as.matrix(varG_per_omics)))) {
    varG_per_omics_se <- matrix(
      NA_real_, nrow = nrow(varG_per_omics), ncol = ncol(varG_per_omics),
      dimnames = dimnames(varG_per_omics)
    )
  } else {
    varG_per_omics_se <- as.matrix(varG_per_omics_se)
  }
  residual_se <- suppressWarnings(as.numeric(res_var_cov_h_ve[["Residual_Var_SE"]]))
  if (length(residual_se) != length(var_residual)) {
    residual_se <- rep(NA_real_, length(var_residual))
  }
  estimation_method <- as.character(
    res_var_cov_h_ve[["Estimation_method"]] %||% "ASReml_REML"
  )[[1L]]

  # Rigor: delta-method heritability SE (asreml::vpredict), threaded into the
  # heritability row's Standard_error. NA where it cannot be computed exactly
  # (e.g. factor-analytic) or for paths that do not supply it.
  heritability_se <- res_var_cov_h_ve[["Heritability_SE"]]
  get_h2_se <- function(j = 1L) {
    if (is.null(heritability_se)) return(NA_real_)
    hs <- as.matrix(heritability_se)
    if (j >= 1L && j <= ncol(hs)) as.numeric(hs[1L, j]) else NA_real_
  }

  ## One or more than one omics but single location/environment
  if(nrow(varG_per_omics)>=1 & ncol(varG_per_omics)==1){
    var_u_omics_list <- list()

    for (i in 1:nrow(varG_per_omics)) {
      var_u_omics_list[[i]] <- varG_per_omics[i, ]
    }

    if(nrow(varG_per_omics)>1){
    variance_components <- data.frame(Components = c(unlist(var_u_omics_list), total_varG, var_residual, hertiability),
                                      Standard_error = c(varG_per_omics_se[, 1L], NA_real_, residual_se[1L], get_h2_se(1L)),
                                      Estimation_method = estimation_method,
                                      row.names = c(paste("genetic_variance", seq_along(var_u_omics_list), sep = "_"), "total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )

    } else {

      if(nrow(varG_per_omics)==1){
        variance_components <- data.frame(Components = c(unlist(var_u_omics_list), var_residual, hertiability),
                                          Standard_error = c(varG_per_omics_se[1L, 1L], residual_se[1L], get_h2_se(1L)),
                                          Estimation_method = estimation_method,
                                          row.names = c("genetic_variance",
                                                        "residual_variance", "heritability"),
                                          stringsAsFactors = FALSE
        )

      }
    }

  } else {
    # Initialize variance_components and var_u_omics_list
    variance_components <- data.frame(Components = numeric(),
                                      Standard_error = numeric(),
                                      Estimation_method = character(),
                                      row.names = character(),
                                      stringsAsFactors = FALSE)

    ## One or more than one omics but multiple location/environment
    if (nrow(varG_per_omics) >= 1 & ncol(varG_per_omics) > 1) {
      for (j in 1:ncol(varG_per_omics)) { ## each column represent environment
        # Initialize var_u_omics_list for each column
        var_u_omics_list <- list()

        # Populate var_u_omics_list with values from varG_per_omics
        for (i in 1:nrow(varG_per_omics)) { ## each row represent each omics
          var_u_omics_list[[i]] <- varG_per_omics[i, j]
        }

        # Construct variance_components_j for each column
        if(nrow(varG_per_omics)>1){
        variance_components_j <- data.frame(Components = c(unlist(var_u_omics_list), total_varG[, j], var_residual[j,], hertiability[, j]),
                                            Standard_error = c(varG_per_omics_se[, j], NA_real_, residual_se[j], get_h2_se(j)),
                                            Estimation_method = estimation_method,
                                            row.names = c(paste(paste("genetic_variance", seq_along(var_u_omics_list), sep = "_"), colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("total_genetic_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("residual_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("heritability", colnames(varG_per_omics)[j], sep = "_")),
                                            stringsAsFactors = FALSE)

        } else {
          if(nrow(varG_per_omics)==1){
            variance_components_j <- data.frame(Components = c(unlist(var_u_omics_list), var_residual[j,], hertiability[, j]),
                                                Standard_error = c(varG_per_omics_se[1L, j], residual_se[j], get_h2_se(j)),
                                                Estimation_method = estimation_method,
                                                row.names = c(paste("genetic_variance",  colnames(varG_per_omics)[j], sep = "_"),
                                                              paste("residual_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                              paste("heritability", colnames(varG_per_omics)[j], sep = "_")),
                                                stringsAsFactors = FALSE)

          }

        }

        # Combine variance_components_j with variance_components
        variance_components <- rbind(variance_components, variance_components_j)
      }

      # Align ASReml MET variance_components with the Bayesian MET pattern
      # (bayes_multitrait_heter_resid.R::bayes_multitrait_env_heter_fit): prepend
      # three explicitly labelled across-environment summary rows ahead of the
      # per-environment rows, computed from the
      # SAME aggregates the Bayes path uses (mean of per-env Vg, mean of
      # per-env Ve, ratio of those means). This is cosmetic for downstream
      # extractors that key off per-env suffixes (heritability_<env> etc.),
      # but it makes the rownames pattern identical across kernel-Bayes
      # MET and ASReml MET so users inspecting variance_components see
      # the same layout regardless of engine.
      per_env_vg <- if (nrow(varG_per_omics) > 1L && !is.null(total_varG)) {
        as.numeric(total_varG[1L, , drop = TRUE])
      } else {
        as.numeric(varG_per_omics[1L, , drop = TRUE])
      }
      per_env_ve <- as.numeric(var_residual)
      mean_vg <- mean(per_env_vg, na.rm = TRUE)
      mean_ve <- mean(per_env_ve, na.rm = TRUE)
      overall_h2 <- if (is.finite(mean_vg + mean_ve) && (mean_vg + mean_ve) > 0) {
        mean_vg / (mean_vg + mean_ve)
      } else {
        NA_real_
      }
      summary_block <- data.frame(
        Components = c(mean_vg, mean_ve, overall_h2),
        Standard_error = NA_real_,
        Estimation_method = estimation_method,
        row.names = c("mean_genetic_variance_across_environments",
                      "mean_residual_variance_across_environments",
                      "heritability_from_mean_variances"),
        stringsAsFactors = FALSE
      )
      variance_components <- rbind(summary_block, variance_components)
    }

  }

  if (all(c("Covariance", "Correlation") %in% names(res_var_cov_h_ve))) {
    return(list(Variance_components = variance_components,
                Covariance = res_var_cov_h_ve[["Covariance"]],
                Correlation = res_var_cov_h_ve[["Correlation"]]))
  } else {
  return(variance_components)

  }

}
