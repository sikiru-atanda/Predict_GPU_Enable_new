#' Title
#'
#' @param res_var_cov_h_ve
#'
#' @return
#' @export
#'
#' @examples
asreml_variance_components <-  function(res_var_cov_h_ve
                                          ){


  varG_per_omics <- res_var_cov_h_ve[["varG_per_omics"]] ## matrix of rows omics and columns env
  total_varG <- res_var_cov_h_ve[["Total_genetic_var"]] ## matrix of rows omics and columns env
  var_residual <- as.matrix(res_var_cov_h_ve[["Residual_Var"]]) # matrix of rows env just single column
  hertiability <- res_var_cov_h_ve[["Heritability"]]  # # matrix of 1 row and columns of envs

  ## One or more than one omics but single location/environment
  if(nrow(varG_per_omics)>=1 & ncol(varG_per_omics)==1){
    var_u_omics_list <- list()

    for (i in 1:nrow(varG_per_omics)) {
      var_u_omics_list[[i]] <- varG_per_omics[i, ]
    }

    if(nrow(varG_per_omics)>1){
    variance_components <- data.frame(Components = c(unlist(var_u_omics_list), total_varG, var_residual, hertiability),
                                      Standard_error = NA,
                                      row.names = c(paste("genetic_variance", seq_along(var_u_omics_list), sep = "_"), "total_genetic_variance",
                                                    "residual_variance", "heritability"),
                                      stringsAsFactors = FALSE
    )

    } else {

      if(nrow(varG_per_omics)==1){
        variance_components <- data.frame(Components = c(unlist(var_u_omics_list), var_residual, hertiability),
                                          Standard_error = NA,
                                          row.names = c("genetic_variance",
                                                        "residual_variance", "heritability"),
                                          stringsAsFactors = FALSE
        )

      }
    }

  } else {
    # Initialize variance_components and var_u_omics_list
    variance_components <- data.frame(Components = character(),
                                      Standard_error = numeric(),
                                      row.names = character(),
                                      stringsAsFactors = FALSE)

    ## One or more than one omics but multiple location/environment
    if (nrow(varG_per_omics) >= 1 & ncol(varG_per_omics) > 1) {
      for (j in 1:ncol(varG_per_omics)) {
        # Initialize var_u_omics_list for each column
        var_u_omics_list <- list()

        # Populate var_u_omics_list with values from varG_per_omics
        for (i in 1:nrow(varG_per_omics)) {
          var_u_omics_list[[i]] <- varG_per_omics[i, j]
        }

        # Construct variance_components_j for each column
        if(nrow(varG_per_omics)>1){
        variance_components_j <- data.frame(Components = c(unlist(var_u_omics_list), total_varG[, j], var_residual[j,], hertiability[, j]),
                                            Standard_error = NA,
                                            row.names = c(paste(paste("genetic_variance", seq_along(var_u_omics_list), sep = "_"), colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("total_genetic_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("residual_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                          paste("heritability", colnames(varG_per_omics)[j], sep = "_")),
                                            stringsAsFactors = FALSE)

        } else {
          if(nrow(varG_per_omics)==1){
            variance_components_j <- data.frame(Components = c(unlist(var_u_omics_list), var_residual[j,], hertiability[, j]),
                                                Standard_error = NA,
                                                row.names = c(paste("genetic_variance",  colnames(varG_per_omics)[j], sep = "_"),
                                                              paste("residual_variance", colnames(varG_per_omics)[j], sep = "_"),
                                                              paste("heritability", colnames(varG_per_omics)[j], sep = "_")),
                                                stringsAsFactors = FALSE)

          }

        }

        # Combine variance_components_j with variance_components
        variance_components <- rbind(variance_components, variance_components_j)
      }
    }

  }

  if (all(c("Covariance", "Correlation") %in% names(res_var_cov_h_ve))) {
    if(nrow(varG_per_omics) == 1){
      names(res_var_cov_h_ve[[grep("Covariance", names(res_var_cov_h_ve), ignore.case = TRUE)]]) <- "Covariance"
      names(res_var_cov_h_ve[[grep("Correlation", names(res_var_cov_h_ve), ignore.case = TRUE)]]) <- "Correlation"

    }

    return(list(Variance_components = variance_components,
                Covariance = res_var_cov_h_ve[[grep("Covariance", names(res_var_cov_h_ve), ignore.case = TRUE)]],
                Correlation = res_var_cov_h_ve[[grep("Correlation", names(res_var_cov_h_ve), ignore.case = TRUE)]]))
  } else {
  return(variance_components)

  }

}
