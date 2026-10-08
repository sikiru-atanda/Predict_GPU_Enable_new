#' Calculate Heritability from ASReml Model with compound symmetry Model (CSM)
#'
#' This function calculates heritability from an ASReml model, particularly focusing on models that include cross-sectional data (CSM). It handles heterogeneity in residuals and can work with multiple genetic components. The function allows for the specification of heterogeneity groups and the inclusion of various kernels or genetic matrices.
#'
#' @param model An `asreml` model object resulting from fitting an ASReml model.
#' @param heter_groups A character string specifying the column in the model's data frame that represents different heterogeneity groups (environments, for example).
#' @param heter_resid A logical value indicating if heterogeneity in residuals should be considered.
#' @param names_in_inv_list A character vector of names identifying the inverse matrices in the model, typically related to different omics data types or genetic components.
#' @param inter_gen_pos An optional parameter to specify the position of interaction terms in genetic models.
#' @param gen_pos An optional parameter to specify the position of genetic terms in the model.
#' @param ... Additional arguments passed to the underlying functions.
#'
#' @return A list containing elements for heritability (`Heritability`), total genetic variance per omics data type (`varG_per_omics`), total genetic variance (`Total_genetic_var`), and residual variance (`Residual_Var`). Each element is tailored to the structure of the input model, especially considering the heterogeneity among groups.
#'
#' @details
#' The function extracts variance components from the specified ASReml model, computing heritability and genetic variances across specified heterogeneity groups. It is designed to work with complex genetic models that may include multiple omics layers or kernels. The function checks for stability in variance components and requires a stable model for accurate computation.
#'
#' @examples
#' # Assuming 'asreml_model' is a fitted ASReml model with proper variance structures:
#' results <- asreml_herit_CSM_new(model = asreml_model,
#'                                 heter_groups = "Environment",
#'                                 heter_resid = TRUE,
#'                                 names_in_inv_list = c("G_inv", "omic1_inv"),
#'                                 inter_gen_pos = 2,
#'                                 gen_pos = 1)
#' @export
asreml_herit_CSM_new <-  function(
    model = mod,
    heter_groups = NULL,
    heter_resid=FALSE,
    names_in_inv_list = NULL,
    inter_gen_pos = NULL,
    gen_pos = NULL,
    ...

){

  msg <- ""

  vc <- asreml_varcomp_table(model)
  ENV <- data.frame(model$mf)[, heter_groups]
  heter_grp = levels(ENV)
  n_heter_grp <- nlevels(ENV)
  VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE), ]
  #VAR_check_Pos = which(vc[, "bound"]=="F" | vc[, "bound"]=="U" |vc[, "bound"] =="?" |vc[, "bound"] =="S" )
  #VAR_check_Pos = which(VAR_check[, "bound"]=="F" | VAR_check[, "bound"]=="U" |VAR_check[, "bound"] =="?" |VAR_check[, "bound"] =="S" )

  VAR_check_Pos = which(VAR_check[, "bound"] =="?" |VAR_check[, "bound"] =="S")

  if(length(VAR_check_Pos)>1) {

    # stop(print(paste(paste("variance component for", as.character(heter_grp[VAR_check_Pos]), collapse = " and "),
    #                    "are unstable, refit the model")), call. = FALSE)

    stop(print(paste(msg, paste(paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = ","),
                                "are unstable, refit the model"))), call. = FALSE)
  } else {
    if(length(VAR_check_Pos)==1) {
      stop(print(paste(msg, paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
                       " is unstable, refix the model")), call. = FALSE)
    }
  }

  if (any(is.na(vc$std.error))) {
    stop(print(paste(msg, "Some variance component are non estimatable. Refix the model")), call. = FALSE)
  }

  VarE  <- vc[grep("!R", rownames(vc)), "component"]
  if (length(VarE) == 1L && n_heter_grp > 1L && !isTRUE(heter_resid)) {
    VarE <- rep(as.numeric(VarE), n_heter_grp)
  }
  if (length(VarE) != length(heter_grp)) {
    stop(print(paste(
      msg,
      "Residual variance components do not match the number of",
      heter_groups,
      "levels."
    )), call. = FALSE)
  }
  names(VarE) = heter_grp
  varG_per_omics = matrix(NA, nrow = length(names_in_inv_list), ncol = length(heter_grp))
  dimnames(varG_per_omics) <- list(names_in_inv_list, heter_grp)

  for (g in seq_along(names_in_inv_list)) {
    VarG = vc[grep(names_in_inv_list[g], rownames(vc)), drop = TRUE,"component"]
    VarG = rep(VarG, n_heter_grp)
    varG_per_omics[g, ] <-  VarG

  }

  #VarG = rep(VarG, n_heter_grp)

  H = matrix(NA, nrow = 1, ncol = length(heter_grp))
  colnames(H) = heter_grp
  Total_varG = matrix(NA, nrow = 1, ncol = length(heter_grp))
  colnames(Total_varG) = heter_grp

  # Rigor: per-environment heritability also reported with a delta-method SE.
  # In the compound-symmetry case the genetic numerator and the (per-env)
  # residual both map directly to variance-component rows, so vpredict gives an
  # exact h2 SE. We reuse the SAME greps that build the point estimate so the
  # vpredict formula matches it exactly; gp_asreml_heritability_se() falls back to
  # NA (never a wrong number) if the components cannot be matched or asreml errors.
  genetic_rows <- unique(unlist(lapply(names_in_inv_list, function(nm)
    rownames(vc)[grep(nm, rownames(vc))])))
  resid_rows <- rownames(vc)[grep("!R", rownames(vc))]
  H_SE = matrix(NA_real_, nrow = 1, ncol = length(heter_grp))
  colnames(H_SE) = heter_grp

  for (k in 1:ncol(H)) {
    varG = c()
    for (v in 1:length(names_in_inv_list)) {
      varG[[v]] = varG_per_omics[v, k]

    }

    Total_varG[, k] <- sum(unlist(varG))
    H[, k] <- sum(unlist(varG))/(sum(unlist(varG))+VarE[k])
    ve_row_k <- if (length(resid_rows) == 1L) resid_rows else resid_rows[k]
    H_SE[, k] <- gp_asreml_heritability_se(model, vg_rownames = genetic_rows,
                                           ve_rowname = ve_row_k)
  }


   return(list(Heritability = H,
                Heritability_SE = H_SE,
                Total_genetic_var = Total_varG,
                varG_per_omics = varG_per_omics,
                Residual_Var = VarE))

}


