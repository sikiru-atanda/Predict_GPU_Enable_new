

#' Build a nested per-omics / per-environment list scaffold
#'
#' Internal helper: allocates a nested list of empty vectors with outer length =
#' number of omics layers (`names_in_inv_list`) and inner length = number of
#' heterogeneous groups (`het_gp`), used to accumulate per-omics, per-env
#' variance / heritability components during ASReml MET output assembly.
#'
#' @param names_in_inv_list Character vector of inverse-relationship object
#'   names (one per omics layer).
#' @param het_gp Character vector of heterogeneous-group (environment) labels.
#'
#' @return A list of length `length(names_in_inv_list)`, each element a list
#'   of length `length(het_gp)`.
#' @export
#'
#' @examples
List_list = function(names_in_inv_list, het_gp){
  All_varGs = vector(mode = "list", length = length(names_in_inv_list))
  for (l in 1:length(names_in_inv_list)) {
    All_varGs[[l]] <- vector(mode ="list", length(het_gp))

  }

  return(All_varGs)
}

#' Resolve response-scale ASReml residual variances
#'
#' Internal helper for the Stage 2 precision-weight contract. ASReml reports a
#' residual scale parameter, while a weighted Gaussian fit has observation-level
#' variance `scale / precision_weight`. When exact weights are supplied with
#' `dispersion = 1`, the reported scale is fixed at one and must not itself be
#' labelled as the response-scale residual variance.
#'
#' @keywords internal
#' @noRd
gp_asreml_residual_variance_contract <- function(model,
                                                 heter_groups,
                                                 heter_grp,
                                                 residual_scale,
                                                 residual_scale_se,
                                                 response = NULL,
                                                 gen_name = NULL) {
  mf <- as.data.frame(model$mf, stringsAsFactors = FALSE)
  env <- as.character(mf[[heter_groups]])
  residual_scale <- stats::setNames(as.numeric(residual_scale), heter_grp)
  residual_scale_se <- stats::setNames(as.numeric(residual_scale_se), heter_grp)

  analysis_row <- !is.na(env) & env %in% heter_grp
  if (!is.null(response) && response %in% names(mf)) {
    observed <- !is.na(mf[[response]])
    if (is.numeric(mf[[response]])) {
      observed <- observed & is.finite(mf[[response]])
    }
    analysis_row <- analysis_row & observed
  }
  if (any(!heter_grp %in% env[analysis_row])) {
    missing_env <- heter_grp[!heter_grp %in% env[analysis_row]]
    stop(
      "ASReml residual-variance summaries have no analyzed observations for environment(s): ",
      paste(missing_env, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  weight_column <- gp_stage2_weight_column()
  weighted <- weight_column %in% names(mf)
  if (!weighted) {
    n_by_env <- vapply(
      heter_grp,
      function(label) sum(analysis_row & env == label),
      integer(1L)
    )
    environment_summary <- data.frame(
      Env = heter_grp,
      N_observations = unname(n_by_env),
      Residual_scale_parameter = unname(residual_scale),
      Mean_precision_weight = NA_real_,
      Residual_variance = unname(residual_scale),
      Minimum_residual_variance = unname(residual_scale),
      Median_residual_variance = unname(residual_scale),
      Maximum_residual_variance = unname(residual_scale),
      Residual_variance_basis = "ASReml_fitted_residual_component",
      stringsAsFactors = FALSE
    )
    return(list(
      residual_variance = residual_scale,
      residual_variance_se = residual_scale_se,
      residual_scale = residual_scale,
      residual_scale_se = residual_scale_se,
      residual_variance_by_environment = environment_summary,
      residual_variance_by_observation = NULL,
      residual_variance_basis = "ASReml_fitted_residual_component",
      weighted = FALSE
    ))
  }

  precision <- suppressWarnings(as.numeric(mf[[weight_column]]))
  invalid_weight <- analysis_row & (!is.finite(precision) | precision <= 0)
  if (any(invalid_weight)) {
    stop(
      "ASReml fitted model contains a non-positive or non-finite Stage 2 precision weight on an analyzed row.",
      call. = FALSE
    )
  }
  scale_for_row <- unname(residual_scale[match(env, heter_grp)])
  residual_variance_for_row <- scale_for_row / precision
  basis <- "ASReml_residual_scale_divided_by_stage2_precision_weight"

  residual_variance <- stats::setNames(
    vapply(
      heter_grp,
      function(label) mean(residual_variance_for_row[analysis_row & env == label]),
      numeric(1L)
    ),
    heter_grp
  )
  mean_inverse_weight <- stats::setNames(
    vapply(
      heter_grp,
      function(label) mean(1 / precision[analysis_row & env == label]),
      numeric(1L)
    ),
    heter_grp
  )
  residual_variance_se <- residual_scale_se * mean_inverse_weight

  analyzed_index <- which(analysis_row)
  observation_detail <- data.frame(
    Observation_index = analyzed_index,
    stringsAsFactors = FALSE
  )
  if (!is.null(gen_name) && gen_name %in% names(mf)) {
    observation_detail[["GID"]] <- as.character(mf[[gen_name]][analyzed_index])
  }
  observation_detail[["Env"]] <- env[analyzed_index]
  observation_detail[["Precision_weight"]] <- precision[analyzed_index]
  observation_detail[["Residual_scale_parameter"]] <- scale_for_row[analyzed_index]
  observation_detail[["Residual_variance"]] <- residual_variance_for_row[analyzed_index]
  observation_detail[["Residual_variance_basis"]] <- basis

  environment_summary <- do.call(
    rbind,
    lapply(heter_grp, function(label) {
      idx <- analysis_row & env == label
      values <- residual_variance_for_row[idx]
      data.frame(
        Env = label,
        N_observations = sum(idx),
        Residual_scale_parameter = unname(residual_scale[[label]]),
        Mean_precision_weight = mean(precision[idx]),
        Residual_variance = mean(values),
        Minimum_residual_variance = min(values),
        Median_residual_variance = stats::median(values),
        Maximum_residual_variance = max(values),
        Residual_variance_basis = basis,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(environment_summary) <- NULL

  list(
    residual_variance = residual_variance,
    residual_variance_se = residual_variance_se,
    residual_scale = residual_scale,
    residual_scale_se = residual_scale_se,
    residual_variance_by_environment = environment_summary,
    residual_variance_by_observation = observation_detail,
    residual_variance_basis = basis,
    weighted = TRUE
  )
}

#' Calculate Heritability and extract Variance component, Covariance Structures from ASReml Model
#'
#' This function calculates heritability and extracts variance-covariance and correlation structures from an ASReml model. It is designed to work with models that may include heterogeneity in residuals and supports the inclusion of multiple genetic or omics components. The function also handles different variance-covariance structures specified by the user.
#'
#' @param model An `asreml` model object resulting from fitting an ASReml model.
#' @param heter_groups A character string specifying the column in the model's data frame that represents different heterogeneity groups (e.g., environments).
#' @param var_cov_str A character string indicating the type of variance-covariance structure to be analyzed (`"us"`, `"fa"`, `"rr"`, `"corgh"`, `"corh"`, or `"corv"`).
#' @param heter_resid A logical indicating whether heterogeneity in residuals is considered.
#' @param names_in_inv_list A character vector of names identifying the inverse matrices in the model, related to different genetic or omics components.
#' @param inter_gen_pos Optional parameter specifying the position of interaction terms in genetic models.
#' @param gen_pos Optional parameter specifying the position of genetic terms in the model.
#' @param response Optional response-column name. Used to exclude prediction-only
#'   rows when summarizing weighted residual variances.
#' @param gen_name Optional genotype identifier column used in the detailed
#'   observation-level residual-variance table.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{Covariance}{A list of variance-covariance matrices for each genetic component.}
#'   \item{Correlation}{A list of correlation matrices for each genetic component.}
#'   \item{Heritability}{A matrix of heritability estimates for each heterogeneity group.}
#'   \item{Total_genetic_var}{Total genetic variance for each heterogeneity group.}
#'   \item{varG_per_omics}{Variance attributed to each genetic component for each heterogeneity group.}
#'   \item{Residual_Var}{Residual variance for each heterogeneity group.}
#' }
#'
#' @details
#' The function primarily focuses on extracting and computing key genetic statistics from a given ASReml model. It supports a range of variance-covariance structures, including but not limited to unstructured (`us`), factor analytic (`fa`), and various correlation structures (`corgh`, `corh`, `corv`). The function checks for stability in variance components and requires a stable model for accurate computations.
#'
#' @examples
#' # Assuming 'asreml_model' is a fitted ASReml model with proper variance structures:
#' results <- asreml_herit_varCov_new(model = asreml_model,
#'                                    heter_groups = "Environment",
#'                                    var_cov_str = "fa2",
#'                                    heter_resid = TRUE,
#'                                    names_in_inv_list = c("G_inv", "omic1_inv"),
#'                                    inter_gen_pos = 2,
#'                                    gen_pos = 1)
#' @export

asreml_herit_varCov_new <-  function(
    model = NULL,
    heter_groups = NULL,
    var_cov_str= NULL,
    heter_resid=FALSE,
    names_in_inv_list = NULL,
    inter_gen_pos = NULL,
    gen_pos = NULL,
    response = NULL,
    gen_name = NULL,
    ...

){

  gp_reject_obsolete_asreml_structure(var_cov_str)
  msg <- ""

  if (!inherits(model, "summary.asreml")) {
    gp_asreml_assert_trustworthy_fit(
      model,
      context = paste0("ASReml MET GBLUP (", as.character(var_cov_str)[1L], ")")
    )
  }

  vc <- asreml_varcomp_table(model)

  ENV <- data.frame(model$mf)[, heter_groups]
  if (!is.factor(ENV)) ENV <- factor(ENV, levels = unique(as.character(ENV)))
  heter_grp <- levels(ENV)
  n_heter_grp <- length(heter_grp)
  if (n_heter_grp < 2L) {
    stop("ASReml MET covariance extraction requires at least two environments.",
         call. = FALSE)
  }

  VarCov <- matrix(NA, ncol = n_heter_grp, nrow = n_heter_grp)
  CORR <- matrix(NA, ncol = n_heter_grp, nrow = n_heter_grp)
  corr_all = vector(mode = 'list', length = length(names_in_inv_list))
  varcov_all = vector(mode = 'list', length = length(names_in_inv_list))
  ########
  extracted_names_from_inv_list <- gsub("_inv", "", names_in_inv_list, ignore.case = TRUE)
  index_names_inv_extracted_omic <- grep("omic", extracted_names_from_inv_list, ignore.case = TRUE)
  index_names_inv_extracted_geno <- grep("gmatrix", extracted_names_from_inv_list, ignore.case = TRUE)

  if (length(index_names_inv_extracted_omic) >1) {
    extracted_names_from_inv_list[index_names_inv_extracted_omic] <- paste("omics", seq_along(index_names_inv_extracted_omic))
  } else {
    if (length(index_names_inv_extracted_omic) == 1) {
      extracted_names_from_inv_list[index_names_inv_extracted_omic] <- "omics"
    }
  }

  if (length(index_names_inv_extracted_geno)==1) {
    extracted_names_from_inv_list[index_names_inv_extracted_geno] <- "geno"
  }
  ##########

  if(length(names_in_inv_list)>1){
  names(varcov_all) <- paste("Genetic_covariance_environments", extracted_names_from_inv_list, sep = "_")
  names(corr_all) <- paste("Genetic_correlation_environments", extracted_names_from_inv_list, sep = "_")
  }
################
  if(length(names_in_inv_list)==1){
    names(varcov_all) <- "Genetic_covariance_environments_by_kernel"
    names(corr_all) <- "Genetic_correlation_environments_by_kernel"

  }
#####################


  for (ca in seq_along(corr_all)) {
    corr_all[[ca]] <-  CORR

    varcov_all[[ca]] <- VarCov

  }

  # B - fixed at a boundary (!GP)
  # ? - liable to change from P to B
  # C - Constrained by user (!VCC) U - unbounded
  # S - Singular Information matrix # S means there is no information in the data for this parameter.
  # F - fixed by user
  # P - positive definite
  # U - unbounded
  ############################################
  #VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE), ]
  VAR_check <- vc[apply(sapply(heter_grp, function(env) grepl(paste0("!", env, "!"), rownames(vc))), 1, any), ]


  #VAR_check_Pos = which(vc[, "bound"]=="F" | vc[, "bound"]=="U" |vc[, "bound"] =="?" |vc[, "bound"] =="S" )
  #VAR_check_Pos = which(VAR_check[, "bound"]=="F" | VAR_check[, "bound"]=="U" |VAR_check[, "bound"] =="?" |VAR_check[, "bound"] =="S" )

  VAR_check_Pos = which(VAR_check[, "bound"] =="?" |VAR_check[, "bound"] =="S" )


  if(length(VAR_check_Pos)>1) {

    # stop(print(paste(paste("variance component for", as.character(heter_grp[VAR_check_Pos]), collapse = " and "),
    #                    "are unstable, refit the model")), call. = FALSE)

    stop(print(paste(msg, paste(paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = ","),
                                "are unstable, refit the model"))), call. = FALSE)
  } else {

    if(length(VAR_check_Pos)==1) {

      stop(print(paste(msg,paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
                       " is unstable, refix the model")), call. = FALSE)

    }

  }

  #### Extract variance and covariance for For factor analytic models
  if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
    ## Check for all variable is positive definitive
    ## Check for this other random term can be present aside the genetic effect
    #VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE),"bound"]

    CheckR <- vc[grep("!R", rownames(vc)), "component"]
    VarG_all = vector("list", length = length(names_in_inv_list))

    ## Extract number of factor(s) specified by users
    N_fa = as.double(substr(var_cov_str, 3, 100))
    Fac = seq(1, N_fa)


    for(l in seq_along(names_in_inv_list)){

      LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
      VarG <- LL[grep("!var", rownames(LL)), drop=FALSE, ]
      VarG <- VarG[, 1]

      # Backward compatibility for fitted objects created before every kernel
      # received the requested environment covariance structure. Those older
      # objects used a single idv() component for kernels after the first.
      if (length(VarG) != length(heter_grp)) {
        idv_components <- suppressWarnings(as.numeric(LL[, "component"]))
        verified_idv_fallback <- length(names_in_inv_list) > 1L && l > 1L &&
          length(VarG) == 0L && nrow(LL) == 1L &&
          length(idv_components) == 1L && is.finite(idv_components)
        if (!isTRUE(verified_idv_fallback)) {
          stop(
            sprintf(
              paste0(
                "Could not reconstruct %s covariance for kernel '%s': ",
                "expected %d environment-specific variances, found %d."
              ),
              var_cov_str,
              names_in_inv_list[l],
              length(heter_grp),
              length(VarG)
            ),
            call. = FALSE
          )
        }
        idv_var <- idv_components[[1L]]
        VarG <- rep(idv_var, length(heter_grp))
        names(VarG) <- heter_grp
        varcov_all[[l]] <- diag(VarG, nrow = length(heter_grp))
        dimnames(varcov_all[[l]]) <- list(heter_grp, heter_grp)
        corr_all[[l]] <- diag(length(heter_grp))
        dimnames(corr_all[[l]]) <- list(heter_grp, heter_grp)
        VarG_all[[l]] <- VarG
        next
      }

      names(VarG) <- heter_grp

      All_Fac = vector(mode = 'list', length = length(Fac))

      for (i in 1:length(Fac)) {
        TT = LL[grep(paste0("!fa", Fac[i]), rownames(LL)), ]
        All_Fac[[i]] = TT[grep(names_in_inv_list[l], rownames(TT)), "component"]
      }

      All_Fac = do.call(cbind, All_Fac)
      if (!is.matrix(All_Fac) || nrow(All_Fac) != length(heter_grp) ||
          ncol(All_Fac) != length(Fac) || anyNA(All_Fac) || any(!is.finite(All_Fac))) {
        stop(
          sprintf(
            "Could not reconstruct %s loadings for kernel '%s': expected a %d x %d finite loading matrix.",
            var_cov_str,
            names_in_inv_list[l],
            length(heter_grp),
            length(Fac)
          ),
          call. = FALSE
        )
      }
      rownames(All_Fac) <- heter_grp

      varcov_all[[l]] <- gp_factor_analytic_covariance(All_Fac, VarG)
      dimnames(varcov_all[[l]]) <- list(heter_grp, heter_grp)
      corr_all[[l]] <- stats::cov2cor(varcov_all[[l]])
      dimnames(corr_all[[l]]) <- list(heter_grp, heter_grp)

      VarG_all[[l]] = VarG

    }

  } else { ## End of factor analytic model

    ##### For US
    if (var_cov_str=="us"){

      VarCovRaw_All = vector(mode = 'list', length = length(names_in_inv_list))
      VarCovRaw <- vc[grep(paste0("!", heter_groups), rownames(vc)), drop=FALSE, "component"]
      All_varGs = vector(mode = "list", length = length(names_in_inv_list))


      VarG_all = List_list(names_in_inv_list = names_in_inv_list, het_gp = heter_grp)

      for(l in 1:length(names_in_inv_list)){
        LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
        for (G in 1:length(heter_grp)) {
          VarG_all[[l]][[G]] = LL[grep(heter_grp[G], rownames(LL)), drop=FALSE,]
        }

      }

      for (g in 1:length(names_in_inv_list)) {
            VarCovRaw_All[[g]] <- VarCovRaw[grep(names_in_inv_list[g], rownames(VarCovRaw)), "component"]

        varcov_all[[g]] <- gp_unstructured_covariance_from_parameters(
          VarCovRaw_All[[g]], n = n_heter_grp, labels = heter_grp
        )
        CORR <- stats::cov2cor(varcov_all[[g]])
        corr_all[[g]] <-  CORR
        dimnames(varcov_all[[g]]) <-  list(heter_grp , heter_grp )
        dimnames(corr_all[[g]]) <-  list(heter_grp , heter_grp )

      }

    }
    # End of us

    #### corgh

    if (var_cov_str=="corgh"){
          var_corr_all = vector("list", length = length(names_in_inv_list))
          CheckR <- vc[grep("!R", rownames(vc)), "component"]
          VarG_all = vector("list", length = length(names_in_inv_list))

      for(l in 1:length(names_in_inv_list)){
        LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
        ######
        #CORGH
        var_corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]
        VarG <- LL[grep(paste0(heter_groups, "_"), rownames(LL)), "component"]
        names(VarG) = heter_grp

        corr_all[[l]] <- gp_correlation_matrix_from_parameters(
          var_corr[, 1], n = n_heter_grp, labels = heter_grp
        )
        varcov_all[[l]] <- gp_covariance_from_correlation(
          VarG, corr_all[[l]], labels = heter_grp
        )
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }
    } ## end of corgh

    ### corh
    if (var_cov_str=="corh"){
        var_corr_all = vector("list", length = length(names_in_inv_list))
        CheckR <- vc[grep("!R", rownames(vc)), "component"]
        VarG_all = vector("list", length = length(names_in_inv_list))

      for(l in 1:length(names_in_inv_list)){
        LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
        ######
        var_corr <- LL[grep(".cor", rownames(LL)),drop=TRUE, "component"]
        VarG <- LL[-grep(".cor", rownames(LL)),drop=TRUE, "component"]
        corr_all[[l]] <- gp_correlation_matrix_from_parameters(
          var_corr, n = n_heter_grp, compound_symmetry = TRUE, labels = heter_grp
        )
        varcov_all[[l]] <- gp_covariance_from_correlation(
          VarG, corr_all[[l]], labels = heter_grp
        )
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }

    } # End of corh

    ## corv
    if (var_cov_str=="corv"){
          var_corr_all = vector("list", length = length(names_in_inv_list))
          CheckR <- vc[grep("!R", rownames(vc)), "component"]
          VarG_all = vector("list", length = length(names_in_inv_list))

      for(l in 1:length(names_in_inv_list)){
        LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
        ######
        var_corr <- LL[grep(".cor", rownames(LL)),drop=TRUE, "component"]
        VarG <- LL[grep(paste0(heter_groups, "!var"), rownames(LL)), drop = TRUE,"component"]
        VarG = rep( VarG, n_heter_grp)
        corr_all[[l]] <- gp_correlation_matrix_from_parameters(
          var_corr, n = n_heter_grp, compound_symmetry = TRUE, labels = heter_grp
        )
        varcov_all[[l]] <- gp_covariance_from_correlation(
          VarG, corr_all[[l]], labels = heter_grp
        )
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }
    }### End corv
 }

  ### Calculate genomic Heritability for Each Location for any of the variance-covariance structure
  residual_rows <- grep("!R", rownames(vc))
  VarE <- suppressWarnings(as.numeric(vc[residual_rows, "component"]))
  if(isTRUE(heter_resid)){
    if (length(VarE) != length(heter_grp)) {
      stop(
        sprintf(
          "Expected one ASReml residual variance for each of %d environments; found %d.",
          length(heter_grp), length(VarE)
        ),
        call. = FALSE
      )
    }
    names(VarE) <- heter_grp

  } else {
    if (length(VarE) < 1L || !is.finite(VarE[[1L]])) {
      stop("ASReml did not return a finite residual variance.", call. = FALSE)
    }
    VarE <- rep(VarE[[1L]], length(heter_grp))

    names(VarE) <-  heter_grp
  }
  if (any(!is.finite(VarE)) || any(VarE < 0)) {
    stop("ASReml residual variances must be finite and non-negative.", call. = FALSE)
  }
  vc_se <- if ("std.error" %in% names(vc)) {
    suppressWarnings(as.numeric(vc[["std.error"]]))
  } else {
    rep(NA_real_, nrow(vc))
  }
  residual_scale_se <- vc_se[residual_rows]
  if (!isTRUE(heter_resid) && length(residual_scale_se)) {
    residual_scale_se <- rep(residual_scale_se[[1L]], length(heter_grp))
  }
  if (length(residual_scale_se) != length(heter_grp)) {
    residual_scale_se <- rep(NA_real_, length(heter_grp))
  }
  names(residual_scale_se) <- heter_grp

  residual_contract <- gp_asreml_residual_variance_contract(
    model = model,
    heter_groups = heter_groups,
    heter_grp = heter_grp,
    residual_scale = VarE,
    residual_scale_se = residual_scale_se,
    response = response,
    gen_name = gen_name
  )
  VarE <- residual_contract$residual_variance
  residual_se <- residual_contract$residual_variance_se

  #VarG = diag(VarCov)
  H = matrix(NA, nrow = 1, ncol = length(heter_grp))
  colnames(H) = heter_grp
  varG_per_omics = matrix(NA, nrow = length(names_in_inv_list), ncol = length(heter_grp))
  dimnames(varG_per_omics) <- list(names_in_inv_list, heter_grp)
  Total_varG = matrix(NA, nrow = 1, ncol = length(heter_grp))
  colnames(Total_varG) = heter_grp

  for (va in 1:length(names_in_inv_list)) {
    varG_per_omics[va, ] <- diag(varcov_all[[va]])

  }

  # Rigor: per-environment heritability with a delta-method SE where it is exact.
  # The per-env genetic variance comes from the diagonal of the reconstructed
  # variance-covariance matrix. For structures whose diagonal entries ARE fitted
  # variance-component rows (e.g. corgh, us) the SE is exact; for factor-analytic
  # structures the diagonal is Lambda Lambda' + psi, which matches no single row,
  # so gp_asreml_h2_se_for_env() self-checks and returns NA rather than a wrong SE.
  H_SE = matrix(NA_real_, nrow = 1, ncol = length(heter_grp))
  colnames(H_SE) = heter_grp

  for (k in 1:ncol(H)) {
    varG = c()
    for (v in 1:length(names_in_inv_list)) {
      varG[[v]] = varG_per_omics[v, k]

    }
### since omics can be more than 1
    Total_varG[, k] <- sum(unlist(varG))
    ve_k <- VarE[min(k, length(VarE))]
    denominator <- sum(unlist(varG)) + ve_k
    H[, k] <- if (is.finite(denominator) && denominator > 0) {
      sum(unlist(varG)) / denominator
    } else {
      NA_real_
    }

    H_SE[, k] <- gp_asreml_h2_se_for_env(model, vc = vc,
                                         names_in_inv_list = names_in_inv_list,
                                         env_label = heter_grp[k],
                                         target_vg = Total_varG[, k],
                                         target_ve = as.numeric(ve_k))

  }

  total_genetic_covariance <- Reduce(`+`, varcov_all)
  dimnames(total_genetic_covariance) <- list(heter_grp, heter_grp)
  total_diagnostics <- relationship_matrix_diagnostics(
    total_genetic_covariance,
    require_names = TRUE
  )
  if (!isTRUE(total_diagnostics$valid)) {
    stop(
      paste(
        "ASReml fitted genetic covariance matrix is invalid:",
        paste(total_diagnostics$issues, collapse = "; ")
      ),
      call. = FALSE
    )
  }
  total_genetic_correlation <- stats::cov2cor(total_genetic_covariance)
  dimnames(total_genetic_correlation) <- list(heter_grp, heter_grp)

  residual_covariance <- diag(VarE, nrow = length(VarE))
  dimnames(residual_covariance) <- list(heter_grp, heter_grp)
  residual_correlation <- diag(length(VarE))
  dimnames(residual_correlation) <- list(heter_grp, heter_grp)

  varG_per_omics_se <- matrix(
    NA_real_,
    nrow = nrow(varG_per_omics),
    ncol = ncol(varG_per_omics),
    dimnames = dimnames(varG_per_omics)
  )
  for (source_idx in seq_along(names_in_inv_list)) {
    kernel_rows <- grepl(names_in_inv_list[[source_idx]], rownames(vc), fixed = TRUE)
    excluded_rows <- grepl("\\.cor$|!cor$|!fa[0-9]+$", rownames(vc))
    candidates <- which(kernel_rows & !excluded_rows)
    for (env_idx in seq_along(heter_grp)) {
      target <- varG_per_omics[source_idx, env_idx]
      matching <- candidates[
        is.finite(vc_se[candidates]) &
          abs(suppressWarnings(as.numeric(vc[candidates, "component"])) - target) <=
            1e-8 + 1e-5 * abs(target)
      ]
      if (length(matching) == 1L) {
        varG_per_omics_se[source_idx, env_idx] <- vc_se[matching]
      }
    }
  }
  ### This is to create variance-covariance and correlation matrix
  ## for any number of omics that the user provide
  return (list(Covariance = varcov_all,
                 Correlation = corr_all,
                 Genetic_covariance_environments = total_genetic_covariance,
                 Genetic_correlation_environments = total_genetic_correlation,
                 Residual_covariance_environments = residual_covariance,
                 Residual_correlation_environments = residual_correlation,
                 Heritability = H,
                 Heritability_SE = H_SE,
                 Total_genetic_var = Total_varG,
                 varG_per_omics = varG_per_omics,
                 varG_per_omics_SE = varG_per_omics_se,
                 Residual_Var = VarE,
                 Residual_Var_SE = residual_se,
                 Residual_scale_parameter = residual_contract$residual_scale,
                 Residual_scale_parameter_SE = residual_contract$residual_scale_se,
                 Residual_variance_by_environment = residual_contract$residual_variance_by_environment,
                 Residual_variance_by_observation = residual_contract$residual_variance_by_observation,
                 Residual_variance_basis = residual_contract$residual_variance_basis,
                 Estimation_method = paste0("ASReml_REML_", var_cov_str)))


}


