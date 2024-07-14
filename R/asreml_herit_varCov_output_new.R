

#' Title
#'
#' @param names_in_inv_list
#' @param het_gp
#'
#' @return
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

#' Calculate Heritability and extract Variance component, Covariance Structures from ASReml Model
#'
#' This function calculates heritability and extracts variance-covariance and correlation structures from an ASReml model. It is designed to work with models that may include heterogeneity in residuals and supports the inclusion of multiple genetic or omics components. The function also handles different variance-covariance structures specified by the user.
#'
#' @param model An `asreml` model object resulting from fitting an ASReml model.
#' @param heter_groups A character string specifying the column in the model's data frame that represents different heterogeneity groups (e.g., environments).
#' @param var_cov_str A character string indicating the type of variance-covariance structure to be analyzed (`"us"`, `"fa"`, `"rr"`, `"corgh"`, `"corgv"`, `"corh"`, or `"corv"`).
#' @param heter_resid A logical indicating whether heterogeneity in residuals is considered.
#' @param names_in_inv_list A character vector of names identifying the inverse matrices in the model, related to different genetic or omics components.
#' @param inter_gen_pos Optional parameter specifying the position of interaction terms in genetic models.
#' @param gen_pos Optional parameter specifying the position of genetic terms in the model.
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
#' The function primarily focuses on extracting and computing key genetic statistics from a given ASReml model. It supports a range of variance-covariance structures, including but not limited to unstructured (`us`), factor analytic (`fa`), and various correlation structures (`corgh`, `corgv`, `corh`, `corv`). The function checks for stability in variance components and requires a stable model for accurate computations.
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
    ...

){

  msg <- "\n==================================================\n"

  vc <- asreml::summary.asreml(model)$varcomp

  ENV <- data.frame(model$mf)[, heter_groups]
  heter_grp = levels(ENV)
  n_heter_grp <- nlevels(ENV)

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
  names(varcov_all) <- paste("covariance", extracted_names_from_inv_list, sep = "_")
  names(corr_all) <- paste("correlation", extracted_names_from_inv_list, sep = "_")
  }
################
  if(length(names_in_inv_list)==1){
    names(varcov_all) <-"covariance"
    names(corr_all) <- "correlation"

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
      names(VarG) <- heter_grp

      All_Fac = vector(mode = 'list', length = length(Fac))

      for (i in 1:length(Fac)) {
        TT = LL[grep(paste0("!fa", Fac[i]), rownames(LL)), ]
        All_Fac[[i]] = TT[grep(names_in_inv_list[l], rownames(TT)), "component"]
      }

      All_Fac = do.call(cbind, All_Fac)

      varcov_all[[l]] <- All_Fac %*% t(All_Fac) + diag(VarG)
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

        a <- 1
        for (r in 1:n_heter_grp) {
          for (c in 1:r) {
            varcov_all[[g]][r, c] <-VarCovRaw_All[[g]][a]
            varcov_all[[g]][c, r] <- VarCovRaw_All[[g]][a]
            a <- a + 1
          }
        }

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

        a <- 1
        for (r in 1:n_heter_grp) {
          for (c in 1:r) {
            if (r == c) {
              corr_all[[l]][r, c] <- 1
            } else {
              corr_all[[l]][r, c] <-  var_corr[a, 1]
              corr_all[[l]][c, r] <- var_corr[a, 1]
              a <- a + 1

            }
          }
        }

        varcov_all[[l]] <- diag(sqrt(VarG)) %*% corr_all[[l]] %*% diag(sqrt(VarG))
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }
    } ## end of corgh

    ## start corgv
    if (var_cov_str=="corgv"){
        var_corr_all = vector("list", length = length(names_in_inv_list))
        CheckR <- vc[grep("!R", rownames(vc)), "component"]
        VarG_all = vector("list", length = length(names_in_inv_list))

      for(l in 1:length(names_in_inv_list)){
          LL = vc[grep(names_in_inv_list[l], rownames(vc)), drop=FALSE,]
          ######
          var_corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]
          VarG <- LL[grep(paste0(heter_groups, "!var"), rownames(LL)), drop = TRUE,"component"]
          VarG = rep( VarG, n_heter_grp)

        a <- 1
        for (r in 1:n_heter_grp) {
          for (c in 1:r) {
            if (r == c) {
              corr_all[[l]][r, c] <- 1
            } else {
              corr_all[[l]][r, c] <-  var_corr[a, 1]
              corr_all[[l]][c, r] <- var_corr[a, 1]
              a <- a + 1

            }
          }
        }

        varcov_all[[l]] <- diag(sqrt(VarG)) %*% corr_all[[l]] %*% diag(sqrt(VarG))
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }

    } ### End corgv

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
        var_corr = rep( var_corr, n_heter_grp)
        corr_all[[l]][1:n_heter_grp, 1:n_heter_grp] <-  var_corr
        varcov_all[[l]] <- diag(sqrt(VarG)) %*% corr_all[[l]] %*% diag(sqrt(VarG))
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
        var_corr = rep( var_corr, n_heter_grp)
        corr_all[[l]][1:n_heter_grp, 1:n_heter_grp] <-  var_corr
        diag(corr_all[[l]]) <- 1
        varcov_all[[l]] <- diag(sqrt(VarG)) %*% corr_all[[l]] %*% diag(sqrt(VarG))
        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)
        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)
        VarG_all[[l]] = VarG

      }
    }### End corv
 }

  ### Calculate genomic Heritability for Each Location for any of the variance-covariance structure
  VarE = vc[grep("!R", rownames(vc)), "component"]
  if(isTRUE(heter_resid)){
  names(VarE) = heter_grp

  } else {
    VarE <- array(rep(VarE[1], length(heter_grp)))

    names(VarE) <-  heter_grp
  }

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


  for (k in 1:ncol(H)) {
    varG = c()
    for (v in 1:length(names_in_inv_list)) {
      varG[[v]] = varG_per_omics[v, k]

    }
### since omics can be more than 1
    Total_varG[, k] <- sum(unlist(varG))
    if(isTRUE(heter_resid) & !is.null(inter_gen_pos)){

      H[, k] <- sum(unlist(varG))/(sum(unlist(varG))+VarE[k])

    } else {

      if((isFALSE(heter_resid) | is.null(heter_resid)) & is.null(inter_gen_pos)){

        H[, k] <- sum(unlist(varG))/(sum(unlist(varG))+VarE)
      }

    }

  }

  ### This is to create variance-covariance and correlation matrix
  ## for any number of omics that the user provide
  return (list(Covariance = varcov_all,
                 Correlation = corr_all,
                 Heritability = H,
                 Total_genetic_var = Total_varG,
                 varG_per_omics = varG_per_omics,
                 Residual_Var = VarE))


}


