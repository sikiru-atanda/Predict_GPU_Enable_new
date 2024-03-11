
#' Title
#'
#' @param model
#' @param heter_groups
#' @param heter_resid
#' @param G_list
#' @param inter_gen_pos
#' @param gen_pos
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
asreml_herit_CSM <-  function(
    model = mod,
    heter_groups = NULL,
    heter_resid=FALSE,
    G_list = NULL,
    inter_gen_pos = NULL,
    gen_pos = NULL,
    ...

){

  msg <- sprintf("==================================================\n")

  vc <- asreml::summary.asreml(model)$varcomp

  ENV <- data.frame(model$mf)[, heter_groups]

  heter_grp = levels(ENV)

  n_heter_grp <- nlevels(ENV)

  VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE), ]

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

      stop(print(paste(msg, paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
                 " is unstable, refix the model")), call. = FALSE)

    }

  }

  if (any(is.na(vc$std.error))) {
    stop(print(paste(msg, "Some variance component are non estimatable. Refix the model")), call. = FALSE)
  }

  VarE  <- vc[grep("!R", rownames(vc)), "component"]
  names(VarE) = heter_grp

  varG_per_omics = matrix(NA, nrow = length(G_list), ncol = length(heter_grp))

  dimnames(varG_per_omics) <- list(unlist(G_list), heter_grp)

  for (g in 1:length(G_list)) {

  VarG = vc[grep(paste(G_list[[g]], "inv", sep = "_"), rownames(vc)), drop = TRUE,"component"]

  VarG = rep(VarG, n_heter_grp)

  varG_per_omics[g, ] <-  VarG

  }

  #VarG = rep(VarG, n_heter_grp)

  H = matrix(NA, nrow = 1, ncol = length(heter_grp))

  colnames(H) = heter_grp

  Total_varG = matrix(NA, nrow = 1, ncol = length(heter_grp))

  colnames(Total_varG) = colnames(heter_grp)


  for (k in 1:ncol(H)) {

    varG = c()

    for (v in 1:length(G_list)) {

      #varG[[v]] = VarG_All[[v]][k, 1]
      varG[[v]] = varG_per_omics[v, k]

    }

    Total_varG[, k] <- sum(unlist(varG))

    if(isTRUE(heter_resid) & !is.null(inter_gen_pos)){


      H[, k] <- sum(unlist(varG))/(sum(unlist(varG))+VarE[k])

    } else {

      if((isFALSE(heter_resid) | is.null(heter_resid)) & is.null(inter_gen_pos)){

        H[, k] <- sum(unlist(varG))/(sum(unlist(varG))+VarE)
      }

    }

  }


  Result = list(Heritability = H,
                Total_genetic_var = Total_varG,
                varG_per_omics = varG_per_omics,
                Residual_Var = VarE)

names(Result) <- c("Heritability",
                   "Total_genetic_var",
                    "varG_per_omics",
                     "Residual_Var")

  return(Result)

}


