
#' Title
#'
#' @param model
#' @param heter_groups
#' @param var_cov_str
#' @param heter_resid
#' @param G_list
#' @param inter_gen_pos
#' @param gen_pos
#' @param ...
#'
#' @return
#'
#'
#' @examples
#'

asreml_herit_varCov <-  function(
    model = NULL,
    heter_groups = NULL,
    var_cov_str= NULL,
    heter_resid=FALSE,
    G_list = NULL,
    inter_gen_pos = NULL,
    gen_pos = NULL,
    ...

){

  msg <- "\n==================================================\n"

  vc <- asreml::summary.asreml(model)$varcomp

  #heter_grp <- as.character(unique(data.frame(model$mf)[, heter_groups]))

  ENV <- data.frame(model$mf)[, heter_groups]

  heter_grp = levels(ENV)

  #n_heter_grp <- length(heter_grp)

  n_heter_grp <- nlevels(ENV)

  VarCov <- matrix(NA, ncol = n_heter_grp, nrow = n_heter_grp)

  CORR <- matrix(NA, ncol = n_heter_grp, nrow = n_heter_grp)

  corr_all = vector(mode = 'list', length = length(G_list))

  varcov_all = vector(mode = 'list', length = length(G_list))

  names(varcov_all) <- paste(G_list, "covariance", sep = "_")

  names(corr_all) <- paste(G_list, "correlation", sep = "_")

  for (ca in 1:length(corr_all)) {
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

      stop(print(paste(msg,paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
                       " is unstable, refix the model")), call. = FALSE)

    }

  }

  #### This needs to be turn on
  # if (any(is.na(vc$std.error))) {
  #   stop(print(paste(msg, "Some variance component are non estimatable. Refix the model")), call. = FALSE)
  # }

  #################################################################


  #### Extract variance and covariance for For factor analytic models
  if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
    ## Check for all variable is positive definitive
    ## Check for this other random term can be present aside the genetic effect
    #VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE),"bound"]


    CheckR <- vc[grep("!R", rownames(vc)), "component"]

    VarG_all = vector("list", length = length(G_list))

    ## Extract number of factor(s) specified by users
    N_fa = as.double(substr(var_cov_str, 3, 100))

    Fac = seq(1, N_fa)


    for(l in 1:length(G_list)){

      LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]

      VarG <- LL[grep("!var", rownames(LL)), drop=FALSE, ]

      VarG <- VarG[, 1]

      names(VarG) <- heter_grp

      #FA_All = vector(mode = 'list', length = length(G_list))

      #for (fa in 1:length(G_list)){

      All_Fac = vector(mode = 'list', length = length(Fac))

      for (i in 1:length(Fac)) {

        TT = LL[grep(paste0("!fa", Fac[i]), rownames(LL)), ]


        All_Fac[[i]] = TT[grep(paste(G_list[[l]],"_inv", sep = ""), rownames(TT)), "component"]


      }


      All_Fac = do.call(cbind, All_Fac)




      varcov_all[[l]] <- All_Fac %*% t(All_Fac) + diag(VarG)
      dimnames(varcov_all[[l]]) <- list(heter_grp, heter_grp)
      corr_all[[l]] <- stats::cov2cor(varcov_all[[l]])
      dimnames(corr_all[[l]]) <- list(heter_grp, heter_grp)

      VarG_all[[l]] = VarG

    }

    # names(varcov_all) <- paste(G_list, "covariance", sep = "_")
    # names(corr_all) <- paste(G_list, "correlation", sep = "_")


  } else { ## End of factor analytic model

    ##### For US
    if (var_cov_str=="us"){

      VarCovRaw_All = vector(mode = 'list', length = length(G_list))

      VarCovRaw <- vc[grep(paste0("!", heter_groups), rownames(vc)), drop=FALSE, "component"]

      All_varGs = vector(mode = "list", length = length(G_list))

      List_list = function(G_list, het_gp){

        All_varGs = vector(mode = "list", length = length(G_list))

        for (l in 1:length(G_list)) {

          All_varGs[[l]] <- vector(mode ="list", length(het_gp))

        }

        return(All_varGs)
      }


      VarG_all = List_list(G_list = G_list, het_gp = heter_grp)


      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]

        for (G in 1:length(heter_grp)) {

          VarG_all[[l]][[G]] = LL[grep(heter_grp[G], rownames(LL)), drop=FALSE,]
        }


      }



      for (g in 1:length(G_list)) {


        VarCovRaw_All[[g]] <- VarCovRaw[grep(paste(G_list[[g]],"_inv", sep = ""), rownames(VarCovRaw)), "component"]

        #}

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

      # dimnames(varcov_all[[g]]) <-  list(heter_grp , heter_grp )
      # dimnames(corr_all[[g]]) <-  list(heter_grp , heter_grp )

    }
    # End of us

    #### corgh

    if (var_cov_str=="corgh"){

      var_corr_all = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_all = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        var_corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]


        VarG <- LL[grep(paste0(heter_groups, "_"), rownames(LL)), "component"]

        names(VarG) = heter_grp
        #for (caa in 1:length(G_list)) {


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
        #} # End of Corgh
      }

      # names(varcov_all) <- paste(G_list, "covariance", sep = "_")
      # names(corr_all) <- paste(G_list, "correlation", sep = "_")

    } ## end of corgh

    if (var_cov_str=="corgv"){

      var_corr_all = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_all = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        var_corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]


        VarG <- LL[grep(paste0(heter_groups, "!var"), rownames(LL)), drop = TRUE,"component"]

        VarG = rep( VarG, n_heter_grp)
        #names(VarG) = heter_grp
        #for (caa in 1:length(G_list)) {


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
        #} # End of Corgh
      }

      # names(varcov_all) <- paste(G_list, "covariance", sep = "_")
      # names(corr_all) <- paste(G_list, "correlation", sep = "_")

    } ### End corgv

    ### corh

    if (var_cov_str=="corh"){

      var_corr_all = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_all = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        var_corr <- LL[grep(".cor", rownames(LL)),drop=TRUE, "component"]


        VarG <- LL[-grep(".cor", rownames(LL)),drop=TRUE, "component"]

        var_corr = rep( var_corr, n_heter_grp)


        corr_all[[l]][1:n_heter_grp, 1:n_heter_grp] <-  var_corr





        varcov_all[[l]] <- diag(sqrt(VarG)) %*% corr_all[[l]] %*% diag(sqrt(VarG))

        dimnames(varcov_all[[l]]) <-  list(heter_grp, heter_grp)

        dimnames(corr_all[[l]]) <-  list(heter_grp, heter_grp)

        VarG_all[[l]] = VarG
        #} # End of corv
      }

      # names(varcov_all) <- paste(G_list, "covariance", sep = "_")
      # names(corr_all) <- paste(G_list, "correlation", sep = "_")
    }

    ## end corh


    if (var_cov_str=="corv"){

      var_corr_all = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_all = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

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
        #} # End of corv
      }

      # names(varcov_all) <- paste(G_list, "covariance", sep = "_")
      # names(corr_all) <- paste(G_list, "correlation", sep = "_")
    }

    ### End corv


  }

  ### Calculate genomic Heritability for Each Location for any of the variance-covariance structure
  VarE = vc[grep("!R", rownames(vc)), "component"]

  names(VarE) = heter_grp

  #VarG = diag(VarCov)
  H = matrix(NA, nrow = 1, ncol = length(heter_grp))

  colnames(H) = heter_grp

  varG_per_omics = matrix(NA, nrow = length(G_list), ncol = length(heter_grp))

  dimnames(varG_per_omics) <- list(unlist(G_list), heter_grp)

  Total_varG = matrix(NA, nrow = 1, ncol = length(heter_grp))

  colnames(Total_varG) = heter_grp



  for (va in 1:length(G_list)) {

    #varG_per_omics[va, ] <- VarG_all[[va]]

    varG_per_omics[va, ] <- diag(varcov_all[[va]])


  }


  for (k in 1:ncol(H)) {

    varG = c()

    for (v in 1:length(G_list)) {

      #varG[[v]] = VarG_all[[v]][k, 1]
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

  ### This is to create variance-covariance and correlation matrix
  ## for every any number of G that the user provide


      Result <- list(Covariance = varcov_all,
                    Correlation = corr_all,
                    Heritability = H,
                    Total_genetic_var = Total_varG,
                    varG_per_omics = varG_per_omics,
                    Residual_Var = VarE)

      # names(Result) <-  c("Covariance",
      #                     "Correlation",
      #                     "Heritability",
      #                     "Total_genetic_var",
      #                     "varG_per_omics",
      #                     "Residual_Var")


  return(Result)

}


