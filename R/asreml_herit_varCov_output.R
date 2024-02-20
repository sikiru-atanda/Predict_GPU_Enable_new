
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

asreml_herit_varCovraw <-  function(
    model = NULL,
    heter_groups = NULL,
    var_cov_str= NULL,
    heter_resid=FALSE,
    G_list = NULL,
    inter_gen_pos = NULL,
    gen_pos = NULL,
    ...

){

  msg <- sprintf("==================================================\n")

  vc <- asreml::summary.asreml(model)$varcomp

  #Heter.Grp <- as.character(unique(data.frame(model$mf)[, heter_groups]))

  ENV <- data.frame(model$mf)[, heter_groups]

  Heter.Grp = levels(ENV)

  #N.Heter.Grp <- length(Heter.Grp)

  N.Heter.Grp <- nlevels(ENV)

  VarCov <- matrix(NA, ncol = N.Heter.Grp, nrow = N.Heter.Grp)

  CORR <- matrix(NA, ncol = N.Heter.Grp, nrow = N.Heter.Grp)

  CORR_ALL = vector(mode = 'list', length = length(G_list))

  VarCov_All = vector(mode = 'list', length = length(G_list))

  for (ca in 1:length(CORR_ALL)) {
    CORR_ALL[[ca]] <-  CORR

    VarCov_All[[ca]] <- VarCov

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

    # stop(print(paste(paste("variance component for", as.character(Heter.Grp[VAR_check_Pos]), collapse = " and "),
    #                    "are unstable, refit the model")), call. = FALSE)

    stop(print(paste(msg, paste(paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = ","),
                     "are unstable, refit the model"))), call. = FALSE)
  } else {

    if(length(VAR_check_Pos)==1) {

      stop(print(paste(msg,paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
                 " is unstable, refix the model")), call. = FALSE)

    }

  }

  if (any(is.na(vc$std.error))) {
    stop(print(paste(msg, "Some variance component are non estimatable. Refix the model")), call. = FALSE)
  }

  #################################################################


  #### Extract variance and covariance for For factor analytic models
  if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
    ## Check for all variable is positive definitive
    ## Check for this other random term can be present aside the genetic effect
    #VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE),"bound"]


    CheckR <- vc[grep("!R", rownames(vc)), "component"]

    VarG_All = vector("list", length = length(G_list))

    ## Extract number of factor(s) specified by users
    N_fa = as.double(substr(var_cov_str, 3, 100))

    Fac = seq(1, N_fa)


    for(l in 1:length(G_list)){

      LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]

      VarG <- LL[grep("!var", rownames(LL)), drop=FALSE, ]

      VarG <- VarG[, 1]

      names(VarG) <- Heter.Grp

      #FA_All = vector(mode = 'list', length = length(G_list))

      #for (fa in 1:length(G_list)){

      All_Fac = vector(mode = 'list', length = length(Fac))

      for (i in 1:length(Fac)) {

        TT = LL[grep(paste0("!fa", Fac[i]), rownames(LL)), ]


        All_Fac[[i]] = TT[grep(paste(G_list[[l]],"_inv", sep = ""), rownames(TT)), "component"]


      }


      All_Fac = do.call(cbind, All_Fac)




      VarCov_All[[l]] <- All_Fac %*% t(All_Fac) + diag(VarG)
      dimnames(VarCov_All[[l]]) <- list(Heter.Grp, Heter.Grp)
      CORR_ALL[[l]] <- stats::cov2cor(VarCov_All[[l]])
      dimnames(CORR_ALL[[l]]) <- list(Heter.Grp, Heter.Grp)

      VarG_All[[l]] = VarG

    }

    names(VarCov_All) = G_list

    names(CORR_ALL) = G_list


  } else { ## End of factor analytic model

    ##### For US
    if (var_cov_str=="us"){

      ## Check for this other random term can be present aside the genetic effect
      #VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE),"bound"]

      # VAR_check <- vc[grep(paste0("!", heter_groups), rownames(vc), value = FALSE), ]
      #
      # #VAR_check_Pos = which(vc[, "bound"]=="F" | vc[, "bound"]=="U" |vc[, "bound"] =="?" |vc[, "bound"] =="S" )
      # VAR_check_Pos = which(VAR_check[, "bound"]=="F" | VAR_check[, "bound"]=="U" |VAR_check[, "bound"] =="?" |VAR_check[, "bound"] =="S" )
      #
      #
      # if(length(VAR_check_Pos)>1) {
      #
      #   # stop(print(paste(paste("variance component for", as.character(Heter.Grp[VAR_check_Pos]), collapse = " and "),
      #   #                    "are unstable, refit the model")), call. = FALSE)
      #
      #   stop(print(paste(paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = ","),
      #                    "are unstable, refit the model")), call. = FALSE)
      # } else {
      #
      #   if(length(VAR_check_Pos)==1) {
      #
      #     stop(print(paste("variance component for", as.character(rownames(VAR_check)[VAR_check_Pos]), collapse = " "),
      #                  " is unstable, refix the model"), call. = FALSE)
      #
      #   }
      #
      # }

      # het_gp = paste(paste(heter_groups, "_", sep = ""),
      #                paste(Heter.Grp, Heter.Grp, sep = ":"), sep = "")



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


     VarG_All = List_list(G_list = G_list, het_gp = Heter.Grp)


      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]

        for (G in 1:length(Heter.Grp)) {

          VarG_All[[l]][[G]] = LL[grep(Heter.Grp[G], rownames(LL)), drop=FALSE,]
        }


      }



     for (g in 1:length(G_list)) {


       VarCovRaw_All[[g]] <- VarCovRaw[grep(paste(G_list[[g]],"_inv", sep = ""), rownames(VarCovRaw)), "component"]

       #}

       a <- 1
       for (r in 1:N.Heter.Grp) {
         for (c in 1:r) {
           VarCov_All[[g]][r, c] <-VarCovRaw_All[[g]][a]
           VarCov_All[[g]][c, r] <- VarCovRaw_All[[g]][a]
           a <- a + 1
         }
       }

       CORR <- stats::cov2cor(VarCov_All[[g]])

       CORR_ALL[[g]] <-  CORR

     }

     dimnames(VarCov_All[[g]]) <-  list(Heter.Grp , Heter.Grp )

     dimnames(CORR_ALL[[g]]) <-  list(Heter.Grp , Heter.Grp )

    }
     # End of us

     #### corgh

    if (var_cov_str=="corgh"){

      Var.corr_All = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_All = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        Var.corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]


        VarG <- LL[grep(paste0(heter_groups, "_"), rownames(LL)), "component"]

        names(VarG) = Heter.Grp
        #for (caa in 1:length(G_list)) {


        a <- 1
        for (r in 1:N.Heter.Grp) {
          for (c in 1:r) {
            if (r == c) {
              CORR_ALL[[l]][r, c] <- 1
            } else {
              CORR_ALL[[l]][r, c] <-  Var.corr[a, 1]
              CORR_ALL[[l]][c, r] <- Var.corr[a, 1]
              a <- a + 1

            }
          }
        }

        VarCov_All[[l]] <- diag(sqrt(VarG)) %*% CORR_ALL[[l]] %*% diag(sqrt(VarG))

        dimnames(VarCov_All[[l]]) <-  list(Heter.Grp, Heter.Grp)

        dimnames(CORR_ALL[[l]]) <-  list(Heter.Grp, Heter.Grp)

        VarG_All[[l]] = VarG
        #} # End of Corgh
      }

      names(VarCov_All) = G_list

      names(CORR_ALL) = G_list

    } ## end of corgh

    if (var_cov_str=="corgv"){

      Var.corr_All = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_All = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        Var.corr <- LL[grep(".cor", rownames(LL)),drop=FALSE, "component"]


        VarG <- LL[grep(paste0(heter_groups, "!var"), rownames(LL)), drop = TRUE,"component"]

        VarG = rep( VarG, N.Heter.Grp)
        #names(VarG) = Heter.Grp
        #for (caa in 1:length(G_list)) {


        a <- 1
        for (r in 1:N.Heter.Grp) {
          for (c in 1:r) {
            if (r == c) {
              CORR_ALL[[l]][r, c] <- 1
            } else {
              CORR_ALL[[l]][r, c] <-  Var.corr[a, 1]
              CORR_ALL[[l]][c, r] <- Var.corr[a, 1]
              a <- a + 1

            }
          }
        }

        VarCov_All[[l]] <- diag(sqrt(VarG)) %*% CORR_ALL[[l]] %*% diag(sqrt(VarG))

        dimnames(VarCov_All[[l]]) <-  list(Heter.Grp, Heter.Grp)

        dimnames(CORR_ALL[[l]]) <-  list(Heter.Grp, Heter.Grp)

        VarG_All[[l]] = VarG
        #} # End of Corgh
      }

      names(VarCov_All) = G_list

      names(CORR_ALL) = G_list

    } ### End corgv

### corh

    if (var_cov_str=="corh"){

      Var.corr_All = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_All = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        Var.corr <- LL[grep(".cor", rownames(LL)),drop=TRUE, "component"]


        VarG <- LL[-grep(".cor", rownames(LL)),drop=TRUE, "component"]

        Var.corr = rep( Var.corr, N.Heter.Grp)


        CORR_ALL[[l]][1:N.Heter.Grp, 1:N.Heter.Grp] <-  Var.corr





        VarCov_All[[l]] <- diag(sqrt(VarG)) %*% CORR_ALL[[l]] %*% diag(sqrt(VarG))

        dimnames(VarCov_All[[l]]) <-  list(Heter.Grp, Heter.Grp)

        dimnames(CORR_ALL[[l]]) <-  list(Heter.Grp, Heter.Grp)

        VarG_All[[l]] = VarG
        #} # End of corv
      }

      names(VarCov_All) = G_list

      names(CORR_ALL) = G_list
    }

    ## end corh


    if (var_cov_str=="corv"){

      Var.corr_All = vector("list", length = length(G_list))

      CheckR <- vc[grep("!R", rownames(vc)), "component"]

      VarG_All = vector("list", length = length(G_list))

      for(l in 1:length(G_list)){

        LL = vc[grep(paste(G_list[l], "inv", sep = "_"), rownames(vc)), drop=FALSE,]
        ######
        ## Ideal for the extraction extrapolated from Johan package Agriutilities
        #CORGH

        Var.corr <- LL[grep(".cor", rownames(LL)),drop=TRUE, "component"]


        VarG <- LL[grep(paste0(heter_groups, "!var"), rownames(LL)), drop = TRUE,"component"]

        VarG = rep( VarG, N.Heter.Grp)

        Var.corr = rep( Var.corr, N.Heter.Grp)


        CORR_ALL[[l]][1:N.Heter.Grp, 1:N.Heter.Grp] <-  Var.corr

        diag(CORR_ALL[[l]]) <- 1



        VarCov_All[[l]] <- diag(sqrt(VarG)) %*% CORR_ALL[[l]] %*% diag(sqrt(VarG))

        dimnames(VarCov_All[[l]]) <-  list(Heter.Grp, Heter.Grp)

        dimnames(CORR_ALL[[l]]) <-  list(Heter.Grp, Heter.Grp)

        VarG_All[[l]] = VarG
        #} # End of corv
      }

      names(VarCov_All) = G_list

      names(CORR_ALL) = G_list
    }

    ### End corv


  }



  ### Calculate genomic Heritability for Each Location for any of the variance-covariance structure
  VarE = vc[grep("!R", rownames(vc)), "component"]

  names(VarE) = Heter.Grp

  #varG_matrix = matrix(NA,  nrow = length(G_list), ncol = length(Heter.Grp))

  # varG_matrix = matrix(NA,  1, ncol = length(Heter.Grp))
  #
  # colnames(Heter.Grp) <-  Heter.Grp
  #
  # #dimnames(varG_matrix) <- list(unlist(G_list), Heter.Grp)
  #
  #
  # for (h in 1:length(G_list)) {
  #
  #   varG_matrix[h, ]  = diag(VarCov_All[[h]])
  #
  # }

  #VarG = diag(VarCov)
  H = matrix(NA, nrow = 1, ncol = length(Heter.Grp))

  colnames(H) = Heter.Grp

  varG_per_omics = matrix(NA, nrow = length(G_list), ncol = length(Heter.Grp))

  dimnames(varG_per_omics) <- list(unlist(G_list), Heter.Grp)

  Total_varG = matrix(NA, nrow = 1, ncol = length(Heter.Grp))

  colnames(Total_varG) = colnames(Heter.Grp)



  for (va in 1:length(G_list)) {

    #varG_per_omics[va, ] <- VarG_All[[va]]

    varG_per_omics[va, ] <- diag(VarCov_All[[va]])


  }





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

  ### This is to create variance-covariance and correlation matrix
  ## for every any number of G that the user provide
 if(length(G_list)>1) {

   ### To get things setup assign the element of
   ## G_list to global enviornment

   for (va in 1:length(G_list)) {

     assign(G_list[[va]], G_list[[va]])


   }

   #### Start the looping
  for (va in 1:length(G_list)) {

  assign(paste(G_list[[va]], "varcov", sep = "_"),  VarCov_All[[va]])

  assign(paste(G_list[[va]], "corr", sep = "_"),  CORR_ALL[[va]])

  }

   if (exists("G") & (exists("omic1") | exists("omic2") | exists("omic3"))){
     if(exists("omic1")){
       omic_varCor = omic1_varcov
       omic_corr = omic1_corr

     } else if (exists("omic2")){

       omic_varCor = omic2_varcov
       omic_corr = omic2_corr
     } else {

       if(exists("omic3")){
         omic_varCor = omic3_varcov
         omic_corr = omic3_corr
       }
     }

   Result = list(Covariance1 = G_varcov,
                 Covariance2 = omic_varCor,
                 Correlation1 = G_corr,
                 Correlation2 = omic_corr,
                 Heritability = H,
                 Total_genetic_var = Total_varG,
                 varG_per_omics = varG_per_omics,
                 Residual_Var = VarE)

   names(Result) <-  c("Covariance1",
                       "Covariance2",
                       "Correlation1",
                       "Correlation2",
                       "Heritability",
                       "Total_genetic_var",
                       "varG_per_omics",
                       "Residual_Var")


   } else if (exists("G") & (exists("omic1") & exists("omic2"))){

     if(exists("omic1")){
       omic1_varCor = omic1_varcov
       omic1_corr = omic1_corr
     }


       if(exists("omic2")){
         omic2_varCor = omic2_varcov
         omic2_corr = omic2_corr
       }


     Result = list(Covariance1 = G_varcov,
                   Covariance2 = omic1_varCor,
                   Covariance3 = omic2_varCor,
                   Correlation1 = G_corr,
                   Correlation2 = omic1_corr,
                   Correlation3 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

        names(Result) <- c("Covariance1",
                           "Covariance2",
                           "Covariance3",
                           "Correlation1",
                           "Correlation2",
                           "Correlation3",
                           "Heritability",
                           "Total_genetic_var",
                           "varG_per_omics",
                           "Residual_Var")


   } else if (exists("G") & (exists("omic1") & exists("omic3"))){

     if(exists("omic1")){
       omic1_varCor = omic1_varcov
       omic1_corr = omic1_corr
     }


     if(exists("omic3")){
       omic2_varCor = omic3_varcov
       omic2_corr = omic3_corr
     }


     Result = list(Covariance1 = G_varcov,
                   Covariance2 = omic1_varCor,
                   Covariance3 = omic2_varCor,
                   Correlation1 = G_corr,
                   Correlation2 = omic1_corr,
                   Correlation3 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

     names(Result) <- c("Covariance1",
                       "Covariance2",
                       "Covariance3",
                       "Correlation1",
                       "Correlation2",
                       "Correlation3",
                       "Heritability",
                       "Total_genetic_var",
                       "varG_per_omics",
                       "Residual_Var")

   } else if (exists("G") & (exists("omic2") & exists("omic3"))){

     if(exists("omic2")){
       omic1_varCor = omic2_varcov
       omic1_corr = omic2_corr
     }


     if(exists("omic3")){
       omic2_varCor = omic3_varcov
       omic2_corr = omic3_corr
     }


     Result = list(Covariance1 = G_varcov,
                   Covariance2 = omic1_varCor,
                   Covariance3 = omic2_varCor,
                   Correlation1 = G_corr,
                   Correlation2 = omic1_corr,
                   Correlation3 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

     names(Result) <- c("Covariance1",
                        "Covariance2",
                        "Covariance3",
                        "Correlation1",
                        "Correlation2",
                        "Correlation3",
                        "Heritability",
                        "Total_genetic_var",
                        "varG_per_omics",
                        "Residual_Var")

   } else if (exists("omic2") & exists("omic3")){


     if(exists("omic2")){
       omic1_varCor = omic2_varcov
       omic1_corr = omic2_corr
     }


     if(exists("omic3")){
       omic2_varCor = omic3_varcov
       omic2_corr = omic3_corr
     }




     Result = list(Covariance1 = omic1_varCor,
                   Covariance2 = omic2_varCor,
                   Correlation1 = omic1_corr,
                   Correlation2 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

     names(Result) <- c("Covariance1",
                        "Covariance2",
                        "Correlation1",
                        "Correlation2",
                        "Heritability",
                        "Total_genetic_var",
                        "varG_per_omics",
                        "Residual_Var")


   } else if ((exists("omic1") & exists("omic2"))){

     if(exists("omic1")){
       omic1_varCor = omic1_varcov
       omic1_corr = omic1_corr
     }


     if(exists("omic2")){
       omic2_varCor = omic2_varcov
       omic2_corr = omic2_corr
     }


     Result = list(Covariance1 = omic1_varCor,
                   Covariance2 = omic2_varCor,
                   Correlation1 = omic1_corr,
                   Correlation2 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

     names(Result) <- c("Covariance1",
                        "Covariance2",
                        "Correlation1",
                        "Correlation2",
                        "Heritability",
                        "Total_genetic_var",
                        "varG_per_omics",
                        "Residual_Var")


   } else if ((exists("omic1") & exists("omic3"))){

     if(exists("omic1")){
       omic1_varCor = omic1_varcov
       omic1_corr = omic1_corr
     }


     if(exists("omic3")){
       omic2_varCor = omic3_varcov
       omic2_corr = omic3_corr
     }


     Result = list( Covariance1 = omic1_varCor,
                   Covariance2 = omic2_varCor,
                   Correlation1 = omic1_corr,
                   Correlation2 = omic2_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

names(Result) <-  c("Covariance1",
                    "Covariance2",
                    "Correlation1",
                    "Correlation2",
                    "Heritability",
                    "Total_genetic_var",
                    "varG_per_omics",
                    "Residual_Var")

   } else if (exists("omic3") & (exists("omic1") & exists("omic2"))){


     if(exists("omic1")){
       omic1_varCor = omic1_varcov
       omic1_corr = omic1_corr
     }

     if(exists("omic2")){
       omic2_varCor = omic2_varcov
       omic2_corr = omic2_corr
     }


     if(exists("omic3")){
       omic3_varCor = omic3_varcov
       omic3_corr = omic3_corr
     }


     Result = list(Covariance1 = omic1_varCor,
                   Covariance2 = omic2_varCor,
                   Covariance3 = omic3_varCor,
                   Correlation1 = omic1_corr,
                   Correlation2 = omic2_corr,
                   Correlation3 = omic3_corr,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   varG_per_omics = varG_per_omics,
                   Residual_Var = VarE)

  names(Result) <- c("Covariance1",
                   "Covariance2",
                   "Covariance3",
                   "Correlation1",
                   "Correlation2",
                   "Correlation3",
                   "Heritability",
                   "Total_genetic_var",
                   "varG_per_omics",
                   "Residual_Var")


   } else {


     if (exists("G") & (exists("omic1") & (exists("omic2") & exists("omic3")))){


       if(exists("omic1")){
         omic1_varCor = omic1_varcov
         omic1_corr = omic1_corr
       }

       if(exists("omic2")){
         omic2_varCor = omic2_varcov
         omic2_corr = omic2_corr
       }


       if(exists("omic3")){
         omic3_varCor = omic3_varcov
         omic3_corr = omic3_corr
       }


       Result = list( Covariance1 = G_varcov,
                     Covariance2 = omic1_varCor,
                     Covariance3 = omic2_varCor,
                     Covariance4 = omic3_varCor,
                     Correlation1 = G_corr,
                     Correlation2 = omic1_corr,
                     Correlation3 = omic2_corr,
                     Correlation4 = omic2_corr,
                     Heritability = H,
                     Total_genetic_var = Total_varG,
                     varG_per_omics = varG_per_omics,
                     Residual_Var = VarE)

  names(Result) <- c("Covariance1",
                    "Covariance2",
                    "Covariance3",
                    "Covariance4",
                    "Correlation1",
                    "Correlation2",
                    "Correlation3",
                    "Correlation4",
                    "Heritability",
                    "Total_genetic_var",
                    "varG_per_omics",
                    "Residual_Var")
     }

   }


 } else {

   if(length(G_list)==1) {

     VarCov_All =  unlist(VarCov_All)
     CORR_ALL = unlist(CORR_ALL)

     Result = list(Covariance = VarCov_All,
                   Correlation = CORR_ALL,
                   Heritability = H,
                   Total_genetic_var = Total_varG,
                   Residual_Var = VarE)

  names(Result) <- c("Covariance",
                   "Correlation",
                   "Heritability",
                   "Total_genetic_var",
                   "Residual_Var")


   }

 }

  #### Results


  return(Result)

}


#Res= VarCov_Her(model= mod.ref, heter_groups=heter_groups, var_cov_str= var_cov_str)


