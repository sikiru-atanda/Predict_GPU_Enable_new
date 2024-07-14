#' Title
#'
#' @param mod_asreml
#' @param pheno_data
#' @param heter_groups
#' @param gen_name
#' @param var_cov_str
#' @param heter_resid
#' @param gkernel
#' @param gmatrix
#' @param omic1_kernel
#' @param omic2_kernel
#' @param omic3_kernel
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
asreml_mod_output <- function(
    mod_asreml = NULL,
    pheno_data = NULL,
    gkernel=NULL,
    gmatrix = NULL,
    omic1_kernel=NULL,
    omic2_kernel=NULL,
    omic3_kernel=NULL,
    heter_groups = NULL,
    gen_name = NULL,
    var_cov_str = NULL,
    heter_resid = NULL,
    pworkspace= 1e15,
    #workspace = 1e08,
    maxit = 50,
    ...
)
{


  msg <- "\n==================================================\n"

  if (!is.null(gkernel)){
    gmatrix = gkernel
  }
  mod = mod_asreml$model
  #mod = mod[[1]]
  ### if every variance component is stable the update will not run
  ## by default in asreml so it safe to keep it
  mod = asreml::update.asreml(mod)
  str.mod = mod_asreml$str.mod
  gen_pos = mod_asreml$gen_pos
  inter_gen_pos = mod_asreml$inter_gen_pos
  G_list = mod_asreml$G_list
  rand_term = mod_asreml$rand_term
  #############################
  if(!is.null(gmatrix)){
    gid_name <- rownames(gmatrix)
    g_retain <-  gmatrix ## To store back in the database for reuse
  } else if (!is.null(gkernel)){
    gid_name <- rownames(gkernel)
    g_retain <-  gkernel ## To store back in the database for reuse
  } else if (!is.null(omic1_kernel)){
    gid_name <- rownames(omic1_kernel)
  } else if (!is.null(omic2_kernel)){
    gid_name <- rownames(omic2_kernel)
  } else {
    if (!is.null(omic3_kernel)){
      gid_name <- rownames(omic3_kernel)
    }
  }
######
  ### This is important because omic_kernel are compromised when there is more than  one environment
  ## So we need to original matrix to store back.
  if (!is.null(omic1_kernel)){
    omic1_retain <- omic1_kernel
  }

  if (!is.null(omic2_kernel)){
    omic2_retain <- omic2_kernel
  }

  if (!is.null(omic3_kernel)){
    omic3_retain <- omic3_kernel
  }
  ###
  ### !is.null(var_cov_str) & is.null(inter_gen_pos) incase user provide var_cov_str
  ## while the data is not MT in nature
  if(!is.null(var_cov_str) & is.null(inter_gen_pos)){
    var_cov_str = NULL
    heter_groups = NULL
    heter_resid = NULL
  }

  #########################
  ### Set up parameter for predict function
  # asreml::asreml.options(trace=FALSE, workspace=workspace,
  #                        pworkspace = pworkspace, maxit = maxit)

  BLUP <- summary(mod, coef=TRUE)$coef.random

  colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Std_error"

  #heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

  if (!is.null(heter_groups)){
    ### It possible the user provide the heter_groups while it actually a single environment,
    ## This will check and turn it off
    if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
      heter_groups = NULL
    } else{
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
      }
    }

  }
  #ENV_Ids = as.character(unique(pheno_data[[heter_groups]]))
  ##############################################
  #### Extract Breeding values/genetic effect estimate for all omics
  ##############################################
  BV_All = vector(mode = 'list', length = length(G_list))

  names(BV_All) <- unlist(G_list)

  ##
  for (b in 1:length(G_list)) {


    BV_All[[b]] <- BLUP[grep(paste(paste(G_list[[b]], '_inv', sep = ""),"\\)", sep = ""),rownames(BLUP)),]

  }
 ###
  ### For variance structure extraction
  if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

    if(isTRUE(grepl("fa", var_cov_str))){
      ## Extract the number of factors
      #N_fa = substr(var_cov_str, 3, 100)

      for (bb in 1:length(G_list)) {
        BV_All[[bb]] <- BV_All[[bb]][!rownames(BV_All[[bb]])%in%rownames(BV_All[[bb]][grep('Comp',rownames(BV_All[[bb]])),]), ]

        BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

        BV_All[[bb]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,3])

      }

    } else {

      if (var_cov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {


        for (bb in 1:length(G_list)) {

          BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

          BV_All[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])

        }


      }


    } ## End

    #################################
    if(!is.null(inter_gen_pos)){
      for (bb in 1:length(G_list)) {

        BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]

        BV_All[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(BV_All[[bb]][, gen_name])))

        BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]

        colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")

        BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "Std_error"]^2
        rownames( BV_All[[bb]]) = NULL


      }

    } else {

      stop(message(paste(msg, "No interaction term")), call. = FALSE)

      # if(is.null(inter_gen_pos)){
      # for (bb in 1:length(G_list)) {
      #
      #   BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]
      #
      #   BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]
      #
      #   colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
      #
      #   BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "std.error"]^2
      #
      # }
      #
      # }
    }

    #######################################

    Res= asreml_herit_varCov(model= mod,
                             heter_groups= heter_groups,
                             var_cov_str= var_cov_str,
                             heter_resid= heter_resid,
                             G_list = G_list,
                             inter_gen_pos = inter_gen_pos,
                             gen_pos = gen_pos)



    #VA = Res$Genetic_Var
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in 1:length(G_list)){

      if(length(G_list)>1){
        VA = Res$varG_per_omics[[bb]]

      } else {
        if(length(G_list)==1){

          VA = unlist(Res$Total_genetic_var)
        }

      }

      if(length(VA)< length(heter_grp)){

        BV_All[[bb]][, "Reliability"] = NA

        message(paste( insight::print_color("WARNINGS\n", "blue"),
                       insight::print_color(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.")), "blue")))
        #message(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.")))
      } else{

        if(length(VA) == length(heter_grp)){

          BV_All[[bb]][, "Reliability"] <- NA
          #BV$Reliability = NA
          for (i in 1:length(VA)) {

            BV_All[[bb]][, "Reliability"] <- ifelse(BV_All[[bb]][, heter_groups]%in% heter_grp[i],
                                                    round(1 - BV_All[[bb]][, "PEV"]/VA[i],6), BV_All[[bb]][, "Reliability"])

          }
        }


      }

    }


 } else {
    ## Problem
    if(is.null(var_cov_str) & is.null(inter_gen_pos) ){


      for (bb in 1:length(G_list)){
        BV_All[[bb]] <- as.data.frame(BV_All[[bb]])
        BV_All[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])
        rownames( BV_All[[bb]]) = NULL
        BV_All[[bb]] =  BV_All[[bb]][, c(4, 1:2)]
        colnames(BV_All[[bb]])[1:2] <- c(gen_name, "BLUP")
      }

      vc <- summary(mod)$varcomp
      #VAR_check_Pos <- which(vc$bound=='F' | vc$bound=='U' | vc$bound=="?" | vc$bound=="S")
      VAR_check_Pos <- which(vc$bound=="?" | vc$bound=="S")
      if(length(VAR_check_Pos)>1) {


        # stop(message(paste( insight::print_color("STOP\n", "red"),
        #                insight::print_color(paste(msg, paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
        #                                                      " is unstable, refix the model")), "red"))), call. = FALSE)
        #
        stop(print(paste(msg, paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
                                    " is unstable, refix the model"))), call. = FALSE)



      } else {

        if(length(VAR_check_Pos)==1) {

          stop(print(paste(msg, paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" " ),
                                      " is unstable, refix the model"))), call. = FALSE)



        }

      }

      VarG_All <- vector("list", length = length(G_list))
      names(VarG_All) <- as.character(unlist(G_list))
      VarG <- vc[grep(paste("^vm\\(", gen_name, sep = ""), rownames(vc)), drop=FALSE, ]

      for (g in 1:length(G_list)) {

        VarG_All[[g]] <-VarG[grep(paste(G_list[[g]],"_inv", sep = ""), rownames(VarG)), "component"]

      }

      VE <- vc[grep("!R", row.names(vc)), "component"]

      varG_matrix = unlist(VarG_All)

      H = matrix(NA, nrow = 1, ncol = length(VE))




      #BV$PEV <- BV$std.error^2
      #BV_All$PEV <- BV_All$std.error^2
      for (bb in 1:length(G_list)){

        BV_All[[bb]][, "PEV"] <- BV_All[[bb]][, "Std_error"]^2

        BV_All[[bb]][, "Reliability"] <- round(1 - BV_All[[bb]][, "PEV"]/as.double((VarG_All[bb])),6)

        if(length(VE)>1){

          for (i in 1:length(VE)) {
            H[, i] <- sum(unlist(VarG_All))/(sum(unlist(VarG_All))+VE[i])

          }

        } else {

          if (length(VE)==1){

            #BV$Reliability <- round(1 - BV$PEV/VA,6)

            H <- sum(unlist(VarG_All))/(sum(unlist(VarG_All))+VE)
          }
        }

      }

      VarG_All <-  as.matrix(unlist(VarG_All))
      colnames(VarG_All) <- "Variance"
      Total_genetic_var <- sum(VarG_All[, 1])
      #VarG_All <- cbind(VarG_All, Standard_error = NA)
      Res_Va_Ve_H2_COV_COR = list(Heritability = H,
                                  varG_per_omics = VarG_All,
                                  Total_genetic_var = Total_genetic_var,
                                  Residual_Var = VE)
    } ## End


  } #### End var_Covar

  ## when variance_covariance structure is not defined by the user and CS is used
  ### For variance structure extraction
  if(is.null(var_cov_str) & !is.null(inter_gen_pos)){


    for (bb in 1:length(G_list)) {

      BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

      BV_All[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])

      rownames( BV_All[[bb]]) = NULL
    }



    #################################
    if(!is.null(inter_gen_pos)){
      for (bb in 1:length(G_list)) {

        BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]

        BV_All[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(BV_All[[bb]][, gen_name])))

        BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]

        colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")

        BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "Std_error"]^2
        rownames( BV_All[[bb]]) = NULL

      }

    } else {

      stop(print(paste(msg, 'No interaction term')), call. = FALSE)

    }

    #######################################

    Res = asreml_herit_CSM(model= mod,
                           heter_groups= heter_groups,
                           heter_resid= heter_resid,
                           G_list = G_list,
                           inter_gen_pos = inter_gen_pos,
                           gen_pos = gen_pos)


    VE <-  Res$Residual_Var
    H <- Res$Heritability
    varG_matrix = Res$varG_per_omics
    #VA = Res$Genetic_Var
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in 1:length(G_list)){

      if(length(G_list)>1){
        VA = Res$varG_per_omics[bb, ]

      } else {
        if(length(G_list)==1){

          VA = Res$Total_genetic_var[1, ]
        }

      }

      if(length(VA)< length(heter_grp)){

        print(paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive"))
      } else{

        if(length(VA) == length(heter_grp)){

          BV_All[[bb]][, "Reliability"] = NA
          #BV$Reliability = NA
          for (i in 1:length(VA)) {

            BV_All[[bb]][, "Reliability"] <- ifelse(BV_All[[bb]][, heter_groups]%in% heter_grp[i],
                                                    round(1 - BV_All[[bb]][, "PEV"]/VA[i],6), BV_All[[bb]][, "Reliability"])

          }
        } else {

          stop(print(paste(msg, 'Genetic variance is missing')), call. = FALSE)

        }


      }

    }


  } ## End CS


  ### Predicted Values
  pred_value <- asreml::predict.asreml(mod, classify=gen_name, sed=FALSE)$pvals
  pred_value = pred_value[, -ncol(pred_value)] ### Remove status
  colnames(pred_value)[colnames(pred_value)%in%c("predicted.value", "std.error")] <- c("Predicted_value", "Std_error")
  pred_value[, "PEV"] <-  pred_value[, "Std_error"]^2

  pred_value[, "Reliability"] <-  NA

  if(!is.null(heter_groups)){
    if(var(pred_value$Predicted_value)==0){
      ### Check if the the across
      message(paste( insight::print_color("WARNING\n", "blue"),
                     insight::print_color(paste(msg, paste(paste('The average prediction across', heter_groups), paste('is a constant value.\n \t Check the model to change', heter_groups), 'to fixed term ')), "blue")))

    }
  }
  gc()
  # if (is.null(heter_groups) & is.null(var_cov_str)){inter_gen_pos= NULL}
  # if(length(gen_pos) == length(rand_term)){inter_gen_pos= NULL}
  # if(!is.null(inter_gen_pos)){
  #
  #   pred_heter_groups <- asreml::predict.asreml(mod, classify= rand_term[[inter_gen_pos]], sed=FALSE)$pvals
  #   pred_heter_groups =  pred_heter_groups[, -ncol(pred_heter_groups)] ### Remove status
  #   colnames(pred_heter_groups)[colnames(pred_heter_groups)%in%c("predicted.value", "std.error")] <- c("Predicted_value", "Std_error")
  #   pred_heter_groups[, "PEV"] <-  pred_heter_groups[, "Std_error"]^2
  #   pred_heter_groups[, "Reliability"] <- NA
  #
  #   name_across_env <-  paste("Across", paste(heter_groups, "Predicted_value", sep = "_"), sep = "_")
  #
  # }
  ### Results
  ############################# Coffient cal

  ### Unlist each breeding value for multi-omics
  # if(length(G_list)>1){
  #   for (b in 1:length(G_list)) {
  #
  #     assign(paste("ebv", G_list[[b]], sep = "_"),  BV_All[[b]])
  #
  #   }
  #
  # } else {
  #   if(length(G_list)==1){
  #     BV_All = BV_All[[1]]
  #   }
  #
  # }

  ### To get things setup assign the element of
  ## G_list to global enviornment

  # for (va in 1:length(G_list)) {
  #
  #   assign(G_list[[va]], G_list[[va]])
  #
  #
  # }
  ### Initialize step to calculate coefficient for each omics

  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    if(!is.null(gmatrix)){
      gmatrix <- Zg%*%gmatrix%*%t(Zg)

      gmatrix <- grm_kernel_precheck(gmatrix)
    }

    if(!is.null(gkernel)){
      gmatrix <- Zg%*%gkernel%*%t(Zg)

      gmatrix <- grm_kernel_precheck(gmatrix)
    }

    if(!is.null(omic1_kernel)){
      omic1_kernel <- Zg%*%omic1_kernel%*%t(Zg)

      omic1_kernel <- grm_kernel_precheck(omic1_kernel)

    }

    if(!is.null(omic2_kernel)){
      omic2_kernel <- Zg%*%omic2_kernel%*%t(Zg)

      omic2_kernel <- grm_kernel_precheck(omic2_kernel)
    }

    if(!is.null(omic3_kernel)){
      omic3_kernel <- Zg%*%omic3_kernel%*%t(Zg)

      omic3_kernel <- grm_kernel_precheck(omic3_kernel)

    }

  }
  #}
  if(length(G_list)>1){

    if(length(intersect(c("G", "omic1"),unlist(G_list)))==length(G_list)){

      #### Check for the Generalized and unique Inverse
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }
        ##
        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }
      }

    } else if(length(intersect(c("G", "omic2"),unlist(G_list)))==length(G_list)){

      #### Check for the Generalized and unique Inverse
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }
        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
      }


    } else if(length(intersect(c("G", "omic3"),unlist(G_list)))==length(G_list)){

      #### Check for the Generalized and unique Inverse
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }
        ##
        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }

    } else if(length(intersect(c("omic1", "omic2"),unlist(G_list)))==length(G_list)){

      #### Check for the Generalized and unique Inverse
      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }
        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
      }

    } else if(length(intersect(c("omic1", "omic3"),unlist(G_list)))==length(G_list)){

      ##
      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }
        ##
        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }


    } else if(length(intersect(c("omic2", "omic3"),unlist(G_list)))==length(G_list)){

     #####
      for (n in 1:length(G_list)) {

        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
        ##
        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }

    } else if(length(intersect(c("G","omic1", "omic2"),unlist(G_list)))==length(G_list)){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }

        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }
        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
      }
      ##

    } else if(length(intersect(c("G","omic1", "omic3"),unlist(G_list)))==length(G_list)){

      # ebv1 = ebv_G
      # ebv2 = ebv_omic1
      # ebv3 = ebv_omic3
      ###
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }

        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }
        ##
        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }

    } else if(length(intersect(c("G","omic2", "omic3"),unlist(G_list)))==length(G_list)){

      # ebv1 = ebv_G
      # ebv2 = ebv_omic2
      # ebv3 = ebv_omic3

      ###
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
          coeffRaw_G <- colMeans(coeffRaw_G)
        }

        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
        ##

        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }


    } else if(length(intersect(c("omic1","omic2", "omic3"),unlist(G_list)))==length(G_list)){

      ###
      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }

        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
        ##

        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }

    } else {

      if(length(intersect(c("G","omic1", "omic2", "omic3"),unlist(G_list)))==length(G_list)){

        ###
        for (n in 1:length(G_list)) {

          if("G"==names(BV_All)[n]){
            coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]][, "BLUP"]
            coeffRaw_G <- colMeans(coeffRaw_G)
          }


          if("omic1"==names(BV_All)[n]){
            coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]][, "BLUP"]
            coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
          }

          ##
          if("omic2"==names(BV_All)[n]){
            coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]][, "BLUP"]
            coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
          }
          ##

          if("omic3"==names(BV_All)[n]){
            coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]][, "BLUP"]
            coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
          }
        }

      }

    }

  } else{

    if(!is.null(gmatrix)){


      coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[1]][, "BLUP"]
      coeffRaw_G <- colMeans(coeffRaw_G)

    } else if(!is.null(gkernel)){


      coeffRaw_G = solve(t(gkernel)*gkernel)*t(gkernel)*BV_All[[1]][, "BLUP"]
      coeffRaw_G <- colMeans(coeffRaw_G)

    } else if(!is.null(omic1_kernel)){

      coeffRaw_G = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[1]][, "BLUP"]
      coeffRaw_G <- colMeans(coeffRaw_G)

    } else if(!is.null(omic2_kernel)){

      coeffRaw_G = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[1]][, "BLUP"]
      coeffRaw_G <- colMeans(coeffRaw_G)

    } else if(!is.null(omic3_kernel)){

      coeffRaw_G = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[1]][, "BLUP"]
      coeffRaw_G <- colMeans(coeffRaw_G)
    }

  }



  if(!is.null(heter_groups)){

    if(exists("coeffRaw_G")){
      coeffRaw_G <- data.frame(x_variables = rep(gid_name , length(unique(heter_grp))),
                        Env = rep(unique(heter_grp), each= length(gid_name)),
                        coeff = coeffRaw_G,
                        stringsAsFactors = FALSE)
    names(coeffRaw_G)[2] <- heter_groups

    }
    ###
    if(exists("coeffRaw_omic1")){
      coeffRaw_omic1 <- data.frame(x_variables = rep(gid_name , length(unique(heter_grp))),
                               Env = rep(unique(heter_grp), each= length(gid_name)),
                               coeff = coeffRaw_omic1,
                               stringsAsFactors = FALSE)
      names(coeffRaw_omic1)[2] <- heter_groups

    }
    ####
    if(exists("coeffRaw_omic2")){
      coeffRaw_omic2 <- data.frame(x_variables = rep(gid_name , length(unique(heter_grp))),
                                   Env = rep(unique(heter_grp), each= length(gid_name)),
                                   coeff = coeffRaw_omic2,
                                   stringsAsFactors = FALSE)
      names(coeffRaw_omic2)[2] <- heter_groups

    }
    ###
    if(exists("coeffRaw_omic3")){
      coeffRaw_omic3 <- data.frame(x_variables = rep(gid_name , length(unique(heter_grp))),
                                   Env = rep(unique(heter_grp), each= length(gid_name)),
                                   coeff = coeffRaw_omic3,
                                   stringsAsFactors = FALSE)
      names(coeffRaw_omic3)[2] <- heter_groups

    }
  } else {

    if(exists("coeffRaw_G")){
    coeffRaw_G <- data.frame(x_variables = gid_name,
                        coeff = coeffRaw_G,
                        stringsAsFactors = FALSE)

    }

    if(exists("coeffRaw_omic1")){
      coeffRaw_omic1 <- data.frame(x_variables = gid_name,
                               coeff = coeffRaw_omic1,
                               stringsAsFactors = FALSE)

    }

    if(exists("coeffRaw_omic2")){
      coeffRaw_omic2 <- data.frame(x_variables = gid_name,
                                   coeff = coeffRaw_omic2,
                                   stringsAsFactors = FALSE)

    }

    if(exists("coeffRaw_omic3")){
      coeffRaw_omic3 <- data.frame(x_variables = gid_name,
                                   coeff = coeffRaw_omic3,
                                   stringsAsFactors = FALSE)

    }

  }

  ##### End

    ####
    if(!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      #     pred_value[, "Reliability"] <-   round(1 - pred_value[, "PEV"]/Res$Total_genetic_var,6)
      #     pred_value[, "Reliability"] = ifelse( pred_value[, "Reliability"]<0, "Alias",  pred_value[, "Reliability"])
      # ####
      #     pred_heter_groups[, "Reliability"] <-   round(1 - pred_heter_groups[, "PEV"]/Res$Total_genetic_var,6)
      #     pred_heter_groups[, "Reliability"] = ifelse(pred_heter_groups[, "Reliability"]<0, "Alias",  pred_heter_groups[, "Reliability"])

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = BV_All[[1]],
                    pred_value =pred_value,
                    Coefficients = coeffRaw_G,
                    list(Geno_model_ready = g_retain)
                    )

      rm(g_retain)

    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){


      Result <- list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = BV_All[[1]],
                    pred_value =pred_value,
                    Coefficients = coeffRaw_G,
                    list(Omic_model_ready = omic1_retain))


      rm(omic1_retain, omic1_kernel)

    } else if (is.null(gmatrix)  & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){


      Result = list(call = str.mod,
                    mod = mod,
                    Estimated_breeding_value = BV_All[[1]],
                    pred_value =pred_value,
                    Coefficients = coeffRaw_G,
                    list(Omic_model_ready = omic2_retain))


      rm(omic2_retain, omic2_kernel)

    } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value=BV_All[[1]],
                    pred_value =pred_value,
                    Coefficients = coeffRaw_G,
                    list(Omic_model_ready = omic3_retain))


      rm(omic3_retain, omic3_kernel)

    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

      }

      #pred_value[, "Reliability"] <- round(1 - pred_value[, "PEV"]/VA[i],6)


      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1= ebv_G,
                               Estimated_breeding_value_2 = ebv_omic1),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 =coeffRaw_G,
                                        Coefficients_2  = coeffRaw_omic1),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic_model_ready = omic1_retain
                    ))


      rm(ebv_G, ebv_omic1, g_retain, omic1_retain, gmatrix, omic1_kernel)

    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(EBV_1 = ebv_G,
                                 EBV_2 = ebv_omic2),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_G,
                                        Coefficients_2 = coeffRaw_omic2),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic_model_ready = omic2_retain
                    ))


      rm(ebv_G, ebv_omic2, g_retain, omic2_retain, gmatrix, omic2_kernel)

    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_G,
                                                    Estimated_breeding_value_2 = ebv_omic3),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_G,
                                        Coefficients_2 = coeffRaw_omic3),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic_model_ready = omic3_retain
                    ))



      rm(ebv_G, ebv_omic3, g_retain, omic3_retain, gmatrix, omic3_kernel)
    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_G,
                                                    Estimated_breeding_value_2 = ebv_omic1,
                                                    Estimated_breeding_value_3 = ebv_omic2),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_G,
                                        Coefficients_2 = coeffRaw_omic1,
                                        Coefficients_3 = coeffRaw_omic2),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic1_kernel_model_ready = omic1_retain,
                                      Omic2_kernel_model_ready = omic2_retain
                    ))


      rm(ebv_G, ebv_omic1, ebv_omic2,
         g_retain, omic1_retain, omic2_retain,
         gmatrix, omic1_kernel, omic2_kernel)

    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_G,
                                                    Estimated_breeding_value_2 = ebv_omic1,
                                                    Estimated_breeding_value_3 = ebv_omic3),
                    pred_value =pred_value,

                    Coefficients = list(Coefficients_1 = coeffRaw_G,
                                        Coefficients_2 = coeffRaw_omic1,
                                        Coefficients_3 = coeffRaw_omic3),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic1_kernel_model_ready = omic1_retain,
                                      Omic2_kernel_model_ready = omic3_retain
                    ))


      rm(ebv_G, ebv_omic1, ebv_omic3,
         g_retain, omic1_retain, omic3_retain,
         gmatrix, omic1_kernel, omic3_kernel)

    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_G,
                                                    Estimated_breeding_value_2 = ebv_omic2,
                                                    Estimated_breeding_value_3 = ebv_omic3),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_G,
                                        Coefficients_2 = coeffRaw_omic2,
                                        Coefficients_3 = coeffRaw_omic3),
                    omics_geno = list(Geno_model_ready = g_retain,
                                      Omic1_kernel_model_ready = omic2_retain,
                                      Omic2_kernel_model_ready = omic3_retain
                    ))

      rm(ebv_G, ebv_omic3, ebv_omic2,
         g_retain, omic3_retain, omic2_retain,
         gmatrix, omic3_kernel, omic2_kernel)

    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

      }
      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_omic1,
                                                    Estimated_breeding_value_2 = ebv_omic2),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_omic1,
                                        Coefficients_2 = coeffRaw_omic2),
                    omics_omics = list(Omic1_kernel_model_ready = omic1_retain,
                                       Omic2_kernel_model_ready = omic2_retain
                    ))

      rm(ebv_omic1, ebv_omic2,
         omic1_retain, omic2_retain,
          omic1_kernel, omic2_kernel)
    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_omic1,
                                                    Estimated_breeding_value_2 = ebv_omic3),
                    pred_value =pred_value,
                    coefficients = list(Coefficients_1 = coeffRaw_omic1,
                                        Coefficients_2 = coeffRaw_omic3),
                    omics_omics = list(Omic1_kernel_model_ready = omic1_retain,
                                       Omic2_kernel_model_ready = omic3_retain
                    ))

      rm(ebv_omic1, ebv_omic3,
         omic1_retain, omic3_retain,
         omic1_kernel, omic3_kernel)

    } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_omic2,
                               Estimated_breeding_value_2 = ebv_omic3),
                    pred_value =pred_value,
                    coefficients = list(Coefficients_1 = coeffRaw_omic2,
                                        Coefficients_2 = coeffRaw_omic3),
                    omics_omics = list(Omic1_kernel_model_ready = omic2_retain,
                                       Omic2_kernel_model_ready = omic3_retain
                    ))

      rm(ebv_omic3, ebv_omic2,
         omic3_retain, omic2_retain,
         omic3_kernel, omic2_kernel)

    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_omic1,
                                                    Estimated_breeding_value_2 = ebv_omic2,
                                                    Estimated_breeding_value_3 = ebv_omic3),
                    pred_value =pred_value,
                    Coefficients = list(Coefficients_1 = coeffRaw_omic1,
                                        Coefficients_2 = coeffRaw_omic2,
                                        Coefficients_3 = coeffRaw_omic3),
                    omics_omics = list(Omic1_kernel_model_ready = omic1_retain,
                                       Omic2_kernel_model_ready = omic2_retain,
                                       Omic3_kernel_model_ready = omic3_retain
                    ))


      rm(ebv_omic3, ebv_omic2,ebv_omic1,
         omic3_retain, omic2_retain,omic1_retain,
         omic3_kernel, omic2_kernel, omic1_kernel)
    } else {

      if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

        for (n in 1:length(G_list)) {

          if("G"==names(BV_All)[n]){
            ebv_G = BV_All[[n]]
          }

          if("omic1"==names(BV_All)[n]){
            ebv_omic1 = BV_All[[n]]
          }

          if("omic2"==names(BV_All)[n]){
            ebv_omic2 = BV_All[[n]]
          }

          if("omic3"==names(BV_All)[n]){
            ebv_omic3 = BV_All[[n]]
          }

        }

        Result = list(call=str.mod,
                      mod=mod,
                      Estimated_breeding_value = list(Estimated_breeding_value_1 = ebv_G,
                                                      Estimated_breeding_value_2 = ebv_omic1,
                                                      Estimated_breeding_value_3 = ebv_omic2,
                                                      Estimated_breeding_value_4 = ebv_omic3),
                      pred_value =pred_value,
                      coefficients = list(Coefficients_1 = coeffRaw_G,
                                          Coefficients_2 = coeffRaw_omic1,
                                          Coefficients_3 = coeffRaw_omic2,
                                          Coefficients_4 = coeffRaw_omic3),
                      omics_geno =   list(Geno_model_ready = g_retain,
                                          Omic1_kernel_model_ready = omic1_retain,
                                          Omic2_kernel_model_ready = omic2_retain,
                                          Omic3_kernel_model_ready = omic3_retain
                      ))



        rm(ebv_G, ebv_omic3, ebv_omic2,ebv_omic1,
           g_retain, omic3_retain, omic2_retain,omic1_retain,
           gmatrix, omic3_kernel, omic2_kernel, omic1_kernel)

      }
    }

if(length(Result)==0){stop("Result miss")}
  names(Result) <- c(
    "Model_structure",
    "Asreml_model",
    "Estimated_breeding_value",
    "Predicted_value",
    "Coefficients",
    "M_matrix_model_ready"

  )

    if(exists("pred_heter_groups")){
      ### Add the across prediction after the Predicted_value for each environment
      aft <-  which(names(Result)=="Predicted_value")
      Result <- append(Result, list(pred_heter_groups = pred_heter_groups), after = aft)
      NN <-  paste("Across", paste(heter_groups, "Predicted_value", sep = "_"), sep = "_")
      names(Result)[(aft+1)] <- NN
    }


  if(!exists("Res_Va_Ve_H2_COV_COR")){
    ### Add the across prediction after the Predicted_value for each environment
    aft <-  which(names(Result)=="Coefficients")

    MM <- rbind(Res$varG_per_omics,
               total_genetic_variance = as.numeric(Res$Total_genetic_var),
               heritability = as.numeric(Res$Heritability),
               residual_variance = Res$Residual_Var)
    #View(MM)
    #MM <- cbind(MM, Standard_error = NA)
    #colnames(MM) <- c("Components", "Standard_error")

    Naming <- unlist(G_list)
    Naming <- gsub("G", "geno", Naming)
    Naming <- paste(Naming, "variance", sep = "_")
    rownames(MM)[1:length(G_list)] <- Naming
    Result <- append(Result, list(Res_Va_Ve_H2_COV_COR = MM), after = aft)
    names(Result)[(aft+1)] <-  "Variance_components"
    if ("Covariance" %in% names(Res) && "Correlation" %in% names(Res)) {
      Result <- append(Result, list(Covariance = Res$Covariance))
      Result <- append(Result, list(Correlation = Res$Correlation))
    }
  } else {
    aft <-  which(names(Result)=="Coefficients")

    MM <- rbind(Res_Va_Ve_H2_COV_COR$varG_per_omics,
               total_genetic_variance = Res_Va_Ve_H2_COV_COR$Total_genetic_var,
               heritability = Res_Va_Ve_H2_COV_COR$Heritability,
               residual_variance = Res_Va_Ve_H2_COV_COR$Residual_Var)

    MM <- cbind(MM, Standard_error = NA)
    colnames(MM) <- c("Components", "Standard_error")

    Naming <- unlist(G_list)
    Naming <- gsub("G", "geno", Naming)
    Naming <- paste(Naming, "variance", sep = "_")
    rownames(MM)[1:length(G_list)] <- Naming

    Result <- append(Result, list(Res_Va_Ve_H2_COV_COR = MM), after = aft)
    names(Result)[(aft+1)] <-  "Variance_components"

  }
  # if(exists("G_inv")) rm(G_inv)
  # if(exists("omic1_inv")) rm(omic1_inv)
  # if(exists("omic2_inv")) rm(omic2_inv)
  # if(exists("omic3_inv")) rm(omic3_inv)
  # ifelse(names(Result) %in% c("Predicted_value", "Variance_components"),
  #        {
  #          Result$Predicted_value <- as.data.frame(Result$Predicted_value)
  #          Result$Variance_components <- as.data.frame(Result$Variance_components)
  #        },
  #        NULL)



  return(Result)
}
