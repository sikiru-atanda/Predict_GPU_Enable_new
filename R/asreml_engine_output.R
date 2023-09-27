#IND = paste(paste0('vm(', gen_name), collapse = ',', 'G_inv)')
# Initiate Process to extract BLUPs/BV from the model

#' Title
#'
#' @param mod_asreml
#' @param pheno_data
#' @param heter_groups
#' @param gen_name
#' @param VarCov_str
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
         VarCov_str = NULL,
         heter_resid = NULL,
         pworkspace= 1e15,
         #workspace = 1e08,
         maxit = 50,
         ...
         )
{

  msg <- sprintf("==================================================\n")

  if (!is.null(gkernel)){
    gmatrix = gkernel
  }
  mod = mod_asreml$model
  mod = mod[[1]]
  ### if every variance component is stable the update will not run
  ## by default in asreml so it safe to keep it
  mod = asreml::update.asreml(mod)
  str.mod = mod_asreml$str.mod
  Gen_pos = mod_asreml$Gen_pos
  Inter_Gen_pos = mod_asreml$Inter_Gen_pos
  G_list = mod_asreml$G_list
  rand_term = mod_asreml$rand_term

  ### Set up parameter for predict function
  # asreml::asreml.options(trace=FALSE, workspace=workspace,
  #                        pworkspace = pworkspace, maxit = maxit)

  BLUP <- summary(mod, coef=TRUE)$coef.random

#Heter.Grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

if (!is.null(heter_groups)){
  Heter.Grp <- as.character(unique(pheno_data[, heter_groups]))
}
#ENV_Ids = as.character(unique(pheno_data[, heter_groups]))
#### Extract Breeding values/genetic effect estimate for all omics
BV_All = vector(mode = 'list', length = length(G_list))

names(BV_All) <- unlist(G_list)

##
for (b in 1:length(G_list)) {


    BV_All[[b]] <- BLUP[grep(paste(paste(G_list[[b]], '_inv', sep = ""),"\\)", sep = ""),rownames(BLUP)),]

}

### For variance structure extraction
if(!is.null(VarCov_str) & !is.null(Inter_Gen_pos)){

  if(isTRUE(grepl("fa", VarCov_str))){
    ## Extract the number of factors
    #N_fa = substr(VarCov_str, 3, 100)

    for (bb in 1:length(G_list)) {
      BV_All[[bb]] <- BV_All[[bb]][!rownames(BV_All[[bb]])%in%rownames(BV_All[[bb]][grep('Comp',rownames(BV_All[[bb]])),]), ]

      BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

      BV_All[[bb]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,3])

    }

  } else {

    if(VarCov_str=="us" |VarCov_str=="corgh" | VarCov_str=="corgv" | VarCov_str=="corh" | VarCov_str=="corv") {


      for (bb in 1:length(G_list)) {

        BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

        BV_All[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])

      }


    }


  } ## End

  #################################
  if(!is.null(Inter_Gen_pos)){
  for (bb in 1:length(G_list)) {

    BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]

    BV_All[[bb]][, heter_groups] <- rep(Heter.Grp, each=length(unique(BV_All[[bb]][, gen_name])))

    BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]

    colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")

    BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "std.error"]^2

  }

  } else {

    stop(print(paste(msg, "No interaction term")), call. = FALSE)

    # if(is.null(Inter_Gen_pos)){
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
                   VarCov_str= VarCov_str,
                   heter_resid= heter_resid,
                   G_list = G_list,
                   Inter_Gen_pos = Inter_Gen_pos,
                   Gen_pos = Gen_pos)



  #VA = Res$Genetic_Var
  Heter.Grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

  for (bb in 1:length(G_list)){

    if(length(G_list)>1){
    VA = Res$varG_per_omics[[bb]]

    } else {
      if(length(G_list)==1){

        VA = unlist(Res$Total_genetic_var)
      }

    }

    if(length(VA)< length(Heter.Grp)){

      message(paste( insight::print_color("WARNINGS\n", "blue"),
                     insight::print_color(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.")), "blue")))
      #message(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.")))
    } else{

      if(length(VA) == length(Heter.Grp)){

        BV_All[[bb]][, "Reliability"] = NA
        #BV$Reliability = NA
        for (i in 1:length(VA)) {

          BV_All[[bb]][, "Reliability"] <- ifelse(BV_All[[bb]][, heter_groups]%in% Heter.Grp[i],
                                                  round(1 - BV_All[[bb]][, "PEV"]/VA[i],6), BV_All[[bb]][, "Reliability"])

        }
      }


    }

  }


} else {
  ## Problem
  if(is.null(VarCov_str) & is.null(Inter_Gen_pos)){


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

    VarG_All = vector("list", length = length(G_list))
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

      BV_All[[bb]][, "PEV"] <- BV_All[[bb]][, "std.error"]^2

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

  } ## End

} #### End var_Covar

## when variance_covariance structure is not defined by the user and CS is used
### For variance structure extraction
if(is.null(VarCov_str) & !is.null(Inter_Gen_pos)){


  for (bb in 1:length(G_list)) {

    BV_All[[bb]] <- as.data.frame(BV_All[[bb]])

    BV_All[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(BV_All[[bb]]), "\\)_", 3)[,2])

    rownames( BV_All[[bb]]) = NULL
  }



  #################################
  if(!is.null(Inter_Gen_pos)){
    for (bb in 1:length(G_list)) {

      BV_All[[bb]] <- BV_All[[bb]][, c(4, 1:2)]

      BV_All[[bb]][, heter_groups] <- rep(Heter.Grp, each=length(unique(BV_All[[bb]][, gen_name])))

      BV_All[[bb]] <- BV_All[[bb]][, c(1, 4, 2:3)]

      colnames(BV_All[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")

      BV_All[[bb]][, "PEV"] <-  BV_All[[bb]][, "std.error"]^2

    }

  } else {

    stop(print(paste(msg, 'No interaction term')), call. = FALSE)

  }

  #######################################

  Res = asreml_herit_CSM(model= mod,
                         heter_groups= heter_groups,
                         heter_resid= heter_resid,
                         G_list = G_list,
                         Inter_Gen_pos = Inter_Gen_pos,
                         Gen_pos = Gen_pos)


  VE <-  Res$Residual_Var
  H <- Res$Heritability
  varG_matrix = Res$varG_per_omics
  #VA = Res$Genetic_Var
  Heter.Grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

  for (bb in 1:length(G_list)){

    if(length(G_list)>1){
      VA = Res$varG_per_omics[bb, ]

    } else {
      if(length(G_list)==1){

        VA = Res$Total_genetic_var[1, ]
      }

    }

    if(length(VA)< length(Heter.Grp)){

      print(paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive"))
    } else{

      if(length(VA) == length(Heter.Grp)){

        BV_All[[bb]][, "Reliability"] = NA
        #BV$Reliability = NA
        for (i in 1:length(VA)) {

          BV_All[[bb]][, "Reliability"] <- ifelse(BV_All[[bb]][, heter_groups]%in% Heter.Grp[i],
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

if(var(pred_value$predicted.value)==0){
  ### Check if the the across

  message(paste( insight::print_color("WARNING\n", "blue"),
                 insight::print_color(paste(msg, paste(paste('The average prediction across', heter_groups), paste('is a constant value.\n \t Check the model to change', heter_groups), 'to fixed term ')), "blue")))

}
gc()
# if (is.null(heter_groups) & is.null(VarCov_str)){Inter_Gen_pos= NULL}
# if(length(Gen_pos) == length(rand_term)){Inter_Gen_pos= NULL}
if(!is.null(Inter_Gen_pos)){

  pred_heter_groups <- asreml::predict.asreml(mod, classify= rand_term[[Inter_Gen_pos]], sed=FALSE)$pvals

}
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
  #if(exists("G") & exists("omic1")){

    #ebv1 = ebv_G
    #ebv2 = ebv_omic1

    #### Check for the Generalized and unique Inverse
    for (n in 1:length(G_list)) {

    if("G"==names(BV_All)[n]){
    coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
    coeffRaw_G <- colMeans(coeffRaw_G)
    }
    ##
      if("omic1"==names(BV_All)[n]){
    coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
    coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }
}

  } else if(length(intersect(c("G", "omic2"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_G
    #ebv2 = ebv_omic2
    ###
    #### Check for the Generalized and unique Inverse
    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
        coeffRaw_G <- colMeans(coeffRaw_G)
      }
      ##
      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
      }
    }


  } else if(length(intersect(c("G", "omic3"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_G
    #ebv2 = ebv_omic3
    ##
    #### Check for the Generalized and unique Inverse
    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
        coeffRaw_G <- colMeans(coeffRaw_G)
      }
      ##
      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
      }
    }

  } else if(length(intersect(c("omic1", "omic2"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_omic1
    #ebv2 = ebv_omic2
    ##
    #### Check for the Generalized and unique Inverse
    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }
      ##
      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
      }
    }

  } else if(length(intersect(c("omic1", "omic3"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_omic1
    #ebv2 = ebv_omic3
    ##
    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }
      ##
      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
      }
    }


  } else if(length(intersect(c("omic2", "omic3"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_omic2
    #ebv2 = ebv_omic3

    for (n in 1:length(G_list)) {

      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
      }
      ##
      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
      }
    }

  } else if(length(intersect(c("G","omic1", "omic2"),unlist(G_list)))==length(G_list)){

    #ebv1 = ebv_G
    #ebv2 = ebv_omic1
    #ebv3 = ebv_omic2



    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
        coeffRaw_G <- colMeans(coeffRaw_G)
      }

      if("omic1"==names(BV_All)[n]){
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }
      ##
      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
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
        coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
        coeffRaw_G <- colMeans(coeffRaw_G)
      }

      if("omic1"==names(BV_All)[n]){
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }
      ##
      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
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
        coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
        coeffRaw_G <- colMeans(coeffRaw_G)
      }

      ##
      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
      }
      ##

      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
      }
    }


  } else if(length(intersect(c("omic1","omic2", "omic3"),unlist(G_list)))==length(G_list)){

    # ebv1 = ebv_omic1
    # ebv2 = ebv_omic2
    # ebv3 = ebv_omic3

    ###
    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
      }

      ##
      if("omic2"==names(BV_All)[n]){
        coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
      }
      ##

      if("omic3"==names(BV_All)[n]){
        coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
        coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
      }
    }

  } else {

    if(length(intersect(c("G","omic1", "omic2", "omic3"),unlist(G_list)))==length(G_list)){

      # ebv1 = ebv_G
      # ebv2 = ebv_omic1
      # ebv3 = ebv_omic2
      # ebv4 = ebv_omic3

      ###
      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[n]]$BLUP
          coeffRaw_G <- colMeans(coeffRaw_G)
        }


        if("omic1"==names(BV_All)[n]){
          coeffRaw_omic1 = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[n]]$BLUP
          coeffRaw_omic1 <- colMeans(coeffRaw_omic1)
        }

        ##
        if("omic2"==names(BV_All)[n]){
          coeffRaw_omic2 = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[n]]$BLUP
          coeffRaw_omic2 <- colMeans(coeffRaw_omic2)
        }
        ##

        if("omic3"==names(BV_All)[n]){
          coeffRaw_omic3 = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[n]]$BLUP
          coeffRaw_omic3 <- colMeans(coeffRaw_omic3)
        }
      }

    }

  }

} else{

  if(!is.null(gmatrix)){


    coeffRaw_G = solve(t(gmatrix)*gmatrix)*t(gmatrix)*BV_All[[1]]$BLUP

  } else if(!is.null(gkernel)){


    coeffRaw_G = solve(t(gkernel)*gkernel)*t(gkernel)*BV_All[[1]]$BLUP

  } else if(!is.null(omic1_kernel)){

    coeffRaw_G = solve(t(omic1_kernel)*omic1_kernel)*t(omic1_kernel)*BV_All[[1]]$BLUP

  } else if(!is.null(omic2_kernel)){

    coeffRaw_G = solve(t(omic2_kernel)*omic2_kernel)*t(omic2_kernel)*BV_All[[1]]$BLUP

  } else if(!is.null(omic3_kernel)){

    coeffRaw_G = solve(t(omic3_kernel)*omic3_kernel)*t(omic3_kernel)*BV_All[[1]]$BLUP
  }

}



##### End

if (!is.null(VarCov_str) & length(Inter_Gen_pos)!=0){
  ####
  if(!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV=BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_G,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR= Res$Correlation,
                  COV=Res$Covariance)
  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV = BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic1,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR= Res$Correlation,
                  COV=Res$Covariance)
  } else if (is.null(gmatrix)  & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call = str.mod,
                  mod = mod,
                  EBV = BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR= Res$Correlation,
                  COV=Res$Covariance)

  } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV=BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR= Res$Correlation,
                  COV=Res$Covariance)

  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

    rm(ebv_G, ebv_omic1)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic2,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

    rm(ebv_G, ebv_omic2)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

    rm(ebv_G, ebv_omic3)
  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  EBV_3 = ebv_omic2,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  coefficients_3 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2,
                  CORR3= Res$Correlation3,
                  COV3=Res$Covariance3)

    rm(ebv_G, ebv_omic1, ebv_omic2)

  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2,
                  CORR3= Res$Correlation3,
                  COV3=Res$Covariance3)

    rm(ebv_G, ebv_omic1, ebv_omic3)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic2,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic2,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2,
                  CORR3= Res$Correlation3,
                  COV3=Res$Covariance3)

    rm(ebv_G, ebv_omic2, ebv_omic3)

  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }
    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic1,
                  EBV_2 = ebv_omic2,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

rm(ebv_omic1, ebv_omic2)
  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic1,
                  EBV_2 = ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

rm(ebv_omic1, ebv_omic3)

  } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic2,
                  EBV_2 = ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic2,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2)

    rm(ebv_omic2, ebv_omic3)

  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_omic1,
                  EBV_2 =  ebv_omic2,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic2,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =Res$Heritability,
                  varG_per_omics = Res$varG_per_omics,
                  Va = Res$Total_genetic_var,
                  Ve = Res$Residual_Var,
                  CORR1= Res$Correlation1,
                  COV1=Res$Covariance1,
                  CORR2= Res$Correlation2,
                  COV2=Res$Covariance2,
                  CORR3= Res$Correlation3,
                  COV3=Res$Covariance3)

    rm(ebv_omic1, ebv_omic2, ebv_omic3)
  } else {

    if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic1,
                    EBV_3 = ebv_omic2,
                    EBV_4 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic1,
                    coefficients_3 = coeffRaw_omic2,
                    coefficients_3 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =Res$Heritability,
                    varG_per_omics = Res$varG_per_omics,
                    Va = Res$Total_genetic_var,
                    Ve = Res$Residual_Var,
                    CORR1= Res$Correlation1,
                    COV1=Res$Covariance1,
                    CORR2= Res$Correlation2,
                    COV2=Res$Covariance2,
                    CORR3= Res$Correlation3,
                    COV3=Res$Covariance3,
                    CORR4= Res$Correlation4,
                    COV4=Res$Covariance4)

      rm(ebv_G, ebv_omic1, ebv_omic2, ebv_omic3)

    }
  }
} else if (is.null(VarCov_str) & length(Inter_Gen_pos)==0){

  if(!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV= BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_G,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)
  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV = BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic1,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)
  } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    Result = list(call = str.mod,
                  mod = mod,
                  EBV = BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

  } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    Result = list(call=str.mod,
                  mod=mod,
                  EBV=BV_All[[1]]$BLUP,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic1)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic2,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic2)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic3)

  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  EBV_3 = ebv_omic2,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  coefficients_3 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic1, ebv_omic2)

  } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic1,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic1,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic1, ebv_omic3)

  } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("G"==names(BV_All)[n]){
        ebv_G = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_G,
                  EBV_2 =  ebv_omic2,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_G,
                  coefficients_2 = coeffRaw_omic2,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_G, ebv_omic2, ebv_omic3)

  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic1,
                  EBV_2 = ebv_omic2,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic2,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_omic1, ebv_omic2)

  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic1,
                  EBV_2 = ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

rm(ebv_omic1, ebv_omic3)

  } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 =  ebv_omic2,
                  EBV_2 = ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic2,
                  coefficients_2 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_omic2, ebv_omic3)

  } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

    for (n in 1:length(G_list)) {

      if("omic1"==names(BV_All)[n]){
        ebv_omic1 = BV_All[[n]]$BLUP
      }

      if("omic2"==names(BV_All)[n]){
        ebv_omic2 = BV_All[[n]]$BLUP
      }

      if("omic3"==names(BV_All)[n]){
        ebv_omic3 = BV_All[[n]]$BLUP
      }

    }

    Result = list(call=str.mod,
                  mod=mod,
                  EBV_1 = ebv_omic1,
                  EBV_2 =  ebv_omic2,
                  EBV_3 = ebv_omic3,
                  pred_value =pred_value,
                  #pred_heter_groups = pred_heter_groups,
                  coefficients_1 = coeffRaw_omic1,
                  coefficients_2 = coeffRaw_omic2,
                  coefficients_3 = coeffRaw_omic3,
                  BIC=summary(mod)$bic,
                  AIC=summary(mod)$aic,
                  H2 =H,
                  Va = varG_matrix,
                  Ve = VE)

    rm(ebv_omic1, ebv_omic2, ebv_omic3)

  } else {

    if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic1,
                    EBV_3 = ebv_omic2,
                    EBV_4 = ebv_omic3,
                    pred_value =pred_value,
                    #pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic1,
                    coefficients_3 = coeffRaw_omic2,
                    coefficients_3 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_G, ebv_omic1, ebv_omic2, ebv_omic3)

    }
  }
    # Result = list(call=str.mod,
    #               mod=mod,
    #               ebv=BV_All,
    #               pred_value =pred_value,
    #               BIC=summary(mod)$bic,
    #               AIC=summary(mod)$aic,
    #               H2 =H,
    #               Va = varG_matrix,
    #               Ve = VE)


  } else{

  #if (is.null(VarCov_str) & length(rand_inter_Pos)!=0){
  if ((is.null(VarCov_str) & length(Inter_Gen_pos)!=0) | (is.null(VarCov_str) & length(Inter_Gen_pos)==0)){

    if(!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      Result = list(call=str.mod,
                    mod=mod,
                    EBV=BV_All[[1]]$BLUP,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients = coeffRaw_G,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      Result = list(call=str.mod,
                    mod=mod,
                    EBV = BV_All[[1]]$BLUP,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients = coeffRaw_omic1,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
    } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      Result = list(call = str.mod,
                    mod = mod,
                    EBV = BV_All[[1]]$BLUP,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients = coeffRaw_omic2,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

    } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      Result = list(call=str.mod,
                    mod=mod,
                    EBV=BV_All[[1]]$BLUP,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }


      }


      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic1,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic1,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
rm(ebv_G, ebv_omic1)
    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }


      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic2,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic2,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
      rm(ebv_G, ebv_omic2)

    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }


      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_G, ebv_omic3)

    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic1,
                    EBV_3 = ebv_omic2,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic1,
                    coefficients_3 = coeffRaw_omic2,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
rm(ebv_G, ebv_omic1, ebv_omic2)
    } else if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic1,
                    EBV_3 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic1,
                    coefficients_3 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_G, ebv_omic1, ebv_omic3)

    } else if (!is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("G"==names(BV_All)[n]){
          ebv_G = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }
      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_G,
                    EBV_2 =  ebv_omic2,
                    EBV_3 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_G,
                    coefficients_2 = coeffRaw_omic2,
                    coefficients_3 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)
rm(ebv_G, ebv_omic2, ebv_omic3)

    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 =  ebv_omic1,
                    EBV_2 = ebv_omic2,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_omic1,
                    coefficients_2 = coeffRaw_omic2,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

rm(ebv_omic1, ebv_omic2)
    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 =  ebv_omic1,
                    EBV_2 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_omic1,
                    coefficients_2 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_omic1, ebv_omic3)

    } else if (is.null(gmatrix) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 =  ebv_omic2,
                    EBV_2 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_omic2,
                    coefficients_2 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_omic2, ebv_omic3)

    } else if (is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

      for (n in 1:length(G_list)) {

        if("omic1"==names(BV_All)[n]){
          ebv_omic1 = BV_All[[n]]$BLUP
        }


        if("omic2"==names(BV_All)[n]){
          ebv_omic2 = BV_All[[n]]$BLUP
        }

        if("omic3"==names(BV_All)[n]){
          ebv_omic3 = BV_All[[n]]$BLUP
        }

      }

      Result = list(call=str.mod,
                    mod=mod,
                    EBV_1 = ebv_omic1,
                    EBV_2 =  ebv_omic2,
                    EBV_3 = ebv_omic3,
                    pred_value =pred_value,
                    pred_heter_groups = pred_heter_groups,
                    coefficients_1 = coeffRaw_omic1,
                    coefficients_2 = coeffRaw_omic2,
                    coefficients_3 = coeffRaw_omic3,
                    BIC=summary(mod)$bic,
                    AIC=summary(mod)$aic,
                    H2 =H,
                    Va = varG_matrix,
                    Ve = VE)

      rm(ebv_omic1, ebv_omic2, ebv_omic3)

    } else {

      if (!is.null(gmatrix) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){

        for (n in 1:length(G_list)) {

          if("G"==names(BV_All)[n]){
            ebv_G = BV_All[[n]]$BLUP
          }

          if("omic1"==names(BV_All)[n]){
            ebv_omic1 = BV_All[[n]]$BLUP
          }


          if("omic2"==names(BV_All)[n]){
            ebv_omic2 = BV_All[[n]]$BLUP
          }

          if("omic3"==names(BV_All)[n]){
            ebv_omic3 = BV_All[[n]]$BLUP
          }

        }


        Result = list(call=str.mod,
                      mod=mod,
                      EBV_1 = ebv_G,
                      EBV_2 =  ebv_omic1,
                      EBV_3 = ebv_omic2,
                      EBV_4 = ebv_omic3,
                      pred_value =pred_value,
                      pred_heter_groups = pred_heter_groups,
                      coefficients_1 = coeffRaw_G,
                      coefficients_2 = coeffRaw_omic1,
                      coefficients_3 = coeffRaw_omic2,
                      coefficients_3 = coeffRaw_omic3,
                      BIC=summary(mod)$bic,
                      AIC=summary(mod)$aic,
                      H2 =H,
                      Va = varG_matrix,
                      Ve = VE)

        rm(ebv_G, ebv_omic1, ebv_omic2, ebv_omic3)

      }
    }
    # Result = list(call=str.mod,
    #               mod=mod,
    #               ebv=BV_All,
    #               pred_value =pred_value,
    #               pred_heter_groups = pred_heter_groups,
    #               BIC=summary(mod)$bic,
    #               AIC=summary(mod)$aic,
    #               H2 =H,
    #               Va = varG_matrix,
    #               Ve = VE)



  }

} ## Result Ends

return(Result)
}
