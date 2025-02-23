
#' Remove Variables from Global Environment
#'
#' This function removes specified variables from the global environment.
#'
#' @param var_names A character vector of variable names to remove from the global environment.
#'
#' @return None. The specified variables are removed from the global environment if they exist.
#' @examples
#' \dontrun{
#'   var1 <- 1
#'   var2 <- 2
#'   remove_from_global(c("var1", "var2"))
#' }
#' @export

remove_from_global <- function(var_names) {
  for (var_name in var_names) {
    if(exists(var_name, envir = .GlobalEnv)) {
      rm(list = var_name, envir = .GlobalEnv)
      #print(paste("Object", var_name, "removed from global environment."))
    } else {
      #print(paste("Object", var_name, "not found in global environment."))
    }
  }
}

#' Calculate Coefficients for ASReml Models
#'
#' This function calculates coefficients for ASReml models given the genetic matrix and estimated breeding values.
#'
#' @param gmatrix A matrix of genetic data.
#' @param ebv A numeric vector of estimated breeding values.
#' @param heter_groups A factor or vector specifying the heterogeneity groups.
#' @param heter_grp A factor or vector specifying the unique heterogeneity groups.
#' @param gid_name A character vector of genotype IDs.
#'
#' @return A data frame containing the coefficients for the ASReml models.
#' @examples
#' \dontrun{
#'   gmatrix <- matrix(rnorm(100), nrow=10)
#'   ebv <- rnorm(10)
#'   heter_groups <- factor(rep(1:2, each=5))
#'   heter_grp <- unique(heter_groups)
#'   gid_name <- as.character(1:10)
#'   coefficients <- cal_coeff_asreml(gmatrix, ebv, heter_groups, heter_grp, gid_name)
#'   print(coefficients)
#' }
#' @export
cal_coeff_asreml <- function(gmatrix,
                             ebv,
                             heter_groups,
                             heter_grp,
                             gid_name){

  coeff = solve(t(gmatrix)*gmatrix)*t(gmatrix)*ebv
  coeff <- colMeans(coeff)


  if(!is.null(heter_groups)){
    coeff <- data.frame(x_variables = rep(gid_name , length(unique(heter_grp))),
                        Env = rep(unique(heter_grp), each= length(gid_name)),
                        coeff = coeff,
                        stringsAsFactors = FALSE)
    names(coeff)[2] <- heter_groups



  }else {
    coeff <-  data.frame(x_variables = rownames(gmatrix),
                         coeff = coeff,
                         stringsAsFactors = FALSE)

  }

  return(coeff)
}

#' Extract and Process ASReml Model Output
#'
#' This function extracts and processes the output from ASReml models, including coefficients, predicted values, and variance components.
#'
#' @param mod_asreml A list containing the fitted ASReml model and related information.
#' @param pheno_data A data frame containing the phenotypic data.
#' @param response A string specifying the response variable in the phenotypic data.
#' @param gmatrix A matrix of genetic data.
#' @param omic1_kernel An optional kernel matrix for the first omic data.
#' @param omic2_kernel An optional kernel matrix for the second omic data.
#' @param omic3_kernel An optional kernel matrix for the third omic data.
#' @param omics_kernel_label A list of labels for the omic kernel matrices.
#' @param heter_groups A factor or vector specifying the heterogeneity groups.
#' @param gen_name A string specifying the column name for the genotypic data.
#' @param var_cov_str An optional string specifying the variance-covariance structure.
#' @param heter_resid An optional string specifying the heterogeneity of residuals.
#' @param pworkspace A numeric value specifying the workspace for ASReml. Default is 1e15.
#' @param maxit An integer specifying the maximum number of iterations for ASReml. Default is 50.
#' @param ... Additional arguments to be passed to underlying functions.
#'
#' @return A list containing the processed ASReml model output, including coefficients, predicted values, and variance components.
#' @examples
#' \dontrun{
#'   # Mock preparation of ASReml model components
#'   mod_asreml <- list(
#'     model = "asreml_model",
#'     str_mod = "model_string",
#'     gen_pos = "genetic_position",
#'     inter_gen_pos = "interaction_genetic_position",
#'     names_in_inv_list = c("GID"),
#'     rand_term = "random_term"
#'   )
#'   pheno_data <- data.frame(
#'     GID = 1:10,
#'     Trait = rnorm(10),
#'     Group = factor(rep(1:2, each = 5))
#'   )
#'   response <- "Trait"
#'   gmatrix <- matrix(rnorm(100), nrow = 10)
#'   result <- asreml_mod_output_new(
#'     mod_asreml = mod_asreml,
#'     pheno_data = pheno_data,
#'     response = response,
#'     gmatrix = gmatrix
#'   )
#'   print(result)
#' }
#' @export
#'

asreml_mod_output_new <- function(
    mod_asreml = NULL,
    pheno_data = NULL,
    response = NULL,
    gmatrix = NULL,
    omic1_kernel=NULL,
    omic2_kernel=NULL,
    omic3_kernel=NULL,
    omics_kernel_label = list(omic1_kernel = NULL,
                             omic2_kernel = NULL,
                             omic3_kernel = NULL),
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

  ##### Print lable
  if(inherits(omics_kernel_label,'list')){
    if(!all(sapply(omics_kernel_label, function(x){ is.null(x)}))!=FALSE){
      label <-  which(sapply(omics_kernel_label, function(x) !is.null(x)))
      print_lable <-  omics_kernel_label[label]
    } else {
      print_lable <-  NULL
    }

  }else {
    if(inherits(omics_kernel_label,"character")){
      print_lable <-  omics_kernel_label
    }
    if(is.null(omics_kernel_label)){
      print_lable <-  NULL
    }
  }
  ######

  mod <-  mod_asreml[["model"]]
  #mod = mod[[1]]
  ### if every variance component is stable the update will not run
  ## by default in asreml so it safe to keep it
  #mod <-  asreml::update.asreml(mod)
  str_mod <-  mod_asreml[["str_mod"]]
  gen_pos <-  mod_asreml[["gen_pos"]]
  inter_gen_pos <-  mod_asreml[["inter_gen_pos"]]
  names_in_inv_list <-  mod_asreml[["names_in_inv_list"]]
  rand_term <-  mod_asreml[["rand_term"]]

  tst <- NULL
  yy <- as.double(unique(data.frame(mod$mf)[, response]))


  sik_yna <- as.data.frame(mod$mf)

  tst_NA <- which(is.na(sik_yna[, response]))

  yNA <- sik_yna[, response]

  # Check if response column contains NAs and if all unique genotypes match the number of rows
  if (any(is.na(sik_yna[, response])) && length(unique(as.character(sik_yna[, gen_name]))) == nrow(sik_yna)) {

    tstt <- which(is.na(sik_yna[, response]))

    if (length(tstt) != 0) {
      train_gid <- unique(as.character(sik_yna[-tstt, gen_name]))
      unique_gid <- unique(as.character(sik_yna[, gen_name]))
      train_test_label_unique <- ifelse(unique_gid%in%train_gid, "Train", "Test")
      #train_test_label <- ifelse(is.na(sik_yna[, response]), "Test", "Train")
    }

  } else {
    if (nrow(sik_yna) > length(unique(sik_yna[, gen_name]))) {

      tstt <- which(is.na(sik_yna[, response]))

      if (length(tstt) != 0) {
        train_gid_across_env <- unique(as.character(sik_yna[-tstt, gen_name]))
        unique_gid_across_env <- unique(as.character(sik_yna[, gen_name]))
        train_test_label_unique <- ifelse(unique_gid_across_env%in%train_gid_across_env, "Train", "Test")
        train_test_label <- ifelse(is.na(sik_yna[, gen_name]), "Test", "Train")
      }
    }
  }


  #############################
  ### !is.null(var_cov_str) & is.null(inter_gen_pos) incase user provide var_cov_str
  ## while the data is not MT in nature
  if(!is.null(var_cov_str) & is.null(inter_gen_pos)){
    var_cov_str <-  NULL
    heter_groups <-  NULL
    heter_resid <-  NULL
  }

  #########################
  ### Set up parameter for predict function
  # asreml::asreml.options(trace=FALSE, workspace=workspace,
  #                        pworkspace = pworkspace, maxit = maxit)

  BLUP <- summary(mod, coef=TRUE)$coef.random

  colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

  #heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

  if (!is.null(heter_groups)){
    ### It possible the user provide the heter_groups while it actually a single environment,
    ## This will check and turn it off
    if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
      heter_groups <-  NULL
    } else{
      if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
        heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
        all_envs_for_met <- as.character(pheno_data[[heter_groups]])
      }
    }

  }
  #ENV_Ids = as.character(unique(pheno_data[[heter_groups]]))
  ##############################################
  #### Extract Breeding values/genetic effect estimate for all omics
  ##############################################
  estimated_breeding_value_list <- list()

  ##
  for (b in seq_along(names_in_inv_list)) {
    estimated_breeding_value_list[[names_in_inv_list[b]]] <- BLUP[grep(paste(names_in_inv_list[b],"\\)", sep = ""),rownames(BLUP)),]

  }
  ###
  ### For variance structure extraction
  if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

    if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
      ## Extract the number of factors
      #N_fa = substr(var_cov_str, 3, 100)

      for (bb in seq_along(names_in_inv_list)) {
        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][!rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])%in%rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]][grep('Comp',rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]])),]), ]

        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

        estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,3])

      }

    } else {

      if (var_cov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {

        for (bb in seq_along(names_in_inv_list)) {

          estimated_breeding_value_list[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list[[names_in_inv_list[bb]]])

          estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]), "\\)_", 3)[,2])

        }


      }


    } ## End

    #################################
    if(!is.null(inter_gen_pos)){
      for (bb in seq_along(names_in_inv_list)) {

        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(4, 1:2)]

        estimated_breeding_value_list[[names_in_inv_list[bb]]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list[[names_in_inv_list[bb]]][, gen_name])))

        estimated_breeding_value_list[[names_in_inv_list[bb]]] <- estimated_breeding_value_list[[names_in_inv_list[bb]]][, c(1, 4, 2:3)]

        colnames(estimated_breeding_value_list[[names_in_inv_list[bb]]])[1:3] <- c(gen_name, heter_groups, "BLUP")

        estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Standard_error"]^2
        rownames(estimated_breeding_value_list[[names_in_inv_list[bb]]]) = NULL


      }

    }
    #######################################

    res_herit_varCov <-  tryCatch({

    Res_Va_Ve_H2_COV_COR <-  asreml_herit_varCov_new(model= mod,
                                                 heter_groups= heter_groups,
                                                 var_cov_str= var_cov_str,
                                                 heter_resid= heter_resid,
                                                 names_in_inv_list = names_in_inv_list,
                                                 inter_gen_pos = inter_gen_pos,
                                                 gen_pos = gen_pos)

    Res_Va_Ve_H2_COV_COR

      },
    error = function(e) {
      # Handle the error, you can print a message or take other actions
      cat(paste("Error in processing the variance components:", conditionMessage(e), "\n"))
      return(NULL)  # Return NULL or an appropriate value to indicate the failure
    }
    )

    #VA = Res$Genetic_Var
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in seq_along(names_in_inv_list)){

      if(!is.null(res_herit_varCov)){

        Res_Va_Ve_H2_COV_COR <-  res_herit_varCov
      if(length(names_in_inv_list)>1){
        VA <-  Res_Va_Ve_H2_COV_COR[["varG_per_omics"]][[bb]]
      } else {
        if(length(names_in_inv_list)==1){
          VA <-  unlist(Res_Va_Ve_H2_COV_COR[["Total_genetic_var"]])
        }

      }
      } else{
        VA <- NULL
      }

    if(is.null(VA)){
      if(length(VA)< length(heter_grp)){
        estimated_breeding_value_list[[bb]][, "Reliability"] = NA
        message(paste( insight::print_color("WARNINGS\n", "blue"),
                       insight::print_color(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.\n")), "blue")))
        #message(paste(msg,paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.")))
      } else{
        if(length(VA) == length(heter_grp)){
          estimated_breeding_value_list[[bb]][, "Reliability"] <- NA
          #BV$Reliability = NA
          for (i in 1:length(VA)) {
            estimated_breeding_value_list[[bb]][, "Reliability"] <- ifelse(estimated_breeding_value_list[[bb]][, heter_groups]%in% heter_grp[i],
                                                    round(1 - (estimated_breeding_value_list[[bb]][, "Prediction_error_variance"]/VA[i]),6), estimated_breeding_value_list[[bb]][, "Reliability"])

          }
        }


      }
##########
} else {

    estimated_breeding_value_list[[bb]][, "Reliability"] = NA
    message(paste( insight::print_color("WARNINGS\n", "blue"),
                   insight::print_color(paste(msg,paste("There is a problem with the variance component, reliability cannot be estimated.\n")), "blue")))

}

    }


  } else {
    ## Problem
    if(is.null(var_cov_str) & is.null(inter_gen_pos) ){
      for (bb in seq_along(names_in_inv_list)){
        estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
        estimated_breeding_value_list[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[bb]]), "\\)_", 3)[,2])
        rownames(estimated_breeding_value_list[[bb]]) <-  NULL
        estimated_breeding_value_list[[bb]] <-   estimated_breeding_value_list[[bb]][, c(4, 1:2)]
        colnames(estimated_breeding_value_list[[bb]])[1:2] <- c(gen_name, "BLUP")
      }

      vc <- summary(mod)$varcomp

      #VAR_check_Pos <- which(vc$bound=='F' | vc$bound=='U' | vc$bound=="?" | vc$bound=="S")

      res_comp_check <-  tryCatch({
        VAR_check_Pos <- which(vc$bound=="?" | vc$bound=="S")
        if(length(VAR_check_Pos)>1) {

          # stop(message(paste( insight::print_color("STOP\n", "red"),
          #                insight::print_color(paste(msg, paste(paste("variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
          #                                                      " is unstable, refix the model")), "red"))), call. = FALSE)
          #
          stop(print(paste(msg, paste(paste("Variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" and" ),
                                      " is unstable, refix the model.\n"))), call. = FALSE)

        } else {

          if(length(VAR_check_Pos)==1) {

            stop(print(paste(msg, paste(paste("Variance component for", as.character(rownames(vc)[VAR_check_Pos]), collapse =" " ),
                                        " is unstable, refix the model.\n"))), call. = FALSE)

          }

        }
           vc
        },
        error = function(e) {
          # Handle the error, you can print a message or take other actions
          cat(paste("Variance component  for", conditionMessage(e), "\n"))
          return(NULL)  # Return NULL or an appropriate value to indicate the failure
        }
      )

      #### Extration of the genetic variance for all omics to calcuate heritability.
      VarG_All <- vector("list", length = length(names_in_inv_list))
      names(VarG_All) <- names_in_inv_list
      VarG <- vc[grep(paste("^vm\\(", gen_name, sep = ""), rownames(vc)), drop=FALSE, ]

      for (g in seq_along(names_in_inv_list)) {

        VarG_All[[g]] <-VarG[grep(names_in_inv_list[g], rownames(VarG)), "component"]

      }
#### Residual variance
      VE <- vc[grep("!R", row.names(vc)), "component"]
      varG_matrix <-  unlist(VarG_All)
      H <-  matrix(NA, nrow = 1, ncol = length(VE))

      #####
      for (bb in seq_along(names_in_inv_list)){
        estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <- estimated_breeding_value_list[[bb]][, "Standard_error"]^2
        estimated_breeding_value_list[[bb]][, "Reliability"] <- round(1 - (estimated_breeding_value_list[[bb]][, "Prediction_error_variance"]/as.double((VarG_All[bb]))),6)

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

      #colnames(VarG_All) <- "Variance"
      Total_genetic_var <- sum(VarG_All[, 1])
#### Single location
      Res_Va_Ve_H2_COV_COR <-  list(Heritability = H,
                                    varG_per_omics = VarG_All,
                                    Total_genetic_var = Total_genetic_var,
                                    Residual_Var = VE)
    } ## End


  } #### End var_Covar

  ## when variance_covariance structure is not defined by the user and CS is used
  ### For variance structure extraction
  if(is.null(var_cov_str) & !is.null(inter_gen_pos)){
    for (bb in seq_along(names_in_inv_list)) {
      estimated_breeding_value_list[[bb]] <- as.data.frame(estimated_breeding_value_list[[bb]])
      estimated_breeding_value_list[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list[[bb]]), "\\)_", 3)[,2])
      rownames(estimated_breeding_value_list[[bb]]) <-  NULL
    }
    #################################
    if(!is.null(inter_gen_pos)){
      for (bb in seq_along(names_in_inv_list)) {
        estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(4, 1:2)]
        estimated_breeding_value_list[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list[[bb]][, gen_name])))
        estimated_breeding_value_list[[bb]] <- estimated_breeding_value_list[[bb]][, c(1, 4, 2:3)]
        colnames(estimated_breeding_value_list[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
        estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[bb]][, "Standard_error"]^2
        rownames(estimated_breeding_value_list[[bb]]) = NULL
      }
    }

    #######################################
    Res_Va_Ve_H2_COV_COR <-  asreml_herit_CSM_new(model= mod,
                                                 heter_groups= heter_groups,
                                                 heter_resid= heter_resid,
                                                 names_in_inv_list = names_in_inv_list,
                                                 inter_gen_pos = inter_gen_pos,
                                                 gen_pos = gen_pos)

    varG_matrix <- Res_Va_Ve_H2_COV_COR[["varG_per_omics"]]
    ##
    heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    for (bb in seq_along(names_in_inv_list)){
      if(length(names_in_inv_list)>1){
        VA <-  Res_Va_Ve_H2_COV_COR[["varG_per_omics"]][bb, ]
      } else {
        if(length(names_in_inv_list)==1){
          VA <-  Res_Va_Ve_H2_COV_COR[["Total_genetic_var"]][1, ]
        }
      }

      if(length(VA)< length(heter_grp)){
        print(paste("Reliability cannot be estimated. Not all varaince components for", heter_groups, "are postive definitive.\n"))
      } else{
        if(length(VA) == length(heter_grp)){
          estimated_breeding_value_list[[bb]][, "Reliability"] <- NA
          ###
          for (i in 1:length(VA)) {
            estimated_breeding_value_list[[bb]][, "Reliability"] <- ifelse(estimated_breeding_value_list[[bb]][, heter_groups]%in% heter_grp[i],
                                                    round(1 - (estimated_breeding_value_list[[bb]][, "Prediction_error_variance"]/VA[i]),6), estimated_breeding_value_list[[bb]][, "Reliability"])
          }
        } else {
          stop(print(paste(msg, 'Genetic variance is missing')), call. = FALSE)
        }
      }
    }
  } ## End CS
####################
  ##### For MET analysis
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    ### Extract all environments in MET
    all_envs_for_met <-  as.character(pheno_data[,heter_groups])

    #### End MET
  } else{ ##  start single enviornment results
    ##################################

    Zg <- NULL
    all_envs_for_met <- NULL
  } ##  End single enviornment results
  #############################
  ### Predicted Values
  tst <- which(is.na(yy))

  gc() ## to free memory

  result_single_loc <-tryCatch({
  predicted_value <- asreml::predict.asreml(mod, classify=gen_name, sed=FALSE)$pvals
  predicted_value <-  predicted_value[, -ncol(predicted_value)] ### Remove status
  colnames(predicted_value)[colnames(predicted_value)%in%c("predicted.value", "std.error")] <- c("Predicted_value", "Standard_error")
  predicted_value[, "Prediction_error_variance"] <-  predicted_value[, "Standard_error"]^2
  #predicted_value[, "Reliability"] <-  NA
  predicted_value[, "Train_Test_Label"] <-  train_test_label_unique
  predicted_value

    },
  error = function(e) {
    # Handle the error, you can print a message or take other actions
    cat(paste("Error in prediction for single location GBLUP model. BLUPs value will be provided:", conditionMessage(e), "\n"))
    return(NULL)  # Return NULL or an appropriate value to indicate the failure
  }
  )

  if(is.null(result_single_loc)){
    BLUP <- summary(mod, coef=TRUE)$coef.random

    colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

    #heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

    if (!is.null(heter_groups)){
      ### It possible the user provide the heter_groups while it actually a single environment,
      ## This will check and turn it off
      if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
        heter_groups <-  NULL
        inter_gen_pos <-  NULL
      } else{
        if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
          heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
          all_envs_for_met <- as.character(pheno_data[[heter_groups]])
        }
      }

    }

    ##############################################
    #### Extract Breeding values/genetic effect estimate for all omics
    ##############################################
    estimated_breeding_value_listt <- list()

    ##
    for (b in seq_along(names_in_inv_list)) {
      estimated_breeding_value_listt[[names_in_inv_list[b]]] <- BLUP[grep(paste(names_in_inv_list[b],"\\)", sep = ""),rownames(BLUP)),]

    }
    ###
    ### For variance structure extraction
    if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

      if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
        ## Extract the number of factors
        #N_fa = substr(var_cov_str, 3, 100)

        for (bb in seq_along(names_in_inv_list)) {
          estimated_breeding_value_listt[[names_in_inv_list[bb]]] <- estimated_breeding_value_listt[[names_in_inv_list[bb]]][!rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]])%in%rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]][grep('Comp',rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]])),]), ]

          estimated_breeding_value_listt[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_listt[[names_in_inv_list[bb]]])

          estimated_breeding_value_listt[[names_in_inv_list[bb]]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]]), "\\)_", 3)[,3])

        }

      } else {

        if (var_cov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {

          for (bb in seq_along(names_in_inv_list)) {

            estimated_breeding_value_listt[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_listt[[names_in_inv_list[bb]]])

            estimated_breeding_value_listt[[names_in_inv_list[bb]]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]]), "\\)_", 3)[,2])

          }


        }


      } ## End

      #################################
      if(!is.null(inter_gen_pos)){
        for (bb in seq_along(names_in_inv_list)) {

          estimated_breeding_value_listt[[names_in_inv_list[bb]]] <- estimated_breeding_value_listt[[names_in_inv_list[bb]]][, c(4, 1:2)]

          estimated_breeding_value_listt[[names_in_inv_list[bb]]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_listt[[names_in_inv_list[bb]]][, gen_name])))

          estimated_breeding_value_listt[[names_in_inv_list[bb]]] <- estimated_breeding_value_listt[[names_in_inv_list[bb]]][, c(1, 4, 2:3)]

          colnames(estimated_breeding_value_listt[[names_in_inv_list[bb]]])[1:3] <- c(gen_name, heter_groups, "BLUP")

          estimated_breeding_value_listt[[names_in_inv_list[bb]]][, "Prediction_error_variance"] <-  estimated_breeding_value_listt[[names_in_inv_list[bb]]][, "Standard_error"]^2
          rownames(estimated_breeding_value_listt[[names_in_inv_list[bb]]]) <-  NULL


        }

      }

      # Bind all data frames into a single data frame
      combined_df <- dplyr::bind_rows(estimated_breeding_value_listt)

    }

    ### For compound symmetric
    if(is.null(var_cov_str) & !is.null(inter_gen_pos)){
      for (bb in seq_along(names_in_inv_list)) {
        estimated_breeding_value_listt[[bb]] <- as.data.frame(estimated_breeding_value_listt[[bb]])
        estimated_breeding_value_listt[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_listt[[bb]]), "\\)_", 3)[,2])
        estimated_breeding_value_listt[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_listt[[bb]][, "Standard_error"]^2
        rownames(estimated_breeding_value_listt[[bb]]) <-  NULL
      }
      #################################
      if(!is.null(inter_gen_pos)){
        for (bb in seq_along(names_in_inv_list)) {
          estimated_breeding_value_listt[[bb]] <- estimated_breeding_value_listt[[bb]][, c(4, 1:2)]
          estimated_breeding_value_listt[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_listt[[bb]][, gen_name])))
          estimated_breeding_value_listt[[bb]] <- estimated_breeding_value_listt[[bb]][, c(1, 4, 2:3)]
          colnames(estimated_breeding_value_listt[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
          estimated_breeding_value_listt[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_listt[[bb]][, "Standard_error"]^2
          rownames(estimated_breeding_value_listt[[bb]]) <- NULL
        }
      }

      # Bind all data frames into a single data frame
      combined_df <- dplyr::bind_rows(estimated_breeding_value_listt)

      }


    #}

    if(is.null(var_cov_str) & is.null(inter_gen_pos) ){
      for (bb in seq_along(names_in_inv_list)){
        estimated_breeding_value_listt[[bb]] <- as.data.frame(estimated_breeding_value_listt[[bb]])
        estimated_breeding_value_listt[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_listt[[bb]]), "\\)_", 3)[,2])
        rownames(estimated_breeding_value_listt[[bb]]) <-  NULL

        estimated_breeding_value_listt[[bb]] <-   estimated_breeding_value_listt[[bb]][, c(4, 1:2)]
        colnames(estimated_breeding_value_listt[[bb]])[1:2] <- c(gen_name, "BLUP")
      }

      combined_df <- dplyr::bind_rows(estimated_breeding_value_listt)

    }

    # Aggregate total BLUP and standard error across all environments
    total_summarized_blup <- combined_df |>
      dplyr::ungroup() |>
      #dplyr::group_by(!!rlang::sym(gen_name)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(gen_name))) |>
      dplyr::summarise(
        BLUP = mean(BLUP, na.rm = TRUE),
        Standard_error = mean(Standard_error, na.rm = TRUE),
        #Standard_error = sqrt(sum(Standard_error^2, na.rm = TRUE)),
        .groups = 'drop'
      )

    # Calculate the Prediction Error Variance
    summarized_blup <- total_summarized_blup |>
      dplyr::mutate(Prediction_error_variance = Standard_error^2)
    summarized_blup <- as.data.frame(summarized_blup)
    predicted_value <- summarized_blup
    predicted_value[, "Train_Test_Label"] <- train_test_label_unique
  }
  ## sik

  if(!is.null(Zg)){
    genotype_means <- pheno_data |>
      dplyr::group_by(!!rlang::sym(gen_name)) |>
      dplyr::summarise(mean_value = mean(!!rlang::sym(response), na.rm = TRUE))
    genotype_means <-  as.data.frame(genotype_means)
    # Reorder genotype_means based on name in predicted_value
    genotype_means <- genotype_means[match(predicted_value[, gen_name], genotype_means[[gen_name]]), ]
    yy <- genotype_means[, "mean_value"]
  }

  if(length(tst_NA)!=0){
    if(!is.null(result_single_loc)){
      residual_value <- data.frame(name =  predicted_value[-tst_NA, gen_name],
                                   #Env = all_envs_for_met[tst],
                                   Predicted_value = predicted_value[-tst_NA, "Predicted_value"],
                                   Residual_value = (yNA[-tst_NA] - predicted_value[-tst_NA, "Predicted_value"]),
                                   stringsAsFactors = FALSE)

      colnames(residual_value)[1] <- gen_name

    } else{

    residual_value <- NULL

    }

  } else {
    if(!is.null(result_single_loc)){
    residual_value <- data.frame(name = unique(as.character(predicted_value[, gen_name])),
                                 #Env = all_envs_for_met,
                                 Predicted_value = predicted_value[, "Predicted_value"],
                                 Residual_value = (yNA - predicted_value[, "Predicted_value"]),
                                 stringsAsFactors = FALSE)
    colnames(residual_value)[1] <- gen_name

    } else {
      residual_value <- NULL
    }

  }


  if(!is.null(heter_groups)){
    if(var(predicted_value[, "Predicted_value"])==0){
      ### Check if the the across
      message(paste( insight::print_color("WARNING\n", "blue"),
                     insight::print_color(paste(msg, paste(paste('The average prediction across', heter_groups), paste('is a constant value.\n \t Check the model to change', heter_groups), 'to fixed term.\n ')), "blue")))
    }
  }
  gc()
  ################

  if (is.null(heter_groups) & is.null(var_cov_str)) inter_gen_pos <-  NULL
  if(length(gen_pos) == length(rand_term)) inter_gen_pos <- NULL

  if(!is.null(inter_gen_pos)){

   result_met <-  tryCatch({
    #pred_heter_groups <- asreml::predict.asreml(mod, classify= rand_term[[inter_gen_pos]], vcov = TRUE, aliased = T)
    across_env_predicted_value <- asreml::predict.asreml(mod, classify= rand_term[[inter_gen_pos]], sed=FALSE)$pvals
    across_env_predicted_value =  across_env_predicted_value[, -ncol(across_env_predicted_value)] ### Remove status
    colnames(across_env_predicted_value)[colnames(across_env_predicted_value)%in%c("predicted.value", "std.error")] <- c("Predicted_value", "Standard_error")
    across_env_predicted_value[, "Prediction_error_variance"] <-  across_env_predicted_value[, "Standard_error"]^2
    #across_env_predicted_value[, "Reliability"] <- NA

    across_env_predicted_value[, "Train_Test_Label"] <- train_test_label

    #name_across_env <-  paste("Across", paste(heter_groups, "Predicted_value", sep = "_"), sep = "_")
    if(length(tst)!=0){
      residual_value_met <- data.frame(name =  across_env_predicted_value[, gen_name][-tst_NA],
                                       Env = all_envs_for_met[-tst_NA],
                                       Predicted_value = across_env_predicted_value[, "Predicted_value"][-tst_NA],
                                       Residual_value = (yNA[-tst_NA] - across_env_predicted_value[, "Predicted_value"][-tst_NA]),
                                       stringsAsFactors = FALSE)
      colnames(residual_value_met)[1:2] <- c(gen_name, heter_groups)
    } else {
      residual_value_met <- data.frame(name = across_env_predicted_value[, gen_name],
                                       Env = all_envs_for_met,
                                       Predicted_value = across_env_predicted_value[, "Predicted_value"],
                                       Residual_value = (yNA - across_env_predicted_value[, "Predicted_value"]),
                                       #Residual_value = (yy - across_env_predicted_value[, "Predicted_value"]),
                                       stringsAsFactors = FALSE)
      colnames(residual_value_met)[1:2] <- c(gen_name, heter_groups)
    }
    list(across_env_predicted_value = across_env_predicted_value,
         residual_value_met = residual_value_met)
      },
    error = function(e) {
      # Handle the error, you can print a message or take other actions
      #cat(paste("Error in across", heter_groups, "prediction:", conditionMessage(e), "\n"))
      cat(paste("Error in multi-environment GS prediction. BLUP will be provided:", conditionMessage(e), "\n"))
      return(NULL)  # Return NULL or an appropriate value to indicate the failure
    }
    )

    if(is.null(result_met)) {
      BLUP <- summary(mod, coef=TRUE)$coef.random

      colnames(BLUP)[colnames(BLUP)%in%"std.error"] <- "Standard_error"

      #heter_grp <- as.character(unique(data.frame(mod$mf)[, heter_groups]))

      if (!is.null(heter_groups)){
        ### It possible the user provide the heter_groups while it actually a single environment,
        ## This will check and turn it off
        if(length(pheno_data[,gen_name])==length(unique(pheno_data[,gen_name]))){
          heter_groups <-  NULL
          inter_gen_pos <-  NULL
        } else{
          if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
            heter_grp <- as.character(unique(pheno_data[[heter_groups]]))
            all_envs_for_met <- as.character(pheno_data[[heter_groups]])
          }
        }

      }

      ##############################################
      #### Extract Breeding values/genetic effect estimate for all omics
      ##############################################
      estimated_breeding_value_list_t <- list()

      ##
      for (b in seq_along(names_in_inv_list)) {
        estimated_breeding_value_list_t[[names_in_inv_list[b]]] <- BLUP[grep(paste(names_in_inv_list[b],"\\)", sep = ""),rownames(BLUP)),]

      }
      ###
      ### For variance structure extraction
      if(!is.null(var_cov_str) & !is.null(inter_gen_pos)){

        if(isTRUE(grepl("fa", var_cov_str)) | isTRUE(grepl("rr", var_cov_str))){
          ## Extract the number of factors
          #N_fa = substr(var_cov_str, 3, 100)

          for (bb in seq_along(names_in_inv_list)) {
            estimated_breeding_value_list_t[[names_in_inv_list[bb]]] <- estimated_breeding_value_list_t[[names_in_inv_list[bb]]][!rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]])%in%rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]][grep('Comp',rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]])),]), ]

            estimated_breeding_value_list_t[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list_t[[names_in_inv_list[bb]]])

            estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, gen_name]<-as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]]), "\\)_", 3)[,3])

          }

        } else {

          if (var_cov_str %in% c("us", "corgh", "corgv", "corh", "corv")) {

            for (bb in seq_along(names_in_inv_list)) {

              estimated_breeding_value_list_t[[names_in_inv_list[bb]]] <- as.data.frame(estimated_breeding_value_list_t[[names_in_inv_list[bb]]])

              estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]]), "\\)_", 3)[,2])

            }


          }


        } ## End

        #################################
        if(!is.null(inter_gen_pos)){
          for (bb in seq_along(names_in_inv_list)) {

            estimated_breeding_value_list_t[[names_in_inv_list[bb]]] <- estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, c(4, 1:2)]

            estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, gen_name])))

            estimated_breeding_value_list_t[[names_in_inv_list[bb]]] <- estimated_breeding_value_list_t[[names_in_inv_list[bb]]][, c(1, 4, 2:3)]

            colnames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]])[1:3] <- c(gen_name, heter_groups, "BLUP")

            #estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[names_in_inv_list[bb]]][, "Standard_error"]^2
            rownames(estimated_breeding_value_list_t[[names_in_inv_list[bb]]]) <-  NULL


          }

        }

        # Bind all data frames into a single data frame
        combined_df <- dplyr::bind_rows(estimated_breeding_value_list_t)

      }

      ### For compound symmetric
      if(is.null(var_cov_str) & !is.null(inter_gen_pos)){

        for (bb in seq_along(names_in_inv_list)) {
          estimated_breeding_value_list_t[[bb]] <- as.data.frame(estimated_breeding_value_list_t[[bb]])
          estimated_breeding_value_list_t[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list_t[[bb]]), "\\)_", 3)[,2])
          rownames(estimated_breeding_value_list_t[[bb]]) <-  NULL
        }
        #################################
        if(!is.null(inter_gen_pos)){
          for (bb in seq_along(names_in_inv_list)) {
            estimated_breeding_value_list_t[[bb]] <- estimated_breeding_value_list_t[[bb]][, c(4, 1:2)]
            estimated_breeding_value_list_t[[bb]][, heter_groups] <- rep(heter_grp, each=length(unique(estimated_breeding_value_list_t[[bb]][, gen_name])))
            estimated_breeding_value_list_t[[bb]] <- estimated_breeding_value_list_t[[bb]][, c(1, 4, 2:3)]
            colnames(estimated_breeding_value_list_t[[bb]])[1:3] <- c(gen_name, heter_groups, "BLUP")
            #estimated_breeding_value_list[[bb]][, "Prediction_error_variance"] <-  estimated_breeding_value_list[[bb]][, "Standard_error"]^2
            rownames(estimated_breeding_value_list_t[[bb]]) <- NULL
          }
        }

        # Bind all data frames into a single data frame
        combined_df <- dplyr::bind_rows(estimated_breeding_value_list_t)

      }

      if(is.null(var_cov_str) & is.null(inter_gen_pos)){
        for (bb in seq_along(names_in_inv_list)){
          estimated_breeding_value_list_t[[bb]] <- as.data.frame(estimated_breeding_value_list_t[[bb]])
          estimated_breeding_value_list_t[[bb]][, gen_name] <- as.character(stringr::str_split_fixed(rownames(estimated_breeding_value_list_t[[bb]]), "\\)_", 3)[,2])
          rownames(estimated_breeding_value_list_t[[bb]]) <-  NULL

          estimated_breeding_value_list_t[[bb]] <-   estimated_breeding_value_list_t[[bb]][, c(4, 1:2)]
          colnames(estimated_breeding_value_list_t[[bb]])[1:2] <- c(gen_name, "BLUP")
        }

        combined_df <- dplyr::bind_rows(estimated_breeding_value_list_t)


      }

      env_summarized_blup <- combined_df |>
        dplyr::group_by(!!rlang::sym(gen_name), !!rlang::sym(heter_groups))|>
        dplyr::summarise(
          BLUP = mean(BLUP, na.rm = TRUE),
          #BLUP = sum(BLUP, na.rm = TRUE),
          Standard_error = mean(Standard_error, na.rm = TRUE),
          #Standard_error = sqrt(sum(Standard_error^2, na.rm = TRUE)),
          .groups = 'drop'  # Ensure the resulting data frame is not grouped
        )

      final_env_summarized_blup <- env_summarized_blup |>
        dplyr::mutate(Prediction_error_variance = Standard_error^2)
      summarized_blup <- as.data.frame(final_env_summarized_blup)

      across_env_predicted_value <- final_env_summarized_blup
      residual_value_met <- NULL
    } else{
      across_env_predicted_value <- result_met[["across_env_predicted_value"]]
      across_env_predicted_value[, "Train_Test_Label"] <- train_test_label
      residual_value_met <- result_met[["residual_value_met"]]
    }

  }

  ### Initialize step to calculate coefficient for each omics
### Gather all the variables in the function some might be null
  datasets <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
  dataset_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
  datasets_index <- which(!sapply(datasets, is.null)) ## identify null variable

  datasets <-  datasets[datasets_index] ## remove null variable
  dataset_names <- dataset_names[datasets_index]
  coefficients_list <-   list()
  m_matrix_model_ready_list <-   list()
  sum_estimated_breeding_value <-  0

  ## names_in_inv_list used to fit asreml model contains names of the kernel/gmatrix with inv this will identify it and remove it
  ## The idea is to make sure the order of the name is the same here and in names_in_in_list
  extracted_names_from_inv_list <- gsub("_inv", "", names_in_inv_list, ignore.case = TRUE)
  index_names_inv_extracted_omic <- grep("omic", extracted_names_from_inv_list, ignore.case = TRUE)
  index_names_inv_extracted_geno <- grep("gmatrix", extracted_names_from_inv_list, ignore.case = TRUE)

  if (length(index_names_inv_extracted_omic) > 1) {
    extracted_names_from_inv_list[index_names_inv_extracted_omic] <- paste(extracted_names_from_inv_list[index_names_inv_extracted_omic], "kernel", sep = "_")
  } else {
    if (length(index_names_inv_extracted_omic) == 1) {
      ## When ist is length 1, we just want omic_kernel irrespective of wither omic1 or omic2 or omic is provided
      extracted_names_from_inv_list[index_names_inv_extracted_omic] <- paste("omic", "kernel", sep = "_")
      dataset_names[index_names_inv_extracted_omic] <- paste("omic", "kernel", sep = "_")
    }
  }

  if (length(index_names_inv_extracted_geno)==1) {
    extracted_names_from_inv_list[index_names_inv_extracted_geno] <- "gmatrix"
  }


if(length(dataset_names)==length(extracted_names_from_inv_list)) {
  # Reorder dataset_names based on the order of name in extracted_names_from_inv_list
  dataset_names <- dataset_names[match(extracted_names_from_inv_list, dataset_names)]

} else {
  stop(paste(msg, "Different names were provided."), call. = FALSE)
}

  ## Intialize step to calculate EBV across environment for each omics
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    if (!is.null(dataset)) {
      if(!is.null(Zg)){
        ZgZg <- Zg%*%dataset%*%t(Zg)
        suppressMessages({
          ZgZg <- grm_kernel_precheck(ZgZg)
        })

      }
      ## This function calcuate EBV for each omics
      coefficients_list[[paste("coefficient",dataset_names[i], sep = "_")]] <- cal_coeff_asreml(gmatrix = if(!is.null(Zg)) ZgZg else dataset,
                                                                                                ebv = estimated_breeding_value_list[[i]][, "BLUP"],
                                                                                                heter_groups = heter_groups,
                                                                                                heter_grp = all_envs_for_met,
                                                                                                gid_name = rownames(dataset))

      sum_estimated_breeding_value <- sum_estimated_breeding_value + estimated_breeding_value_list[[i]][["BLUP"]]

      #colnames(estimated_breeding_value_list[[i]])[which("BLUP"%in%colnames(estimated_breeding_value_list[[i]]))] <- "sik" # Estimated_breeding_value

      if(dataset_names[i]=="gmatrix"){
        m_matrix_model_ready_list[[paste(gsub("gmatrix", "geno", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
      } else{
        m_matrix_model_ready_list[[paste(gsub("_kernel", "", dataset_names[i]), "model_ready", sep = "_")]] <- dataset
      }
    }


  }
#### process the sum_ebv
  combined_dff <- dplyr::bind_rows(estimated_breeding_value_list)
  if(length(datasets)>1){
    summarized_blup_use <- combined_dff |>
      dplyr::group_by(!!rlang::sym(gen_name))|>
      dplyr::summarise(
        #BLUP = sum(BLUP, na.rm = TRUE),
        BLUP = mean(BLUP, na.rm = TRUE),
        Standard_error = mean(Standard_error, na.rm = TRUE),
        #Standard_error = sqrt(sum(Standard_error^2, na.rm = TRUE)),
        #Reliability = NA,
        .groups = 'drop'  # Ensure the resulting data frame is not grouped
      )

    # Calculate the Prediction Error Variance
    sum_ebv <- summarized_blup_use |>
      dplyr::mutate(Prediction_error_variance = Standard_error^2,
                    Reliability = NA)

  } else{
    if(length(datasets)==1){
      sum_ebv <- dplyr::bind_rows(estimated_breeding_value_list)

    }
  }
  ### process of sum_ebv ends
  variance_components <-asreml_variance_components(res_var_cov_h_ve = Res_Va_Ve_H2_COV_COR)

  if(inherits(variance_components, "list")){
    covariance <-  variance_components[["Covariance"]]
    correlation <-  variance_components[["Correlation"]]
    variance_components <-  variance_components[["Variance_components"]]


} else{
  if(inherits(variance_components, "data.frame")){
    covariance <-  NULL
    correlation <-  NULL

  }
}

#######################################################
  if(!"geno_model_ready" %in%names(m_matrix_model_ready_list)){
    if(length(m_matrix_model_ready_list)>1){
      names(m_matrix_model_ready_list) <- paste(paste("Omics", seq_along(m_matrix_model_ready_list), sep = ""), "model_ready", sep = "_")
      names(coefficients_list) <- paste(paste("Omics", seq_along(coefficients_list), sep = ""), "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- paste(paste("Omics", seq_along(estimated_breeding_value_list), sep = ""), "estimated_breeding_value", sep = "_")
    } else {
      names(m_matrix_model_ready_list) <- paste("Omics", "model_ready", sep = "_")
      names(coefficients_list) <- "Coefficient" # paste("Omics", "coefficient", sep = "_")
      names(estimated_breeding_value_list) <- "Estimated_breeding_value" # paste("Omics", "estimated_breeding_value", sep = "_")

    }
  } else if(length(grep("omic", names(m_matrix_model_ready_list)))==0 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
    geno_index <- grep("geno", names(m_matrix_model_ready_list))
    names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
    names(coefficients_list)[geno_index] <- "Coefficient"  #paste("Geno", "coefficient", sep = "_")
    names(estimated_breeding_value_list)[geno_index] <- "Estimated_breeding_value"  # paste("Geno", "estimated_breeding_value", sep = "_")

  } else{
    if(length(grep("omic", names(m_matrix_model_ready_list)))>=1 & "geno_model_ready" %in%names(m_matrix_model_ready_list)) {
      omic_index <- grep("omic", names(m_matrix_model_ready_list))
      geno_index <- grep("geno", names(m_matrix_model_ready_list))
      if(length(omic_index)>1){
        names(m_matrix_model_ready_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste(paste("Omics", seq_along(omic_index), sep = ""), "estimated_breeding_value", sep = "_")
      } else {
        names(m_matrix_model_ready_list)[omic_index] <- paste("Omics", "model_ready", sep = "_")
        names(coefficients_list)[omic_index] <- paste("Omics", "coefficient", sep = "_")
        names(estimated_breeding_value_list)[omic_index] <- paste("Omics", "estimated_breeding_value", sep = "_")

      }
      ####
      names(m_matrix_model_ready_list)[geno_index] <- paste("Geno", "model_ready", sep = "_")
      names(coefficients_list)[geno_index] <- paste("Geno", "coefficient", sep = "_")
      names(estimated_breeding_value_list)[geno_index] <- paste("Geno", "estimated_breeding_value", sep = "_")

    }

  }
  ######
  if(!is.null(print_lable)){

    print_lable <- unlist(print_lable)

    if(!is.null(geno_data)){
      print_lable <- c("Genomic", print_lable)
    }
    if(length(print_lable)>length(m_matrix_model_ready_list) | length(print_lable)<length(m_matrix_model_ready_list)){
      message(insight::print_color(paste(msg,paste("More than two Omics lable were provided. Default name was applied.\n")), "blue"))



    } else {
      if(length(print_lable)==length(m_matrix_model_ready_list)){
        rownames(variance_components)[1:length(print_lable)] <- paste(print_lable, "variance", sep="_")
        ######
        names(coefficients_list) <-  paste(print_lable, "coefficient", sep="_")
        names(estimated_breeding_value_list) <- paste(print_lable, "estimated_breeding_value", sep="_")
        #names(Res$Genetic_variance) = c("Genomic_variance",  paste(print_lable, "variance", sep="_"))
        names(m_matrix_model_ready_list) <- paste(print_lable, "model_ready", sep="_")
        ####
      }
    }

  } else {
    if(is.null(gmatrix)){
    message(insight::print_color(paste(msg,paste("Omics lable was not provided. Default name was applied.\n")), "blue"))
}

  }
#####################################################################
if(exists("sum_ebv")) sum_ebv <- as.data.frame(sum_ebv)
  if(is.null(Zg)){

    res <- list(Coefficients = coefficients_list,
                Asreml_model = mod,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Predicted_value =  predicted_value,
                Residual_value = residual_value,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list
    )

  } else {
    if(!is.null(Zg) & !is.null(result_met)){
      if(is.null(correlation) & is.null(covariance)){
    res <- list(Coefficients = coefficients_list,
                Asreml_model = mod,
                Estimated_breeding_value = estimated_breeding_value_list,
                Total_estimated_breeding_value = sum_ebv,
                Total_Predicted_value =  predicted_value,
                Predicted_value = across_env_predicted_value,
                Residual_value = residual_value,
                Residual_value_MET = residual_value_met,
                Variance_components = variance_components,
                M_matrix_model_ready =  m_matrix_model_ready_list)

      } else {
        if(!is.null(correlation) & !is.null(covariance)){
          res <- list(Coefficients = coefficients_list,
                      Asreml_model = mod,
                      Estimated_breeding_value = estimated_breeding_value_list,
                      Total_estimated_breeding_value = sum_ebv,
                      Total_Predicted_value =  predicted_value,
                      Predicted_value = across_env_predicted_value,
                      Correlation = correlation,
                      Covariance = covariance,
                      Residual_value = residual_value,
                      Residual_value_MET = residual_value_met,
                      Variance_components = variance_components,
                      M_matrix_model_ready =  m_matrix_model_ready_list)

        }
      }

    } else{
      if(!is.null(Zg) & is.null(result_met)){
        if(is.null(correlation) & is.null(covariance)){
        res <- list(Coefficients = coefficients_list,
                    Asreml_model = mod,
                    Estimated_breeding_value = estimated_breeding_value_list,
                    Total_estimated_breeding_value = sum_ebv,
                    Predicted_value =  predicted_value,
                    Residual_value = residual_value,
                    Variance_components = variance_components,
                    M_matrix_model_ready =  m_matrix_model_ready_list)

        } else{
          if(!is.null(correlation) & !is.null(covariance)){
            res <- list(Coefficients = coefficients_list,
                        Asreml_model = mod,
                        Estimated_breeding_value = estimated_breeding_value_list,
                        Total_estimated_breeding_value = sum_ebv,
                        Predicted_value =  predicted_value,
                        Residual_value = residual_value,
                        Correlation = correlation,
                        Covariance = covariance,
                        Variance_components = variance_components,
                        M_matrix_model_ready =  m_matrix_model_ready_list)
          }

        }
      }
    }
  }

  remove_from_global(names_in_inv_list)

    return(res)
}
