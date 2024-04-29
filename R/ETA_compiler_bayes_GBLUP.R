#' Compile ETA for Bayesian Genomic Prediction Models
#'
#' This function compiles the ETA components for Bayesian genomic prediction models,
#' including RKHS, and GBLUP_BRR, based on the provided fixed and random model terms,
#' genomic or other omics relationship matrices, and phenotypic data.
#'
#' @param fixed A formula or a list of formulas specifying the fixed effects.
#' @param random A formula or a list of formulas specifying the random effects.
#' @param GS_model A character string specifying the genomic selection model to use.
#' @param fixed_term_model_bayesian Character string specifying the model for fixed terms. Default is NULL.
#' @param rand_term_model_bayesian Character vector specifying the models for each random term. Default is NULL.
#' @param pheno_data A data.frame containing the phenotypic data.
#' @param gmatrix A numeric matrix representing the genomic relationship matrix.
#' @param gkernel A numeric matrix representing a genomic kernel for RKHS models. Default is NULL.
#' @param omic1_kernel A numeric matrix representing an omics-based kernel. Default is NULL.
#' @param omic2_kernel Same as `omic1_kernel`. Default is NULL.
#' @param omic3_kernel Same as `omic1_kernel`. Default is NULL.
#' @param gen_name A character string specifying the column name in `pheno_data` that contains the genotype identifiers.
#' @param heter_groups A character string specifying the column name in `pheno_data` for heterogeneous groups. Default is NULL.
#' @return A list containing the compiled ETA components, the modified phenotypic data, and names of ETA elements.
#' @examples
#' # Assuming pheno_data is your phenotypic dataset, gmatrix is the genomic relationship matrix:
#' result <- ETA_compiler_bayes_GBLUP(fixed = ~ fixed_effect,
#'                                    random = ~ random_effect,
#'                                    GS_model = "RKHS",
#'                                    pheno_data = pheno_data,
#'                                    gmatrix = gmatrix,
#'                                    gen_name = "GenotypeID")
#' @export

ETA_compiler_bayes_GBLUP <- function(
    fixed = NULL,
    random = NULL,
    GS_model = NULL,
    fixed_term_model_bayesian = NULL,
    rand_term_model_bayesian = NULL,
    pheno_data = NULL,
    gmatrix= NULL,
    gkernel = NULL,
    omic1_kernel = NULL,
    omic2_kernel = NULL,
    omic3_kernel = NULL,
    gen_name = NULL,
    heter_groups = NULL,
    ...
){
#browser()
  ### Create empty list for ETA

  #rm(ZE, ZEZE, Zg, K1, K2, ETA)

  non_gen_inter_test <-  NULL

  inter_gen_pos_mod <-  NULL

  non_gen_pos_mod <-  NULL

  ETA <- list()

  msg <- sprintf("==================================================\n")
  ### Get the random terms. Both no interaction and interaction terms if present in the random terms
  rand_terms <- random_terms(random = random,
                             pheno_data = pheno_data)

  ### Assign
  rand_model <- random_term_model(rand_terms = rand_terms,
                                  rand_terms_model_bayesian = rand_term_model_bayesian,
                                  GS_model = GS_model,
                                  gen_name = gen_name)

  if(!is.null(fixed)){
    ETA <-  ETA_compiler_fixed_term(fixed = fixed,
                                    pheno_data = pheno_data,
                                    fixed_term_model_bayesian = fixed_term_model_bayesian)


  }
  ################################################
  #### Check for Interaction and and non-interaction term
  ## No interaction term
  rand_terms_no_inter <- rand_terms[!grepl(":", rand_terms)]

  ## Interaction term
  check_rand_inter <- rand_terms[grepl(":", rand_terms)]

  ## Start with the No interaction terms
  if(length(rand_terms_no_inter)!=0){
    ## Check if gen_name is present and store the position
    gen_pos_mod=  match(gen_name, rand_terms)
    if(length(gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect is missing"))), call. = FALSE)}
    if(length(gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
    ## Extract other terms from the rand_terms_no_inter  expect the gen_name
    rand_terms_no_inter_no_gen <- rand_terms_no_inter[!rand_terms_no_inter%in%gen_name]
    ## Get the position of other terms (No interaction) in the random that is not gen_name
    non_gen_pos_mod <- match(rand_terms_no_inter_no_gen, rand_terms)

  } else {
    non_gen_pos_mod <-  NULL
  }

  ######## Initialize step For model adjustment for the random term with interaction

  if(length(check_rand_inter)!=0){

    test_present_of_geno <-  grep(gen_name, check_rand_inter, value = TRUE)

    if(length(test_present_of_geno)!=0 | !is.na(test_present_of_geno)){

      inter_gen_pos_mod <-  match(test_present_of_geno, rand_terms)

      non_gen_inter_test <- check_rand_inter[!check_rand_inter%in%test_present_of_geno]

    } else {

      non_gen_inter_test <-  NULL

      inter_gen_pos_mod <-  NULL
    }

    if(!is.null(non_gen_inter_test)){

      inter_non_gen_pos_mod <- match(non_gen_inter_test, rand_terms)

    }

  }

  if (length(inter_gen_pos_mod)>=1){
    if (length(pheno_data[,gen_name]) ==length(unique(pheno_data[,gen_name]))){
      stop(print(paste(msg, "Phenotypic data contain single environment but you specify multi-environment analysis.")), call. = FALSE)
    }

  }
  #########
  #### When the genotype are present in more than one environment/location
  if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
    ### incidence matrix for main eff. of the genotypes
    Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)

    if(!is.null(heter_groups)){
      ZE <- model.matrix(~factor(pheno_data[,heter_groups])-1)
      ZEZE <-tcrossprod(ZE)

    }

  } else {

    Zg <-  NULL
    ZE <-  NULL
    ZEZE <- NULL
  }

  rand_model_copy <- rand_model
  #####
  #####################################################
  for (ra in 1:length(rand_terms)){

    rand_model <- rand_model_copy[ra]

    ### start when no genotype
    if(is.null(non_gen_pos_mod)){
      if(length(non_gen_pos_mod)!=0){
        if((ra == non_gen_pos_mod & !is.null(ZE))){
          ETA[[length(ETA) + 1]] <- list(X= ZE,
                                         model=rand_model,
                                         saveEffects=TRUE)

        }

      }

    } ### End

    datasets <- list(gmatrix, omic1_kernel, omic2_kernel, omic3_kernel)
    dataset_names <- c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
    datasets_index <- which(!sapply(datasets, is.null))
    datasets <-  datasets[datasets_index]
    dataset_names <- dataset_names[datasets_index]
    ETA_element_name <- character()

    for (i in seq_along(datasets)) {
      dataset <- datasets[[i]]
      if (!is.null(dataset)) { ## this seems redundant but useful
        if(rand_model== "RKHS"){
          if((ra == gen_pos_mod ) & is.null(Zg)){
            ETA[[length(ETA) + 1]] <- list(K= as.matrix(dataset),
                                           model = rand_model,
                                           saveEffects = TRUE)
          } else if((ra == gen_pos_mod) & !is.null(Zg)){
            K1 <- Zg%*%as.matrix(dataset)%*%t(Zg)
            ETA[[length(ETA) + 1]] <- list(K= K1,
                                           model = rand_model,
                                           saveEffects = TRUE)
          } else {

            if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
              if((ra == inter_gen_pos_mod & !is.null(Zg)) & !is.null(ZEZE)){

                K1 <- Zg%*%as.matrix(dataset)%*%t(Zg)

                K2<-K1*ZEZE

                ETA[[length(ETA) + 1]] <- list(K= K2,
                                               model=rand_model,
                                               saveEffects=TRUE)

              }

            }
          }
          ETA_element_name <- c(ETA_element_name, dataset_names[i])
        } else{

          if(rand_model== "BRR"){
            if((ra == gen_pos_mod ) & is.null(Zg)){
              ETA[[length(ETA) + 1]] <- list(X= as.matrix(dataset),
                                             model = rand_model,
                                             saveEffects = TRUE)
            } else if((ra == gen_pos_mod) & !is.null(Zg)){
              K1 <- Zg%*%as.matrix(dataset)%*%t(Zg)
              ETA[[length(ETA) + 1]] <- list(X= K1,
                                             model = rand_model,
                                             saveEffects = TRUE)
            } else {

              if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
                if((ra == inter_gen_pos_mod & !is.null(Zg)) & !is.null(ZEZE)){

                  K1 <- Zg%*%as.matrix(dataset)%*%t(Zg)

                  K2<-K1*ZEZE

                  ETA[[length(ETA) + 1]] <- list(X= K2,
                                                 model=rand_model,
                                                 saveEffects=TRUE)

                }

              }
            }
            ETA_element_name <- c(ETA_element_name, dataset_names[i])
          }

        }
      }
    }
  }
  output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)

  return(output)

}


# ETA_compiler_bayes_GBLUP <- function(
#     fixed = NULL,
#     random = NULL,
#     GS_model = NULL,
#     fixed_term_model_bayesian = NULL,
#     rand_term_model_bayesian = NULL,
#     pheno_data = NULL,
#     gmatrix= NULL,
#     gkernel = NULL,
#     omic1_kernel = NULL,
#     omic2_kernel = NULL,
#     omic3_kernel = NULL,
#     gen_name = NULL,
#     heter_groups = NULL,
#     ...
# ){
#
#   ### Create empty list for ETA
#
#   #rm(ZE, ZEZE, Zg, K1, K2, ETA)
#   ETA = list()
#
#   msg <- sprintf("==================================================\n")
#   ### Get the random terms. Both no interaction and interaction terms if present in the random terms
#   rand_terms <- random_terms(random = random,
#                              object = pheno_data)
#
#   ### Assign
#   rand_model <- random_term_model(rand_terms = rand_terms,
#                                   rand_terms_model_bayesian = rand_term_model_bayesian,
#                                   GS_model = GS_model,
#                                   gen_name = gen_name)
#
#   if(!is.null(fixed)){
#     ETA <-  ETA_compiler_fixed_term(fixed = fixed,
#                                     pheno_data = pheno_data,
#                                     fixed_term_model_bayesian = fixed_term_model_bayesian)
#
#
#   }
#  ################################################
#   #### Check for Interaction and and non-interaction term
#   ## No interaction term
#   rand_terms_no_inter = rand_terms[!grepl(":", rand_terms)]
#
#   ## Interaction term
#   check_rand_inter = rand_terms[grepl(":", rand_terms)]
#
#   ## Start with the No interaction terms
#   if(length(rand_terms_no_inter)!=0){
#     ## Check if gen_name is present and store the position
#     gen_pos_mod=  match(gen_name, rand_terms)
#     if(length(gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect is missing"))), call. = FALSE)}
#     if(length(gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
#     ## Extract other terms from the rand_terms_no_inter  expect the gen_name
#     rand_terms_no_inter_no_gen = rand_terms_no_inter[!rand_terms_no_inter%in%gen_name]
#     ## Get the position of other terms (No interaction) in the random that is not gen_name
#     non_gen_pos_mod = match(rand_terms_no_inter_no_gen, rand_terms)
#
#   }
#
#   ######## Initialize step For model adjustment for the random term with interaction
#
#   if(length(check_rand_inter)!=0){
#
#     test_present_of_geno = grep(gen_name, check_rand_inter, value = TRUE)
#
#     if(length(test_present_of_geno)!=0 | !is.na(test_present_of_geno)){
#
#       inter_gen_pos_mod = match(test_present_of_geno, rand_terms)
#
#       non_gen_inter_test =  check_rand_inter[!check_rand_inter%in%test_present_of_geno]
#
#     } else {
#
#       non_gen_inter_test = NULL
#
#       inter_gen_pos_mod = NULL
#     }
#
#     if(!is.null(non_gen_inter_test)){
#
#       inter_non_gen_pos_mod = match(non_gen_inter_test, rand_terms)
#
#     }
#
#   }
#
#   #########
#   #### When the genotype are present in more than one environment/location
#   if(length(pheno_data[,gen_name])>length(unique(pheno_data[,gen_name]))){
#     ### incidence matrix for main eff. of the genotypes
#     Zg<-stats::model.matrix(~factor(pheno_data[,gen_name])-1)
#
#     if(!is.null(heter_groups)){
#       ZE <- model.matrix(~factor(pheno_data[,heter_groups])-1)
#       ZEZE<-tcrossprod(ZE)
#
#     }
#
#   }
#
#   rand_model_copy = rand_model
#   #####
#   #####################################################
#   for (ra in 1:length(rand_terms)){
#
#     rand_model = rand_model_copy[ra]
#
#
#     ### start when no genotype
#     if(exists("non_gen_pos_mod")){
#       if(length(non_gen_pos_mod)!=0){
#         if((ra == non_gen_pos_mod & exists("ZE"))){
#           ETA[[length(ETA) + 1]] <- list(X= ZE,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#         }
#
#       }
#
#     } ### End
#
#     if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod ) & !exists("Zg")){
#           ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("gmatrix")
#           } else if((ra == gen_pos_mod) & exists("Zg")){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix")
#             } else if(((ra == gen_pos_mod & exists("Zg")))){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix")
#                 }
#
#               }
#
#
#
#               }
#
#
#
#             }
#
#
#           }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel")
#                 }
#
#
#
#                 }
#
#               }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel")
#                 }
#
#                 if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                   if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel")
#                   }
#
#                 }
#
#               }
#             }
#
#             ETA_element_name = c("gkernel")
#           }
#
#
#         }
#
#       }
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(omic1_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if(((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE"))){
#
#                   K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic1_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#       ETA_element_name = c("omic1_kernel")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(omic2_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic2_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel")
#           } else if(((ra == gen_pos_mod & exists("Zg")))){
#
#             K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel")
#           } else {
#             if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel")
#             } else {
#               if((length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod))){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic2_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#
#       ETA_element_name = c("omic2_kernel")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(omic3_kernel)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic3_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic3_kernel")
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic3_kernel")
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("omic3_kernel")
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       }
#
#
#       ETA_element_name = c("omic3_kernel")
#       #If user provide geno_data and omic1_data
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#       #ETA_element_name = c("geno_data", "omic2_data")
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic2_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic2_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic3_data")
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic3_kernel")
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K=  K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic ,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic3_kernel")
#
#                 }
#
#               }
#
#             }
#
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic3_kernel")
#
#             } else {
#
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                   K2_omic <-K1_omic*ZEZE
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic3_kernel")
#                 }
#
#
#
#               }
#
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel")
#               } else {
#
#                 if((ra == gen_pos_mod & exists("Zg"))){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic3_kernel")
#                 }
#
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic3_kernel%*%t(Zg)
#                     K2_omic <-K1_omic*ZEZE
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic3_kernel")
#                   }
#
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#       } ## End
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic2_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if(((ra == gen_pos_mod & !exists("Zg")))){
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#       #ETA_element_name = c("omic1_data", "omic2_data")
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#
#       #ETA_element_name = c("omic1_data", "omic3_data")
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic2_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic2_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2<-K1*ZEZE
#
#               K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_1<-K1_1*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic2_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic2_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_1 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_1<-K1_1*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         }
#
#
#       }
#
#
#
#
#       #ETA_element_name = c("omic2_data", "omic3_data")
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic1_data", "omic2_data")
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(rand_model== "RKHS"){
#
#         if(!is.null(gkernel)){
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(gkernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#         } else {
#
#           if(!is.null(gmatrix)){
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(gmatrix),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic1_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#           }
#
#           ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#         }
#
#       } else {
#
#         if(rand_model== "BRR"){
#           if(!is.null(gkernel)){
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(gkernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#           } else {
#
#             if(!is.null(gmatrix)){
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(gmatrix),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_kernel),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_kernel),
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#             }
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#           }
#         }
#
#       }
#       #ETA_element_name = c("geno_data", "omic1_data", "omic3_data")
#       ###
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((!is.null(omic1_kernel) &  is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#     } else if((!is.null(gmatrix) | !is.null(gkernel)) & ((is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_omic<-K1_omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_3omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_omic<-K1_omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_3omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_3omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_omic<-K1_omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_3omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_3omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic2_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#       #ETA_element_name = c("geno_data", "omic2_data", "omic3_data")
#
#     } else if((is.null(gmatrix) & is.null(gkernel)) & ((!is.null(omic1_kernel) &  !is.null(omic2_kernel)) & !is.null(omic3_kernel))){
#
#
#       if(rand_model== "RKHS"){
#         if((ra == gen_pos_mod & !exists("Zg"))){
#
#           ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#         } else if((ra == gen_pos_mod & exists("Zg"))){
#
#           K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#           K1_2 <- Zg%*%omic2_kernel%*%t(Zg)
#
#           K1_3 <- Zg%*%omic3_kernel%*%t(Zg)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_2,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#           ETA[[length(ETA) + 1]] <- list(K= K1_3,
#                                          model=rand_model,
#                                          saveEffects=TRUE)
#
#
#           ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#         } else {
#           if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#             if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#               K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K2 <- K1*ZEZE
#
#               K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K2_1<- K1_1*ZEZE
#
#               K1_2 <- Zg%*%omic3_kernel%*%t(Zg)
#
#               K2_2 <- K1_2*ZEZE
#
#               ETA[[length(ETA) + 1]] <- list(K= K2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K2_2,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#             }
#
#           }
#         }
#
#
#       } else {
#
#         if(rand_model== "BRR"){
#
#           if((ra == gen_pos_mod & !exists("Zg"))){
#
#             ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_2 <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3 <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_2,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(X= K1_3,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K2 <- K1*ZEZE
#
#                 K1_1 <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1<- K1_1*ZEZE
#
#                 K1_2 <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_2 <- K1_2*ZEZE
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K2_2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("omic1_kernel", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#         }
#
#
#       }
#
#
#     } else {
#
#       if(!is.null(gmatrix)){
#         if(rand_model== "RKHS"){
#           if((ra == gen_pos_mod & !exists("Zg"))){
#             ETA[[length(ETA) + 1]] <- list(K= gmatrix,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#           } else if((ra == gen_pos_mod & exists("Zg"))){
#
#             K1 <- Zg%*%gmatrix%*%t(Zg)
#
#             K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#             K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#             K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_2omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                            model=rand_model,
#                                            saveEffects=TRUE)
#
#             ETA_element_name = c("gmatrix", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#           } else {
#             if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#               if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                 K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                 K2<-K1*ZEZE
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_1omic<-K1_omic*ZEZE
#                 ###
#                 K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K2_1omic<-K1_2omic*ZEZE
#                 ###
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 K2_2omic<-K1_3omic*ZEZE
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K1_1omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_1omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(K= K2_2omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#               }
#
#             }
#           }
#
#
#         } else {
#
#           if(rand_model== "BRR"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(X= gmatrix,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gmatrix%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_2omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gmatrix", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gmatrix%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K1_1omic<-K1_omic*ZEZE
#                   ###
#                   K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_1omic<-K1_2omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_2omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K1_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(X= K2_2omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gmatrix", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#
#           }
#
#
#         }
#
#
#       } else {
#
#         if(!is.null(gkernel)){
#
#           if(rand_model== "RKHS"){
#             if((ra == gen_pos_mod & !exists("Zg"))){
#               ETA[[length(ETA) + 1]] <- list(K= gkernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic1_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic2_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= omic3_kernel,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#             } else if((ra == gen_pos_mod & exists("Zg"))){
#
#               K1 <- Zg%*%gkernel%*%t(Zg)
#
#               K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#               K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#               K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_2omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA[[length(ETA) + 1]] <- list(K= K1_3omic,
#                                              model=rand_model,
#                                              saveEffects=TRUE)
#
#               ETA_element_name = c("gkernel", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#             } else {
#               if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                 if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                   K1 <- Zg%*%gkernel%*%t(Zg)
#
#                   K2<-K1*ZEZE
#
#                   K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                   K1_1omic<-K1_omic*ZEZE
#                   ###
#                   K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                   K2_1omic<-K1_2omic*ZEZE
#                   ###
#                   K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                   K2_2omic<-K1_3omic*ZEZE
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K1_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_1omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA[[length(ETA) + 1]] <- list(K= K2_2omic,
#                                                  model=rand_model,
#                                                  saveEffects=TRUE)
#
#                   ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                 }
#
#               }
#             }
#
#           } else {
#
#             if(rand_model== "BRR"){
#               if((ra == gen_pos_mod & !exists("Zg"))){
#                 ETA[[length(ETA) + 1]] <- list(X= gkernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic1_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic2_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= omic3_kernel,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#
#               } else if((ra == gen_pos_mod & exists("Zg"))){
#
#                 K1 <- Zg%*%gkernel%*%t(Zg)
#
#                 K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                 K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                 K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_2omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA[[length(ETA) + 1]] <- list(X= K1_3omic,
#                                                model=rand_model,
#                                                saveEffects=TRUE)
#
#                 ETA_element_name = c("gkernel", "omic1_kernel",  "omic2_kernel", "omic3_kernel")
#
#               } else {
#                 if(length(inter_gen_pos_mod)!=0 | !is.na(inter_gen_pos_mod)){
#                   if((ra == inter_gen_pos_mod & exists("Zg")) & exists("ZEZE")){
#
#                     K1 <- Zg%*%gkernel%*%t(Zg)
#
#                     K2<-K1*ZEZE
#
#                     K1_omic <- Zg%*%omic1_kernel%*%t(Zg)
#
#                     K1_1omic<-K1_omic*ZEZE
#                     ###
#                     K1_2omic <- Zg%*%omic2_kernel%*%t(Zg)
#
#                     K2_1omic<-K1_2omic*ZEZE
#                     ###
#                     K1_3omic <- Zg%*%omic3_kernel%*%t(Zg)
#
#                     K2_2omic<-K1_3omic*ZEZE
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K1_1omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_1omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA[[length(ETA) + 1]] <- list(X= K2_2omic,
#                                                    model=rand_model,
#                                                    saveEffects=TRUE)
#
#                     ETA_element_name = c("gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")
#                   }
#
#                 }
#               }
#             }
#
#
#           }
#
#
#         }
#
#       }
#
#     }
#
#   }
#
#   output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)
#
#   names(output) <- c("ETA", "pheno_data", "ETA_element_name")
#
#   return(output)
#
#
#
# }
#
#
