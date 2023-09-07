

#' Title
#' This function deal with fixed terms defined by the user and assign model to each fixed term
#'
#' @param fixed fixed terms defined by the user
#' @param fixed_term_model_bayesian model for the fixed term and the default is FIXED
#' @param pheno_data phenotypic data
#'
#' @return
#' @export
#'
#' @examples

ETA_compiler_fixed_term <- function(fixed = NULL,
                         fixed_term_model_bayesian = NULL,
                         pheno_data = NULL
                         ){


  ### Initialize steps for compiling the Fixed terms
  ## Start with creating empty list for ETA compilation

  ETA = list()

  #if(!is.null(fixed)){
    fixed_term_no_inter <- fixed_terms(fixed = fixed, pheno_data = pheno_data)
    fixed_model <- fixed_term_model(fixed_term_no_inter,
                                    fixed_term_model_bayesian)
    #### Fit Fixed terms in ETA
    if(length(fixed_term_no_inter)!=0){
      for (ET in 1:length(fixed_term_no_inter)) {


        ETA[[ET]] <- list(X=stats::model.matrix(~factor(pheno_data[, fixed_term_no_inter[ET]])-1),
                          model=fixed_model)

      }

    }
  #}

  return(ETA)

}


#' Title
#'
#' @param fixed
#' @param random
#' @param GS_model
#' @param fixed_term_model_bayesian
#' @param rand_term_model_bayesian
#' @param pheno_data
#' @param geno_data
#' @param omic1_data
#' @param omic2_data
#' @param omic3_data
#' @param gen_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

ETA_compiler_bayes <- function(
         fixed = NULL,
         random = NULL,
         GS_model = NULL,
         fixed_term_model_bayesian = NULL,
         rand_term_model_bayesian = NULL,
         pheno_data = NULL,
         geno_data = NULL,
         omic1_data = NULL,
         omic2_data = NULL,
         omic3_data = NULL,
         gen_name = NULL,
         ...
){

  ### Create empty list for ETA
  ETA = list()

  msg <- sprintf("==================================================\n")
  ### Get the random terms. Here no interaction terms in the random effect
  rand_term_no_inter <- random_terms(random = random,
                                     pheno_data = pheno_data)

  ### Assign
  rand_model <- random_term_model(rand_terms = rand_term_no_inter,
                                  gen_name = gen_name,
                                  rand_terms_model_bayesian = rand_term_model_bayesian,
                                  GS_model = GS_model,
                                  message = message)

  if(!is.null(fixed)){
   ETA <-  ETA_compiler_fixed_term(fixed = fixed,
                                   pheno_data = pheno_data,
                                   fixed_term_model_bayesian = fixed_term_model_bayesian)

   #len_ETA = length(ETA)
}
  # } else {
  #
  #   len_ETA = 0
  # }

   if(length(rand_model)==1 && length(rand_term_no_inter)==1){

     if(!is.null(geno_data) & ((is.null(omic1_data) &  is.null(omic2_data)) & is.null(omic3_data))){

     ETA[[length(ETA) + 1]] <- list(X= as.matrix(geno_data),
                                 model=rand_model,
                                 saveEffects=TRUE)

     ETA_element_name = c("geno_data")
     } else if(is.null(geno_data) & ((!is.null(omic1_data) &  is.null(omic2_data)) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic1_data),
                                  model=rand_model,
                                  saveEffects=TRUE)
      ETA_element_name = c("omic1_data")
     } else if ((is.null(geno_data) & is.null(omic1_data)) &  (!is.null(omic2_data) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic2_data),
                                  model=rand_model,
                                  saveEffects=TRUE)

      ETA_element_name = c("omic2_data")

     } else if((is.null(geno_data) & is.null(omic1_data)) &  (is.null(omic2_data) & !is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= as.matrix(omic3_data),
                                  model=rand_model,
                                  saveEffects=TRUE)

       ETA_element_name = c("omic3_data")
       #If user provide geno_data and omic1_data

     } else if ((!is.null(geno_data) & !is.null(omic1_data)) &  (is.null(omic2_data) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA_element_name = c("geno_data", "omic1_data")

     } else if((!is.null(geno_data) & !is.null(omic2_data)) &  (is.null(omic1_data) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                  model=rand_model,
                                  saveEffects=TRUE)

      ETA_element_name = c("geno_data", "omic2_data")

     } else if ((!is.null(geno_data) & is.null(omic1_data)) &  (is.null(omic2_data) & !is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                  model=rand_model,
                                  saveEffects=TRUE)

       ETA_element_name = c("geno_data", "omic3_data")

     } else if ((is.null(geno_data) & !is.null(omic1_data)) &  (!is.null(omic2_data) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                      model=rand_model,
                                      saveEffects=TRUE)


       ETA_element_name = c("omic1_data", "omic2_data")
     } else if((is.null(geno_data) & !is.null(omic1_data)) &  (is.null(omic2_data) & !is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA_element_name = c("omic1_data", "omic3_data")
     } else if((is.null(geno_data) & is.null(omic1_data)) &  (!is.null(omic2_data) & !is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA_element_name = c("omic2_data", "omic3_data")
     } else if ((((!is.null(geno_data) & !is.null(omic1_data)) &  !is.null(omic2_data)) & is.null(omic3_data))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA_element_name = c("geno_data", "omic1_data", "omic2_data")

     } else if(((((!is.null(geno_data) & !is.null(omic1_data)) &  is.null(omic2_data)) & !is.null(omic3_data)))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA_element_name = c("geno_data", "omic1_data", "omic3_data")
       ###
     } else if(((((!is.null(geno_data) & is.null(omic1_data)) &  !is.null(omic2_data)) & !is.null(omic3_data)))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                      model=rand_model,
                                      saveEffects=TRUE)


       ETA_element_name = c("geno_data", "omic2_data", "omic3_data")

       } else if(((((is.null(geno_data) & !is.null(omic1_data)) &  !is.null(omic2_data)) & !is.null(omic3_data)))){

         ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                    model=rand_model,
                                    saveEffects=TRUE)


         ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                        model=rand_model,
                                        saveEffects=TRUE)

         ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                        model=rand_model,
                                        saveEffects=TRUE)


         ETA_element_name = c("omic1_data", "omic2_data", "omic3_data")


     } else {

       if(((((!is.null(geno_data) & !is.null(omic1_data)) &  !is.null(omic2_data)) & !is.null(omic3_data)))){

       ETA[[length(ETA) + 1]] <- list(X= geno_data,
                                  model=rand_model,
                                  saveEffects=TRUE)


       ETA[[length(ETA) + 1]] <- list(X= omic1_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA[[length(ETA) + 1]] <- list(X= omic2_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       ETA[[length(ETA) + 1]] <- list(X= omic3_data,
                                      model=rand_model,
                                      saveEffects=TRUE)

       }

       ETA_element_name = c("geno_data", "omic1_data", "omic2_data", "omic3_data")
     }

   } else {

     stop(paste(msg, 'This works only for one random effect'))


   }

  output <- list(ETA= ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)

  names(output) <- c("ETA", "pheno_data", "ETA_element_name")

  return(output)



}


