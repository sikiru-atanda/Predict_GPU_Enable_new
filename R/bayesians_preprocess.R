

#' Title
#'
#' @param nIter
#' @param burnIn
#' @param thin
#'
#' @return
#' @export
#'
#' @examples

## This function check if the following parameters are provided by the user,
# if not assign default values are provided,however the user is aware

##Check
# 1. Check if user provide number of iteration. If not default value is provided.
#    More so, it check if the value provided by the user is too small, it alert the user
#    on the implication and suggest optimal value.
# 2. Check if user provide number of burning. If not default value is provided.
#    More so, it check if the value provided by the user is too small, it alert the user
#    on the implication and suggest optimal value.
# 3. Check if user provide number of thinning. If not default value is provided.
#    More so, it check if the value provided by the user is too small, it alert the user
#    on the implication and suggest optimal value.
#' Title
#'
#' @param nIter the number of iteration for MCMC chain
#' @param burnIn  the number of burning
#' @param thin    the number of thinning
#' @param message
#'
#' @return
#' @export
#'
#' @examples
bayes_parameter_check <- function(
         nIter = NULL,
         burnIn = NULL,
         thin = NULL,
         message = TRUE
         )
  {

  msg <- sprintf("==================================================\n")

  #if(is.null(nIter) || nIter< 16000){
  if(is.null(nIter) ){
    if(isTRUE(message)){

      message(insight::print_color(paste(msg,paste("Number of iteration is missing. Default value of 26000 was assigned. \n Check if this appropriate for your data.")), "blue"))

    }
      nIter <- 200 # 26000

  } else {
    nIter< 16000
    message(paste( insight::print_color("WARNINGS\n", "red"),
                   insight::print_color(paste(msg,paste("Number of iteration is provided is less than 16000 which we consider optimal. \n Check if this appropriate for your data.")), "red")))

  }

  #if(is.null(burnIn) || burnIn < 5000){
  if(is.null(burnIn)){

    if(isTRUE(message)){

      message(insight::print_color(paste(msg,paste("Number of burn-in is missing. Default value of 1600 was assigned. \n Check if this appropriate for your data.")), "blue"))

    }
      burnIn <- 50 # 5000

  } else {

    burnIn< 1600
    message(paste( insight::print_color("WARNINGS\n", "red"),
                   insight::print_color(paste(msg,paste("Number of burn-in provided is less than 1600 which we consider optimal. \n Check if this appropriate for your data.")), "red")))


  }

  #if(is.null(thin) || thin<10){
  if(is.null(thin)){
    if(message){

      message(insight::print_color(paste(msg,paste("Number of thining is missing. Default value of 10 was assigned. \n Check if this appropriate for your data.")), "blue"))

    }
      thin <- 1 # 10

  }

  bayes_para = list(nIter = nIter, burnIn = burnIn, thin = thin)

  names(bayes_para) = c("nIter", "burnIn", "thin")

  class(bayes_para) <- "Bayesian Parameters"
  attr(bayes_para, which = 'complete_parameters') <- TRUE


  return(bayes_para)
}

#' Title
#'
#' @param fixed
#'
#' @return
#' @export
#'
#' @examples

### Adjust for fixed terms supply by the users. It is expected the user provide it
## in formula with each term separated by +.
## The output will be vector of each term
#' Title
#'
#' @param fixed fixed effect, only formula is accepted
#' @param pheno_data phenotypic data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fixed_terms <- function(
    fixed  = NULL,
    pheno_data = NULL,
    ...){

  ## Check details for the rand_fix_check.
  ## In General it check all conditions for required for fixed and random terms
  ## are fulfilled. It will have to pass through this check and ensure it pass the attribute
  fixed <- rand_fix_check(rand_fix_term = fixed,
                           pheno_data= pheno_data)

  if(attr(fixed, "cleared")!="pass" & class(fixed)!="forumla") {

    msg <- sprintf("==================================================\n")
    stop(print(paste(msg,paste(pheno_data, 'is not class formula.', sep = ""))), call. = FALSE)
  }
  # msg <- sprintf("==================================================\n")
  # ### Check if the random term provided are present in the phenotypic data
  #
  # if(!all(sapply(all.vars(fixed), function(x, pheno_data) x%in%names(pheno_data),  pheno_data))) {
  #   stop(print(paste(msg,"All variables indicated in argument 'random' should be present in phenotypic data")), call. = FALSE)
  # }

   fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
   ## Remove any space in the terms
   fixed_term = gsub(" ", "", fixed_term)

   return(fixed_term)

}


## Adjust for the interaction terms. Here single environment was in mind.
##
#'
#' @param random random effect, it must be in formula
#' @param pheno_data phenotypic data.
#' @param ..
#'
#' @return
#' @export
#'
#' @examples
random_terms <- function(random = NULL,
                         pheno_data = NULL,
                         ...){

  msg <- sprintf("==================================================\n")

  random <- rand_fix_check(rand_fix_term = random,
                           pheno_data= pheno_data)

  if(attr(random, "cleared")!="pass" && class(random)!="forumla") {

    stop(print(paste(msg,paste('the random term', 'is not formular.', sep = ""))), call. = FALSE)
  }

  # msg <- sprintf("==================================================\n")
  # ### Check if the random term provided are present in the phenotypic data
  # if(!all(sapply(all.vars(random), function(x, pheno_data) x%in%names(pheno_data),  pheno_data))) {
  #   stop(print(paste(msg,"All variables indicated in argument 'random' should be present in phenotypic data")), call. = FALSE)
  # }
  rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  ## Remove any space in the terms
  rand_term = gsub(" ", "", rand_term)

  return(rand_term)
}


#' Title
#'
#' @param fixed_term all the fixed term, which is an output from fixed_terms function
#' @param fixed_term_model_bayesian default is FIXED
#'
#' @return
#' @export
#'
#' @examples

### Check if user provide model for fixed term. Which is typically FIXED
fixed_term_model <- function(fixed_term = NULL,
                             fixed_term_model_bayesian = NULL,
                             message = TRUE){

  msg <- sprintf("==================================================\n")
  if(is.null(fixed_term_model_bayesian)){
    if(isTRUE(message)){
    warning(paste(msg, "The model for fixed term(s) is missing. We fix it for you"))
}
    fixed_term_model_bayesian = 'FIXED'
  } else{

    if(fixed_term_model_bayesian!='FIXED'){
      if(isTRUE(message)){
      warning(paste(msg, "The model for fixed term(s) should be equal to FIXED. We fix it for you"))
}
    }

    if(length(fixed_term_model_bayesian)>1){

      fixed_term_model_bayesian = 'FIXED'
    }
  }

  return(fixed_term_model_bayesian)
}





## Check the random term provided by the users and provide one in scenario where
## user provide one or less than the random terms.
# Make adjustment if more than the random terms. User should be let aware of the
# implications
#' Title
#'
#' @param rand_terms all the random terms in the formula. This is an output from random_terms function
#' @param GS_model model defined by the user. It is expected this will be for the genotype effect estimate
#' @param gen_name Column name containing the individuals/genotype in the phenotypic data file
#' @param rand_terms_model_bayesian model for all the random terms including the genotype effect.
#'                                  the total number of the model must be equal the number of random term.
#' @param message If the user want to print out the message/warning info or not
#'
#' @return
#' @export
#'
#' @examples
random_term_model <- function(rand_terms = NULL,
                              GS_model = NULL,
                              gen_name = NULL,
                              rand_terms_model_bayesian = NULL,
                              message = TRUE){


  msg <- sprintf("==================================================\n")

  #####
  #### Check for Interaction and and non-interaction term
  ## No interaction term
  rand_terms_No_Inter <-  rand_terms[!grepl(":", rand_terms)]

  ## Interaction term
  Check_rand_Inter <-  rand_terms[grepl(":", rand_terms)]

  ## Start with the No interaction terms
  if(length(rand_terms_No_Inter)!=0){
    ## Check if gen_name is present and store the position
    Gen_pos_mod =  match(gen_name, rand_terms)
    if(length(Gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect sik is missing"))), call. = FALSE)}
    if(length(Gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
    ## Extract other terms from the rand_terms_No_Inter  expect the gen_name
    rand_terms_No_Inter_No_Gen = rand_terms_No_Inter[!rand_terms_No_Inter%in%gen_name]
    ## Get the position of other terms (No interaction) in the random that is not gen_name
    Non_Gen_pos_mod = match(rand_terms_No_Inter_No_Gen, rand_terms)

  }

  ######## Initialize step For model adjustment for the random term with interaction

  if(length(Check_rand_Inter)!=0){

    TestPresentofGeno = grep(gen_name, Check_rand_Inter, value = TRUE)

    if(length(TestPresentofGeno)!=0 | !is.na(TestPresentofGeno)){

      Inter_Gen_pos_mod = match(TestPresentofGeno, rand_terms)

      Non_Gen_Inter_Test =  Check_rand_Inter[!Check_rand_Inter%in%TestPresentofGeno]

    } else {

      Non_Gen_Inter_Test = NULL

      Inter_Gen_pos_mod = NULL
    }

    if(!is.null(Non_Gen_Inter_Test)){

      Inter_Non_Gen_pos_mod = match(Non_Gen_Inter_Test, rand_terms)

    }

  } ## End

  ### Each random term will have a specific model for parameter estimate.
  ## It is expected the user will provide model for each random term,
  ## In a scenario where the model provided is not equal to the number of random terms
  ## or no model was provided at all. This next line of code take care of it.

  ### If user provide only the GS_model and random terms model is missing
  ### Especially when there is more than one random term

  if(is.null(rand_terms_model_bayesian) & !is.null(GS_model)){
    #warning(paste(msg, "The model for ranom term(s) is missing. We fix it for you"))

    ## Check to see which GS_model(s) the user provide
    #if(!is.null(GS_model)){
    ### Here the GS_model is of length rand_terms
    if(length(GS_model)==length(rand_terms)){
      ### Check GS_model to be consistent with models present in the engine
      mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_GS_model)==length(rand_terms)){
        rand_terms_model_bayesian = GS_model
        #rand_mod_copy = GS_model
      }
      ### if GS_model is less than the length of rand_terms
    } else{
      mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_GS_model)!=0){
        #rand_mod_copy = mod_present_in_GS_model
        if(message){
          warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
        }
        mod_len = length(rand_terms) - length(mod_present_in_GS_model)

        rand_terms_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
        rand_terms_model_bayesian[Gen_pos_mod] <- GS_model[1]

        ## For other terms in random effect aside gen_name
        if(exists("NonNon_Gen_pos_mod")){

          rand_terms_model_bayesian[Non_Gen_pos_mod] = "BRR"

        } ## end

        #### Random interaction terms
        if(length(Check_rand_Inter)!=0){
          ### For gen_name part of interaction
          if(!is.na(Inter_Gen_pos_mod) | length(Inter_Gen_pos_mod)!=0){

            ### This account for multi-kernel
            rand_terms_model_bayesian[Inter_Gen_pos_mod] = "RKHS"


          }

          ### For Non gen_name part of interaction
          if(length(Non_Gen_Inter_Test)!=0){

            rand_terms_model_bayesian[Non_Gen_Inter_Test] = "BRR"

          }

        } ## End of random interaction terms


      } else {
        stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)
      }
    }

    #}

  } else if (!is.null(rand_terms_model_bayesian) & is.null(GS_model)){
    # if(length(rand_terms_model_bayesian)!= length(rand_terms)){
    #   warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
    # }
    ### This will be used only for the genetic effect assuming the user only provide that
    if(length(rand_terms_model_bayesian)!= length(rand_terms)){
      if(isTRUE(message)){
        warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
      }
      ### Check GS_model to be consistent with models present in the engine
      mod_present_in_model_bayesian = rand_terms_model_bayesian[rand_terms_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_model_bayesian)==0){

        stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)

      } else {

        mod_len = (length(rand_terms) - length(mod_present_in_model_bayesian))
        rand_terms_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
        #rand_mod_copy = rand_terms_model_bayesian

      }

      rand_terms_model_bayesian[Gen_pos_mod] <- GS_model[1]

      ## For other terms in random effect aside gen_name
      if(exists("NonNon_Gen_pos_mod")){

        rand_terms_model_bayesian[Non_Gen_pos_mod] = "BRR"

      } ## end

      #### Random interaction terms
      if(length(Check_rand_Inter)!=0){
        ### For gen_name part of interaction
        if(!is.na(Inter_Gen_pos_mod) | length(Inter_Gen_pos_mod)!=0){

          ### This account for multi-kernel
          rand_terms_model_bayesian[Inter_Gen_pos_mod] = "RKHS"


        }

        ### For Non gen_name part of interaction
        if(length(Non_Gen_Inter_Test)!=0){

          rand_terms_model_bayesian[Non_Gen_Inter_Test] = "BRR"

        }

      } ## End of random interaction terms

    }

  }else {

    ### When user provide provide both GS_model and rand_terms_model_bayesian
    if(!is.null(rand_terms_model_bayesian) & !is.null(GS_model)){
 ### Check which GS_model the user supply, it has it to match the available models
      mod_present_in_GS_model <- GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_GS_model)==0){
        stop(print(paste(msg,'Provided appropiate name for the baysian model in GS_model.')), call. = FALSE)
      }
      ### Check which rand_terms_model_bayesian the user supply, it has it to match the available models
      mod_present_in_model_bayesian <- rand_terms_model_bayesian[rand_terms_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_model_bayesian)==0){
        stop(print(paste(msg,'Provided appropiate name for the baysian model in rand_terms_model_bayesian.')), call. = FALSE)
      }
      ### Check if the models term define is greater than the random term
      if(length(rand_terms_model_bayesian)>length(rand_terms)){
        stop(print(paste(msg,'The number of model is greater than the random terms')), call. = FALSE)
      }
      ##
      if(length(mod_present_in_GS_model)!=0 & length(mod_present_in_model_bayesian)!=0){

        ## Combined all the models
        rand_terms_model_bayesian  = c(GS_model, rand_terms_model_bayesian)

        if(length(rand_terms_model_bayesian)!= length(rand_terms)){

          if(length(rand_terms_model_bayesian)>length(rand_terms)){
            #stop(print(paste(msg,'The number of model is greater than the random terms')), call. = FALSE)
            rand_terms_model_bayesian = rand_terms_model_bayesian[length(rand_terms)]
            }

          if(length(rand_terms_model_bayesian)<length(rand_terms)){
          mod_len = (length(rand_terms) - length(rand_terms_model_bayesian))
          rand_terms_model_bayesian = c(rand_terms_model_bayesian,  rep("BRR",  mod_len))

          }

        }

        rand_terms_model_bayesian[Gen_pos_mod] = GS_model[1]
      }
      ###
      ###  rand_terms_model_bayesian[Gen_pos_mod] = GS_model[1]
      ## For other terms in random effect aside gen_name
      if(exists("NonNon_Gen_pos_mod")){

        rand_terms_model_bayesian[Non_Gen_pos_mod] = "BRR"

      } ## end

      #### Random interaction terms
      if(length(Check_rand_Inter)!=0){
        ### For gen_name part of interaction
        if(!is.na(Inter_Gen_pos_mod) | length(Inter_Gen_pos_mod)!=0){

          ### This account for multi-kernel
          rand_terms_model_bayesian[Inter_Gen_pos_mod] = "RKHS"


        }

        ### For Non gen_name part of interaction
        if(length(Non_Gen_Inter_Test)!=0){

          rand_terms_model_bayesian[Non_Gen_Inter_Test] = "BRR"

        }

      } ## End of random interaction terms



    } ## Sik


  } #ENd
  ####




  #} ## End else
  return(rand_terms_model_bayesian)

}


# ## Check the random term provided by the users and provide one in scenario where
# ## user provide one or less than the random terms.
# # Make adjustment if more than the random terms. User should be let aware of the
# # implications
# random_term_model <- function(rand_terms = NULL,
#                               GS_model = NULL,
#                               gen_name = NULL,
#                               rand_term_model_bayesian = NULL,
#                               message = TRUE)
#   {
#
#   msg <- sprintf("==================================================\n")
#   ### Each random term will have a specific model for parameter estimate.
#   ## It is expected the user will provide model for each random term,
#   ## In a scenario where the model provided is not equal to the number of random terms
#   ## or no model was provided at all. This next line of code take care of it.
#
#   ### If user provide only the GS_model and random terms model is missing
#   ### Especially when there is more than one random term
#   if(is.null(rand_term_model_bayesian) && !is.null(GS_model)){
#     #warning(paste(msg, "The model for ranom term(s) is missing. We fix it for you"))
#
#     ## Check to see which GS_model(s) the user provide
#     #if(!is.null(GS_model)){
#     ### Here the GS_model is of length rand_terms
#     if(length(GS_model)==length(rand_terms)){
#       ### Check GS_model to be consistent with models present in the engine
#       mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
#       if(length(mod_present_in_GS_model)==length(rand_terms)){
#         rand_term_model_bayesian = GS_model
#         #rand_mod_copy = GS_model
#       }
#       ### if GS_model is less than the length of rand_terms
#     } else{
#       mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
#       if(length(mod_present_in_GS_model)!=0){
#         #rand_mod_copy = mod_present_in_GS_model
#         if(message){
#         warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
#         }
#          mod_len = length(rand_terms) - length(mod_present_in_GS_model)
#
#         rand_term_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
#       } else {
#         stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)
#       }
#     }
#
#     #}
#
#   } else {
#     ### If user provide only rand_term_model_bayesian and GS_model is missing
#     ### Especially when there is more than one random term
#     if(!is.null(rand_term_model_bayesian) && is.null(GS_model)){
#       # if(length(rand_term_model_bayesian)!= length(rand_terms)){
#       #   warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
#       # }
#       ### This will be used only for the genetic effect assuming the user only provide that
#       if(length(rand_term_model_bayesian)!= length(rand_terms)){
#         if(message){
#         warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
#         }
#           ### Check GS_model to be consistent with models present in the engine
#         mod_present_in_model_bayesian = rand_term_model_bayesian[rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
#         if(length(mod_present_in_model_bayesian)==0){ stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)
#
#         } else {
#
#           mod_len = length(rand_terms) - length(mod_present_in_model_bayesian)
#           rand_term_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
#           #rand_mod_copy = rand_term_model_bayesian
#
#         }
#
#       }
#
#     }
#
#     ### If user provide GS_model and rand_term_model_bayesian
#     if(!is.null(rand_term_model_bayesian) && !is.null(GS_model)){
#       mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
#       if(length(mod_present_in_GS_model)==0){stop(print(paste(msg,'Provided appropiate name for the baysian model in GS_model.')), call. = FALSE)}
#
#       mod_present_in_model_bayesian = rand_term_model_bayesian[rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
#       if(length(mod_present_in_model_bayesian)==0){ stop(print(paste(msg,'Provided appropiate name for the baysian model in rand_term_model_bayesian.')), call. = FALSE)}
#
#       if(length(mod_present_in_GS_model)!=0 && length(mod_present_in_model_bayesian)!=0){
#
#         rand_term_model_bayesian  = c(GS_model, rand_term_model_bayesian)
#
#         if(length(rand_term_model_bayesian)!= length(rand_terms)){
#
#           if(length(rand_term_model_bayesian)>length(rand_terms)){
#             stop(print(paste(msg,'The number of model is greater than the random terms')), call. = FALSE)
#           }
#
#           mod_len = length(rand_terms) - length(rand_term_model_bayesian)
#           rand_term_model_bayesian = c(rand_term_model_bayesian,  rep("BRR",  mod_len))
#
#         }
#       }
#
#
#     }
#
#     ### if user didn't specific any GS_model the default rrBLUP will be used
#     if(is.null(rand_term_model_bayesian) && is.null(GS_model)){
#       rand_term_model_bayesian = array(dim = length(rand_terms), "BRR")
#     }
#
#
#   }
#
#   ###
#   # For Didier processing the first element of the vector will be used
#   # as pointer shiny production because that is the on that represent the model
#   # for genetic effect.
#
#   return(rand_term_model_bayesian)
#
# }
#
