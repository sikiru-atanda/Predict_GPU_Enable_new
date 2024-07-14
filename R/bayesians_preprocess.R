
#' Check and Set Default Parameters for Bayesian Model Execution
#'
#' Validates and assigns default values to the key parameters of a Bayesian model:
#' number of iterations (`nIter`), burn-in period (`burnIn`), and thinning interval (`thin`).
#' It optionally displays messages for missing parameters or parameters set below recommended thresholds.
#'
#' @param nIter An optional integer specifying the total number of MCMC iterations to perform.
#' Default is set to 200 if not provided, but a message suggests 26000 as optimal.
#' @param burnIn An optional integer specifying the number of initial MCMC iterations to discard (burn-in).
#' Default is set to 50 if not provided, but a message suggests 1600 as optimal.
#' @param thin An optional integer specifying the thinning interval for MCMC sampling.
#' Default is set to 1 if not provided, but a message suggests 10 as optimal.
#' @param message A logical indicating whether to print messages about missing parameters or those set below the recommended thresholds.
#'
#' @return A list of class "Bayesian Parameters" with named elements `nIter`, `burnIn`, and `thin`,
#' each containing the respective parameter value. The list also has an attribute 'complete_parameters' set to TRUE.
#'
#' @details
#' This function is particularly useful in preparing Bayesian model parameters, ensuring that they
#' meet minimum recommended values for adequate model convergence and accuracy. It helps users identify
#' and rectify potentially inadequate parameter settings before model execution.
#'
#' @examples
#' bayes_params <- bayes_parameter_check(
#'   nIter = 25000,
#'   burnIn = 5000,
#'   thin = 5,
#'   message = TRUE
#' )
#' print(bayes_params)
#'
#' @export
#'
bayes_parameter_check <- function(
         nIter = NULL,
         burnIn = NULL,
         thin = NULL,
         message = TRUE
         )
  {

  msg <- "\n ==================================================\n"

  if(is.null(nIter) ){
    if(isTRUE(message)){

      message(insight::print_color(paste(msg,paste("\n Number of iteration is missing. Default value of 26000 was assigned. \n Check if this appropriate for your data.\n ")), "blue"))

    }
      nIter <- 200 # 26000

  } else {

    if(nIter< 16000){
    message(paste( insight::print_color("\nWARNINGS\n", "red"),
                   insight::print_color(paste(msg,paste("Number of iteration is provided is less than 16000 which we consider optimal. \n Check if this appropriate for your data.")), "red")))

    }
  }

  if(is.null(burnIn)){

    if(isTRUE(message)){

      message(insight::print_color(paste(msg,paste("Number of burn-in is missing. Default value of 1600 was assigned. \n Check if this appropriate for your data.")), "blue"))

    }
      burnIn <- 50 # 5000

  } else {

    if(burnIn< 1600){
    message(paste( insight::print_color("\n WARNINGS\n", "red"),
                   insight::print_color(paste(msg,paste("Number of burn-in provided is less than 1600 which we consider optimal. \n Check if this appropriate for your data.\n ")), "red")))

    }

  }

  if(is.null(thin)){
    if(isTRUE(message)){

      message(insight::print_color(paste(msg,paste("\n Number of thining is missing. Default value of 10 was assigned. \n Check if this appropriate for your data.\n ")), "blue"))

    }
      thin <- 1 # 10

  }

  bayes_para = list(nIter = nIter, burnIn = burnIn, thin = thin)

  names(bayes_para) = c("nIter", "burnIn", "thin")

  class(bayes_para) <- "Bayesian Parameters"
  attr(bayes_para, which = 'complete_parameters') <- TRUE


  return(bayes_para)
}


#' Process and Validate Fixed Effect Terms for Modeling
#'
#' This function checks and validates the fixed effect terms specified for a statistical model,
#' ensuring they conform to required conditions. It leverages an auxiliary function `rand_fix_check`
#' to perform the validation and extracts the terms from a formula object.
#'
#' @param fixed A formula specifying the fixed effects to be included in the model.
#' @param pheno_data A data frame containing the phenotypic data, which is used by the `rand_fix_check` function to validate the fixed terms.
#' @param ... Additional arguments passed to the `rand_fix_check` function.
#'
#' @return A character vector of fixed effect terms that have been validated and extracted from the input formula.
#'
#' @details
#' The function is designed to preprocess and validate fixed effect terms before their use in statistical models,
#' particularly those involving genetic or phenotypic data analysis. It ensures that the input terms meet the
#' necessary conditions for modeling and are properly formatted. The validation process involves checking that
#' the fixed effects are appropriately specified as a formula and pass through a custom check (`rand_fix_check`),
#' which may involve various validations such as ensuring the presence of the terms in the provided phenotypic data.
#'
#' @examples
#' # Assuming pheno_data is your dataset and you have a formula for fixed effects:
#' fixed_effects_formula <- ~ Trait1 + Trait2
#' fixed_terms <- fixed_terms(fixed = fixed_effects_formula, pheno_data = pheno_data)
#' print(fixed_terms)
#'
#' @seealso \code{\link{rand_fix_check}} for details on the validation process.
#' @export
#'
fixed_terms <- function(
    fixed  = NULL,
    pheno_data = NULL,
    ...){

  ## Check details for the rand_fix_check.
  ## In General it check all conditions for required for fixed and random terms
  ## are fulfilled. It will have to pass through this check and ensure it pass the attribute
  fixed <- rand_fix_check(rand_fix_term = fixed,
                          term_type = "fixed",
                           pheno_data= pheno_data)

  if(attr(fixed, "cleared")!="pass" & !inherits(fixed, "formula")) {

    msg <- "\n ==================================================\n"
    stop(msg, "The fixed term is not a class of type 'formula'. Example: ", "fixed = ~ X + Y")
  }

   fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
   ## Remove any space in the terms
   fixed_term = gsub(" ", "", fixed_term)

   return(fixed_term)

}


#' Validate and Process Random Effect Terms
#'
#' This function processes and validates random effect terms for  modeling,
#' ensuring they conform to required conditions. It leverages `rand_fix_check` for validation.
#'
#' @param random A formula specifying the random effects to be included in the model.
#' @param pheno_data A data frame containing the phenotypic data used for validation.
#' @param ... Additional arguments passed to the `rand_fix_check` function.
#'
#' @return The validated `random` term as a formula object, with an attribute "cleared" indicating successful validation.
#'
#' @details
#' Ensures the random effect terms are properly specified and present in the phenotypic data,
#' converting relevant variables to factors where necessary. This function is integral to the
#' preprocessing steps required for accurate and effective statistical analysis in genetic studies.
#'
#' @examples
#' pheno_data <- data.frame(Trait1 = rnorm(100), Trait2 = rnorm(100), Genotype = as.factor(rep(1:10, each = 10)))
#' random_effects_formula <- ~ Genotype
#' validated_random <- random_terms(random = random_effects_formula, pheno_data = pheno_data)
#' print(validated_random)
#'
#' @export
random_terms <- function(random = NULL,
                         pheno_data = NULL,
                         ...){

  msg <- "\n ==================================================\n"

  random <- rand_fix_check(rand_fix_term = random,
                           term_type  = "random",
                           pheno_data= pheno_data)

  if(attr(random, "cleared")!="pass" && !inherits(random, "formula")) {

    stop(msg, "\n The random term is not a class of type 'formula'. Example: ", "fixed = ~ X + Y \n")
  }

  rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  ## Remove any space in the terms
  rand_term = gsub(" ", "", rand_term)

  return(rand_term)
}


#' Validate and Set Model Type for Fixed Terms in Bayesian Analysis
#'
#' This function checks the specified model type for fixed terms in a Bayesian analysis context,
#' ensuring it is correctly set to 'FIXED'. If not set or incorrectly specified, it automatically corrects it.
#'
#' @param fixed_term Unused parameter placeholder for potential future use, allowing for flexibility in function calls.
#' @param fixed_term_model_bayesian Character string specifying the model type for fixed terms. Expected to be 'FIXED'.
#' @param message Logical indicating whether to display warning messages when the model type is missing or incorrectly specified.
#'
#' @return A character string 'FIXED', confirming the model type for fixed terms is correctly set.
#'
#' @details
#' In the context of Bayesian analysis, especially within certain statistical or machine learning frameworks,
#' fixed terms may require specific model type settings. This function ensures that the model type for fixed terms
#' is explicitly set to 'FIXED', correcting it if necessary and optionally warning the user of any adjustments made.
#' This process helps maintain consistency and prevent errors in model specification.
#'
#' @examples
#' fixed_term_model_bayesian <- fixed_term_model(fixed_term_model_bayesian = "INCORRECT", message = TRUE)
#' print(fixed_term_model_bayesian)
#'
#' @export

### Check if user provide model for fixed term. Which is typically FIXED
fixed_term_model <- function(fixed_term = NULL,
                             fixed_term_model_bayesian = NULL,
                             message = TRUE){

  msg <- "\n ==================================================\n"
  if(is.null(fixed_term_model_bayesian)){
    if(isTRUE(message)){
    warning(msg, "The model for fixed term(s) is missing. We fix it for you.")
}
    fixed_term_model_bayesian = 'FIXED'
  } else{

    if(fixed_term_model_bayesian!='FIXED'){
      if(isTRUE(message)){
      warning(msg, "The model for fixed term(s) should be equal to FIXED. We fix it for you.")
}
    }

    if(length(fixed_term_model_bayesian)>1){

      fixed_term_model_bayesian <- 'FIXED'
    }
  }

  return(fixed_term_model_bayesian)
}


#' Validate and Adjust Model Specification for Random Terms
#'
#' This function evaluates and adjusts the model specifications for random terms in a Bayesian analysis setting,
#' ensuring they align with the specified general and random term models. It handles various scenarios including
#' interactions, the presence of genetic names, and checks for consistency with available Bayesian models.
#'
#' @param rand_terms A character vector specifying the random effect terms to be included in the model.
#' @param GS_model A character string or vector specifying the general statistical model(s) to be applied,
#' such as "BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS".
#' @param gen_name A character string specifying the name of the genetic effect in the model.
#' @param rand_terms_model_bayesian A character vector specifying the model(s) for each random term.
#' If not specified or inconsistent, default adjustments are made based on `GS_model` and available terms.
#' @param message A logical indicating whether to display warning messages for missing or incorrect model specifications.
#'
#' @return A character vector with the adjusted or validated model specifications for each random term.
#'
#' @details
#' The function ensures that the model specifications for random terms are correctly set for Bayesian analysis.
#' It verifies the presence and appropriateness of the genetic effect name, handles interaction terms,
#' and adjusts model specifications based on the general model (`GS_model`) provided. If `GS_model` is specified
#' but `rand_terms_model_bayesian` is missing or incomplete, the function attempts to intelligently assign models
#' to each term, defaulting to "BRR" where necessary and using "RKHS" for multi-kernel interactions.
#'
#' @examples
#' rand_terms <- c("Genotype", "Genotype:Environment")
#' GS_model <- "BRR"
#' gen_name <- "Genotype"
#' rand_terms_model_bayesian <- c("BRR", "RKHS")
#' adjusted_models <- random_term_model(rand_terms = rand_terms,
#'                                      GS_model = GS_model,
#'                                      gen_name = gen_name,
#'                                      rand_terms_model_bayesian = rand_terms_model_bayesian,
#'                                      message = TRUE)
#' print(adjusted_models)
#'
#' @export
random_term_model <- function(rand_terms = NULL,
                              GS_model = NULL,
                              gen_name = NULL,
                              rand_terms_model_bayesian = NULL,
                              message = TRUE){


  msg <- "\n ==================================================\n"

  #####
  #### Check for Interaction and and non-interaction term
  ## No interaction term
  rand_terms_no_inter <-  rand_terms[!grepl(":", rand_terms)]

  ## Interaction term
  check_rand_inter <-  rand_terms[grepl(":", rand_terms)]

  non_gen_pos_mod <- NULL

  ## Start with the No interaction terms
  if(length(rand_terms_no_inter)!=0){
    ## Check if gen_name is present and store the position
    gen_pos_mod <-   match(gen_name, rand_terms)
    if(length(gen_pos_mod)==0){stop(print(paste(message(msg), paste(gen_name, "effect is missing"))), call. = FALSE)}
    if(length(gen_pos_mod)>1){stop(print(paste(message(msg), paste(gen_name, "effect should not be greater than 1"))), call. = FALSE)}
    ## Extract other terms from the rand_terms_no_inter  expect the gen_name
    rand_terms_no_inter_no_gen <-  rand_terms_no_inter[!rand_terms_no_inter%in%gen_name]
    ## Get the position of other terms (No interaction) in the random that is not gen_name
    non_gen_pos_mod <- match(rand_terms_no_inter_no_gen, rand_terms)

  }

  ######## Initialize step For model adjustment for the random term with interaction

  if(length(check_rand_inter)!=0){

    test_present_of_geno <-  grep(gen_name, check_rand_inter, value = TRUE)

    if(length(test_present_of_geno)!=0 | !is.na(test_present_of_geno)){

      inter_gen_pos_mod <-  match(test_present_of_geno, rand_terms)

      non_gen_inter_test <-   check_rand_inter[!check_rand_inter%in%test_present_of_geno]

    } else {

      non_gen_inter_test <-  NULL

      inter_gen_pos_mod <-  NULL
    }

    if(!is.null(non_gen_inter_test)){

      inter_non_gen_pos_mod <-  match(non_gen_inter_test, rand_terms)

    }

  } ## End

  ### Each random term will have a specific model for parameter estimate.
  ## It is expected the user will provide model for each random term,
  ## In a scenario where the model provided is not equal to the number of random terms
  ## or no model was provided at all. This next line of code take care of it.

  ### If user provide only the GS_model and random terms model is missing
  ### Especially when there is more than one random term

  valid_models <- c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")
  if(is.null(rand_terms_model_bayesian) & !is.null(GS_model)){
    #warning(paste(msg, "The model for ranom term(s) is missing. We fix it for you"))

    ## Check to see which GS_model(s) the user provide
    #if(!is.null(GS_model)){
    ### Here the GS_model is of length rand_terms
    if(length(GS_model)==length(rand_terms)){
      ### Check GS_model to be consistent with models present in the engine
      mod_present_in_GS_model = GS_model[GS_model%in% valid_models]
      if(length(mod_present_in_GS_model)==length(rand_terms)){
        rand_terms_model_bayesian = GS_model
        #rand_mod_copy = GS_model
      }
      ### if GS_model is less than the length of rand_terms
    } else{
      mod_present_in_GS_model = GS_model[GS_model%in% valid_models]
      if(length(mod_present_in_GS_model)!=0){
        #rand_mod_copy = mod_present_in_GS_model
        if(message){
          warning(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed.")
        }
        mod_len = length(rand_terms) - length(mod_present_in_GS_model)

        rand_terms_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
        rand_terms_model_bayesian[gen_pos_mod] <- GS_model[1]

        ## For other terms in random effect aside gen_name
        if(!is.null(non_gen_pos_mod)){

          rand_terms_model_bayesian[non_gen_pos_mod] <-  "BRR"

        } ## end

        #### Random interaction terms
        if(length(check_rand_inter)!=0){
          ### For gen_name part of interaction
          if(!is.na(inter_gen_pos_mod) | length(inter_gen_pos_mod)!=0){

            ### This account for multi-kernel
            rand_terms_model_bayesian[inter_gen_pos_mod] <- "RKHS"


          }

          ### For Non gen_name part of interaction
          if(length(non_gen_inter_test)!=0){

            rand_terms_model_bayesian[non_gen_inter_test] <- "BRR"

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
        warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed."))
      }
      ### Check GS_model to be consistent with models present in the engine
      mod_present_in_model_bayesian = rand_terms_model_bayesian[rand_terms_model_bayesian%in% valid_models]
      if(length(mod_present_in_model_bayesian)==0){

        stop(msg,'Provided appropiate name for the Baysian model.\n ')

      } else {

        mod_len <-  (length(rand_terms) - length(mod_present_in_model_bayesian))
        rand_terms_model_bayesian <-  c(mod_present_in_GS_model,  rep("BRR",  mod_len))
        #rand_mod_copy = rand_terms_model_bayesian

      }

      rand_terms_model_bayesian[gen_pos_mod] <- GS_model[1]

      ## For other terms in random effect aside gen_name
      if((!is.null(non_gen_pos_mod) & length(non_gen_pos_mod)==0) & !is.na(non_gen_pos_mod)){

        rand_terms_model_bayesian[non_gen_pos_mod] <-  "BRR"

      } ## end

      #### Random interaction terms
      if(length(check_rand_inter)!=0){
        ### For gen_name part of interaction
        if(!is.na(inter_gen_pos_mod) | length(inter_gen_pos_mod)!=0){

          ### This account for multi-kernel
          rand_terms_model_bayesian[inter_gen_pos_mod] <-  "RKHS"


        }

        ### For Non gen_name part of interaction
        if(length(non_gen_inter_test)!=0){

          rand_terms_model_bayesian[non_gen_inter_test] <-  "BRR"

        }

      } ## End of random interaction terms

    }

  }else {

    ### When user provide provide both GS_model and rand_terms_model_bayesian
    if(!is.null(rand_terms_model_bayesian) & !is.null(GS_model)){
 ### Check which GS_model the user supply, it has it to match the available models
      mod_present_in_GS_model <- GS_model[GS_model%in% valid_models]
      if(length(mod_present_in_GS_model)==0){
        stop(print(paste(msg,'Provided appropiate name for the baysian model in GS_model.')), call. = FALSE)
      }
      ### Check which rand_terms_model_bayesian the user supply, it has it to match the available models
      mod_present_in_model_bayesian <- rand_terms_model_bayesian[rand_terms_model_bayesian%in% valid_models]
      if(length(mod_present_in_model_bayesian)==0){
        stop(print(paste(msg,'Provided appropiate name for the baysian model in rand_terms_model_bayesian.')), call. = FALSE)
      }
      ### Check if the models term define is greater than the random term
      if(length(rand_terms_model_bayesian)>length(rand_terms)){
        stop(print(paste(msg,'The number of model is greater than the random terms.')), call. = FALSE)
      }
      ##
      if(length(mod_present_in_GS_model)!=0 & length(mod_present_in_model_bayesian)!=0){

        ## Combined all the models
        rand_terms_model_bayesian  <- c(GS_model, rand_terms_model_bayesian)

        if(length(rand_terms_model_bayesian)!= length(rand_terms)){

          if(length(rand_terms_model_bayesian)>length(rand_terms)){
            #stop(print(paste(msg,'The number of model is greater than the random terms')), call. = FALSE)
            rand_terms_model_bayesian <-  rand_terms_model_bayesian[length(rand_terms)]
            }

          if(length(rand_terms_model_bayesian)<length(rand_terms)){
          mod_len <-  (length(rand_terms) - length(rand_terms_model_bayesian))
          rand_terms_model_bayesian <-  c(rand_terms_model_bayesian,  rep("BRR",  mod_len))

          }

        }

        rand_terms_model_bayesian[gen_pos_mod] <-  GS_model[1]
      }
      ###
      ###  rand_terms_model_bayesian[gen_pos_mod] = GS_model[1]
      ## For other terms in random effect aside gen_name
      if((!is.null(non_gen_pos_mod) & length(non_gen_pos_mod)==0) & !is.na(non_gen_pos_mod)){

        rand_terms_model_bayesian[non_gen_pos_mod] <-  "BRR"

      } ## end

      #### Random interaction terms
      if(length(check_rand_inter)!=0){
        ### For gen_name part of interaction
        if(!is.na(inter_gen_pos_mod) | length(inter_gen_pos_mod)!=0){

          ### This account for multi-kernel
          rand_terms_model_bayesian[inter_gen_pos_mod] <-  "RKHS"


        }

        ### For Non gen_name part of interaction
        if(length(non_gen_inter_test)!=0){

          rand_terms_model_bayesian[non_gen_inter_test] <-  "BRR"

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
#   msg <- "\n==================================================\n"
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
