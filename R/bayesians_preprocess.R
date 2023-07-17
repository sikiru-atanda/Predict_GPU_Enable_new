

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
# if not assign a default value, but let the user be aware if it
#' Title
#'
#' @param nIter
#' @param burnIn
#' @param thin
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
    if(message){
    message(paste(msg, "Number of iteration is missing. Default value of 16000 was assigned. \n Check if this appropriate for your data"))
    }
      nIter = 200

  }

  #if(is.null(burnIn) || burnIn < 5000){
  if(is.null(burnIn)){

    if(message){
      message(paste(msg, "Number of burn-in is missing. Default value of 6000 was assigned. \n Check if this appropriate for your data"))
    }
      burnIn = 30

  }

  #if(is.null(thin) || thin<10){
  if(is.null(thin)){
    if(message){
    message(paste(msg, "Number of thining is missing. Default value of 10 was assigned. \n Check if this appropriate for your data"))
    }
      thin = 5

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

### Adjust for fixed terms.
#' Title
#'
#' @param fixed
#' @param object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fixed_terms <- function(
    fixed  = NULL,
    object = NULL,
    ...){

  fixed <- rand_fix_check(rand_fix_term = fixed,
                           object= object)

  if(attr(fixed, "cleared")!="pass" & class(fixed)!="forumla") {

    msg <- sprintf("==================================================\n")
    stop(print(paste(msg,paste(object, 'is not class formula.', sep = ""))), call. = FALSE)
  }
  # msg <- sprintf("==================================================\n")
  # ### Check if the random term provided are present in the phenotypic data
  #
  # if(!all(sapply(all.vars(fixed), function(x, object) x%in%names(object),  object))) {
  #   stop(print(paste(msg,"All variables indicated in argument 'random' should be present in phenotypic data")), call. = FALSE)
  # }

   fixed_term <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
   ## Remove any space in the terms
   fixed_term = gsub(" ", "", fixed_term)

   return(fixed_term)

}


## Adjust for the interaction terms. Here single environment was in mind.
## Multiple environment will be build ontop of it
#' Title
#'
#' @param random
#' @param object
#' @param ..
#'
#' @return
#' @export
#'
#' @examples
random_terms <- function(random = NULL,
                         object = NULL,
                         ...){

  msg <- sprintf("==================================================\n")

  random <- rand_fix_check(rand_fix_term = random,
                           object= object)

  if(attr(random, "cleared")!="pass" && class(random)!="forumla") {

    stop(print(paste(msg,paste('the random term', 'is not formular.', sep = ""))), call. = FALSE)
  }

  # msg <- sprintf("==================================================\n")
  # ### Check if the random term provided are present in the phenotypic data
  # if(!all(sapply(all.vars(random), function(x, object) x%in%names(object),  object))) {
  #   stop(print(paste(msg,"All variables indicated in argument 'random' should be present in phenotypic data")), call. = FALSE)
  # }
  rand_term <- strsplit(as.character(random[2]), split = "[+]")[[1]]
  ## Remove any space in the terms
  rand_term = gsub(" ", "", rand_term)

  return(rand_term)
}


#' Title
#'
#' @param fixed_term
#' @param fixed_term_model_bayesian
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
    if(message){
    warning(paste(msg, "The model for fixed term(s) is missing. We fix it for you"))
}
    fixed_term_model_bayesian = 'FIXED'
  } else{

    if(fixed_term_model_bayesian!='FIXED'){
      if(message){
      warning(paste(msg, "The model for fixed term(s) should be equal to fix. We fix it for you"))
}
    }

    if(length(fixed_term_model_bayesian)>1){

      fixed_term_model_bayesian = 'FIXED'
    }
  }

  return(fixed_term_model_bayesian)
}


#' Title
#'
#' @param rand_terms
#' @param GS_model
#' @param rand_term_model_bayesian
#'
#' @return
#' @export
#'
#' @examples

## Check the random term provided by the users and provide one in scenario where
## user provide one or less than the random terms.
# Make adjustment if more than the random terms. User should be let aware of the
# implications
random_term_model <- function(rand_terms = NULL,
                              GS_model = NULL,
                              rand_term_model_bayesian = NULL,
                              message = TRUE)
  {

  msg <- sprintf("==================================================\n")
  ### Each random term will have a specific model for parameter estimate.
  ## It is expected the user will provide model for each random term,
  ## In a scenario where the model provided is not equal to the number of random terms
  ## or no model was provided at all. This next line of code take care of it.

  ### If user provide only the GS_model and random terms model is missing
  ### Especially when there is more than one random term
  if(is.null(rand_term_model_bayesian) && !is.null(GS_model)){
    #warning(paste(msg, "The model for ranom term(s) is missing. We fix it for you"))

    ## Check to see which GS_model(s) the user provide
    #if(!is.null(GS_model)){
    ### Here the GS_model is of length rand_terms
    if(length(GS_model)==length(rand_terms)){
      ### Check GS_model to be consistent with models present in the engine
      mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_GS_model)==length(rand_terms)){
        rand_term_model_bayesian = GS_model
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

        rand_term_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
      } else {
        stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)
      }
    }

    #}

  } else {
    ### If user provide only rand_term_model_bayesian and GS_model is missing
    ### Especially when there is more than one random term
    if(!is.null(rand_term_model_bayesian) && is.null(GS_model)){
      # if(length(rand_term_model_bayesian)!= length(rand_terms)){
      #   warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
      # }
      ### This will be used only for the genetic effect assuming the user only provide that
      if(length(rand_term_model_bayesian)!= length(rand_terms)){
        if(message){
        warning(paste(msg, "The model for random term(s) should be equal to the total number of random term.\n\t Default model (BRR) was assigned, provide desired models if needed"))
        }
          ### Check GS_model to be consistent with models present in the engine
        mod_present_in_model_bayesian = rand_term_model_bayesian[rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
        if(length(mod_present_in_model_bayesian)==0){ stop(print(paste(msg,'Provided appropiate name for the Baysian model.')), call. = FALSE)

        } else {

          mod_len = length(rand_terms) - length(mod_present_in_model_bayesian)
          rand_term_model_bayesian = c(mod_present_in_GS_model,  rep("BRR",  mod_len))
          #rand_mod_copy = rand_term_model_bayesian

        }

      }

    }

    ### If user provide GS_model and rand_term_model_bayesian
    if(!is.null(rand_term_model_bayesian) && !is.null(GS_model)){
      mod_present_in_GS_model = GS_model[GS_model%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_GS_model)==0){stop(print(paste(msg,'Provided appropiate name for the baysian model in GS_model.')), call. = FALSE)}

      mod_present_in_model_bayesian = rand_term_model_bayesian[rand_term_model_bayesian%in%c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")]
      if(length(mod_present_in_model_bayesian)==0){ stop(print(paste(msg,'Provided appropiate name for the baysian model in rand_term_model_bayesian.')), call. = FALSE)}

      if(length(mod_present_in_GS_model)!=0 && length(mod_present_in_model_bayesian)!=0){

        rand_term_model_bayesian  = c(GS_model, rand_term_model_bayesian)

        if(length(rand_term_model_bayesian)!= length(rand_terms)){

          if(length(rand_term_model_bayesian)>length(rand_terms)){
            stop(print(paste(msg,'The number of model is greater than the random terms')), call. = FALSE)
          }

          mod_len = length(rand_terms) - length(rand_term_model_bayesian)
          rand_term_model_bayesian = c(rand_term_model_bayesian,  rep("BRR",  mod_len))

        }
      }


    }

    ### if user didn't specific any GS_model the default rrBLUP will be used
    if(is.null(rand_term_model_bayesian) && is.null(GS_model)){
      rand_term_model_bayesian = array(dim = length(rand_terms), "BRR")
    }


  }

  ###
  # For Didier processing the first element of the vector will be used
  # as pointer shiny production because that is the on that represent the model
  # for genetic effect.

  return(rand_term_model_bayesian)

}

