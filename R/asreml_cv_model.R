


#' Title
#'
#' @param pheno_dataa
#' @param y
#' @param code_asr_fit_cv
#' @param tst
#'
#' @return
#' @export
#'
#' @examples
asreml_cv_model <- function(pheno_dataa = NULL,
                            response = NULL,
                            gen_name = NULL,
                            #heter_groups = NULL,
                            asreml_models_prep_cv = NULL,
                            tst = NULL
                             ){

  asreml::asreml.options(trace=FALSE)
  names_in_inv_list <-  asreml_models_prep_cv[["names_in_inv_list"]]
  code_asr_fit_cv <-  asreml_models_prep_cv[["code_asr_fit"]]
  pheno_dataa[tst, response] <- NA
  #gen_tst <- pheno_dataa[tst, gen_name]
  code_asr_fit_cv[4] <- 'na.action=list(x="include",y="include"),data=pheno_dataa)'
  inv_list <- asreml_models_prep_cv[["inv_list"]]
  ####
  #code.asr[1] <- paste('mod<-', code.asr[1], sep='')
  code_asr_fit_cv[1] <- paste('mod_cv<-', code_asr_fit_cv[1], sep='')
  str_mod_cv <- paste(code_asr_fit_cv[1],code_asr_fit_cv[2],code_asr_fit_cv[3],code_asr_fit_cv[4],sep=',')

  ## This is useful because asreml want the inv_object in the current environment
  ## though present in the global environment. so this call it from the global
  ## to the current environment
  # Loop through each name in the list
  #cat(names_in_inv_list)
  # if(is.null(names_in_inv_list)| length(names_in_inv_list)==0){
  #   message("it is null")
  # }
  #
  # # Assign the object to the current environment with the same name
  #
  # for (name in names_in_inv_list) {
  #   message(names_in_inv_list)
  #   # Check if the object exists in the global environment
  #   #if (exists(name, envir = .GlobalEnv)) {
  #   if (name %in% ls(envir = .GlobalEnv)) {
  #     # Get the object from the global environment
  #     obj <- get(name, envir = .GlobalEnv)
  #
  #     # Assign the object to the current environment with the same name
  #     assign(name, obj, envir = environment())
  #     # Print a confirmation message
  #     message(paste("Object", name, "has been copied to the current environment."))
  #
  #   } else {
  #     # Print a message if the object does not exist in the global environment
  #     message(paste("Object", name, "not found in the global environment."))
  #   }
  # }

  ##### THis is important for asreml inorder to update the model if need be
  for (i in seq_along(inv_list)) {
    assign(names(inv_list)[i], inv_list[[i]], envir = .GlobalEnv)
  }
  # Loop through each variable name in the list

  ## Calls the current environment for evaluation
  eval(parse(text=str_mod_cv), envir=environment())

  for (i in 1:3) {
    if (!mod_cv$converge) {

      eval(parse(text='mod_cv<-asreml::update.asreml(mod_cv)'))
    }else {
      # Exit the loop if model has converged
      break
    }
  }

  # # Attempt to evaluate in the current environment
  # tryCatch({
  #   eval(parse(text = str_mod_cv), envir = environment())
  # }, error = function(e) {
  #   # If an error occurs, check for the missing object in the global environment
  #   message("Error occurred: ", e$message)
  #   message("Attempting to execute in the global environment...")
  #   eval(parse(text = str_mod_cv), envir = .GlobalEnv)
  # })



  output <- list(model_cv = mod_cv,
                 asreml_models_prep_cv = asreml_models_prep_cv
                 #gen_tst = gen_tst
                 )


  return(output)

}
