
#' Cross-Validation for ASReml Models
#'
#' This function performs cross-validation for ASReml models using the provided phenotypic data.
#'
#' @param pheno_dataa A data frame containing the phenotypic data.
#' @param response A string specifying the response variable in the phenotypic data.
#' @param gen_name A string specifying the column name for the genotypic data.
#' @param heter_groups A factor or vector specifying the heterogeneity groups.
#' @param asreml_models_prep_cv A list containing the prepared ASReml model components for cross-validation.
#' @param tst An integer vector specifying the indices of the test set in the phenotypic data.
#'
#' @return A list containing the fitted cross-validation model, the prepared ASReml model components, and the modified phenotypic data.
#' \itemize{
#'   \item{model_cv}{The fitted ASReml cross-validation model.}
#'   \item{asreml_models_prep_cv}{The prepared ASReml model components.}
#'   \item{pheno_dataa}{The modified phenotypic data with missing values for the test set.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Toy example phenotypic data
#'   pheno_data <- data.frame(
#'     GID = 1:10,
#'     Trait = rnorm(10),
#'     Group = factor(rep(1:2, each = 5))
#'   )
#'
#'   # Mock preparation of ASReml model components for cross-validation
#'   asreml_models_prep_cv <- list(
#'     names_in_inv_list = c("GID"),
#'     code_asr_fit = c("Trait ~ 1 + Group", "random = ~ GID", "residual = ~ units", ""),
#'     inv_list = list(GID = diag(10))
#'   )
#'
#'   # Indices of the test set
#'   tst_indices <- c(1, 3, 5, 7, 9)
#'
#'   # Perform cross-validation for ASReml models
#'   result <- asreml_cv_model(
#'     pheno_dataa = pheno_data,
#'     response = "Trait",
#'     gen_name = "GID",
#'     heter_groups = pheno_data$Group,
#'     asreml_models_prep_cv = asreml_models_prep_cv,
#'     tst = tst_indices
#'   )
#'   print(result)
#' }
#'
#' @export

asreml_cv_model <- function(pheno_dataa = NULL,
                            response = NULL,
                            gen_name = NULL,
                            heter_groups = NULL,
                            asreml_models_prep_cv = NULL,
                            tst = NULL
                             ){

  GIDs <- as.character(pheno_dataa[[gen_name]])
  asreml::asreml.options(trace=FALSE)
  names_in_inv_list <-  asreml_models_prep_cv[["names_in_inv_list"]]
  code_asr_fit_cv <-  asreml_models_prep_cv[["code_asr_fit"]]
  #pheno_dataa[as.character(pheno_dataa[[gen_name]])%in%GIDs[tst], response] <- NA
  pheno_dataa[tst, response] <- NA
  #pheno_dataa[tst, response] <- NA
  #gen_tst <- pheno_dataa[tst, gen_name]
  code_asr_fit_cv[1] <-  gsub("trait", response, code_asr_fit_cv[1])
  code_asr_fit_cv[4] <- 'na.action=list(x="include",y="include"),data=pheno_dataa)'
  inv_list <- asreml_models_prep_cv[["inv_list"]]
  ####
  #code.asr[1] <- paste('mod<-', code.asr[1], sep='')
  code_asr_fit_cv[1] <- paste('mod_cv<-', code_asr_fit_cv[1], sep='')
  str_mod_cv <- paste(code_asr_fit_cv[1],code_asr_fit_cv[2],code_asr_fit_cv[3],code_asr_fit_cv[4],sep=',')

  ##### THis is important for asreml inorder to update the model if need be
  for (i in seq_along(inv_list)) {
    assign(names(inv_list)[i], inv_list[[i]], envir = .GlobalEnv)
  }
  # Loop through each variable name in the list

  result_model <-  tryCatch(
    {
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
      mod_cv
    },
  error = function(e) {
    # Handle the error, you can print a message or take other actions
    cat(paste("Asreml model fail:", conditionMessage(e), "\n"))
    return(NULL)  # Return NULL or an appropriate value to indicate the failure
  }
  )

if(!is.null(result_model)){

  output <- list(model_cv = result_model,
                 asreml_models_prep_cv = asreml_models_prep_cv,
                 pheno_dataa = pheno_dataa
                 #gen_tst = gen_tst
                 )

}else{
  output <- list(model_cv = NULL,
                 asreml_models_prep_cv = NULL,
                 pheno_dataa = NULL
                 #gen_tst = gen_tst
  )
}


  return(output)

}
