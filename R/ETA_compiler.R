#' Compile ETA for Fixed Terms in Bayesian Models
#'
#' This function prepares the ETA components for fixed terms in Bayesian genomic prediction models.
#' It uses the specified fixed effects from the phenotypic data and compiles them into a format suitable
#' for inclusion in a Bayesian analysis. This is particularly useful for preparing data for models
#' such as GBLUP, RKHS, and BRR.
#'
#' @param fixed A formula specifying the fixed effects to be included in the model.
#' @param fixed_term_model_bayesian Optionally specify the model for the fixed terms in the Bayesian framework.
#'        This is usually "FIXED" but can be left NULL for automatic handling.
#' @return A list of lists where each inner list represents an ETA component for a fixed term.
#'         Each component contains the model matrix (`X`) and the model type (`model`) for that fixed term.
#' @examples
#' \dontrun{
#' # Assuming pheno_data is your phenotypic dataset and fixed_effect is your fixed effect formula:
#' result <- ETA_compiler_fixed_term(fixed = ~fixed_effect1 + fixed_effect2,
#'                                   pheno_data = pheno_data)
#' }
#' @export

gp_bayes_normalize_fixed_argument <- function(fixed = NULL,
                                              fixed_term_model_bayesian = NULL) {
  if (is.null(fixed) || inherits(fixed, "formula")) {
    return(fixed)
  }
  if (length(fixed) == 1L &&
      length(fixed_term_model_bayesian) == 1L &&
      identical(as.character(fixed), as.character(fixed_term_model_bayesian))) {
    return(NULL)
  }
  fixed
}

#' Compile ETA for Fixed Terms in Bayesian Models
#'
#' This function prepares model-matrix ETA components for fixed terms in
#' Bayesian genomic prediction models.
#'
#' @param fixed A formula specifying fixed effects.
#' @param fixed_term_model_bayesian Optional model label for fixed terms.
#' @param pheno_data Phenotypic data containing variables in `fixed`.
#' @return A list of ETA components for fixed terms.
#' @export
ETA_compiler_fixed_term <- function(fixed = NULL,
                         fixed_term_model_bayesian = NULL,
                         pheno_data = NULL
                         ){

  ### Initialize steps for compiling the Fixed terms
  ## Start with creating empty list for ETA compilation

  ETA = list()
  fixed <- gp_bayes_normalize_fixed_argument(
    fixed = fixed,
    fixed_term_model_bayesian = fixed_term_model_bayesian
  )

  if (is.null(fixed)) {
    return(ETA)
  }

    fixed_term_no_inter <- fixed_terms(fixed = fixed, pheno_data = pheno_data)
    fixed_term_no_inter <- setdiff(fixed_term_no_inter, c("1", "0"))
    fixed_model <- fixed_term_model(fixed_term_no_inter,
                                    fixed_term_model_bayesian)
    #### Fit Fixed terms in ETA
    if(length(fixed_term_no_inter)!=0){
      for (ET in 1:length(fixed_term_no_inter)) {


        fixed_value <- pheno_data[[fixed_term_no_inter[ET]]]
        fixed_design <- stats::model.matrix(
          ~ .fixed_value,
          data = data.frame(.fixed_value = fixed_value)
        )
        fixed_design <- fixed_design[
          , colnames(fixed_design) != "(Intercept)", drop = FALSE
        ]
        informative <- if (ncol(fixed_design)) {
          vapply(seq_len(ncol(fixed_design)), function(j) {
            x <- fixed_design[, j]
            x <- as.numeric(x)
            any(is.finite(x)) && stats::sd(x, na.rm = TRUE) > sqrt(.Machine$double.eps)
          }, logical(1L))
        } else {
          logical()
        }
        fixed_design <- fixed_design[, informative, drop = FALSE]
        if (!ncol(fixed_design)) {
          next
        }
        ETA[[length(ETA) + 1L]] <- list(
          X = fixed_design,
          model = fixed_model[[min(ET, length(fixed_model))]]
        )

      }

    }

  return(ETA)

}


#' Compile ETA for random and fixed term combined for Bayesian Genomic Prediction Models
#'
#' This function compiles the effects to be estimated (ETA) for Bayesian genomic prediction models, including fixed and random effects, using genomic and omics data. It prepares the model matrices for fixed effects and the relationship matrices for random effects to be used in Bayesian analysis, such as GBLUP, RKHS, and Bayesian Regression.
#'
#' @param fixed Formula specifying the fixed effects to be included in the model.
#' @param random Formula specifying the random effects to be included in the model.
#' @param GS_model A character string specifying the genomic selection model to be used (e.g., "BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS").
#' @param fixed_term_model_bayesian Optionally, specify the model for the fixed terms in the Bayesian framework.
#' @param rand_term_model_bayesian Optionally, specify the model for the random terms in the Bayesian framework.
#' @param pheno_data A data frame containing the phenotypic data, including columns for all variables specified in the `fixed` and `random` formulas.
#' @param geno_data Genomic relationship matrix or genotype data.
#' @param omic1_data First omics data matrix.
#' @param omic2_data Second omics data matrix.
#' @param omic3_data Third omics data matrix.
#' @param gen_name Name of the genotype column in `pheno_data`.
#' @param scaling Logical; if `TRUE`, centre and scale the marker / omics data
#'   before compiling the random-effects ETA components.
#' @param ... Reserved for future extensions; currently ignored.
#' @return A list containing the compiled ETA components, the modified phenotypic data, and names of the ETA elements. Each element in the ETA list represents a model component (fixed or random effect) along with its corresponding model matrix (`X`) or relationship matrix (`K`), and the specified model type.
#' @examples
#' \dontrun{
#' # Assuming pheno_data is your phenotypic dataset and you have genomic and omics data:
#' ETA_result <- ETA_compiler_bayes(
#'   fixed = ~ fixed_effect,
#'   random = ~ (1|gen_name),
#'   GS_model = "BRR",
#'   pheno_data = pheno_data,
#'   geno_data = geno_matrix,
#'   omic1_data = omic1_matrix,
#'   # Additional omics data can be added as needed
#'   gen_name = "GenotypeID"
#' )
#' }
#' @export
#'
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
    scaling = TRUE,
    ...
) {
  ETA <- list()
  msg <- ""
  fixed <- gp_bayes_normalize_fixed_argument(
    fixed = fixed,
    fixed_term_model_bayesian = fixed_term_model_bayesian
  )

  rand_term_no_inter <- random_terms(random = random, pheno_data = pheno_data)

  rand_model <- random_term_model(
    rand_terms = rand_term_no_inter,
    gen_name = gen_name,
    rand_terms_model_bayesian = rand_term_model_bayesian,
    GS_model = GS_model,
    message = TRUE
  )

  if (!is.null(fixed)) {
    ETA <- ETA_compiler_fixed_term(
      fixed = fixed,
      pheno_data = pheno_data,
      fixed_term_model_bayesian = fixed_term_model_bayesian
    )
  }

  if (length(rand_model) != 1 || length(rand_term_no_inter) != 1) {
    stop(paste(msg, 'This works only for one random effect'))
  }

  datasets <- list(geno_data, omic1_data, omic2_data, omic3_data)
  dataset_names <- c("geno_data", "omic1_data", "omic2_data", "omic3_data")
  datasets_index <- which(!sapply(datasets, is.null))
  datasets <-  datasets[datasets_index]
  dataset_names <- dataset_names[datasets_index]
  ETA_element_name <- character()

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    if("BRR" %in% rand_model){
      if(!is.null(dataset)) {
        if(isTRUE(scaling)){
        #dataset <- scale(dataset, center = TRUE, scale = TRUE)
          dataset <- scale(dataset, center = TRUE, scale = apply(dataset, 2, function(x) max(sd(x), 1e-6)))
          cols_with_na <- which(colSums(is.na(dataset)) > 0)
          if(length(cols_with_na)!=0){
          dataset <- dataset[, -cols_with_na]

          }
          #colnames(dataset)[cols_with_na]
        }
      }
    }
    if (!is.null(dataset)) {  ## this seems redundant but necessary
      ETA[[length(ETA) + 1]] <- list(X = as.matrix(dataset), model = rand_model, saveEffects = TRUE)
      ETA_element_name <- c(ETA_element_name, dataset_names[i])
    }
  }

  output <- list(ETA = ETA, pheno_data = pheno_data, ETA_element_name = ETA_element_name)
  names(output) <- c("ETA", "pheno_data", "ETA_element_name")
  return(output)
}
