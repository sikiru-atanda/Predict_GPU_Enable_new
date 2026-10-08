
#' Validate and Process Random/Fixed Effect Terms Against Phenotypic Data
#'
#' This function checks the validity of specified random or fixed effect terms against a provided phenotypic dataset. It ensures that the terms are correctly specified as a formula and that all variables in the formula are present in the dataset. It also attempts to convert variables to factors if they are not already.
#'
#' @param pheno_data A data frame, data.table, or an object of class "phenotype" containing the phenotypic data against which the terms are validated.
#' @param rand_fix_term A formula specifying the random or fixed effect terms to be included in the model.
#' @param term_type character indicating if the term is random or fixed
#' @param ... Additional arguments for future extensions or compatibility.
#'
#' @return The validated and potentially modified `rand_fix_term` as a formula object, with an attribute "cleared" set to "pass" indicating successful validation.
#'
#' @details
#' The function is crucial for ensuring the integrity of the model specification before executing statistical analyses, particularly in genetic studies. It performs several checks:
#' - Verifies that `rand_fix_term` is a formula.
#' - Ensures that all variables specified in `rand_fix_term` are present in `pheno_data`.
#' - Converts variables in `pheno_data` referenced by `rand_fix_term` to factors, if they are not already.
#' This preprocessing helps prevent common errors in model specification and data analysis workflows.
#'
#' @examples
#' pheno_data <- data.frame(Trait1 = rnorm(100), Trait2 = rnorm(100), Genotype = factor(rep(1:10, each = 10)))
#' rand_fix_term <- ~ Trait1 + Trait2 + Genotype
#' validated_term <- rand_fix_check(pheno_data = pheno_data, rand_fix_term = rand_fix_term)
#' print(validated_term)
#'
#' @export

rand_fix_check <- function(pheno_data = NULL,
                           rand_fix_term = NULL,
                           term_type = NULL,
                            ...){


  msg <- ""

  if (!inherits(rand_fix_term, 'formula')) {
    stop(msg, "The ", term_type, " term is not a class of type 'formula'. Example: ", term_type, " = ~ X + Y")
  }

  term_vars <- all.vars(rand_fix_term)
  if (!all(term_vars %in% names(pheno_data))) {
    stop(msg, "All variables indicated in the ", term_type, " term should be present in phenotypic data.")
  }

    ## Check all variables in pheno_data term are factors

  non_factor_vars <- term_vars[!sapply(pheno_data[term_vars], is.factor)]
  if (length(non_factor_vars) > 0) {
    pheno_data[non_factor_vars] <- lapply(pheno_data[non_factor_vars], factor)
  }

  attr(rand_fix_term, "cleared") <- "pass"
  return(rand_fix_term)


}
