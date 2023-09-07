
#' Title
#' This function do the following:
#'
#' 1. Check if the random and fixed term are in formula. If not it will stop
#' 2. Check if the ranom and fixed term are in class factor in the pheno_data,
#'    if not we fix it for the user
#' 3. The output will be class formula and attribute pass
#'
#' @param pheno_data phenotypic data
#' @param ...
#' @param rand_fix_term this can be random or fixed term and expected to be class formula
#'
#' @return
#' @export
#'
#' @examples
rand_fix_check <- function(pheno_data = NULL,
                      rand_fix_term = NULL,
                      ...){
  msg <- sprintf("==================================================\n")

  ## pheno_data is the phenotypic data

  # if(attr(pheno_data, "cleared")!="pass" & class(pheno_data)!=c("data.table", "data.frame", "phenotype")) {
  #
  #   stop(print(paste(msg,paste(pheno_data, 'is not class phenotype.', sep = ""))), call. = FALSE)
  # }
    if (inherits(rand_fix_term, what = 'formula')) {
      rand_fix_term <- stats::as.formula(rand_fix_term)
    } else {
      stop(print(paste(msg,paste('The term defined is not class formula. \n \t Example:\n \t \t pheno_data = ~ X or pheno_data = ~ X + Y'))), call. = FALSE)
    }

    ## Check all variables in pheno_data term are factors
### ..x was because it is an pheno_data of data.table
    if(!all(sapply(all.vars(rand_fix_term), function(x, pheno_data) is.factor(pheno_data[, x]),  pheno_data))) {
      # stop("All variables indicated in argument 'pheno_data' should be factors")

    } else {

      rand = all.vars(rand_fix_term)
      for (r in 1:length(rand)) {

        pheno_data[, rand[r]] = as.factor(pheno_data[, rand[r]])

      }
    }

  rm(pheno_data)
  # Assign appropriate class.
  class(rand_fix_term) <- c('formula')

  attr(rand_fix_term, "cleared") <- "pass"


  return(rand_fix_term)



}
