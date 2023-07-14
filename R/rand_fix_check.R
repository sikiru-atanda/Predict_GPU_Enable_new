
#' Title
#'
#' @param object
#' @param ...
#' @param rand_fix_term
#'
#' @return
#' @export
#'
#' @examples
rand_fix_check <- function(object = NULL,
                      rand_fix_term = NULL,
                      ...){
  msg <- sprintf("==================================================\n")

  # if(attr(object, "cleared")!="pass" & class(object)!=c("data.table", "data.frame", "phenotype")) {
  #
  #   stop(print(paste(msg,paste(object, 'is not class phenotype.', sep = ""))), call. = FALSE)
  # }
    if (inherits(rand_fix_term, what = 'formula')) {
      rand_fix_term <- stats::as.formula(rand_fix_term)
    } else {
      stop(print(paste(msg,paste('The term defined is not class formula. \n \t Example:\n \t \t object = ~ X or object = ~ X + Y'))), call. = FALSE)
    }

    ## Check all variables in object term are factors
### ..x was because it is an object of data.table
    if(!all(sapply(all.vars(rand_fix_term), function(x, object) is.factor(object[, x]),  object))) {
      # stop("All variables indicated in argument 'object' should be factors")

    } else {

      rand = all.vars(rand_fix_term)
      for (r in 1:length(rand)) {

        object[, rand[r]] = as.factor(object[, rand[r]])

      }
    }

  rm(object)
  # Assign appropriate class.
  class(rand_fix_term) <- c('formula')

  attr(rand_fix_term, "cleared") <- "pass"


  return(rand_fix_term)



}
