#' Subset ellipsis based on a functions formals.
#'
#' @description
#' Sometimes we want the ellipsis to bear the input for more than one internal function.
#' Therefore we can split it according to the formals of each function.
#'
#' @param ellipsis The ellipsis (i.e., ...) (default = \code{NULL}).
#' @param fun The function to check formals and collect respective input from ... (default = \code{NULL}).
#'
#' @return A list with the ... corresponding to the function.
#'
#' @keywords internal

ellipsis_subset <- function(ellipsis = NULL, fun = NULL){

  # Collect function formals (arguments values and names).
  vals <- formals(fun)
  args <- formalArgs(fun)

  # Run across vals and get first value is multiple are used (e.g., c(1,2)).
  # Logic, in the case c() is in the input, it will be x[[1]].
  vals <- lapply(vals, function(x) if (is.call(x)) return(x[[2]]) else x)

  # Identify arguments passed using ellipsis.
  ellipsis_subset <- names(ellipsis)

  # Identify which arguments from ellipsis match the arguments from function.
  vars2get <- match(args, ellipsis_subset)

  # Remove eventual NAs.
  vars2get <- na.exclude(vars2get)

  # Replace function's formals with values provided with ellipsis.
  vals[args %in% ellipsis_subset] <- ellipsis[vars2get]

  return(vals)
}
