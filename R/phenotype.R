#' Prepare phenotypic data: coerce column types and check NAs
#'
#' Casts named columns to numeric or factor, ensures specified response columns
#' contain no missing values, and returns the cleaned phenotype frame.
#'
#' @param object A phenotype data frame.
#' @param to.numeric Character vector of column names to coerce to numeric.
#' @param to.factor Character vector of column names to coerce to factor.
#' @param noNAon Character vector of column names that must not contain `NA`;
#'   `phenotype()` errors if any are missing.
#' @param message Logical; if `TRUE`, print summary messages about coercions
#'   and NA checks.
#'
#' @return A phenotype data frame with the requested coercions applied.
#' @export
#'
#' @examples
#'
#' pheno <- data.frame(
#'   GID = paste0("g", 1:4),
#'   Rep = c(1, 1, 2, 2),
#'   IncBlock = c("1", "2", "1", "2"),
#'   ht6 = c(12.1, 11.4, 13.0, 10.8),
#'   weight = c(NA, 2.3, 2.8, 2.1)
#' )
#'
#' pheno_checked <- pheno |>
#'   phenotype(
#'     to.factor = c("Rep"),
#'     to.numeric = c("IncBlock"),
#'     noNAon = c("ht6", "weight"))
#'
#' pheno_checked |> str()
#'
#' pheno_checked |> class()
#'

phenotype <- function(
    object = NULL,
    to.numeric = NULL,
    to.factor = NULL,
    noNAon = NULL,
    message = TRUE)
{

  # Traps -------------------------------------------------------------------

  # TODO check what inputs to allow here.
  # Check class of input.
  # if(!inherits(x = object, what = "data.frame"))
  #   stop("'phenotype' must be of class 'data.frame'", call. = FALSE)

  #object = object[[1]]

  if(!inherits(x = object, what = "data.frame"))
    stop("'phenotype' must be of class 'data.frame'", call. = FALSE)

  # Body --------------------------------------------------------------------

  # Transform into data.table.
  if (!data.table::is.data.table(object))
    object <- data.table::as.data.table(object)

  # Report.
  if (message) {
    #message(col_blue("\nChecking phenotypic data frame."))
    message("\nChecking phenotypic data frame.")
  }

  # Structure ---------------------------------------------------------------

  # Mutate variables to requested class.
  if (!is.null(to.factor)) {
    factor_values <- lapply(as.data.frame(object[, to.factor, with = FALSE]), FUN = factor)
    for (nm in to.factor) {
      data.table::set(object, j = nm, value = factor_values[[nm]])
    }
  }

  if (!is.null(to.numeric)) {
    numeric_values <- lapply(
      as.data.frame(object[, to.numeric, with = FALSE]),
      function(x) as.numeric(as.character(x))
    )
    for (nm in to.numeric) {
      data.table::set(object, j = nm, value = numeric_values[[nm]])
    }
  }

  # Remove missing ----------------------------------------------------------

  # Remove NA on noNAon.
  if (!is.null(noNAon)){

    # TODO choose which one to use and drop the other ones. Probably use data.table as it is more efficient.
    # In case object is a data.table.
    if(data.table::is.data.table(object)){
      not.na.cases <- stats::complete.cases(as.data.frame(object[, noNAon, with = FALSE]))
    } # else {
    #   # In case object is data.frame, tibble, etc.
    #   not.na.cases <- complete.cases(object[, noNAon])
    # }

    # Remove if needed.
    if (any(!not.na.cases)){
      # Apply removal.
      object <- object[not.na.cases,]

      # Report change.
      if (message){
        message("A total of ", sum(!not.na.cases),
                " sample(s) were removed due to missing values in ",
                paste0(noNAon, collapse = " and/or "), ".")
      }
    }

    # GC.
    rm(not.na.cases)
  }

  # Remove any ghost levels (once just in case).
  object <- droplevels(object)

  # Assign appropriate class.
  class(object) <- c("data.table", "data.frame", "phenotype")

  attr(object, "data.table") <- "passed"

  if(isTRUE(attr(object, "passed"))) print('ok')

  # Return.
  return(object)
}
