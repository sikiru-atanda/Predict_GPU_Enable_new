#' Prepares phenotypic data
#' @export
#'
#' @examples
#'
#' # Messing with the data.
#' nassau <- asreml::nassau
#' nassau$weight <- 1 ; nassau$weight[1] <- NA
#' nassau$Rep <- as.numeric(nassau$Rep)
#'
#' nassau. <- nassau |>
#'   phenotype(
#'     to.factor = c("Rep"),
#'     to.numeric = c("IncBlock"),
#'     noNAon = c("ht6", "weight"))
#'
#' nassau. |> str()
#'
#' nassau. |> class()
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

  if(!inherits(x = object, what = "data.table"))
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
  if (!is.null(to.factor))
    object[, to.factor] <-
      lapply(object[, ..to.factor, drop = FALSE], FUN = factor)

  if (!is.null(to.numeric))
    object[, to.numeric] <-
      lapply(object[, ..to.numeric, drop = FALSE],
             function(x) as.numeric(as.character(x)))

  # Remove missing ----------------------------------------------------------

  # Remove NA on noNAon.
  if (!is.null(noNAon)){

    # TODO choose which one to use and drop the other ones. Probably use data.table as it is more efficient.
    # In case object is a data.table.
    if(data.table::is.data.table(object)){
      not.na.cases <- complete.cases(object[, ..noNAon])
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
