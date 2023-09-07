#' Prepares kernel matrix
#'
#' If an object of class \code{genotype} is provided, calculates the
#' relationship matrix based on a requested method and perform diagnostics.
#' If the object is an user-provided GRM matrix, do basic checking and proceed
#' to diagnostics. If the object is of class \code{grm}, proceeds to
#' diagnostics.
#'
#' @export
#' @examples
#' snp <- as.matrix(asreml::nassau.snp) + 1
#' class(x = snp) <- c(class(x = snp), "genotype")
#'
#' snp |>
#'   prepare::grm(bend = T) |>
#'   head(c(5, 5))
#'
#' asreml::nassau.grm |>
#'   prepare::grm() |>
#'   head(c(5, 5))
#'

# TODO add kernel calculations from Sikiru.

grm <- function(object = NULL, message = TRUE, ...){

  # Traps -------------------------------------------------------------------

  # Check class of input.
  # If it is of class "genotype".
  if (inherits(x = object, what = "genotype"))
    obj.class <- "genotype"

  # If it is of class "kernel".
  else if (inherits(x = object, what = "kernel"))
    obj.class <- "kernel"

  # If it is of class "matrix" (non-checked square kernel).
  else if (
    (nrow(object) == ncol(object)) & inherits(x = object, what = "matrix"))
    obj.class <- "user.kernel"
  # TODO we might want to do some more checking here.

  # If none of the above.
  else
    stop("'object' must be of class 'genotype' or an user provided kernel",
         call. = FALSE)


  # Collect ellipsis -------------------------------------------------------

  # Logic: any arguments to the ASRgenomics functions can be passed here via the
  # ellipsis (...). However, we have to distinguish them as some ASRgenomics
  # functions do not have ellipsis and would fail if incorrect input is parsed.
  # Each object (e.g. G.matrix) will contain the arguments with respect to that
  # function parsed via '...'.

  # Assign.
  ellipsis <- list(...)

  # Split into functions.
  G.matrix <- ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::G.matrix)
  G.tuneup <- ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::G.tuneup)
  kinship.diagnostics <- ellipsis.subset(
    ellipsis = ellipsis, fun = ASRgenomics::kinship.diagnostics)

  # Get K ------------------------------------------------------------------

  # Report K matrix computation.
  if (obj.class == "genotype"){

    if (message) {
      message(col_blue("\nObtaining relationship matrix (K)."))
    }

    # Check if round is needed based on the fractional parts.
    if (!all(na.omit(c(object) %% 1 == 0)))
    {
      if (message){
        message("Converting genotypes to integers.")
      }

      # Round numbers.
      object <- round(x = object, digits = 0)
    }

    # Get kinship matrix.

    if (message) {
      message(paste0("Computing kernel matrix with ", G.matrix$method), "'s method.")
    }

    ASRgenomics:::silent_(
      code =
        object <- ASRgenomics::G.matrix(
          M = object,
          method = G.matrix$method, # User.
          na.string = NA, # It will be NA after qc.filtering.
          sparseform = FALSE, # Fixed.
          digits = 8 # Fixed.
        )$G
    )

  } else {
    if (message){
      message("Using provided kernel matrix.")
    }
  }


  # Tuning G with user-defined parameters (if any).
  # G.tuneup has to be called so changes are cumulative.
  if (G.tuneup$blend | G.tuneup$bend | G.tuneup$align){

    if (message) message(col_blue("\nImplementing tune-up on kernel matrix."))

    # TODO check if we are tunning up with A, probably not.
    # # Adapt A if necessary.
    # if (!is.null(G.tuneup$A)){
    #
    #   tmpGA <- ASRgenomics::match.G2A(
    #     A = G.tuneup$A, # User.
    #     G = object, # Computed.
    #     clean = TRUE, # Fixed.
    #     ord = TRUE, # Fixed.
    #     mism = FALSE) # Fixed.
    #
    #   # Assign objects.
    #   G.tuneup$A <- tmpGA$Aclean
    #   object <- tmpGA$Gclean
    # }

    object <- ASRgenomics::G.tuneup(
      G = object, # Computed.
      A = G.tuneup$A, # User.
      blend = G.tuneup$blend, # User.
      pblend = G.tuneup$pblend, # User.
      bend = G.tuneup$bend, # User.
      eig.tol = G.tuneup$eig.tol, # User.
      align = G.tuneup$align, # User.
      rcn = FALSE, # Fixed.
      digits = 8, # Fixed.
      sparseform = FALSE, # Fixed.
      determinant = FALSE, # Fixed.
      message = message # Fixed.
    )$Gb
  }


  # Diagnosing K ------------------------------------------------------------

  # Report K diagnosis.
  if (message) message(col_blue("\nPerforming diagnostics on kernel matrix."))

  # Checking if everything is OK with K.
  tmp.obj <- ASRgenomics::kinship.diagnostics(
    K = object, # Computed/provided.
    diagonal.thr.large = kinship.diagnostics$diagonal.thr.large, # User.
    diagonal.thr.small = kinship.diagnostics$diagonal.thr.small, # User.
    duplicate.thr = kinship.diagnostics$duplicate.thr, # User.
    clean.diagonal = kinship.diagnostics$clean.diagonal, # User.
    clean.duplicate = kinship.diagnostics$clean.duplicate, # User.
    message = message
  )

  # Clean up K if requested by user.
  if (!is.null(tmp.obj$clean.kinship)){

    if (nrow(tmp.obj$clean.kinship) == 0){
      stop("No genotypes remain after cleaning based on diagonal values and/or duplicates.")
    }

    # Reporting clean up.
    if (message) message("Cleaning kernel matrix based on diagonal values and/or duplicates.")

    # Return clean K if changes were made.
    object <- tmp.obj$clean.kinship
  }

  # Removing objects no longer required.
  rm(tmp.obj)

  # -------------------------------------------------------------------------

  # Assign appropriate class.
  class(object) <-c("matrix", "array", "kernel")

  # Return.
  return(object)
}
