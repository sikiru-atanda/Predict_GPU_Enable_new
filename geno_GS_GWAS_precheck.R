#' Prepares genotype / genomic data
#'
#' @param object
#' @param rename_markers
#' @param message
#' @param map_data
#' @param ...
#'
#' @export
#' @examples
#' snp <- as.matrix(asreml::nassau.snp) + 1
#'
#' snp |>
#'   prepare::genotype() |>
#'   head(c(5, 5))
#'
#' snp |>
#'   prepare::genotype(maf = 0.05, ind.callrate = 0.3, pruning.thr = 0.95) |>
#'   head(c(5, 5))
#'
#' snp |>
#'   prepare::genotype(message = F) |>
#'   class()
#'

genotype <- function(
    object = NULL,
    rename_markers = FALSE,
    message = TRUE,
    map_data = NULL,
    ...){

  # Traps -------------------------------------------------------------------

  # TODO check what inputs to allow here.
  # TODO check which traps we want here.

  # Check class of input.
  if(!inherits(x = object, what = "matrix"))
    stop("'object' must be of class 'matrix'", call. = FALSE)

  # Check row and column names in object.
  if(is.null(rownames(object))){
    stop("Individual names not assigned to rows of \'object'.")
  }

  if(is.null(colnames(object))){
    stop("Marker names not assigned to columns of \'object'.")
  }

  # Collect ellipsis -------------------------------------------------------

  # See logic on grm.R.

  # Assign.
  ellipsis <- list(...)

  # Overwrite ASRgenomics defaults for this specific package.
  if (is.null(ellipsis$pruning.thr)){
    ellipsis$pruning.thr = FALSE
  }

  # Split into functions.
  qc_filtering <- ellipsis_subset(
    ellipsis = ellipsis, fun = ASRgenomics::qc.filtering)
  snp_pruning <- ellipsis_subset(
    ellipsis = ellipsis, fun = ASRgenomics::snp.pruning)

  # Fix marker names if necessary ------------------------------------------

  # Check if fix is requested.
  if (rename_markers){

    # This means it will look for - * / % and +
    math_operators <- "\\-|\\*|\\/|\\%|\\+|\\:"

    if (message){
      message(col_blue("\nVerifying if marker names need to be modified."))
    }

    # Check if there are math operators in the colnames of object.
    if (any(grepl(pattern = math_operators, x = colnames(object)))){

      if(message){
        message("Marker names in \'object' contain mathematical operators",
                "(e.g., + - * / % :), which will be replaced by \'_'.")
      }

      # Fix the names.
      colnames(object) <-
        sapply(X = colnames(object), FUN = gsub,
               pattern = math_operators, replacement = "_", USE.NAMES = FALSE)
    }

    if (!is.null(map_data)){

      # Check if there are math operators in the map_data$marker.
      if (any(grepl(pattern = math_operators, map_data$marker))){

        if(message){
          message("Marker names in \'map_data' contain mathematical operators",
                  "(e.g., + - * / % :), which will be replaced by \'_'.")
        }

        # Fix the names.
        map_data$marker <-
          sapply(X = map_data$marker, FUN = gsub,
                 pattern = math_operators, replacement = "_")
      }
    }
  }

  # Marker quality control --------------------------------------------------

  # Use qc.filtering for object quality control.
  # Perform quality control.
  if (message) {
    message("\nQuality control of data frame \'object' in progress.")
  }

  tmp_qc <- ASRgenomics::qc.filtering(
    M = object, # User.
    map = map_data, # User.
    marker = "marker", # Fixed.
    chrom = "chrom", # Fixed.
    pos = "pos", # Fixed.
    base = qc_filtering$base, # User.
    ref = qc_filtering$ref, # User.
    marker.callrate = qc_filtering$marker.callrate, # User.
    ind.callrate = qc_filtering$ind.callrate, # User.
    maf = qc_filtering$maf, # User.
    heterozygosity = qc_filtering$heterozygosity, # User.
    Fis = qc_filtering$Fis, # User.
    na.string = qc_filtering$na.string, # User.
    impute = qc_filtering$impute, # User.
    message = message, # User.
    Mrecode = FALSE,
    plots = FALSE, # User.
    digits = 2 # Fixed.
  )

  #### Meta data
  marker_callrate <- qc_filtering$marker.callrate
  ind_callrate <-  qc_filtering$ind.callrate
  maf <-  qc_filtering$maf
  heterozygosity <-  qc_filtering$heterozygosity
  Fis <-  qc_filtering$Fis
  Impute_method <- "mean"
  # Collect output (object).
  object <- tmp_qc$M.clean

  # Collect output (map_data) if it was provided by the user.
  if (!is.null(map_data)){
    map_data <- tmp_qc$map
  }

  # GC.
  rm(tmp_qc)


  # Pruning ----------------------------------------------------------------

  # Logic. Prune as required.
  # The map is not exported because it is not needed from the following code.
  # The map can be matched with object later on.
  # isFALSE is required here because the value can be a number or FALSE (cannot be NULL).

  if (!isFALSE(ellipsis$pruning.thr)){
    object <- ASRgenomics::snp.pruning(
      M = object,
      map = map_data,
      marker = "marker",
      chrom = "chrom",
      pos = "pos",
      method = "correlation",  # CONFLICT (maybe change in ASRgenomics in the future).
      criteria = snp_pruning$criteria,
      pruning.thr = snp_pruning$pruning.thr,
      by.chrom = snp_pruning$by.chrom,
      window.n = snp_pruning$window.n,
      overlap.n = snp_pruning$overlap.n,
      iterations = snp_pruning$iterations,
      seed = snp_pruning$seed,
      message = message
    )$Mpruned
  }

  # Check object for missing values ----------------------------------------

  if (message) {
    message(col_blue("\nChecking data frame \'object' for missing data."))
  }

  # Check if there is any missing values in object (required for the choice of Q.method).
  any.miss.object <- any(is.na(object))

  if (message & any.miss.object) {
    message("Missing data present in \'object'. This will impact in processing time of downstream analyses.")
    message("Use \'impute = TRUE' to request mean-based imputation.")
  }


  # Finalize ----------------------------------------------------------------

  # Assign appropriate class.
  class(object) <-c("matrix", "array", "genotype")

  marker_meta_data = list(marker_callrate, ind_callrate,
                    maf, heterozygosity, Impute_method)

  names(marker_meta_data) <- c("marker_callrate", "ind_callrate",
                           "maf", "heterozygosity", "Impute_method")

  output <-  list(object, map_data, marker_meta_data)

  names(output) <-  c("marker_matrix", "map", "marker_meta_data")


  # Return.
  return(output)
}
