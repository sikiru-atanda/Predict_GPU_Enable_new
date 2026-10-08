
#' Normalize testing-set identifiers
#'
#' Validates and normalizes a vector, one-column data frame, or one-column
#' matrix of testing identifiers.
#'
#' @param test_set Testing-set identifiers as a vector, one-column data frame,
#'   or one-column matrix.
#'
#' @return A unique vector of testing identifiers, or \code{NULL}.
#' @export
check_test_set <- function(test_set = NULL
) {
  msg <- ""
  if (is.null(test_set)) return(NULL)

  if (!is.data.frame(test_set) && !is.matrix(test_set) && !is.vector(test_set)) {
    stop(paste(msg, 'The testing set should be a dataframe, matrix, or a vector.'), call. = FALSE)
  }

  if (is.data.frame(test_set) || is.matrix(test_set)) {
    if(ncol(test_set)>1){
      stop(paste(msg, 'The testing set should be a dataframe, matrix with single column.'), call. = FALSE)
    }
    test_set <- unique(test_set[, 1])
  } else {
    test_set <- unique(test_set)
  }

  if (length(test_set) == 0) {
    return(NULL)
  }

  return(test_set)
}

#' Match phenotype IDs with genomic or omic data
#'
#' Aligns genomic, omic, relationship, or kernel data to phenotype records.
#' Phenotyped IDs missing from genomic/omic data stop the workflow. Extra
#' genomic/omic IDs are treated as inferred testing IDs unless an explicit
#' test set is configured to override inferred testing IDs.
#'
#' @param object_geno Matrix-like genomic/omic data, relationship matrix, or
#'   kernel matrix.
#' @param object_pheno Phenotype data frame.
#' @param gen_name Name of the genotype or sample identifier column in
#'   \code{object_pheno}.
#' @param train_set Optional training-set identifiers.
#' @param test_set Optional testing-set identifiers.
#' @param heter_groups Optional environment or grouping column used for ordering.
#' @param low_call_rate_inds_removed Optional IDs removed during low-call-rate
#'   QC.
#' @param test_set_overrides_inferred Logical. When \code{TRUE}, an explicit
#'   \code{test_set} is used as supplied and genomic/omic-only individuals are
#'   not added automatically.
#' @param message Logical indicating whether progress messages should be
#'   printed.
#'
#' @return A list containing \code{geno_pheno_match_data} and, when present,
#'   \code{test_data}.
#' @export
pheno_geno_match <- function(object_geno = NULL,
                             object_pheno = NULL,
                             gen_name = NULL,
                             train_set = NULL,
                             test_set = NULL,
                             heter_groups = NULL,
                             low_call_rate_inds_removed = NULL,
                             test_set_overrides_inferred = FALSE,
                             message = TRUE) {
  msg <- ""

  if (is.null(object_geno)) {
    stop(paste(msg, "Genotypic/omic data cannot be empty."), call. = FALSE)
  }
  if (is.null(object_pheno) || is.null(gen_name) || !gen_name %in% names(object_pheno)) {
    stop(paste(msg, "Phenotypic data and a valid genotype ID column are required."), call. = FALSE)
  }
  if (is.null(rownames(object_geno))) {
    stop(paste(msg, "Genotypic/omic data must have row names matching individual IDs."), call. = FALSE)
  }

  # Check if genotypes are consistent across all environments
  if(!is.null(heter_groups)){
    object_pheno <- object_pheno |>
      dplyr::arrange(!!rlang::sym(heter_groups), !!rlang::sym(gen_name))

  }

  qc_removed_ids <- check_test_set(low_call_rate_inds_removed)
  if(!is.null(qc_removed_ids)){
    qc_removed_ids <- as.character(qc_removed_ids)
    pheno_ids_before_qc_drop <- as.character(object_pheno[[gen_name]])
    removed_pheno_ids <- unique(pheno_ids_before_qc_drop[pheno_ids_before_qc_drop %in% qc_removed_ids])
    if (length(removed_pheno_ids) > 0 && isTRUE(message)) {
      message(insight::print_color(
        paste(
          msg,
          "Removing phenotypic records for individuals removed during genomic/omic QC:",
          length(removed_pheno_ids),
          "IDs:",
          paste(utils::head(removed_pheno_ids, 20), collapse = ", ")
        ),
        "blue"
      ))
    }
    object_pheno <- object_pheno[!pheno_ids_before_qc_drop %in% qc_removed_ids, , drop = FALSE]
    test_set <- check_test_set(setdiff(check_test_set(test_set), qc_removed_ids))
    train_set <- check_test_set(setdiff(check_test_set(train_set), qc_removed_ids))

  }
  ID_pheno <- as.character(unique(object_pheno[[gen_name]]))
  ID_pheno <- ID_pheno[!is.na(ID_pheno)]

  geno_ids <- as.character(rownames(object_geno))
  is_kernel <- is.matrix(object_geno) &&
    nrow(object_geno) == ncol(object_geno) &&
    !is.null(colnames(object_geno)) &&
    setequal(as.character(rownames(object_geno)), as.character(colnames(object_geno)))

  missing_pheno_ids <- setdiff(ID_pheno, geno_ids)
  if (length(missing_pheno_ids) > 0) {
    stop(
      paste(
        msg,
        "Not all individuals with phenotypic records have genotypic/omic records.",
        "Missing IDs:",
        paste(utils::head(missing_pheno_ids, 20), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (isTRUE(is_kernel)) {
    missing_kernel_cols <- setdiff(ID_pheno, as.character(colnames(object_geno)))
    if (length(missing_kernel_cols) > 0) {
      stop(
        paste(
          msg,
          "Not all individuals with phenotypic records are present in the kernel columns.",
          "Missing IDs:",
          paste(utils::head(missing_kernel_cols, 20), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  explicit_test_set <- check_test_set(test_set)
  train_set <- check_test_set(train_set)
  if (!is.null(train_set)) {
    missing_train_ids <- setdiff(train_set, ID_pheno)
    if (length(missing_train_ids) > 0) {
      stop(
        paste(
          msg,
          "All train_set IDs must be present in the phenotypic data.",
          "Unknown IDs:",
          paste(utils::head(missing_train_ids, 20), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  test_from_train <- if (!is.null(train_set)) setdiff(ID_pheno, train_set) else NULL
  genomic_only_test_set <- setdiff(geno_ids, ID_pheno)
  if (length(genomic_only_test_set) > 0) {
    if (isTRUE(message)) {
      message(insight::print_color(paste(msg, 'Not all individuals with genotypic/omic records have phenotypic records.'), "blue"))
    }
  }

  infer_genomic_only_test_set <- !isTRUE(test_set_overrides_inferred) || is.null(explicit_test_set)
  resolved_test_set <- check_test_set(unique(c(
    explicit_test_set,
    test_from_train,
    if (isTRUE(infer_genomic_only_test_set)) genomic_only_test_set else NULL
  )))
  if (!is.null(resolved_test_set)) {
    missing_test_ids <- setdiff(resolved_test_set, geno_ids)
    if (length(missing_test_ids) > 0) {
      stop(
        paste(
          msg,
          "All testing IDs must have genotypic/omic records.",
          "Missing IDs:",
          paste(utils::head(missing_test_ids, 20), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  model_ids <- unique(c(
    ID_pheno,
    if (isTRUE(infer_genomic_only_test_set)) genomic_only_test_set else NULL,
    setdiff(resolved_test_set, ID_pheno)
  ))
  model_ids <- model_ids[model_ids %in% geno_ids]

  if (isTRUE(is_kernel)) {
    object_geno <- object_geno[
      match(model_ids, rownames(object_geno)),
      match(model_ids, colnames(object_geno)),
      drop = FALSE
    ]
  } else {
    object_geno <- object_geno[
      match(model_ids, rownames(object_geno)),
      ,
      drop = FALSE
    ]
  }

  attr(object_geno, "cleared") <- "model_ready_use"
  out <- list(geno_pheno_match_data = object_geno)
  if (!is.null(resolved_test_set)) {
    out[["test_data"]] <- resolved_test_set
  }
  out

}

