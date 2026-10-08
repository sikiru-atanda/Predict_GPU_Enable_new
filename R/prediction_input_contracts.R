#' Describe PredictProR input standards through one interface
#'
#' @param task One of `"all"`, `"general"`, `"hybrid"`, `"multi_trait"`,
#'   or `"multi_environment"`.
#'
#' @return The selected input standard, or a named list of every standard when
#'   `task = "all"`.
#' @export
prediction_input_standard <- function(task = c(
                                        "all", "general", "hybrid",
                                        "multi_trait", "multi_environment"
                                      )) {
  task <- match.arg(task)
  standards <- list(
    general = general_prediction_data_standard(),
    hybrid = hybrid_data_standard(),
    multi_trait = multi_trait_data_standard(),
    multi_environment = met_data_standard()
  )
  if (identical(task, "all")) standards else standards[[task]]
}

gp_prediction_input_validator <- function(task) {
  switch(
    task,
    general = validate_general_input_standard,
    hybrid = validate_hybrid_input_standard,
    multi_trait = validate_multi_trait_input_standard,
    multi_environment = validate_met_input_standard,
    stop("Unknown prediction input task.", call. = FALSE)
  )
}

gp_prediction_input_infer_task <- function(args) {
  true_flag <- function(name) isTRUE(args[[name]])
  if (any(vapply(
    c("hybrid_asreml", "hybrid_bayes", "hybrid_gp", "hybrid_ml", "hybrid_dl"),
    true_flag,
    logical(1)
  )) || any(c("female_parent", "male_parent", "female_gmatrix", "male_gmatrix") %in% names(args))) {
    return("hybrid")
  }
  if (any(vapply(
    c("multi_trait_asreml", "multi_trait_bayes", "multi_trait_gp", "multi_trait_ml", "multi_trait_dl"),
    true_flag,
    logical(1)
  )) || length(args$response %||% character()) > 1L) {
    return("multi_trait")
  }
  if (isTRUE(args$met_ml_dl)) {
    return("multi_environment")
  }
  "general"
}

#' Validate PredictProR input through one interface
#'
#' Dispatches to the same workflow-specific validator used by
#' `model_execute()`. With `task = "auto"`, hybrid flags take precedence,
#' followed by multi-trait flags or multiple responses, MET ML/DL, and then
#' the general contract.
#'
#' @param task One of `"auto"`, `"general"`, `"hybrid"`, `"multi_trait"`,
#'   or `"multi_environment"`.
#' @param ... Named arguments accepted by the selected workflow validator.
#'
#' @return Invisibly returns the validator result with the resolved `task` and
#'   validator name.
#' @rdname prediction_input_standard
#' @export
validate_prediction_input <- function(task = c(
                                        "auto", "general", "hybrid",
                                        "multi_trait", "multi_environment"
                                      ),
                                      ...) {
  task <- match.arg(task)
  args <- list(...)
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop("All validation arguments must be named.", call. = FALSE)
  }
  resolved_task <- if (identical(task, "auto")) gp_prediction_input_infer_task(args) else task
  validator <- gp_prediction_input_validator(resolved_task)
  accepted <- names(formals(validator))
  unknown <- setdiff(names(args), accepted)
  if (length(unknown)) {
    stop(
      sprintf(
        "Arguments not valid for the %s input contract: %s.",
        resolved_task,
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  result <- do.call(validator, args)
  result <- result %||% list(valid = TRUE)
  result$task <- resolved_task
  result$validator <- switch(
    resolved_task,
    general = "validate_general_input_standard",
    hybrid = "validate_hybrid_input_standard",
    multi_trait = "validate_multi_trait_input_standard",
    multi_environment = "validate_met_input_standard"
  )
  invisible(result)
}
