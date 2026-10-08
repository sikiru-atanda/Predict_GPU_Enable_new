# Run-scoped preprocessing reuse for composite model_execute() workflows.
#
# A normal model_execute() call already prepares phenotype, genomic, omic and
# kernel inputs once before it dispatches model/fold tasks.  Composite
# workflows which call model_execute() once per protected route must carry one
# cache environment across those child calls or each child will repeat QC,
# imputation, LD pruning and kernel repair.  The cache is deliberately private,
# process-local and immutable-by-contract after each stage is stored.

gp_shared_preprocess_cache_new <- function(models = character(),
                                           requires_asreml = FALSE,
                                           requires_final_prediction = FALSE) {
  cache <- new.env(parent = emptyenv())
  cache$schema_version <- 1L
  cache$values <- new.env(parent = emptyenv())
  cache$models <- unique(as.character(models %||% character()))
  cache$requires_asreml <- isTRUE(requires_asreml)
  cache$requires_final_prediction <- isTRUE(requires_final_prediction)
  cache$builds <- integer()
  cache$hits <- integer()
  class(cache) <- c("predictpror_preprocess_cache", "environment")
  cache
}

gp_shared_preprocess_cache_validate <- function(cache) {
  if (is.null(cache)) return(invisible(NULL))
  valid <- is.environment(cache) &&
    inherits(cache, "predictpror_preprocess_cache") &&
    identical(cache$schema_version, 1L) &&
    is.environment(cache$values)
  if (!isTRUE(valid)) {
    stop(
      "The internal PredictProR preprocessing cache is invalid or incompatible.",
      call. = FALSE
    )
  }
  invisible(cache)
}

gp_shared_preprocess_cache_has <- function(cache, stage) {
  if (is.null(cache)) return(FALSE)
  gp_shared_preprocess_cache_validate(cache)
  exists(as.character(stage)[[1L]], envir = cache$values, inherits = FALSE)
}

gp_shared_preprocess_cache_get <- function(cache, stage) {
  if (!gp_shared_preprocess_cache_has(cache, stage)) return(NULL)
  stage <- as.character(stage)[[1L]]
  value <- get(stage, envir = cache$values, inherits = FALSE)
  count <- if (stage %in% names(cache$hits)) cache$hits[[stage]] else 0L
  cache$hits[[stage]] <- as.integer(count) + 1L
  value
}

gp_shared_preprocess_cache_set <- function(cache, stage, value) {
  if (is.null(cache)) return(invisible(value))
  gp_shared_preprocess_cache_validate(cache)
  stage <- as.character(stage)[[1L]]
  if (exists(stage, envir = cache$values, inherits = FALSE)) {
    stop(
      "PredictProR attempted to replace immutable shared preprocessing stage `",
      stage,
      "`.",
      call. = FALSE
    )
  }
  assign(stage, value, envir = cache$values)
  count <- if (stage %in% names(cache$builds)) cache$builds[[stage]] else 0L
  cache$builds[[stage]] <- as.integer(count) + 1L
  invisible(value)
}

gp_shared_preprocess_requested_models <- function(ctx) {
  shared_models <- ctx$predictpror_shared_models %||%
    ctx$.predictpror_shared_models %||%
    character()
  models <- unique(c(
    as.character(ctx$GS_model %||% character()),
    as.character(ctx$GS_model_cv %||% character()),
    as.character(shared_models)
  ))
  models <- models[!is.na(models) & nzchar(trimws(models))]
  gp_canonicalize_supported_model_names(models)
}

gp_shared_preprocess_summary <- function(cache) {
  if (is.null(cache)) {
    return(list(
      enabled = FALSE,
      models = character(),
      stages_built = integer(),
      stage_reuse_hits = integer()
    ))
  }
  gp_shared_preprocess_cache_validate(cache)
  list(
    enabled = TRUE,
    schema_version = cache$schema_version,
    models = cache$models,
    requires_asreml_order = isTRUE(cache$requires_asreml),
    requires_final_prediction = isTRUE(cache$requires_final_prediction),
    stages_built = cache$builds,
    stage_reuse_hits = cache$hits
  )
}
