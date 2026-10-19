# Preserve the estimator settings in Run_metadata as well as model_parameters.
gp_regularization_metadata <- function(x) {
  found <- list()
  visit <- function(one, path = character(), depth = 0L) {
    if (!is.list(one) || is.data.frame(one) || depth > 6L) return(invisible(NULL))
    p <- one$model_parameters
    if (is.data.frame(p) && all(c("stat", "summary") %in% names(p))) {
      p <- p[grepl("^gp_krr_lambda", p$stat), , drop = FALSE]
      for (i in seq_len(nrow(p))) {
        key <- paste(c(path, as.character(p$stat[i])), collapse = ".")
        found[[key]] <<- as.character(p$summary[i])
      }
      return(invisible(NULL))
    }
    for (nm in names(one)) visit(one[[nm]], c(path, nm), depth + 1L)
  }
  visit(x)
  found
}
