# Skip tests that need a Python backend when its interpreter is missing or
# lacks the packages (e.g. a bare system python3 on Linux). The test suite
# never installs them (setup-python-install.R).
skip_if_python_lacks <- function(py, modules, what) {
  testthat::skip_if(is.null(py) || !nzchar(py) || !file.exists(py), paste("No", what, "Python runtime"))
  probe <- tryCatch(PredictProR:::gp_python_probe(py, modules = modules), error = function(e) NULL)
  lacking <- modules[!vapply(modules, function(m) isTRUE(probe[[m]]), logical(1))]
  testthat::skip_if(length(lacking) > 0L,
                    paste0(what, " Python runtime lacks ", paste(lacking, collapse = "/")))
  invisible(py)
}

skip_if_no_ml_python <- function() {
  py <- tryCatch(PredictProR:::gp_detect_ml_python(), error = function(e) NULL)
  skip_if_python_lacks(py, c("numpy", "sklearn"), "ML")
}

skip_if_no_dl_python <- function() {
  py <- tryCatch(PredictProR:::gp_detect_dl_python(), error = function(e) NULL)
  skip_if_python_lacks(py, c("numpy", "pandas", "torch"), "deep-learning")
}

skip_if_no_gp_python <- function(py = tryCatch(PredictProR:::gp_detect_gp_python(), error = function(e) NULL)) {
  skip_if_python_lacks(py, c("numpy", "torch", "gpytorch"), "Gaussian-process")
}
