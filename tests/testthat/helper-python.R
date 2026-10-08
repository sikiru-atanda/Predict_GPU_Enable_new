# Skip tests that need the Python ML bridge when no interpreter with
# numpy + scikit-learn is available (e.g. a bare system python3 on Linux).
skip_if_no_ml_python <- function() {
  py <- tryCatch(PredictProR:::gp_detect_ml_python(), error = function(e) NULL)
  testthat::skip_if(is.null(py) || !nzchar(py) || !file.exists(py), "No ML Python runtime with numpy/scikit-learn")
  probe <- tryCatch(PredictProR:::gp_python_probe(py, modules = c("numpy", "sklearn")), error = function(e) NULL)
  testthat::skip_if(
    !isTRUE(probe[["numpy"]]) || !isTRUE(probe[["sklearn"]]),
    "ML Python runtime lacks numpy/scikit-learn"
  )
  invisible(py)
}
