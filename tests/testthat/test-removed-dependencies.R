# Packages and code PredictProR no longer uses. mgcv was imported only for
# the GAM helpers (gam_new(), ND_modes_cv(), ...), which nothing called;
# they were removed in 0.20.207. Every import is loaded by the parallel
# policy's worker-memory check. The older random-forest feature-selection
# engine went in 0.20.210: feature selection is feature_scoring.R.

test_that("mgcv is neither imported nor called", {
  desc <- test_path("..", "..", "DESCRIPTION")
  files <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(!file.exists(desc) || !length(files), "package sources not available")
  imports <- read.dcf(desc, fields = "Imports")[1, "Imports"]
  imports <- trimws(sub("[(].*$", "", strsplit(imports, ",")[[1]]))
  expect_false("mgcv" %in% imports)
  hits <- unlist(lapply(files, function(f) {
    code <- sub("#.*$", "", readLines(f, warn = FALSE))
    idx <- grep("mgcv::", code, fixed = TRUE)
    if (length(idx)) paste0(basename(f), ":", idx) else character()
  }))
  expect_identical(hits, character())
})

test_that("the older random-forest feature-selection engine is not part of the package", {
  # it ranked by how often a marker's importance was above zero and chose the
  # set size by RMSE on the rows it was fitted to
  ns <- asNamespace("PredictProR")
  for (fn in c("feature_selection_rf", "cv_with_replication_stability",
               "sequential_feature_selection", "set_parallel_plan_rf")) {
    expect_false(exists(fn, envir = ns, inherits = FALSE), info = fn)
  }
})
