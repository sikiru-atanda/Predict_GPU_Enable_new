# The multi-trait ASReml routes checked for asreml with installed.packages(),
# which reads every installed package (0.36-0.39 s in a session's first fit
# on the Windows test machine, 14.6 s on a cold disk, again in every worker).
# Every route now checks with system.file() and loads the namespace, without
# attaching it.

test_that("no package code scans installed.packages()", {
  files <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(!length(files), "package sources not available")
  hits <- unlist(lapply(files, function(f) {
    code <- sub("#.*$", "", readLines(f, warn = FALSE))
    if (any(grepl("installed.packages(", code, fixed = TRUE))) basename(f)
  }))
  expect_identical(hits, NULL)
})

test_that("every ASReml route stops with the install message when the engine is missing", {
  expect_error(PredictProR:::gp_multitrait_asreml_gaussian_model(pheno_object = NULL, engine = "notapackage"),
               "install asreml-R to use multi-trait ASReml-R")
  expect_error(PredictProR:::gp_multitrait_asreml_gaussian_cv(pheno_object = NULL, engine = "notapackage"),
               "install asreml-R to use multi-trait ASReml-R CV")
  expect_error(PredictProR:::gp_multitrait_asreml_mtmet_gaussian_model(pheno_object = NULL, engine = "notapackage"),
               "install asreml-R to use grouped MT-MET GBLUP")
  expect_error(suppressMessages(capture.output(
    PredictProR::asreml_utilis_new(pheno_data = data.frame(GID = "a", y = 1), gen_name = "GID", response = "y",
                                   engine = "notapackage"))),
    "install asreml-R")
})

test_that("the engine check loads asreml without attaching it", {
  skip_if_not_installed("asreml")
  attached <- "package:asreml" %in% search()
  if (attached) {
    detach("package:asreml")
    withr::defer(suppressPackageStartupMessages(library(asreml)))
  }
  expect_invisible(PredictProR:::gp_asreml_require_engine("asreml", "install asreml-R"))
  expect_true("asreml" %in% loadedNamespaces())
  expect_false("package:asreml" %in% search())
})
