# Deep-learning epochs. model_execute() trained every DL model for 10 epochs by
# default, while the internal DL functions defaulted to 30 (hybrid: 20). A fixed
# epoch count also ignores data size: a 100-line panel gets ~2 optimizer steps
# per epoch, a 5000-line panel ~60. epochs = "auto" (one default everywhere)
# lets the Python trainer ask for ~200 steps, kept within 20-100 epochs.

test_that("deep-learning epochs default to \"auto\" everywhere they are declared", {
  expect_identical(eval(formals(PredictProR::model_execute)$epochs), "auto")
  expect_identical(PredictProR:::normalize_common(list())$epochs, "auto")
  for (fn in c("deep_learning_model", "gp_met_deep_learning_model",
               "gp_multitrait_dl_gaussian_model", "gp_multitrait_dl_gaussian_cv",
               "gp_hybrid_dl_fit_response", "gp_hybrid_dl_gaussian_model", "gp_hybrid_dl_gaussian_cv")) {
    f <- get(fn, envir = asNamespace("PredictProR"))
    expect_identical(eval(formals(f)$epochs), "auto", info = fn)
  }
  expect_identical(PredictProR:::gp_met_dl_model_params("mlp", list())$epochs, "auto")
  # an explicit number is kept
  expect_identical(PredictProR:::gp_dl_bridge_fit_params("mlp", list(epochs = 12))$epochs, 12L)
})

test_that("epochs accepts \"auto\" or a positive whole number", {
  arg <- PredictProR:::gp_dl_epochs_arg
  expect_identical(arg("auto"), "auto")
  expect_identical(arg(" AUTO "), "auto")
  expect_identical(arg(NULL), "auto")
  expect_identical(arg(40), 40L)
  for (bad in list(0, -5, 2.5, "lots", c(10, 20), NA)) {
    expect_error(arg(bad), "positive whole number")
  }
})

test_that("the Python trainer sizes epochs = \"auto\" by the steps per epoch", {
  py <- skip_if_no_dl_python()
  script <- withr::local_tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    sprintf("sys.path.insert(0, r'%s')", system.file("python", package = "PredictProR")),
    "import dl_models as d",
    # steps per epoch 1, 2 (100 lines), 4 (250), 7 (500), 25 (2000); explicit 12
    "print([d._resolve_epochs('auto', s) for s in (1, 2, 4, 7, 25)], d._resolve_epochs(12, 2))"
  ), script)
  out <- suppressWarnings(system2(py, script, stdout = TRUE, stderr = FALSE))
  expect_identical(trimws(utils::tail(out, 1L)), "[100, 100, 50, 29, 20] 12")
})
