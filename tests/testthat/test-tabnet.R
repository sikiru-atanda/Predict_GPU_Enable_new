# TabNet did not generalize on marker data (test r about 0 at 10-200 epochs):
# sparsemax kept about 2 of 400 markers per step, the output used only the
# masked inputs, and the sparsity penalty was detached (no gradient). It is now
# the reference architecture (input batch norm, shared/step GLU blocks,
# decision/attention split) with dense softmax attention, the reference sizes
# and its own learning rate.

test_that("TabNet fit parameters use the reference defaults and their own learning rate", {
  p <- PredictProR:::gp_dl_bridge_fit_params("tabnet", list())
  expect_identical(p$tabnet_steps, 3L)
  expect_identical(p$tabnet_feature_dim, 8L)
  expect_identical(p$tabnet_output_dim, 8L)
  expect_equal(p$tabnet_gamma, 1.3)
  expect_equal(p$tabnet_lambda_sparse, 1e-3)
  expect_equal(p$learning_rate, 2e-2)
  expect_equal(PredictProR:::gp_dl_bridge_fit_params("tabnet", list(learning_rate = 1e-3))$learning_rate, 1e-3)
  expect_equal(PredictProR:::gp_dl_bridge_fit_params("tabnet", list(tabnet_learning_rate = 5e-3))$learning_rate, 5e-3)
  expect_equal(eval(formals(PredictProR::model_execute)$tabnet_steps), 3L)
  expect_equal(eval(formals(PredictProR::model_execute)$tabnet_lambda_sparse), 1e-3)
})

test_that("tabnet_learning_rate has one default and is forwarded on the TabNet routes", {
  canon_default <- PredictProR:::gp_dl_bridge_fit_params("tabnet", list())$learning_rate
  expect_equal(eval(formals(PredictProR::model_execute)$tabnet_learning_rate), canon_default)
  for (fn in c("deep_learning_model", "gp_met_deep_learning_model")) {
    f <- get(fn, envir = asNamespace("PredictProR"))
    expect_equal(eval(formals(f)$tabnet_learning_rate), canon_default, info = fn)
  }
  # MET ML/DL cross-validation builds its fit parameters from the CV list
  met <- PredictProR:::gp_met_dl_model_params
  expect_equal(met("tabnet", list(tabnet_learning_rate = 5e-3))$tabnet_learning_rate, 5e-3)
  expect_equal(met("tabnet", list())$tabnet_learning_rate, canon_default)
  # single-environment cross-validation
  seen <- NULL
  local_mocked_bindings(deep_learning_model = function(..., tabnet_learning_rate) {
    seen <<- tabnet_learning_rate
    numeric(0)
  })
  PredictProR::predict_with_model("tabnet", y = 1:4, omics_data = matrix(0, 4, 2), tst = 1L,
                                  additional_params = list(tabnet_learning_rate = 5e-3))
  expect_equal(seen, 5e-3)
})

test_that("TabNet uses dense attention and its sparsity penalty has a gradient", {
  py <- skip_if_no_dl_python()
  py_dir <- system.file("python", package = "PredictProR")
  expect_true(file.exists(file.path(py_dir, "dl_models.py")),
              info = paste0("python dir: '", py_dir, "'; namespace path: ",
                            getNamespaceInfo("PredictProR", "path"), "; wd: ", getwd()))
  script <- withr::local_tempfile(fileext = ".py")
  writeLines(c(
    "import sys, torch",
    sprintf("sys.path.insert(0, r'%s')", py_dir),
    "import dl_models as d",
    "torch.manual_seed(1)",
    "m = d.TabNet(in_dim=40, out_dim=1, task='regression')",
    "out = m(torch.randn(16, 40))",
    "pen = m.extra_loss()",
    "print(m.mask_type, tuple(out.shape), bool(torch.is_tensor(pen) and pen.requires_grad))"
  ), script)
  err <- withr::local_tempfile(fileext = ".txt")
  out <- suppressWarnings(system2(py, script, stdout = TRUE, stderr = err))
  expect_identical(trimws(utils::tail(out, 1L)), "softmax (16, 1) True",
                   info = paste(c(py, out, utils::tail(readLines(err, warn = FALSE), 20L)), collapse = "\n"))
})
