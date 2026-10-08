if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

# Regression guard: on small panels the validation split holds only a handful
# of lines, its loss rarely beats the first epochs, and restoring that "best"
# epoch returned a model predicting about the mean for every line (TabNet,
# TabAttention, GPNet). Best-epoch selection now waits for a warm-up and keeps
# the final epoch when the selected one is near-constant.

test_that("best-epoch selection has a warm-up and rejects near-constant epochs", {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  testthat::skip_if(is.null(py) || !file.exists(py), "No DL Python runtime found")
  py_dir <- system.file("python", package = "PredictProR")
  testthat::skip_if(!file.exists(file.path(py_dir, "dl_models.py")), "dl_models.py not installed")

  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, numpy as np, torch",
    sprintf("sys.path.insert(0, r'%s')", normalizePath(py_dir, winslash = "/")),
    "import dl_models as d",
    "print('WARMUP', d._best_epoch_warmup(5), d._best_epoch_warmup(30), d._best_epoch_warmup(100), d._best_epoch_warmup(400))",
    "y = np.linspace(-1, 1, 50)",
    "print('CONST', d._near_constant(np.full(50, 0.3) + 1e-6 * y, y), d._near_constant(0.5 * y, y), d._near_constant(y, np.zeros(50)))",
    "rng = np.random.default_rng(3)",
    "X = rng.integers(0, 3, size=(70, 120)).astype('float32')",
    "b = np.zeros(120); b[:10] = rng.normal(size=10)",
    "yy = (X @ b + rng.normal(scale=0.5, size=70)).astype('float32'); yy = (yy - yy.mean()) / yy.std()",
    "for mt in ['tabnet', 'gp_dkl']:",
    "    torch.manual_seed(1)",
    "    m, h = d.fit_model(X, yy, model_type=mt, task='regression', epochs=60, batch_size=16, device='cpu', random_seed=7, validation_split=0.2)",
    "    p = m.predict_mean(X) if hasattr(m, 'predict_mean') else d._regression_train_predictions(m, X, torch.device('cpu'))",
    "    print('SPREAD', mt, 'ok' if not d._near_constant(p, yy) else 'constant')"
  ), script)

  out <- suppressWarnings(system2(py, script, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(out, collapse = "\n"))
  expect_true("WARMUP 5 10 25 100" %in% out, info = paste(out, collapse = "\n"))
  expect_true("CONST True False False" %in% out, info = paste(out, collapse = "\n"))
  expect_true("SPREAD tabnet ok" %in% out, info = paste(out, collapse = "\n"))
  if (any(grepl("^SPREAD gp_dkl", out))) {
    expect_true("SPREAD gp_dkl ok" %in% out, info = paste(out, collapse = "\n"))
  }
})
