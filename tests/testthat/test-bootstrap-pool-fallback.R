# If the Python bootstrap process pool cannot start (Windows "OSError: [Errno
# 22]" under memory pressure lost a whole mice RandomForest prediction), the
# refits run in-process from the same pre-drawn indices: identical results.
test_that("bootstrap falls back to in-process refits with identical results when the pool fails", {
  py <- skip_if_no_ml_python()
  bridge_dir <- dirname(PredictProR:::gp_ml_bridge_script_path())
  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, numpy as np, concurrent.futures as cf",
    sprintf("sys.path.insert(0, r'%s')", normalizePath(bridge_dir, winslash = "/")),
    "import ml_bridge",
    "rng = np.random.default_rng(1)",
    "X = rng.normal(size=(80, 20)); y = X[:, 0] + rng.normal(size=80); Xp = rng.normal(size=(15, 20))",
    "args = dict(model='ridge', task='gaussian', X_train=X, y_train=y, X_pred=Xp, n_bootstrap=6, seed=3, params={})",
    "serial = ml_bridge.bootstrap_fit_predict_payload(n_jobs=1, **args)['bootstrap']",
    "class FailingPool:",
    "    def __init__(self, *a, **k): raise OSError(22, 'Invalid argument')",
    "cf.ProcessPoolExecutor = FailingPool",
    "fallback = ml_bridge.bootstrap_fit_predict_payload(n_jobs=4, **args)['bootstrap']",
    "print('IDENTICAL' if np.array_equal(serial, fallback) and fallback.shape == (6, 15) else 'DIFFERENT')"
  ), script)
  out <- suppressWarnings(system2(py, shQuote(script), stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("IDENTICAL", out)), info = paste(out, collapse = "\n"))
  expect_true(any(grepl("running the refits in-process", out)))
})
