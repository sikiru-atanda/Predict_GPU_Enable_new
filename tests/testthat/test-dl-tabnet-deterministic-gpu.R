if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

# Regression guard: with deterministic = TRUE on a GPU, TabNet's sparsemax used
# torch.cumsum (no deterministic CUDA kernel), every fit fell back to the CPU,
# and the fallback left the class-weighted loss on CUDA, so binary and
# multiclass TabNet failed ("Expected all tensors to be on the same device").

test_that("TabNet trains on a deterministic GPU and the CPU fallback moves the loss", {
  py <- PredictProR:::gp_preferred_python(purpose = "dl")
  testthat::skip_if(is.null(py) || !file.exists(py), "No DL Python runtime found")
  py_dir <- system.file("python", package = "PredictProR")
  testthat::skip_if(!file.exists(file.path(py_dir, "dl_models.py")), "dl_models.py not installed")

  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, numpy as np, torch",
    sprintf("sys.path.insert(0, r'%s')", normalizePath(py_dir, winslash = "/")),
    "import dl_models",
    "if not torch.cuda.is_available():",
    "    print('NO_CUDA'); sys.exit(0)",
    "fallbacks = []",
    "orig = dl_models._safe_to_device",
    "def spy(model, dev):",
    "    if dev.type == 'cpu': fallbacks.append(1)",
    "    return orig(model, dev)",
    "dl_models._safe_to_device = spy",
    "torch.use_deterministic_algorithms(True)",
    "rng = np.random.default_rng(1); X = rng.normal(size=(48, 40)).astype('float32')",
    "ys = {'regression': rng.normal(size=48).astype('float32'),",
    "      'binary': (rng.random(48) > 0.5).astype('float32'),",
    "      'multiclass': rng.integers(0, 3, 48).astype('int64')}",
    "for task, y in ys.items():",
    "    fallbacks.clear()",
    "    dl_models.fit_model(X, y, model_type='tabnet', task=task, epochs=1, batch_size=16, device='cuda', auto_class_weights=True)",
    "    print('FIT', task, 'fallback' if fallbacks else 'gpu')",
    "real = dl_models.TabNet.forward",
    "for task in ['binary', 'multiclass']:",
    "    calls = {'n': 0}",
    "    def bad(self, x, _r=real):",
    "        calls['n'] += 1",
    "        if calls['n'] == 1 and x.is_cuda: raise RuntimeError('forced')",
    "        return _r(self, x)",
    "    dl_models.TabNet.forward = bad",
    "    dl_models.fit_model(X, ys[task], model_type='tabnet', task=task, epochs=1, batch_size=16, device='cuda', auto_class_weights=True)",
    "    dl_models.TabNet.forward = real",
    "    print('FALLBACK', task, 'ok')"
  ), script)

  out <- withr::with_envvar(
    c(CUBLAS_WORKSPACE_CONFIG = ":4096:8"),
    suppressWarnings(system2(py, script, stdout = TRUE, stderr = TRUE))
  )
  testthat::skip_if(any(out == "NO_CUDA"), "No CUDA device")
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(out, collapse = "\n"))
  for (task in c("regression", "binary", "multiclass")) {
    expect_true(paste("FIT", task, "gpu") %in% out, info = paste(out, collapse = "\n"))
  }
  expect_true(all(c("FALLBACK binary ok", "FALLBACK multiclass ok") %in% out),
              info = paste(out, collapse = "\n"))
})
