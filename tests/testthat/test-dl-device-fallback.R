# Deep-learning device fallbacks. A fit moved to the CPU without a word when
# device = "cuda" found no CUDA GPU, when moving the model to the GPU failed,
# or on any error during training or prediction on the GPU; nothing recorded
# which device a fit used. The Python trainer now reports every such switch
# (one "PREDICTPRO_DL_DEVICE {json}" line per bridge command), R warns, and
# model_execute() records the devices in Run_metadata.

device_report_line <- function(fits = list(cuda = 2L), notes = list()) {
  paste0(PredictProR:::gp_dl_device_report_prefix,
         jsonlite::toJSON(list(fits = fits, notes = notes), auto_unbox = TRUE))
}

oom_note <- list(kind = "training", model = "ft_transformer", count = 2L,
                 message = "training on cuda failed (OutOfMemoryError: CUDA out of memory); the fit continued on the CPU.")

test_that("device is NULL (automatic) or a device PyTorch can use", {
  arg <- PredictProR:::gp_dl_device_arg
  expect_null(arg(NULL))
  expect_null(arg(NA))
  expect_null(arg(" Auto "))
  expect_null(arg(""))
  expect_identical(arg("CUDA"), "cuda")
  expect_identical(arg("gpu"), "cuda")
  expect_identical(arg("cuda:1"), "cuda:1")
  expect_identical(arg("mps"), "mps")
  expect_identical(arg("cpu"), "cpu")
  for (bad in list("tpu", "cuda:x", c("cpu", "cuda"), 1)) {
    expect_error(arg(bad), "`device` must be NULL")
  }
  expect_null(PredictProR:::normalize_common(list(device = "auto"))$device)
  # model_execute() checks it before reading any data
  expect_error(PredictProR::model_execute(device = "tpu"), "`device` must be NULL")
})

test_that("a device fallback reported by the bridge is a warning outside a run", {
  withr::local_envvar(c(PREDICTPRO_DL_DEVICE_LOG = NA))
  out <- c("PredictProR deep learning: ft_transformer: training on cuda failed", "noise",
           device_report_line(list(cuda = 1L, cpu = 2L), list(oom_note)))
  expect_warning(
    PredictProR:::gp_dl_device_report_record(out),
    "TabTransformer, 2 fits: training on cuda failed.*Devices used: cuda: 1 fit; cpu: 2 fits"
  )
  # fits on the expected device alone are not worth a warning
  expect_no_warning(PredictProR:::gp_dl_device_report_record(device_report_line()))
  expect_no_warning(PredictProR:::gp_dl_device_report_record(c("no report", "PREDICTPRO_DL_DEVICE {broken")))
})

test_that("gp_dl_run_cli passes the bridge's device report on", {
  withr::local_envvar(c(PREDICTPRO_DL_DEVICE_LOG = NA))
  fake_bridge <- function(...) c("training ...", device_report_line(list(cpu = 2L), list(oom_note)))
  expect_warning(
    PredictProR:::gp_dl_run_cli("fit-predict", character(), python_bin = "python", run_fn = fake_bridge),
    "TabTransformer, 2 fits"
  )
})

test_that("inside a run, every process's reports reach Run_metadata and one warning at the end", {
  withr::local_envvar(c(PREDICTPRO_DL_DEVICE_LOG = NA))
  log <- PredictProR:::gp_dl_device_log_start()
  withr::defer(unlink(log$dir, recursive = TRUE, force = TRUE))
  expect_true(log$owner)
  expect_true(dir.exists(log$dir))
  expect_identical(Sys.getenv("PREDICTPRO_DL_DEVICE_LOG"), log$dir)
  # a model_execute() call inside another does not take the log over
  expect_false(PredictProR:::gp_dl_device_log_start()$owner)
  # gp_parallel_runtime_env() hands the folder to the workers
  expect_identical(PredictProR:::gp_parallel_runtime_env()[["PREDICTPRO_DL_DEVICE_LOG"]], log$dir)

  # this session, and a worker (another process writes its own file)
  expect_no_warning(PredictProR:::gp_dl_device_report_record(device_report_line(list(cuda = 3L))))
  writeLines(sub("^PREDICTPRO_DL_DEVICE ", "", device_report_line(list(cpu = 2L), list(oom_note))),
             file.path(log$dir, "dl-device-99999.jsonl"))

  fields <- PredictProR:::gp_dl_device_metadata_fields()
  expect_identical(fields$dl_device, "cuda: 3 fits; cpu: 2 fits")
  expect_match(fields$dl_device_notes, "^TabTransformer, 2 fits: training on cuda failed")
  local_mocked_bindings(gp_python_runtime_status = function(...) list(), .package = "PredictProR")
  meta <- PredictProR:::gp_runtime_metadata()
  expect_identical(meta$value[meta$key == "dl_device"], "cuda: 3 fits; cpu: 2 fits")
  expect_match(meta$value[meta$key == "dl_device_notes"], "OutOfMemoryError")

  expect_warning(PredictProR:::gp_dl_device_log_finish(log),
                 "TabTransformer, 2 fits.*Recorded in Run_metadata")
  expect_false(dir.exists(log$dir))
  expect_identical(Sys.getenv("PREDICTPRO_DL_DEVICE_LOG"), "")
  # no deep-learning fits: no Run_metadata fields
  expect_identical(PredictProR:::gp_dl_device_metadata_fields(), list())
})

test_that("a run whose fits all used the expected device ends without a warning", {
  withr::local_envvar(c(PREDICTPRO_DL_DEVICE_LOG = NA))
  log <- PredictProR:::gp_dl_device_log_start()
  PredictProR:::gp_dl_device_report_record(device_report_line(list(cpu = 4L)))
  expect_identical(PredictProR:::gp_dl_device_metadata_fields()$dl_device_notes, "none")
  expect_no_warning(PredictProR:::gp_dl_device_log_finish(log))
})

run_dl_device_script <- function(py, body) {
  py_dir <- system.file("python", package = "PredictProR")
  script <- withr::local_tempfile(fileext = ".py", .local_envir = parent.frame())
  writeLines(c(
    "import importlib.util, json, sys",
    "import numpy as np, torch, torch.nn as nn",
    sprintf("spec = importlib.util.spec_from_file_location('dlm', r'%s')", file.path(py_dir, "dl_models.py")),
    "d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)",
    "rng = np.random.default_rng(1)",
    "X = rng.normal(size=(60, 30)).astype('float32'); y = X[:, 0] + 0.1 * rng.normal(size=60)",
    "fast = dict(epochs=3, batch_size=16, compile_model=False, task='regression')",
    "res = {}",
    body,
    "print('RESULT ' + json.dumps(res))"
  ), script)
  err <- withr::local_tempfile(fileext = ".txt", .local_envir = parent.frame())
  out <- suppressWarnings(system2(py, script, stdout = TRUE, stderr = err))
  line <- grep("^RESULT ", out, value = TRUE)
  if (!length(line)) {
    stop(paste(c(out, readLines(err, warn = FALSE)), collapse = "\n"))
  }
  list(result = jsonlite::fromJSON(sub("^RESULT ", "", line[[1L]]), simplifyVector = FALSE),
       stderr = readLines(err, warn = FALSE))
}

test_that("device = \"cuda\" without a CUDA GPU trains on the CPU and says so", {
  py <- skip_if_no_dl_python()
  run <- run_dl_device_script(py, c(
    "torch.cuda.is_available = lambda: False",
    "m, _ = d.fit_model(X, y, model_type='mlp', device='cuda', **fast)",
    "res['trained'] = m._trained_device",
    "res['report'] = d.device_report()",
    "try:",
    "    d.fit_model(X, y, model_type='mlp', device='tpu', **fast)",
    "except ValueError as e:",
    "    res['bad_device'] = str(e)"
  ))
  res <- run$result
  expect_identical(res$trained, "cpu")
  expect_identical(res$report$fits$cpu, 1L)
  expect_length(res$report$notes, 1L)
  expect_identical(res$report$notes[[1L]]$kind, "cuda_unavailable")
  expect_match(res$report$notes[[1L]]$message, "device = \"cuda\" was given, but PyTorch .* sees no CUDA GPU")
  expect_match(res$bad_device, "device must be")
  expect_true(any(grepl("PredictProR deep learning: mlp: device = \"cuda\" was given", run$stderr)))
})

test_that("an NVIDIA GPU that a CPU-only PyTorch cannot use is reported, unless the GPUs are hidden", {
  py <- skip_if_no_dl_python()
  run <- run_dl_device_script(py, c(
    "import os",
    "os.environ.pop('CUDA_VISIBLE_DEVICES', None); os.environ.pop('PREDICTPRO_TORCH_DEVICE', None)",
    "torch.cuda.is_available = lambda: False",
    "d._mps_available = lambda: False",
    "torch.version.cuda = None",
    "d._deps._driver_cuda_version = lambda: (12, 4)",
    "d.fit_model(X, y, model_type='mlp', **fast)",
    "d.fit_model(X, y, model_type='mlp', **fast)",
    "res['gpu'] = d.device_report()",
    "d._GPU_UNUSED = None",
    "os.environ['CUDA_VISIBLE_DEVICES'] = '-1'",
    "d.fit_model(X, y, model_type='mlp', **fast)",
    "res['hidden'] = d.device_report()",
    "d._GPU_UNUSED = None; del os.environ['CUDA_VISIBLE_DEVICES']",
    "d._deps._driver_cuda_version = lambda: None",
    "d.fit_model(X, y, model_type='mlp', **fast)",
    "res['no_gpu'] = d.device_report()",
    "d._GPU_UNUSED = None; d._deps._driver_cuda_version = lambda: (12, 4)",
    "d.fit_model(X, y, model_type='mlp', device='cpu', **fast)",
    "res['cpu_asked'] = d.device_report()"
  ))
  res <- run$result
  expect_identical(res$gpu$fits$cpu, 2L)
  expect_identical(res$gpu$notes[[1L]]$kind, "gpu_unused")
  expect_identical(res$gpu$notes[[1L]]$count, 2L)
  expect_match(res$gpu$notes[[1L]]$message, "is a CPU-only build")
  # checked with nvidia-smi once per process, printed once
  expect_identical(sum(grepl("CPU-only build", run$stderr)), 1L)
  for (case in c("hidden", "no_gpu", "cpu_asked")) {
    expect_identical(res[[case]]$fits$cpu, 1L, info = case)
    expect_length(res[[case]]$notes, 0L)
  }
})

test_that("a CUDA error moves the fit to the CPU, a failing compiled model trains uncompiled, and both are reported", {
  py <- skip_if_no_dl_python()
  cuda <- suppressWarnings(system2(py, c("-c", shQuote("import torch; print(torch.cuda.is_available())")),
                                   stdout = TRUE, stderr = FALSE))
  skip_if(!identical(utils::tail(cuda, 1L), "True"), "no CUDA GPU for this PyTorch")
  run <- run_dl_device_script(py, c(
    # the 3rd training step on the GPU runs out of memory
    "real_build = d.build_model",
    "def failing_build(*a, **k):",
    "    net = real_build(*a, **k); fwd = net.forward; calls = [0]",
    "    def forward(x):",
    "        if x.is_cuda and net.training:",
    "            calls[0] += 1",
    "            if calls[0] == 3:",
    "                raise torch.cuda.OutOfMemoryError('CUDA out of memory (test)')",
    "        return fwd(x)",
    "    net.forward = forward",
    "    return net",
    "d.build_model = failing_build",
    "out = d.fit_predict_payload('mlp', X, y, X[:5], params=dict(epochs=12, batch_size=16, compile_model=False))",
    "res['oom'] = {'n': len(out['predictions']), 'finite': bool(np.isfinite(out['predictions']).all()),",
    "              'epochs': len(out['history']['train_loss']), 'report': d.device_report()}",
    "d.build_model = real_build",
    # torch.compile that fails: the same GPU, uncompiled
    "class BrokenCompile(nn.Module):",
    "    def __init__(self, inner):",
    "        super().__init__(); self._orig_mod = inner",
    "    def forward(self, x):",
    "        raise RuntimeError('BackendCompilerFailed (test)')",
    "real_compile = d._maybe_compile",
    "d._maybe_compile = lambda model, dev, compile_model=True: BrokenCompile(model)",
    "m, h = d.fit_model(X, y, model_type='mlp', epochs=12, batch_size=16, task='regression')",
    "d._maybe_compile = real_compile",
    "res['compile'] = {'trained': m._trained_device, 'cls': type(m).__name__, 'epochs': len(h['train_loss']),",
    "                  'report': d.device_report()}",
    # prediction on the GPU fails: predicted on the CPU
    "fwd = m.forward",
    "def pfwd(x):",
    "    if x.is_cuda and not m.training:",
    "        raise RuntimeError('CUDA error: an illegal memory access was encountered (test)')",
    "    return fwd(x)",
    "m.forward = pfwd",
    "p = d.predict(m, X[:5], device=m._trained_device)",
    "res['predict'] = {'n': len(p), 'report': d.device_report()}",
    # no error: one cuda fit, no notes
    "d.fit_model(X, y, model_type='mlp', **fast)",
    "res['normal'] = d.device_report()"
  ))
  res <- run$result
  expect_identical(res$oom$n, 5L)
  expect_true(res$oom$finite)
  expect_identical(res$oom$epochs, 12L)
  expect_identical(names(res$oom$report$fits), "cpu")
  expect_identical(res$oom$report$notes[[1L]]$kind, "training")
  expect_match(res$oom$report$notes[[1L]]$message, "training on cuda.* failed \\(OutOfMemoryError: CUDA out of memory \\(test\\)\\); the fit continued on the CPU")

  expect_match(res$compile$trained, "^cuda")
  expect_identical(res$compile$cls, "MLPNet")
  expect_identical(res$compile$epochs, 12L)
  expect_identical(res$compile$report$notes[[1L]]$kind, "compile")

  expect_identical(res$predict$n, 5L)
  expect_identical(res$predict$report$notes[[1L]]$kind, "prediction")

  expect_identical(names(res$normal$fits), "cuda")
  expect_length(res$normal$notes, 0L)
})

test_that("the bridge's device report reaches R as a warning", {
  py <- skip_if_no_dl_python()
  skip_if_python_lacks(py, c("numpy", "pandas", "torch", "sklearn"), "deep-learning")
  withr::local_envvar(c(PREDICTPRO_DL_DEVICE_LOG = NA, CUDA_VISIBLE_DEVICES = "-1"))
  set.seed(1)
  X <- matrix(rnorm(40 * 12), 40, 12)
  y <- X[, 1] + rnorm(40, sd = 0.1)
  expect_warning(
    res <- PredictProR:::gp_dl_bridge_fit_predict_raw(
      "mlp", X[1:30, ], y[1:30], X[31:40, ],
      dl_args = list(device = "cuda", epochs = 2L, compile_model = FALSE)
    ),
    "DenseNeuralNet, 1 fit: device = \"cuda\" was given, but PyTorch .* sees no CUDA GPU.*Devices used: cpu: 1 fit"
  )
  expect_length(res$predictions, 10L)
})
