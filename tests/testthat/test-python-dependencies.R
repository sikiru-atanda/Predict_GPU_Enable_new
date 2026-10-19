# PredictProR installs the Python packages its models need on first use
# (inst/python/predictpror_deps.py, driven from R by R/python_dependencies.R).
# The R tests mock the installer; the Python tests install a locally built
# wheel into a fresh virtualenv (no network), race three processes for the
# environment lock, and check failure reporting and the PyTorch build choice.

test_that("requested models map to the Python packages they need", {
  plan <- PredictProR:::gp_python_dependency_groups(c("Xgboost", "CatBoost", "tabnet", "gp_dkl", "GP", "GBLUP"))
  expect_setequal(plan$ml, c("tabular", "sklearn", "xgboost", "catboost"))
  expect_setequal(plan$dl, c("tabular", "torch", "sklearn", "gpytorch"))
  expect_setequal(plan$gp, c("tabular", "torch", "gp"))
  # Lasso is fitted in R (glmnet); Bayesian/ASReml models need no Python
  expect_length(PredictProR:::gp_python_dependency_groups(c("Lasso", "BayesA", "GBLUP")), 0L)
  # hybrid_gp with "auto" falls back to the R solver, so it installs nothing
  expect_null(PredictProR:::gp_python_dependency_groups("GP", gp_backend = "auto", hybrid_gp = TRUE)$gp)
  expect_setequal(PredictProR:::gp_python_dependency_groups("GP", gp_backend = "python", hybrid_gp = TRUE)$gp,
                  c("tabular", "torch", "gp"))
})

test_that("verified groups (or a superset) skip the check, also in parallel workers", {
  withr::local_envvar(c(PREDICTPRO_PYTHON_DEPS_OK = NA))
  py <- tempfile("python-", fileext = ".exe")
  file.create(py)
  expect_false(PredictProR:::gp_python_deps_verified(py, c("tabular", "torch")))
  PredictProR:::gp_python_mark_verified(py, c("tabular", "torch", "sklearn", "gpytorch"))
  expect_true(PredictProR:::gp_python_deps_verified(py, c("torch", "tabular")))
  expect_false(PredictProR:::gp_python_deps_verified(py, c("tabular", "gp")))
  # workers receive the list with the rest of the runtime environment
  expect_true("PREDICTPRO_PYTHON_DEPS_OK" %in% names(PredictProR:::gp_parallel_runtime_env()))
})

test_that("missing packages are installed once, then the interpreter counts as verified", {
  withr::local_envvar(c(PREDICTPRO_PYTHON_DEPS_OK = NA, PREDICTPRO_AUTO_INSTALL = NA))
  py <- normalizePath(tempfile("python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  file.create(py)
  calls <- character()
  local_mocked_bindings(
    gp_python_deps_run = function(python_bin, command, groups, stream = FALSE, extra_args = character(), timeout = 0) {
      calls <<- c(calls, command)
      if (identical(command, "check")) {
        list(status = "ok", missing = c("numpy>=1.24", "scikit-learn"), install_problem = NULL)
      } else {
        list(status = "ok", installed = c("numpy>=1.24", "scikit-learn"))
      }
    },
    .package = "PredictProR"
  )
  expect_message(
    out <- PredictProR:::gp_python_ensure_dependencies("ml", c("tabular", "sklearn"), python_bin = py),
    "installing Python packages for machine-learning models"
  )
  expect_identical(out, normalizePath(py, winslash = "/", mustWork = TRUE))
  expect_identical(calls, c("check", "ensure"))
  PredictProR:::gp_python_ensure_dependencies("ml", c("tabular", "sklearn"), python_bin = py)
  expect_identical(calls, c("check", "ensure"))   # verified: no second check
})

test_that("Run_metadata names the interpreter each Python model family of the run used", {
  ml_py <- normalizePath(tempfile("ml-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  local_mocked_bindings(
    gp_python_ensure_dependencies = function(purpose, groups, python_bin = NULL, quiet = FALSE) ml_py,
    gp_python_runtime_status = function(...) list(),
    .package = "PredictProR"
  )
  withr::defer(PredictProR:::gp_python_prepare_runtimes("BayesA"))
  PredictProR:::gp_python_prepare_runtimes(c("RandomForest", "BayesA"))
  meta <- PredictProR:::gp_runtime_metadata(include_accelerator = FALSE)
  expect_identical(meta$value[meta$key == "python_ml"], ml_py)
  expect_false(any(c("python_dl", "python_gp") %in% meta$key))
  # the next run starts from an empty record
  PredictProR:::gp_python_prepare_runtimes("BayesA")
  meta <- PredictProR:::gp_runtime_metadata(include_accelerator = FALSE)
  expect_false("python_ml" %in% meta$key)
})

test_that("an interpreter that cannot take installs is replaced by a PredictProR virtualenv", {
  withr::local_envvar(c(PREDICTPRO_PYTHON_DEPS_OK = NA, PREDICTPRO_ML_PYTHON = NA, PREDICTPRO_AUTO_INSTALL = NA))
  base <- normalizePath(tempfile("system-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  venv <- normalizePath(tempfile("venv-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  file.create(c(base, venv))
  seen <- character()
  local_mocked_bindings(
    gp_python_deps_run = function(python_bin, command, groups, stream = FALSE, extra_args = character(), timeout = 0) {
      seen <<- c(seen, paste(basename(python_bin), command))
      if (identical(python_bin, base)) {
        return(list(status = "ok", missing = "scikit-learn", install_problem = "externally_managed"))
      }
      if (identical(command, "check")) list(status = "ok", missing = "scikit-learn", install_problem = NULL)
      else list(status = "ok")
    },
    gp_python_managed_env = function(base_python, purpose, reason) {
      expect_identical(reason, "externally_managed")
      venv
    },
    .package = "PredictProR"
  )
  out <- suppressMessages(PredictProR:::gp_python_ensure_dependencies("ml", c("tabular", "sklearn"), python_bin = base))
  expect_identical(out, normalizePath(venv, winslash = "/", mustWork = TRUE))
  expect_identical(Sys.getenv("PREDICTPRO_ML_PYTHON"), venv)   # the rest of the session uses it
  expect_identical(seen, c(paste(basename(base), "check"), paste(basename(venv), "check"),
                           paste(basename(venv), "ensure")))
})

test_that("a PredictProR environment that cannot take installs either stops, it is not replaced again", {
  withr::local_envvar(c(PREDICTPRO_PYTHON_DEPS_OK = NA, PREDICTPRO_ML_PYTHON = NA, PREDICTPRO_AUTO_INSTALL = NA))
  base <- normalizePath(tempfile("system-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  venv <- normalizePath(tempfile("venv-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  file.create(c(base, venv))
  made <- 0L
  local_mocked_bindings(
    gp_python_deps_run = function(python_bin, command, groups, stream = FALSE, extra_args = character(), timeout = 0) {
      list(status = "ok", missing = "scikit-learn", install_problem = "not_writable")
    },
    gp_python_managed_env = function(base_python, purpose, reason) {
      made <<- made + 1L
      venv
    },
    .package = "PredictProR"
  )
  expect_error(
    suppressMessages(PredictProR:::gp_python_ensure_dependencies("ml", c("tabular", "sklearn"), python_bin = base)),
    "PredictProR Python environment .* cannot be written to"
  )
  expect_identical(made, 1L)
})

test_that("interpreter paths keep a virtualenv's python link instead of resolving it to the base Python", {
  skip_on_os("windows")
  root <- withr::local_tempdir()
  base <- file.path(root, "base-python3")
  file.create(base)
  dir.create(file.path(root, "venv", "bin"), recursive = TRUE)
  link <- file.path(root, "venv", "bin", "python")
  skip_if_not(file.symlink(base, link), "cannot create symlinks here")
  kept <- file.path(normalizePath(file.path(root, "venv", "bin")), "python")
  expect_identical(PredictProR:::gp_python_path(link, mustWork = TRUE), kept)
  expect_identical(normalizePath(link), normalizePath(base))   # what normalizePath() would have used
  withr::local_envvar(c(PREDICTPRO_ML_PYTHON = link, PREDICTPRO_DL_PYTHON = link, PREDICTPRO_GP_PYTHON = link))
  expect_identical(PredictProR:::gp_detect_ml_python_uncached(), kept)
  expect_identical(PredictProR:::gp_detect_dl_python(), kept)
  expect_identical(PredictProR:::gp_detect_gp_python(), kept)
})

test_that("a venv module without ensurepip (Debian/Ubuntu without python3-venv) makes the environment without pip", {
  skip_on_os("windows")
  root <- withr::local_tempdir()
  fake <- file.path(root, "python3")
  writeLines(c(
    "#!/bin/sh",
    "# stands in for a Debian python3 whose venv module lacks ensurepip",
    "if [ \"$1\" = \"--version\" ]; then echo 'Python 3.10.12'; exit 0; fi",
    "for last; do :; done",
    "mkdir -p \"$last/bin\"",
    "printf '#!/bin/sh\\necho Python 3.10.12\\n' > \"$last/bin/python\"",
    "chmod +x \"$last/bin/python\"",
    "case \" $* \" in *' --without-pip '*) exit 0 ;; esac",
    "echo 'The virtual environment was not created successfully because ensurepip is not available.' >&2",
    "exit 1"
  ), fake)
  Sys.chmod(fake, "0755")
  env_dir <- file.path(root, "env")
  made <- PredictProR:::gp_python_make_venv(fake, env_dir)
  expect_true(made$ok, info = paste(made$output, collapse = "\n"))
  expect_true(any(grepl("ensurepip", made$output)))
  expect_true(file.exists(file.path(env_dir, "bin", "python")))

  # any other venv failure is reported, and the half-made directory removed
  writeLines(c("#!/bin/sh", "if [ \"$1\" = \"--version\" ]; then exit 0; fi",
               "for last; do :; done", "mkdir -p \"$last\"", "echo 'disk full' >&2", "exit 1"), fake)
  failed <- PredictProR:::gp_python_make_venv(fake, env_dir)
  expect_false(failed$ok)
  expect_true(any(grepl("disk full", failed$output)))
  expect_false(dir.exists(env_dir))
})

test_that("PREDICTPRO_AUTO_INSTALL = false stops with the packages that are missing", {
  withr::local_envvar(c(PREDICTPRO_PYTHON_DEPS_OK = NA, PREDICTPRO_AUTO_INSTALL = "false"))
  py <- tempfile("python-", fileext = ".exe")
  file.create(py)
  local_mocked_bindings(
    gp_python_deps_run = function(python_bin, command, groups, stream = FALSE, extra_args = character(), timeout = 0) {
      if (!identical(command, "check")) stop("must not install")
      list(status = "ok", missing = c("torch", "pandas>=1.5"), install_problem = NULL)
    },
    .package = "PredictProR"
  )
  expect_error(
    PredictProR:::gp_python_ensure_dependencies("dl", c("tabular", "torch"), python_bin = py),
    "missing from .*torch, pandas>=1.5.*PREDICTPRO_AUTO_INSTALL is false"
  )
})

test_that("an explicitly configured GP interpreter takes precedence", {
  py <- normalizePath(tempfile("gp-python-", fileext = ".exe"), winslash = "/", mustWork = FALSE)
  file.create(py)
  withr::local_envvar(c(PREDICTPRO_GP_PYTHON = py))
  local_mocked_bindings(gp_python_candidates_bridge = function() stop("must not score candidates"),
                        .package = "PredictProR")
  expect_identical(PredictProR:::gp_detect_gp_python(), normalizePath(py, winslash = "/", mustWork = TRUE))
})

test_that("true prediction stops with the causes when no model produced a result", {
  results <- list(
    list(GS_model = "Xgboost", res_model_output = NULL, model_error = "ML Python bridge failed: no numpy"),
    list(GS_model = "RandomForest", res_model_output = NULL, model_error = "")
  )
  expect_error(
    PredictProR:::gp_finalize_model_execute_results(results = results, GS_model = "Xgboost"),
    "No model produced a result.*Xgboost: ML Python bridge failed: no numpy.*RandomForest: no result"
  )
})

# ---- Python installer -------------------------------------------------------

deps_python <- function() {
  py <- tryCatch(PredictProR:::gp_detect_dl_python(), error = function(e) NULL)
  if (is.null(py) || !nzchar(py) || !file.exists(py)) {
    py <- tryCatch(PredictProR:::gp_detect_ml_python(), error = function(e) NULL)
  }
  skip_if(is.null(py) || !nzchar(py) || !file.exists(py), "no Python interpreter available")
  py
}

run_py <- function(py, code, env = character()) {
  script <- tempfile(fileext = ".py")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "import sys",
    sprintf("sys.path.insert(0, r'%s')", system.file("python", package = "PredictProR")),
    "import predictpror_deps as d",
    code
  ), script)
  out <- withr::with_envvar(env, suppressWarnings(system2(py, PredictProR:::gp_quote_system_args(script), stdout = TRUE, stderr = TRUE)))
  out
}

test_that("the PyTorch build follows the NVIDIA driver, with CPU as the fallback", {
  py <- deps_python()
  out <- run_py(py, c(
    "def tags(drv, **kw):",
    "    d._driver_cuda_version = lambda: drv",
    "    return ','.join(label for label, _ in d.torch_candidates(**kw))",
    "import platform; platform.system = lambda: 'Linux'",
    "print('A', tags((13, 2)))",
    "print('B', tags((12, 4)))",
    "print('C', tags(None))",
    "print('D', tags((13, 2), prefer_gpu=False))",
    "print('E', tags((13, 2), cuda='cu126'))",
    "import os; os.environ['PREDICTPRO_TORCH_DEVICE'] = 'cpu'; print('F', tags((13, 2)))",
    "os.environ.pop('PREDICTPRO_TORCH_DEVICE'); os.environ['PYTORCH_INDEX_URL'] = 'https://mirror/whl'",
    "print('G', d.torch_candidates())"
  ), env = c(PYTORCH_INDEX_URL = NA, PREDICTPRO_TORCH_DEVICE = NA))
  expect_true("A cu130,cu129,cu128,cpu" %in% out)
  expect_true("B cu124,cu121,cu118,cpu" %in% out)
  expect_true("C cpu" %in% out)
  expect_true("D cpu" %in% out)
  expect_true("E cu126,cpu" %in% out)
  expect_true("F cpu" %in% out)
  expect_true("G [('forced', 'https://mirror/whl')]" %in% out)
})

test_that("install problems are recognised before pip runs", {
  py <- deps_python()
  out <- run_py(py, c(
    "import os, tempfile, sysconfig",
    "print('ok', d.install_problem())",
    # PEP 668 is checked first; the system Python of e.g. Ubuntu 24.04 is externally managed
    "real = d._dir_writable, d._externally_managed",
    "d._dir_writable = lambda p: False; d._externally_managed = lambda: False",
    "print('ro', d.install_problem())",
    "d._dir_writable, d._externally_managed = real",
    "stdlib = tempfile.mkdtemp(); open(os.path.join(stdlib, 'EXTERNALLY-MANAGED'), 'w').close()",
    "sysconfig.get_path = lambda name, *a, **k: stdlib",
    "sys.base_prefix = sys.prefix",
    "print('pep668', d.install_problem())"
  ), env = c(PIP_BREAK_SYSTEM_PACKAGES = NA))
  expect_true("ro not_writable" %in% out)
  expect_true("pep668 externally_managed" %in% out)
})

test_that("a directory the user may not write to is reported at once, not retried", {
  skip_on_cran()
  py <- deps_python()
  dir <- withr::local_tempdir()
  if (.Platform$OS.type == "windows") {
    # an ACL deny: os.access() still says writable, which made tempfile loop
    user <- Sys.getenv("USERNAME")
    system2("icacls", c(shQuote(dir, type = "cmd"), "/deny", paste0(user, ":(OI)(CI)(WD,AD,WEA,WA)")),
            stdout = FALSE, stderr = FALSE)
    withr::defer(system2("icacls", c(shQuote(dir, type = "cmd"), "/remove:d", user, "/T"),
                         stdout = FALSE, stderr = FALSE))
  } else {
    Sys.chmod(dir, "0555")
    withr::defer(Sys.chmod(dir, "0755"))
    skip_if(file.create(file.path(dir, "probe")), "running as a user who can write anywhere")
  }
  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    sprintf("sys.path.insert(0, r'%s')", system.file("python", package = "PredictProR")),
    "import predictpror_deps as d",
    sprintf("print('WRITABLE', d._dir_writable(r'%s'))", dir)
  ), script)
  t0 <- Sys.time()
  out <- suppressWarnings(system2(py, PredictProR:::gp_quote_system_args(script), stdout = TRUE,
                                  stderr = TRUE, timeout = 60))
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 60)
  expect_true("WRITABLE False" %in% out, info = paste(out, collapse = "\n"))
})

test_that("missing packages install offline into a fresh virtualenv, once, under concurrent requests", {
  skip_on_cran()
  skip_if_not_installed("processx")
  base <- deps_python()
  root <- withr::local_tempdir()
  venv <- file.path(root, "venv")
  made <- suppressWarnings(system2(base, PredictProR:::gp_quote_system_args(c("-m", "venv", venv)), stdout = TRUE, stderr = TRUE))
  py <- PredictProR:::gp_python_venv_python(venv)
  skip_if(!file.exists(py), paste("could not create a virtualenv:", paste(made, collapse = " ")))

  # a local wheel stands in for a package index
  wheels <- file.path(root, "wheels")
  dir.create(wheels)
  build <- run_py(base, c(
    "import zipfile, hashlib, base64, os",
    sprintf("out = r'%s'", wheels),
    "files = {'predictpror_dummy/__init__.py': 'VALUE = 42\\n',",
    "  'predictpror_dummy-0.1.dist-info/METADATA': 'Metadata-Version: 2.1\\nName: predictpror-dummy\\nVersion: 0.1\\n',",
    "  'predictpror_dummy-0.1.dist-info/WHEEL': 'Wheel-Version: 1.0\\nGenerator: test\\nRoot-Is-Purelib: true\\nTag: py3-none-any\\n'}",
    "rec = []",
    "with zipfile.ZipFile(os.path.join(out, 'predictpror_dummy-0.1-py3-none-any.whl'), 'w') as z:",
    "    for name, text in files.items():",
    "        data = text.encode(); z.writestr(name, data)",
    "        h = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b'=').decode()",
    "        rec.append(f'{name},sha256={h},{len(data)}')",
    "    rec.append('predictpror_dummy-0.1.dist-info/RECORD,,')",
    "    z.writestr('predictpror_dummy-0.1.dist-info/RECORD', '\\n'.join(rec) + '\\n')",
    "print('built')"
  ))
  expect_true("built" %in% build)

  script <- file.path(root, "ensure_dummy.py")
  writeLines(c(
    "import sys",
    sprintf("sys.path.insert(0, r'%s')", system.file("python", package = "PredictProR")),
    "import predictpror_deps as d",
    "d.GROUPS['dummy'] = [('predictpror_dummy', 'predictpror-dummy')]",
    "r = d.ensure_groups(['dummy'], verbose=False)",
    "import predictpror_dummy",   # importable in the process that installed it
    "print('RESULT', len(r['installed']), predictpror_dummy.VALUE)"
  ), script)
  offline <- c(PIP_NO_INDEX = "1", PIP_FIND_LINKS = wheels, PREDICTPRO_AUTO_INSTALL = NA)
  procs <- withr::with_envvar(offline, lapply(1:3, function(i) {
    processx::process$new(py, script, stdout = "|", stderr = "|")
  }))
  for (p in procs) p$wait(timeout = 180000)
  outs <- vapply(procs, function(p) paste(p$read_all_output_lines(), collapse = "\n"), character(1))
  errs <- vapply(procs, function(p) paste(p$read_all_error_lines(), collapse = "\n"), character(1))
  expect_true(all(vapply(procs, function(p) p$get_exit_status(), integer(1)) == 0L), info = paste(errs, collapse = "\n---\n"))
  expect_true(all(grepl("RESULT [01] 42", outs)), info = paste(outs, collapse = " | "))
  # the lock serialises the processes: one installs, the others find it installed
  expect_identical(sum(grepl("RESULT 1 42", outs)), 1L, info = paste(outs, collapse = " | "))

  # an unreachable index is reported as a network problem, with pip's output
  net <- run_py(py, c(
    "d.GROUPS['other'] = [('predictpror_other', 'predictpror-other')]",
    "try:",
    "    d.ensure_groups(['other'], verbose=False)",
    "except d.DependencyError as e:",
    "    print('REASON', e.reason); print(str(e))"
  ), env = c(PIP_NO_INDEX = NA, PIP_FIND_LINKS = NA, PIP_INDEX_URL = "http://127.0.0.1:9/simple",
             PIP_RETRIES = "0", PIP_TIMEOUT = "3", PREDICTPRO_AUTO_INSTALL = NA))
  expect_true("REASON network" %in% net, info = paste(net, collapse = "\n"))
  expect_true(any(grepl("PIP_INDEX_URL", net)))

  # no pip and no ensurepip: the installer tries pip's zipapp and says where it looked
  nopip <- run_py(py, c(
    "d._pip_works = lambda: False",
    "try:",
    "    d._ensure_pip()",
    "except d.DependencyError as e:",
    "    print('REASON', e.reason); print(str(e))"
  ), env = c(PREDICTPRO_PIP_BOOTSTRAP_URL = "http://127.0.0.1:9/pip.pyz"))
  expect_true("REASON pip_missing" %in% nopip, info = paste(nopip, collapse = "\n"))
  expect_true(any(grepl("127.0.0.1:9/pip.pyz", nopip, fixed = TRUE)))

  # installing turned off
  off <- run_py(py, c(
    "d.GROUPS['other'] = [('predictpror_other', 'predictpror-other')]",
    "try:",
    "    d.ensure_groups(['other'])",
    "except d.DependencyError as e:",
    "    print('REASON', e.reason)"
  ), env = c(PREDICTPRO_AUTO_INSTALL = "false"))
  expect_true("REASON auto_install_disabled" %in% off)
})
