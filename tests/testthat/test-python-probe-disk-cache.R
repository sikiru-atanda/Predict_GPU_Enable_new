# Checking an interpreter starts it and imports torch (0.2-3 s each). An
# interactive library(PredictProR) checked 12 interpreters three times and
# listed conda environments three times: 23 s per fresh session on a Windows
# test machine. Checks are kept on disk until the interpreter or its
# site-packages change, so the next session starts no Python.

# A fake interpreter whose check is answered by a mock (counted), and a
# private cache file instead of the user's.
local_probe_cache <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  pc <- new.env()
  pc$dir <- normalizePath(dir, winslash = "/")
  pc$file <- file.path(pc$dir, "python-discovery.rds")
  pc$site <- file.path(pc$dir, "site-packages")
  dir.create(pc$site)
  pc$py <- file.path(pc$dir, "python.exe")
  writeLines("fake interpreter", pc$py)
  pc$runs <- 0L
  pc$site_dirs <- pc$site
  pc$answer <- function() {
    sprintf('{"numpy": true, "torch": false, "_site_dirs": [%s]}',
            paste0('"', pc$site_dirs, '"', collapse = ", "))
  }
  local_mocked_bindings(
    gp_python_discovery_cache_file = function() pc$file,
    gp_python_probe_run = function(py_bin, py_code) {
      pc$runs <- pc$runs + 1L
      pc$answer()
    },
    .package = "PredictProR",
    .env = env
  )
  withr::local_envvar(PREDICTPRO_PROBE_CACHE = "true", .local_envir = env)
  cache <- PredictProR:::PredictProR_runtime_cache
  saved <- mget(ls(cache), envir = cache)
  withr::defer({
    rm(list = ls(cache), envir = cache)
    list2env(saved, envir = cache)
  }, envir = env)
  new_session()
  pc
}

# What a new R session starts with: nothing checked yet in memory.
new_session <- function() {
  cache <- PredictProR:::PredictProR_runtime_cache
  keys <- grep("^(python_probe|preferred_python|ml_python|conda_candidates)::", ls(cache), value = TRUE)
  rm(list = keys, envir = cache)
}

probe <- function(pc) PredictProR:::gp_python_probe(pc$py, modules = c("numpy", "torch"))

test_that("an interpreter checked in one session is not started again in the next", {
  pc <- local_probe_cache()
  first <- probe(pc)
  expect_identical(pc$runs, 1L)
  expect_true(isTRUE(first$numpy))
  expect_false(isTRUE(first$torch))
  expect_null(first[["_site_dirs"]])
  expect_true(file.exists(pc$file))

  new_session()
  again <- probe(pc)
  expect_identical(pc$runs, 1L)
  expect_identical(again, first)
})

test_that("installing a package into the interpreter makes the next session check it again", {
  pc <- local_probe_cache()
  probe(pc)
  # pip and conda add or rename a dist-info folder, which changes the
  # modification time of site-packages
  Sys.setFileTime(pc$site, Sys.time() + 5)
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("a user site-packages folder created later makes the next session check again", {
  pc <- local_probe_cache()
  user_site <- file.path(pc$dir, "user-site-packages")
  pc$site_dirs <- c(pc$site, user_site)
  probe(pc)
  dir.create(user_site)   # first pip install --user
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("a replaced interpreter is checked again", {
  pc <- local_probe_cache()
  probe(pc)
  writeLines("another interpreter build", pc$py)
  Sys.setFileTime(pc$py, Sys.time() + 5)
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("checks older than 30 days are made again", {
  pc <- local_probe_cache()
  probe(pc)
  entries <- readRDS(pc$file)
  entries <- lapply(entries, function(e) { e$saved <- Sys.time() - 31 * 86400; e })
  saveRDS(entries, pc$file)
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("a check made by another probe script is not reused", {
  pc <- local_probe_cache()
  probe(pc)
  entries <- readRDS(pc$file)
  entries <- lapply(entries, function(e) { e$extra$code <- "print('older PredictProR')"; e })
  saveRDS(entries, pc$file)
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("a failed check is not kept", {
  pc <- local_probe_cache()
  pc$answer <- function() character()
  expect_null(probe(pc))
  new_session()
  expect_null(probe(pc))
  expect_identical(pc$runs, 2L)
})

# e.g. the Windows Store alias python3.exe, a candidate on every Windows
# machine, which fails each time it is started
test_that("an interpreter that fails is started once per session", {
  pc <- local_probe_cache()
  pc$answer <- function() character()
  expect_null(probe(pc))
  expect_null(probe(pc))
  expect_identical(pc$runs, 1L)
})

test_that("an unreadable cache file is ignored and replaced", {
  pc <- local_probe_cache()
  writeLines("not an rds file", pc$file)
  expect_true(isTRUE(probe(pc)$numpy))
  new_session()
  probe(pc)
  expect_identical(pc$runs, 1L)
})

test_that("PREDICTPRO_PROBE_CACHE=false keeps nothing on disk", {
  pc <- local_probe_cache()
  withr::local_envvar(PREDICTPRO_PROBE_CACHE = "false")
  probe(pc)
  expect_false(file.exists(pc$file))
  new_session()
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("parallel workers follow the session's PREDICTPRO_PROBE_CACHE", {
  withr::local_envvar(PREDICTPRO_PROBE_CACHE = "false")
  expect_identical(PredictProR:::gp_parallel_runtime_env()[["PREDICTPRO_PROBE_CACHE"]], "false")
})

test_that("forgetting an interpreter's checks after installing into it also drops them from disk", {
  pc <- local_probe_cache()
  probe(pc)
  PredictProR:::gp_python_forget_cached_probes(pc$py)
  probe(pc)
  expect_identical(pc$runs, 2L)
})

test_that("conda environments are listed once per session and reused until one is added", {
  skip_if_not_installed("reticulate")
  pc <- local_probe_cache()
  registry <- file.path(pc$dir, "environments.txt")
  writeLines(file.path(pc$dir, "envs", "predictgp"), registry)
  env_py <- file.path(pc$dir, "envs", "predictgp", "python.exe")
  dir.create(dirname(env_py), recursive = TRUE)
  writeLines("", env_py)
  listed <- 0L
  local_mocked_bindings(
    conda_list = function(...) {
      listed <<- listed + 1L
      data.frame(name = "predictgp", python = env_py, stringsAsFactors = FALSE)
    },
    .package = "reticulate"
  )
  local_mocked_bindings(gp_conda_registry_files = function() registry, .package = "PredictProR")

  envs <- c("predictgp", "predictdl")
  expect_true(env_py %in% PredictProR:::gp_conda_python_candidates(envs))
  PredictProR:::gp_conda_python_candidates(envs)
  expect_identical(listed, 1L)

  new_session()
  expect_true(env_py %in% PredictProR:::gp_conda_python_candidates(envs))
  expect_identical(listed, 1L)

  # conda create registers the new environment
  cat(file.path(pc$dir, "envs", "predictdl"), "\n", file = registry, append = TRUE)
  Sys.setFileTime(registry, Sys.time() + 5)
  new_session()
  PredictProR:::gp_conda_python_candidates(envs)
  expect_identical(listed, 2L)
})

test_that("an interactive attach after the first starts no Python and lists no conda environments", {
  skip_if_not_installed("reticulate")
  pc <- local_probe_cache()
  withr::local_envvar(
    PREDICTPRO_PYTHON = NA, PREDICTPRO_GP_PYTHON = NA, PREDICTPRO_DL_PYTHON = NA,
    PREDICTPRO_PYTHON_DL = NA, RETICULATE_PYTHON = NA, `_R_CHECK_PACKAGE_NAME_` = NA,
    WORKON_HOME = file.path(pc$dir, "virtualenvs")
  )
  venv_py <- file.path(pc$dir, "virtualenvs", "predictdl", PredictProR:::gp_python_executable_relpath())
  dir.create(dirname(venv_py), recursive = TRUE)
  writeLines("", venv_py)
  withr::local_dir(pc$dir)
  listed <- 0L
  local_mocked_bindings(
    conda_list = function(...) {
      listed <<- listed + 1L
      data.frame(name = character(), python = character(), stringsAsFactors = FALSE)
    },
    .package = "reticulate"
  )
  local_mocked_bindings(
    gp_conda_registry_files = function() file.path(pc$dir, "environments.txt"),
    gp_system_python_candidates = function() character(),
    .package = "PredictProR"
  )
  attach_interactive <- PredictProR:::.onAttach
  environment(attach_interactive) <- list2env(list(interactive = function() TRUE),
                                              parent = asNamespace("PredictProR"))
  attach_message <- function() {
    msg <- character()
    withCallingHandlers(attach_interactive(NULL, "PredictProR"),
                        packageStartupMessage = function(m) {
                          msg <<- c(msg, conditionMessage(m))
                          invokeRestart("muffleMessage")
                        })
    msg
  }

  first <- attach_message()
  runs_first <- pc$runs
  expect_gt(runs_first, 0L)
  expect_identical(listed, 1L)
  expect_match(paste(first, collapse = "\n"), "Configured DL Python", fixed = TRUE)

  new_session()
  second <- attach_message()
  expect_identical(pc$runs, runs_first)
  expect_identical(listed, 1L)
  expect_identical(second, first)
})
