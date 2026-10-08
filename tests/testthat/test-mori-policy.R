test_that("mori policy is enabled by default and can be disabled", {
  old <- Sys.getenv("GP_USE_MORI", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("GP_USE_MORI") else Sys.setenv(GP_USE_MORI = old)
  }, add = TRUE)
  Sys.unsetenv("GP_USE_MORI")

  expect_true(PredictProR:::gp_mori_enabled())

  Sys.setenv(GP_USE_MORI = "false")
  expect_false(PredictProR:::gp_mori_enabled())
  expect_false(PredictProR:::gp_mori_should_share(matrix(0, nrow = 100, ncol = 100), backend = "mirai"))
})

test_that("mori backend eligibility follows backend policy", {
  expect_true(PredictProR:::gp_mori_backend_eligible("mirai", NULL))
  expect_true(PredictProR:::gp_mori_backend_eligible("base_parallel", NULL))
  expect_true(PredictProR:::gp_mori_backend_eligible("foreach", NULL))
  expect_true(PredictProR:::gp_mori_backend_eligible("future", future::multisession))
  expect_false(PredictProR:::gp_mori_backend_eligible("future", future::multicore))
  expect_false(PredictProR:::gp_mori_backend_eligible("sequential", NULL))
})

test_that("shared globals can be rebound into task closure environment", {
  external_big <- matrix(1, nrow = 10, ncol = 10)
  run_task <- local({
    external_big <- matrix(2, nrow = 2, ncol = 2)
    function(i) c(sum = sum(external_big), rows = nrow(external_big), idx = i)
  })

  rebound <- PredictProR:::gp_sp_bind_shared_globals(run_task, list(external_big = external_big))
  out <- rebound(7L)

  expect_equal(unname(out[["sum"]]), 100)
  expect_equal(unname(out[["rows"]]), 10)
  expect_equal(unname(out[["idx"]]), 7)
})

test_that("sp prepare globals falls back safely when mori is unavailable or disabled", {
  old_use <- Sys.getenv("GP_USE_MORI", unset = NA_character_)
  old_thr <- Sys.getenv("GP_MORI_MIN_MB", unset = NA_character_)
  on.exit({
    if (is.na(old_use)) Sys.unsetenv("GP_USE_MORI") else Sys.setenv(GP_USE_MORI = old_use)
    if (is.na(old_thr)) Sys.unsetenv("GP_MORI_MIN_MB") else Sys.setenv(GP_MORI_MIN_MB = old_thr)
  }, add = TRUE)

  Sys.setenv(GP_USE_MORI = "false", GP_MORI_MIN_MB = "0")
  globals <- list(big = matrix(rnorm(1000), nrow = 100))
  decision <- list(backend = "mirai", plan = NULL)
  prep <- PredictProR:::gp_sp_prepare_globals(globals, decision)

  expect_true(is.list(prep$globals))
  expect_true(is.matrix(prep$globals$big))
  expect_true(is.data.frame(prep$mori_summary))
  expect_false(prep$mori_summary$shared[[1]])
  expect_match(prep$mori_summary$reason[[1]], "^gate_|^score_")
})

test_that("mori decision engine accounts for size, workers, and task reuse", {
  old_use <- Sys.getenv("GP_USE_MORI", unset = NA_character_)
  old_thr <- Sys.getenv("GP_MORI_MIN_MB", unset = NA_character_)
  on.exit({
    if (is.na(old_use)) Sys.unsetenv("GP_USE_MORI") else Sys.setenv(GP_USE_MORI = old_use)
    if (is.na(old_thr)) Sys.unsetenv("GP_MORI_MIN_MB") else Sys.setenv(GP_MORI_MIN_MB = old_thr)
  }, add = TRUE)

  Sys.setenv(GP_USE_MORI = "true", GP_MORI_MIN_MB = "1")
  skip_if_not(PredictProR:::gp_mori_available(), "mori not installed")

  small <- matrix(1, nrow = 20, ncol = 20)
  large <- matrix(1, nrow = 2000, ncol = 2000)

  dec_small <- PredictProR:::gp_mori_decide_one(small, backend = "future", plan = future::multisession, workers = 2L, n_tasks = 2L, name = "small")
  dec_large <- PredictProR:::gp_mori_decide_one(large, backend = "future", plan = future::multisession, workers = 2L, n_tasks = 8L, name = "large")

  expect_false(dec_small$use_mori)
  expect_true(dec_large$use_mori)
  expect_gt(dec_large$score, dec_small$score)
  expect_gt(dec_large$copy_volume_mb, dec_small$copy_volume_mb)
})

test_that("mori can mark eligible large objects when installed", {
  skip_if_not_installed("mori")

  old_use <- Sys.getenv("GP_USE_MORI", unset = NA_character_)
  old_thr <- Sys.getenv("GP_MORI_MIN_MB", unset = NA_character_)
  on.exit({
    if (is.na(old_use)) Sys.unsetenv("GP_USE_MORI") else Sys.setenv(GP_USE_MORI = old_use)
    if (is.na(old_thr)) Sys.unsetenv("GP_MORI_MIN_MB") else Sys.setenv(GP_MORI_MIN_MB = old_thr)
  }, add = TRUE)

  Sys.setenv(GP_USE_MORI = "true", GP_MORI_MIN_MB = "0")
  globals <- list(big = matrix(rnorm(2000), nrow = 100))
  decision <- list(backend = "mirai", plan = NULL)
  prep <- PredictProR:::gp_sp_prepare_globals(globals, decision)

  expect_true(is.data.frame(prep$mori_summary))
  expect_true(prep$mori_summary$shared[[1]])
})
