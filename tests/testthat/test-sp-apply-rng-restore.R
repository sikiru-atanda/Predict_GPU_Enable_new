test_that("sequential sp_apply leaves the caller's RNG kind and seed unchanged", {
  # Sequential tasks run in the caller's session and each installs an
  # L'Ecuyer-CMRG stream seed into .Random.seed. That switched the session's
  # RNG kind and never restored it, so after a model_execute() call the
  # user's own set.seed() code produced different numbers.
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  out <- suppressMessages(PredictProR:::sp_apply(
    X = as.list(1:3),
    FUN = function(x) stats::runif(1),
    decision = list(backend = "sequential"),
    seed = 123L
  ))

  expect_length(out, 3L)
  expect_identical(RNGkind(), kind_before)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), seed_before)
  # and the caller's stream continues exactly as if sp_apply had not run
  expect_identical(stats::runif(2), {
    assign(".Random.seed", seed_before, envir = .GlobalEnv); stats::runif(2)
  })
})
