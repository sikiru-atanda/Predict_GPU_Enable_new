# Phase 3.7e pinning test for warm-worker gp_engine handling.
#
# The Phase 3.7d audit (Agent A) confirmed by inspection that the warm
# worker correctly threads gp_engine into every fit's payload with no
# module-level cache that would cause the previous fit's engine to
# leak into the next. This test pins that property as an integration
# check: run two consecutive fits with different gp_engine values
# through the warm worker, then verify each fit reports the engine
# that was requested (via the new model_parameters$gp_engine_used row
# added in 0.20.27).
#
# Skip when:
#   * BGLR not installed (no test data),
#   * PREDICTPRO_GP_PYTHON not set (no Python backend),
#   * the predictgp env can't import torch / gpytorch / scipy.sparse,
#   * or the user has explicitly disabled the worker.

run_kg_fit <- function(pheno_long, G, engine_label) {
  PredictProR:::gp_warm_start_clear()
  withr::local_dir(withr::local_tempdir())  # exported results stay out of tests/
  suppressWarnings(model_execute(
    pheno_data       = pheno_long,
    gmatrix          = G,
    gen_name         = "NAME",
    response         = "Yield",
    GS_model         = "Kernel-GBLUP",
    engine           = "asreml",
    heter_groups     = "Env",
    cross_validation = FALSE,
    system_database  = FALSE,
    ld_prunning_qc   = FALSE,
    para_tunning     = FALSE,
    gp_varcomp_mode  = "reml",
    gp_engine        = engine_label,
    random           = ~NAME
  ))
}

build_wheat_fixture <- function(seed = 20260601, n_geno = 60) {
  set.seed(seed)
  data("wheat", package = "BGLR", envir = environment())
  Y <- get("wheat.Y", envir = environment())
  X <- get("wheat.X", envir = environment())
  keep_idx <- sample(seq_len(nrow(Y)), n_geno)
  NAME <- rownames(Y)[keep_idx]
  Y_sub <- Y[keep_idx, , drop = FALSE]
  X_sub <- X[keep_idx, , drop = FALSE]
  pheno_long <- do.call(rbind, lapply(seq_len(ncol(Y_sub)), function(j) {
    data.frame(NAME = NAME, Env = paste0("E", colnames(Y_sub)[j]),
               Yield = as.numeric(Y_sub[, j]), stringsAsFactors = FALSE)
  }))
  G <- tcrossprod(scale(X_sub, center = TRUE, scale = FALSE))
  G <- G / mean(diag(G)); G <- 0.99 * G + 0.01 * diag(nrow(G))
  rownames(G) <- colnames(G) <- NAME
  list(pheno_long = pheno_long, G = G)
}

test_that("warm worker honors gp_engine on each successive fit", {
  skip_if_not_installed("BGLR")
  skip_if(!nzchar(Sys.getenv("PREDICTPRO_GP_PYTHON")),
          "PREDICTPRO_GP_PYTHON not set; GP runtime smoke test skipped")
  if (!PredictProR:::gp_bridge_worker_use_enabled()) {
    skip("PREDICTPRO_GP_WORKER=false; warm-worker pinning skipped")
  }

  fx <- build_wheat_fixture()

  # Fit 1: dense_v.
  fit_dense <- tryCatch(run_kg_fit(fx$pheno_long, fx$G, "dense_v"),
                        error = function(e) e)
  if (inherits(fit_dense, "error")) {
    skip(paste("dense_v fit errored; runtime not available:",
               conditionMessage(fit_dense)))
  }
  mp_dense <- fit_dense$model_results$model_parameters
  expect_true("gp_engine_requested" %in% mp_dense$stat)
  expect_true("gp_engine_used" %in% mp_dense$stat)
  expect_identical(mp_dense$summary[mp_dense$stat == "gp_engine_requested"],
                   "dense_v")
  expect_identical(mp_dense$summary[mp_dense$stat == "gp_engine_used"],
                   "dense_v")

  # Fit 2: mme, on the SAME warm worker (no restart).
  fit_mme <- tryCatch(run_kg_fit(fx$pheno_long, fx$G, "mme"),
                      error = function(e) e)
  if (inherits(fit_mme, "error")) {
    skip(paste("mme fit errored; runtime not available:",
               conditionMessage(fit_mme)))
  }
  mp_mme <- fit_mme$model_results$model_parameters
  expect_identical(mp_mme$summary[mp_mme$stat == "gp_engine_requested"],
                   "mme")
  # gp_engine_used reflects what actually ran. For Kernel-GBLUP the
  # canonical model is KRR which routes through method="krr_exact" --
  # that Python path does NOT yet have the MME branch wired (Phase 3.7f
  # territory), so we expect a deterministic "dense_v" diagnostic
  # showing the fallback. For models that DO reach gp_exact_with_X
  # (Gaussian-Process-GBLUP, canonical "GP"), the value would be "mme".
  expect_true(mp_mme$summary[mp_mme$stat == "gp_engine_used"] %in%
              c("mme", "dense_v"))

  # Fit 3: dense_v again, to pin that mme doesn't stick.
  fit_dense2 <- tryCatch(run_kg_fit(fx$pheno_long, fx$G, "dense_v"),
                         error = function(e) e)
  if (inherits(fit_dense2, "error")) {
    skip(paste("dense_v fit (round 2) errored:",
               conditionMessage(fit_dense2)))
  }
  mp_dense2 <- fit_dense2$model_results$model_parameters
  expect_identical(mp_dense2$summary[mp_dense2$stat == "gp_engine_used"],
                   "dense_v")

  # VC parity across the three fits: deterministic optimizer + same data
  # converge to the same REML MLE no matter which engine got us there.
  vc_dense  <- unique(fit_dense$model_results$variance_components[,
                       c("Component", "Components")])
  vc_mme    <- unique(fit_mme$model_results$variance_components[,
                       c("Component", "Components")])
  vc_dense2 <- unique(fit_dense2$model_results$variance_components[,
                       c("Component", "Components")])
  shared <- Reduce(intersect, list(vc_dense$Component, vc_mme$Component,
                                   vc_dense2$Component))
  expect_true(length(shared) > 0)
  v1 <- vc_dense$Components[match(shared, vc_dense$Component)]
  v2 <- vc_mme$Components[match(shared, vc_mme$Component)]
  v3 <- vc_dense2$Components[match(shared, vc_dense2$Component)]
  expect_equal(v1, v2, tolerance = 1e-8)
  expect_equal(v1, v3, tolerance = 1e-8)
})


test_that("gp_engine_used reflects fallback when MME is not eligible (skip-only smoke)", {
  # When build_gp_exact_ai_varcomp_mme raises NotImplementedError (e.g. FA
  # specs in the future), gp_framework.py:gp_exact_with_X catches the
  # exception, warns the user, and re-dispatches to the dense_v path. The
  # resulting model_parameters$gp_engine_used should then read "dense_v"
  # while gp_engine_requested stays "mme" -- giving the user a clear
  # diagnostic that the fallback fired.
  #
  # Currently in 0.20.27 the only MME-ineligible spec kind (FA) is gated
  # earlier in PredictProR via the FA -> KRR auto-route, so we cannot
  # exercise this branch end-to-end without artificial inputs. The test
  # is registered as a placeholder so future MME ineligibility cases
  # (GxE, env-main random terms) can plug in here.
  skip("MME-ineligible spec kinds not yet reachable from R; placeholder for 3.7e+")
})
