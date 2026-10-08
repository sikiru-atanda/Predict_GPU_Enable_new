# Pins the cross-validation parity for the kernel-Bayesian MET heter_resid
# fix (commit 35a3bcb covered true-prediction only). When all four conditions
# hold -- gaussian, heter_resid=TRUE, kernel-Bayes (RKHS/GBLUP_BRR/BRR), MET
# panel (heter_groups labels >= 2 envs with repeated genotype rows) -- CV
# folds must take the multitrait kernel-prior path, not BGLR::BGLR with
# native `groups`.

# ---- gp_bayes_should_route_met_kernel_cv -------------------------------------

test_that("gp_bayes_should_route_met_kernel_cv true for kernel-Bayes MET heter_resid", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3","g4"), each = 2),
    Env = rep(c("E1","E2"), times = 4),
    Yield = rnorm(8),
    stringsAsFactors = FALSE
  )
  expect_true(gp_bayes_should_route_met_kernel_cv(
    GS_model = "RKHS", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env"))
  expect_true(gp_bayes_should_route_met_kernel_cv(
    GS_model = "GBLUP_BRR", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env"))
  expect_true(gp_bayes_should_route_met_kernel_cv(
    GS_model = "BRR", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env"))
})

test_that("gp_bayes_should_route_met_kernel_cv false when heter_resid is FALSE", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3"), each = 2),
    Env = rep(c("E1","E2"), times = 3),
    Yield = rnorm(6),
    stringsAsFactors = FALSE
  )
  expect_false(gp_bayes_should_route_met_kernel_cv(
    GS_model = "RKHS", heter_resid = FALSE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env"))
})

test_that("gp_bayes_should_route_met_kernel_cv false for single-environment panels", {
  pheno <- data.frame(
    GID = c("g1","g2","g3"),
    Env = "E1",
    Yield = rnorm(3),
    stringsAsFactors = FALSE
  )
  expect_false(gp_bayes_should_route_met_kernel_cv(
    GS_model = "RKHS", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env"))
})

test_that("gp_bayes_should_route_met_kernel_cv false for marker-Bayesian models", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3"), each = 2),
    Env = rep(c("E1","E2"), times = 3),
    Yield = rnorm(6),
    stringsAsFactors = FALSE
  )
  for (m in c("BayesA","BayesB","BayesC","BL")) {
    expect_false(gp_bayes_should_route_met_kernel_cv(
      GS_model = m, heter_resid = TRUE, pheno_data = pheno,
      gen_name = "GID", heter_groups = "Env"))
  }
})

test_that("gp_bayes_should_route_met_kernel_cv false for non-Gaussian families", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3"), each = 2),
    Env = rep(c("E1","E2"), times = 3),
    Yield = rnorm(6),
    stringsAsFactors = FALSE
  )
  expect_false(gp_bayes_should_route_met_kernel_cv(
    GS_model = "RKHS", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env",
    response_family = "ordinal"))
  expect_false(gp_bayes_should_route_met_kernel_cv(
    GS_model = "RKHS", heter_resid = TRUE, pheno_data = pheno,
    gen_name = "GID", heter_groups = "Env",
    response_family = "binary"))
})

# ---- bayes_mod_cv routing on met_kernel_cv_meta ------------------------------

test_that("bayes_mod_cv dispatches to met_kernel_cv_predict when meta is supplied", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3","g4"), each = 2),
    Env = rep(c("E1","E2"), times = 4),
    Yield = c(1.1, 1.4, 2.0, 1.9, 3.1, 3.4, 4.2, 4.0),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- c("g1","g2","g3","g4")
  tst <- c(2L, 6L) # mask two rows

  seen <- list(called = FALSE, n_tst = NA_integer_, meta_gs = NA_character_)
  local_mocked_bindings(
    bayes_mod_cv_met_kernel_predict = function(y, tst, met_kernel_cv_meta) {
      seen$called <<- TRUE
      seen$n_tst <<- length(tst)
      seen$meta_gs <<- met_kernel_cv_meta$GS_model
      # Return placeholder predictions of the right length
      rep(0.0, length(tst))
    },
    .package = "PredictProR"
  )

  y <- pheno$Yield; y[tst] <- NA_real_
  meta <- list(
    pheno_data = pheno, response = "Yield", gen_name = "GID",
    heter_groups = "Env", kernels = list(K = K), GS_model = "BRR",
    bayes_para = list(nIter = 200, burnIn = 50, thin = 5)
  )

  preds <- bayes_mod_cv(
    y = y, ETA = list(), weights = NULL,
    bayes_para = meta$bayes_para,
    tst = tst, bayes_model = "BRR", bayes_trait = "Yield",
    response_family = "gaussian", groups = NULL,
    met_kernel_cv_meta = meta
  )

  expect_true(seen$called)
  expect_identical(seen$n_tst, length(tst))
  expect_identical(seen$meta_gs, "BRR")
  expect_length(preds, length(tst))
})

test_that("bayes_mod_cv uses the legacy BGLR path when met_kernel_cv_meta is NULL", {
  # Guards against accidentally always routing through the multitrait path
  # (which would break single-env or marker-Bayes CV fold fits).
  seen_native <- list(called = FALSE)
  local_mocked_bindings(
    gp_bayes_bglr_fit = function(...) {
      seen_native$called <<- TRUE
      list(yHat = rep(0.0, 6), mu = 0, ETA = list(), levels = NULL, probs = NULL)
    },
    bayes_mod_cv_met_kernel_predict = function(...) {
      stop("met_kernel_cv_predict should NOT be called when met_kernel_cv_meta is NULL")
    },
    .package = "PredictProR"
  )

  y <- c(1, 2, NA, 4, 5, NA)
  bayes_mod_cv(
    y = y, ETA = list(), weights = NULL,
    bayes_para = list(nIter = 200, burnIn = 50, thin = 5),
    tst = c(3L, 6L), bayes_model = "RKHS", bayes_trait = "y",
    response_family = "gaussian", groups = NULL,
    met_kernel_cv_meta = NULL
  )
  expect_true(seen_native$called)
})

# ---- end-to-end: model_prep_bayes_cv attaches the meta -----------------------

test_that("model_prep_bayes_cv attaches met_kernel_cv_meta for kernel-Bayes MET heter_resid", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3","g4"), each = 2),
    Env = rep(c("E1","E2"), times = 4),
    Yield = rnorm(8),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- c("g1","g2","g3","g4")
  prep <- model_prep_bayes_cv(
    random = ~ GID + GID:Env, GS_model_cv = "GBLUP_BRR",
    response = "Yield", gen_name = "GID", pheno_data = pheno,
    gmatrix = K, heter_groups = "Env", heter_resid = TRUE,
    nIter = 200, burnIn = 50, thin = 5,
    cross_validation = TRUE, response_family = "gaussian"
  )
  expect_true("met_kernel_cv_meta" %in% names(prep[["GBLUP_BRR"]]))
  meta <- prep[["GBLUP_BRR"]][["met_kernel_cv_meta"]]
  expect_identical(meta$GS_model, "BRR")  # collapsed alias
  expect_identical(meta$response, "Yield")
  expect_identical(meta$heter_groups, "Env")
  expect_identical(meta$gen_name, "GID")
  expect_true(length(meta$kernels) >= 1L)
})

test_that("model_prep_bayes_cv omits met_kernel_cv_meta when heter_resid is FALSE", {
  pheno <- data.frame(
    GID = rep(c("g1","g2","g3","g4"), each = 2),
    Env = rep(c("E1","E2"), times = 4),
    Yield = rnorm(8),
    stringsAsFactors = FALSE
  )
  K <- diag(4); rownames(K) <- colnames(K) <- c("g1","g2","g3","g4")
  prep <- model_prep_bayes_cv(
    random = ~ GID + GID:Env, GS_model_cv = "RKHS",
    response = "Yield", gen_name = "GID", pheno_data = pheno,
    gmatrix = K, heter_groups = "Env", heter_resid = FALSE,
    nIter = 200, burnIn = 50, thin = 5,
    cross_validation = TRUE, response_family = "gaussian"
  )
  expect_false("met_kernel_cv_meta" %in% names(prep[["RKHS"]]))
})
