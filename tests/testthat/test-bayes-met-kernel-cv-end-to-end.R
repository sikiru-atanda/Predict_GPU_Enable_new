# End-to-end test for kernel-Bayes MET CV (Phase 3.20 regression guard).
# The existing test-bayes-met-kernel-cv-routing.R only mocks
# bayes_mod_cv_met_kernel_predict and bayes_mod_cv. This test invokes the
# real BGLR::Multitrait sampler with a small synthetic dataset so that the
# (gid, env) key lookup and per-fold response handling are actually exercised.

test_that("bayes_mod_cv_met_kernel_predict returns non-NA preds for single-trait MET CV", {
  skip_if_not_installed("BGLR")

  set.seed(1)
  n_gid <- 12
  envs <- c("E1","E2","E3")
  pheno <- expand.grid(GID = paste0("g", seq_len(n_gid)),
                       Env = envs, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
  pheno$Yield <- with(pheno, as.numeric(factor(GID)) * 0.2 + rnorm(nrow(pheno), 0, 0.5))
  K <- diag(n_gid); rownames(K) <- colnames(K) <- paste0("g", seq_len(n_gid))

  meta <- list(
    pheno_data = pheno,
    response = "Yield",
    gen_name = "GID",
    heter_groups = "Env",
    kernels = list(gmatrix = K),
    GS_model = "RKHS",
    bayes_para = list(nIter = 200, burnIn = 50, thin = 2),
    omics_kernel_label = NULL
  )

  tst <- c(1L, 7L, 15L, 22L, 30L)
  y <- pheno$Yield
  y[tst] <- NA_real_

  preds <- suppressMessages(suppressWarnings(
    PredictProR:::bayes_mod_cv_met_kernel_predict(
      y = y, tst = tst, met_kernel_cv_meta = meta)
  ))

  expect_length(preds, length(tst))
  expect_true(sum(!is.na(preds)) == length(tst),
              info = paste("Expected all non-NA, got", sum(is.na(preds)), "NA"))
  expect_true(is.numeric(preds))
})

test_that("bayes_mod_cv single-trait-from-multi-trait override prevents recursive-indexing crash", {
  # Phase 3.18 regression guard: meta$response carries the full multi-trait
  # response vector from prep, dispatcher injects bayes_trait per-fold, the
  # override must collapse meta$response to that one trait before calling
  # bayes_mod_cv_met_kernel_predict.
  skip_if_not_installed("BGLR")

  set.seed(2)
  n_gid <- 10
  envs <- c("E1","E2")
  pheno <- expand.grid(GID = paste0("g", seq_len(n_gid)),
                       Env = envs, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
  pheno$Yield <- rnorm(nrow(pheno))
  pheno$BLUE  <- rnorm(nrow(pheno))
  pheno$deBLUP <- rnorm(nrow(pheno))
  K <- diag(n_gid); rownames(K) <- colnames(K) <- paste0("g", seq_len(n_gid))

  meta <- list(
    pheno_data = pheno,
    response = c("Yield","BLUE","deBLUP"),  # multi-trait, length 3 from prep
    gen_name = "GID",
    heter_groups = "Env",
    kernels = list(gmatrix = K),
    GS_model = "BRR",
    bayes_para = list(nIter = 200, burnIn = 50, thin = 2),
    omics_kernel_label = NULL
  )

  y <- pheno$Yield; tst <- c(2L, 5L, 11L, 17L)
  y[tst] <- NA_real_

  # Calling bayes_mod_cv (not the predictor directly) -- this is the dispatcher
  # entry point and must apply the Phase 3.18 override.
  preds <- suppressMessages(suppressWarnings(
    bayes_mod_cv(y = y, ETA = list(), weights = NULL,
                 bayes_para = meta$bayes_para, tst = tst,
                 bayes_model = "BRR", bayes_trait = "Yield",
                 response_family = "gaussian", groups = NULL,
                 met_kernel_cv_meta = meta)
  ))

  expect_length(preds, length(tst))
  expect_true(all(!is.na(preds)),
              info = paste("All-NA predictions (recursive-indexing or downstream bug). NAs:",
                           sum(is.na(preds))))
})
