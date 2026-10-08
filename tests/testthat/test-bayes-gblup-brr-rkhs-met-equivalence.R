# Pins the documented MET equivalence between GBLUP_BRR and RKHS:
# both are routed through bayes_multitrait_env_heter_fit -> BGLR::Multitrait
# with model = "RKHS", so under X = eigen-sqrt(K) they produce numerically
# identical fits when given the same RNG seed. Also pins that the per-fit
# BGLR save-prefix is unique enough to avoid collision when both run in the
# same session.

skip_if_no_bglr <- function() {
  testthat::skip_if_not_installed("BGLR")
}

make_tiny_met <- function() {
  set.seed(1L)
  gid <- paste0("g", 1:8)
  envs <- c("E1", "E2")
  pheno <- expand.grid(GID = gid, Env = envs, stringsAsFactors = FALSE)
  pheno$Yield <- rnorm(nrow(pheno))
  K <- diag(length(gid))
  rownames(K) <- colnames(K) <- gid
  list(pheno = pheno, K = K)
}

test_that("MET heter_resid GBLUP_BRR and RKHS produce identical fits under the same seed", {
  skip_if_no_bglr()
  dat <- make_tiny_met()

  fit_one <- function(gs, seed) {
    set.seed(seed)
    suppressMessages(suppressWarnings(
      bayes_finalize_RKHS_GBLUPBRR(
        random = ~ GID + GID:Env,
        GS_model = gs,
        response = "Yield",
        pheno_data = dat$pheno,
        gmatrix = dat$K,
        gen_name = "GID",
        heter_groups = "Env",
        heter_resid = TRUE,
        nIter = 300,
        burnIn = 100,
        thin = 5
      )
    ))$bglr_multitrait_model
  }

  a <- fit_one("GBLUP_BRR", 42L)
  b <- fit_one("RKHS",      42L)

  # mu and ETAHat must match to numerical precision: same kernel, same y,
  # same priors, same MCMC RNG state at BGLR entry -> same chain.
  expect_equal(as.double(a$mu), as.double(b$mu), tolerance = 1e-10)
  expect_equal(as.vector(a$ETAHat), as.vector(b$ETAHat), tolerance = 1e-10)
  expect_equal(diag(a$resCov$R), diag(b$resCov$R), tolerance = 1e-10)
})

test_that("MET heter_resid GBLUP_BRR responds to R's set.seed (different seeds -> different draws)", {
  # Guards against accidentally freezing the chain (a regression where every
  # call returns the same numbers regardless of seed would mean an upstream
  # seed-override is masking BGLR's RNG dependence).
  skip_if_no_bglr()
  dat <- make_tiny_met()

  fit_one <- function(seed) {
    set.seed(seed)
    suppressMessages(suppressWarnings(
      bayes_finalize_RKHS_GBLUPBRR(
        random = ~ GID + GID:Env, GS_model = "GBLUP_BRR",
        response = "Yield", pheno_data = dat$pheno, gmatrix = dat$K,
        gen_name = "GID", heter_groups = "Env", heter_resid = TRUE,
        nIter = 300, burnIn = 100, thin = 5
      )
    ))$bglr_multitrait_model
  }

  a <- fit_one(42L)
  b <- fit_one(999L)

  # At least some difference somewhere in mu / ETAHat / resCov$R; the chain
  # is not constant-by-construction.
  any_diff <- max(abs(as.double(a$mu) - as.double(b$mu))) +
              max(abs(as.vector(a$ETAHat) - as.vector(b$ETAHat))) +
              max(abs(diag(a$resCov$R) - diag(b$resCov$R)))
  expect_gt(any_diff, 0)
})

test_that("Sequential GBLUP_BRR + RKHS MET fits do not collide on BGLR save_prefix", {
  skip_if_no_bglr()
  dat <- make_tiny_met()

  before <- list.files(tempdir(), pattern = "^BGLR_Multitrait_", full.names = FALSE)
  for (gs in c("GBLUP_BRR", "RKHS")) {
    set.seed(7L)
    suppressMessages(suppressWarnings(
      bayes_finalize_RKHS_GBLUPBRR(
        random = ~ GID + GID:Env, GS_model = gs,
        response = "Yield", pheno_data = dat$pheno, gmatrix = dat$K,
        gen_name = "GID", heter_groups = "Env", heter_resid = TRUE,
        nIter = 200, burnIn = 50, thin = 5
      )
    ))
  }
  after <- list.files(tempdir(), pattern = "^BGLR_Multitrait_", full.names = FALSE)
  new_files <- setdiff(after, before)

  # The fitter is expected to clean up its temp files via on.exit + unlink, so
  # the post-set may be empty here; what matters is that the prefix builder
  # generates distinct names (timestamp + PID + random int + GS_model label).
  # Test the builder directly:
  p1 <- gp_bglr_save_prefix(model_name = "BGLR_Multitrait_GBLUP_BRR_UN", response = "Yield")
  p2 <- gp_bglr_save_prefix(model_name = "BGLR_Multitrait_RKHS_UN",      response = "Yield")
  expect_false(identical(p1, p2))
  expect_true(grepl("GBLUP_BRR", basename(p1), fixed = TRUE))
  expect_true(grepl("_RKHS_",    basename(p2), fixed = TRUE))
})
