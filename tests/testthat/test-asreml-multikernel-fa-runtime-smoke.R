if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

local_asreml_fa_runtime <- function() {
  skip_if_not_installed("asreml")
  status <- tryCatch(asreml::asreml.license.status(), error = function(e) NULL)
  if (is.null(status) || is.null(status$status) || is.na(status$status) ||
      as.integer(status$status) < 0L) {
    msg <- if (!is.null(status$statusMessage) && nzchar(status$statusMessage)) {
      status$statusMessage
    } else {
      "ASReml-R license unavailable"
    }
    skip(msg)
  }
}

mk_asreml_multikernel_fa_data <- function(n = 42L, p = 18L) {
  set.seed(814)
  gids <- paste0("g", seq_len(n))
  envs <- paste0("E", seq_len(3L))

  make_kernel <- function(seed) {
    set.seed(seed)
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    x <- scale(x)
    k <- tcrossprod(x) / ncol(x)
    k <- k / mean(diag(k))
    diag(k) <- diag(k) + 0.15
    dimnames(k) <- list(gids, gids)
    k
  }

  k1 <- make_kernel(815)
  k2 <- make_kernel(816)
  u1 <- as.numeric(crossprod(chol(k1), rnorm(n)))
  u2 <- as.numeric(crossprod(chol(k2), rnorm(n)))
  lambda1 <- c(1.00, 0.72, 0.48)
  lambda2 <- c(0.35, 0.75, 1.05)

  pheno <- expand.grid(
    GID = gids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  g <- match(pheno$GID, gids)
  e <- match(pheno$Env, envs)
  pheno$Yield <- c(0.0, 0.5, -0.3)[e] +
    lambda1[e] * u1[g] + lambda2[e] * u2[g] + rnorm(nrow(pheno), sd = 0.35)
  pheno$Yield[c(8L, 39L, 91L)] <- NA_real_
  pheno$GID <- factor(pheno$GID, levels = gids)
  pheno$Env <- factor(pheno$Env, levels = envs)

  list(pheno = pheno, genomic = k1, omic = k2)
}

test_that("licensed ASReml-R fits and reconstructs every kernel in an FA1 MET", {
  local_asreml_fa_runtime()
  dat <- mk_asreml_multikernel_fa_data()

  fit <- PredictProR:::asreml_utilis_new(
    fixed = ~ Env,
    random = ~ GID + GID:Env,
    GS_model = "GBLUP",
    response = "Yield",
    pheno_data = dat$pheno,
    gmatrix = dat$genomic,
    kernel_list = list(omic = dat$omic),
    inverse = TRUE,
    gen_name = "GID",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "fa1",
    engine = "asreml",
    maxit = 100L
  )

  fa_terms <- gregexpr("fa\\(Env,1\\):vm\\(GID,", fit$str_mod)[[1L]]
  expect_length(fa_terms[fa_terms > 0L], 2L)
  expect_true(isTRUE(fit$model$converge))

  vc <- PredictProR:::asreml_herit_varCov_new(
    model = fit$model,
    heter_groups = "Env",
    var_cov_str = "fa1",
    heter_resid = TRUE,
    names_in_inv_list = fit$names_in_inv_list,
    inter_gen_pos = fit$inter_gen_pos,
    gen_pos = fit$gen_pos
  )

  expect_length(vc$Covariance, 2L)
  expect_equal(dim(vc$varG_per_omics), c(2L, 3L))
  expect_true(all(is.finite(vc$varG_per_omics)))
  expect_true(all(is.finite(vc$Heritability)))
  expect_true(all(vc$Heritability >= 0 & vc$Heritability <= 1))

  for (covariance in vc$Covariance) {
    expect_equal(covariance, t(covariance), tolerance = 1e-8)
    expect_true(min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  }
})
