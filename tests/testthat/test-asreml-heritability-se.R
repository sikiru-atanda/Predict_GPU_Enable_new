# Rigor patch (ASReml finding #3): heritability was reported as a point estimate
# with no standard error. gp_asreml_heritability_se() supplies a delta-method SE
# via asreml::vpredict for the cases where the genetic and residual variances map
# directly to variance-component rows (compound symmetry, corgh diagonal). It must
# (a) match a direct vpredict call exactly, and (b) fall back to NA (never error,
# never a wrong number) when the requested components cannot be matched.

resolve_fn <- function(name) {
  if (exists(name, mode = "function")) return(get(name, mode = "function"))
  ns <- tryCatch(asNamespace("PredictProR"), error = function(e) NULL)
  if (!is.null(ns) && exists(name, envir = ns, mode = "function")) {
    return(get(name, envir = ns, mode = "function"))
  }
  NULL
}

asreml_ready <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) return(FALSE)
  tryCatch(grepl("No error", asreml::asreml.license.status(quiet = TRUE)$rv_text),
           error = function(e) FALSE)
}

# A small, well-conditioned single-env GBLUP that converges quickly.
fit_toy_gblup <- function() {
  inv_builder <- resolve_fn("compute_inverse_and_sparse")
  if (is.null(inv_builder)) return(NULL)
  set.seed(7L); n <- 40L; p <- 250L
  X <- matrix(stats::rnorm(n * p), n, p)
  rownames(X) <- paste0("g", seq_len(n))
  K <- tcrossprod(scale(X)) / p
  K <- K / mean(diag(K)); diag(K) <- diag(K) + 0.02
  g <- as.numeric(crossprod(chol(K), stats::rnorm(n)))
  ph <- data.frame(GID = factor(rownames(X)),
                   y = 5 + 2 * scale(g)[, 1] + stats::rnorm(n, sd = 1.5))
  Ginv <<- inv_builder(kernel = K, epsilon = 1e-6, inverse = TRUE)
  mod <- tryCatch(asreml::asreml(fixed = y ~ 1, random = ~ vm(GID, Ginv),
                                 data = ph, trace = FALSE),
                  error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  for (i in 1:3) if (!isTRUE(mod$converge)) mod <- asreml::update.asreml(mod) else break
  if (!isTRUE(mod$converge)) return(NULL)
  mod
}

test_that("gp_asreml_heritability_se matches a direct vpredict delta-method SE", {
  skip_if_not(asreml_ready(), "asreml + valid license required")
  fn <- resolve_fn("gp_asreml_heritability_se")
  skip_if(is.null(fn), "gp_asreml_heritability_se not available")
  mod <- fit_toy_gblup()
  skip_if(is.null(mod), "toy GBLUP did not converge")

  vc_rn <- rownames(summary(mod)$varcomp)
  vg <- grep("vm\\(GID", vc_rn, value = TRUE)
  ve <- grep("!R$", vc_rn, value = TRUE)
  expect_length(vg, 1L); expect_length(ve, 1L)

  se <- fn(mod, vg_rownames = vg, ve_rowname = ve)
  ig <- match(vg, vc_rn); ie <- match(ve, vc_rn)
  direct <- asreml::vpredict(mod, stats::as.formula(sprintf("h2 ~ V%d/(V%d+V%d)", ig, ig, ie)))$SE[1]

  expect_true(is.finite(se) && se > 0)
  expect_equal(se, direct, tolerance = 1e-6)
})

# A small 2-env corgh MET fit (heterogeneous genetic variances + correlation).
fit_toy_corgh <- function() {
  inv_builder <- resolve_fn("compute_inverse_and_sparse")
  if (is.null(inv_builder)) return(NULL)
  set.seed(11L); n <- 45L; p <- 300L
  X <- matrix(stats::rnorm(n * p), n, p); rownames(X) <- paste0("g", seq_len(n))
  K <- tcrossprod(scale(X)) / p; K <- K / mean(diag(K)); diag(K) <- diag(K) + 0.02
  g <- crossprod(chol(K), matrix(stats::rnorm(2 * n), n, 2))
  ph <- rbind(
    data.frame(GID = rownames(X), Env = "E1", y = 5 + 2 * scale(g[, 1])[, 1] + stats::rnorm(n, sd = 2.2)),
    data.frame(GID = rownames(X), Env = "E2", y = 6 + 2 * scale(g[, 2])[, 1] + stats::rnorm(n, sd = 2.5)))
  ph$GID <- factor(ph$GID); ph$Env <- factor(ph$Env); ph <- ph[order(ph$Env, ph$GID), ]
  Ginv <<- inv_builder(kernel = K, epsilon = 1e-6, inverse = TRUE)
  mod <- tryCatch(asreml::asreml(fixed = y ~ Env, random = ~ corgh(Env):vm(GID, Ginv),
                                 residual = ~ dsum(~ units | Env), data = ph, trace = FALSE),
                  error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  for (i in 1:3) if (!isTRUE(mod$converge)) mod <- asreml::update.asreml(mod) else break
  if (!isTRUE(mod$converge)) return(NULL)
  mod
}

test_that("gp_asreml_h2_se_for_env: exact SE for a corgh diagonal, NA when self-check fails", {
  skip_if_not(asreml_ready(), "asreml + valid license required")
  fn <- resolve_fn("gp_asreml_h2_se_for_env")
  skip_if(is.null(fn), "gp_asreml_h2_se_for_env not available")
  mod <- fit_toy_corgh()
  skip_if(is.null(mod), "toy corgh did not converge")

  vc <- summary(mod)$varcomp; rn <- rownames(vc)
  g_row <- grep("vm\\(GID.*E1$", rn, value = TRUE)
  e_row <- grep("E1!R$", rn, value = TRUE)
  expect_length(g_row, 1L); expect_length(e_row, 1L)
  vg <- vc[g_row, "component"]; ve <- vc[e_row, "component"]

  # Correct targets -> SE matches a direct vpredict call exactly.
  se <- fn(mod, vc = vc, names_in_inv_list = "Ginv", env_label = "E1",
           target_vg = vg, target_ve = ve)
  ig <- match(g_row, rn); ie <- match(e_row, rn)
  direct <- asreml::vpredict(mod, stats::as.formula(sprintf("h2 ~ V%d/(V%d+V%d)", ig, ig, ie)))$SE[1]
  expect_true(is.finite(se) && se > 0)
  expect_equal(se, direct, tolerance = 1e-6)

  # Self-check guard: a target genetic variance that does NOT match the row
  # (as happens for a factor-analytic diagonal Lambda Lambda' + psi) -> NA,
  # never a wrong SE.
  expect_true(is.na(fn(mod, vc = vc, names_in_inv_list = "Ginv", env_label = "E1",
                       target_vg = vg * 1.5, target_ve = ve)))
  expect_true(is.na(fn(mod, vc = vc, names_in_inv_list = "Ginv", env_label = "E1",
                       target_vg = vg, target_ve = ve * 2)))
})

test_that("gp_asreml_heritability_se returns NA (no error) when components are unmatched", {
  skip_if_not(asreml_ready(), "asreml + valid license required")
  fn <- resolve_fn("gp_asreml_heritability_se")
  skip_if(is.null(fn), "gp_asreml_heritability_se not available")
  mod <- fit_toy_gblup()
  skip_if(is.null(mod), "toy GBLUP did not converge")

  expect_true(is.na(fn(mod, vg_rownames = "no_such_component", ve_rowname = "units!R")))
  expect_true(is.na(fn(mod, vg_rownames = grep("vm\\(GID", rownames(summary(mod)$varcomp), value = TRUE),
                       ve_rowname = "no_such_resid")))
  # A non-model object must not error.
  expect_true(is.na(fn(list(), vg_rownames = "x", ve_rowname = "y")))
})
