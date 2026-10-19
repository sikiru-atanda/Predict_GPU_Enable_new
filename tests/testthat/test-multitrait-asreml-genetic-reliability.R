# 0.20.221: multi-trait ASReml reliability is genetic, as on the single-trait route and the GP
# engine: 1 - Var(ghat - g) / sum_k sigma2_g(trait, k) * K_k,ii. It was 1 - total SE^2 / sigma2_g
# (the total SE includes the trait mean), optimistic by about 0.06 in the validation study.

asreml_ready <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) return(FALSE)
  tryCatch(grepl("No error", asreml::asreml.license.status(quiet = TRUE)$rv_text),
           error = function(e) FALSE)
}

# Three correlated traits on a marker kernel, some records missing per trait.
mt_rel_panel <- function(n = 120, p = 400, traits = c("T1", "T2", "T3")) {
  withr::local_seed(21)
  M <- matrix(stats::rbinom(n * p, 2, 0.3), n, p)
  Z <- scale(M); K <- tcrossprod(Z) / p
  ids <- sprintf("L%03d", seq_len(n)); dimnames(K) <- list(ids, ids)
  G0 <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3)
  B <- matrix(stats::rnorm(20 * 3), 20, 3) %*% chol(G0)
  g <- scale(Z[, 1:20] %*% B)
  y <- g + matrix(stats::rnorm(3 * n), n) + matrix(c(10, 5, 2), n, 3, byrow = TRUE)
  ph <- data.frame(GID = ids, y, stringsAsFactors = FALSE)
  names(ph)[-1] <- traits
  ph$T1[111:120] <- NA; ph$T2[101:120] <- NA; ph$T3[c(1:5, 116:120)] <- NA
  list(K = K, ph = ph, traits = traits, ids = ids)
}

test_that("multi-trait genetic variances become a kernels x traits matrix", {
  specs <- data.frame(Kernel = c("A", "D"), inverse_name = c("i1", "i2"))
  us <- function(inv) paste0("Trait:vm(GID, ", inv, ")!Trait_", c("T1:T1", "T2:T1", "T2:T2"))
  raw <- data.frame(component = c(1, 0.5, 2, 3, -0.2, 4, 9, 0.1, 8),
                    row.names = c(us("i1"), us("i2"), paste0("units:Trait!Trait_", c("T1:T1", "T2:T1", "T2:T2"))))
  m <- gp_multitrait_asreml_genetic_variance_matrix(raw, c("T1", "T2"), specs, "us")
  expect_equal(m, matrix(c(1, 3, 2, 4), 2, dimnames = list(c("A", "D"), c("T1", "T2"))))
  expect_equal(gp_multitrait_asreml_genetic_variance_matrix(raw, c("T1", "T2"), specs[1, ], "us")[1, ],
               c(T1 = 1, T2 = 2))
  fa <- data.frame(component = c(0.1, 0.2, 1, 0.5),
                   row.names = paste0("fa(Trait, 1):vm(GID, i1)!", c("T1!var", "T2!var", "T1!fa1", "T2!fa1")))
  expect_equal(gp_multitrait_asreml_genetic_variance_matrix(fa, c("T1", "T2"), specs[1, ], "fa1")[1, ],
               c(T1 = 1.1, T2 = 0.45))
  expect_error(gp_multitrait_asreml_genetic_variance_matrix(raw[-3, , drop = FALSE], c("T1", "T2"), specs, "us"),
               "Expected exactly one")
})

test_that("three-trait ASReml genetic PEV equals a mixed-model-equations reference", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_rel_panel()
  long <- gp_multitrait_asreml_long_data(d$ph, d$traits, "GID")
  long$Trait <- factor(long$Trait, levels = d$traits)
  long$GID <- factor(long$GID, levels = d$ids)
  Kinv <- compute_inverse_and_sparse(kernel = d$K, epsilon = 1e-6, inverse = TRUE)
  assign("mt_rel_Kinv", Kinv, envir = .GlobalEnv)
  withr::defer(rm("mt_rel_Kinv", envir = .GlobalEnv))
  fit <- function(...) suppressWarnings(asreml::asreml(
    fixed = Observed_value ~ Trait, random = ~ us(Trait):vm(GID, mt_rel_Kinv),
    residual = ~ units:us(Trait), data = long, asmv = "Trait",
    na.action = list(x = "include", y = "include"), workspace = 1e8, trace = FALSE, ...
  ))
  # predict() works at the last iteration's parameters, and this fit is near the boundary
  # (they still move about 1% per update): refit with every parameter fixed at the
  # estimates, so the reference built from the reported parameters describes the same MME
  est <- fit()$vparameters
  tab <- fit(start.values = TRUE)$vparameters.table
  tab$Value <- unname(est[as.character(tab$Component)]); tab$Constraint <- "F"
  mod <- fit(G.param = tab, R.param = tab)
  vp <- mod$vparameters
  expect_equal(vp, est)
  G0 <- mt_rel_matrix(vp, "vm(GID", d$traits); R0 <- mt_rel_matrix(vp, "units:Trait!Trait", d$traits)
  pred <- asreml_run_with_model_data_context(mod, long, asreml::predict.asreml(
    mod, classify = "Trait:GID", sed = FALSE)$pvals)
  pred <- as.data.frame(pred)
  pred$Trait <- as.character(pred$Trait); pred$GID <- as.character(pred$GID)
  out <- asreml_apply_genetic_reliability(
    pred_df = pred, mod = mod, pheno_data = long, gen_name = "GID", classify = "Trait:GID",
    varG_per_omics = matrix(diag(G0), 1, dimnames = list(NULL, d$traits)),
    names_in_inv_list = "mt_rel_Kinv", heter_groups = "Trait"
  )
  Kf <- d$K + diag(1e-6, nrow(d$K))
  pev_ref <- mt_rel_reference_pev(d$ph, d$traits, Kf, G0, R0)
  prior_ref <- rep(diag(G0), each = nrow(Kf)) * rep(diag(Kf), length(d$traits))
  idx <- match(paste(rep(d$traits, each = nrow(Kf)), rep(d$ids, length(d$traits))),
               paste(out$Trait, out$GID))
  expect_false(anyNA(idx))
  expect_equal(out$Reliability_variance_input[idx], unname(pev_ref), tolerance = 1e-8)
  expect_equal(out$Reliability_reference_variance[idx], unname(prior_ref), tolerance = 1e-10)
  expect_equal(out$Reliability[idx], unname(pmax(0, pmin(1, 1 - pev_ref / prior_ref))),
               tolerance = 1e-8)
  expect_match(out$Reliability_basis[1], "^genetic_prior_variance: ")
  # the genetic PEV leaves out the trait means, so it is below the total SE^2
  expect_true(all(out$Reliability_variance_input < out$std.error^2))
})

test_that("the multi-trait ASReml route reports genetic reliability for three traits", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_rel_panel()
  res <- suppressWarnings(gp_multitrait_asreml_gaussian_model(
    d$ph, d$traits, "GID", gmatrix = d$K, trait_covariance = "us"
  ))
  expect_identical(res$asreml_prediction_mode, "predict")
  pl <- res$predicted_values
  expect_true(all(grepl("^genetic_prior_variance: ", pl$Reliability_basis)))
  expect_true(all(is.finite(pl$Reliability)))
  expect_equal(pl$Reliability,
               pmax(0, pmin(1, 1 - pl$Reliability_variance_input / pl$Reliability_reference_variance)))
  vc <- res$variance_components
  s2g <- vapply(d$traits, function(t) {
    vc$Components[vc$Component == "genetic_variance" & vc$Trait %in% t]
  }, numeric(1))
  Kf <- d$K + diag(1e-6, nrow(d$K))
  expect_equal(pl$Reliability_reference_variance,
               unname(s2g[pl$Trait] * diag(Kf)[as.character(pl$GID)]), tolerance = 1e-10)
  # against the reference at the reported parameters: predict() uses the last iteration's
  # parameters, so the route's default convergence leaves a small gap (about 2% of PEV and
  # 0.01 in reliability for three traits; 0.006 for six); the exact check is the test above
  vp <- res$Asreml_model$vparameters
  G0 <- mt_rel_matrix(vp, "vm(GID", d$traits); R0 <- mt_rel_matrix(vp, "units:Trait!Trait", d$traits)
  pev_ref <- mt_rel_reference_pev(d$ph, d$traits, Kf, G0, R0)
  idx <- match(paste(rep(d$traits, each = nrow(Kf)), rep(d$ids, length(d$traits))),
               paste(pl$Trait, pl$GID))
  rel_ref <- 1 - pev_ref / (rep(diag(G0), each = nrow(Kf)) * rep(diag(Kf), length(d$traits)))
  expect_equal(pl$Reliability_variance_input[idx], unname(pev_ref), tolerance = 0.05)
  expect_lt(max(abs(pl$Reliability[idx] - rel_ref)), 0.02)
})
