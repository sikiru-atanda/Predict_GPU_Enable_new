# 0.20.221: single-environment multi-trait ASReml accepts a factor-analytic genetic structure,
# var_cov_str = "fa1", "fa2", ...: fa(Trait, k):vm(<gen>, Kinv) with an unstructured residual.
# With many traits it has far fewer parameters than us (6 traits, fa1: 12 against 21).

asreml_ready <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) return(FALSE)
  tryCatch(grepl("No error", asreml::asreml.license.status(quiet = TRUE)$rv_text),
           error = function(e) FALSE)
}

# Four traits driven by one common factor plus trait-specific effects, on a full-sib panel.
mt_fa_panel <- function(n_fam = 30, fam_size = 5, p = 400, traits = paste0("T", 1:4)) {
  withr::local_seed(3)
  n <- n_fam * fam_size; nt <- length(traits)
  fam <- rep(seq_len(n_fam), each = fam_size)
  P <- matrix(stats::rbinom(2 * n_fam * p, 1, 0.3), 2 * n_fam, p)
  M <- t(vapply(seq_len(n), function(i) {
    f <- fam[i]
    P[2 * f - 1, ] + ifelse(stats::runif(p) < 0.5, P[2 * f, ], P[2 * f - 1, ])
  }, numeric(p)))
  Z <- scale(M); Z[is.na(Z)] <- 0
  K <- tcrossprod(Z) / p
  ids <- sprintf("L%03d", seq_len(n)); dimnames(K) <- list(ids, ids)
  g <- outer(drop(Z[, 1:40] %*% stats::rnorm(40)), seq(0.9, 0.6, length.out = nt)) +
    Z[, 41:80] %*% matrix(stats::rnorm(40 * nt, sd = 0.3), 40)
  y <- scale(g) * 0.8 + matrix(stats::rnorm(n * nt), n) * 0.7
  ph <- data.frame(GID = ids, y, stringsAsFactors = FALSE)
  names(ph)[-1] <- traits
  for (j in seq_len(nt)) ph[sample(n, 15), traits[j]] <- NA
  ph[(n - 9):n, traits] <- NA
  list(K = K, ph = ph, traits = traits, ids = ids)
}

test_that("the multi-trait ASReml structure accepts fa1, fa2, ... and checks the rank", {
  expect_identical(gp_multitrait_asreml_trait_covariance("FA2"), "fa2")
  expect_identical(gp_multitrait_asreml_trait_covariance("us"), "us")
  expect_identical(gp_multitrait_asreml_fa_rank("fa3"), 3L)
  expect_identical(gp_multitrait_asreml_fa_rank("us"), NA_integer_)
  expect_identical(gp_multitrait_asreml_residual_covariance("fa2"), "us")
  expect_identical(gp_multitrait_asreml_residual_covariance("corgh"), "corgh")
  expect_error(gp_multitrait_asreml_trait_covariance("fa0"), "fa1, fa2")
  expect_error(gp_multitrait_asreml_trait_covariance("ar1"), "fa1, fa2")
  # identifiable ranks (Ledermann bound): n(k + 1) - k(k - 1)/2 <= n(n + 1)/2
  expect_identical(vapply(c(2L, 3L, 4L, 5L, 6L, 8L, 10L, 12L), gp_asreml_fa_max_rank,
                          integer(1)), c(0L, 1L, 1L, 2L, 3L, 4L, 6L, 7L))
  expect_error(gp_asreml_check_fa_rank(1L, 2L, "traits"),
               "not identifiable for 2 traits: it has 4 parameters, more than the 3")
  expect_error(gp_asreml_check_fa_rank(2L, 4L, "traits"), "Use rank 1 or less")
  expect_true(gp_asreml_check_fa_rank(2L, 5L, "traits"))
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel(n_fam = 6, p = 100, traits = c("T1", "T2"))
  expect_error(
    gp_multitrait_asreml_gaussian_model(d$ph, d$traits, "GID", gmatrix = d$K, trait_covariance = "fa1"),
    "not identifiable for 2 traits"
  )
})

test_that("FA loadings and specific variances are read by name into Lambda Lambda' + Psi", {
  nm <- c(paste0("fa(Trait, 2):vm(GID, inv)!", rep(c("A", "B", "C"), 3), "!",
                 rep(c("var", "fa1", "fa2"), each = 3)), "units:Trait!R")
  raw <- data.frame(component = c(0.1, 0.2, 0.3, 1, 0.5, -0.4, 0, 0.6, 0.2, 1), row.names = nm)
  fa <- gp_multitrait_asreml_fa_parameters(raw, c("A", "B", "C"), 1:9, 2L)
  L <- matrix(c(1, 0.5, -0.4, 0, 0.6, 0.2), 3)
  expect_equal(unname(fa$covariance), tcrossprod(L) + diag(c(0.1, 0.2, 0.3)))
  expect_equal(fa$index["B", "fa2"], 8L)
  expect_error(gp_multitrait_asreml_fa_parameters(raw, c("A", "B", "D"), 1:9, 2L), "trait `D`")
})

test_that("FA genetic PEV equals a mixed-model-equations reference with G0 = Lambda Lambda' + Psi", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel()
  long <- gp_multitrait_asreml_long_data(d$ph, d$traits, "GID")
  long$Trait <- factor(long$Trait, levels = d$traits)
  long$GID <- factor(long$GID, levels = d$ids)
  Kinv <- compute_inverse_and_sparse(kernel = d$K, epsilon = 1e-6, inverse = TRUE)
  assign("mt_fa_Kinv", Kinv, envir = .GlobalEnv)
  withr::defer(rm("mt_fa_Kinv", envir = .GlobalEnv))
  fit <- function(...) suppressWarnings(asreml::asreml(
    fixed = Observed_value ~ Trait, random = ~ fa(Trait, 1):vm(GID, mt_fa_Kinv),
    residual = ~ units:us(Trait), data = long, asmv = "Trait",
    na.action = list(x = "include", y = "include"), workspace = 1e8, trace = FALSE, ...
  ))
  # every parameter fixed at its estimate, so predict() and the reference use the same values
  est <- fit()$vparameters
  tab <- fit(start.values = TRUE)$vparameters.table
  tab$Value <- unname(est[as.character(tab$Component)]); tab$Constraint <- "F"
  mod <- fit(G.param = tab, R.param = tab)
  raw <- summary(mod)$varcomp
  fa <- gp_multitrait_asreml_fa_parameters(raw, d$traits, grep("fa(Trait", rownames(raw), fixed = TRUE), 1L)
  R0 <- mt_rel_matrix(mod$vparameters, "units:Trait!Trait", d$traits)
  pred <- asreml_run_with_model_data_context(mod, long, asreml::predict.asreml(
    mod, classify = "Trait:GID", sed = FALSE)$pvals)
  pred <- as.data.frame(pred)
  pred$Trait <- as.character(pred$Trait); pred$GID <- as.character(pred$GID)
  out <- asreml_apply_genetic_reliability(
    pred_df = pred, mod = mod, pheno_data = long, gen_name = "GID", classify = "Trait:GID",
    varG_per_omics = matrix(diag(fa$covariance), 1, dimnames = list(NULL, d$traits)),
    names_in_inv_list = "mt_fa_Kinv", heter_groups = "Trait"
  )
  Kf <- d$K + diag(1e-6, nrow(d$K))
  pev_ref <- mt_rel_reference_pev(d$ph, d$traits, Kf, unname(fa$covariance), R0)
  idx <- match(paste(rep(d$traits, each = nrow(Kf)), rep(d$ids, length(d$traits))),
               paste(out$Trait, out$GID))
  expect_false(anyNA(idx))
  expect_equal(out$Reliability_variance_input[idx], unname(pev_ref), tolerance = 1e-8)
  expect_equal(out$Reliability_reference_variance[idx],
               unname(rep(diag(fa$covariance), each = nrow(Kf)) * rep(diag(Kf), length(d$traits))),
               tolerance = 1e-10)
})

test_that("the multi-trait ASReml route fits fa1 and reports its derived covariances", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel()
  res <- suppressWarnings(gp_multitrait_asreml_gaussian_model(
    d$ph, d$traits, "GID", gmatrix = d$K, trait_covariance = "fa1"
  ))
  expect_true(isTRUE(res$Asreml_model$converge))
  expect_identical(res$asreml_prediction_mode, "predict")
  expect_identical(names(res$Asreml_model$G.param), "fa(Trait, 1):vm(GID, pp_mt_inv_1)")
  raw <- summary(res$Asreml_model)$varcomp
  fa <- gp_multitrait_asreml_fa_parameters(raw, d$traits, grep("fa(Trait", rownames(raw), fixed = TRUE), 1L)
  expect_equal(unname(res$Genetic_covariance_traits), unname(fa$covariance), tolerance = 1e-10)
  vc <- res$variance_components
  g <- vc[vc$Component == "genetic_variance", ]
  expect_equal(g$Components[match(d$traits, g$Trait)], unname(diag(fa$covariance)), tolerance = 1e-10)
  expect_true(all(is.finite(g$Standard_error) & g$Standard_error > 0))
  expect_equal(sum(vc$Component == "genetic_covariance"), choose(length(d$traits), 2))
  expect_equal(sum(vc$Component == "residual_variance"), length(d$traits))
  expect_equal(sum(vc$Component == "residual_covariance"), choose(length(d$traits), 2))
  pl <- res$predicted_values
  expect_true(all(grepl("^genetic_prior_variance: ", pl$Reliability_basis)))
  Kf <- d$K + diag(1e-6, nrow(d$K))
  expect_equal(pl$Reliability_reference_variance,
               unname(diag(fa$covariance)[pl$Trait] * diag(Kf)[as.character(pl$GID)]),
               tolerance = 1e-10)
})

test_that("model_execute runs single-environment multi-trait ASReml with var_cov_str = 'fa1'", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel()
  withr::local_dir(withr::local_tempdir())
  out <- suppressWarnings(PredictProR::model_execute(
    pheno_data = d$ph, gmatrix = d$K, response = d$traits, gen_name = "GID",
    GS_model = "GBLUP", engine = "asreml", response_family = "gaussian",
    multi_trait_asreml = TRUE, var_cov_str = "fa1", message = FALSE
  ))
  expect_true(PredictProR::validate_prediction_output(out$model_results)$valid)
  preds <- out$model_results$predicted_values
  expect_setequal(unique(preds$Trait), d$traits)
  expect_true(all(grepl("^genetic_prior_variance: ", preds$Reliability_basis)))
  expect_true(any(grepl("independent_fa1_trait_covariance",
                        unlist(out$model_results$model_parameters), fixed = TRUE)))
})

test_that("multi-trait ASReml cross-validation runs with fa1 and loses no fold", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel()
  failed <- character()
  cv <- withCallingHandlers(
    gp_multitrait_asreml_gaussian_cv(d$ph, d$traits, "GID", gmatrix = d$K, nfolds = 3L,
                                     trait_covariance = "fa1"),
    warning = function(w) {
      if (grepl("fit failed", conditionMessage(w))) failed <<- c(failed, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(failed, 0L)
  pred <- cv$cv_results[[1L]]$ypred_cv_Reps_all
  expect_setequal(unique(pred$fold), 1:3)
  expect_setequal(unique(pred$Trait), d$traits)
  expect_true(all(is.finite(pred$Predicted_value)))
})

test_that("single-trait MET ASReml refuses an FA rank the environments cannot identify", {
  skip_if_not(asreml_ready(), "ASReml-R not installed or not licensed")
  d <- mt_fa_panel(n_fam = 6, p = 100, traits = "T1")
  ph <- do.call(rbind, lapply(c("E1", "E2", "E3"), function(e) {
    data.frame(GID = d$ph$GID, Env = e, Yield = d$ph$T1 + rnorm(nrow(d$ph), sd = 0.1))
  }))
  withr::local_dir(withr::local_tempdir())
  run <- function(vcs) PredictProR::model_execute(
    pheno_data = ph, gmatrix = d$K, response = "Yield", gen_name = "GID", GS_model = "GBLUP",
    engine = "asreml", fixed = ~Env, random = ~GID + GID:Env, heter_groups = "Env",
    heter_resid = TRUE, var_cov_str = vcs, message = FALSE
  )
  # fa2 over 3 environments: 8 parameters for 6 covariance entries; stopped before any fit
  expect_error(run("fa2"), "not identifiable for 3 environments: it has 8 parameters, more than the 6")
})

test_that("multi-trait ASReml retries a fit once at 4gb when it runs out of workspace", {
  # 0.20.221: at about 1,000 lines and 6 traits the default 1e8-word (800 MB) workspace is
  # too small and every multi-trait fit failed.
  calls <- list()
  fake <- function(...) {
    a <- list(...); calls[[length(calls) + 1L]] <<- a$workspace
    if (length(calls) == 1L) stop("Error   : Insufficient workspace; current workspace is approx 800.000 Mb")
    structure(list(workspace = a$workspace), class = "fake_fit")
  }
  fit <- gp_multitrait_asreml_fit(fake, list(workspace = 1e8), "Multi-trait ASReml-R fit failed: ")
  expect_equal(calls[[1L]], 1e8); expect_identical(calls[[2L]], "4gb")
  expect_identical(fit$workspace, "4gb")
  # other errors: no retry, the route's message
  calls <- list()
  other <- function(...) { calls[[length(calls) + 1L]] <<- 1; stop("singular") }
  expect_error(gp_multitrait_asreml_fit(other, list(workspace = 1e8), "Multi-trait ASReml-R fit failed: "),
               "Multi-trait ASReml-R fit failed: singular")
  expect_length(calls, 1L)
  # already at 4gb: no retry
  expect_null(gp_multitrait_asreml_retry_workspace("4gb"))
  expect_null(gp_multitrait_asreml_retry_workspace(6e8))
  expect_identical(gp_multitrait_asreml_retry_workspace(1e8), "4gb")
  expect_identical(gp_multitrait_asreml_retry_workspace(NULL), "4gb")
  # retry also fails: both sizes in the message
  always <- function(...) stop("Insufficient workspace available when reordering matrices")
  expect_error(gp_multitrait_asreml_fit(always, list(workspace = 1e8), "Grouped MT-MET ASReml GBLUP fit failed: "),
               "at the supplied workspace and at 4gb")
  # ASReml's other wording for the same shortage (fa2, 6 traits, 1,000 lines)
  expect_true(gp_multitrait_asreml_workspace_error("Error   : Failed ordering equations. Increase workspace"))
  expect_false(gp_multitrait_asreml_workspace_error("Error   : singular"))
  # the prediction retry still targets 1gb
  expect_identical(gp_multitrait_asreml_retry_pworkspace(1e6), "1gb")
  expect_null(gp_multitrait_asreml_retry_pworkspace("2gb"))
})

test_that("multi-trait ASReml defaults to fa2 from six traits and us below", {
  # 0.20.221: capacity study, 300 lines with 10% trait gaps: us converged conveniently up to 4
  # traits, fa2 up to 10.
  expect_identical(gp_multitrait_asreml_trait_covariance(NULL, 5L), "us")
  expect_identical(gp_multitrait_asreml_trait_covariance(NULL, 6L), "fa2")
  expect_identical(gp_multitrait_asreml_trait_covariance(NULL, 12L), "fa2")
  expect_identical(gp_multitrait_asreml_trait_covariance(NULL), "us")
  # an explicit choice is never replaced
  expect_identical(gp_multitrait_asreml_trait_covariance("us", 8L), "us")
  expect_identical(gp_multitrait_asreml_trait_covariance("fa1", 8L), "fa1")
})

test_that("model_execute passes the six-trait fa2 default to the multi-trait ASReml model", {
  seen <- NULL
  local_mocked_bindings(gp_multitrait_asreml_gaussian_model = function(..., trait_covariance) {
    seen <<- trait_covariance
    stop("stop after routing")
  })
  d <- mt_fa_panel(n_fam = 6, p = 100, traits = paste0("T", 1:6))
  withr::local_dir(withr::local_tempdir())
  run <- function(ph, traits) try(suppressWarnings(PredictProR::model_execute(
    pheno_data = ph, gmatrix = d$K, response = traits, gen_name = "GID", GS_model = "GBLUP",
    engine = "asreml", response_family = "gaussian", multi_trait_asreml = TRUE, message = FALSE
  )), silent = TRUE)
  run(d$ph, d$traits)
  expect_identical(seen, "fa2")
  seen <- NULL
  run(d$ph[, c("GID", "T1", "T2", "T3")], c("T1", "T2", "T3"))
  expect_identical(seen, "us")
})
