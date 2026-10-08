# Bayesian MET (GBLUP_BRR / RKHS) on unbalanced data. The G2F benchmark
# (2,311 hybrids, 11,270 records in 10 environments) lost the whole
# prediction twice: the record-level kernel was built with an alphabetical
# incidence matrix that did not fit a kernel holding genotypes without
# records, and the output table paired one-per-genotype names with
# one-per-record effects ("differing number of rows: 2311, 11270").

test_that("record-level kernel matches records to kernel rows by genotype ID", {
  ids <- c("G3", "G1", "G4", "G2")                       # not alphabetical; G4 has no record
  K <- matrix(c(4, 1, 0, 2,
                1, 5, 1, 0,
                0, 1, 6, 1,
                2, 0, 1, 7), 4, 4, dimnames = list(ids, ids))
  rec <- c("G1", "G2", "G1", "G3", "G2")
  Z <- outer(rec, ids, "==") * 1
  expected <- Z %*% K %*% t(Z)
  out <- PredictProR:::gp_bayes_record_kernel(rec, K)
  expect_equal(unname(out), unname(expected))
  expect_identical(dim(out), c(5L, 5L))
  expect_error(PredictProR:::gp_bayes_record_kernel(c("G1", "G9"), K), "not in the kernel")
  expect_error(PredictProR:::gp_bayes_record_kernel(rec, unname(K)), "no genotype row names")
})

# The records x records eigen-decomposition (BRR design, BGLR RKHS setup) took
# about an hour at 11,270 G2F records. Its non-zero eigenpairs come from the
# genotype kernel scaled by record counts (per environment block for G x E).
test_that("record-kernel eigenpairs reproduce Z K Z' and its environment-blocked form", {
  set.seed(4)
  ids <- c("G5", "G2", "G9", "G1", "G7", "G3")               # G3 has no records
  M <- matrix(stats::rnorm(6 * 20), 6, 20)
  K <- tcrossprod(M) / 20 + diag(0.1, 6); dimnames(K) <- list(ids, ids)
  rec <- c("G1", "G2", "G2", "G5", "G9", "G1", "G7", "G9", "G9", "G2")
  env <- c("E1", "E1", "E2", "E1", "E2", "E2", "E1", "E1", "E2", "E2")
  K1 <- PredictProR:::gp_bayes_record_kernel(rec, K)
  e1 <- PredictProR:::gp_bayes_record_eigen(rec, K)
  expect_equal(unname(e1$vectors %*% diag(e1$values) %*% t(e1$vectors)), unname(K1), tolerance = 1e-10)
  expect_equal(crossprod(e1$vectors), diag(length(e1$values)), tolerance = 1e-10)
  expect_equal(e1$values, sort(eigen(K1, symmetric = TRUE)$values, decreasing = TRUE)[seq_along(e1$values)],
               tolerance = 1e-10)
  K2 <- K1 * outer(env, env, "==")
  e2 <- PredictProR:::gp_bayes_record_eigen(rec, K, groups = env)
  expect_equal(unname(e2$vectors %*% diag(e2$values) %*% t(e2$vectors)), unname(K2), tolerance = 1e-10)
  expect_equal(crossprod(e2$vectors), diag(length(e2$values)), tolerance = 1e-10)
  X <- PredictProR:::gp_bayes_brr_design_from_eigen(e2)
  expect_equal(unname(tcrossprod(X)), unname(K2), tolerance = 1e-8)
  expect_error(PredictProR:::gp_bayes_record_eigen(c("G1", "G8"), K), "not in the kernel")
})

test_that("precomputed RKHS eigenpairs are dropped for weighted fits", {
  eta <- list(list(K = diag(2), V = diag(2), d = c(1, 1), model = "RKHS"), list(X = diag(2), model = "BRR"))
  out <- PredictProR:::gp_bayes_drop_precomputed_eigen(eta)
  expect_null(out[[1]]$V); expect_null(out[[1]]$d); expect_false(is.null(out[[1]]$K))
  expect_identical(out[[2]], eta[[2]])
})

test_that("RKHS coefficient table labels record-level effects with each record's genotype and group", {
  mod <- list(model = list(ETA = list(list(u = c(0.1, 0.2, 0.3), SD.u = c(1, 1, 1)))))
  out <- PredictProR:::cal_coeff_ebv_pev_rel_RHKS_glub(
    mod = mod, gen_name = "GID", var_u = 1, var_E = 1,
    hetero = c("E2", "E1", "E2"), heter_groups = "Env", gid_name = c("A", "B", "B"), eta_index = 1L
  )
  ebv <- out$Estimated_breeding_value
  expect_identical(ebv$GID, c("A", "B", "B"))
  expect_identical(ebv$Env, c("E2", "E1", "E2"))
  expect_equal(ebv$Estimated_breeding_value, c(0.1, 0.2, 0.3))
})

test_that("unbalanced MET GBLUP_BRR and RKHS predict the full grid with aligned labels", {
  skip_on_cran()
  skip_if_not_installed("BGLR")
  set.seed(21)
  n <- 30; p <- 150
  ids <- sprintf("G%02d", seq_len(n))
  X <- matrix(stats::rbinom(n * p, 2, 0.4), n, p, dimnames = list(ids, sprintf("m%03d", seq_len(p))))
  g <- as.numeric(scale(X[, 1:10] %*% stats::rnorm(10)))
  ph <- expand.grid(Env = c("E1", "E2", "E3"), GID = ids, stringsAsFactors = FALSE)[, c("GID", "Env")]
  ph <- ph[sample(nrow(ph), 70), ]                       # unbalanced, shuffled; not a multiple of n
  ph$Yield <- g[match(ph$GID, ids)] + stats::rnorm(nrow(ph), sd = 0.5)
  test <- unique(ph$GID)[1:6]
  ph$Yield[ph$GID %in% test] <- NA
  withr::local_dir(withr::local_tempdir())
  for (m in c("GBLUP_BRR", "RKHS")) {
    res <- suppressWarnings(model_execute(
      pheno_data = ph, geno_data = X, gen_name = "GID", response = "Yield",
      fixed = ~1, random = ~ GID + GID:Env, heter_groups = "Env", GS_model = m,
      gmatrix_method = "VanRaden", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
      nIter = 300L, burnIn = 100L, thin = 5L, random_state = 7L, parallel_mode = "sequential",
      system_database = FALSE, message = FALSE, verbose = FALSE
    ))
    find_pred <- function(x, d = 0) {
      if (is.data.frame(x) && all(c("Predicted_value", "GID") %in% names(x))) return(x)
      if (is.list(x) && d < 6) for (e in x) { r <- find_pred(e, d + 1); if (!is.null(r)) return(r) }
      NULL
    }
    pv <- find_pred(res)
    expect_false(is.null(pv), info = m)
    # full line x environment grid (0.20.179), covering every input record
    expect_identical(nrow(pv), length(unique(ph$GID)) * length(unique(ph$Env)), info = m)
    key_out <- paste(pv$GID, pv$Env)
    expect_true(all(paste(ph$GID, ph$Env) %in% key_out), info = m)
    obs <- ph$Yield[match(key_out, paste(ph$GID, ph$Env))]
    train <- pv$Train_Test_Label == "Train"
    expect_equal(pv$Observed_value[train], obs[train], info = m)
  }
})
