# Three pinning behaviours after the MET ML/DL input-contract hardening:
#   1. Auto-promote a Yang GRM build when met_ml_dl=TRUE and the user
#      supplied geno_data but no kernel and no gmatrix_method.
#   2. Surface dropped tasks (warning + dropped_tasks attribute) when CV
#      tasks come back with all-NA predictions, instead of silently
#      "completing with 0 valid results".
#   3. Multi-kernel acceptance: kernel_list with arbitrary extra kernels
#      flows through the guardrail unchanged (no clobbering by auto-promote).

# ---- gp_auto_promote_met_ml_dl_kernel_build ----------------------------------

test_that("auto-promotes gmatrix_method='Yang' for MET ML/DL with only geno_data", {
  ctx <- list(
    met_ml_dl = TRUE,
    geno_data = matrix(0, nrow = 4, ncol = 8,
                       dimnames = list(c("g1","g2","g3","g4"), paste0("snp", 1:8))),
    gmatrix = NULL, gkernel = NULL, kernel_list = NULL,
    omic1_kernel = NULL, omic2_kernel = NULL, omic3_kernel = NULL,
    gmatrix_method = NULL
  )
  expect_message(
    gp_auto_promote_met_ml_dl_kernel_build(ctx),
    "MET ML/DL"
  )
  out <- suppressMessages(gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_identical(out$gmatrix_method, "Yang")
  expect_true(length(out$model_input_guardrail_notes) >= 1L)
  expect_true(any(grepl("auto-set gmatrix_method",
                        out$model_input_guardrail_notes, fixed = TRUE)))
})

test_that("auto-promote uses the ploidy-aware VanRaden GRM for polyploid dosage", {
  G <- matrix(c(0, 1, 2, 3, 4, 2, 1, 3), nrow = 4, ncol = 8,
              dimnames = list(c("g1","g2","g3","g4"), paste0("snp", 1:8)))
  ctx <- list(met_ml_dl = TRUE, geno_data = G, gmatrix_method = NULL, ploidy = 4L)
  out <- suppressMessages(gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_identical(out$gmatrix_method, "VanRaden")
  # ploidy carried by the matrix attribute, ploidy = "auto"
  attr(G, "ploidy") <- 6L
  ctx <- list(met_ml_dl = TRUE, geno_data = G, gmatrix_method = NULL, ploidy = "auto")
  expect_message(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx), "VanRaden GRM")
  expect_identical(out$gmatrix_method, "VanRaden")
  # diploid stays Yang
  ctx$ploidy <- 2L
  expect_identical(suppressMessages(gp_auto_promote_met_ml_dl_kernel_build(ctx))$gmatrix_method, "Yang")
})

test_that("does not auto-promote when user already supplied gmatrix_method", {
  ctx <- list(
    met_ml_dl = TRUE,
    geno_data = matrix(0, nrow = 4, ncol = 8,
                       dimnames = list(c("g1","g2","g3","g4"), paste0("snp", 1:8))),
    gmatrix_method = "VanRaden"
  )
  expect_silent(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_identical(out$gmatrix_method, "VanRaden")  # user choice preserved
})

test_that("does not auto-promote when a kernel is already supplied", {
  K <- diag(4); rownames(K) <- colnames(K) <- paste0("g", 1:4)
  for (slot in c("gmatrix", "gkernel", "omic1_kernel", "omic2_kernel", "omic3_kernel")) {
    ctx <- list(met_ml_dl = TRUE, geno_data = matrix(0, 4, 8))
    ctx[[slot]] <- K
    expect_silent(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx))
    expect_null(out$gmatrix_method)
  }
})

test_that("does not auto-promote when kernel_list carries kernels", {
  K <- diag(4); rownames(K) <- colnames(K) <- paste0("g", 1:4)
  ctx <- list(
    met_ml_dl = TRUE,
    geno_data = matrix(0, 4, 8),
    kernel_list = list(K_aux = K)
  )
  expect_silent(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_null(out$gmatrix_method)
})

test_that("does not auto-promote when met_ml_dl is FALSE", {
  ctx <- list(met_ml_dl = FALSE, geno_data = matrix(0, 4, 8))
  expect_silent(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_null(out$gmatrix_method)
})

test_that("does not auto-promote when there is no geno_data to build from", {
  ctx <- list(met_ml_dl = TRUE,
              geno_data = NULL, train_geno_data = NULL, test_geno_data = NULL)
  expect_silent(out <- gp_auto_promote_met_ml_dl_kernel_build(ctx))
  expect_null(out$gmatrix_method)
})

test_that("gp_apply_model_input_guardrails plumbs the auto-promote through", {
  ctx <- list(
    met_ml_dl = TRUE,
    geno_data = matrix(0, 4, 8,
                       dimnames = list(c("g1","g2","g3","g4"), paste0("snp", 1:8))),
    pheno_data = data.frame(GID = c("g1","g1","g2","g2","g3","g3","g4","g4"),
                            Env = rep(c("E1","E2"), 4),
                            Yield = rnorm(8)),
    gen_name = "GID", heter_groups = "Env", response = "Yield",
    response_family = "gaussian"
  )
  out <- suppressMessages(gp_apply_model_input_guardrails(ctx))
  expect_identical(out$gmatrix_method, "Yang")
})

test_that("MET auto-promotion defers to a dedicated hybrid workflow", {
  ctx <- list(
    met_ml_dl = FALSE,
    hybrid_ml = TRUE,
    pheno_data = data.frame(
      HybridID = rep(c("h1", "h2"), each = 2L),
      Env = rep(c("E1", "E2"), times = 2L),
      Yield = c(1, 2, 3, 4)
    ),
    gen_name = "HybridID",
    heter_groups = "Env",
    response = "Yield",
    GS_model = "RandomForest",
    GS_model_cv = character()
  )

  out <- suppressMessages(gp_auto_promote_met_ml_dl_flag(ctx))

  expect_false(isTRUE(out$met_ml_dl))
  expect_true(any(grepl(
    "DEFERRED to a dedicated multi-trait or hybrid workflow",
    out$model_input_guardrail_notes,
    fixed = TRUE
  )))
})

test_that("MET auto-promotion rejects unsupported mixed model lists before dispatch", {
  ctx <- list(
    met_ml_dl = FALSE,
    pheno_data = data.frame(
      GID = rep(c("g1", "g2"), each = 2L),
      Env = rep(c("E1", "E2"), times = 2L),
      Yield = c(1, 2, 3, 4)
    ),
    gen_name = "GID",
    heter_groups = "Env",
    response = "Yield",
    GS_model = c("Xgboost", "Lasso")
  )

  expect_error(
    gp_auto_promote_met_ml_dl_flag(ctx),
    "Unsupported single-trait multi-environment model.*Lasso.*stopped before fitting"
  )
})

# ---- gp_filter_crossval_results: dropped tasks surfaced ----------------------

test_that("gp_filter_crossval_results warns and records collapsed tasks", {
  # Two tasks. One has all-NA predictions (collapsed); the other is fine.
  results <- list(
    list(
      trait = "Yield", model = "Xgboost", model_canonical = "Xgboost", rep = 1L,
      eval_metrics_reps = data.frame(accuracy = NA_real_),
      ypred_cv_Reps_all = data.frame(yhat = rep(NA_real_, 5L),
                                     cv_role = c("train","train","test","test","test"))
    ),
    list(
      trait = "Yield", model = "RKHS", model_canonical = "RKHS", rep = 1L,
      eval_metrics_reps = data.frame(accuracy = 0.45),
      ypred_cv_Reps_all = data.frame(yhat = c(NA, NA, 1.0, 2.0, 1.5),
                                     cv_role = c("train","train","test","test","test"))
    )
  )
  expect_warning(
    out <- gp_filter_crossval_results(results),
    "all folds produced NA predictions"
  )
  expect_length(out, 1L)
  dropped <- attr(out, "dropped_tasks")
  expect_true(is.list(dropped))
  expect_true("Xgboost/Yield rep=1" %in% dropped$collapsed_all_na_predictions)
  expect_length(dropped$nonfinite_metric, 0L)
})

test_that("gp_filter_crossval_results warns separately for nonfinite metrics", {
  results <- list(
    list(
      trait = "Yield", model = "RandomForest", model_canonical = "RandomForest", rep = 1L,
      eval_metrics_reps = data.frame(accuracy = NA_real_),
      ypred_cv_Reps_all = data.frame(yhat = c(NA, NA, 1.0, 2.0),
                                     cv_role = c("train","train","test","test"))
    )
  )
  expect_warning(
    out <- gp_filter_crossval_results(results),
    "non-finite"
  )
  expect_length(out, 0L)
  dropped <- attr(out, "dropped_tasks")
  expect_true("RandomForest/Yield rep=1" %in% dropped$nonfinite_metric)
})

test_that("gp_filter_crossval_results retains valid CV when one class metric is undefined", {
  results <- list(
    list(
      trait = "Class", model = "BRR", model_canonical = "BRR", rep = 1L,
      eval_metrics_reps = data.frame(
        Rep = 1L, accuracy = 0.5, balanced_accuracy = 0.5,
        mcc = NA_real_, log_loss = 0.8
      ),
      ypred_cv_Reps_all = data.frame(
        yhat = c(0.8, 0.7, 0.9, 0.6), cv_role = rep("test", 4L)
      )
    )
  )
  expect_silent(out <- gp_filter_crossval_results(results))
  expect_length(out, 1L)
  expect_true(is.na(out[[1L]]$eval_metrics_reps$mcc[[1L]]))
  expect_length(attr(out, "dropped_tasks")$nonfinite_metric, 0L)
})

test_that("gp_filter_crossval_results is silent when all tasks pass", {
  results <- list(
    list(
      trait = "Yield", model = "RKHS", model_canonical = "RKHS", rep = 1L,
      eval_metrics_reps = data.frame(accuracy = 0.5),
      ypred_cv_Reps_all = data.frame(yhat = c(NA, NA, 1.0, 2.0, 1.5),
                                     cv_role = c("train","train","test","test","test"))
    )
  )
  expect_silent(out <- gp_filter_crossval_results(results))
  expect_length(out, 1L)
  dropped <- attr(out, "dropped_tasks")
  expect_length(dropped$collapsed_all_na_predictions, 0L)
  expect_length(dropped$nonfinite_metric, 0L)
})
