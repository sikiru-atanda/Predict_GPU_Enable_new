# Feature scoring rigor (0.20.218, spec Part B): records-only fold scoring, fixed-effect and
# environment adjustment, AIR RandomForest importance, ties at the k cut-off, the selection remark
# and the leakage-free CV default.

test_that("rankings use raw scores even when no score is positive", {
  meta <- gp_feature_metadata_table(
    trait = "Y", predictors = c("b", "a", "c"),
    source_block = c(b = "merged", a = "merged", c = "merged"),
    raw_score = c(-0.1, -0.5, -0.2), scoring_model = "RandomForest", seed = 1L, n_train = 10L
  )
  expect_equal(meta$predictor[order(meta$rank)], c("b", "c", "a"))
  expect_equal(meta$normalized_score[meta$predictor == "b"], -0.2)   # -0.1 / max|raw| = -0.1 / 0.5
})

test_that("non-negative scores keep their previous ranks and normalized scores", {
  meta <- gp_feature_metadata_table(
    trait = "Y", predictors = c("m1", "m2", "m3"),
    source_block = c(m1 = "merged", m2 = "merged", m3 = "merged"),
    raw_score = c(0.2, 0.8, 0.4), scoring_model = "Ridge_Regression", seed = 1L, n_train = 10L
  )
  expect_equal(meta$predictor[order(meta$rank)], c("m2", "m3", "m1"))
  expect_equal(meta$normalized_score[meta$predictor == "m3"], 0.5)
})

.tie_meta <- function() {
  data.frame(
    trait = "Y", predictor = c("a", "b", "c", "d"), source_block = "merged", rank = 1:4,
    raw_score = c(3, 2, 2, 1), normalized_score = c(1, 2 / 3, 2 / 3, 1 / 3),
    stringsAsFactors = FALSE
  )
}

test_that("a tie at the top-k cut-off warns once per run and is recorded", {
  gp_feature_reset_tie_warning()
  meta <- .tie_meta()
  expect_warning(sel <- gp_feature_selected_metadata(meta, "Y", 2L), "share the score at the top-2 cut-off")
  expect_equal(sel$predictor, c("a", "b"))
  expect_true(isTRUE(attr(sel, "tied_at_cutoff")))
  expect_no_warning(sel2 <- gp_feature_selected_metadata(meta, "Y", 2L))
  expect_true(isTRUE(attr(sel2, "tied_at_cutoff")))
  x <- data.frame(a = 1:3, b = 1:3, c = 1:3, d = 1:3)
  expect_no_warning(top <- feature_select_top_k(x, meta, trait = "Y", k = 2))
  expect_true(isTRUE(attr(top, "tied_at_cutoff")))
})

test_that("no tie, k = all, or k at least the number of predictors gives no tie warning", {
  gp_feature_reset_tie_warning()
  meta <- .tie_meta()
  expect_no_warning(sel <- gp_feature_selected_metadata(meta, "Y", 3L))
  expect_false(isTRUE(attr(sel, "tied_at_cutoff")))
  expect_no_warning(gp_feature_selected_metadata(meta, "Y", "all"))
  expect_no_warning(gp_feature_selected_metadata(meta, "Y", 10L))
})

# 30 lines in two environments; environment E2's records are the held-out part of the fold (as in
# CV0). Lines with training records in E1 also have held-out records in E2.
.leak_data <- function() {
  set.seed(2)
  X <- matrix(rbinom(30 * 20, 2, 0.4), 30,
              dimnames = list(sprintf("g%02d", 1:30), sprintf("m%02d", 1:20)))
  ph <- data.frame(GID = rep(rownames(X), 2), Env = rep(c("E1", "E2"), each = 30),
                   stringsAsFactors = FALSE)
  ph$Yield <- X[ph$GID, 1] + rnorm(60, sd = 0.3)
  list(X = X, ph = ph, test = which(ph$Env == "E2"))
}

test_that("fold-internal scoring uses only the fold's training records", {
  d <- .leak_data()
  view <- function(p) gp_feature_fold_view(
    predictor_data = d$X, pheno_data = p, response = "Yield", gen_name = "GID",
    test_rows = d$test, k = 5, scoring_model = "Ridge_Regression", seed = 1L
  )
  v1 <- view(d$ph)
  ph2 <- d$ph
  ph2$Yield[d$test] <- d$X[ph2$GID[d$test], 7] * 5 + rnorm(30, sd = 0.3)   # change only held-out rows
  v2 <- view(ph2)
  expect_equal(unique(v1$metadata$n_train), 30L)
  expect_identical(v1$selected_predictors, v2$selected_predictors)
  expect_identical(v1$metadata$raw_score, v2$metadata$raw_score)
})

test_that("the specialised routes' fold view uses only training records too", {
  d <- .leak_data()
  v <- gp_feature_trait_view(
    predictor_data = d$X, pheno_data = d$ph, response = "Yield", gen_name = "GID", k = 5,
    selection_mode = "fold_internal", test_rows = d$test, scoring_model = "Ridge_Regression", seed = 1L
  )
  expect_equal(unique(v$metadata$n_train), 30L)
})

# Unbalanced two-environment trial: E1 has lines 1-30, E2 lines 11-40, and E2's mean is 5 higher.
# Marker m05 is 0 on the lines only in E1 and 2 on the lines only in E2, so it tracks the
# environment; m01 is the real QTL.
.met_confounded <- function() {
  set.seed(31)
  lines <- sprintf("L%02d", 1:40)
  X <- matrix(rbinom(40 * 20, 2, 0.4), 40, dimnames = list(lines, sprintf("m%02d", 1:20)))
  X[1:10, "m05"] <- 0
  X[31:40, "m05"] <- 2
  ph <- rbind(data.frame(GID = lines[1:30], Env = "E1"), data.frame(GID = lines[11:40], Env = "E2"))
  ph$Yield <- ifelse(ph$Env == "E2", 5, 0) + 1.5 * X[ph$GID, "m01"] + rnorm(nrow(ph), sd = 0.5)
  list(X = X, ph = ph)
}
.rank_of <- function(meta, p) meta$rank[meta$predictor == p]

test_that("environment adjustment removes a marker that only tracks the environment (Ridge)", {
  d <- .met_confounded()
  raw <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression")
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  heter_groups = "Env")
  expect_lte(.rank_of(raw, "m05"), 2L)
  expect_equal(.rank_of(adj, "m01"), 1L)
  expect_gt(.rank_of(adj, "m05"), 3L)
  expect_identical(unique(adj$covariates_adjusted), "Env")
})

test_that("a numeric covariate in fixed is adjusted for as well", {
  d <- .met_confounded()
  d$ph$Shade <- as.numeric(d$ph$Env == "E2") * 10 + rnorm(nrow(d$ph), sd = 0.01)   # tracks Env
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ Shade)
  expect_gt(.rank_of(adj, "m05"), 3L)
  expect_identical(unique(adj$covariates_adjusted), "Shade")
})

test_that("ridge scoring with covariates equals ridge with unpenalized covariate effects", {
  d <- .met_confounded()
  ph <- d$ph
  x <- d$X[ph$GID, ]
  meta <- feature_score_predictors(d$X, ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                   heter_groups = "Env")
  C <- cbind(1, as.numeric(ph$Env == "E2"))
  qc <- qr(C)
  xm <- x - qr.fitted(qc, x)
  Z0 <- sweep(x, 2, apply(xm, 2, stats::sd), "/")   # standardized, not residualized: C unpenalized
  ys <- ph$Yield / stats::sd(ph$Yield - qr.fitted(qc, ph$Yield))
  A <- cbind(C, Z0)
  beta <- solve(crossprod(A) + diag(c(0, 0, rep(1, ncol(Z0)))), crossprod(A, ys))[-(1:2)]
  expect_equal(meta$raw_score[match(colnames(x), meta$predictor)], abs(as.numeric(beta)),
               tolerance = 1e-8)
})

test_that("environment adjustment works for BayesB scoring", {
  skip_if_not_installed("BGLR")
  d <- .met_confounded()
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "BayesB",
                                  heter_groups = "Env", bayes_nIter = 1200L, bayes_burnIn = 400L,
                                  bayes_thin = 4L, seed = 5L)
  expect_equal(.rank_of(adj, "m01"), 1L)
  expect_gt(.rank_of(adj, "m05"), 3L)
})

test_that("fixed = ~1 and no environment give the scores of the unadjusted call", {
  d <- .met_confounded()
  a <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression")
  b <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                fixed = ~1, heter_groups = NULL)
  expect_identical(a, b)
  expect_identical(unique(a$covariates_adjusted), "")
})

test_that("records with a missing covariate are dropped and counted", {
  d <- .met_confounded()
  d$ph$Env[c(1, 2)] <- NA
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  heter_groups = "Env")
  expect_equal(unique(adj$n_records_dropped), 2L)
  expect_equal(unique(adj$n_train), nrow(d$ph) - 2L)
})

test_that("a fixed term with one level among the scored records is dropped, not an error", {
  d <- .met_confounded()
  one_env <- d$ph[d$ph$Env == "E1", ]
  adj <- feature_score_predictors(d$X, one_env, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ Env)
  expect_identical(unique(adj$covariates_adjusted), "")
})

test_that("a fixed term naming a missing column stops with the column's name", {
  d <- .met_confounded()
  expect_error(
    feature_score_predictors(d$X, d$ph, "Yield", "GID", fixed = ~ Block),
    "column\\(s\\) not in the phenotype data: Block"
  )
})

test_that("a marker monomorphic among the scored records scores exactly 0 under adjustment", {
  d <- .met_confounded()
  d$X[, "m07"] <- 1
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  heter_groups = "Env")
  expect_identical(adj$raw_score[adj$predictor == "m07"], 0)
})

test_that("a marker constant within each environment scores exactly 0 under adjustment", {
  set.seed(32)
  lines <- sprintf("L%02d", 1:40)
  X <- matrix(rbinom(40 * 20, 2, 0.4), 40, dimnames = list(lines, sprintf("m%02d", 1:20)))
  X[1:20, "m08"] <- 0
  X[21:40, "m08"] <- 2
  ph <- rbind(data.frame(GID = lines[1:20], Env = "E1"), data.frame(GID = lines[21:40], Env = "E2"))
  ph$Yield <- ifelse(ph$Env == "E2", 5, 0) + 1.5 * X[ph$GID, "m01"] + rnorm(nrow(ph), sd = 0.5)
  adj <- feature_score_predictors(X, ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  heter_groups = "Env")
  expect_identical(adj$raw_score[adj$predictor == "m08"], 0)
})

test_that("a single-level factor made inside the formula is dropped, not an error", {
  d <- .met_confounded()
  one_env <- d$ph[d$ph$Env == "E1", ]
  one_env$EnvCode <- 1
  adj <- feature_score_predictors(d$X, one_env, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ factor(EnvCode))
  base <- feature_score_predictors(d$X, one_env, "Yield", "GID", scoring_model = "Ridge_Regression")
  expect_identical(unique(adj$covariates_adjusted), "")
  expect_identical(adj, base)
})

test_that("a covariate with an infinite value is dropped and counted", {
  d <- .met_confounded()
  d$ph$Shade <- as.numeric(d$ph$Env == "E2") * 10 + rnorm(nrow(d$ph), sd = 0.01)
  d$ph$Shade[c(1, 2)] <- c(Inf, -Inf)
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ Shade)
  expect_equal(unique(adj$n_records_dropped), 2L)
  expect_equal(unique(adj$n_train), nrow(d$ph) - 2L)
})

test_that("RandomForest scoring asks for AIR and keeps its sign", {
  captured <- NULL
  local_mocked_bindings(
    gp_py_ml_feature_importance = function(model_type, X_train, y_train, response_family = "gaussian",
                                           model_params = list(), class_levels = NULL) {
      captured <<- list(params = model_params, ncol = ncol(X_train))
      c(-0.3, 0.1, 0.05)
    },
    .package = "PredictProR"
  )
  x <- data.frame(m1 = 1:6, m2 = 6:1, m3 = c(1, 1, 0, 0, 1, 0))
  rownames(x) <- paste0("g", 1:6)
  ph <- data.frame(ID = rownames(x), Y = c(1, 2, 3, 4, 5, 6))
  meta <- feature_score_predictors(x, ph, "Y", "ID", scoring_model = "RandomForest", ntree = 20L)
  expect_identical(captured$params$importance, "air")
  expect_equal(captured$params$air_unranked, 0L)
  expect_equal(meta$predictor[order(meta$rank)], c("m2", "m3", "m1"))   # -0.3 ranks last
  expect_identical(unique(meta$scoring_engine), "sklearn RandomForest AIR")
})

test_that("non-Gaussian RandomForest scoring passes covariates as unranked columns", {
  captured <- NULL
  local_mocked_bindings(
    gp_py_ml_feature_importance = function(model_type, X_train, y_train, response_family = "gaussian",
                                           model_params = list(), class_levels = NULL) {
      captured <<- list(params = model_params, ncol = ncol(X_train))
      c(0.3, 0.1, 0.2)
    },
    .package = "PredictProR"
  )
  x <- data.frame(m1 = c(0, 1, 2, 0, 1, 2, 0, 1), m2 = c(2, 2, 1, 1, 0, 0, 1, 2), m3 = c(1, 0, 1, 0, 1, 0, 1, 0))
  rownames(x) <- paste0("g", 1:8)
  ph <- data.frame(ID = rep(rownames(x), 2), Env = rep(c("A", "B"), each = 8),
                   Y = factor(rep(c("lo", "hi"), 8)))
  meta <- feature_score_predictors(x, ph, "Y", "ID", scoring_model = "RandomForest",
                                   response_family = "binary", heter_groups = "Env", ntree = 20L)
  expect_equal(captured$ncol, 4L)                 # 3 markers + 1 environment indicator
  expect_equal(captured$params$air_unranked, 1L)
  expect_equal(nrow(meta), 3L)
})

test_that("AIR removes impurity importance's bias towards a continuous noise predictor", {
  skip_on_cran()
  skip_if_no_ml_python()
  noise_rank <- function(seed, importance) {
    set.seed(seed)
    snp <- matrix(rbinom(200 * 30, 2, 0.3), 200, 30)
    X <- cbind(snp, rnorm(200))
    colnames(X) <- c(sprintf("s%02d", 1:30), "noise")
    y <- snp[, 1] + snp[, 2] + snp[, 3] + rnorm(200)
    imp <- gp_py_ml_feature_importance(
      "randomforest", X, y, "gaussian",
      list(ntree = 300L, n_jobs = 1L, random_state = seed, importance = importance)
    )
    match("noise", colnames(X)[order(-imp)])
  }
  mdi <- vapply(1:10, noise_rank, integer(1), importance = "impurity")
  air <- vapply(1:10, noise_rank, integer(1), importance = "air")
  # spec B3 evidence (120 data sets): MDI puts the noise predictor in the top 4 in every data set,
  # AIR in about 22% (a null SNP: about 3%). AIR's ranks are spread out (about 38% in the bottom 4),
  # so a median over few data sets says little; the top-4 share is the measure.
  expect_gte(sum(mdi <= 4L), 8L)
  expect_lte(sum(air <= 4L), 5L)
})

test_that("AIR does not penalise markers in LD blocks under a global null", {
  skip_on_cran()
  skip_if_no_ml_python()
  gap <- function(seed) {
    set.seed(seed)
    n <- 200
    block <- do.call(cbind, lapply(1:5, function(b) {
      base <- rbinom(n, 2, 0.3)
      vapply(1:6, function(k) {
        v <- base
        redraw <- runif(n) < 0.05
        v[redraw] <- rbinom(sum(redraw), 2, 0.3)
        v
      }, numeric(n))
    }))
    unlinked <- matrix(rbinom(n * 30, 2, 0.3), n, 30)
    X <- cbind(block, unlinked)
    colnames(X) <- c(sprintf("b%02d", 1:30), sprintf("u%02d", 1:30))
    y <- rnorm(n)
    imp <- gp_py_ml_feature_importance(
      "randomforest", X, y, "gaussian",
      list(ntree = 200L, n_jobs = 1L, random_state = seed, importance = "air")
    )
    mean(imp[1:30]) - mean(imp[31:60])
  }
  gaps <- vapply(1:10, gap, numeric(1))
  expect_gt(mean(gaps), -0.5)
})

test_that("the CV engine's fold-internal scoring adjusts for the run's environments", {
  d <- .met_confounded()
  ph <- d$ph
  local_mocked_bindings(
    predict_with_model = function(model, y, omics_data = NULL, tst, additional_params = list(), ...) {
      rep(mean(as.numeric(y[-tst]), na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )
  res <- models_execute_crossval(
    pheno_data = ph,
    params = list(
      response = "Yield", gen_name = "GID", random_state = 3L, replication = 1L,
      cross_validation_meth = "K-Folds", nfolds = 3L, eval_metrics = "mean_squared_error",
      response_family = "gaussian", GS_model_cv = "RandomForest",
      ml_dat_res = list(merged_data = list(merge_data = d$X[ph$GID, ]), omic_count = 1L),
      feature_score_metadata = feature_score_predictors(d$X, ph, "Yield", "GID"),
      feature_k_grid = 5L, feature_scoring_model = "Ridge_Regression",
      feature_scoring_cv = "fold_internal", heter_groups = "Env", fixed = ~1,
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )
  fold_meta <- res[[1]]$feature_selection_metadata
  expect_true("covariates_adjusted" %in% names(fold_meta))
  expect_gt(nrow(fold_meta), 0L)
  expect_true(all(fold_meta$covariates_adjusted == "Env"))
})

test_that("model_execute passes fixed and the environment column to scoring", {
  d <- .met_confounded()
  captured <- NULL
  local_mocked_bindings(
    feature_score_predictors = function(...) {
      captured <<- list(...)
      stop("captured scoring call")
    },
    gp_python_prepare_runtimes = function(...) invisible(NULL),
    .package = "PredictProR"
  )
  withr::local_dir(withr::local_tempdir())
  expect_error(suppressWarnings(model_execute(
    pheno_data = d$ph, geno_data = d$X, gen_name = "GID", response = "Yield", heter_groups = "Env",
    fixed = ~ Env, random = ~ GID, GS_model = "GBLUP_BRR", feature_scoring = TRUE, feature_k = 5L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, gmatrix_method = "VanRaden",
    parallel_mode = "sequential", system_database = FALSE, message = FALSE
  )), "captured scoring call")
  expect_identical(captured$heter_groups, "Env")
  expect_true(inherits(captured$fixed, "formula"))
})

test_that("the selection kind is response, explicit or none", {
  meta <- data.frame(trait = "Y", predictor = paste0("m", 1:10), rank = 1:10, raw_score = 10:1)
  expect_identical(gp_feature_selection_kind(meta, 5L, NULL, "Y"), "response")
  expect_identical(gp_feature_selection_kind(meta, "all", NULL, "Y"), "none")
  expect_identical(gp_feature_selection_kind(meta, 10L, NULL, "Y"), "none")       # k >= predictors
  expect_identical(gp_feature_selection_kind(NULL, NULL, c("m1", "m2"), "Y"), "explicit")
  expect_identical(gp_feature_selection_kind(NULL, NULL, NULL, "Y"), "none")
})

test_that("the remark is added to variance components and predictions only", {
  res <- list(
    predicted_values = data.frame(GID = "a", Predicted_value = 1),
    variance_components = data.frame(Component = "Genetic", Estimate = 1),
    model_parameters = data.frame(stat = "k", summary = "5"),
    model_results_by_trait = list(Y = list(Variance_components = data.frame(x = 1))),
    cv_results_raw = list(variance_components = data.frame(x = 1))
  )
  out <- gp_feature_add_selection_remark(res, "response")
  expect_true("Feature_selection_remark" %in% names(out$predicted_values))
  expect_true("Feature_selection_remark" %in% names(out$variance_components))
  expect_true("Feature_selection_remark" %in% names(out$model_results_by_trait$Y$Variance_components))
  expect_false("Feature_selection_remark" %in% names(out$model_parameters))
  expect_false("Feature_selection_remark" %in% names(out$cv_results_raw$variance_components))
  expect_match(out$variance_components$Feature_selection_remark, "overestimated")
  expect_identical(gp_feature_add_selection_remark(res, "none"), res)
})

test_that("a true-prediction run with response-based selection carries the remark", {
  skip_if_not_installed("glmnet")
  set.seed(5)
  G <- matrix(rbinom(80 * 40, 2, 0.4), 80, dimnames = list(sprintf("g%02d", 1:80), sprintf("m%02d", 1:40)))
  ph <- data.frame(GID = rownames(G), Yield = drop(G[, 1:8] %*% c(2, -2, 1, 1, -1, 1, -1, 0.5)) + rnorm(80, sd = 0.5))
  ph$Yield[1:10] <- NA
  withr::local_dir(withr::local_tempdir())
  run <- function(...) suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "Lasso", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    gmatrix_method = "VanRaden", parallel_mode = "sequential", system_database = FALSE,
    message = FALSE, verbose = FALSE, ...
  )))
  sel <- run(feature_scoring = TRUE, feature_k = 10L)
  expect_true("Feature_selection_remark" %in% names(sel$model_results$predicted_values))
  vc <- Filter(is.data.frame, list(sel$model_results$variance_components, sel$model_results$Variance_components))
  expect_true(length(vc) > 0L)
  expect_true(all(vapply(vc, function(d) "Feature_selection_remark" %in% names(d), logical(1))))
  none <- run()
  expect_false("Feature_selection_remark" %in% names(none$model_results$predicted_values))
})

test_that("the final stage marks only the traits whose best model used a selection", {
  meta <- rbind(
    data.frame(trait = "A", predictor = paste0("m", 1:10), rank = 1:10, raw_score = 10:1),
    data.frame(trait = "B", predictor = paste0("m", 1:10), rank = 1:10, raw_score = 10:1)
  )
  mk <- function() list(
    predicted_values = data.frame(GID = "a", Predicted_value = 1),
    variance_components = data.frame(Component = "Genetic", Estimate = 1)
  )
  fr <- list(
    model_results_by_trait = list(A = mk(), B = mk()),
    model_results_by_model = list(Lasso = list(A = mk(), B = mk())),
    model_results = mk()
  )
  bm <- data.frame(trait = c("A", "B"), model = "Lasso", feature_k = c(5L, NA))
  out <- gp_feature_apply_final_selection_remarks(fr, bm, meta, NULL, NULL)
  has <- function(x) "Feature_selection_remark" %in% names(x)
  expect_true(has(out$model_results_by_trait$A$predicted_values))
  expect_true(has(out$model_results_by_trait$A$variance_components))
  expect_true(has(out$model_results_by_model$Lasso$A$predicted_values))
  expect_false(has(out$model_results_by_trait$B$predicted_values))
  expect_false(has(out$model_results_by_model$Lasso$B$variance_components))
  expect_false(has(out$model_results$predicted_values))   # two traits: left alone
})

test_that("auto resolves to fold_internal for package scores and fixed for a supplied ranking", {
  expect_identical(gp_feature_resolve_cv_policy("auto", ranking_supplied = FALSE), "fold_internal")
  expect_identical(gp_feature_resolve_cv_policy("auto", ranking_supplied = TRUE), "fixed")
  expect_identical(gp_feature_resolve_cv_policy("fixed", ranking_supplied = FALSE), "fixed")
  expect_identical(gp_feature_resolve_cv_policy("both", ranking_supplied = FALSE), "both")
  expect_error(gp_feature_normalize_cv_policy("auto"), "must be resolved")
})

test_that("a model_execute CV with package scoring defaults to fold-internal selection", {
  skip_if_not_installed("glmnet")
  set.seed(5)
  G <- matrix(rbinom(80 * 40, 2, 0.4), 80, dimnames = list(sprintf("g%02d", 1:80), sprintf("m%02d", 1:40)))
  ph <- data.frame(GID = rownames(G), Yield = drop(G[, 1:8] %*% c(2, -2, 1, 1, -1, 1, -1, 0.5)) + rnorm(80, sd = 0.5))
  withr::local_dir(withr::local_tempdir())
  cv <- function(...) suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model_cv = "Lasso", cross_validation = TRUE, cross_validation_meth = "K-Folds", nfolds = 3,
    replication = 1, cv_generate_plots = FALSE, cv_evaluation_only = TRUE,
    eval_metrics = "accuracy", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    gmatrix_method = "VanRaden", parallel_mode = "sequential", system_database = FALSE,
    message = FALSE, verbose = FALSE, feature_k_grid = 10L, ...
  )))
  pkg <- cv(feature_scoring = TRUE)
  expect_identical(unique(pkg$cv_results_raw[[1]]$eval_metrics_reps$feature_scoring_cv), "fold_internal")
  user <- cv(feature_score_metadata = feature_score_predictors(G, ph, "Yield", "GID"))
  expect_identical(unique(user$cv_results_raw[[1]]$eval_metrics_reps$feature_scoring_cv), "fixed")
  expect_true(all(user$cv_results_raw[[1]]$feature_selection_metadata$ranking_source == "user"))
})

# ---- final-review fix wave (0.20.218) ----

test_that("source-filtered selection keeps rank order whatever the metadata order (M-a)", {
  meta <- data.frame(trait = "Y", predictor = c("m01", "m04", "m06", "m07", "m02"),
                     source_block = "geno_data", rank = c(1L, 4L, 3L, 2L, 5L),
                     raw_score = c(5, 2, 3, 4, 1), stringsAsFactors = FALSE)
  x <- matrix(0, 2, 5, dimnames = list(NULL, c("m01", "m02", "m04", "m06", "m07")))
  unsorted <- gp_feature_select_matrix_by_metadata(x, meta, "Y", 4L, source = "geno_data")
  sorted <- gp_feature_select_matrix_by_metadata(x, meta[order(meta$rank), ], "Y", 4L, source = "geno_data")
  expect_identical(colnames(unsorted), c("m01", "m07", "m06", "m04"))
  expect_identical(colnames(unsorted), colnames(sorted))
})

test_that("a QTL or a trait on a tiny scale keeps its ranking under adjustment (M-b)", {
  d <- .met_confounded()
  ref <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  heter_groups = "Env")
  X <- d$X
  X[, "m01"] <- X[, "m01"] * 1e-9
  tiny_x <- feature_score_predictors(X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                     heter_groups = "Env")
  expect_equal(.rank_of(tiny_x, "m01"), 1L)
  expect_equal(tiny_x$raw_score[match(ref$predictor, tiny_x$predictor)], ref$raw_score,
               tolerance = 1e-6)
  ph <- d$ph
  ph$Yield <- ph$Yield * 1e-9
  tiny_y <- feature_score_predictors(d$X, ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                     heter_groups = "Env")
  expect_identical(tiny_y$predictor, ref$predictor)
  expect_equal(tiny_y$raw_score, ref$raw_score, tolerance = 1e-6)
})

test_that("a fixed term base R cannot evaluate is dropped with a warning; the rest is adjusted (M-c)", {
  gp_feature_reset_run_state()
  d <- .met_confounded()
  d$ph$Rep <- factor(rep(1:2, length.out = nrow(d$ph)))
  expect_warning(
    adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                    fixed = ~ Env + at(Env):Rep),
    "at(Env):Rep", fixed = TRUE
  )
  env_only <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                       fixed = ~ Env)
  expect_identical(unique(adj$covariates_adjusted), "Env")
  expect_identical(adj, env_only)
  expect_gt(.rank_of(adj, "m05"), 3L)
})

test_that("a matrix-valued fixed term such as poly() is adjusted for (M-c)", {
  d <- .met_confounded()
  set.seed(33)
  d$ph$Shade <- as.numeric(d$ph$Env == "E2") * 10 + rnorm(nrow(d$ph), sd = 0.01)
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ poly(Shade, 2))
  expect_identical(unique(adj$covariates_adjusted), "poly(Shade, 2)")
  expect_equal(.rank_of(adj, "m01"), 1L)
  expect_gt(.rank_of(adj, "m05"), 2L)         # 1st without adjustment
  d$ph$Shade[3] <- NA
  adj_na <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                     fixed = ~ Shade)
  expect_equal(unique(adj_na$n_records_dropped), 1L)
})

test_that("covariates_adjusted lists only the terms C keeps (M-i)", {
  d <- .met_confounded()
  d$ph$Shade <- 5
  adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                  fixed = ~ Shade)
  base <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression")
  expect_identical(unique(adj$covariates_adjusted), "")
  expect_identical(adj$raw_score, base$raw_score)
  d$ph$Shade <- ifelse(d$ph$Env == "E2", 10, 0)      # aliased with Env
  both <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                   fixed = ~ Env + Shade)
  expect_identical(unique(both$covariates_adjusted), "Env")
})

test_that("BayesB scores markers standardized within the covariates, with C as a FIXED term (M-f)", {
  skip_if_not_installed("BGLR")
  d <- .met_confounded()
  captured <- NULL
  local_mocked_bindings(
    BGLR = function(y, ETA, ...) {
      captured <<- ETA
      list(ETA = list(markers = list(b = rep(0, ncol(ETA$markers$X)))))
    },
    .package = "BGLR"
  )
  feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "BayesB", heter_groups = "Env",
                           bayes_nIter = 10L, bayes_burnIn = 5L, bayes_thin = 1L)
  C <- gp_feature_covariate_design(d$ph, heter_groups = "Env")$C
  expect_equal(unname(captured$markers$X), unname(gp_feature_scale_matrix(d$X[d$ph$GID, ], covariates = C)))
  expect_identical(captured$covariates$model, "FIXED")
  expect_equal(ncol(captured$covariates$X), 1L)
  feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "BayesB",
                           bayes_nIter = 10L, bayes_burnIn = 5L, bayes_thin = 1L)
  expect_equal(unname(captured$markers$X), unname(gp_feature_scale_matrix(d$X[d$ph$GID, ])))
  expect_null(captured$covariates)
})

test_that("Gaussian RandomForest scoring with covariates sends the residualized trait (M-T4a)", {
  d <- .met_confounded()
  captured <- NULL
  local_mocked_bindings(
    gp_py_ml_feature_importance = function(model_type, X_train, y_train, response_family = "gaussian",
                                           model_params = list(), class_levels = NULL) {
      captured <<- list(y = y_train, ncol = ncol(X_train), params = model_params)
      seq_len(ncol(X_train)) / 10
    },
    .package = "PredictProR"
  )
  feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "RandomForest",
                           heter_groups = "Env", ntree = 20L)
  C <- gp_feature_covariate_design(d$ph, heter_groups = "Env")$C
  expect_equal(captured$y, gp_feature_residualize(as.double(d$ph$Yield), C))
  expect_equal(captured$ncol, ncol(d$X))           # no unranked covariate columns
  expect_equal(captured$params$air_unranked, 0L)
})

test_that("given training rows, fold scoring uses exactly those observed rows (I2)", {
  d <- .leak_data()
  test <- 1:10
  train <- 21:60          # rows 11-20 are neither, as the half-sibs Both_New_Parents keeps out
  view <- function(p) gp_feature_fold_view(
    predictor_data = d$X, pheno_data = p, response = "Yield", gen_name = "GID",
    test_rows = test, train_rows = train, k = 5, scoring_model = "Ridge_Regression", seed = 1L
  )
  v1 <- view(d$ph)
  ph2 <- d$ph
  ph2$Yield[11:20] <- d$X[ph2$GID[11:20], 7] * 50          # change only the unassigned rows
  v2 <- view(ph2)
  expect_equal(unique(v1$metadata$n_train), length(train))
  expect_identical(v1$selected_predictors, v2$selected_predictors)
  expect_identical(v1$metadata$raw_score, v2$metadata$raw_score)
  ph2$Yield[25] <- NA
  expect_equal(unique(view(ph2)$metadata$n_train), length(train) - 1L)
  mt <- gp_feature_trait_view(
    predictor_data = d$X, pheno_data = d$ph, response = "Yield", gen_name = "GID", k = 5,
    selection_mode = "fold_internal", test_rows = test, train_rows = train,
    scoring_model = "Ridge_Regression", seed = 1L
  )
  expect_equal(unique(mt$metadata$n_train), length(train))
  expect_identical(mt$selected_predictors, v1$selected_predictors)
})

.hybrid_toy <- function() {
  set.seed(41)
  grid <- expand.grid(Female = paste0("F", 1:4), Male = paste0("M", 1:4), stringsAsFactors = FALSE)
  grid$Hybrid <- paste(grid$Female, grid$Male, sep = "x")
  X <- matrix(rbinom(16 * 12, 2, 0.4), 16, dimnames = list(grid$Hybrid, sprintf("m%02d", 1:12)))
  grid$Yield <- X[, 1] + rnorm(16, sd = 0.3)
  meta <- feature_score_predictors(X, grid, "Yield", "Hybrid")
  list(ph = grid, X = X, meta = meta)
}

test_that("hybrid ML and DL CV score their scenario's training rows (I2)", {
  h <- .hybrid_toy()
  sc <- gp_hybrid_build_cv_scenarios(h$ph, "Yield", "Hybrid", "Female", "Male",
                                     cross_validation_meth = "Hybrid_Both_New_Parents")[[1L]]
  aligned <- gp_hybrid_ml_align_inputs(h$ph, "Yield", "Hybrid", "Female", "Male", h$X)
  rows <- gp_hybrid_cv_scenario_rows(sc, aligned$ph, "Hybrid")
  expect_lt(length(rows$train), nrow(h$ph) - length(rows$test))   # some rows are neither
  captured <- NULL
  local_mocked_bindings(
    gp_feature_trait_view = function(...) {
      captured <<- list(...)
      stop("captured view")
    },
    .package = "PredictProR"
  )
  common <- list(pheno_object = h$ph, response = "Yield", gen_name = "Hybrid",
                 female_parent = "Female", male_parent = "Male", geno_omic_object = h$X,
                 cross_validation_meth = "Hybrid_Both_New_Parents", feature_score_metadata = h$meta,
                 feature_k = 3L, feature_scoring_cv = "fold_internal")
  expect_error(do.call(gp_hybrid_ml_gaussian_cv, c(list(model_type = "Ridge_Regression"), common)),
               "captured view")
  expect_identical(as.integer(captured$train_rows), as.integer(rows$train))
  captured <- NULL
  expect_error(do.call(gp_hybrid_dl_gaussian_cv, c(list(model_type = "mlp"), common)), "captured view")
  expect_identical(as.integer(captured$train_rows), as.integer(rows$train))
})

.single_env <- function() {
  set.seed(7)
  X <- matrix(rbinom(60 * 20, 2, 0.4), 60, dimnames = list(sprintf("g%02d", 1:60), sprintf("m%02d", 1:20)))
  ph <- data.frame(GID = rownames(X), Yield = 2 * X[, 1] - X[, 2] + rnorm(60), stringsAsFactors = FALSE)
  list(X = X, ph = ph)
}
.mean_predictor <- function(model, y, omics_data = NULL, tst, additional_params = list(), ...) {
  rep(mean(as.numeric(y[-tst]), na.rm = TRUE), length(tst))
}
.cv_engine <- function(d, models, k_grid, metadata, mode) {
  models_execute_crossval(
    pheno_data = d$ph,
    params = list(
      response = "Yield", gen_name = "GID", random_state = 3L, replication = 1L,
      cross_validation_meth = "K-Folds", nfolds = 3L, eval_metrics = "mean_squared_error",
      response_family = "gaussian", GS_model_cv = models,
      ml_dat_res = list(merged_data = list(merge_data = d$X[d$ph$GID, ]), omic_count = 1L),
      feature_score_metadata = metadata, feature_k_grid = k_grid,
      feature_scoring_model = "Ridge_Regression", feature_scoring_cv = mode, fixed = ~1,
      parallel_mode = "sequential"
    ),
    verbose = FALSE
  )
}

test_that("fold-internal CV scores each fold once for all models and k values (I1)", {
  d <- .single_env()
  meta <- feature_score_predictors(d$X, d$ph, "Yield", "GID")
  orig <- gp_feature_score_training_records
  n_scored <- 0L
  local_mocked_bindings(
    gp_feature_score_training_records = function(...) {
      n_scored <<- n_scored + 1L
      orig(...)
    },
    predict_with_model = .mean_predictor,
    .package = "PredictProR"
  )
  cached <- .cv_engine(d, c("RandomForest", "Xgboost"), c(5L, 10L), meta, "fold_internal")
  expect_equal(n_scored, 3L)                          # 3 folds, shared by 2 models x 2 k values
  expect_length(ls(gp_feature_fold_cache), 0L)        # emptied when the CV ends
  local_mocked_bindings(gp_feature_fold_cache_get = function(key) NULL, .package = "PredictProR")
  n_scored <- 0L
  uncached <- .cv_engine(d, c("RandomForest", "Xgboost"), c(5L, 10L), meta, "fold_internal")
  expect_equal(n_scored, 12L)
  expect_identical(lapply(cached, `[[`, "feature_selection_metadata"),
                   lapply(uncached, `[[`, "feature_selection_metadata"))
  expect_identical(lapply(cached, `[[`, "eval_metrics_reps"), lapply(uncached, `[[`, "eval_metrics_reps"))
  expect_identical(lapply(cached, `[[`, "ypred_cv_Reps_all"), lapply(uncached, `[[`, "ypred_cv_Reps_all"))
})

test_that("a fold ranking is reused only for identical scoring inputs (I1)", {
  d <- .single_env()
  view <- function(p, seed = 1L, hash = "h1") gp_feature_fold_view(
    predictor_data = d$X, pheno_data = p, response = "Yield", gen_name = "GID", test_rows = 1:20,
    k = 5, scoring_model = "Ridge_Regression", seed = seed, predictor_hash = hash
  )
  gp_feature_fold_cache_clear()
  on.exit(gp_feature_fold_cache_clear(), add = TRUE)
  a <- view(d$ph)
  expect_length(ls(gp_feature_fold_cache), 1L)
  view(d$ph)
  expect_length(ls(gp_feature_fold_cache), 1L)        # same inputs: reused
  ph2 <- d$ph
  ph2$Yield[30] <- ph2$Yield[30] + 1                  # a training record changed
  view(ph2)
  view(d$ph, seed = 2L)
  view(d$ph, hash = "h2")
  expect_length(ls(gp_feature_fold_cache), 4L)
  gp_feature_fold_cache_clear()
  b <- gp_feature_fold_view(
    predictor_data = d$X, pheno_data = d$ph, response = "Yield", gen_name = "GID", test_rows = 1:20,
    k = 5, scoring_model = "Ridge_Regression", seed = 1L
  )
  expect_length(ls(gp_feature_fold_cache), 0L)        # no predictor_hash: not memoised
  expect_identical(a, b)
})

test_that("a tie found in the CV's tasks warns once in the main process (M-d)", {
  d <- .single_env()
  tied <- feature_score_predictors(d$X, d$ph, "Yield", "GID")
  tied$raw_score <- 1
  orig <- gp_execute_partitioned_tasks
  local_mocked_bindings(
    predict_with_model = .mean_predictor,
    # as on a parallel worker: the task's warning and flag stay in the worker
    gp_execute_partitioned_tasks = function(...) {
      res <- suppressWarnings(orig(...))
      gp_feature_reset_tie_warning()
      res
    },
    .package = "PredictProR"
  )
  gp_feature_reset_tie_warning()
  expect_warning(.cv_engine(d, "RandomForest", 5L, tied, "fixed"), "share the score at the top-5 cut-off")
  recs <- list(list(feature_selection_metadata = data.frame(
    trait = "Y", predictor = c("a", "b", "c"), rank = 1:3, raw_score = c(2, 1, 1),
    replication = 1L, fold = 1L, selection_mode = "fold_internal", feature_k = 2L,
    tied_at_cutoff = TRUE
  )))
  expect_no_warning(gp_feature_warn_returned_ties(recs))   # already warned in this run
  gp_feature_reset_tie_warning()
  expect_warning(gp_feature_warn_returned_ties(recs), "2 predictors share the score at the top-2")
  expect_no_warning(gp_feature_warn_returned_ties(recs))
})

test_that("a run leaves no tie flag, warned terms or memoised rankings behind (M-d)", {
  skip_if_not_installed("glmnet")
  set.seed(5)
  G <- matrix(rbinom(80 * 40, 2, 0.4), 80, dimnames = list(sprintf("g%02d", 1:80), sprintf("m%02d", 1:40)))
  ph <- data.frame(GID = rownames(G), Yield = drop(G[, 1:8] %*% c(2, -2, 1, 1, -1, 1, -1, 0.5)) + rnorm(80, sd = 0.5))
  local_mocked_bindings(
    feature_score_predictors = function(predictor_data, pheno_data, response, gen_name, ...) {
      gp_feature_metadata_table(
        trait = response[[1L]], predictors = colnames(predictor_data),
        source_block = gp_feature_source_block(colnames(predictor_data)),
        raw_score = rep(1, ncol(predictor_data)), scoring_model = "Ridge_Regression", seed = 1L,
        n_train = nrow(pheno_data)
      )
    },
    .package = "PredictProR"
  )
  withr::local_dir(withr::local_tempdir())
  gp_feature_reset_tie_warning()
  suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "Lasso", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    gmatrix_method = "VanRaden", parallel_mode = "sequential", system_database = FALSE,
    message = FALSE, verbose = FALSE, feature_scoring = TRUE, feature_k = 10L
  )))
  expect_false(isTRUE(gp_feature_state$tie_warned))
  expect_length(gp_feature_state$fixed_terms_warned, 0L)
  expect_length(ls(gp_feature_fold_cache), 0L)
})

test_that("a supplied ranking gets its own remark kind and wording (R-2)", {
  meta <- data.frame(trait = "Y", predictor = paste0("m", 1:10), rank = 1:10, raw_score = 10:1,
                     ranking_source = "user")
  expect_identical(gp_feature_selection_kind(meta, 5L, NULL, "Y"), "supplied_ranking")
  expect_identical(gp_feature_selection_kind(meta, "all", NULL, "Y"), "none")
  meta$ranking_source <- "package"
  expect_identical(gp_feature_selection_kind(meta, 5L, NULL, "Y"), "response")
  txt <- gp_feature_selection_remark_text("supplied_ranking")
  expect_match(txt, "selected with a supplied ranking (feature_score_metadata)", fixed = TRUE)
  expect_match(txt, "If that ranking was computed using these phenotypes", fixed = TRUE)
  out <- gp_feature_add_selection_remark(list(variance_components = data.frame(x = 1)), "supplied_ranking")
  expect_identical(out$variance_components$Feature_selection_remark, txt)
})

test_that("the remark marks genetic correlation and covariance tables (M-e)", {
  cor_mat <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(c("A", "B"), c("A", "B")))
  res <- list(
    Genetic_correlation_traits = cor_mat,
    Genetic_covariance_traits = data.frame(Trait = c("A", "B"), A = c(1, 0.5), B = c(0.5, 1)),
    Genetic_correlation_environments_SE = cor_mat,
    genetic_correlation = cor_mat,
    gxe_covariance = cor_mat,
    environment_correlation = cor_mat,
    covariance_components = data.frame(Component = "G", Estimate = 1),
    Covariance = list(`GID:Env` = cor_mat),
    model_parameters = data.frame(stat = "k", summary = "5"),
    cv_results_processed = list(environment_correlation = cor_mat)
  )
  out <- gp_feature_add_selection_remark(res, "response")
  for (nm in c("Genetic_correlation_traits", "Genetic_correlation_environments_SE",
               "genetic_correlation", "gxe_covariance", "environment_correlation")) {
    expect_match(attr(out[[nm]], "Feature_selection_remark"), "overestimated", info = nm)
    expect_equal(out[[nm]], cor_mat, ignore_attr = "Feature_selection_remark", info = nm)
  }
  expect_true("Feature_selection_remark" %in% names(out$Genetic_covariance_traits))
  expect_true("Feature_selection_remark" %in% names(out$covariance_components))
  expect_match(attr(out$Covariance$`GID:Env`, "Feature_selection_remark"), "overestimated")
  expect_false("Feature_selection_remark" %in% names(out$model_parameters))
  expect_null(attr(out$cv_results_processed$environment_correlation, "Feature_selection_remark"))
  csv <- gp_feature_remark_export_table(out$genetic_correlation)
  expect_identical(names(csv), c("A", "B", "Feature_selection_remark"))
  expect_identical(rownames(csv), c("A", "B"))
})

.lasso_data <- function() {
  set.seed(5)
  G <- matrix(rbinom(80 * 40, 2, 0.4), 80, dimnames = list(sprintf("g%02d", 1:80), sprintf("m%02d", 1:40)))
  ph <- data.frame(GID = rownames(G), Yield = drop(G[, 1:8] %*% c(2, -2, 1, 1, -1, 1, -1, 0.5)) + rnorm(80, sd = 0.5))
  ph$Yield[1:10] <- NA
  list(G = G, ph = ph)
}
.lasso_tp <- function(d, ...) {
  suppressWarnings(suppressMessages(model_execute(
    pheno_data = d$ph, geno_data = d$G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model = "Lasso", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    gmatrix_method = "VanRaden", parallel_mode = "sequential", message = FALSE, verbose = FALSE, ...
  )))
}
.exported <- function(name) {
  list.files(".", pattern = paste0("^", name, "[.]csv$"), recursive = TRUE, full.names = TRUE)
}

test_that("the exported prediction and variance-component files carry the remark (I3)", {
  skip_if_not_installed("glmnet")
  d <- .lasso_data()
  withr::local_dir(withr::local_tempdir())
  sel <- .lasso_tp(d, feature_scoring = TRUE, feature_k = 10L)        # default output: files
  pv <- .exported("Predicted_Value")
  vc <- .exported("variance_components")
  expect_length(pv, 1L)
  expect_length(vc, 1L)
  pv_tab <- utils::read.csv(pv)
  vc_tab <- utils::read.csv(vc)
  expect_true("Feature_selection_remark" %in% names(pv_tab))
  expect_true("Feature_selection_remark" %in% names(vc_tab))
  expect_identical(unique(pv_tab$Feature_selection_remark),
                   unique(sel$model_results$predicted_values$Feature_selection_remark))
  returned_vc <- sel$model_results$variance_components %||% sel$model_results$Variance_components
  expect_identical(unique(vc_tab$Feature_selection_remark), unique(returned_vc$Feature_selection_remark))
  expect_match(unique(vc_tab$Feature_selection_remark), "with these phenotypes")
})

test_that("runs with a supplied ranking or set say so; runs without selection write no remark (R-2, M-h)", {
  skip_if_not_installed("glmnet")
  d <- .lasso_data()
  withr::local_dir(withr::local_tempdir())
  user <- .lasso_tp(d, feature_score_metadata = feature_score_predictors(d$G, d$ph, "Yield", "GID"),
                    feature_k = 10L)
  expect_match(unique(user$model_results$predicted_values$Feature_selection_remark),
               "selected with a supplied ranking", fixed = TRUE)
  expect_match(unique(utils::read.csv(.exported("Predicted_Value"))$Feature_selection_remark),
               "selected with a supplied ranking", fixed = TRUE)
  withr::local_dir(withr::local_tempdir())
  explicit <- .lasso_tp(d, feature_selected = sprintf("m%02d", 1:8))
  expect_match(unique(explicit$model_results$predicted_values$Feature_selection_remark),
               "user-supplied marker set (feature_selected)", fixed = TRUE)
  withr::local_dir(withr::local_tempdir())
  none <- .lasso_tp(d)
  expect_false("Feature_selection_remark" %in% names(none$model_results$predicted_values))
  expect_false("Feature_selection_remark" %in% names(utils::read.csv(.exported("Predicted_Value"))))
})

test_that("environment adjustment removes an environment-tracking marker from the top (RandomForest, M-h)", {
  skip_on_cran()
  skip_if_no_ml_python()
  d <- .met_confounded()
  ranks <- vapply(1:5, function(s) {
    raw <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "RandomForest",
                                    ntree = 300L, seed = s)
    adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "RandomForest",
                                    ntree = 300L, seed = s, heter_groups = "Env")
    c(raw_m05 = .rank_of(raw, "m05"), adj_m05 = .rank_of(adj, "m05"), adj_m01 = .rank_of(adj, "m01"))
  }, numeric(3))
  expect_true(all(ranks["raw_m05", ] <= 2))         # unadjusted: the environment marker leads
  expect_gte(sum(ranks["adj_m05", ] > 3), 4L)       # adjusted: it no longer does
  expect_gte(sum(ranks["adj_m01", ] == 1), 4L)      # and the QTL does
})

test_that("feature_scoring_cv = 'both' gives fixed and fold-internal CV result sets (M-h)", {
  skip_if_not_installed("glmnet")
  d <- .lasso_data()
  d$ph <- d$ph[!is.na(d$ph$Yield), ]
  withr::local_dir(withr::local_tempdir())
  cv <- suppressWarnings(suppressMessages(model_execute(
    pheno_data = d$ph, geno_data = d$G, gen_name = "GID", response = "Yield", fixed = ~1, random = ~GID,
    GS_model_cv = "Lasso", cross_validation = TRUE, cross_validation_meth = "K-Folds", nfolds = 3,
    replication = 1, cv_generate_plots = FALSE, cv_evaluation_only = TRUE,
    eval_metrics = "accuracy", qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    gmatrix_method = "VanRaden", parallel_mode = "sequential", system_database = FALSE,
    message = FALSE, verbose = FALSE, feature_scoring = TRUE, feature_k_grid = 10L,
    feature_scoring_cv = "both"
  )))
  modes <- vapply(cv$cv_results_raw, function(e) unique(e$eval_metrics_reps$feature_scoring_cv), character(1))
  expect_setequal(unique(modes), c("fixed", "fold_internal"))
  fold_meta <- do.call(rbind, lapply(cv$cv_results_raw[modes == "fold_internal"], function(e) {
    m <- e$feature_selection_metadata
    if (is.data.frame(m)) m[, c("selection_mode", "fold")] else NULL
  }))
  expect_gt(NROW(fold_meta), 0L)
  expect_true(all(fold_meta$selection_mode == "fold_internal"))
  expect_true(all(is.finite(fold_meta$fold)))
})

test_that("a fixed term that fails only on the complete records is dropped, not an error (M-c)", {
  gp_feature_reset_run_state()
  d <- .met_confounded()
  d$ph$Shade <- rep(c(0, 5, 10), length.out = nrow(d$ph))
  d$ph$Block <- factor(rep(c("a", "b"), length.out = nrow(d$ph)))
  d$ph$Block[d$ph$Shade == 5] <- NA        # dropped records leave Shade two distinct values
  expect_warning(
    adj <- feature_score_predictors(d$X, d$ph, "Yield", "GID", scoring_model = "Ridge_Regression",
                                    fixed = ~ Block + poly(Shade, 2)),
    "poly(Shade, 2)", fixed = TRUE
  )
  expect_identical(unique(adj$covariates_adjusted), "Block")
  expect_equal(unique(adj$n_records_dropped), sum(d$ph$Shade == 5))
})
