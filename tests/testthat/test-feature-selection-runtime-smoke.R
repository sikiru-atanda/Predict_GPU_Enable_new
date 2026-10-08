if (!identical(Sys.getenv("PREDICTPROR_RUN_FEATURE_SELECTION_SMOKE"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_FEATURE_SELECTION_SMOKE=1 to run feature-selection runtime tests.")
}

test_that("RandomForest scoring supports binary, ordinal, and nominal multiclass traits", {
  set.seed(8100)
  ids <- paste0("g", seq_len(30))
  x <- matrix(rnorm(30 * 5), nrow = 30, dimnames = list(ids, paste0("m", 1:5)))
  pheno <- data.frame(
    ID = ids,
    Binary = factor(ifelse(x[, "m1"] > 0, "case", "control")),
    Ordinal = ordered(
      cut(x[, "m2"], breaks = c(-Inf, -0.4, 0.4, Inf), labels = c("low", "mid", "high")),
      levels = c("low", "mid", "high")
    ),
    Nominal = factor(c("red", "green", "blue")[(max.col(cbind(x[, "m3"], x[, "m4"], x[, "m5"])))],
      levels = c("red", "green", "blue")),
    stringsAsFactors = FALSE
  )

  scores <- feature_score_predictors(
    predictor_data = x,
    pheno_data = pheno,
    response = c("Binary", "Ordinal", "Nominal"),
    gen_name = "ID",
    scoring_model = "RandomForest",
    response_family = "auto",
    ntree = 25L,
    rf_n_jobs = 1L,
    seed = 8100L
  )

  family_by_trait <- tapply(scores$response_family, scores$trait, unique)
  expect_equal(as.vector(family_by_trait[c("Binary", "Ordinal", "Nominal")]),
               c("binary", "ordinal", "multiclass"))
  expect_true(all(is.finite(scores$raw_score)))
  expect_true(all(tapply(scores$rank, scores$trait, function(z) identical(sort(z), 1:5))))
})

# geno_data must hold dosages (0/1/2 for diploids): genotype QC rejects
# continuous values above 2 without ploidy metadata.
dose <- function(v) findInterval(v, stats::quantile(v, c(1 / 3, 2 / 3)))

test_that("multi-trait ML uses fold-internal trait-specific feature sets", {
  set.seed(8101)
  ids <- paste0("g", seq_len(16))
  # inbred lines: homozygous 0/2 dosages (the default het_threshold = 0.1
  # would remove heterozygous markers)
  x <- ifelse(matrix(rnorm(16 * 6), nrow = 16) > 0, 2, 0)
  dimnames(x) <- list(ids, paste0("m", 1:6))
  storage.mode(x) <- "double"
  pheno <- data.frame(
    ID = ids,
    T1 = 1.2 * x[, "m1"] - 0.2 * x[, "m2"] + rnorm(16, sd = 0.05),
    T2 = 1.1 * x[, "m5"] + 0.2 * x[, "m6"] + rnorm(16, sd = 0.05),
    stringsAsFactors = FALSE
  )

  out <- model_execute(
    pheno_data = pheno,
    geno_data = x,
    response = c("T1", "T2"),
    gen_name = "ID",
    GS_model_cv = "RandomForest",
    response_family = "gaussian",
    multi_trait_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "root_mean_squared_error",
    feature_scoring = TRUE,
    feature_scoring_model = "Ridge_Regression",
    feature_scoring_cv = "fold_internal",
    feature_k = 2L,
    ntree = 25L,
    system_database = TRUE,
    message = FALSE
  )

  meta <- out$cv_results_processed$feature_selection_metadata
  expect_true(is.data.frame(meta))
  expect_equal(sort(unique(meta$trait)), c("T1", "T2"))
  expect_true(all(meta$selection_mode == "fold_internal"))
  expect_equal(sort(unique(meta$fold)), 1:2)
  expect_true(all(tapply(meta$selected, interaction(meta$trait, meta$fold), sum) == 2L))
  expect_equal(length(unique(meta$training_id_hash)), 2L)
})

test_that("hybrid ML uses fold-internal hybrid-level feature sets", {
  set.seed(8102)
  pheno <- data.frame(
    HybridID = c("A_M1", "A_M2", "B_M1", "B_M2", "C_M1", "C_M2", "D_M1", "D_M2"),
    Female = rep(c("A", "B", "C", "D"), each = 2),
    Male = rep(c("M1", "M2"), 4),
    stringsAsFactors = FALSE
  )
  female_effect <- c(A = 0.8, B = 0.3, C = -0.4, D = 0.1)
  male_effect <- c(M1 = 0.5, M2 = -0.2)
  pheno$Yield <- female_effect[pheno$Female] + male_effect[pheno$Male] +
    rnorm(nrow(pheno), sd = 0.04)
  x <- apply(matrix(rnorm(nrow(pheno) * 6), nrow = nrow(pheno)), 2, dose)
  dimnames(x) <- list(pheno$HybridID, paste0("m", 1:6))
  storage.mode(x) <- "double"
  x[, "m1"] <- dose(female_effect[pheno$Female] + rnorm(nrow(pheno), sd = 0.05))
  x[, "m5"] <- dose(male_effect[pheno$Male] + rnorm(nrow(pheno), sd = 0.05))

  out <- model_execute(
    pheno_data = pheno,
    geno_data = x,
    response = "Yield",
    gen_name = "HybridID",
    female_parent = "Female",
    male_parent = "Male",
    GS_model_cv = "RandomForest",
    response_family = "gaussian",
    hybrid_ml = TRUE,
    cross_validation = TRUE,
    cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_Known_Parents",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "root_mean_squared_error",
    feature_scoring = TRUE,
    feature_scoring_model = "Ridge_Regression",
    feature_scoring_cv = "fold_internal",
    feature_k = 2L,
    ntree = 25L,
    system_database = TRUE,
    message = FALSE
  )

  meta <- out$cv_results_processed$feature_selection_metadata
  expect_true(is.data.frame(meta))
  expect_equal(unique(meta$trait), "Yield")
  expect_true(all(meta$selection_mode == "fold_internal"))
  expect_true(all(tapply(meta$selected, meta$fold, sum) == 2L))
  expect_true(all(out$cv_results_processed$hybrid_cv_predictions$Train_Test_Label == "Test"))
})
