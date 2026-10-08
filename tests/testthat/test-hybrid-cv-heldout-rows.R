test_that("hybrid CV predictions contain only the fold's held-out observed hybrids", {
  skip_if_not_installed("BGLR")
  # Each fold refits with its held-out hybrids masked to NA. Selecting every
  # row the refit labels "Test" also picked up the hybrids that were already
  # NA in pheno_data (true-prediction targets), so they were predicted in every
  # fold, inflated n_test_hybrids and bloated the CV prediction table.
  set.seed(42)
  females <- sprintf("F%02d", 1:8); males <- sprintf("M%02d", 1:6)
  parents <- c(females, males)
  geno <- matrix(sample(c(0, 2), length(parents) * 80, replace = TRUE), length(parents), 80,
                 dimnames = list(parents, sprintf("SNP%02d", 1:80)))
  crosses <- expand.grid(Female = females, Male = males, stringsAsFactors = FALSE)
  crosses$HybridID <- sprintf("H%02d", seq_len(nrow(crosses)))
  effect <- stats::rnorm(ncol(geno))
  gca <- as.vector(scale(geno %*% effect))
  names(gca) <- parents
  crosses$Yield <- gca[crosses$Female] + gca[crosses$Male] + stats::rnorm(nrow(crosses), sd = 0.5)
  unobserved <- c("H03", "H17", "H30", "H44")
  crosses$Yield[crosses$HybridID %in% unobserved] <- NA
  pheno <- crosses[, c("HybridID", "Female", "Male", "Yield")]

  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  res <- suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = geno, gen_name = "HybridID",
    female_parent = "Female", male_parent = "Male", response = "Yield", ploidy = 2L,
    hybrid_bayes = TRUE, GS_model_cv = "RKHS", gmatrix_method = "VanRaden",
    cross_validation = TRUE, cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_One_New_Parent", nfolds = 2L, replication = 1L,
    random_state = 1L, nIter = 300L, burnIn = 100L, thin = 5L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    system_database = FALSE, message = FALSE, verbose = FALSE, parallel_mode = "sequential"
  ))
  processed <- res$cv_results_processed
  pred <- processed$hybrid_cv_predictions

  expect_gt(nrow(pred), 0L)
  expect_false(any(pred$HybridID %in% unobserved))
  expect_true(all(is.finite(pred$Observed_value)))
  # One_New_Parent holds each observed hybrid out once as a new female and
  # once as a new male
  n_observed <- sum(!is.na(pheno$Yield))
  expect_equal(sum(processed$hybrid_prediction_counts$n_test_hybrids), nrow(pred))
  expect_equal(nrow(pred), 2L * n_observed)
})
