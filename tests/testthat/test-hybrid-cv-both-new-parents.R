test_that("Both_New_Parents CV trains kernel models without either parent of the held-out cross", {
  skip_if_not_installed("BGLR")
  # The scenario's training set is the hybrids sharing NEITHER parent with the
  # held-out cross (train_ids). ML/DL honoured it, but the ASReml, Bayesian and
  # GP paths masked only the held-out cross and kept training on its parents'
  # other hybrids, so for those models the parents were not new.
  set.seed(7)
  females <- sprintf("F%d", 1:5); males <- sprintf("M%d", 1:4)
  parents <- c(females, males)
  geno <- matrix(sample(c(0, 2), length(parents) * 40, replace = TRUE), length(parents), 40,
                 dimnames = list(parents, sprintf("SNP%02d", 1:40)))
  crosses <- expand.grid(Female = females, Male = males, stringsAsFactors = FALSE)
  crosses$HybridID <- sprintf("H%02d", seq_len(nrow(crosses)))
  crosses$Yield <- stats::rnorm(nrow(crosses))
  pheno <- crosses[, c("HybridID", "Female", "Male", "Yield")]

  fold_checks <- list()
  real_fit <- PredictProR:::gp_hybrid_bayes_gaussian_model
  local_mocked_bindings(
    gp_hybrid_bayes_gaussian_model = function(pheno_object, ...) {
      ph <- as.data.frame(pheno_object)
      train <- ph[!is.na(ph$Yield), , drop = FALSE]
      masked <- ph[is.na(ph$Yield), , drop = FALSE]
      # some masked hybrid must have BOTH parents absent from training
      fold_checks[[length(fold_checks) + 1L]] <<- any(
        !masked$Female %in% train$Female & !masked$Male %in% train$Male
      )
      real_fit(pheno_object = pheno_object, ...)
    },
    .package = "PredictProR"
  )

  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  res <- suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = geno, gen_name = "HybridID",
    female_parent = "Female", male_parent = "Male", response = "Yield", ploidy = 2L,
    hybrid_bayes = TRUE, GS_model_cv = "GBLUP_BRR", gmatrix_method = "VanRaden",
    cross_validation = TRUE, cv_evaluation_only = TRUE,
    cross_validation_meth = "Hybrid_Both_New_Parents", replication = 1L,
    random_state = 1L, nIter = 200L, burnIn = 50L, thin = 5L,
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    system_database = FALSE, message = FALSE, verbose = FALSE, parallel_mode = "sequential"
  ))

  expect_gt(length(fold_checks), 0L)
  expect_true(all(unlist(fold_checks)))
  expect_gt(nrow(res$cv_results_processed$hybrid_cv_predictions), 0L)
})
