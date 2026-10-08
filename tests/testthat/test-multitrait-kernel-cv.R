multitrait_kernel_cv_fixture <- function() {
  data.frame(
    Sample = paste0("G", seq_len(10)),
    Trait1 = seq_len(10) + c(0.1, -0.1),
    Trait2 = c(NA, seq_len(8) / 2, NA),
    Trait3 = c(seq_len(7) / 3, NA, 3.1, 3.4),
    stringsAsFactors = FALSE
  )
}

test_that("joint kernel CV masks every trait by genotype without leakage", {
  ph <- multitrait_kernel_cv_fixture()
  traits <- c("Trait1", "Trait2", "Trait3")
  calls <- list()
  fake_fit <- function(pheno_masked, test_ids, fold_seed, rep_idx, fold_idx) {
    test_rows <- match(test_ids, pheno_masked$Sample)
    expect_true(all(is.na(as.matrix(pheno_masked[test_rows, traits, drop = FALSE]))))
    train_rows <- setdiff(seq_len(nrow(ph)), test_rows)
    expect_equal(
      pheno_masked[train_rows, traits, drop = FALSE],
      ph[train_rows, traits, drop = FALSE]
    )
    calls[[length(calls) + 1L]] <<- list(
      ids = test_ids,
      rep = rep_idx,
      fold = fold_idx,
      seed = fold_seed
    )
    grid <- expand.grid(
      Sample = test_ids,
      Trait = traits,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid$Predicted_value <- seq_len(nrow(grid)) / 10 + rep_idx + fold_idx
    grid$Standard_error <- 0.2
    grid$PEV <- 0.04
    grid
  }

  out <- PredictProR:::gp_multitrait_joint_gaussian_cv(
    pheno_object = ph,
    response = traits,
    gen_name = "Sample",
    model_type = "RKHS",
    fit_fold = fake_fit,
    cross_validation_meth = "Repeated_K-Folds",
    nfolds = 2L,
    replication = 2L,
    eval_metrics = c("mean_squared_error", "bias"),
    random_seed = 707L
  )

  pred <- out$cv_results_processed$multitrait_oof_predictions
  expected_per_rep <- sum(is.finite(as.matrix(ph[, traits, drop = FALSE])))
  expect_equal(nrow(pred), 2L * expected_per_rep)
  expect_false(anyDuplicated(paste(pred$rep, pred$GID, pred$Trait, sep = "\r")) > 0L)
  expect_true(all(pred$Train_Test_Label == "Test"))
  expect_true(all(is.finite(pred$Predicted_value)))
  expect_true(all(pred$PEV == 0.04))
  expect_true(isTRUE(out$cv_results_processed$cv_scheme$all_traits_masked_together))
  expect_equal(length(calls), 4L)

  assignments <- out$cv_results_processed$fold_assignments
  expect_false(anyDuplicated(paste(assignments$rep, assignments$GID, sep = "\r")) > 0L)
  expect_equal(
    sort(unique(assignments$GID)),
    sort(ph$Sample[rowSums(is.finite(as.matrix(ph[, traits, drop = FALSE]))) > 0L])
  )
})

test_that("joint kernel CV rejects missing, duplicate, and non-finite fold predictions", {
  ph <- multitrait_kernel_cv_fixture()
  traits <- c("Trait1", "Trait2", "Trait3")
  make_predictions <- function(test_ids) {
    grid <- expand.grid(
      GID = test_ids,
      Trait = traits,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid$Predicted_value <- seq_len(nrow(grid))
    grid
  }
  run <- function(kind) {
    PredictProR:::gp_multitrait_joint_gaussian_cv(
      pheno_object = ph,
      response = traits,
      gen_name = "Sample",
      model_type = "Gaussian-Process-GBLUP",
      fit_fold = function(pheno_masked, test_ids, fold_seed, rep_idx, fold_idx) {
        pred <- make_predictions(test_ids)
        if (identical(kind, "missing")) pred <- pred[-1L, , drop = FALSE]
        if (identical(kind, "duplicate")) pred <- rbind(pred, pred[1L, , drop = FALSE])
        if (identical(kind, "nonfinite")) pred$Predicted_value[[1L]] <- NA_real_
        pred
      },
      nfolds = 2L,
      replication = 1L,
      random_seed = 19L
    )
  }

  expect_error(run("missing"), "did not return every held-out observed")
  expect_error(run("duplicate"), "duplicate genotype-trait predictions")
  expect_error(run("nonfinite"), "non-finite held-out predictions")
})

test_that("joint kernel CV fold planning fails loudly for traits that are too sparse", {
  ph <- multitrait_kernel_cv_fixture()
  ph$Trait3 <- NA_real_
  ph$Trait3[1:3] <- 1:3
  expect_error(
    PredictProR:::gp_multitrait_joint_cv_fold_plan(
      pheno_data = ph,
      response = c("Trait1", "Trait2", "Trait3"),
      gen_name = "Sample",
      nfolds = 2L
    ),
    "could not create leakage-free genotype folds|Too-sparse traits"
  )
})

multitrait_met_kernel_cv_fixture <- function() {
  grid <- expand.grid(
    Sample = paste0("G", seq_len(6)),
    Env = paste0("E", seq_len(3)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$Trait1 <- seq_len(nrow(grid)) / 3
  grid$Trait2 <- seq_len(nrow(grid)) / 5 + 1
  grid
}

test_that("joint MT-MET kernel CV uses CV0 CV1 and CV2 row semantics", {
  ph <- multitrait_met_kernel_cv_fixture()
  traits <- c("Trait1", "Trait2")
  run_scenario <- function(method) {
    PredictProR:::gp_multitrait_joint_gaussian_cv(
      pheno_object = ph,
      response = traits,
      gen_name = "Sample",
      heter_groups = "Env",
      model_type = "RKHS",
      fit_fold = function(pheno_masked, test_ids, fold_seed, rep_idx, fold_idx) {
        masked <- apply(is.na(pheno_masked[, traits, drop = FALSE]), 1L, all)
        expect_true(any(masked))
        grid <- do.call(rbind, lapply(traits, function(trait) {
          data.frame(
            Sample = pheno_masked$Sample,
            Env = pheno_masked$Env,
            Trait = trait,
            Predicted_value = seq_len(nrow(pheno_masked)) / 10,
            Standard_error = 0.2,
            PEV = 0.04,
            stringsAsFactors = FALSE
          )
        }))
        grid
      },
      cross_validation_meth = method,
      nfolds = if (identical(method, "CV0")) 3L else 2L,
      replication = 1L,
      eval_metrics = "mean_squared_error",
      random_seed = 31L
    )
  }

  cv0 <- run_scenario("CV0")
  cv1 <- run_scenario("CV1")
  cv2 <- run_scenario("CV2")

  expect_identical(cv0$cv_results_processed$cv_scheme$unit, "environment")
  expect_identical(cv1$cv_results_processed$cv_scheme$unit, "genotype")
  expect_identical(cv2$cv_results_processed$cv_scheme$unit, "genotype_environment_cell")
  expect_equal(nrow(cv0$cv_results_processed$multitrait_oof_predictions), nrow(ph) * length(traits))
  expect_equal(nrow(cv1$cv_results_processed$multitrait_oof_predictions), nrow(ph) * length(traits))
  expect_equal(nrow(cv2$cv_results_processed$multitrait_oof_predictions), nrow(ph) * length(traits))

  a0 <- cv0$cv_results_processed$fold_assignments
  a1 <- cv1$cv_results_processed$fold_assignments
  a2 <- cv2$cv_results_processed$fold_assignments
  expect_true(all(vapply(split(a0$fold, a0$Env), function(x) length(unique(x)) == 1L, logical(1L))))
  expect_true(all(vapply(split(a1$fold, a1$GID), function(x) length(unique(x)) == 1L, logical(1L))))
  expect_true(all(vapply(split(a2$fold, a2$GID), function(x) length(unique(x)) > 1L, logical(1L))))
})

test_that("joint MT-MET kernel CV rejects single-environment K-fold names", {
  ph <- multitrait_met_kernel_cv_fixture()
  expect_error(
    PredictProR:::gp_multitrait_joint_gaussian_cv(
      pheno_object = ph,
      response = c("Trait1", "Trait2"),
      gen_name = "Sample",
      heter_groups = "Env",
      model_type = "RKHS",
      fit_fold = function(...) data.frame(),
      cross_validation_meth = "K-Folds",
      nfolds = 2L
    ),
    "MT-MET.*CV0"
  )
})

test_that("joint GP and Bayesian wrappers consume the shared MT-MET scenarios", {
  ph <- multitrait_met_kernel_cv_fixture()
  traits <- c("Trait1", "Trait2")
  K <- diag(6)
  rownames(K) <- colnames(K) <- paste0("G", seq_len(6))

  prediction_grid <- function(masked) {
    do.call(rbind, lapply(traits, function(trait) {
      data.frame(
        GID = masked$Sample,
        Env = masked$Env,
        Trait = trait,
        Predicted_value = seq_len(nrow(masked)) / 10,
        Standard_error = 0.2,
        PEV = 0.04,
        stringsAsFactors = FALSE
      )
    }))
  }

  testthat::local_mocked_bindings(
    gp_multi_trait_met_model = function(pheno_data, ...) {
      list(predictions = prediction_grid(pheno_data))
    },
    bayes_multitrait_joint_fit = function(pheno_data, ...) {
      list(bayes_result = list(predicted_values = prediction_grid(pheno_data)))
    },
    .package = "PredictProR"
  )

  gp <- PredictProR:::gp_multitrait_gp_gaussian_cv(
    model_type = "GP",
    pheno_object = ph,
    response = traits,
    gen_name = "Sample",
    heter_groups = "Env",
    gmatrix = K,
    cross_validation_meth = "CV1",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )
  bayes <- PredictProR:::gp_multitrait_bayes_gaussian_cv(
    model_type = "RKHS",
    pheno_object = ph,
    response = traits,
    gen_name = "Sample",
    heter_groups = "Env",
    kernels = list(genomic = K),
    bayes_para = list(nIter = 100L, burnIn = 20L, thin = 2L),
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )

  expect_identical(gp$cv_results_processed$cv_scheme$method, "cv1")
  expect_identical(bayes$cv_results_processed$cv_scheme$method, "cv2")
  expect_equal(nrow(gp$cv_results_processed$multitrait_oof_predictions), nrow(ph) * length(traits))
  expect_equal(nrow(bayes$cv_results_processed$multitrait_oof_predictions), nrow(ph) * length(traits))
})
