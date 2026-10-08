if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

collision_runtime_data <- function(family) {
  set.seed(20260920)
  ids <- paste0("g", seq_len(30L))
  blocks <- lapply(seq_len(4L), function(i) {
    matrix(rnorm(30L * 5L), 30L, dimnames = list(ids, paste0("x", 1:5)))
  })
  kernels <- lapply(blocks, function(x) tcrossprod(x) / ncol(x) + diag(0.1, nrow(x)))
  names(kernels) <- c("genomic", "RNA-A", "RNA_A", "RNA_A_1")
  y <- switch(family,
    gaussian = blocks[[1L]][, 1L] + blocks[[2L]][, 2L],
    binary = factor(rep(c("no", "yes"), 15L), levels = c("no", "yes")),
    multiclass = factor(rep(c("A", "B", "C"), 10L), levels = c("A", "B", "C")),
    ordinal = ordered(rep(c("low", "medium", "high"), 10L), levels = c("low", "medium", "high"))
  )
  list(pheno = data.frame(GID = ids, Trait = y), kernels = kernels, omics = blocks[3:4])
}

collision_prediction_table <- function(x) {
  if (is.data.frame(x)) {
    if ("Train_Test_Label" %in% names(x) &&
        any(c("Predicted_value", "Predicted_class") %in% names(x))) return(x)
    return(NULL)
  }
  if (is.list(x)) {
    for (item in x) {
      found <- collision_prediction_table(item)
      if (!is.null(found)) return(found)
    }
  }
  NULL
}

for (family in c("gaussian", "binary", "multiclass", "ordinal")) {
  test_that(paste("colliding kernel names survive direct prediction and CV for", family), {
    dat <- collision_runtime_data(family)
    model <- switch(family, gaussian = "Ridge_Regression", ordinal = "RKHS", "RandomForest")
    if (family == "ordinal") {
      skip_if_not_installed("BGLR")
    } else {
      py <- gp_preferred_python(purpose = "ml")
      skip_if(is.null(py) || !file.exists(py), "No ML Python runtime available")
      withr::local_envvar(c(PREDICTPRO_PYTHON = py, PREDICTPRO_ML_BRIDGE_BACKEND = "cli"))
    }
    prepare <- gp_prepare_model_input_objects
    observed <- list()
    local_mocked_bindings(
      gp_prepare_model_input_objects = function(ctx) {
        result <- prepare(ctx)
        observed[[length(observed) + 1L]] <<- result
        result
      },
      .package = "PredictProR"
    )
    args <- list(
      pheno_data = dat$pheno, kernel_list = dat$kernels,
      response = "Trait", gen_name = "GID", fixed = ~ 1, random = ~ GID,
      response_family = family, GS_model = model, GS_model_cv = model,
      random_state = 291L, para_tunning = FALSE, n_bootstrap = 2L,
      nIter = 80L, burnIn = 20L, thin = 2L,
      system_database = TRUE, message = FALSE,
      cross_validation_meth = "Stratified_K-Folds", nfolds = 2L,
      replication = 1L, cv_evaluation_only = TRUE,
      eval_metrics = if (family == "gaussian") "mean_squared_error" else if (family == "ordinal") "mean_absolute_error_class" else "log_loss"
    )
    if (family != "ordinal") {
      args$omic1_data <- dat$omics[[1L]]
      args$omic2_data <- dat$omics[[2L]]
    }
    for (cv in c(FALSE, TRUE)) {
      observed <- list()
      args$cross_validation <- cv
      args$pheno_data <- dat$pheno
      if (!cv) args$pheno_data$Trait[28:30] <- NA
      out <- do.call(PredictProR::model_execute, args)
      expect_gt(length(observed), 0L)
      for (prepared in observed) {
        expect_length(prepared$gmatrix_kernel_model_ready_list, 4L)
        if (family != "ordinal") {
          expect_equal(nrow(prepared$ml_dat_res$kernel_feature_summary), 4L)
          # Both raw omics blocks must remain alongside all four kernel blocks.
          expect_true(all(c("omic1_model_ready", "omic2_model_ready") %in%
                            names(prepared$geno_omic_model_ready_list)))
        }
      }
      if (cv) {
        expect_gt(length(out$cv_results_raw), 0L)
        raw <- out$cv_results_raw[[1L]]
        expect_true(any(raw$ypred_cv_Reps_all$cv_role == "test"))
        expect_false(anyNA(raw$ypred_cv_Reps_all$yhat))
        if (family == "gaussian") {
          expect_true(all(is.finite(raw$ypred_cv_Reps_all$yhat)))
        }
        if (family != "gaussian") {
          prob <- as.matrix(raw$yprob_cv_Reps_all[, levels(dat$pheno$Trait), drop = FALSE])
          expect_true(all(is.finite(prob) & prob >= 0 & prob <= 1))
          expect_equal(unname(rowSums(prob)), rep(1, nrow(prob)), tolerance = 1e-6)
        }
      } else {
        pred <- collision_prediction_table(out)
        expect_s3_class(pred, "data.frame")
        expect_equal(sum(pred$Train_Test_Label == "Test"), 3L)
        if (family == "gaussian") {
          expect_true(all(is.finite(pred$Predicted_value)))
        } else {
          prob <- as.matrix(pred[, grep("^Probability_", names(pred)), drop = FALSE])
          expect_equal(ncol(prob), nlevels(dat$pheno$Trait))
          expect_true(all(is.finite(prob) & prob >= 0 & prob <= 1))
          expect_equal(unname(rowSums(prob)), rep(1, nrow(prob)), tolerance = 1e-6)
        }
      }
    }
  })
}

for (family in c("gaussian", "binary", "multiclass")) {
  test_that(paste("deep learning retains colliding kernels for", family), {
    py <- gp_preferred_python(purpose = "dl")
    skip_if(is.null(py) || !file.exists(py), "No DL Python runtime available")
    withr::local_envvar(c(PREDICTPRO_DL_PYTHON = py, PREDICTPRO_DL_DEVICE = "cpu"))
    dat <- collision_runtime_data(family)
    dat$pheno$Trait[28:30] <- NA
    prepare <- gp_prepare_met_kernel_features
    observed <- NULL
    local_mocked_bindings(
      gp_prepare_met_kernel_features = function(...) {
        observed <<- prepare(...)
        observed
      },
      .package = "PredictProR"
    )
    expect_warning(out <- PredictProR::model_execute(
      pheno_data = dat$pheno, kernel_list = dat$kernels,
      omic1_data = dat$omics[[1L]], omic2_data = dat$omics[[2L]],
      response = "Trait", gen_name = "GID", response_family = family,
      GS_model = "DenseNeuralNet", epochs = 2L, batch_size = 8L,
      mlp_neurons_per_layer = 8L, device = "cpu", n_bootstrap = 2L,
      random_state = 291L, random_seed = 291L,
      para_tunning = FALSE, system_database = TRUE, message = FALSE
    ), NA)
    expect_length(observed$feature_blocks, 4L)
    expect_identical(anyDuplicated(names(observed$feature_table)), 0L)
    pred <- collision_prediction_table(out)
    expect_s3_class(pred, "data.frame")
    expect_equal(sum(pred$Train_Test_Label == "Test"), 3L)
    if (family == "gaussian") {
      expect_true(all(is.finite(pred$Predicted_value)))
    } else {
      prob <- as.matrix(pred[, grep("^Probability_", names(pred)), drop = FALSE])
      expect_equal(ncol(prob), nlevels(dat$pheno$Trait))
      expect_true(all(is.finite(prob) & prob >= 0 & prob <= 1))
      expect_equal(unname(rowSums(prob)), rep(1, nrow(prob)), tolerance = 1e-6)
    }
  })
}
