if (!identical(Sys.getenv("PREDICTPROR_RUN_USER_CONTRACT_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_USER_CONTRACT_TESTS=1 to run user-facing output contract tests.")
}

contract_gaussian_cols <- PredictProR:::gp_gaussian_prediction_columns()

contract_classification_base_cols <- c(
  "GID", "Predicted_class", "Train_Test_Label", "Observed_class",
  "Prediction_confidence", "Classification_uncertainty",
  "Reliability", "Reliability_remarks"
)

contract_model_result_names <- c(
  "model_parameters", "predicted_values", "diagnostic_plots", "variance_components"
)

make_contract_gaussian_data <- function(n = 18L, p = 10L, seed = 20260518L) {
  set.seed(seed)
  X <- matrix(sample(c(0, 1, 2), n * p, replace = TRUE), nrow = n)
  for (j in seq_len(ncol(X))) {
    if (length(unique(X[, j])) < 2L) {
      X[, j] <- (seq_len(n) + j) %% 3L
    }
  }
  rownames(X) <- paste0("g", seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))
  y <- 2 + 0.7 * X[, 1] - 0.4 * X[, 2] + rnorm(n, sd = 0.2)
  y[sample(seq_len(n), 4L)] <- NA_real_
  list(
    pheno = data.frame(GID = rownames(X), Yield = y, stringsAsFactors = FALSE),
    geno = X
  )
}

make_contract_class_data <- function(kind = c("binary", "ordinal", "multiclass"),
                                     n = 24L,
                                     p = 8L,
                                     seed = 20260519L) {
  kind <- match.arg(kind)
  set.seed(seed + match(kind, c("binary", "ordinal", "multiclass")))
  X <- matrix(sample(c(0, 1, 2), n * p, replace = TRUE), nrow = n)
  rownames(X) <- paste0(substr(kind, 1, 1), seq_len(n))
  colnames(X) <- paste0("m", seq_len(p))
  score <- X[, 1] - 0.6 * X[, 2] + 0.25 * X[, 3]
  y <- switch(
    kind,
    binary = factor(ifelse(score > stats::median(score), "yes", "no"), levels = c("no", "yes")),
    ordinal = ordered(
      cut(score, breaks = stats::quantile(score, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE,
          labels = c("low", "medium", "high")),
      levels = c("low", "medium", "high")
    ),
    multiclass = factor(
      ifelse(score < stats::quantile(score, 1 / 3), "A",
             ifelse(score < stats::quantile(score, 2 / 3), "B", "C")),
      levels = c("A", "B", "C")
    )
  )
  y[sample(seq_len(n), 5L)] <- NA
  list(
    pheno = data.frame(GID = rownames(X), Class = y, stringsAsFactors = FALSE),
    geno = X
  )
}

make_contract_met_data <- function(n_gid = 8L, p = 8L, seed = 20260520L) {
  set.seed(seed)
  gids <- paste0("g", seq_len(n_gid))
  envs <- c("E1", "E2")
  X <- matrix(sample(c(0, 1, 2), n_gid * p, replace = TRUE), nrow = n_gid)
  rownames(X) <- gids
  colnames(X) <- paste0("m", seq_len(p))
  pheno <- expand.grid(GID = gids, Env = envs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  gid_i <- match(pheno$GID, gids)
  env_eff <- c(E1 = -0.2, E2 = 0.35)[pheno$Env]
  pheno$Yield <- 3 + 0.6 * X[gid_i, 1] - 0.25 * X[gid_i, 2] + env_eff + rnorm(nrow(pheno), sd = 0.15)
  pheno$Yield[sample(seq_len(nrow(pheno)), 4L)] <- NA_real_
  list(pheno = pheno, geno = X)
}

make_contract_hybrid_data <- function(n_parent = 6L, p = 6L, seed = 20260521L) {
  set.seed(seed)
  females <- paste0("F", seq_len(n_parent))
  males <- paste0("M", seq_len(n_parent))
  hybrids <- data.frame(
    Female = rep(females, each = 2L),
    Male = rep(males, length.out = 2L * n_parent),
    stringsAsFactors = FALSE
  )
  hybrids$HybridID <- paste(hybrids$Female, hybrids$Male, sep = "_")
  parent_ids <- unique(c(females, males))
  geno <- matrix(
    sample(c(0, 1, 2), length(parent_ids) * p, replace = TRUE),
    nrow = length(parent_ids)
  )
  rownames(geno) <- parent_ids
  colnames(geno) <- paste0("m", seq_len(p))
  f_eff <- rowMeans(geno[hybrids$Female, , drop = FALSE])
  m_eff <- rowMeans(geno[hybrids$Male, , drop = FALSE])
  hybrids$Yield <- 5 + 0.7 * f_eff + 0.4 * m_eff + rnorm(nrow(hybrids), sd = 0.1)
  hybrids$Yield[sample(seq_len(nrow(hybrids)), 3L)] <- NA_real_
  hybrid_geno <- (geno[hybrids$Female, , drop = FALSE] + geno[hybrids$Male, , drop = FALSE]) / 2
  rownames(hybrid_geno) <- hybrids$HybridID
  list(pheno = hybrids, geno = hybrid_geno)
}

base_contract_args <- function(dat, response = "Yield", model = "BRR", response_family = "gaussian") {
  list(
    pheno_data = dat$pheno,
    geno_data = dat$geno,
    gen_name = "GID",
    response = response,
    random = ~GID,
    GS_model = model,
    response_family = response_family,
    cross_validation = FALSE,
    system_database = TRUE,
    message = FALSE,
    para_tunning = FALSE,
    gmatrix_method = "Yang",
    kernel_method = "Gaussian_kernel",
    ind_call_rate_threshold = 0,
    snp_call_rate_threshold = 0,
    maf_threshold = 0,
    het_threshold = 1,
    ld_prunning_qc = FALSE,
    impute = "mean",
    nIter = 80L,
    burnIn = 20L,
    thin = 2L,
    n_bootstrap = 2L,
    ntree = 20L,
    maxnodes = 10L,
    nodesize = 1L,
    ncomp = 2L,
    catboost_iterations = 10L,
    catboost_depth = 3L,
    lightgbm_nrounds = 10L,
    lightgbm_num_leaves = 7L,
    epochs = 1L,
    batch_size = 8L,
    validation_split = 0.2,
    compile_model = FALSE,
    deterministic = TRUE,
    random_seed = 11L,
    device = "cpu",
    use_amp = FALSE,
    batch_norm = FALSE,
    mlp_neurons_per_layer = as.integer(c(12L, 6L)),
    mlp_learning_rate = 1e-3
  )
}

run_contract_model <- function(args) {
  do.call(PredictProR::model_execute, args)
}

extract_contract_result <- function(out, trait) {
  testthat::expect_true(is.list(out$model_results_by_model))
  testthat::expect_true(is.list(out$model_results_by_trait))
  testthat::expect_false(trait %in% names(out))
  model_names <- names(out$model_results_by_model)
  testthat::expect_equal(length(model_names), 1L)
  res <- out$model_results_by_model[[model_names[[1L]]]][[trait]]
  testthat::expect_true(is.list(res))
  res
}

expect_one_public_prediction_table <- function(res, optional_names = character()) {
  testthat::expect_identical(
    names(res)[seq_along(contract_model_result_names)],
    contract_model_result_names
  )
  extra_names <- setdiff(names(res), contract_model_result_names)
  if (length(optional_names)) {
    testthat::expect_true(all(optional_names %in% extra_names))
  } else {
    testthat::expect_length(extra_names, 0L)
  }
  testthat::expect_true(is.data.frame(res$predicted_values))
  testthat::expect_false("Predicted_value" %in% names(res))
  testthat::expect_false("diagnostic_tst_plot" %in% names(res))
  testthat::expect_false("Variance_components" %in% names(res))
}

expect_variance_component_contract <- function(res) {
  testthat::expect_true(is.data.frame(res$variance_components))
  testthat::expect_gt(nrow(res$variance_components), 0L)
  vc_names <- names(res$variance_components)
  testthat::expect_true(
    all(c("Component", "Components", "Standard_error") %in% vc_names) ||
      all(c("Trait", "Component", "Components", "Standard_error") %in% vc_names) ||
      all(c("Trait", "Env", "Component", "Components", "Standard_error") %in% vc_names)
  )
  testthat::expect_false("Estimation_method" %in% vc_names)
  testthat::expect_true(any(nzchar(as.character(res$variance_components$Component))))
  testthat::expect_true(any(is.finite(suppressWarnings(as.numeric(res$variance_components$Components)))))
}

expect_gaussian_contract <- function(res, include_env = FALSE, include_trait = FALSE) {
  optional_names <- if (isTRUE(include_env)) {
    c(
      "across_environment_predicted_values",
      "covariance_components",
      "correlation_components"
    )
  } else {
    character()
  }
  expect_one_public_prediction_table(res, optional_names = optional_names)
  expected <- c("GID")
  if (isTRUE(include_env)) expected <- c(expected, "Env")
  if (isTRUE(include_trait)) expected <- c(expected, "Trait")
  expected <- c(expected, contract_gaussian_cols[-1L])
  testthat::expect_identical(names(res$predicted_values), expected)
  train <- as.character(res$predicted_values$Train_Test_Label) == "Train"
  if (any(train, na.rm = TRUE)) {
    testthat::expect_true(all(is.finite(
      suppressWarnings(as.numeric(res$predicted_values$Observed_value[train]))
    )))
  }
  expect_variance_component_contract(res)
}

expect_classification_contract <- function(res) {
  expect_one_public_prediction_table(res)
  testthat::expect_true(all(contract_classification_base_cols %in% names(res$predicted_values)))
  testthat::expect_false("Standard_error" %in% names(res$predicted_values))
  testthat::expect_false("PEV" %in% names(res$predicted_values))
  testthat::expect_true(any(grepl("^Probability_", names(res$predicted_values))))
}

single_env_models_to_check <- function() {
  defaults <- c("BRR", "BayesA", "BayesB", "BayesC", "BL", "GBLUP_BRR", "RKHS",
                "GBLUP", "Kernel-GBLUP", "Gaussian-Process-GBLUP",
                "RandomForest", "Ridge_Regression", "Lasso", "PartialLeastSquare",
                "SupportVectorMachine", "K-NearestNeighbors", "Xgboost", "CatBoost",
                "LightGBM", "Conv1DNet", "TabNet", "TabTransformer", "DenseNeuralNet")
  requested <- Sys.getenv("PREDICTPROR_CONTRACT_MODELS", "")
  if (nzchar(requested)) {
    trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  } else {
    defaults
  }
}

test_that("installed model_execute returns one compact single-environment Gaussian table across model families", {
  dat <- make_contract_gaussian_data()
  failures <- list()
  for (model in single_env_models_to_check()) {
    if (identical(model, "GBLUP") && !suppressWarnings(requireNamespace("asreml", quietly = TRUE))) {
      next
    }
    args <- base_contract_args(dat, model = model)
    out <- tryCatch(run_contract_model(args), error = function(e) e)
    if (inherits(out, "error")) {
      failures[[model]] <- conditionMessage(out)
      next
    }
    err <- tryCatch({
      res <- extract_contract_result(out, "Yield")
      expect_gaussian_contract(res)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) {
      failures[[model]] <- err
    }
  }
  testthat::expect_equal(failures, list())
})

test_that("installed model_execute infers binary, ordinal, and multiclass contracts", {
  binary <- make_contract_class_data("binary")
  out_binary <- run_contract_model(modifyList(
    base_contract_args(binary, response = "Class", model = "RandomForest", response_family = "auto"),
    list(gen_name = "GID")
  ))
  expect_classification_contract(extract_contract_result(out_binary, "Class"))

  ordinal <- make_contract_class_data("ordinal")
  out_ord <- run_contract_model(modifyList(
    base_contract_args(ordinal, response = "Class", model = "BayesA", response_family = "auto"),
    list(gen_name = "GID")
  ))
  expect_classification_contract(extract_contract_result(out_ord, "Class"))

  multiclass <- make_contract_class_data("multiclass")
  out_multi <- run_contract_model(modifyList(
    base_contract_args(multiclass, response = "Class", model = "CatBoost", response_family = "auto"),
    list(gen_name = "GID")
  ))
  expect_classification_contract(extract_contract_result(out_multi, "Class"))
})

test_that("installed MET and hybrid paths use task-specific public schemas", {
  met <- make_contract_met_data()
  out_met <- run_contract_model(modifyList(
    base_contract_args(met, model = "RKHS"),
    list(heter_groups = "Env", heter_resid = TRUE)
  ))
  met_res <- extract_contract_result(out_met, "Yield")
  expect_gaussian_contract(met_res, include_env = TRUE)
  expect_false(any(grepl(
    "not_identifiable$",
    as.character(met_res$variance_components$Component)
  )))
  expect_true(any(grepl(
    "^genetic_variance",
    as.character(met_res$variance_components$Component)
  )))
  expect_true(any(grepl(
    "^residual_variance",
    as.character(met_res$variance_components$Component)
  )))
  expect_true(any(grepl(
    "^heritability",
    as.character(met_res$variance_components$Component)
  )))
  if ("across_environment_predicted_values" %in% names(met_res)) {
    testthat::expect_identical(names(met_res$across_environment_predicted_values), contract_gaussian_cols)
  }

  hybrid <- make_contract_hybrid_data()
  out_hybrid <- run_contract_model(modifyList(
    base_contract_args(hybrid, model = "RandomForest"),
    list(
      gen_name = "HybridID",
      random = ~HybridID,
      hybrid_ml = TRUE,
      hybrid_id = "HybridID",
      female_parent = "Female",
      male_parent = "Male"
    )
  ))
  hres <- extract_contract_result(out_hybrid, "Yield")
  expect_one_public_prediction_table(
    hres,
    optional_names = "prediction_uncertainty_source"
  )
  testthat::expect_identical(
    names(hres$predicted_values),
    c("Hybrid_ID", "Female", "Male", contract_gaussian_cols[-1L])
  )
})

test_that("multi-model MET keeps fitted results aligned with model and trait labels", {
  met <- make_contract_met_data(n_gid = 10L, seed = 20260522L)
  met$pheno$Yield2 <- 1.4 + 0.65 * met$pheno$Yield +
    stats::rnorm(nrow(met$pheno), sd = 0.1)
  met$pheno$Yield2[c(2L, 7L, 12L, 18L)] <- NA_real_

  args <- base_contract_args(met, response = c("Yield", "Yield2"), model = NULL)
  args$GS_model <- c("RandomForest", "RKHS")
  args$heter_groups <- "Env"
  args$heter_resid <- TRUE
  args$met_ml_dl <- TRUE
  args$sequential_models <- "RKHS"
  args$parallel_mode <- "sequential"
  args$num_cores <- 1L

  out <- run_contract_model(args)
  expect_true(all(c("RandomForest", "RKHS") %in% names(out$model_results_by_model)))

  rf <- out$model_results_by_model[["RandomForest"]][["Yield"]]
  rkhs <- out$model_results_by_model[["RKHS"]][["Yield2"]]
  expect_true(is.list(rf))
  expect_true(is.list(rkhs))

  rf_requested <- rf$model_parameters$summary[
    rf$model_parameters$stat == "requested_model"
  ]
  rkhs_model <- rkhs$model_parameters$summary[
    rkhs$model_parameters$stat == "bayesian_kernel_model"
  ]
  expect_identical(as.character(rf_requested[[1L]]), "RandomForest")
  expect_identical(as.character(rkhs_model[[1L]]), "RKHS")

  expect_true(any(grepl(
    "not_identifiable$",
    as.character(rf$variance_components$Component)
  )))
  expect_false(any(grepl(
    "not_identifiable$",
    as.character(rkhs$variance_components$Component)
  )))
  expect_true(any(grepl(
    "^genetic_variance",
    as.character(rkhs$variance_components$Component)
  )))
  expect_true(any(grepl(
    "^residual_variance",
    as.character(rkhs$variance_components$Component)
  )))
})
