if (!identical(Sys.getenv("PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS"), "1")) {
  testthat::skip("Set PREDICTPROR_RUN_RUNTIME_SMOKE_TESTS=1 to run backend runtime smoke tests.")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bayes_heter_toy_data <- function(n_ids = 8L, n_markers = 8L) {
  set.seed(12)
  ids <- paste0("g", seq_len(n_ids))
  envs <- c("E1", "E2")
  pheno <- expand.grid(
    GID = ids,
    Env = envs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  marker <- matrix(
    stats::rnorm(n_ids * n_markers),
    nrow = n_ids,
    dimnames = list(ids, paste0("m", seq_len(n_markers)))
  )
  marker_signal <- scale(marker[, 1] - 0.4 * marker[, 2])[, 1]
  signal <- marker_signal[match(pheno$GID, ids)]
  pheno$Y1 <- 2 + 0.4 * signal + ifelse(pheno$Env == "E2", 0.8, 0) +
    stats::rnorm(nrow(pheno), sd = 0.25)
  pheno$Y2 <- 1 - 0.3 * signal + ifelse(pheno$Env == "E2", -0.4, 0) +
    stats::rnorm(nrow(pheno), sd = 0.25)
  pheno$Y1[c(3L, 12L)] <- NA_real_
  pheno$Y2[c(5L, 14L)] <- NA_real_

  k <- tcrossprod(scale(marker, center = TRUE, scale = FALSE))
  k <- k / mean(diag(k)) + diag(0.05, n_ids)
  rownames(k) <- colnames(k) <- ids
  k2 <- tcrossprod(scale(marker[, seq(2L, n_markers, by = 2L), drop = FALSE], center = TRUE, scale = FALSE))
  k2 <- k2 / mean(diag(k2)) + diag(0.05, n_ids)
  rownames(k2) <- colnames(k2) <- ids
  list(pheno = pheno, marker = marker, gmatrix = k, omic_kernel = k2)
}

bayes_heter_model_key <- function(model) {
  make.names(PredictProR:::gp_public_model_label(model)[1])
}

expect_met_gaussian_contract <- function(pred, test_n) {
  expect_s3_class(pred, "data.frame")
  expect_identical(
    names(pred),
    PredictProR:::gp_gaussian_prediction_columns(include_env = TRUE)
  )
  expect_equal(sum(pred$Train_Test_Label == "Test"), test_n)
  expect_true(any(is.finite(pred$Predicted_value)))
  expect_true(any(is.finite(pred$Standard_error)))
  expect_true(any(is.finite(pred$PEV)))
  finite <- is.finite(pred$PEV) & is.finite(pred$Standard_error)
  expect_equal(pred$PEV[finite], pred$Standard_error[finite] ^ 2, tolerance = 1e-10)
  expect_true(all(pred$Reliability[is.finite(pred$Reliability)] >= 0 &
                    pred$Reliability[is.finite(pred$Reliability)] <= 1))
  covered <- is.finite(pred$lower_bound) & is.finite(pred$upper_bound) & is.finite(pred$Predicted_value)
  expect_true(all(pred$lower_bound[covered] <= pred$Predicted_value[covered]))
  expect_true(all(pred$upper_bound[covered] >= pred$Predicted_value[covered]))
}

expect_across_environment_contract <- function(pred) {
  expect_s3_class(pred, "data.frame")
  expect_identical(names(pred), PredictProR:::gp_gaussian_prediction_columns())
  expect_true(any(is.finite(pred$Predicted_value)))
  expect_true(any(is.finite(pred$Standard_error)))
  expect_true(any(is.finite(pred$PEV)))
  finite <- is.finite(pred$PEV) & is.finite(pred$Standard_error)
  expect_equal(pred$PEV[finite], pred$Standard_error[finite] ^ 2, tolerance = 1e-10)
}

expect_variance_component_contract <- function(vc) {
  expect_s3_class(vc, "data.frame")
  expect_true(nrow(vc) > 0L)
  expect_true(all(c("Component", "Components", "Standard_error") %in% names(vc)))
  expect_true(all(nzchar(vc$Component)))
  expect_true(any(is.finite(vc$Components)))
}

expect_covariance_contract <- function(covariance, correlation, labels) {
  expect_true(is.matrix(covariance))
  expect_true(is.matrix(correlation))
  expect_identical(dim(covariance), c(length(labels), length(labels)))
  expect_identical(dim(correlation), c(length(labels), length(labels)))
  expect_equal(covariance, t(covariance), tolerance = 1e-10)
  expect_equal(correlation, t(correlation), tolerance = 1e-10)
  expect_identical(rownames(covariance), labels)
  expect_identical(colnames(covariance), labels)
  expect_true(min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  expect_equal(unname(diag(correlation)), rep(1, length(labels)), tolerance = 1e-10)
}

expect_covariance_se_contract <- function(covariance_se, correlation_se, labels) {
  expect_true(is.matrix(covariance_se))
  expect_true(is.matrix(correlation_se))
  expect_identical(dim(covariance_se), c(length(labels), length(labels)))
  expect_identical(dim(correlation_se), c(length(labels), length(labels)))
  expect_true(all(is.finite(covariance_se)))
  expect_true(all(is.finite(correlation_se)))
  expect_true(all(covariance_se >= 0))
  expect_true(all(correlation_se >= 0))
  expect_equal(unname(diag(correlation_se)), rep(0, length(labels)), tolerance = 1e-10)
}

fit_bayes_heter_model <- function(model, response = "Y1", data = bayes_heter_toy_data()) {
  PredictProR::model_execute(
    pheno_data = data$pheno,
    gmatrix = data$gmatrix,
    response = response,
    gen_name = "GID",
    random = ~ GID + GID:Env,
    fixed = ~ Env,
    GS_model = model,
    heter_groups = "Env",
    heter_resid = TRUE,
    response_family = "gaussian",
    para_tunning = FALSE,
    nIter = 60L,
    burnIn = 20L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )
}

test_that("GBLUP_BRR and RKHS support heterogeneous residual MET true prediction", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data()

  for (model in c("GBLUP_BRR", "RKHS")) {
    out <- fit_bayes_heter_model(model, response = "Y1", data = dat)
    trait_res <- out$model_results_by_model[[bayes_heter_model_key(model)]][["Y1"]]

    expect_met_gaussian_contract(trait_res$predicted_values, sum(is.na(dat$pheno$Y1)))
    expect_across_environment_contract(trait_res$across_environment_predicted_values)
    expect_variance_component_contract(trait_res$variance_components)
    expect_true(all(is.finite(trait_res$variance_components$Standard_error)))
    expect_covariance_contract(
      trait_res$Genetic_covariance_environments,
      trait_res$Genetic_correlation_environments,
      c("E1", "E2")
    )
    expect_covariance_contract(
      trait_res$Residual_covariance_environments,
      trait_res$Residual_correlation_environments,
      c("E1", "E2")
    )
    expect_covariance_se_contract(
      trait_res$Genetic_covariance_environments_SE,
      trait_res$Genetic_correlation_environments_SE,
      c("E1", "E2")
    )
    expect_covariance_se_contract(
      trait_res$Residual_covariance_environments_SE,
      trait_res$Residual_correlation_environments_SE,
      c("E1", "E2")
    )
  }
})

test_that("heterogeneous residual Bayesian MET works trait by trait for multiple responses", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data()

  for (model in c("GBLUP_BRR", "RKHS")) {
    out <- suppressWarnings(
      fit_bayes_heter_model(model, response = c("Y1", "Y2"), data = dat),
    )
    model_block <- out$model_results_by_model[[bayes_heter_model_key(model)]]
    expect_true(all(c("Y1", "Y2") %in% names(model_block)))

    for (trait in c("Y1", "Y2")) {
      trait_res <- model_block[[trait]]
      expect_met_gaussian_contract(
        trait_res$predicted_values,
        sum(is.na(dat$pheno[[trait]]))
      )
      expect_across_environment_contract(trait_res$across_environment_predicted_values)
      expect_variance_component_contract(trait_res$variance_components)
    }
  }
})

test_that("single-environment Bayesian calls ignore heter controls with a user message", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data()
  single <- dat$pheno[dat$pheno$Env == "E1", c("GID", "Y1"), drop = FALSE]
  single$Y1[c(2L, 7L)] <- NA_real_

  for (model in c("GBLUP_BRR", "RKHS")) {
    expect_message(
      out <- PredictProR::model_execute(
        pheno_data = single,
        gmatrix = dat$gmatrix,
        response = "Y1",
        gen_name = "GID",
        random = ~ GID,
        GS_model = model,
        heter_groups = "Env",
        heter_resid = TRUE,
        response_family = "gaussian",
        para_tunning = FALSE,
        nIter = 60L,
        burnIn = 20L,
        thin = 2L,
        system_database = TRUE,
        message = TRUE
      ),
      "single-environment",
      fixed = TRUE
    )

    trait_res <- out$model_results_by_model[[bayes_heter_model_key(model)]][["Y1"]]
    expect_identical(
      names(trait_res$predicted_values),
      PredictProR:::gp_gaussian_prediction_columns()
    )
  }
})

test_that("joint multi-trait Bayesian BGLR route works for single-environment and MET layouts", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data(n_ids = 7L, n_markers = 6L)

  for (model in c("GBLUP_BRR", "RKHS")) {
    single <- dat$pheno[dat$pheno$Env == "E1", c("GID", "Y1", "Y2"), drop = FALSE]
    out_single <- PredictProR::model_execute(
      pheno_data = single,
      gmatrix = dat$gmatrix,
      response = c("Y1", "Y2"),
      gen_name = "GID",
      GS_model = model,
      multi_trait_bayes = TRUE,
      response_family = "gaussian",
      para_tunning = FALSE,
      nIter = 50L,
      burnIn = 15L,
      thin = 2L,
      system_database = TRUE,
      message = FALSE
    )
    expect_identical(
      names(out_single$model_results$predicted_values),
      PredictProR:::gp_gaussian_prediction_columns(include_trait = TRUE)
    )
    expect_true(all(c("Y1", "Y2") %in% unique(out_single$model_results$predicted_values$Trait)))
    expect_variance_component_contract(out_single$model_results$variance_components)
    expect_true("Trait" %in% names(out_single$model_results$variance_components))
    expect_true(all(is.finite(out_single$model_results$variance_components$Standard_error)))
    expect_covariance_contract(
      out_single$model_results$Genetic_covariance_traits,
      out_single$model_results$Genetic_correlation_traits,
      c("Y1", "Y2")
    )
    expect_covariance_contract(
      out_single$model_results$Residual_covariance_traits,
      out_single$model_results$Residual_correlation_traits,
      c("Y1", "Y2")
    )
    expect_covariance_se_contract(
      out_single$model_results$Genetic_covariance_traits_SE,
      out_single$model_results$Genetic_correlation_traits_SE,
      c("Y1", "Y2")
    )
    expect_covariance_se_contract(
      out_single$model_results$Residual_covariance_traits_SE,
      out_single$model_results$Residual_correlation_traits_SE,
      c("Y1", "Y2")
    )

    out_met <- suppressMessages(PredictProR::model_execute(
      pheno_data = dat$pheno,
      gmatrix = dat$gmatrix,
      response = c("Y1", "Y2"),
      gen_name = "GID",
      GS_model = model,
      multi_trait_bayes = TRUE,
      heter_groups = "Env",
      heter_resid = TRUE,
      response_family = "gaussian",
      para_tunning = FALSE,
      nIter = 50L,
      burnIn = 15L,
      thin = 2L,
      system_database = TRUE,
      message = FALSE
    ))
    expect_identical(
      names(out_met$model_results$predicted_values),
      PredictProR:::gp_gaussian_prediction_columns(include_env = TRUE, include_trait = TRUE)
    )
    expect_s3_class(out_met$model_results$across_environment_predicted_values, "data.frame")
    expect_true(all(c("Y1", "Y2") %in% unique(out_met$model_results$across_environment_predicted_values$Trait)))
    expect_variance_component_contract(out_met$model_results$variance_components)
    expect_true(all(is.finite(out_met$model_results$variance_components$Standard_error)))
    expect_covariance_contract(
      out_met$model_results$Total_genetic_covariance_traits,
      out_met$model_results$Total_genetic_correlation_traits,
      c("Y1", "Y2")
    )
    expect_covariance_contract(
      out_met$model_results$Genetic_covariance_trait_environment,
      out_met$model_results$Genetic_correlation_trait_environment,
      c("Y1", "Y2")
    )
    expect_covariance_se_contract(
      out_met$model_results$Total_genetic_covariance_traits_SE,
      out_met$model_results$Total_genetic_correlation_traits_SE,
      c("Y1", "Y2")
    )
  }
})

test_that("single-trait single-environment marker Bayesian models expose finite uncertainty and variance SE", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data(n_ids = 8L, n_markers = 6L)
  single <- dat$pheno[dat$pheno$Env == "E1", c("GID", "Y1"), drop = FALSE]
  # Marker models take genotypes through QC: homozygous 0/2 dosages of the
  # latent markers (values above 2 need ploidy; >10% heterozygotes are removed).
  marker_dosage <- ifelse(dat$marker > 0, 2, 0)

  for (model in c("BRR", "BayesA", "BayesB", "BayesC", "BL")) {
    set.seed(220 + match(model, c("BRR", "BayesA", "BayesB", "BayesC", "BL")))
    out <- PredictProR::model_execute(
      pheno_data = single,
      geno_data = marker_dosage,
      response = "Y1",
      gen_name = "GID",
      random = ~ GID,
      GS_model = model,
      response_family = "gaussian",
      para_tunning = FALSE,
      nIter = 70L,
      burnIn = 20L,
      thin = 2L,
      system_database = TRUE,
      message = FALSE
    )
    trait_res <- out$model_results_by_model[[bayes_heter_model_key(model)]][["Y1"]]
    pred <- trait_res$predicted_values
    expect_true(all(is.finite(pred$Standard_error)))
    expect_equal(pred$PEV, pred$Standard_error ^ 2, tolerance = 1e-10)
    expect_true(all(pred$Reliability >= 0 & pred$Reliability <= 1))
    expect_variance_component_contract(trait_res$variance_components)
    expect_true(all(is.finite(trait_res$variance_components$Standard_error)))
  }
})

test_that("multi-kernel Bayesian MTMET covariance totals remain PSD", {
  skip_if_not_installed("BGLR")
  dat <- bayes_heter_toy_data(n_ids = 7L, n_markers = 6L)
  out <- suppressMessages(PredictProR::model_execute(
    pheno_data = dat$pheno,
    gmatrix = dat$gmatrix,
    kernel_list = list(omic = dat$omic_kernel),
    response = c("Y1", "Y2"),
    gen_name = "GID",
    GS_model = "RKHS",
    multi_trait_bayes = TRUE,
    heter_groups = "Env",
    response_family = "gaussian",
    para_tunning = FALSE,
    nIter = 55L,
    burnIn = 15L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  ))
  expect_covariance_contract(
    out$model_results$Total_genetic_covariance_traits,
    out$model_results$Total_genetic_correlation_traits,
    c("Y1", "Y2")
  )
  expect_true(all(is.finite(out$model_results$variance_components$Standard_error)))
})

test_that("BGLR wheat data reproduce Bayesian MET and multi-trait uncertainty contracts", {
  skip_if_not_installed("BGLR")
  wheat_env <- new.env(parent = emptyenv())
  utils::data("wheat", package = "BGLR", envir = wheat_env)
  n_lines <- 50L
  ids <- paste0("wheat", seq_len(n_lines))
  wheat_y <- as.matrix(wheat_env$wheat.Y[seq_len(n_lines), seq_len(4L), drop = FALSE])
  wheat_k <- as.matrix(wheat_env$wheat.A[seq_len(n_lines), seq_len(n_lines), drop = FALSE])
  rownames(wheat_k) <- colnames(wheat_k) <- ids

  wheat_met <- expand.grid(
    GID = ids,
    Env = paste0("E", seq_len(4L)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  wheat_met$Yield <- as.numeric(wheat_y)
  wheat_met$Yield[seq(11L, nrow(wheat_met), by = 29L)] <- NA_real_
  met_out <- PredictProR::model_execute(
    pheno_data = wheat_met,
    gmatrix = wheat_k,
    response = "Yield",
    gen_name = "GID",
    random = ~ GID + GID:Env,
    fixed = ~ Env,
    GS_model = "RKHS",
    heter_groups = "Env",
    heter_resid = TRUE,
    response_family = "gaussian",
    para_tunning = FALSE,
    nIter = 55L,
    burnIn = 15L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )
  met_res <- met_out$model_results_by_model[[bayes_heter_model_key("RKHS")]][["Yield"]]
  expect_covariance_contract(
    met_res$Genetic_covariance_environments,
    met_res$Genetic_correlation_environments,
    paste0("E", seq_len(4L))
  )
  expect_true(all(is.finite(met_res$variance_components$Standard_error)))

  wheat_mt <- data.frame(
    GID = ids,
    Trait1 = wheat_y[, 1L],
    Trait2 = wheat_y[, 2L],
    stringsAsFactors = FALSE
  )
  wheat_mt$Trait1[c(7L, 31L)] <- NA_real_
  wheat_mt$Trait2[c(12L, 42L)] <- NA_real_
  mt_out <- PredictProR::model_execute(
    pheno_data = wheat_mt,
    gmatrix = wheat_k,
    response = c("Trait1", "Trait2"),
    gen_name = "GID",
    GS_model = "RKHS",
    multi_trait_bayes = TRUE,
    response_family = "gaussian",
    para_tunning = FALSE,
    nIter = 55L,
    burnIn = 15L,
    thin = 2L,
    system_database = TRUE,
    message = FALSE
  )
  expect_covariance_contract(
    mt_out$model_results$Genetic_covariance_traits,
    mt_out$model_results$Genetic_correlation_traits,
    c("Trait1", "Trait2")
  )
  expect_true(all(is.finite(mt_out$model_results$variance_components$Standard_error)))
})
