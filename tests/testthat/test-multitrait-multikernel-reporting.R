# Multi-trait ASReml multi-kernel reporting contract.
#
# A multi-trait run may fit several independent genetic terms (gmatrix plus any
# extra kernels). Every multi-trait ASReml route must say how many kernels it
# actually fitted and name them, so a two-kernel run cannot be mistaken for a
# one-kernel run. The contract is the same for single-environment and MET, and
# for cross-validation and true prediction:
#   * kernel_configuration : one row per fitted kernel, with its role
#   * model_parameters     : multi_kernel_count / _names / _strategy rows

mtmk_inputs <- function(n = 6L) {
  ids <- paste0("G", seq_len(n))
  G <- diag(n)
  COP <- diag(n) + 0.1
  dimnames(G) <- dimnames(COP) <- list(ids, ids)
  list(ids = ids, G = G, COP = COP)
}

mtmk_param <- function(params, stat) {
  params$summary[params$stat == stat]
}

mtmk_met_fixture <- function() {
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

test_that("kernel configuration names every fitted kernel and its role", {
  cfg <- PredictProR:::gp_multitrait_kernel_configuration(
    list(gmatrix = diag(2), COP = diag(2))
  )
  expect_identical(cfg$Kernel, c("gmatrix", "COP"))
  expect_identical(
    unique(cfg$Kernel_role), "independent_genetic_covariance_term"
  )
})

test_that("single-environment multi-trait ASReml CV reports both kernels", {
  skip_if_not_installed("asreml")  # the CV entry checks for asreml before the mocked fit
  k <- mtmk_inputs()
  ph <- data.frame(
    GID = k$ids,
    Yield = seq_along(k$ids) / 3,
    Height = seq_along(k$ids) / 5 + 2,
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    gp_multitrait_asreml_gaussian_model = function(pheno_object, response,
                                                   gen_name, ...) {
      pred <- expand.grid(
        id_value = as.character(pheno_object[[gen_name]]),
        Trait = response,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      names(pred)[1] <- gen_name
      pred$Predicted_value <- seq_len(nrow(pred)) / 10
      list(predicted_values = pred)
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_multitrait_asreml_gaussian_cv(
    pheno_object = ph,
    response = c("Yield", "Height"),
    gen_name = "GID",
    gmatrix = k$G,
    kernel_list = list(COP = k$COP),
    cross_validation_meth = "K-Folds",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )

  processed <- out$cv_results_processed
  expect_setequal(processed$kernel_configuration$Kernel, c("gmatrix", "COP"))
  expect_identical(
    mtmk_param(processed$model_parameters, "multi_kernel_count"), "2"
  )
  expect_match(
    mtmk_param(processed$model_parameters, "multi_kernel_names"), "COP"
  )
})

test_that("grouped MT-MET ASReml CV reports both kernels in model parameters", {
  k <- mtmk_inputs()
  ph <- expand.grid(
    GID = k$ids,
    Env = paste0("E", seq_len(3)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ph$Yield <- seq_len(nrow(ph)) / 3
  ph$Height <- seq_len(nrow(ph)) / 5 + 2

  testthat::local_mocked_bindings(
    gp_multitrait_asreml_mtmet_gaussian_model = function(pheno_object, response,
                                                         gen_name, heter_groups,
                                                         ...) {
      bundle <- PredictProR:::gp_multitrait_asreml_mtmet_long_data(
        pheno_object,
        response = response,
        gen_name = gen_name,
        heter_groups = heter_groups
      )
      pred <- merge(
        data.frame(
          GID = unique(as.character(pheno_object[[gen_name]])),
          stringsAsFactors = FALSE
        ),
        bundle$group_lookup,
        by = NULL
      )
      names(pred)[names(pred) == "GID"] <- gen_name
      pred$Predicted_value <- seq_len(nrow(pred)) / 10
      list(predicted_values = pred)
    },
    .package = "PredictProR"
  )

  out <- PredictProR:::gp_multitrait_asreml_mtmet_gaussian_cv(
    pheno_object = ph,
    response = c("Yield", "Height"),
    gen_name = "GID",
    heter_groups = "Env",
    gmatrix = k$G,
    kernel_list = list(COP = k$COP),
    cross_validation_meth = "CV1",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )

  processed <- out$cv_results_processed
  expect_identical(
    mtmk_param(processed$model_parameters, "multi_kernel_count"), "2"
  )
  expect_match(
    mtmk_param(processed$model_parameters, "multi_kernel_names"), "COP"
  )
})

test_that("single-environment multi-trait ASReml prediction reports both kernels", {
  skip_if_not_installed("asreml")
  k <- mtmk_inputs()
  set.seed(11)
  ph <- data.frame(
    GID = k$ids,
    Yield = rnorm(length(k$ids)),
    Height = rnorm(length(k$ids)),
    stringsAsFactors = FALSE
  )

  fit <- tryCatch(
    PredictProR:::gp_multitrait_asreml_gaussian_model(
      pheno_object = ph,
      response = c("Yield", "Height"),
      gen_name = "GID",
      gmatrix = k$G,
      kernel_list = list(COP = k$COP),
      trait_covariance = "diag",
      extract_variance = FALSE
    ),
    error = function(e) skip(paste("ASReml fit unavailable:", conditionMessage(e)))
  )

  expect_setequal(fit$kernel_configuration$Kernel, c("gmatrix", "COP"))
  expect_identical(
    mtmk_param(fit$model_parameters, "multi_kernel_count"), "2"
  )
})

test_that("grouped MT-MET ASReml prediction reports both kernels", {
  skip_if_not_installed("asreml")
  k <- mtmk_inputs()
  ph <- expand.grid(
    GID = k$ids,
    Env = paste0("E", seq_len(2)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  set.seed(12)
  ph$Yield <- rnorm(nrow(ph))
  ph$Height <- rnorm(nrow(ph))

  fit <- tryCatch(
    PredictProR:::gp_multitrait_asreml_mtmet_gaussian_model(
      pheno_object = ph,
      response = c("Yield", "Height"),
      gen_name = "GID",
      heter_groups = "Env",
      gmatrix = k$G,
      kernel_list = list(COP = k$COP),
      heter_resid = TRUE,
      extract_variance = FALSE
    ),
    error = function(e) skip(paste("ASReml fit unavailable:", conditionMessage(e)))
  )

  expect_setequal(fit$kernel_configuration$Kernel, c("gmatrix", "COP"))
  expect_identical(
    mtmk_param(fit$model_parameters, "multi_kernel_count"), "2"
  )
})

test_that("multi-trait GP and Bayesian CV report every kernel they fitted", {
  ph <- mtmk_met_fixture()
  traits <- c("Trait1", "Trait2")
  K <- diag(6)
  rownames(K) <- colnames(K) <- paste0("G", seq_len(6))
  COP <- K + 0.1

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
    kernel_list = list(COP = COP),
    cross_validation_meth = "CV1",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )
  expect_setequal(
    gp$cv_results_processed$kernel_configuration$Kernel, c("gmatrix", "COP")
  )
  expect_identical(
    mtmk_param(gp$cv_results_processed$model_parameters, "multi_kernel_count"),
    "2"
  )

  bayes <- PredictProR:::gp_multitrait_bayes_gaussian_cv(
    model_type = "RKHS",
    pheno_object = ph,
    response = traits,
    gen_name = "Sample",
    heter_groups = "Env",
    kernels = list(genomic = K, COP = COP),
    bayes_para = list(nIter = 100L, burnIn = 20L, thin = 2L),
    cross_validation_meth = "CV2",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )
  expect_setequal(
    bayes$cv_results_processed$kernel_configuration$Kernel, c("genomic", "COP")
  )
  expect_identical(
    mtmk_param(bayes$cv_results_processed$model_parameters, "multi_kernel_count"),
    "2"
  )

  # The two engines combine kernels differently and must not claim otherwise.
  # Joint GP sums the bank into one fixed-weight kernel (mt_gp.py builds
  # K_geno += w * K), so per-kernel variances are not separately identified.
  # BGLR gives each kernel its own ETA term, so they are.
  expect_identical(
    mtmk_param(gp$cv_results_processed$model_parameters, "multi_kernel_strategy"),
    "fixed_weight_combined_kernel"
  )
  expect_identical(
    unique(gp$cv_results_processed$kernel_configuration$Kernel_role),
    "fixed_weight_component_of_combined_kernel"
  )
  expect_identical(
    mtmk_param(bayes$cv_results_processed$model_parameters, "multi_kernel_strategy"),
    "independent_bayesian_eta_term_per_kernel"
  )
  expect_identical(
    unique(bayes$cv_results_processed$kernel_configuration$Kernel_role),
    "independent_genetic_covariance_term"
  )
})

test_that("single-environment CV kernel report collapses under feature selection", {
  multi <- PredictProR:::gp_multitrait_asreml_cv_kernel_report(
    kernel_names = c("gmatrix", "COP"),
    trait_covariance = "us",
    feature_active = FALSE
  )
  expect_identical(multi$configuration$Kernel, c("gmatrix", "COP"))
  expect_identical(mtmk_param(multi$parameters, "multi_kernel_count"), "2")

  selected <- PredictProR:::gp_multitrait_asreml_cv_kernel_report(
    kernel_names = c("gmatrix", "COP"),
    trait_covariance = "us",
    feature_active = TRUE
  )
  expect_identical(mtmk_param(selected$parameters, "multi_kernel_count"), "1")
  expect_match(selected$configuration$Kernel_role, "selected_features")
  expect_match(
    mtmk_param(selected$parameters, "multi_kernel_strategy"),
    "replaced by one fold-rebuilt genomic kernel"
  )
})

test_that("GBLUP_BRR and RKHS kernel contract reports its fitted kernels", {
  labels <- c("Yield", "Protein")
  mk <- function(a, b, c) matrix(c(a, b, b, c), 2, 2, dimnames = list(labels, labels))
  draw <- function(a, b, c) rbind(c(a * 0.9, b * 0.9, c * 0.9),
                                  c(a, b, c),
                                  c(a * 1.1, b * 1.1, c * 1.1))
  bundle <- list(
    main = list(
      covariance_by_term = list(k1_G = mk(1.2, 0.2, 0.8), k2_G = mk(0.5, 0.1, 0.4)),
      draws_by_term = list(k1_G = draw(1.2, 0.2, 0.8), k2_G = draw(0.5, 0.1, 0.4))
    )
  )
  out <- PredictProR:::bayes_multitrait_kernel_covariance_contract(
    covariance_bundle = bundle,
    term_kernel = c(k1_G = "k1", k2_G = "k2"),
    kernel_matrices = list(k1 = diag(2), k2 = diag(c(2, 1))),
    labels = labels,
    is_met = FALSE
  )

  expect_setequal(out$kernel_configuration$Kernel, c("k1", "k2"))
  expect_identical(
    unique(out$kernel_configuration$Kernel_role),
    "independent_genetic_covariance_term"
  )
})

test_that("joint GP CV labels the mixture as fitted only when it was fitted", {
  ph <- mtmk_met_fixture()
  traits <- c("Trait1", "Trait2")
  K <- diag(6)
  rownames(K) <- colnames(K) <- paste0("G", seq_len(6))
  COP <- K + 0.1

  single_env <- ph[ph$Env == "E1", , drop = FALSE]

  prediction_grid <- function(masked, with_env) {
    do.call(rbind, lapply(traits, function(trait) {
      out <- data.frame(
        GID = masked$Sample,
        Trait = trait,
        Predicted_value = seq_len(nrow(masked)) / 10,
        Standard_error = 0.2,
        PEV = 0.04,
        stringsAsFactors = FALSE
      )
      if (with_env) out$Env <- masked$Env
      out
    }))
  }

  testthat::local_mocked_bindings(
    gp_multi_trait_model = function(pheno_data, ...) {
      list(predictions = prediction_grid(pheno_data, with_env = FALSE))
    },
    gp_multi_trait_met_model = function(pheno_data, ...) {
      list(predictions = prediction_grid(pheno_data, with_env = TRUE))
    },
    .package = "PredictProR"
  )

  fitted_mix <- PredictProR:::gp_multitrait_gp_gaussian_cv(
    model_type = "GP",
    pheno_object = single_env,
    response = traits,
    gen_name = "Sample",
    heter_groups = NULL,
    gmatrix = K,
    kernel_list = list(COP = COP),
    estimate_kernel_weights = TRUE,
    cross_validation_meth = "K-Folds",
    nfolds = 2L,
    replication = 1L,
    eval_metrics = "mean_squared_error"
  )
  expect_identical(
    mtmk_param(fitted_mix$cv_results_processed$model_parameters,
               "multi_kernel_strategy"),
    "reml_estimated_weight_combined_kernel"
  )

  # MET has no likelihood to profile, so the request must be refused out loud
  # rather than silently reported as a fitted mixture.
  expect_warning(
    met <- PredictProR:::gp_multitrait_gp_gaussian_cv(
      model_type = "GP",
      pheno_object = ph,
      response = traits,
      gen_name = "Sample",
      heter_groups = "Env",
      gmatrix = K,
      kernel_list = list(COP = COP),
      estimate_kernel_weights = TRUE,
      cross_validation_meth = "CV1",
      nfolds = 2L,
      replication = 1L,
      eval_metrics = "mean_squared_error"
    ),
    "not available for joint multi-trait MET GP"
  )
  expect_identical(
    mtmk_param(met$cv_results_processed$model_parameters, "multi_kernel_strategy"),
    "fixed_weight_combined_kernel"
  )
})
