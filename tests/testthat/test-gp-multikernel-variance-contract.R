test_that("single-trait exact GP exposes independently estimated kernel variances", {
  gp_result <- list(
    var_components_summary = data.frame(
      term_group = "main",
      term = "G",
      kernel = c("genomic", "transcriptomic", "total"),
      component_type = "genetic_main",
      estimate = c(0.8, 0.2, 1.0),
      stringsAsFactors = FALSE
    ),
    varcomp = data.frame(
      component = c(
        "vm(GID,genomic)!var",
        "vm(GID,transcriptomic)!var",
        "Residual!R"
      ),
      estimate = c(0.81, 0, 0.35),
      std.error = c(0.12, NA_real_, 0.08),
      bound = c("P", "B", "P"),
      stringsAsFactors = FALSE
    )
  )

  out <- PredictProR:::gp_public_single_kernel_variance_contract(
    gp_result,
    kernel_names = c("genomic", "transcriptomic"),
    model_name = "GP",
    independently_estimated = TRUE
  )

  expect_identical(out$Kernel, c("genomic", "transcriptomic"))
  expect_equal(out$Components, c(0.81, 0))
  expect_equal(out$Standard_error[[1L]], 0.12)
  expect_true(is.na(out$Standard_error[[2L]]))
  expect_true(all(out$Independent_kernel_estimate))
  expect_identical(out$Boundary[[2L]], "bound_at_zero")
})

test_that("KRR does not report fixed allocations as kernel variance estimates", {
  gp_result <- list(
    var_components_summary = data.frame(
      term_group = "main",
      term = "G",
      kernel = c("G1", "G2", "total"),
      component_type = "genetic_main",
      estimate = c(0.5, 0.5, 1.0),
      stringsAsFactors = FALSE
    ),
    varcomp = data.frame(
      component = c("vm(GID,KRR_kernel)!var", "Residual!R"),
      estimate = c(1, 0.4),
      stringsAsFactors = FALSE
    )
  )

  out <- PredictProR:::gp_public_single_kernel_variance_contract(
    gp_result,
    kernel_names = c("genomic", "metabolomic"),
    model_name = "KRR",
    independently_estimated = FALSE
  )

  expect_identical(out$Kernel, c("genomic", "metabolomic"))
  expect_true(all(is.na(out$Components)))
  expect_equal(out$Fixed_weight_allocation, c(0.5, 0.5))
  expect_true(all(out$Component == "kernel_genetic_variance_not_identifiable"))
  expect_true(all(out$Variance_estimand == "not_identifiable_from_fixed_kernel_weights"))
  expect_false(any(out$Independent_kernel_estimate))
  expect_true(all(out$Estimation_method == "KRR_GCV_combined_kernel_only"))
})

test_that("joint GP covariance output matches the weighted kernel bank used for prediction", {
  traits <- c("Yield", "Protein")
  sigma_g <- matrix(c(0.4, 0.1, 0.1, 0.25), 2L,
                    dimnames = list(traits, traits))
  sigma_e <- matrix(c(0.2, 0.02, 0.02, 0.3), 2L,
                    dimnames = list(traits, traits))
  ids <- paste0("G", 1:3)
  k1 <- diag(1, 3L)
  k2 <- diag(2, 3L)
  dimnames(k1) <- dimnames(k2) <- list(ids, ids)
  kernel <- list(Ks = list(A = k1, omic = k2), K = k1, geno_ids = ids)
  source <- list(
    result = list(),
    fit = list(
      Sigma_G_response_scale = sigma_g,
      Sigma_eps_response_scale = sigma_e,
      kernel_weights = c(0.25, 1.5),
      converged = TRUE
    ),
    info = list(trait_levels = traits)
  )

  out <- PredictProR:::gp_public_joint_multikernel_contract(
    source,
    kernel = kernel,
    kernel_weights = c(0.25, 1.5),
    kernel_names = c("genomic", "omic"),
    varcomp_mode = "reml"
  )

  expect_equal(out$result$genetic_covariance, 1.75 * sigma_g)
  expect_equal(out$result$Genetic_covariance_by_kernel$genomic, 0.25 * sigma_g)
  expect_equal(out$result$Genetic_covariance_by_kernel$omic, 1.5 * sigma_g)
  summed <- aggregate(
    Components ~ Trait,
    data = out$result$kernel_variance_components,
    FUN = sum
  )
  expect_equal(
    unname(summed$Components[match(traits, summed$Trait)]),
    unname(1.75 * diag(sigma_g))
  )
  expect_false(any(out$result$kernel_variance_components$Independent_kernel_estimate))
  expect_equal(out$fit$Combined_kernel_genetic_covariance_coefficient, sigma_g)
})

test_that("nonconverged joint GP does not promote covariance as variance evidence", {
  traits <- c("T1", "T2")
  ids <- paste0("G", 1:3)
  k1 <- diag(3L)
  k2 <- diag(3L)
  dimnames(k1) <- dimnames(k2) <- list(ids, ids)
  sigma_g <- diag(c(0.5, 0.3))
  sigma_e <- diag(c(0.2, 0.25))
  dimnames(sigma_g) <- dimnames(sigma_e) <- list(traits, traits)
  source <- list(
    result = list(),
    fit = list(
      Sigma_G_response_scale = sigma_g,
      Sigma_eps_response_scale = sigma_e,
      kernel_weights = c(1, 1),
      converged = FALSE
    ),
    info = list(trait_levels = traits, varcomp_mode = "reml")
  )
  kernel <- list(Ks = list(A = k1, omic = k2), K = k1, geno_ids = ids)

  out <- PredictProR:::gp_public_joint_multikernel_contract(
    source,
    kernel = kernel,
    kernel_weights = c(1, 1),
    kernel_names = c("genomic", "omic"),
    varcomp_mode = "reml"
  )
  public_source <- list(
    gp_result = out$result,
    gp_fit = out$fit,
    gp_info = out$info
  )
  vc <- PredictProR:::gp_public_variance_components(public_source)

  expect_identical(out$result$variance_component_status, "not_promoted_nonconverged")
  expect_true(all(is.na(out$result$kernel_variance_components$Components)))
  expect_true(all(is.na(vc$Components)))
  expect_true(all(grepl("not_promoted_nonconverged$", vc$Component)))
})

test_that("joint GP accepts and records every supplied kernel weight", {
  ids <- paste0("G", 1:4)
  pheno <- data.frame(
    GID = ids,
    T1 = c(1, 2, 3, NA_real_),
    T2 = c(4, 3, 2, NA_real_)
  )
  k1 <- diag(1, 4L)
  k2 <- diag(1, 4L)
  dimnames(k1) <- dimnames(k2) <- list(ids, ids)
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    gp_bridge_direct_execution_enabled = function() TRUE,
    gp_bridge_fit_public_model_direct = function(args, kernel, ...) {
      captured$args <- args
      captured$kernel_names <- names(kernel$Ks)
      traits <- c("T1", "T2")
      pred <- expand.grid(
        GID = ids,
        Trait = traits,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      pred$Prediction <- 0
      pred$SE <- 0.2
      pred$SE_latent <- 0.2
      pred$PEV <- 0.04
      sigma_g <- diag(c(0.4, 0.3))
      sigma_e <- diag(c(0.2, 0.25))
      dimnames(sigma_g) <- dimnames(sigma_e) <- list(traits, traits)
      list(
        result = list(predictions = pred),
        predictions = pred,
        fit = list(
          Sigma_G_response_scale = sigma_g,
          Sigma_eps_response_scale = sigma_e,
          kernel_weights = args$kernel_weights,
          converged = TRUE
        ),
        info = list(trait_levels = traits)
      )
    },
    .package = "PredictProR"
  )

  out <- PredictProR::gp_multi_trait_model(
    pheno_data = pheno,
    gmatrix = k1,
    kernel_list = list(omic = k2),
    kernel_weights = c(gmatrix = 0.2, omic = 0.8),
    response = c("T1", "T2"),
    gen_name = "GID",
    tune_gp = FALSE,
    prediction_blend = "none",
    gp_uncertainty_calibration = "none",
    return_trait_correlations = TRUE
  )

  expect_identical(captured$kernel_names, c("A", "omic"))
  expect_equal(captured$args$kernel_weights, c(0.2, 0.8))
  expect_identical(out$info$kernel_names, c("gmatrix", "omic"))
  expect_equal(unname(out$info$kernel_weights), c(0.2, 0.8))
  expect_identical(names(out$result$Genetic_covariance_by_kernel), c("gmatrix", "omic"))
})
