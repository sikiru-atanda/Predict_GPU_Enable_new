test_that("relationship matrix diagnostics enforce covariance invariants", {
  ids <- paste0("g", 1:3)
  K <- matrix(
    c(1, 0.2, 0.1, 0.2, 1, 0.3, 0.1, 0.3, 1),
    nrow = 3,
    dimnames = list(ids, ids)
  )

  diagnostics <- relationship_matrix_diagnostics(K)
  expect_true(diagnostics$valid)
  expect_gte(diagnostics$minimum_eigenvalue, 0)
  expect_identical(diagnostics$effective_rank, 3L)
  expect_no_error(validate_relationship_matrix(K))

  asymmetric <- K
  asymmetric[1, 2] <- 0.8
  expect_false(relationship_matrix_diagnostics(asymmetric)$valid)
  expect_error(validate_relationship_matrix(asymmetric), "symmetric")

  indefinite <- K
  indefinite[1, 1] <- -1
  expect_false(relationship_matrix_diagnostics(indefinite)$valid)
  expect_error(validate_relationship_matrix(indefinite), "positive semidefinite")
})

test_that("factor analytic covariance is Lambda Lambda-prime plus Psi", {
  loadings <- matrix(
    c(0.8, 0.1, 0.4, -0.2, 0.6, 0.3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("E1", "E2", "E3"), c("F1", "F2"))
  )
  specific <- c(0.2, 0.1, 0.3)
  expected <- loadings %*% t(loadings) + diag(specific)
  dimnames(expected) <- list(rownames(loadings), rownames(loadings))

  actual <- PredictProR:::gp_factor_analytic_covariance(loadings, specific)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_true(relationship_matrix_diagnostics(actual)$valid)
  expect_error(
    PredictProR:::gp_factor_analytic_covariance(loadings, c(0.2, -0.1, 0.3)),
    "cannot be negative"
  )
})

test_that("correlation covariance constructors preserve fitted variances", {
  labels <- c("E1", "E2", "E3")
  variances <- c(1, 2, 3)
  correlation <- PredictProR:::gp_correlation_matrix_from_parameters(
    0.25,
    n = 3,
    compound_symmetry = TRUE,
    labels = labels
  )
  covariance <- PredictProR:::gp_covariance_from_correlation(
    variances,
    correlation,
    labels = labels
  )

  expect_identical(unname(diag(correlation)), rep(1, 3))
  expect_equal(unname(diag(covariance)), variances, tolerance = 1e-12)
  expect_equal(covariance[1, 2], 0.25 * sqrt(2), tolerance = 1e-12)
  expect_error(
    PredictProR:::gp_correlation_matrix_from_parameters(c(0.1, 0.2), n = 3),
    "Expected 3 correlation parameters"
  )
})

test_that("ASReml corh reconstruction does not scale diagonal variances by rho", {
  envs <- c("E1", "E2", "E3")
  kernel_rows <- c(
    "corh(YYY):vm(GID,gmatrix_inv)!cor",
    paste0("corh(YYY):vm(GID,gmatrix_inv)!", envs, "!var")
  )
  residual_rows <- paste0("YYY_", envs, "!R")
  components <- c(0.25, 1, 2, 3, 0.2, 0.3, 0.4)
  varcomp <- data.frame(
    component = components,
    std.error = rep(0.05, length(components)),
    bound = "P",
    row.names = c(kernel_rows, residual_rows),
    stringsAsFactors = FALSE
  )
  model <- structure(
    list(
      varcomp = varcomp,
      mf = data.frame(YYY = factor(rep(envs, each = 2), levels = envs))
    ),
    class = "summary.asreml"
  )

  result <- asreml_herit_varCov_new(
    model = model,
    heter_groups = "YYY",
    var_cov_str = "corh",
    heter_resid = TRUE,
    names_in_inv_list = "gmatrix_inv",
    inter_gen_pos = 2L,
    gen_pos = 1L
  )

  expect_equal(unname(diag(result$Covariance[[1L]])), c(1, 2, 3), tolerance = 1e-12)
  expect_identical(unname(diag(result$Correlation[[1L]])), rep(1, 3))
  expect_equal(as.numeric(result$Total_genetic_var), c(1, 2, 3), tolerance = 1e-12)
})

test_that("kernel hyperparameters cannot create invalid covariance mixtures", {
  x <- matrix(
    c(0, 1, 1, 0, 2, 1, 1, 2, 0),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("m", 1:3))
  )

  expect_error(
    kernel_calculation(x, scaling = FALSE, method = "Gaussian_kernel", theta = 0, message = FALSE),
    "theta must be"
  )
  expect_error(
    kernel_calculation(x, scaling = FALSE, method = "Composite_kernel", alpha = 1.1, message = FALSE),
    "alpha must be"
  )

  weights <- matrix(c(1, -0.2, 1), ncol = 1,
                    dimnames = list(colnames(x), "weight"))
  expect_error(
    PredictProR:::grm_calculation(x, weight = weights, method = "Weighted_VanRaden"),
    "non-negative"
  )
})

test_that("factor analytic terms are built independently for all kernels", {
  pheno <- data.frame(
    GID = factor(rep(c("g1", "g2"), each = 3)),
    Env = factor(rep(c("E1", "E2", "E3"), times = 2)),
    Yield = seq_len(6)
  )
  result <- PredictProR:::random_terms_fit_new(
    random = ~ GID + GID:Env,
    fixed = ~ Env,
    fixed_term = "Env",
    heter_groups = "Env",
    heter_resid = TRUE,
    var_cov_str = "fa1",
    code_asr = c("asreml::asreml(fixed=trait~1+Env)", "random=~", "residual=~"),
    names_in_inv_list = c("gmatrix_inv", "omics_inv"),
    gen_name = "GID",
    pheno_data = pheno
  )

  expect_match(result$code_asr[[2]], "fa\\(Env,1\\):vm\\(GID,gmatrix_inv\\)")
  expect_match(result$code_asr[[2]], "fa\\(Env,1\\):vm\\(GID,omics_inv\\)")
})
