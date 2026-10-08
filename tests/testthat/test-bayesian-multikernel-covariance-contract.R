test_that("joint Bayesian kernel covariance output preserves independently fitted terms", {
  labels <- c("Yield", "Protein")
  main_k1 <- matrix(c(1.2, 0.2, 0.2, 0.8), 2, 2, dimnames = list(labels, labels))
  main_k2 <- matrix(c(0.5, 0.1, 0.1, 0.4), 2, 2, dimnames = list(labels, labels))
  gxe_k1 <- matrix(c(0.3, 0.05, 0.05, 0.2), 2, 2, dimnames = list(labels, labels))
  gxe_k2 <- matrix(c(0.2, 0.03, 0.03, 0.1), 2, 2, dimnames = list(labels, labels))
  draw <- function(a, b, c) rbind(
    c(a * 0.9, b * 0.9, c * 0.9),
    c(a, b, c),
    c(a * 1.1, b * 1.1, c * 1.1)
  )
  bundle <- list(
    main = list(
      covariance_by_term = list(k1_G = main_k1, k2_G = main_k2),
      draws_by_term = list(
        k1_G = draw(1.2, 0.2, 0.8),
        k2_G = draw(0.5, 0.1, 0.4)
      )
    ),
    gxe = list(
      covariance_by_term = list(k1_GxE = gxe_k1, k2_GxE = gxe_k2),
      draws_by_term = list(
        k1_GxE = draw(0.3, 0.05, 0.2),
        k2_GxE = draw(0.2, 0.03, 0.1)
      )
    )
  )
  kernels <- list(k1 = diag(2), k2 = diag(c(2, 1)))
  term_kernel <- c(k1_G = "k1", k1_GxE = "k1", k2_G = "k2", k2_GxE = "k2")

  out <- PredictProR:::bayes_multitrait_kernel_covariance_contract(
    covariance_bundle = bundle,
    term_kernel = term_kernel,
    kernel_matrices = kernels,
    labels = labels,
    is_met = TRUE
  )

  expect_equal(out$Genetic_covariance_by_kernel$k1, main_k1 + gxe_k1)
  expect_equal(out$Genetic_covariance_by_kernel$k2, main_k2 + gxe_k2)
  expect_equal(
    Reduce(`+`, out$Genetic_covariance_by_kernel),
    main_k1 + main_k2 + gxe_k1 + gxe_k2
  )
  expect_setequal(
    unique(out$kernel_variance_components$Component),
    c("kernel_genetic_variance", "kernel_gxe_variance")
  )
  expect_true(all(out$kernel_variance_components$Independent_kernel_estimate))
  expect_true(all(is.finite(out$kernel_variance_components$Standard_error)))
  expect_true(out$kernel_combination$independent_kernel_covariances_estimated)
  expect_identical(
    out$kernel_combination$covariance_basis,
    "sum_of_independently_fitted_BGLR_kernel_covariances"
  )
})

