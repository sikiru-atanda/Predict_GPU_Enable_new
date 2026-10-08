test_that("joint MT-MET GP heritability includes GxE variance", {
  labels <- c("Trait1", "Trait2")
  covariances <- list(
    trait_levels = labels,
    genetic_covariance = structure(
      matrix(c(2, 0.4, 0.4, 3), 2),
      dimnames = list(labels, labels)
    ),
    gxe_covariance = structure(
      matrix(c(1, 0.1, 0.1, 2), 2),
      dimnames = list(labels, labels)
    ),
    residual_covariance = structure(
      matrix(c(3, 0.2, 0.2, 5), 2),
      dimnames = list(labels, labels)
    )
  )

  out <- PredictProR:::gp_variance_components_from_covariances(covariances)
  h2_rows <- out$Component == "heritability"
  h2 <- out$Components[h2_rows]
  h2 <- h2[match(labels, out$Trait[h2_rows])]

  expect_equal(h2, c((2 + 1) / (2 + 1 + 3), (3 + 2) / (3 + 2 + 5)))
})

test_that("joint GP covariance matrices recover matching axis labels", {
  labels <- c("Trait1", "Trait2")
  cov_mat <- matrix(c(2, 0.4, 0.4, 3), 2,
                    dimnames = list(c("1", "2"), labels))
  source <- list(Genetic_covariance_traits = cov_mat)

  out <- PredictProR:::gp_attach_correlation_matrices(list(), source)

  expect_identical(rownames(out$Genetic_covariance_traits), labels)
  expect_identical(colnames(out$Genetic_covariance_traits), labels)
  expect_equal(unname(out$Genetic_covariance_traits), unname(cov_mat))
})
