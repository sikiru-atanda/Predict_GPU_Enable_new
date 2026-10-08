test_that("multi-trait ASReml omits the fixed residual scale from variance components", {
  raw <- data.frame(
    component = c(0.42, 0.40, 0.39, 1, 0.27, 0.29, 0.34),
    std.error = c(0.20, 0.20, 0.20, NA, 0.13, 0.14, 0.15),
    bound = c("?", "?", "?", "F", "P", "P", "P"),
    row.names = c(
      "Trait:vm(GID, gmatrix)!Trait_Yield:Yield",
      "Trait:vm(GID, gmatrix)!Trait_Height:Yield",
      "Trait:vm(GID, gmatrix)!Trait_Height:Height",
      "units:Trait!R",
      "units:Trait!Trait_Yield:Yield",
      "units:Trait!Trait_Height:Yield",
      "units:Trait!Trait_Height:Height"
    ),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    asreml_varcomp_table = function(...) raw,
    .package = "PredictProR"
  )

  vc <- PredictProR:::gp_multitrait_asreml_variance_components(
    model = structure(list(), class = "mock_asreml"),
    response = c("Yield", "Height")
  )

  expect_equal(nrow(vc), 6L)
  expect_false(any(vc$Components == 1 & vc$Component == "residual_covariance"))
  expect_equal(sum(vc$Component == "genetic_variance"), 2L)
  expect_equal(sum(vc$Component == "genetic_covariance"), 1L)
  expect_equal(sum(vc$Component == "residual_variance"), 2L)
  expect_equal(sum(vc$Component == "residual_covariance"), 1L)
})
