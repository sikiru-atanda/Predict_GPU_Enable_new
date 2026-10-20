# A default multi-trait GP run (gp_return_trait_correlations = FALSE) reported
# genetic and residual covariance matrices without correlation tables up to
# 0.20.223; every other engine reports both.
test_that("default multi-trait GP reports correlation tables with its covariances", {
  skip_on_cran()
  skip_if_no_gp_python()
  withr::local_options(PredictProR.output_dir = withr::local_tempdir())
  set.seed(8)
  n <- 40L
  ids <- sprintf("L%02d", seq_len(n))
  geno <- matrix(2L * sample(0:1, n * 150, replace = TRUE), n, dimnames = list(ids, paste0("m", 1:150)))
  g <- drop(scale(geno) %*% rnorm(150, sd = 0.08))
  pheno <- data.frame(GID = ids, Yield = 10 + g + rnorm(n, sd = 0.6), stringsAsFactors = FALSE)
  pheno$Starch <- 0.5 * pheno$Yield + rnorm(n)
  pheno$Yield[1:6] <- NA
  res <- suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = geno, gen_name = "GID", response = c("Yield", "Starch"),
    response_family = "gaussian", multi_trait_gp = TRUE, GS_model = "Gaussian-Process-GBLUP",
    gmatrix_method = "VanRaden", cross_validation = FALSE, random_state = 5L,
    system_database = FALSE, message = FALSE, verbose = FALSE
  ))
  mr <- res$model_results
  cov_m <- unique(mr$covariance_components$Matrix)
  cor_m <- unique(mr$correlation_components$Matrix)
  expect_true(all(sub("covariance", "correlation", cov_m) %in% cor_m))
  audit <- PredictProR:::gp_audit_output_tables(res)
  expect_identical(audit$check[audit$ok %in% FALSE], character())
})
