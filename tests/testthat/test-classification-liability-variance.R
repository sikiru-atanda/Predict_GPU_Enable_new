test_that("RKHS classification variance is reported on the probit liability scale", {
  td <- tempfile("liability-rkhs-")
  dir.create(td)
  variance_file <- file.path(td, "fit_ETA_1_varU.dat")
  writeLines(c("0.5", "1", "1.5", "2"), variance_file)

  pheno <- data.frame(
    GID = paste0("G", 1:4),
    Trait = ordered(c("low", "high", "low", NA), levels = c("low", "high"))
  )
  eta <- list(ETA = list(list(K = diag(2, 4), model = "RKHS")))
  out <- PredictProR:::gp_bayes_classification_liability_variance(
    mod = list(output_files_names = variance_file),
    ETA = eta,
    pheno_data = pheno,
    gen_name = "GID",
    bayes_para = list(nIter = 8L, burnIn = 4L, thin = 2L),
    confidence_level = 0.95
  )

  vc <- out$variance_components
  expect_equal(
    vc$Component,
    c(
      "liability_scale_genetic_variance",
      "liability_scale_residual_variance_fixed",
      "liability_scale_heritability"
    )
  )
  expect_equal(vc$Components[[1L]], mean(c(3, 4)), tolerance = 1e-12)
  expect_equal(vc$Components[[2L]], 1)
  expect_equal(vc$Standard_error[[2L]], 0)
  expect_equal(vc$Components[[3L]], mean(c(3 / 4, 4 / 5)), tolerance = 1e-12)

  intervals <- out$variance_component_intervals
  expect_true(all(intervals$Scale == "latent_liability"))
  expect_true(all(intervals$Link_function == "probit"))
  expect_false(intervals$Residual_variance_estimated[[2L]])
  expect_true(all(is.finite(intervals$Lower_credible_limit)))
  expect_true(all(is.finite(intervals$Upper_credible_limit)))
})

test_that("GBLUP_BRR classification uses its implied kernel covariance", {
  td <- tempfile("liability-brr-")
  dir.create(td)
  variance_file <- file.path(td, "fit_ETA_1_varB.dat")
  writeLines(c("0.25", "0.5", "0.75", "1"), variance_file)

  pheno <- data.frame(GID = paste0("G", 1:3), Trait = c("no", "yes", NA))
  design <- diag(sqrt(3), 3)
  eta <- list(ETA = list(list(X = design, model = "BRR")))
  out <- PredictProR:::gp_bayes_classification_liability_variance(
    mod = list(output_files_names = variance_file),
    ETA = eta,
    pheno_data = pheno,
    gen_name = "GID",
    bayes_para = list(nIter = 8L, burnIn = 4L, thin = 2L)
  )

  genetic <- out$variance_components$Components[
    out$variance_components$Component == "liability_scale_genetic_variance"
  ]
  expect_equal(genetic, mean(c(2.25, 3)), tolerance = 1e-12)
})

test_that("MET liability variance separates genetic main and GxE components", {
  td <- tempfile("liability-met-")
  dir.create(td)
  main_file <- file.path(td, "fit_ETA_1_varB.dat")
  gxe_file <- file.path(td, "fit_ETA_2_varU.dat")
  writeLines(c("0.4", "0.6", "0.8", "1.0"), main_file)
  writeLines(c("0.2", "0.3", "0.4", "0.5"), gxe_file)

  gid <- rep(paste0("G", 1:3), 2)
  env <- rep(c("E1", "E2"), each = 3)
  incidence <- stats::model.matrix(~ factor(gid) - 1)
  base_kernel <- tcrossprod(incidence)
  gxe_kernel <- base_kernel * tcrossprod(stats::model.matrix(~ factor(env) - 1))
  eta <- list(ETA = list(
    list(X = PredictProR:::gp_bayes_brr_design_from_kernel(base_kernel), model = "BRR"),
    list(K = gxe_kernel, model = "RKHS")
  ))
  pheno <- data.frame(
    GID = gid,
    Env = env,
    Trait = ordered(c("low", "medium", "high", "low", "medium", NA),
                    levels = c("low", "medium", "high"))
  )
  out <- PredictProR:::gp_bayes_classification_liability_variance(
    mod = list(output_files_names = c(main_file, gxe_file)),
    ETA = eta,
    pheno_data = pheno,
    gen_name = "GID",
    heter_groups = "Env",
    bayes_para = list(nIter = 8L, burnIn = 4L, thin = 2L)
  )

  vc <- out$variance_components
  expect_true("Env" %in% names(vc))
  expect_setequal(unique(stats::na.omit(vc$Env)), c("E1", "E2"))
  expect_true(all(c(
    "liability_scale_genetic_main_variance",
    "liability_scale_gxe_variance",
    "liability_scale_total_genetic_variance",
    "liability_scale_residual_variance_fixed",
    "liability_scale_heritability"
  ) %in% vc$Component[!is.na(vc$Env)]))
  residual <- vc[vc$Component == "liability_scale_residual_variance_fixed" &
                   !is.na(vc$Env), , drop = FALSE]
  expect_equal(residual$Components, c(1, 1))
  expect_equal(residual$Standard_error, c(0, 0))
})
