# BGLR writes FIXED-term draws to ETA_<k>_b.dat (text). Up to 0.20.223 the
# GBLUP_BRR MET prediction targets looked only for ETA_<k>_b.bin, so the
# environment means of `fixed = ~Env` were left out of every prediction.

test_that("FIXED-term draws are read from BGLR's text file", {
  dir <- withr::local_tempdir()
  dat <- file.path(dir, "runETA_1_b.dat")
  writeLines(c("10 0 0", "11 1 2", "12 2 4", "13 3 6"), dat)
  other <- file.path(dir, "runETA_11_b.dat")
  writeLines("99 99 99", other)
  draws <- PredictProR:::gp_bayes_read_fixed_draws(c(other, dat), eta_index = 1L, n_coef = 3L,
                                                   n_draws = 2L)
  expect_equal(unname(draws), rbind(c(12, 2, 4), c(13, 3, 6)))
  expect_null(PredictProR:::gp_bayes_read_fixed_draws(dat, eta_index = 1L, n_coef = 2L))
  expect_null(PredictProR:::gp_bayes_read_fixed_draws(dat, eta_index = 2L, n_coef = 3L))
})

test_that("MET GBLUP_BRR predictions include the environment fixed effects", {
  skip_if_not_installed("BGLR")
  skip_on_cran()
  withr::local_options(PredictProR.output_dir = withr::local_tempdir())
  set.seed(42)
  n <- 40L
  ids <- sprintf("L%02d", seq_len(n))
  geno <- matrix(2L * sample(0:1, n * 120, replace = TRUE), n, dimnames = list(ids, paste0("m", 1:120)))
  g <- drop(scale(geno) %*% rnorm(120, sd = 0.08))
  pheno <- expand.grid(GID = ids, Env = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
  pheno$Yield <- 20 + c(E1 = 0, E2 = 5, E3 = 10)[pheno$Env] + g[match(pheno$GID, ids)] +
    rnorm(nrow(pheno), sd = 0.5)
  pheno$Yield[pheno$GID %in% ids[1:6]] <- NA_real_
  res <- suppressWarnings(PredictProR::model_execute(
    pheno_data = pheno, geno_data = geno, gen_name = "GID", response = "Yield",
    response_family = "gaussian", heter_groups = "Env", heter_resid = FALSE,
    bayes_kernel_heter_resid = FALSE, met_predict_all_environments = TRUE,
    fixed = ~Env, random = ~GID + GID:Env, GS_model = "GBLUP_BRR",
    gmatrix_method = "VanRaden", bending = TRUE, cross_validation = FALSE,
    nIter = 400L, burnIn = 100L, thin = 5L, random_state = 7L,
    system_database = FALSE, message = FALSE, verbose = FALSE
  ))
  pred <- res$model_results$predicted_values
  tr <- pred$Train_Test_Label == "Train"
  bias <- tapply(pred$Observed_value[tr] - pred$Predicted_value[tr], pred$Env[tr], mean)
  expect_true(all(abs(bias) < 0.5), info = paste(names(bias), round(bias, 3), collapse = " "))
  vc <- res$model_results$variance_components
  expect_false(any(grepl("genetic_variance_genetic_variance", vc$Component)))
  total <- vc$Components[vc$Component == "total_genetic_variance"]
  parts <- vc$Components[grepl("^genetic_variance_", vc$Component) &
                           !grepl("_E[0-9]+$", vc$Component)]
  # the total covers every random genetic term (main + GxE)
  expect_gt(total, max(parts))
})
