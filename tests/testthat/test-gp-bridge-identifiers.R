# Genotype identifiers are opaque text. Up to 0.20.223 the GP bridge parsed
# phenotype and prediction CSVs with default type guessing, so zero-padded IDs
# such as "0001" became 1 and every GP route stopped with "Unknown genotype IDs
# encountered". 0.20.224 reads identifier columns as text in every bridge
# reader.

test_that("GP prediction files keep opaque identifiers", {
  path <- tempfile(fileext = ".csv")
  writeLines(c('Name,Prediction,Prediction_SE', '001,1.2,0.3', 'NA,2.3,', '" A B ",3.4,0.2'), path)
  p <- PredictProR:::gp_bridge_read_predictions(path)
  expect_identical(p$Name, c("001", "NA", " A B "))
  expect_equal(p$Prediction, c(1.2, 2.3, 3.4))
  expect_true(is.na(p$Prediction_SE[2]))
})

test_that("mode imputation counts exact values", {
  x <- matrix(c(0.1 + 0.2, 0.3, 0.3, 1, 1, NA), ncol = 1)
  colnames(x) <- "m1"
  # 0.1 + 0.2 and 0.3 differ in the last digit; table() labels merged them
  expect_equal(unname(PredictProR:::gp_fit_summary_imputation(x[1:5, , drop = FALSE], "mode")), 0.3)
  expect_equal(unname(PredictProR:::gp_fit_summary_imputation(cbind(m = c(2, 0, 2, 0)), "mode")), 0)
})

test_that("GP routes run with zero-padded numeric genotype IDs", {
  skip_on_cran()
  skip_if_no_gp_python()
  withr::local_options(PredictProR.output_dir = withr::local_tempdir())
  set.seed(4)
  n <- 40L
  ids <- sprintf("%04d", seq_len(n))
  geno <- matrix(2L * sample(0:1, n * 150, replace = TRUE), n, dimnames = list(ids, paste0("m", 1:150)))
  g <- drop(scale(geno) %*% rnorm(150, sd = 0.08))
  single <- data.frame(GID = ids, Yield = 10 + g + rnorm(n, sd = 0.6), stringsAsFactors = FALSE)
  single$Yield[1:6] <- NA
  met <- expand.grid(GID = ids, Env = c("E1", "E2"), stringsAsFactors = FALSE)
  met$Yield <- 10 + c(E1 = 0, E2 = 2)[met$Env] + g[match(met$GID, ids)] + rnorm(nrow(met), sd = 0.6)
  met$Yield[met$GID %in% ids[1:6]] <- NA
  multi <- single
  multi$Starch <- 0.5 * multi$Yield + rnorm(n)
  runs <- list(
    list(pheno = single, model = "Gaussian-Process-GBLUP", args = list()),
    list(pheno = single, model = "Kernel-GBLUP", args = list()),
    list(pheno = met, model = "Gaussian-Process-GBLUP",
         args = list(heter_groups = "Env", fixed = ~Env, random = ~GID + GID:Env)),
    list(pheno = multi, model = "Gaussian-Process-GBLUP", args = list(multi_trait_gp = TRUE),
         response = c("Yield", "Starch"))
  )
  for (r in runs) {
    res <- suppressWarnings(do.call(PredictProR::model_execute, c(list(
      pheno_data = r$pheno, geno_data = geno, gen_name = "GID", response = r$response %||% "Yield",
      response_family = "gaussian", GS_model = r$model, gmatrix_method = "VanRaden",
      cross_validation = FALSE, random_state = 3L, system_database = FALSE, message = FALSE,
      verbose = FALSE), r$args)))
    pv <- res$model_results$predicted_values
    expect_true(is.data.frame(pv) && nrow(pv) > 0, info = r$model)
    expect_true(all(pv$GID %in% ids), info = r$model)
    audit <- PredictProR:::gp_audit_output_tables(res)
    expect_identical(audit$check[audit$ok %in% FALSE], character(), info = r$model)
  }
})
