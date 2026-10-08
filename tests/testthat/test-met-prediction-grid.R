# Multi-environment true prediction returns the full line x environment grid
# with every engine, labelled Train / Test / Unobserved (it used to depend on
# the engine: Bayes, GP and ML/DL returned only the input records, ASReml the
# grid with every extra row labelled "Train").

met_pheno <- function() {
  data.frame(
    GID = c("A", "A", "A", "B", "B", "C", "C"),
    Env = c("E1", "E2", "E3", "E2", "E3", "E1", "E2"),
    Loc = c("north", "south", "east", "south", "east", "north", "south"),   # constant within Env
    Group = c("g1", "g1", "g1", "g2", "g2", "g1", "g1"),                    # constant within GID
    Weight = c(2, 2, 2, 2, 2, 2, 2),
    Yield = c(1, 2, 3, 4, 5, NA, NA),
    stringsAsFactors = FALSE
  )
}

test_that("the grid is completed with NA records and filled covariates", {
  out <- PredictProR:::gp_met_complete_grid(met_pheno(), "GID", "Env", "Yield", fixed = ~Env, weights = "Weight")
  expect_identical(nrow(out), 9L)
  added <- out[8:9, ]
  expect_setequal(paste(added$GID, added$Env), c("B E1", "C E3"))
  expect_true(all(is.na(added$Yield)))
  expect_identical(added$Loc[added$Env == "E1"], "north")
  expect_identical(added$Group[added$GID == "B"], "g2")
  expect_identical(added$Weight, c(1, 1))
  # a fixed covariate that varies within environments is not guessed
  ph <- met_pheno(); ph$Rep <- c(1, 2, 1, 2, 1, 2, 1)
  expect_message(same <- PredictProR:::gp_met_complete_grid(ph, "GID", "Env", "Yield", fixed = ~Env + Rep),
                 "observed line x Env combinations only")
  expect_identical(nrow(same), 7L)
  # nothing missing: unchanged
  full <- PredictProR:::gp_met_complete_grid(out, "GID", "Env", "Yield")
  expect_identical(nrow(full), 9L)
})

test_that("rows are labelled Train, Test or Unobserved", {
  ph <- met_pheno()
  pred <- data.frame(GID = rep(c("A", "B", "C"), each = 3), Env = rep(c("E1", "E2", "E3"), 3),
                     Predicted_value = 0, Train_Test_Label = "Train", stringsAsFactors = FALSE)
  keys <- paste(ph$GID, ph$Env, sep = "\r")            # input records before the grid was completed
  out <- PredictProR:::gp_met_relabel_predictions(pred, ph, "GID", "Env", "Yield", input_keys = keys)
  expect_identical(out$Train_Test_Label,
                   c("Train", "Train", "Train", "Unobserved", "Train", "Train", "Test", "Test", "Test"))
  # without added rows every NA record is a user test record
  out0 <- PredictProR:::gp_met_relabel_predictions(pred, ph, "GID", "Env", "Yield")
  expect_identical(out0$Train_Test_Label[4], "Test")
  # multi-trait tables are labelled per trait; an NA record the user supplied
  # (C in E2) is a test record, a combination added to the grid is Unobserved
  ph$Protein <- c(NA, NA, NA, 1, 1, 1, NA)
  mt <- rbind(cbind(pred, Trait = "Yield"), cbind(pred, Trait = "Protein"))
  out_mt <- PredictProR:::gp_met_relabel_predictions(mt, ph, "GID", "Env", c("Yield", "Protein"), input_keys = keys)
  prot <- out_mt[out_mt$Trait == "Protein", "Train_Test_Label"]
  expect_identical(prot, c("Test", "Test", "Test", "Unobserved", "Train", "Train", "Train", "Test", "Unobserved"))
})

# Lines present only in the genotype data are always predicted. They used to
# be dropped whenever pheno_data had any NA, so whether such a line was
# predicted depended on how the other lines were entered.
test_that("genotyped-only lines are predicted alongside NA lines", {
  skip_on_cran()
  set.seed(8)
  ids <- sprintf("L%02d", 1:30)
  G <- matrix(2L * stats::rbinom(30 * 80, 1, 0.5), 30, dimnames = list(ids, sprintf("m%02d", 1:80)))
  ph <- data.frame(GID = ids[1:26], Yield = stats::rnorm(26), stringsAsFactors = FALSE)
  ph$Yield[1:4] <- NA                                        # NA lines ...
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, response = "Yield", gen_name = "GID", random = ~GID,
    GS_model = "GBLUP_BRR", gmatrix_method = "VanRaden", qc_filtering = FALSE, impute = FALSE,
    ld_prunning_qc = FALSE, nIter = 200L, burnIn = 50L, thin = 2L, parallel_mode = "sequential",
    system_database = TRUE, message = FALSE
  )))
  pv <- res$model_results$predicted_values
  expect_setequal(pv$GID, ids)                               # ... and L27-L30 (genotypes only)
  expect_true(all(pv$Train_Test_Label[pv$GID %in% c(ids[1:4], ids[27:30])] == "Test"))
  expect_true(all(is.finite(pv$Predicted_value)))
})

# geno_data + a labelled precomputed omics kernel (omics_kernel_label) made
# the RKHS / GBLUP_BRR output step fail ("object 'geno_data' not found"),
# which returned no predictions.
test_that("labelled omics kernels with geno_data return predictions", {
  skip_on_cran()
  set.seed(12)
  ids <- sprintf("L%02d", 1:30)
  G <- matrix(2L * stats::rbinom(30 * 80, 1, 0.5), 30, dimnames = list(ids, sprintf("m%02d", 1:80)))
  O <- matrix(stats::rnorm(30 * 20), 30, dimnames = list(ids, NULL))
  COP <- tcrossprod(scale(O)) / ncol(O) + diag(1e-4, 30)
  dimnames(COP) <- list(ids, ids)
  ph <- data.frame(GID = ids, Yield = stats::rnorm(30), stringsAsFactors = FALSE)
  ph$Yield[1:5] <- NA
  withr::local_dir(withr::local_tempdir())
  for (m in c("RKHS", "GBLUP_BRR")) {
    res <- suppressWarnings(suppressMessages(model_execute(
      pheno_data = ph, geno_data = G, omic1_kernel = COP,
      omics_kernel_label = list(omic1_kernel = "COP", omic2_kernel = NULL, omic3_kernel = NULL),
      response = "Yield", gen_name = "GID", fixed = ~1, random = ~GID, GS_model = m, gmatrix_method = "VanRaden",
      qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, nIter = 200L, burnIn = 50L, thin = 2L,
      parallel_mode = "sequential", system_database = TRUE, message = FALSE
    )))
    pv <- res$model_results$predicted_values
    expect_true(is.data.frame(pv) && nrow(pv) == 30L, info = m)
  }
})

test_that("intercept-only fixed formulas are accepted by dedicated routes", {
  expect_true(PredictProR:::gp_formula_is_intercept_only(NULL))
  expect_true(PredictProR:::gp_formula_is_intercept_only(~1))
  expect_false(PredictProR:::gp_formula_is_intercept_only(~Env))
  expect_false(PredictProR:::gp_formula_is_intercept_only(~0))
})

# Traits keep their own NA patterns (trait-specific test sets); added
# genotyped-only lines are test lines for every trait.
test_that("traits with different NA lines keep trait-specific test sets", {
  skip_on_cran()
  set.seed(9)
  ids <- sprintf("L%02d", 1:30)
  G <- matrix(2L * stats::rbinom(30 * 80, 1, 0.5), 30, dimnames = list(ids, sprintf("m%02d", 1:80)))
  ph <- data.frame(GID = ids[1:27], Yield = stats::rnorm(27), Protein = stats::rnorm(27), stringsAsFactors = FALSE)
  ph$Yield[1:3] <- NA                                        # Yield test lines L01-L03
  ph$Protein[4:7] <- NA                                      # Protein test lines L04-L07
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, response = c("Yield", "Protein"), gen_name = "GID", random = ~GID,
    GS_model = c("GBLUP_BRR", "GBLUP_BRR"), gmatrix_method = "VanRaden", qc_filtering = FALSE, impute = FALSE,
    ld_prunning_qc = FALSE, nIter = 200L, burnIn = 50L, thin = 2L, parallel_mode = "sequential",
    system_database = TRUE, message = FALSE
  )))
  geno_only <- ids[28:30]
  for (tr in c("Yield", "Protein")) {
    pv <- res$model_results_by_trait[[tr]][[1]]$predicted_values
    expect_setequal(pv$GID, ids)
    test_lines <- c(if (tr == "Yield") ids[1:3] else ids[4:7], geno_only)
    expect_setequal(pv$GID[pv$Train_Test_Label == "Test"], test_lines)
    observed <- ph$GID[!is.na(ph[[tr]])]
    expect_true(all(pv$Train_Test_Label[pv$GID %in% observed] == "Train"), info = tr)
    expect_equal(pv$Observed_value[match(observed, pv$GID)], ph[[tr]][match(observed, ph$GID)], info = tr)
  }
})

test_that("Bayesian MET true prediction returns the labelled grid", {
  skip_on_cran()
  set.seed(4)
  ids <- sprintf("L%02d", 1:24)
  G <- matrix(2L * stats::rbinom(24 * 80, 1, 0.5), 24, dimnames = list(ids, sprintf("m%02d", 1:80)))
  ph <- expand.grid(GID = ids, Env = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
  ph$Yield <- as.numeric(G[match(ph$GID, ids), 1:5] %*% stats::rnorm(5)) + stats::rnorm(nrow(ph))
  removed <- ph[c(2, 30, 55), ]
  ph <- ph[-c(2, 30, 55), ]                                  # three missing line x environment records
  ph$Yield[ph$GID %in% c("L01", "L02")] <- NA                # two test lines
  expected_unobserved <- sum(!removed$GID %in% c("L01", "L02"))
  withr::local_dir(withr::local_tempdir())
  res <- suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, response = "Yield", gen_name = "GID", heter_groups = "Env",
    random = ~ GID + GID:Env, GS_model = "GBLUP_BRR", gmatrix_method = "VanRaden",
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, nIter = 200L, burnIn = 50L, thin = 2L,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE
  )))
  pv <- res$model_results$predicted_values
  expect_identical(nrow(pv), 24L * 3L)
  expect_setequal(unique(pv$Train_Test_Label), c("Train", "Test", "Unobserved"))
  expect_true(all(pv$Train_Test_Label[pv$GID %in% c("L01", "L02")] == "Test"))
  expect_identical(sum(pv$Train_Test_Label == "Unobserved"), expected_unobserved)
  expect_true(all(is.finite(pv$Predicted_value)))
  # opt-out keeps the input records only
  res2 <- suppressWarnings(suppressMessages(model_execute(
    pheno_data = ph, geno_data = G, response = "Yield", gen_name = "GID", heter_groups = "Env",
    random = ~ GID + GID:Env, GS_model = "GBLUP_BRR", gmatrix_method = "VanRaden",
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE, nIter = 200L, burnIn = 50L, thin = 2L,
    parallel_mode = "sequential", system_database = TRUE, message = FALSE,
    met_predict_all_environments = FALSE
  )))
  expect_identical(nrow(res2$model_results$predicted_values), nrow(ph))
})
