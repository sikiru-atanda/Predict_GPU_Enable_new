# ASReml cross-validation used to catch a fold's fit error, print it with cat()
# and return NULL: the error never reached the CV engine, a failed fold left
# NA predictions without a warning, and when every fold failed (for example the
# licence server refusing ASReml) model_execute() returned NULL with only
# generic warnings instead of stopping with the ASReml error.

licence_error <- "License error: License server does not support this product (-18)"

test_that("asreml_cv_model() passes a fit error on instead of returning NULL", {
  local_mocked_bindings(asreml_apply_options = function(...) invisible(NULL))
  ph <- data.frame(GID = paste0("g", 1:6), Yield = stats::rnorm(6))
  prep <- list(
    names_in_inv_list = character(), inv_list = list(),
    # a stand-in fitter in place of asreml::asreml(fixed = trait ~ 1, ...)
    code_asr_fit = c(sprintf("(function(...) stop('%s'))(fixed = trait ~ 1", licence_error),
                     "random = ~GID", "residual = ~units", "")
  )
  expect_error(
    PredictProR:::asreml_cv_model(pheno_dataa = ph, response = "Yield", gen_name = "GID",
                                  asreml_models_prep_cv = prep, tst = 1:2),
    "ASReml fit failed.*License error"
  )
})

asreml_cv_fold_data <- function() {
  set.seed(7)
  ids <- sprintf("L%02d", 1:30)
  X <- matrix(sample(0:2, 30 * 80, replace = TRUE), 30, 80, dimnames = list(ids, sprintf("m%02d", 1:80)))
  ph <- expand.grid(GID = ids, Env = c("E1", "E2"), stringsAsFactors = FALSE)
  ph$Yield <- 5 + as.numeric(scale(X[, 1:5]) %*% rnorm(5))[match(ph$GID, ids)] + rnorm(nrow(ph), sd = 0.5)
  list(pheno = ph, geno = X)
}

run_asreml_cv <- function(dat, mode = "sequential") {
  withr::local_dir(withr::local_tempdir())
  model_execute(
    pheno_data = dat$pheno, geno_data = dat$geno, gen_name = "GID", response = "Yield",
    # corgh: fa1 over these 2 environments is not identifiable (4 parameters for 3 entries)
    heter_groups = "Env", heter_resid = TRUE, var_cov_str = "corgh", fixed = ~Env, random = ~GID + GID:Env,
    GS_model_cv = "GBLUP", engine = "asreml", gmatrix_method = "VanRaden",
    cross_validation = TRUE, cross_validation_meth = "CV1", nfolds = 3, replication = 1,
    random_state = 11L, cv_evaluation_only = TRUE, cv_generate_plots = FALSE,
    eval_metrics = c("accuracy", "root_mean_squared_error"),
    qc_filtering = FALSE, impute = FALSE, ld_prunning_qc = FALSE,
    parallel_mode = mode, system_database = FALSE, message = FALSE, verbose = FALSE
  )
}

test_that("ASReml failing in every CV fold stops with the ASReml error", {
  skip_if_not_installed("asreml")
  local_mocked_bindings(asreml_mod_cv = function(...) stop("ASReml fit failed: ", licence_error, call. = FALSE))
  expect_error(suppressWarnings(run_asreml_cv(asreml_cv_fold_data())),
               "no result.*GBLUP.*License error")
})

test_that("the fold error survives fold-level CV tasks (split and merged)", {
  skip_if_not_installed("asreml")
  skip_if(parallel::detectCores() < 3L, "needs more than one usable core for fold tasks")
  # "auto" splits the CV into one task per fold; on this tiny data the cost
  # model runs them in-process, so the mock applies and the merge path is used
  withr::local_envvar(GP_CV_FOLD_TASKS = "true", GP_PAR_WORKER_OVERHEAD_GB = "0.5")
  local_mocked_bindings(asreml_mod_cv = function(...) stop("ASReml fit failed: ", licence_error, call. = FALSE))
  expect_error(suppressWarnings(run_asreml_cv(asreml_cv_fold_data(), mode = "auto")),
               "no result.*GBLUP.*License error")
})

test_that("ASReml failing in one CV fold warns with the error and keeps the other folds", {
  skip_if_not_installed("asreml")
  calls <- 0L
  local_mocked_bindings(asreml_mod_cv = function(..., tst) {
    calls <<- calls + 1L
    if (calls == 1L) stop("ASReml fit failed: ", licence_error, call. = FALSE)
    rep(5, length(tst))
  })
  warnings_seen <- character()
  out <- withCallingHandlers(
    run_asreml_cv(asreml_cv_fold_data()),
    warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("GBLUP \\(ASReml\\) failed.*fold.*License error", warnings_seen)))
  raw <- out$cv_results_raw[[1L]]$ypred_cv_Reps_all
  test_rows <- raw$cv_role == "test"
  expect_true(any(is.na(raw$yhat[test_rows])))      # the failed fold
  expect_true(any(!is.na(raw$yhat[test_rows])))     # the folds that worked
})
