test_that("k-fold splitters respect stratified and unstratified semantics", {
  ph <- data.frame(
    GID = paste0("g", seq_len(12)),
    Yield = c(rep(0, 9), rep(1, 3)),
    stringsAsFactors = FALSE
  )

  stratified <- kfolds_stratified_un(
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    nfolds = 3,
    random_state = 2,
    sampling_method = "stratified",
    replication = 1,
    message = FALSE
  )[[1]]

  expect_equal(length(stratified), nrow(ph))
  expect_setequal(unique(stratified), 1:3)
  stratified_tab <- table(ph$Yield, stratified)
  expect_equal(as.integer(stratified_tab["0", ]), rep(3L, 3))
  expect_equal(as.integer(stratified_tab["1", ]), rep(1L, 3))

  unstratified <- kfolds_stratified_un(
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    nfolds = 3,
    random_state = 2,
    sampling_method = "unstratified",
    replication = 1,
    message = FALSE
  )[[1]]

  expect_equal(length(unstratified), nrow(ph))
  expect_equal(as.integer(table(unstratified)), rep(4L, 3))
  unstratified_tab <- table(ph$Yield, unstratified)
  expect_false(all(as.integer(unstratified_tab["1", ]) == 1L))

  legacy_unstratified <- CV_nfolds(
    pheno_data = ph,
    response = "Yield",
    gen_name = "GID",
    nFolds = 3,
    random_state = 2,
    method = "unstratified",
    replication = 1
  )[[1]]
  expect_equal(as.integer(table(legacy_unstratified)), rep(4L, 3))
})

test_that("holdout splitters return valid stratified and unstratified test indices", {
  ph <- data.frame(
    GID = paste0("g", seq_len(20)),
    Yield = c(rep(0, 10), rep(1, 10)),
    stringsAsFactors = FALSE
  )

  stratified <- hold_out_stratified_and_un(
    pheno_data = ph,
    gen_name = "GID",
    response = "Yield",
    test_size = 0.2,
    random_state = 11,
    replication = 1,
    sampling_method = "stratified",
    message = FALSE
  )[[1]]

  expect_equal(length(stratified), 4L)
  expect_equal(as.integer(table(ph$Yield[stratified])), c(2L, 2L))
  expect_true(all(stratified %in% seq_len(nrow(ph))))
  expect_equal(length(unique(stratified)), length(stratified))

  unstratified <- hold_out_stratified_and_un(
    pheno_data = ph,
    gen_name = "GID",
    response = "Yield",
    test_size = 0.2,
    random_state = 11,
    replication = 1,
    sampling_method = "unstratified",
    message = FALSE
  )[[1]]

  expect_equal(length(unstratified), 4L)
  expect_true(all(unstratified %in% seq_len(nrow(ph))))
  expect_equal(length(unique(unstratified)), length(unstratified))

  legacy_unstratified <- train_test_split(
    pheno_data = ph,
    gen_name = "GID",
    response = "Yield",
    test_size = 0.2,
    random_state = 11,
    replication = 1,
    method = "unstratified",
    message = FALSE
  )[[1]]
  expect_equal(length(legacy_unstratified), 4L)
  expect_true(all(legacy_unstratified %in% seq_len(nrow(ph))))
})

test_that("CV method normalization is shared across public and canonical paths", {
  allowed <- c(
    "Hold_Out",
    "Stratified_K-Folds",
    "Leave_one_Out",
    "CV2",
    "Repeated_CV2"
  )

  expect_identical(
    normalize_cv_token(c("Hold-Out", "stratified k folds", "cv2")),
    c("hold_out", "stratified_k_folds", "cv2")
  )
  expect_identical(gp_normalize_cv_method_name("hold-out", allowed), "Hold_Out")
  expect_identical(gp_normalize_cv_method_name("stratified k folds", allowed), "Stratified_K-Folds")
  expect_identical(gp_normalize_cv_method_name("leave-one-out", allowed), "Leave_one_Out")
  expect_identical(gp_resolve_cv_sampling_method("K-Folds"), "unstratified")
  expect_identical(gp_resolve_cv_sampling_method("Stratified_K-Folds"), "stratified")
  expect_identical(gp_resolve_cv_sampling_method("leave-one-out"), "unstratified")
})

test_that("CV0 CV1 and CV2 multi-environment assignments preserve their invariants", {
  met <- expand.grid(
    GID = paste0("g", 1:4),
    Env = paste0("E", 1:3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  met <- rbind(
    met,
    data.frame(GID = "g5", Env = "E1", stringsAsFactors = FALSE),
    data.frame(GID = NA_character_, Env = "E1", stringsAsFactors = FALSE)
  )
  met$Yield <- seq_len(nrow(met))
  keep <- stats::complete.cases(met[, c("GID", "Env")])

  cv0 <- CV0_CV1_CV2_for_multi_environment(
    pheno_data = met,
    gen_name = "GID",
    heter_groups = "Env",
    CV = 0,
    nfolds = 3,
    random_state = 5,
    replication = 1,
    message = FALSE
  )[[1]]

  expect_equal(length(cv0), nrow(met))
  expect_true(is.na(cv0[!keep]))
  expect_true(all(tapply(cv0[keep], met$Env[keep], function(x) length(unique(x))) == 1L))
  expect_setequal(unique(cv0[keep]), 1:3)

  cv1 <- CV0_CV1_CV2_for_multi_environment(
    pheno_data = met,
    gen_name = "GID",
    heter_groups = "Env",
    CV = 1,
    nfolds = 3,
    random_state = 5,
    replication = 1,
    message = FALSE
  )[[1]]

  expect_equal(length(cv1), nrow(met))
  expect_true(all(tapply(cv1[keep], met$GID[keep], function(x) length(unique(x))) == 1L))

  cv2 <- CV0_CV1_CV2_for_multi_environment(
    pheno_data = met,
    gen_name = "GID",
    heter_groups = "Env",
    CV = 2,
    nfolds = 3,
    random_state = 5,
    replication = 1,
    message = FALSE
  )[[1]]

  expect_equal(length(cv2), nrow(met))
  expect_true(is.na(cv2[!keep]))
  expect_equal(cv2[met$GID == "g5" & keep], 0L)
  multi_gid <- setdiff(unique(met$GID[keep]), "g5")
  expect_true(all(vapply(multi_gid, function(gid) {
    length(unique(cv2[met$GID == gid & keep])) > 1L
  }, logical(1))))
  expect_setequal(sort(unique(cv2[cv2 > 0])), 1:3)

  met_replicated <- rbind(
    met[keep, , drop = FALSE],
    met[which(met$GID == "g1" & met$Env == "E1")[1L], , drop = FALSE]
  )
  cv2_replicated <- CV0_CV1_CV2_for_multi_environment(
    pheno_data = met_replicated,
    gen_name = "GID",
    heter_groups = "Env",
    CV = 2,
    nfolds = 3,
    random_state = 5,
    replication = 1,
    message = FALSE
  )[[1]]
  g1e1 <- met_replicated$GID == "g1" & met_replicated$Env == "E1"
  expect_length(unique(cv2_replicated[g1e1]), 1L)

  disconnected <- data.frame(
    GID = paste0("g", 1:3),
    Env = paste0("E", 1:3),
    Yield = 1:3,
    stringsAsFactors = FALSE
  )
  expect_error(
    CV0_CV1_CV2_for_multi_environment(
      pheno_data = disconnected,
      gen_name = "GID",
      heter_groups = "Env",
      CV = 2,
      nfolds = 2,
      random_state = 5,
      replication = 1,
      message = FALSE
    ),
    "at least one genotype observed in at least two environments"
  )
})

test_that("cross-validation execution resolves methods and ignores CV2 train-only fold zero", {
  ph <- data.frame(
    GID = paste0("g", seq_len(8)),
    Env = rep(c("E1", "E2"), 4),
    Yield = seq_len(8),
    stringsAsFactors = FALSE
  )
  seen_sampling <- character()

  local_mocked_bindings(
    stop_daemons = function(...) NULL,
    gp_detect_python = function(...) NULL,
    gp_init_python_once = function(...) NULL,
    gp_lowrank_supported_models = function() character(),
    gp_parallel_extract_model_policy_params = function(...) list(),
    gp_execute_partitioned_tasks = function(model_ids, run_task, ...) {
      lapply(seq_along(model_ids), run_task)
    },
    gp_filter_crossval_results = function(results) results,
    kfolds_stratified_un = function(pheno_data, sampling_method, ...) {
      seen_sampling <<- c(seen_sampling, sampling_method)
      list(rep(1:2, length.out = nrow(pheno_data)))
    },
    CV0_CV1_CV2_for_multi_environment = function(pheno_data, ...) {
      list(c(1L, 0L, 2L, 1L, 0L, 2L, 1L, 2L))
    },
    predict_with_model = function(model = NULL, y = NULL, tst = NULL, additional_params = NULL, ...) {
      rep(mean(y, na.rm = TRUE), length(tst))
    },
    .package = "PredictProR"
  )

  base_params <- list(
    response = "Yield",
    gen_name = "GID",
    test_size = 0.25,
    random_state = 1,
    replication = 1,
    GS_model_cv = "RandomForest",
    nfolds = 2,
    eval_metrics = "mean_squared_error",
    response_family = "gaussian",
    parallel_mode = "sequential"
  )

  PredictProR:::models_execute_crossval(
    pheno_data = ph,
    params = modifyList(base_params, list(cross_validation_meth = "K-Folds")),
    verbose = FALSE
  )
  PredictProR:::models_execute_crossval(
    pheno_data = ph,
    params = modifyList(base_params, list(cross_validation_meth = "Stratified_K-Folds")),
    verbose = FALSE
  )

  expect_identical(seen_sampling, c("unstratified", "stratified"))

  cv2_res <- PredictProR:::models_execute_crossval(
    pheno_data = ph,
    params = modifyList(base_params, list(cross_validation_meth = "CV2", heter_groups = "Env")),
    verbose = FALSE
  )

  expect_equal(cv2_res[[1]]$cv_info$n_test, 6L)
  expect_equal(cv2_res[[1]]$cv_info$n_train, 2L)
  expect_true(all(cv2_res[[1]]$ypred_cv_Reps_all$cv_role[c(2, 5)] == "train"))

  loo_res <- PredictProR:::models_execute_crossval(
    pheno_data = ph,
    params = modifyList(base_params, list(cross_validation_meth = "leave-one-out")),
    verbose = FALSE
  )

  expect_equal(loo_res[[1]]$cv_info$token, "leave_one_out")
  expect_equal(loo_res[[1]]$cv_info$nfolds, nrow(ph))
  expect_equal(loo_res[[1]]$cv_info$n_test, nrow(ph))
})
