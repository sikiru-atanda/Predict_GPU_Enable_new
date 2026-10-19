# With cross-validation followed by the refit of the best models,
# Run_metadata used to describe only the refit (often one sequential task), so
# Run_metadata.csv said "sequential" for a CV that ran on mirai workers.

md <- function(...) {
  kv <- list(...)
  data.frame(key = names(kv), value = vapply(kv, as.character, character(1)), stringsAsFactors = FALSE)
}

test_that("the CV stage's backend, workers and decision are added as cv_ rows after the refit's policy rows", {
  refit <- md(context = "model_execute", active_python = "/py", dl_device = "cuda: 3 fits",
              policy_backend = "sequential", policy_workers = 1, policy_decision_reason = "no_parallel_eligible_tasks",
              task_unit = "cv_best_model_trait", cv_best_models_retained = "Yield:BayesC")
  cv <- md(context = "cross_validation", active_python = "/py", python_ml = "/py", runtime_num_gpus = 1,
           dl_device = "cuda: 2 fits", policy_backend = "mirai", policy_workers = 4,
           policy_decision_reason = "user_forced_mirai", policy_hidden_threads = 20, policy_parallel_count = 25,
           task_unit = "trait_model_replication", cv_generate_plots = FALSE, metric_for_ranking = "accuracy",
           cv_best_models_retained = "should not replace the refit's")
  out <- PredictProR:::gp_runtime_metadata_add_cv_stage(refit, cv)
  val <- function(k) out$value[out$key == k]

  expect_identical(val("policy_stage"), "best_model_refit")
  expect_identical(val("policy_backend"), "sequential")
  expect_identical(val("cv_policy_backend"), "mirai")
  expect_identical(val("cv_policy_workers"), "4")
  expect_identical(val("cv_policy_decision_reason"), "user_forced_mirai")
  expect_identical(val("cv_policy_parallel_count"), "25")
  expect_identical(val("cv_context"), "cross_validation")
  expect_identical(val("cv_task_unit"), "trait_model_replication")
  expect_identical(val("cv_metric_for_ranking"), "accuracy")
  expect_identical(val("cv_generate_plots"), "FALSE")
  expect_identical(val("cv_best_models_retained"), "Yield:BayesC")
  # session-wide fields are not repeated
  expect_false(any(c("cv_active_python", "cv_python_ml", "cv_runtime_num_gpus", "cv_dl_device",
                     "cv_policy_hidden_threads") %in% out$key))
  expect_identical(val("dl_device"), "cuda: 3 fits")
  expect_false(anyDuplicated(out$key) > 0L)
  # order: policy_stage, the refit's policy rows, then the CV rows, then the rest
  expect_identical(
    out$key,
    c("context", "active_python", "dl_device", "policy_stage", "policy_backend", "policy_workers",
      "policy_decision_reason", "cv_context", "cv_policy_backend", "cv_policy_workers",
      "cv_policy_decision_reason", "cv_policy_parallel_count", "cv_task_unit", "cv_generate_plots",
      "cv_metric_for_ranking", "task_unit", "cv_best_models_retained")
  )
})

test_that("without a CV stage the run metadata is unchanged", {
  refit <- md(context = "model_execute", policy_backend = "sequential")
  expect_identical(PredictProR:::gp_runtime_metadata_add_cv_stage(refit, NULL), refit)
  expect_identical(PredictProR:::gp_runtime_metadata_add_cv_stage(refit, refit[0, ]), refit)
  no_policy <- md(context = "model_execute", task_unit = "model_trait")
  out <- PredictProR:::gp_runtime_metadata_add_cv_stage(no_policy, md(policy_backend = "mirai"))
  expect_identical(out$key, c("context", "task_unit", "cv_policy_backend"))
})
