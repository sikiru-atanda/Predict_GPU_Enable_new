test_that("shared owner produces a compact direct-map reference", {
  skip_if_not_installed("mori", minimum_version = "0.2.0")

  x <- matrix(seq_len(200), nrow = 20)
  owner <- parallel_shared_create(x)
  reference <- parallel_shared_reference(owner)
  mapped <- parallel_shared_map(reference)

  expect_s3_class(owner, "predictpror_shared_owner")
  expect_s3_class(reference, "predictpror_shared_reference")
  expect_true(mori::is_shared(owner$data))
  expect_true(mori::is_shared(mapped))
  expect_identical(mori::shared_name(mapped), reference$name)
  expect_equal(as.matrix(mapped), x)
  expect_false("data" %in% names(reference))
  expect_lt(length(serialize(reference, NULL)), 4096L)
  expect_true(reference$immutable)
  expect_true(reference$direct_final_worker_map)

  status <- parallel_shared_status(reference)
  expect_true(status$accessible)
  expect_false(status$owner_held_here)
  expect_identical(status$owner_pid, Sys.getpid())
  expect_error(
    parallel_shared_map(owner),
    "Do not send the shared owner to a worker",
    fixed = TRUE
  )
})

test_that("shared-data API rejects unsupported and invented handles", {
  skip_if_not_installed("mori", minimum_version = "0.2.0")

  expect_error(parallel_shared_create(NULL), "cannot be NULL", fixed = TRUE)
  expect_error(
    parallel_shared_create(new.env(parent = emptyenv())),
    "not a data type supported",
    fixed = TRUE
  )
  expect_error(
    parallel_shared_reference(list(name = "invented")),
    "must be created by parallel_shared_create",
    fixed = TRUE
  )
  expect_error(
    parallel_shared_map(list(name = "invented")),
    "must be created by parallel_shared_reference",
    fixed = TRUE
  )
})

fake_shared_reference <- function(size_gib = 0.25) {
  bytes <- size_gib * 1024^3
  structure(list(
    name = "Local\\predictpror_test_reference",
    bytes = bytes,
    size_mib = bytes / 1024^2,
    created_pid = Sys.getpid(),
    created_at = "2026-09-03T00:00:00-0500",
    immutable = TRUE,
    direct_final_worker_map = TRUE,
    schema_version = 1L
  ), class = c("predictpror_shared_reference", "list"))
}

test_that("batch plan requires an explicit data ownership mode", {
  expect_error(
    parallel_batch_plan(
      objects = list(x = numeric(10)), jobs = 2, workers_per_job = 1
    ),
    "data_mode must be explicitly selected",
    fixed = TRUE
  )
  expect_error(
    parallel_batch_plan(
      objects = list(x = numeric(10)), jobs = 2, workers_per_job = 1,
      data_mode = "automatic"
    ),
    "should be one of",
    fixed = TRUE
  )
})

test_that("batch projections distinguish independent and shared ownership", {
  skip_if_not_installed("mori", minimum_version = "0.2.0")
  object <- fake_shared_reference(0.25)
  make_plan <- function(mode) {
    parallel_batch_plan(
      objects = list(genotypes = object),
      jobs = 4L,
      workers_per_job = 2L,
      tasks_per_job = 4L,
      data_mode = mode,
      memory_budget_gb = 20,
      total_worker_limit = 8L,
      payload_multiplier = 1,
      worker_overhead_gb = 0.1
    )
  }

  independent <- make_plan("independent")
  per_job <- make_plan("per_job_shared")
  supervisor <- make_plan("supervisor_shared")

  expect_s3_class(supervisor, "predictpror_batch_plan")
  expect_equal(independent$decision$projected_batch_memory_gb, 4.2,
               tolerance = 1e-12)
  expect_equal(per_job$decision$projected_batch_memory_gb, 3.2,
               tolerance = 1e-12)
  expect_equal(supervisor$decision$projected_batch_memory_gb, 1.8,
               tolerance = 1e-12)
  expect_gt(independent$decision$projected_batch_memory_gb,
            per_job$decision$projected_batch_memory_gb)
  expect_gt(per_job$decision$projected_batch_memory_gb,
            supervisor$decision$projected_batch_memory_gb)
  expect_identical(
    supervisor$job_environment$data_handoff,
    rep("reference_only_final_worker_maps", 4L)
  )
})

test_that("batch plan queues jobs rather than oversubscribing worker slots", {
  skip_if_not_installed("mori", minimum_version = "0.2.0")
  plan <- parallel_batch_plan(
    objects = list(x = fake_shared_reference(0.01)),
    jobs = 4L,
    workers_per_job = 2L,
    tasks_per_job = 8L,
    data_mode = "supervisor_shared",
    memory_budget_gb = 20,
    total_worker_limit = 4L,
    payload_multiplier = 1,
    worker_overhead_gb = 0.1
  )

  expect_identical(plan$decision$jobs_admitted, 2L)
  expect_identical(plan$decision$jobs_queued, 2L)
  expect_identical(plan$decision$total_worker_slots_admitted, 4L)
  expect_lte(plan$decision$oversubscription_ratio, 1)
  expect_match(plan$decision$decision_reason, "jobs_capped_by_total_limit")
  expect_equal(nrow(plan$job_environment), 2L)
  expect_true(all(plan$job_environment$GP_BLAS_THREADS == "1"))
})

test_that("batch plan reduces workers for memory and can refuse admission", {
  object <- fake_shared_reference(0.25)
  reduced <- parallel_batch_plan(
    objects = list(x = object),
    jobs = 1L,
    workers_per_job = 4L,
    tasks_per_job = 8L,
    data_mode = "independent",
    memory_budget_gb = 1.6,
    total_worker_limit = 8L,
    payload_multiplier = 1,
    worker_overhead_gb = 0.1
  )
  expect_identical(reduced$decision$workers_per_job_selected, 3L)
  expect_true(reduced$decision$can_run)
  expect_match(reduced$decision$decision_reason, "workers_capped_by_memory")
  expect_lte(reduced$decision$projected_batch_memory_gb, 1.6)

  skip_if_not_installed("mori", minimum_version = "0.2.0")
  refused <- parallel_batch_plan(
    objects = list(x = object),
    jobs = 2L,
    workers_per_job = 2L,
    tasks_per_job = 4L,
    data_mode = "supervisor_shared",
    memory_budget_gb = 0.7,
    total_worker_limit = 4L,
    payload_multiplier = 1,
    worker_overhead_gb = 0.1
  )
  expect_false(refused$decision$can_run)
  expect_identical(refused$decision$jobs_admitted, 0L)
  expect_identical(refused$decision$jobs_queued, 2L)
  expect_equal(nrow(refused$job_environment), 0L)
  expect_match(refused$decision$decision_reason, "jobs_capped_by_memory")
})

test_that("batch planning is inspectable and does not mutate runtime settings", {
  before <- Sys.getenv(c("GP_PAR_MEMORY_BUDGET_GB", "GP_BLAS_THREADS", "GP_USE_MORI"),
                       unset = NA_character_)
  plan <- parallel_batch_plan(
    objects = list(x = numeric(100)),
    jobs = 2L,
    workers_per_job = 1L,
    data_mode = "independent",
    memory_budget_gb = 2,
    total_worker_limit = 2L,
    payload_multiplier = 1,
    worker_overhead_gb = 0.1
  )
  after <- Sys.getenv(names(before), unset = NA_character_)

  expect_identical(after, before)
  expect_named(plan, c("inventory", "decision", "job_environment"))
  expect_named(
    plan$job_environment,
    c("job_id", "GP_PAR_MEMORY_BUDGET_GB", "GP_BLAS_THREADS", "GP_USE_MORI",
      "data_handoff")
  )
  expect_output(print(plan), "jobs: 2 admitted / 2 requested")
})
