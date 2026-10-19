# A task closure sends its whole frame to each PSOCK worker, and R writes an
# object once for every name or list element that holds it. The CV's frame
# held params and five of its elements under their own names as well
# (52.3 MB, 17 MB of it twice). gp_share_repeated_objects() sends a copy in
# which a repeat is a reference to the name holding it, filled in again on
# the worker before the task first runs.

.frame_with_repeats <- function() {
  local({
    set.seed(3)
    big <- matrix(stats::rnorm(200 * 200), 200, 200)
    other <- matrix(stats::rnorm(100 * 100), 100, 100)
    alias <- big
    params <- list(small = 1:3, big = big, nested = list(k = big, o = other, z = "z"), df = data.frame(a = 1:5))
    counter <- 0
    run <- function(i) {
      counter <<- counter + 1
      list(i = i, sum_big = sum(params$big), sum_nested = sum(params$nested$k), sum_other = sum(params$nested$o),
           alias_same = identical(alias, big), names = names(params), nested_names = names(params$nested),
           counter = counter)
    }
    init <- function() exists("params", inherits = TRUE)
    # without the test file's source records, so sizes are the data's
    run <- utils::removeSource(run)
    init <- utils::removeSource(init)
    list(run = run, init = init)
  })
}
# as an installed package sends it: under load_all() a package function keeps
# its file's text (the shared copy's entry is one)
.mb <- function(x) {
  if (is.function(x)) x <- utils::removeSource(x)
  length(serialize(x, NULL)) / 1024^2
}
# a closure over `data` whose frames end at the global environment
.closure_over <- function(data, body) {
  env <- new.env(parent = globalenv())
  env$data <- data
  eval(body, env)
}

test_that("a repeated object is sent once", {
  f <- .frame_with_repeats()
  shared <- gp_share_repeated_objects(f)
  big_mb <- .mb(get("big", environment(f$run)))
  other_mb <- .mb(get("other", environment(f$run)))
  # the frame's data, serialised in one stream as a worker gets it (the
  # function's other frames and source records depend on the session: in
  # the full Linux suite the whole copy measured 0.26 MB more than alone)
  frame_data_mb <- function(fun) {
    env <- environment(fun)
    if (is.function(get0("f", envir = env, inherits = FALSE))) env <- environment(env$f)  # the copy's entry
    vals <- mget(ls(env, all.names = TRUE), envir = env)
    .mb(Filter(Negate(is.function), vals))
  }
  expect_gt(frame_data_mb(f$run), 3.5 * big_mb)    # big, alias, params$big, params$nested$k
  # the three repeats of big and the one of other are gone
  expect_lt(frame_data_mb(shared$run), big_mb + other_mb + 0.02)
  expect_gt(frame_data_mb(shared$run), big_mb + other_mb - 0.01)
})

test_that("the task gives the same results after a round trip to a worker", {
  f <- .frame_with_repeats()
  expected <- f$run(7)
  f2 <- .frame_with_repeats()
  shared <- gp_share_repeated_objects(f2)
  # as a PSOCK worker gets it: serialized together, the frame once
  arrived <- unserialize(serialize(shared, NULL))
  expect_true(arrived$init())
  expect_identical(arrived$run(7), expected)
  expect_identical(arrived$run(8)$counter, 2)    # <<- still updates the task's frame
})

test_that("a repeated list holding a repeat is filled in after what it holds", {
  # `z` repeats `a` and is found before the repeat of `b` inside `a`
  f <- local({
    b <- matrix(stats::rnorm(5000), 50)
    a <- list(b, "small")
    z <- a
    list(run = function() list(z_is_a = identical(z, a), inner = identical(z[[1]], b), b = b))
  })
  expected <- f$run()
  arrived <- unserialize(serialize(gp_share_repeated_objects(f), NULL))
  expect_identical(arrived$run(), expected)
  expect_true(arrived$run()$z_is_a && arrived$run()$inner)
})

test_that("the caller's frame is left as it was", {
  f <- .frame_with_repeats()
  env <- environment(f$run)
  before <- as.list(env, all.names = TRUE)
  shared <- gp_share_repeated_objects(f)
  invisible(unserialize(serialize(shared, NULL))$run(1))
  after <- as.list(env, all.names = TRUE)
  expect_identical(names(after), names(before))
  for (nm in setdiff(names(before), c("run", "init"))) expect_identical(after[[nm]], before[[nm]], info = nm)
  expect_identical(f$run(1)$counter, 1)
})

test_that("the policy's size measure follows functions into their frames, each frame once", {
  big <- matrix(0, 1000, 1000)                                  # ~7.6 MB
  big_gb <- as.numeric(utils::object.size(big)) / 1024^3
  # a function holding the task, as gp_routed_task() and the shared copy's entry do
  inner <- .closure_over(big, quote(function(i) sum(data[i, ])))
  outer <- .closure_over(NULL, quote(function(i) run(i)))
  assign("run", inner, envir = environment(outer))
  expect_equal(gp_parallel_closure_payload_gb(outer), big_gb, tolerance = 0.01)
  # two functions of one frame: the frame counts once
  two <- .closure_over(big, quote(function(i) g()))
  assign("g", eval(quote(function() data), environment(two)), envir = environment(two))
  expect_equal(gp_parallel_closure_payload_gb(two), big_gb, tolerance = 0.01)
  # through the copy's entry: its data once, without the repeats
  f <- .frame_with_repeats()
  data_gb <- sum(vapply(c("big", "other"), function(nm) {
    as.numeric(utils::object.size(get(nm, envir = environment(f$run))))
  }, 0)) / 1024^3
  shared_gb <- gp_parallel_closure_payload_gb(gp_share_repeated_objects(f)$run)
  expect_gt(shared_gb, data_gb)
  expect_lt(shared_gb, data_gb * 1.1)
  expect_gt(gp_parallel_closure_payload_gb(f$run), 3.5 * data_gb * 0.8)
})

test_that("functions outside a closure frame are returned as they are", {
  expect_identical(gp_share_repeated_objects(list(run = stats::median))$run, stats::median)
  expect_identical(gp_share_repeated_objects(list(run = NULL, init = NULL)), list(run = NULL, init = NULL))
})
