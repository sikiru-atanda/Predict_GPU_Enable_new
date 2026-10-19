# ASReml finds a relationship inverse in vm(GID, name) by its name in the
# global environment. The multi-trait and hybrid routes named those objects
# with the process id (pp_mtmet_inv_1_6072), so the variance-component and
# FA tables of two runs never matched. The names are now fixed, and an
# object of the user's with the same name is put back afterwards.

.small_kernels <- function(n_kernels, n = 4) {
  ids <- paste0("L", seq_len(n))
  lapply(seq_len(n_kernels), function(k) {
    m <- diag(n) + 0.1 * k
    dimnames(m) <- list(ids, ids)
    m
  })
}

test_that("multi-trait ASReml inverses are named without the process id", {
  kernels <- stats::setNames(.small_kernels(2), c("G", "K"))
  specs <- gp_multitrait_asreml_assign_kernel_inverses(kernels, "pp_mtmet_inv", inverse = TRUE, epsilon = 1e-6)
  restore <- attr(specs, "restore")
  expect_true(is.function(restore))
  expect_identical(specs$inverse_name, c("pp_mtmet_inv_1", "pp_mtmet_inv_2"))
  expect_true(all(vapply(specs$inverse_name, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
  restore()
  expect_false(any(vapply(specs$inverse_name, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
})

test_that("with ten or more kernels no inverse name contains another", {
  # components are matched to kernels with grepl(name, fixed = TRUE)
  kernels <- stats::setNames(.small_kernels(10), paste0("K", 1:10))
  specs <- gp_multitrait_asreml_assign_kernel_inverses(kernels, "pp_mt_inv", inverse = FALSE, epsilon = 1e-6)
  on.exit(attr(specs, "restore")())
  nm <- specs$inverse_name
  expect_identical(nm[c(1, 10)], c("pp_mt_inv_01", "pp_mt_inv_10"))
  contained <- outer(nm, nm, Vectorize(function(a, b) a != b && grepl(a, b, fixed = TRUE)))
  expect_false(any(contained))
})

test_that("a user's object with an inverse's name is put back", {
  assign("pp_mtmet_inv_1", "the user's object", envir = .GlobalEnv)
  on.exit(if (exists("pp_mtmet_inv_1", envir = .GlobalEnv, inherits = FALSE)) rm("pp_mtmet_inv_1", envir = .GlobalEnv))
  specs <- gp_multitrait_asreml_assign_kernel_inverses(list(G = .small_kernels(1)[[1]]), "pp_mtmet_inv",
                                                       inverse = TRUE, epsilon = 1e-6)
  expect_false(is.character(get("pp_mtmet_inv_1", envir = .GlobalEnv)))
  attr(specs, "restore")()
  expect_identical(get("pp_mtmet_inv_1", envir = .GlobalEnv), "the user's object")
})

test_that("global objects for ASReml are assigned and put back as they were", {
  assign("female_gca_inv", 1:3, envir = .GlobalEnv)
  on.exit(suppressWarnings(rm("female_gca_inv", "male_gca_inv", envir = .GlobalEnv)))
  restore <- gp_asreml_assign_globals(list(female_gca_inv = "f", male_gca_inv = "m"))
  expect_identical(get("female_gca_inv", envir = .GlobalEnv), "f")
  expect_identical(get("male_gca_inv", envir = .GlobalEnv), "m")
  restore()
  expect_identical(get("female_gca_inv", envir = .GlobalEnv), 1:3)
  expect_false(exists("male_gca_inv", envir = .GlobalEnv, inherits = FALSE))
})

test_that("no ASReml route names its global objects with the process id", {
  files <- test_path("..", "..", "R", c("multi_trait_asreml_helpers.R", "hybrid_asreml_helpers.R"))
  skip_if(!all(file.exists(files)), "package sources not available")
  code <- unlist(lapply(files, function(f) sub("#.*$", "", readLines(f, warn = FALSE))))
  expect_false(any(grepl("Sys.getpid()", code, fixed = TRUE)))
})
