local_package_sources <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  r_dir <- file.path(root, "R")
  if (!dir.exists(r_dir)) {
    testthat::skip("raw R sources are not available in installed-package test context")
  }
  list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
}

call_head_name <- function(x) {
  if (is.symbol(x)) {
    as.character(x)
  } else {
    character()
  }
}

collect_top_level_function_defs <- function(exprs_by_file) {
  defs <- list()
  for (file in names(exprs_by_file)) {
    for (expr in exprs_by_file[[file]]) {
      if (!is.call(expr)) {
        next
      }
      op <- call_head_name(expr[[1]])
      if (length(op) != 1L || !(op %in% c("<-", "="))) {
        next
      }
      lhs <- expr[[2]]
      rhs <- expr[[3]]
      if (is.symbol(lhs) &&
          is.call(rhs) &&
          identical(call_head_name(rhs[[1]]), "function")) {
        sr <- attr(expr, "srcref")
        fmls <- names(as.list(rhs[[2]]))
        defs[[as.character(lhs)]] <- list(
          formals = fmls,
          has_dots = "..." %in% fmls,
          file = basename(file),
          line = if (!is.null(sr)) sr[[1]] else NA_integer_
        )
      }
    }
  }
  defs
}

find_unknown_local_call_args <- function(x, file, defs, out = character()) {
  if (!is.call(x)) {
    return(out)
  }

  fname <- call_head_name(x[[1]])
  if (length(fname) == 1L) {
    def <- defs[[fname]]
    if (!is.null(def) && !isTRUE(def$has_dots)) {
      arg_names <- names(as.list(x)[-1])
      if (is.null(arg_names)) {
        arg_names <- rep("", length(as.list(x)) - 1L)
      }
      unknown <- setdiff(arg_names[nzchar(arg_names)], def$formals)
      if (length(unknown)) {
        sr <- attr(x, "srcref")
        line <- if (!is.null(sr)) sr[[1]] else NA_integer_
        out <- c(
          out,
          sprintf(
            "%s:%s calls %s() with unknown argument(s): %s; signature at %s:%s",
            basename(file),
            line,
            fname,
            paste(unknown, collapse = ", "),
            def$file,
            def$line
          )
        )
      }
    }
  }

  for (i in seq_along(x)[-1]) {
    out <- find_unknown_local_call_args(x[[i]], file, defs, out)
  }
  out
}

test_that("package-local calls use current formal argument names", {
  files <- local_package_sources()
  exprs_by_file <- setNames(lapply(files, parse, keep.source = TRUE), files)
  defs <- collect_top_level_function_defs(exprs_by_file)

  issues <- character()
  for (file in names(exprs_by_file)) {
    for (expr in exprs_by_file[[file]]) {
      issues <- find_unknown_local_call_args(expr, file, defs, issues)
    }
  }

  expect_equal(issues, character())
})

test_that("results_handling call sites use run_metadata consistently", {
  expect_true("run_metadata" %in% names(formals(results_handling)))

  sources <- paste(
    vapply(
      local_package_sources(),
      function(file) paste(readLines(file, warn = FALSE), collapse = "\n"),
      character(1)
    ),
    collapse = "\n"
  )

  expect_false(grepl("runmeta_data|runmeta\\s*=", sources))
})
