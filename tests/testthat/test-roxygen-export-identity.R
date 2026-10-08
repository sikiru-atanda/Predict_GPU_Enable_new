local_roxygen_source_root <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  if (!file.exists(file.path(root, "NAMESPACE")) || !dir.exists(file.path(root, "R"))) {
    testthat::skip("raw package sources are not available in installed-package test context")
  }
  root
}

local_namespace_exports <- function(root) {
  ns <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  export_lines <- grep("^export\\(", ns, value = TRUE)
  unique(sub("^export\\(([^)]+)\\).*", "\\1", export_lines))
}

local_top_level_function_defs <- function(root) {
  files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
  defs <- list()
  for (file in files) {
    src <- readLines(file, warn = FALSE)
    hit <- grep("^\\s*[A-Za-z_.][A-Za-z0-9._]*\\s*(<-|=)\\s*function\\b", src)
    for (line in hit) {
      name <- sub("^\\s*([A-Za-z_.][A-Za-z0-9._]*)\\s*(<-|=)\\s*function\\b.*$", "\\1", src[[line]])
      defs[[name]] <- list(
        file = file,
        line = line
      )
    }
  }
  defs
}

local_has_attached_roxygen_block <- function(def) {
  if (is.na(def$line) || def$line <= 1L) {
    return(FALSE)
  }
  src <- readLines(def$file, warn = FALSE)
  line <- def$line - 1L
  while (line >= 1L && !nzchar(trimws(src[[line]]))) {
    line <- line - 1L
  }
  line >= 1L && grepl("^#'", src[[line]])
}

test_that("exported functions have directly attached roxygen identity blocks", {
  root <- local_roxygen_source_root()
  exported <- local_namespace_exports(root)
  defs <- local_top_level_function_defs(root)

  missing_defs <- setdiff(exported, names(defs))
  expect_equal(missing_defs, character())

  missing_roxygen <- exported[!vapply(exported, function(name) {
    local_has_attached_roxygen_block(defs[[name]])
  }, logical(1))]

  expect_equal(missing_roxygen, character())
})

test_that("pheno_geno_match is exported for public calls", {
  root <- local_roxygen_source_root()
  expect_true("pheno_geno_match" %in% local_namespace_exports(root))
})
