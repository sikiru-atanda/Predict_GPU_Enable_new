#!/usr/bin/env Rscript

# Development-time reference benchmark for PredictProR's ploidy-aware dosage
# contract. AGHmatrix and polyBreedR are deliberately not package dependencies.

args <- commandArgs(trailingOnly = TRUE)
allow_download <- "--allow-download" %in% args ||
  identical(Sys.getenv("PREDICTPROR_ALLOW_POLYPLOID_SOURCE_DOWNLOAD"), "1")
output_arg <- grep("^--output=", args, value = TRUE)
output_file <- if (length(output_arg)) {
  sub("^--output=", "", output_arg[[1L]])
} else {
  file.path("tools", "outputs", "polyploid_reference_benchmark.csv")
}
tolerance <- 1e-10

benchmark_lib <- normalizePath(
  ".r-polyploid-benchmark-lib", winslash = "/", mustWork = FALSE
)
if (dir.exists(benchmark_lib)) .libPaths(c(benchmark_lib, .libPaths()))

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
} else if (!requireNamespace("PredictProR", quietly = TRUE)) {
  stop("Install PredictProR or run this script from its source checkout.", call. = FALSE)
}

pp_fun <- function(name) {
  ns <- asNamespace("PredictProR")
  get(name, envir = ns, inherits = FALSE)
}
grm <- pp_fun("grm_calculation")
parse_vcf <- pp_fun("gp_parse_vcf_alt_dosage")

metric_row <- function(comparator, dataset, ploidy, observed, reference,
                       comparator_version, note = "") {
  observed <- unclass(as.matrix(observed))
  reference <- unclass(as.matrix(reference))
  delta <- observed - reference
  max_abs <- max(abs(delta))
  rmse <- sqrt(mean(delta^2))
  correlation <- if (length(observed) > 1L) {
    suppressWarnings(stats::cor(c(observed), c(reference)))
  } else {
    NA_real_
  }
  data.frame(
    comparator = comparator,
    comparator_version = comparator_version,
    dataset = dataset,
    ploidy = as.integer(ploidy),
    n_individuals = nrow(observed),
    n_markers = NA_integer_,
    max_abs_difference = max_abs,
    rmse_difference = rmse,
    correlation = correlation,
    tolerance = tolerance,
    passed = is.finite(max_abs) && max_abs <= tolerance,
    note = note,
    stringsAsFactors = FALSE
  )
}

synthetic_dosage <- function(ploidy, seed, n = 24L, m = 10L) {
  set.seed(seed)
  x <- matrix(
    sample.int(ploidy + 1L, n * m, replace = TRUE) - 1L,
    nrow = n,
    dimnames = list(paste0("id", seq_len(n)), paste0("m", seq_len(m)))
  )
  storage.mode(x) <- "double"
  x
}

results <- list()
append_result <- function(row, n_markers) {
  row$n_markers <- as.integer(n_markers)
  results[[length(results) + 1L]] <<- row
}

if (!requireNamespace("AGHmatrix", quietly = TRUE)) {
  stop(
    "AGHmatrix is required only for this benchmark. Install it in ",
    ".r-polyploid-benchmark-lib; do not add it to PredictProR dependencies.",
    call. = FALSE
  )
}

agh_version <- as.character(utils::packageVersion("AGHmatrix"))
for (ploidy in c(2L, 4L, 6L)) {
  x <- synthetic_dosage(ploidy, 800L + ploidy)
  observed <- grm(x, method = "VanRaden", ploidy = ploidy, backend = "r")
  invisible(utils::capture.output(reference <- AGHmatrix::Gmatrix(
    SNPmatrix = x, method = "VanRaden", ploidy = ploidy,
    ploidy.correction = TRUE, maf = 0, verify.posdef = FALSE
  )))
  append_result(
    metric_row(
      "AGHmatrix::Gmatrix", paste0("synthetic_P", ploidy), ploidy,
      observed, reference, agh_version,
      "Parametric ploidy correction; identical ALT-dosage orientation."
    ),
    ncol(x)
  )
}

utils::data("snp.sol", package = "AGHmatrix", envir = environment())
observed <- grm(snp.sol, method = "VanRaden", ploidy = 4L, backend = "r")
invisible(utils::capture.output(reference <- AGHmatrix::Gmatrix(
  SNPmatrix = snp.sol, method = "VanRaden", ploidy = 4L,
  ploidy.correction = TRUE, maf = 0, verify.posdef = FALSE
)))
append_result(
  metric_row(
    "AGHmatrix::Gmatrix", "AGHmatrix_snp.sol", 4L,
    observed, reference, agh_version,
    "Public autotetraploid potato dosage data (571 x 3895)."
  ),
  ncol(snp.sol)
)

missing_x <- synthetic_dosage(4L, 844L)
missing_x[cbind(c(1L, 5L, 10L, 15L), c(1L, 3L, 5L, 7L))] <- NA_real_
imputed <- PredictProR::impute_genotypes_native(
  missing_x, ploidy = 4L, method = "mean"
)$snps_matrix
observed <- grm(imputed, method = "VanRaden", ploidy = 4L, backend = "r")
invisible(utils::capture.output(reference <- AGHmatrix::Gmatrix(
  SNPmatrix = missing_x, method = "VanRaden", ploidy = 4L,
  ploidy.correction = TRUE, impute.method = "mean", maf = 0,
  verify.posdef = FALSE
)))
append_result(
  metric_row(
    "AGHmatrix::Gmatrix", "synthetic_P4_missing_mean", 4L,
    observed, reference, agh_version,
    "PredictProR native marker-mean dosage imputation versus AGHmatrix mean imputation."
  ),
  ncol(missing_x)
)

polybreed_sha <- "91f35d689cba55b89997698b638d97e48dc6885d"
polybreed_version <- paste0("source@", substr(polybreed_sha, 1L, 12L))
polybreed_env <- NULL
if (requireNamespace("polyBreedR", quietly = TRUE)) {
  polybreed_g <- getExportedValue("polyBreedR", "G_mat")
  polybreed_gt <- getExportedValue("polyBreedR", "GT2DS")
  polybreed_version <- as.character(utils::packageVersion("polyBreedR"))
} else if (isTRUE(allow_download)) {
  polybreed_env <- new.env(parent = baseenv())
  for (source_file in c("G_mat.R", "GT2DS.R")) {
    local_file <- tempfile(fileext = ".R")
    source_url <- sprintf(
      "https://raw.githubusercontent.com/jendelman/polyBreedR/%s/R/%s",
      polybreed_sha, source_file
    )
    utils::download.file(source_url, local_file, quiet = TRUE, mode = "wb")
    sys.source(local_file, envir = polybreed_env)
  }
  polybreed_g <- get("G_mat", envir = polybreed_env, inherits = FALSE)
  polybreed_gt <- get("GT2DS", envir = polybreed_env, inherits = FALSE)
} else {
  stop(
    "polyBreedR is unavailable. Install it for benchmarking or rerun with ",
    "--allow-download to execute the pinned public G_mat.R and GT2DS.R sources.",
    call. = FALSE
  )
}

for (ploidy in c(2L, 4L, 6L)) {
  x <- synthetic_dosage(ploidy, 900L + ploidy)
  observed <- grm(x, method = "VanRaden", ploidy = ploidy, backend = "r")
  reference <- polybreed_g(t(x), ploidy = ploidy, method = "VR1")
  append_result(
    metric_row(
      "polyBreedR::G_mat", paste0("synthetic_P", ploidy), ploidy,
      observed, reference, polybreed_version,
      "VR1; polyBreedR input transposed to markers by individuals."
    ),
    ncol(x)
  )
}

observed <- grm(snp.sol, method = "VanRaden", ploidy = 4L, backend = "r")
reference <- polybreed_g(t(snp.sol), ploidy = 4L, method = "VR1")
append_result(
  metric_row(
    "polyBreedR::G_mat", "AGHmatrix_snp.sol", 4L,
    observed, reference, polybreed_version,
    "Public autotetraploid potato dosage data (571 x 3895)."
  ),
  ncol(snp.sol)
)

observed <- grm(imputed, method = "VanRaden", ploidy = 4L, backend = "r")
reference <- polybreed_g(t(missing_x), ploidy = 4L, method = "VR1")
append_result(
  metric_row(
    "polyBreedR::G_mat", "synthetic_P4_missing_mean", 4L,
    observed, reference, polybreed_version,
    "Both paths replace missing marker dosage by the marker mean before VR1."
  ),
  ncol(missing_x)
)

for (ploidy in c(2L, 3L, 4L, 6L)) {
  calls <- vapply(0:ploidy, function(dosage) {
    paste(c(rep("0", ploidy - dosage), rep("1", dosage)), collapse = "/")
  }, character(1L))
  samples <- paste0("s", 0:ploidy)
  vcf_row <- as.data.frame(as.list(c(
    `#CHROM` = "1", POS = "10", ID = "m1", REF = "A", ALT = "G",
    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
    stats::setNames(calls, samples)
  )), check.names = FALSE, stringsAsFactors = FALSE)
  observed_dosage <- as.numeric(parse_vcf(vcf_row, ploidy = ploidy)[1L, ])
  reference_dosage <- as.numeric(polybreed_gt(calls, diploidize = FALSE, n.core = 1L))
  append_result(
    metric_row(
      "polyBreedR::GT2DS", paste0("VCF_GT_P", ploidy), ploidy,
      matrix(observed_dosage, nrow = 1L),
      matrix(reference_dosage, nrow = 1L),
      polybreed_version,
      "Complete unphased biallelic GT calls covering every dosage class."
    ),
    1L
  )
}

results <- do.call(rbind, results)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output_file, row.names = FALSE, na = "")
print(results, row.names = FALSE)

if (!all(results$passed)) {
  stop("At least one public polyploid reference comparison exceeded tolerance.", call. = FALSE)
}
cat("\nAll ", nrow(results), " public-reference comparisons passed at tolerance ",
    format(tolerance, scientific = TRUE), ".\n", sep = "")
