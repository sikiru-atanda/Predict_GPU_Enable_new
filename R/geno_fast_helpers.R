predictpror_native_symbol_available <- function(symbol) {
  !is.null(tryCatch(
    getNativeSymbolInfo(symbol, PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

geno_qc_cpp_available <- function() {
  predictpror_native_symbol_available("predictpror_geno_qc_metrics")
}

geno_impute_summary_cpp_available <- function() {
  predictpror_native_symbol_available("predictpror_geno_impute_summary")
}

geno_impute_knn_cpp_available <- function() {
  predictpror_native_symbol_available("predictpror_geno_impute_knn")
}

geno_qc_metrics_cpp <- function(geno_data, ploidy = 2L) {
  if (!is.matrix(geno_data)) {
    geno_data <- as.matrix(geno_data)
  }
  ploidy <- gp_validate_ploidy(ploidy, allow_auto = FALSE)
  .Call("predictpror_geno_qc_metrics", geno_data, ploidy, PACKAGE = "PredictProR")
}

geno_impute_summary_cpp <- function(geno_data, method = c("mean", "median")) {
  method <- match.arg(method)
  if (!is.matrix(geno_data)) {
    geno_data <- as.matrix(geno_data)
  }
  .Call("predictpror_geno_impute_summary", geno_data, method, PACKAGE = "PredictProR")
}

geno_impute_knn_cpp <- function(geno_data, k = 5L, max_donors = NULL) {
  if (!is.matrix(geno_data)) {
    geno_data <- as.matrix(geno_data)
  }
  max_features <- suppressWarnings(as.integer(
    Sys.getenv("PREDICTPROR_KNN_MAX_DISTANCE_FEATURES", "2048")
  ))
  if (is.na(max_features)) {
    max_features <- 2048L
  }
  if (is.null(max_donors)) {
    max_donors <- suppressWarnings(as.integer(
      Sys.getenv("PREDICTPROR_KNN_MAX_DONORS", "512")
    ))
  }
  if (is.na(max_donors)) {
    max_donors <- 512L
  }
  .Call(
    "predictpror_geno_impute_knn",
    geno_data,
    as.integer(k),
    as.integer(max_features),
    as.integer(max_donors),
    PACKAGE = "PredictProR"
  )
}

geno_qc_metrics_r <- function(geno_data, ploidy = 2L) {
  if (!is.matrix(geno_data)) {
    geno_data <- as.matrix(geno_data)
  }
  list(
    monomorphic = as.integer(which(apply(geno_data, 2L, function(x) length(table(x)) <= 1L))),
    marker_missing_rate = colMeans(is.na(geno_data)),
    individual_missing_rate = rowMeans(is.na(geno_data)),
    maf = {
      phat <- colMeans(geno_data, na.rm = TRUE) / ploidy
      ifelse(phat < 0.5, phat, 1 - phat)
    },
    heterozygosity = apply(geno_data, 2L, function(x) {
      sum(x > 0 & x < ploidy, na.rm = TRUE) / length(x)
    })
  )
}

geno_qc_metrics <- function(geno_data, ploidy = "auto") {
  ploidy <- gp_resolve_matrix_ploidy(geno_data, ploidy)
  gp_validate_alt_dosage(geno_data, ploidy, hard_calls = FALSE)
  if (geno_qc_cpp_available()) {
    return(geno_qc_metrics_cpp(geno_data, ploidy = ploidy))
  }
  geno_qc_metrics_r(geno_data, ploidy = ploidy)
}
