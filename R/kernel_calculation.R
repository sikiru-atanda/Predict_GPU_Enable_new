gp_kernel_dense_cpp_available <- function() {
  !is.null(tryCatch(
    getNativeSymbolInfo("predictpror_kernel_dense", PACKAGE = "PredictProR"),
    error = function(e) NULL
  ))
}

gp_kernel_backend <- function(backend = NULL, method = NULL) {
  if (is.null(backend)) {
    backend <- Sys.getenv("PREDICTPRO_KERNEL_BACKEND", "auto")
  }
  backend <- as.character(backend)
  if (!identical(length(backend), 1L) || !(backend %in% c("auto", "cpp", "r"))) {
    stop("backend must be one of: auto, cpp, r.", call. = FALSE)
  }
  if (identical(backend, "auto") && gp_kernel_dense_cpp_available()) {
    method <- if (is.null(method)) "" else as.character(method)
    native_auto_methods <- c(
      "Matern_kernel",
      "Matern12_kernel",
      "Matern32_kernel",
      "Matern52_kernel",
      "Laplacian_kernel",
      "RationalQuadratic_kernel"
    )
    if (method %in% native_auto_methods) {
      return("cpp")
    }
  }
  if (identical(backend, "auto")) {
    return("r")
  }
  if (identical(backend, "cpp") && !gp_kernel_dense_cpp_available()) {
    stop("PredictProR native dense kernel backend is not available.", call. = FALSE)
  }
  backend
}

#' Calculate Kernel Matrix for Omics Data
#'
#' This function computes a dense kernel matrix for omics data using Gaussian,
#' linear, composite, polynomial, Matern-family, Laplacian, and rational
#' quadratic kernels. By default, PredictProR uses the
#' fastest dense route for the selected method unless a backend is forced.
#'
#' @param M_matrix_clean A numeric matrix representing the cleaned omics or M matrix.
#' @param scaling Logical, indicating if the input matrix should be scaled. Defaults to TRUE.
#' @param centering Logical, indicating if the input matrix should be centered when `scaling` is FALSE.
#' @param theta Numeric, the theta parameter for the Gaussian kernel. Defaults to 1 if not provided.
#' @param alpha Numeric, the mixing parameter for the Composite kernel. Defaults to 0.5.
#' @param gamma Positive numeric shape parameter for the rational quadratic
#' kernel.
#' @param smoothness_parameter Positive numeric smoothness parameter for
#' `"Matern_kernel"`.
#' @param length_scale Optional positive numeric bandwidth/length scale for
#'   distance kernels. The default `NULL` uses the median nonzero Euclidean
#'   distance after the requested scaling/centering, avoiding dimension-driven
#'   near-identity kernels.
#' @param method Character vector specifying one or more kernel calculation
#' methods. Supported methods include `"Gaussian_kernel"`, `"Linear_kernel"`,
#' `"Composite_kernel"`, `"Poly2_kernel"`, `"Poly3_kernel"`,
#' `"Poly4_kernel"`, `"Matern_kernel"`, `"Matern12_kernel"`,
#' `"Matern32_kernel"`, `"Matern52_kernel"`, `"Laplacian_kernel"`, and
#' `"RationalQuadratic_kernel"`. Short aliases without `"_kernel"` are also
#' accepted.
#' @param message Logical, indicating if messages should be printed. Defaults to TRUE.
#' @param scale Optional legacy alias for `scaling`, kept for pipeline callers
#'   that already pass `scale =`.
#' @param center Optional legacy alias for `centering`.
#' @param backend Backend selector. One of `"auto"`, `"cpp"`, or `"r"`. `NULL` uses
#'   `PREDICTPRO_KERNEL_BACKEND` or `"auto"`. The automatic route is
#'   performance-aware and keeps BLAS-backed dense methods on the R path in
#'   this phase.
#' @param ... Additional arguments passed to the kernel calculation function.
#'
#' @return Returns a numeric matrix for one method, or a named list of matrices
#' when multiple methods are requested.
#'
#' @examples
#' M_matrix <- matrix(rnorm(100), ncol = 10)
#' kernel_matrix <- kernel_calculation(M_matrix_clean = M_matrix, method = "Gaussian_kernel")
#'
#' @export
kernel_calculation <- function(
    M_matrix_clean = NULL,
    scaling = TRUE,
    centering = FALSE,
    theta = NULL,
    alpha = 0.5,
    gamma = 1,
    smoothness_parameter = 1.5,
    length_scale = NULL,
    method = NULL,
    message = TRUE,
    scale = NULL,
    center = NULL,
    backend = NULL,
    ...) {

  msg <- ""

  if (!is.null(scale)) {
    scaling <- scale
  }
  if (!is.null(center)) {
    centering <- center
  }

  if (is.null(theta)) {
    theta <- 1
  }
  if (is.null(M_matrix_clean)) {
    stop(paste(msg, "object omics/M_matrix is missing."), call. = FALSE)
  }
  input_ploidy <- attr(M_matrix_clean, "ploidy", exact = TRUE)

  kernel_method_avaliable <- gp_kernel_supported_methods()
  method <- gp_normalize_kernel_methods(method)
  if (is.null(method) || !length(method) || anyNA(method) || !all(method %in% kernel_method_avaliable)) {
    stop(
      "Invalid kernel method. Choose from: ",
      paste(kernel_method_avaliable, collapse = ", "),
      call. = FALSE
    )
  }
  method <- unique(method)
  scalar_or_default <- function(x, default) {
    if (is.null(x)) {
      return(default)
    }
    as.numeric(x[[1L]])
  }
  length_scale_auto <- is.null(length_scale)
  length_scale <- if (isTRUE(length_scale_auto)) NA_real_ else scalar_or_default(length_scale, NA_real_)
  smoothness_parameter <- scalar_or_default(smoothness_parameter, 1.5)
  gamma <- scalar_or_default(gamma, 1)
  theta <- scalar_or_default(theta, 1)
  alpha <- scalar_or_default(alpha, 0.5)
  if (any(method %in% c("Gaussian_kernel", "Composite_kernel")) &&
      (!is.finite(theta) || theta <= 0)) {
    stop(paste(msg, "theta must be a positive finite value for Gaussian kernels."), call. = FALSE)
  }
  if ("Composite_kernel" %in% method &&
      (!is.finite(alpha) || alpha < 0 || alpha > 1)) {
    stop(paste(msg, "alpha must be a finite value in [0, 1] for Composite_kernel."), call. = FALSE)
  }
  length_scale_methods <- c(
    "Matern_kernel",
    "Matern12_kernel",
    "Matern32_kernel",
    "Matern52_kernel",
    "Laplacian_kernel",
    "RationalQuadratic_kernel"
  )
  if (!isTRUE(length_scale_auto) && any(method %in% length_scale_methods) &&
      (!is.finite(length_scale) || length_scale <= 0)) {
    stop(paste(msg, "length_scale must be a positive finite value."), call. = FALSE)
  }
  if ("Matern_kernel" %in% method && (!is.finite(smoothness_parameter) || smoothness_parameter <= 0)) {
    stop(paste(msg, "smoothness_parameter must be a positive finite value."), call. = FALSE)
  }
  if ("RationalQuadratic_kernel" %in% method && (!is.finite(gamma) || gamma <= 0)) {
    stop(paste(msg, "gamma must be a positive finite value."), call. = FALSE)
  }

  if (!inherits(M_matrix_clean, "matrix")) {
    M_matrix_clean <- as.matrix(M_matrix_clean)
    if (isTRUE(message)) {
      base::message(insight::print_color(paste(msg, "M_matrix is not class matrix. We fix it."), "blue"))
    }
  }
  if (!is.numeric(M_matrix_clean)) {
    storage.mode(M_matrix_clean) <- "double"
  }
  if (anyNA(M_matrix_clean)) {
    stop(paste(msg, "Missing value is not expected."), call. = FALSE)
  }

  if (isFALSE(scaling) && isTRUE(message)) {
    base::message(insight::print_color(
      paste(msg, "If data is not previously scaled, it is recommended you scale the data."),
      "blue"
    ))
  }
  if (isTRUE(scaling)) {
    M_matrix_clean <- scale(x = M_matrix_clean, center = TRUE, scale = TRUE)
  } else if (isTRUE(centering)) {
    M_matrix_clean <- scale(x = M_matrix_clean, center = TRUE, scale = FALSE)
  }

  if (anyNA(M_matrix_clean) || any(!is.finite(M_matrix_clean))) {
    bad_cols <- colnames(M_matrix_clean)[colSums(is.na(M_matrix_clean) | !is.finite(M_matrix_clean)) > 0]
    bad_msg <- if (length(bad_cols)) paste(utils::head(bad_cols, 10), collapse = ", ") else "unknown columns"
    stop(
      paste(msg, "Kernel input contains non-finite values after scaling; check constant columns:", bad_msg),
      call. = FALSE
    )
  }
  storage.mode(M_matrix_clean) <- "double"

  if (isTRUE(length_scale_auto)) {
    if (any(method %in% length_scale_methods)) {
      distance_values <- as.numeric(stats::dist(M_matrix_clean))
      distance_values <- distance_values[is.finite(distance_values) & distance_values > 0]
      if (!length(distance_values)) {
        stop(paste(msg, "Cannot infer a distance-kernel length scale from identical rows."), call. = FALSE)
      }
      length_scale <- stats::median(distance_values)
    } else {
      # The native kernel call accepts one common parameter bundle even when a
      # distance kernel was not requested; keep that unused slot finite.
      length_scale <- 1
    }
  }

  linear_cache <- NULL
  gaussian_cache <- NULL
  dist2_cache <- NULL
  dist_cache <- NULL
  l1_cache <- NULL
  linear_kernel <- function(x) {
    if (is.null(linear_cache)) {
      linear_cache <<- tcrossprod(x) / ncol(x)
    }
    linear_cache
  }
  gaussian_kernel <- function(x, theta) {
    if (is.null(gaussian_cache)) {
      dist2 <- squared_distance(x)
      median_dist2 <- stats::median(dist2)
      if (!is.finite(median_dist2) || median_dist2 <= 0) {
        stop(paste(msg, "Gaussian kernel median distance must be positive."), call. = FALSE)
      }
      gaussian_cache <<- exp(-theta * dist2 / median_dist2)
    }
    gaussian_cache
  }
  squared_distance <- function(x) {
    if (is.null(dist2_cache)) {
      dist2_cache <<- as.matrix(stats::dist(x))^2
    }
    dist2_cache
  }
  euclidean_distance <- function(x) {
    if (is.null(dist_cache)) {
      dist_cache <<- sqrt(squared_distance(x))
    }
    dist_cache
  }
  manhattan_distance <- function(x) {
    if (is.null(l1_cache)) {
      l1_cache <<- as.matrix(stats::dist(x, method = "manhattan"))
    }
    l1_cache
  }
  matern_general_kernel <- function(x, smoothness, scale_value) {
    d <- euclidean_distance(x)
    z <- sqrt(2 * smoothness) * d / scale_value
    out <- matrix(1, nrow = nrow(d), ncol = ncol(d), dimnames = dimnames(d))
    nonzero <- z > 0
    if (any(nonzero)) {
      out[nonzero] <- (2^(1 - smoothness) / base::gamma(smoothness)) *
        z[nonzero]^smoothness *
        exp(-z[nonzero]) *
        besselK(z[nonzero], nu = smoothness, expon.scaled = TRUE)
    }
    out
  }
  matern12_kernel <- function(x, scale_value) {
    exp(-euclidean_distance(x) / scale_value)
  }
  matern32_kernel <- function(x, scale_value) {
    d <- euclidean_distance(x) / scale_value
    z <- sqrt(3) * d
    (1 + z) * exp(-z)
  }
  matern52_kernel <- function(x, scale_value) {
    d <- euclidean_distance(x) / scale_value
    z <- sqrt(5) * d
    (1 + z + 5 * d^2 / 3) * exp(-z)
  }
  laplacian_kernel <- function(x, scale_value) {
    exp(-manhattan_distance(x) / scale_value)
  }
  rational_quadratic_kernel <- function(x, scale_value, rq_alpha) {
    (1 + squared_distance(x) / (2 * rq_alpha * scale_value^2))^(-rq_alpha)
  }

  calculate_one_kernel <- function(method_one) {
    selected_backend <- gp_kernel_backend(backend, method_one)
    degree <- switch(
      method_one,
      Poly2_kernel = 2L,
      Poly3_kernel = 3L,
      Poly4_kernel = 4L,
      1L
    )
    if (identical(selected_backend, "cpp")) {
      KRM <- .Call(
        "predictpror_kernel_dense",
        M_matrix_clean,
        as.character(method_one),
        as.numeric(theta),
        as.numeric(alpha),
        as.integer(degree),
        as.numeric(length_scale),
        as.numeric(smoothness_parameter),
        as.numeric(gamma),
        PACKAGE = "PredictProR"
      )
      rownames(KRM) <- rownames(M_matrix_clean)
      colnames(KRM) <- rownames(M_matrix_clean)
      if (!is.null(input_ploidy)) attr(KRM, "ploidy") <- input_ploidy
      return(KRM)
    }

    KRM <- switch(
      method_one,
      Gaussian_kernel = gaussian_kernel(M_matrix_clean, theta),
      Linear_kernel = linear_kernel(M_matrix_clean),
      Composite_kernel = {
        alpha * linear_kernel(M_matrix_clean) + (1 - alpha) * gaussian_kernel(M_matrix_clean, theta)
      },
      Poly2_kernel = (linear_kernel(M_matrix_clean) + 1)^2,
      Poly3_kernel = (linear_kernel(M_matrix_clean) + 1)^3,
      Poly4_kernel = (linear_kernel(M_matrix_clean) + 1)^4,
      Matern_kernel = matern_general_kernel(M_matrix_clean, smoothness_parameter, length_scale),
      Matern12_kernel = matern12_kernel(M_matrix_clean, length_scale),
      Matern32_kernel = matern32_kernel(M_matrix_clean, length_scale),
      Matern52_kernel = matern52_kernel(M_matrix_clean, length_scale),
      Laplacian_kernel = laplacian_kernel(M_matrix_clean, length_scale),
      RationalQuadratic_kernel = rational_quadratic_kernel(M_matrix_clean, length_scale, gamma),
      stop(paste(msg, "Select method to calculate kernel relationship matrix"), call. = FALSE)
    )

    KRM <- as.matrix(KRM)
    rownames(KRM) <- rownames(M_matrix_clean)
    colnames(KRM) <- rownames(M_matrix_clean)
    if (!is.null(input_ploidy)) attr(KRM, "ploidy") <- input_ploidy
    KRM
  }

  if (length(method) > 1L) {
    return(gp_named_method_list(lapply(method, calculate_one_kernel), method))
  }

  calculate_one_kernel(method)
}
