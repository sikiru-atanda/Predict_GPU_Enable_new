#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#ifndef FCONE
#define FCONE
#endif

#include <cmath>
#include <cstring>
#include <exception>
#include <new>
#include <vector>
#include <algorithm>

enum KernelMethod {
  KERNEL_GAUSSIAN,
  KERNEL_LINEAR,
  KERNEL_COMPOSITE,
  KERNEL_POLY,
  KERNEL_MATERN,
  KERNEL_MATERN12,
  KERNEL_MATERN32,
  KERNEL_MATERN52,
  KERNEL_LAPLACIAN,
  KERNEL_RATIONAL_QUADRATIC
};

static double scalar_finite_numeric(SEXP x, const char* name) {
  if (!Rf_isNumeric(x) || Rf_length(x) != 1) {
    Rf_error("%s must be a finite numeric scalar", name);
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, REALSXP));
  double value = REAL(coerced)[0];
  UNPROTECT(1);
  if (!R_FINITE(value)) {
    Rf_error("%s must be a finite numeric scalar", name);
  }
  return value;
}

static int scalar_integer(SEXP x, const char* name) {
  if (!Rf_isInteger(x) || Rf_length(x) != 1) {
    Rf_error("%s must be an integer scalar", name);
  }

  int value = INTEGER(x)[0];
  if (value == NA_INTEGER) {
    Rf_error("%s must be an integer scalar", name);
  }
  return value;
}

static KernelMethod parse_method(SEXP method_sexp) {
  if (!Rf_isString(method_sexp) || Rf_length(method_sexp) != 1 ||
      STRING_ELT(method_sexp, 0) == NA_STRING) {
    Rf_error("method must be a string scalar");
  }

  const char* method = CHAR(STRING_ELT(method_sexp, 0));
  if (std::strcmp(method, "Gaussian_kernel") == 0) {
    return KERNEL_GAUSSIAN;
  }
  if (std::strcmp(method, "Linear_kernel") == 0) {
    return KERNEL_LINEAR;
  }
  if (std::strcmp(method, "Composite_kernel") == 0) {
    return KERNEL_COMPOSITE;
  }
  if (std::strcmp(method, "Poly2_kernel") == 0 ||
      std::strcmp(method, "Poly3_kernel") == 0 ||
      std::strcmp(method, "Poly4_kernel") == 0) {
    return KERNEL_POLY;
  }
  if (std::strcmp(method, "Matern_kernel") == 0) {
    return KERNEL_MATERN;
  }
  if (std::strcmp(method, "Matern12_kernel") == 0) {
    return KERNEL_MATERN12;
  }
  if (std::strcmp(method, "Matern32_kernel") == 0) {
    return KERNEL_MATERN32;
  }
  if (std::strcmp(method, "Matern52_kernel") == 0) {
    return KERNEL_MATERN52;
  }
  if (std::strcmp(method, "Laplacian_kernel") == 0) {
    return KERNEL_LAPLACIAN;
  }
  if (std::strcmp(method, "RationalQuadratic_kernel") == 0) {
    return KERNEL_RATIONAL_QUADRATIC;
  }

  Rf_error("unsupported kernel method");
  return KERNEL_LINEAR;
}

// Gram matrix G = X X' (upper triangle) with one BLAS call; x is column-major n x p.
static void gram_upper(const double* x, int n, int p, double* g) {
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)("U", "N", &n, &p, &one, x, &n, &zero, g, &n FCONE FCONE);
}

// Squared Euclidean distance from the Gram matrix: |a|^2 + |b|^2 - 2 a.b.
static inline double gram_squared_distance(const double* g, int n, int row, int col) {
  if (row == col) return 0.0;
  double d2 = g[row + row * static_cast<R_xlen_t>(n)] + g[col + col * static_cast<R_xlen_t>(n)] -
    2.0 * g[row + col * static_cast<R_xlen_t>(n)];
  return d2 > 0.0 ? d2 : 0.0;
}

// L1 distance between two rows stored contiguously (row-major copy).
static inline double contiguous_l1_distance(const double* xr, int p, int row_a, int row_b) {
  const double* a = xr + static_cast<R_xlen_t>(row_a) * p;
  const double* b = xr + static_cast<R_xlen_t>(row_b) * p;
  double value = 0.0;
  for (int col = 0; col < p; ++col) {
    value += std::fabs(a[col] - b[col]);
  }
  return value;
}

static double matern_general(double distance,
                             double smoothness,
                             double length_scale) {
  if (distance <= 0.0) {
    return 1.0;
  }
  double z = std::sqrt(2.0 * smoothness) * distance / length_scale;
  double scaled_bessel = bessel_k(z, smoothness, 1.0);
  double value = std::pow(2.0, 1.0 - smoothness) /
    std::tgamma(smoothness) *
    std::pow(z, smoothness) *
    scaled_bessel;
  if (!R_FINITE(value)) {
    Rf_error("Matern kernel produced a non-finite value");
  }
  return value;
}

// Median of all n x n squared distances (diagonal zeros and both triangles),
// from the n(n-1)/2 upper-triangle values: the full sorted multiset is n zeros
// followed by every off-diagonal value twice. Linear-time selection.
static double median_distance_squared(const double* g, int n) {
  std::vector<double> off;
  off.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2);
  for (int col = 1; col < n; ++col) {
    for (int row = 0; row < col; ++row) {
      off.push_back(gram_squared_distance(g, n, row, col));
    }
  }
  const std::size_t total = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
  auto full_at = [&](std::size_t pos) -> double {
    if (pos < static_cast<std::size_t>(n)) {
      // zeros come first, but only while no off-diagonal value is negative (never)
      return 0.0;
    }
    std::size_t k = (pos - static_cast<std::size_t>(n)) / 2;
    std::nth_element(off.begin(), off.begin() + k, off.end());
    return off[k];
  };
  const std::size_t mid = total / 2;
  if (total % 2 == 0) {
    return 0.5 * (full_at(mid - 1) + full_at(mid));
  }
  return full_at(mid);
}

static SEXP predictpror_kernel_dense_impl(SEXP xSEXP,
                                          SEXP methodSEXP,
                                          SEXP thetaSEXP,
                                          SEXP alphaSEXP,
                                          SEXP degreeSEXP,
                                          SEXP lengthScaleSEXP,
                                          SEXP smoothnessSEXP,
                                          SEXP gammaSEXP) {

  if (!Rf_isReal(xSEXP) || !Rf_isMatrix(xSEXP)) {
    Rf_error("x must be a numeric matrix");
  }

  SEXP dim = Rf_getAttrib(xSEXP, R_DimSymbol);
  if (!Rf_isInteger(dim) || Rf_length(dim) != 2) {
    Rf_error("x must be a matrix with two dimensions");
  }

  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (n <= 0 || p <= 0) {
    Rf_error("x must have positive dimensions");
  }

  KernelMethod method = parse_method(methodSEXP);
  double theta = 1.0;
  double alpha = 0.5;
  double length_scale = 1.0;
  double smoothness = 1.5;
  double rq_alpha = 1.0;
  int degree = 1;
  if (method == KERNEL_GAUSSIAN || method == KERNEL_COMPOSITE) {
    theta = scalar_finite_numeric(thetaSEXP, "theta");
  }
  if (method == KERNEL_COMPOSITE) {
    alpha = scalar_finite_numeric(alphaSEXP, "alpha");
  }
  if (method == KERNEL_POLY) {
    degree = scalar_integer(degreeSEXP, "degree");
  }
  bool needs_length_scale =
    method == KERNEL_MATERN ||
    method == KERNEL_MATERN12 ||
    method == KERNEL_MATERN32 ||
    method == KERNEL_MATERN52 ||
    method == KERNEL_LAPLACIAN ||
    method == KERNEL_RATIONAL_QUADRATIC;
  if (needs_length_scale) {
    length_scale = scalar_finite_numeric(lengthScaleSEXP, "length_scale");
    if (length_scale <= 0.0) {
      Rf_error("length_scale must be positive");
    }
  }
  if (method == KERNEL_MATERN) {
    smoothness = scalar_finite_numeric(smoothnessSEXP, "smoothness_parameter");
    if (smoothness <= 0.0) {
      Rf_error("smoothness_parameter must be positive");
    }
  }
  if (method == KERNEL_RATIONAL_QUADRATIC) {
    rq_alpha = scalar_finite_numeric(gammaSEXP, "gamma");
    if (rq_alpha <= 0.0) {
      Rf_error("gamma must be positive");
    }
  }

  const double* x = REAL(xSEXP);
  for (R_xlen_t idx = 0; idx < XLENGTH(xSEXP); ++idx) {
    if (!R_FINITE(x[idx])) {
      Rf_error("x must contain only finite values");
    }
  }

  SEXP outSEXP = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(outSEXP);

  // Laplacian needs L1 distances (no BLAS form): rows copied contiguously so
  // each pair reads two cache-friendly blocks. Every other kernel works from
  // the Gram matrix X X', computed once with BLAS into `out` (upper triangle)
  // and turned into the kernel in place. Scratch buffers use R_alloc so an
  // interrupt (R_CheckUserInterrupt) cannot leak them.
  double* xr = nullptr;
  if (method == KERNEL_LAPLACIAN) {
    xr = reinterpret_cast<double*>(
      R_alloc(static_cast<std::size_t>(n) * static_cast<std::size_t>(p), sizeof(double)));
    for (int col = 0; col < p; ++col) {
      for (int row = 0; row < n; ++row) {
        xr[static_cast<std::size_t>(row) * p + col] = x[row + static_cast<R_xlen_t>(col) * n];
      }
    }
  } else {
    gram_upper(x, n, p, out);
  }

  double median_d2 = 0.0;
  bool needs_gaussian = method == KERNEL_GAUSSIAN || method == KERNEL_COMPOSITE;
  if (needs_gaussian) {
    median_d2 = median_distance_squared(out, n);
    if (!R_FINITE(median_d2) || median_d2 <= 0.0) {
      UNPROTECT(1);
      Rf_error("median squared distance must be positive and finite");
    }
  }

  // Squared norms are needed after the diagonal is overwritten.
  double* norms = nullptr;
  if (method != KERNEL_LAPLACIAN) {
    norms = reinterpret_cast<double*>(R_alloc(static_cast<std::size_t>(n), sizeof(double)));
    for (int i = 0; i < n; ++i) norms[i] = out[i + static_cast<R_xlen_t>(i) * n];
  }
  auto squared_distance = [&](int row, int col) -> double {
    if (row == col) return 0.0;
    double d2 = norms[row] + norms[col] - 2.0 * out[row + static_cast<R_xlen_t>(col) * n];
    return d2 > 0.0 ? d2 : 0.0;
  };

  // The Gram entry (row, col), row <= col, is read and then replaced by the
  // kernel value; the mirrored write goes to the strict lower triangle, which
  // holds no Gram data, and diagonal Gram values were saved in `norms`.
  for (int col = n - 1; col >= 0; --col) {
    if ((col & 63) == 0) {
      R_CheckUserInterrupt();
    }

    for (int row = 0; row <= col; ++row) {
      double gram = (method == KERNEL_LAPLACIAN) ? 0.0 : out[row + static_cast<R_xlen_t>(col) * n];
      double linear = gram / static_cast<double>(p);
      double value = linear;

      if (method == KERNEL_GAUSSIAN || method == KERNEL_COMPOSITE) {
        double d2 = squared_distance(row, col);
        double gaussian = std::exp(-theta * d2 / median_d2);
        value = gaussian;
        if (method == KERNEL_COMPOSITE) {
          value = alpha * linear + (1.0 - alpha) * gaussian;
        }
      } else if (method == KERNEL_POLY) {
        value = std::pow(linear + 1.0, static_cast<double>(degree));
      } else if (method == KERNEL_MATERN ||
                 method == KERNEL_MATERN12 ||
                 method == KERNEL_MATERN32 ||
                 method == KERNEL_MATERN52 ||
                 method == KERNEL_RATIONAL_QUADRATIC) {
        double d2 = squared_distance(row, col);
        if (method == KERNEL_RATIONAL_QUADRATIC) {
          value = std::pow(
            1.0 + d2 / (2.0 * rq_alpha * length_scale * length_scale),
            -rq_alpha
          );
        } else {
          double d = std::sqrt(d2);
          if (method == KERNEL_MATERN) {
            value = matern_general(d, smoothness, length_scale);
          } else if (method == KERNEL_MATERN12) {
            value = std::exp(-d / length_scale);
          } else if (method == KERNEL_MATERN32) {
            double z = std::sqrt(3.0) * d / length_scale;
            value = (1.0 + z) * std::exp(-z);
          } else {
            double d_scaled = d / length_scale;
            double z = std::sqrt(5.0) * d_scaled;
            value = (1.0 + z + 5.0 * d_scaled * d_scaled / 3.0) * std::exp(-z);
          }
        }
      } else if (method == KERNEL_LAPLACIAN) {
        double d1 = contiguous_l1_distance(xr, p, row, col);
        value = std::exp(-d1 / length_scale);
      }

      out[row + static_cast<R_xlen_t>(col) * n] = value;
      out[col + static_cast<R_xlen_t>(row) * n] = value;
    }
  }

  UNPROTECT(1);
  return outSEXP;
}

extern "C" SEXP predictpror_kernel_dense(SEXP xSEXP,
                                         SEXP methodSEXP,
                                         SEXP thetaSEXP,
                                         SEXP alphaSEXP,
                                         SEXP degreeSEXP,
                                         SEXP lengthScaleSEXP,
                                         SEXP smoothnessSEXP,
                                         SEXP gammaSEXP) {
  try {
    return predictpror_kernel_dense_impl(
      xSEXP,
      methodSEXP,
      thetaSEXP,
      alphaSEXP,
      degreeSEXP,
      lengthScaleSEXP,
      smoothnessSEXP,
      gammaSEXP
    );
  } catch (const std::bad_alloc&) {
    Rf_error("PredictProR native dense kernel backend ran out of memory.");
  } catch (const std::exception& ex) {
    Rf_error("PredictProR native dense kernel backend failed: %s", ex.what());
  } catch (...) {
    Rf_error("PredictProR native dense kernel backend failed.");
  }
  return R_NilValue;
}

// Median of the positive pairwise Euclidean distances between rows (the
// automatic length scale of the distance kernels); same definition as
// median(dist(x)[dist(x) > 0]) but from one BLAS Gram product.
static SEXP predictpror_median_positive_distance_impl(SEXP xSEXP) {
  if (!Rf_isReal(xSEXP) || !Rf_isMatrix(xSEXP)) {
    Rf_error("x must be a numeric matrix");
  }
  SEXP dim = Rf_getAttrib(xSEXP, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (n < 2 || p <= 0) {
    return Rf_ScalarReal(NA_REAL);
  }
  const double* x = REAL(xSEXP);
  // R_alloc: released by R even when R_CheckUserInterrupt jumps out
  const std::size_t nn = static_cast<std::size_t>(n);
  double* g = reinterpret_cast<double*>(R_alloc(nn * nn, sizeof(double)));
  gram_upper(x, n, p, g);
  double* d = reinterpret_cast<double*>(R_alloc(nn * (nn - 1) / 2, sizeof(double)));
  std::size_t count = 0;
  // Identical rows give |a|^2 + |b|^2 - 2 a.b = 0 only up to rounding, which an
  // optimized BLAS may leave a few ulps above zero; such values count as zero
  // (stats::dist on the rows gives exactly zero).
  const double rel_tol = 1e-10;
  for (int col = 1; col < n; ++col) {
    if ((col & 255) == 0) R_CheckUserInterrupt();
    for (int row = 0; row < col; ++row) {
      double d2 = gram_squared_distance(g, n, row, col);
      double scale = g[row + row * static_cast<R_xlen_t>(n)] + g[col + col * static_cast<R_xlen_t>(n)];
      if (R_FINITE(d2) && d2 > rel_tol * scale) d[count++] = std::sqrt(d2);
    }
  }
  if (count == 0) {
    return Rf_ScalarReal(NA_REAL);
  }
  std::size_t mid = count / 2;
  std::nth_element(d, d + mid, d + count);
  double upper = d[mid];
  if (count % 2 == 0) {
    double lower = *std::max_element(d, d + mid);
    return Rf_ScalarReal(0.5 * (lower + upper));
  }
  return Rf_ScalarReal(upper);
}

extern "C" SEXP predictpror_median_positive_distance(SEXP xSEXP) {
  try {
    return predictpror_median_positive_distance_impl(xSEXP);
  } catch (const std::bad_alloc&) {
    Rf_error("PredictProR native median distance ran out of memory.");
  } catch (...) {
    Rf_error("PredictProR native median distance failed.");
  }
  return R_NilValue;
}
