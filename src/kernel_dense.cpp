#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

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

static double row_dot(const double* x, int n, int p, int row_a, int row_b) {
  double value = 0.0;
  for (int col = 0; col < p; ++col) {
    value += x[row_a + col * n] * x[row_b + col * n];
  }
  return value;
}

static double row_squared_distance(const double* x, int n, int p,
                                   int row_a, int row_b) {
  double value = 0.0;
  for (int col = 0; col < p; ++col) {
    double diff = x[row_a + col * n] - x[row_b + col * n];
    value += diff * diff;
  }
  return value;
}

static double row_l1_distance(const double* x, int n, int p,
                              int row_a, int row_b) {
  double value = 0.0;
  for (int col = 0; col < p; ++col) {
    value += std::fabs(x[row_a + col * n] - x[row_b + col * n]);
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

static double median_distance_squared(const double* x, int n, int p) {
  std::vector<double> distances;
  distances.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));

  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      double value = row_squared_distance(x, n, p, row, col);
      if (!R_FINITE(value)) {
        Rf_error("x must contain only finite values");
      }
      distances.push_back(value);
    }
  }

  std::sort(distances.begin(), distances.end());
  std::size_t count = distances.size();
  std::size_t mid = count / 2;
  if (count % 2 == 0) {
    return 0.5 * (distances[mid - 1] + distances[mid]);
  }
  return distances[mid];
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

  double median_d2 = 0.0;
  bool needs_gaussian = method == KERNEL_GAUSSIAN || method == KERNEL_COMPOSITE;
  if (needs_gaussian) {
    median_d2 = median_distance_squared(x, n, p);
    if (!R_FINITE(median_d2) || median_d2 <= 0.0) {
      Rf_error("median squared distance must be positive and finite");
    }
  }

  SEXP outSEXP = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(outSEXP);

  for (int col = 0; col < n; ++col) {
    if ((col & 63) == 0) {
      R_CheckUserInterrupt();
    }

    for (int row = 0; row <= col; ++row) {
      double linear = row_dot(x, n, p, row, col) / static_cast<double>(p);
      double value = linear;

      if (method == KERNEL_GAUSSIAN || method == KERNEL_COMPOSITE) {
        double d2 = row_squared_distance(x, n, p, row, col);
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
        double d2 = row_squared_distance(x, n, p, row, col);
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
        double d1 = row_l1_distance(x, n, p, row, col);
        value = std::exp(-d1 / length_scale);
      }

      out[row + col * n] = value;
      out[col + row * n] = value;
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
