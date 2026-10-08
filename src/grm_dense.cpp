#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <cmath>
#include <cstring>
#include <exception>
#include <new>
#include <vector>

static const char* grm_method_string(SEXP x) {
  if (TYPEOF(x) != STRSXP || Rf_length(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("method must be a string scalar");
  }
  return CHAR(STRING_ELT(x, 0));
}

static bool is_hard_call_genotype(double value) {
  return value == 0.0 || value == 1.0 || value == 2.0;
}

static double dominance_vitezica_value(double genotype, double p, double q) {
  if (genotype == 0.0) {
    return -2.0 * p * p;
  }
  if (genotype == 1.0) {
    return 2.0 * p * q;
  }
  return -2.0 * q * q;
}

static double dominance_su_value(double genotype, double p, double q) {
  const double expected_heterozygosity = 2.0 * p * q;
  if (genotype == 1.0) {
    return 1.0 - expected_heterozygosity;
  }
  return -expected_heterozygosity;
}

static SEXP predictpror_grm_dense_impl(SEXP genoSEXP,
                                       SEXP methodSEXP,
                                       SEXP weightsSEXP,
                                       SEXP paramsSEXP) {
  (void)paramsSEXP;

  if (!Rf_isReal(genoSEXP) || !Rf_isMatrix(genoSEXP)) {
    Rf_error("geno must be a numeric matrix");
  }

  SEXP dimSEXP = Rf_getAttrib(genoSEXP, R_DimSymbol);
  if (TYPEOF(dimSEXP) != INTSXP || Rf_length(dimSEXP) != 2) {
    Rf_error("geno must have matrix dimensions");
  }

  const int n = INTEGER(dimSEXP)[0];
  const int p = INTEGER(dimSEXP)[1];
  if (n <= 0 || p <= 0) {
    Rf_error("geno dimensions must be positive");
  }

  const char* method = grm_method_string(methodSEXP);
  const bool is_vanraden = std::strcmp(method, "VanRaden") == 0;
  const bool is_weighted = std::strcmp(method, "Weighted_VanRaden") == 0;
  const bool is_yang = std::strcmp(method, "Yang") == 0;
  const bool is_epistasis = std::strcmp(method, "Epistasis") == 0;
  const bool is_dominance_vitezica =
      std::strcmp(method, "Dominance") == 0 ||
      std::strcmp(method, "Dominance_Vitezica") == 0;
  const bool is_dominance_su =
      std::strcmp(method, "Dominance_Su") == 0 ||
      std::strcmp(method, "Dominance_Heterozygosity") == 0;
  const bool is_dominance = is_dominance_vitezica || is_dominance_su;
  if (!is_vanraden && !is_weighted && !is_yang && !is_epistasis && !is_dominance) {
    Rf_error("unsupported GRM method");
  }

  const double* geno = REAL(genoSEXP);
  const double* weights = NULL;
  if (weightsSEXP != R_NilValue) {
    if (!Rf_isReal(weightsSEXP) || Rf_length(weightsSEXP) < p) {
      Rf_error("weights must be a numeric vector with length at least ncol(geno)");
    }
    weights = REAL(weightsSEXP);
    for (int k = 0; k < p; ++k) {
      if (!R_FINITE(weights[k])) {
        Rf_error("weights must be finite");
      }
    }
  }
  if (is_weighted && weights == NULL) {
    Rf_error("weights are required for Weighted_VanRaden");
  }

  std::vector<double> freq(p);
  double denom = 0.0;
  for (int k = 0; k < p; ++k) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double value = geno[i + n * k];
      if (!R_FINITE(value)) {
        Rf_error("geno must contain only finite values");
      }
      if (!is_hard_call_genotype(value)) {
        Rf_error("geno must be coded as 0, 1, 2");
      }
      sum += value;
    }
    const double f = sum / (2.0 * static_cast<double>(n));
    freq[k] = f;
    denom += 2.0 * f * (1.0 - f);
  }
  if (!R_FINITE(denom) || denom <= 0.0) {
    Rf_error("GRM denominator must be positive and finite");
  }

  SEXP outSEXP = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(outSEXP);

  if (is_yang) {
    std::vector<double> inv_var(p);
    for (int k = 0; k < p; ++k) {
      const double var = 2.0 * freq[k] * (1.0 - freq[k]);
      if (!R_FINITE(var) || var <= 0.0) {
        Rf_error("Yang method requires markers with positive finite variance");
      }
      inv_var[k] = 1.0 / var;
    }

    const double inv_marker_count = 1.0 / static_cast<double>(p);
    for (int i = 0; i < n; ++i) {
      if ((i & 63) == 0) {
        R_CheckUserInterrupt();
      }

      double diag_sum = 0.0;
      for (int k = 0; k < p; ++k) {
        const double g = geno[i + n * k];
        const double f = freq[k];
        diag_sum += (g * g - g * (1.0 + 2.0 * f) + 2.0 * f * f) * inv_var[k];
      }
      out[i + n * i] = 1.0 + inv_marker_count * diag_sum;

      for (int j = i + 1; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < p; ++k) {
          const double centered_i = geno[i + n * k] - 2.0 * freq[k];
          const double centered_j = geno[j + n * k] - 2.0 * freq[k];
          sum += centered_i * centered_j * inv_var[k];
        }
        const double value = inv_marker_count * sum;
        out[i + n * j] = value;
        out[j + n * i] = value;
      }
    }
  } else if (is_dominance) {
    double dominance_denom = 0.0;
    for (int k = 0; k < p; ++k) {
      const double f = freq[k];
      const double q = 1.0 - f;
      const double expected_heterozygosity = 2.0 * f * q;
      dominance_denom += is_dominance_vitezica
        ? expected_heterozygosity * expected_heterozygosity
        : expected_heterozygosity * (1.0 - expected_heterozygosity);
    }
    if (!R_FINITE(dominance_denom) || dominance_denom <= 0.0) {
      Rf_error("dominance GRM denominator must be positive and finite");
    }

    for (int i = 0; i < n; ++i) {
      if ((i & 63) == 0) {
        R_CheckUserInterrupt();
      }

      for (int j = i; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < p; ++k) {
          const double f = freq[k];
          const double q = 1.0 - f;
          const double dominance_i = is_dominance_vitezica
            ? dominance_vitezica_value(geno[i + n * k], f, q)
            : dominance_su_value(geno[i + n * k], f, q);
          const double dominance_j = is_dominance_vitezica
            ? dominance_vitezica_value(geno[j + n * k], f, q)
            : dominance_su_value(geno[j + n * k], f, q);
          sum += dominance_i * dominance_j;
        }
        const double value = sum / dominance_denom;
        out[i + n * j] = value;
        out[j + n * i] = value;
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      if ((i & 63) == 0) {
        R_CheckUserInterrupt();
      }

      for (int j = i; j < n; ++j) {
        double sum = 0.0;
        for (int k = 0; k < p; ++k) {
          const double centered_i = geno[i + n * k] - 2.0 * freq[k];
          const double centered_j = geno[j + n * k] - 2.0 * freq[k];
          const double weight = is_weighted ? weights[k] : 1.0;
          sum += centered_i * centered_j * weight;
        }
        double value = sum / denom;
        if (is_epistasis) {
          value *= value;
        }
        out[i + n * j] = value;
        out[j + n * i] = value;
      }
    }
  }

  UNPROTECT(1);
  return outSEXP;
}

extern "C" SEXP predictpror_grm_dense(SEXP genoSEXP,
                                      SEXP methodSEXP,
                                      SEXP weightsSEXP,
                                      SEXP paramsSEXP) {
  try {
    return predictpror_grm_dense_impl(genoSEXP, methodSEXP, weightsSEXP, paramsSEXP);
  } catch (const std::bad_alloc&) {
    Rf_error("PredictProR native dense GRM backend ran out of memory.");
  } catch (const std::exception& ex) {
    Rf_error("PredictProR native dense GRM backend failed: %s", ex.what());
  } catch (...) {
    Rf_error("PredictProR native dense GRM backend failed.");
  }
  return R_NilValue;
}
