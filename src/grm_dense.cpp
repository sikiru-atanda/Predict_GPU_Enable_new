#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif

#include <cmath>
#include <cstring>
#include <exception>
#include <new>

// Dense genomic relationship matrices. Each method is a product G = c Z Z'
// of a marker covariate matrix Z (n lines x p markers), built here in one pass
// over the dosage, and formed with one BLAS dsyrk call: a triple loop over
// pairs of lines took 143 s for 1,000 x 10,000 where dsyrk takes 5-8 s. The
// additive methods take any ploidy; Yang and the dominance codings are diploid.
// Every step is the R backend's, in its order (sums in long double as R's
// sum() and colSums() take them, dsyrk with alpha 1, then the division), so
// the two backends give the same numbers.

enum GrmMethod {
  GRM_VANRADEN,
  GRM_WEIGHTED,
  GRM_EPISTASIS,
  GRM_YANG,
  GRM_DOMINANCE_VITEZICA,
  GRM_DOMINANCE_SU
};

static GrmMethod parse_grm_method(SEXP x) {
  if (TYPEOF(x) != STRSXP || Rf_length(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("method must be a string scalar");
  }
  const char* method = CHAR(STRING_ELT(x, 0));
  if (std::strcmp(method, "VanRaden") == 0) return GRM_VANRADEN;
  if (std::strcmp(method, "Weighted_VanRaden") == 0) return GRM_WEIGHTED;
  if (std::strcmp(method, "Epistasis") == 0) return GRM_EPISTASIS;
  if (std::strcmp(method, "Yang") == 0) return GRM_YANG;
  if (std::strcmp(method, "Dominance") == 0 || std::strcmp(method, "Dominance_Vitezica") == 0) {
    return GRM_DOMINANCE_VITEZICA;
  }
  if (std::strcmp(method, "Dominance_Su") == 0 || std::strcmp(method, "Dominance_Heterozygosity") == 0) {
    return GRM_DOMINANCE_SU;
  }
  Rf_error("unsupported GRM method");
  return GRM_VANRADEN;
}

// ploidy from params$ploidy (default 2)
static int grm_ploidy(SEXP params) {
  if (TYPEOF(params) != VECSXP) return 2;
  SEXP names = Rf_getAttrib(params, R_NamesSymbol);
  for (R_xlen_t i = 0; i < Rf_xlength(params); ++i) {
    if (names != R_NilValue && std::strcmp(CHAR(STRING_ELT(names, i)), "ploidy") == 0) {
      SEXP v = VECTOR_ELT(params, i);
      if (!Rf_isNumeric(v) || Rf_length(v) != 1) Rf_error("ploidy must be one positive integer");
      const int ploidy = Rf_asInteger(v);
      if (ploidy == NA_INTEGER || ploidy < 1) Rf_error("ploidy must be one positive integer");
      return ploidy;
    }
  }
  return 2;
}

static SEXP predictpror_grm_dense_impl(SEXP genoSEXP, SEXP methodSEXP, SEXP weightsSEXP, SEXP paramsSEXP) {
  if (!Rf_isReal(genoSEXP) || !Rf_isMatrix(genoSEXP)) {
    Rf_error("geno must be a numeric matrix");
  }
  SEXP dimSEXP = Rf_getAttrib(genoSEXP, R_DimSymbol);
  const int n = INTEGER(dimSEXP)[0];
  const int p = INTEGER(dimSEXP)[1];
  if (n <= 0 || p <= 0) Rf_error("geno dimensions must be positive");
  const R_xlen_t nn = static_cast<R_xlen_t>(n);

  const GrmMethod method = parse_grm_method(methodSEXP);
  const int ploidy = grm_ploidy(paramsSEXP);
  const bool diploid_only = method == GRM_YANG || method == GRM_DOMINANCE_VITEZICA || method == GRM_DOMINANCE_SU;
  if (diploid_only && ploidy != 2) {
    Rf_error("this GRM method is diploid-only");
  }

  const double* geno = REAL(genoSEXP);
  const double* weights = NULL;
  if (method == GRM_WEIGHTED) {
    if (weightsSEXP == R_NilValue || !Rf_isReal(weightsSEXP) || Rf_xlength(weightsSEXP) < p) {
      Rf_error("weights must be a numeric vector with one value per marker");
    }
    weights = REAL(weightsSEXP);
    for (int k = 0; k < p; ++k) {
      if (!R_FINITE(weights[k]) || weights[k] < 0.0) Rf_error("weights must be finite and non-negative");
    }
  }

  // R-managed buffers: freed on return, on an R error and on an interrupt
  // (a std::vector's destructor never runs when Rf_error longjmps out)
  const void* vmax = vmaxget();
  double* mean = reinterpret_cast<double*>(R_alloc(p, sizeof(double)));
  double* freq = reinterpret_cast<double*>(R_alloc(p, sizeof(double)));

  // column means as R's colMeans() takes them (long double accumulation)
  for (int k = 0; k < p; ++k) {
    const double* col = geno + nn * k;
    long double sum = 0.0L;
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(col[i])) Rf_error("geno must contain only finite values");
      sum += col[i];
    }
    mean[k] = static_cast<double>(sum / n);
  }

  // the denominator, as the R backend sums it
  long double denom_sum = 0.0L;
  for (int k = 0; k < p; ++k) {
    const double f = mean[k] / ploidy;
    freq[k] = f;
    const double q = 1.0 - f;
    if (method == GRM_YANG && !(R_FINITE(2.0 * f * q) && 2.0 * f * q > 0.0)) {
      Rf_error("Yang method requires markers with positive finite variance");
    }
    double term;
    switch (method) {
      case GRM_DOMINANCE_VITEZICA: term = (2.0 * f * q) * (2.0 * f * q); break;
      case GRM_DOMINANCE_SU: term = 2.0 * f * q * (1.0 - 2.0 * f * q); break;
      default: term = ploidy * f * q; break;
    }
    denom_sum += term;
  }
  const double denom = static_cast<double>(denom_sum);
  if (method != GRM_YANG && (!R_FINITE(denom) || denom <= 0.0)) {
    Rf_error("GRM denominator must be positive and finite");
  }

  // the covariates, column by column
  SEXP outSEXP = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(outSEXP);
  long double* yang_diag = NULL;
  if (method == GRM_YANG) {
    yang_diag = reinterpret_cast<long double*>(R_alloc(n, sizeof(long double)));
    for (int i = 0; i < n; ++i) yang_diag[i] = 0.0L;
  }
  const void* vmax_z = vmaxget();
  double* z = reinterpret_cast<double*>(R_alloc(nn * p, sizeof(double)));
  for (int k = 0; k < p; ++k) {
    const double* col = geno + nn * k;
    double* zc = z + nn * k;
    const double f = freq[k];
    const double q = 1.0 - f;
    switch (method) {
      case GRM_VANRADEN:
      case GRM_EPISTASIS:
        for (int i = 0; i < n; ++i) zc[i] = col[i] - mean[k];
        break;
      case GRM_WEIGHTED: {
        const double s = std::sqrt(weights[k]);
        for (int i = 0; i < n; ++i) zc[i] = (col[i] - mean[k]) * s;
        break;
      }
      case GRM_YANG: {
        const double inv_var = 1.0 / (2.0 * f * q);
        const double s = std::sqrt(inv_var);
        for (int i = 0; i < n; ++i) {
          const double g = col[i];
          zc[i] = (g - mean[k]) * s;
          const double term = (g * g - g * (1.0 + 2.0 * f) + 2.0 * f * f) * inv_var;
          yang_diag[i] += term;
        }
        break;
      }
      case GRM_DOMINANCE_VITEZICA: {
        const double v0 = -2.0 * f * f, v1 = 2.0 * f * q, v2 = -2.0 * q * q;
        for (int i = 0; i < n; ++i) {
          const double g = col[i];
          zc[i] = g == 0.0 ? v0 : (g == 1.0 ? v1 : (g == 2.0 ? v2 : 0.0));
        }
        break;
      }
      case GRM_DOMINANCE_SU: {
        const double pq2 = 2.0 * f * q;
        for (int i = 0; i < n; ++i) {
          const double g = col[i];
          zc[i] = (g == 0.0 || g == 2.0) ? -pq2 : (g == 1.0 ? 1.0 - pq2 : 0.0);
        }
        break;
      }
    }
    if ((k & 1023) == 0) R_CheckUserInterrupt();
  }

  const double alpha = 1.0, beta = 0.0;
  F77_CALL(dsyrk)("U", "N", &n, &p, &alpha, z, &n, &beta, out, &n FCONE FCONE);
  vmaxset(vmax_z);  // the covariates are no longer needed

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) out[j + nn * i] = out[i + nn * j];
  }
  const R_xlen_t total = nn * nn;
  if (method == GRM_YANG) {
    const double inv_p = 1.0 / p;
    for (R_xlen_t t = 0; t < total; ++t) out[t] = inv_p * out[t];
    for (int i = 0; i < n; ++i) out[i + nn * i] = 1.0 + inv_p * static_cast<double>(yang_diag[i]);
  } else {
    for (R_xlen_t t = 0; t < total; ++t) out[t] = out[t] / denom;
  }
  if (method == GRM_EPISTASIS) {
    for (R_xlen_t t = 0; t < total; ++t) out[t] *= out[t];
  }

  vmaxset(vmax);
  UNPROTECT(1);
  return outSEXP;
}

extern "C" SEXP predictpror_grm_dense(SEXP genoSEXP, SEXP methodSEXP, SEXP weightsSEXP, SEXP paramsSEXP) {
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
