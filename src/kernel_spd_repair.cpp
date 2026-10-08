#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

static double scalar_real(SEXP x, double fallback) {
  if (Rf_length(x) < 1) {
    return fallback;
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, REALSXP));
  double value = REAL(coerced)[0];
  UNPROTECT(1);
  return R_FINITE(value) ? value : fallback;
}

static bool scalar_logical(SEXP x, bool fallback) {
  if (Rf_length(x) < 1) {
    return fallback;
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, LGLSXP));
  int value = LOGICAL(coerced)[0];
  UNPROTECT(1);
  if (value == NA_LOGICAL) {
    return fallback;
  }
  return value == TRUE;
}

static double median_positive_diagonal(const double* x, int n) {
  std::vector<double> diag;
  diag.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double value = x[i + i * n];
    if (R_FINITE(value) && value > 0.0) {
      diag.push_back(value);
    }
  }
  if (diag.empty()) {
    return 1.0;
  }
  std::sort(diag.begin(), diag.end());
  std::size_t mid = diag.size() / 2;
  if (diag.size() % 2 == 0) {
    return 0.5 * (diag[mid - 1] + diag[mid]);
  }
  return diag[mid];
}

extern "C" SEXP predictpror_kernel_spd_repair(SEXP x_sexp,
                                              SEXP min_eigen_sexp,
                                              SEXP keep_diag_sexp) {
  if (!Rf_isReal(x_sexp) || !Rf_isMatrix(x_sexp)) {
    Rf_error("x must be a numeric matrix");
  }

  SEXP dim = Rf_getAttrib(x_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (n != p) {
    Rf_error("x must be square");
  }
  if (n < 1) {
    Rf_error("x must have at least one row");
  }

  double min_eigen = scalar_real(min_eigen_sexp, 1e-8);
  if (!R_FINITE(min_eigen) || min_eigen <= 0.0) {
    min_eigen = 1e-8;
  }
  bool keep_diag = scalar_logical(keep_diag_sexp, true);

  const double* x = REAL(x_sexp);
  double diag_target = median_positive_diagonal(x, n);

  std::vector<double> original_diag(static_cast<std::size_t>(n));
  std::vector<double> a(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    double di = x[i + i * n];
    if (!R_FINITE(di) || di <= 0.0) {
      di = diag_target;
    }
    original_diag[static_cast<std::size_t>(i)] = di;
  }

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) {
      double xij = x[i + j * n];
      double xji = x[j + i * n];
      if (!R_FINITE(xij)) {
        xij = (i == j) ? original_diag[static_cast<std::size_t>(i)] : 0.0;
      }
      if (!R_FINITE(xji)) {
        xji = (i == j) ? original_diag[static_cast<std::size_t>(i)] : 0.0;
      }
      double value = (i == j) ? original_diag[static_cast<std::size_t>(i)] : 0.5 * (xij + xji);
      a[i + j * n] = value;
      a[j + i * n] = value;
    }
  }

  std::vector<double> eigenvalues(static_cast<std::size_t>(n));
  int lda = n;
  int lwork = -1;
  int info = 0;
  double work_query = 0.0;
  char jobz = 'V';
  char uplo = 'U';

  F77_CALL(dsyev)(&jobz, &uplo, &n, a.data(), &lda, eigenvalues.data(),
                  &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    Rf_error("LAPACK dsyev workspace query failed with code %d", info);
  }

  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(static_cast<std::size_t>(lwork));
  F77_CALL(dsyev)(&jobz, &uplo, &n, a.data(), &lda, eigenvalues.data(),
                  work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    Rf_error("LAPACK dsyev failed with code %d", info);
  }

  double min_before = eigenvalues[0];
  double max_before = eigenvalues[0];
  int adjusted_count = 0;
  for (int k = 0; k < n; ++k) {
    min_before = std::min(min_before, eigenvalues[static_cast<std::size_t>(k)]);
    max_before = std::max(max_before, eigenvalues[static_cast<std::size_t>(k)]);
    if (eigenvalues[static_cast<std::size_t>(k)] < min_eigen) {
      eigenvalues[static_cast<std::size_t>(k)] = min_eigen;
      ++adjusted_count;
    }
  }

  SEXP out_matrix = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_matrix);
  std::fill(out, out + static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) {
      double value = 0.0;
      for (int k = 0; k < n; ++k) {
        double qik = a[i + k * n];
        double qjk = a[j + k * n];
        value += qik * eigenvalues[static_cast<std::size_t>(k)] * qjk;
      }
      out[i + j * n] = value;
      out[j + i * n] = value;
    }
  }

  bool diag_rescaled = false;
  if (keep_diag) {
    std::vector<double> scale(static_cast<std::size_t>(n), 1.0);
    for (int i = 0; i < n; ++i) {
      double repaired_diag = out[i + i * n];
      if (R_FINITE(repaired_diag) && repaired_diag > 0.0) {
        scale[static_cast<std::size_t>(i)] =
          std::sqrt(original_diag[static_cast<std::size_t>(i)] / repaired_diag);
      }
    }
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i <= j; ++i) {
        double value = out[i + j * n] *
          scale[static_cast<std::size_t>(i)] *
          scale[static_cast<std::size_t>(j)];
        out[i + j * n] = value;
        out[j + i * n] = value;
      }
    }
    diag_rescaled = true;
  }

  double min_after = min_eigen;
  if (adjusted_count == 0) {
    min_after = min_before;
  }

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 8));
  SET_STRING_ELT(names, 0, Rf_mkChar("matrix"));
  SET_STRING_ELT(names, 1, Rf_mkChar("method"));
  SET_STRING_ELT(names, 2, Rf_mkChar("min_eigen_before"));
  SET_STRING_ELT(names, 3, Rf_mkChar("max_eigen_before"));
  SET_STRING_ELT(names, 4, Rf_mkChar("min_eigen_after"));
  SET_STRING_ELT(names, 5, Rf_mkChar("floor_value"));
  SET_STRING_ELT(names, 6, Rf_mkChar("adjusted_count"));
  SET_STRING_ELT(names, 7, Rf_mkChar("diag_rescaled"));

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(result, 0, out_matrix);
  SET_VECTOR_ELT(result, 1, Rf_mkString("nearPD_cpp"));
  SET_VECTOR_ELT(result, 2, Rf_ScalarReal(min_before));
  SET_VECTOR_ELT(result, 3, Rf_ScalarReal(max_before));
  SET_VECTOR_ELT(result, 4, Rf_ScalarReal(min_after));
  SET_VECTOR_ELT(result, 5, Rf_ScalarReal(min_eigen));
  SET_VECTOR_ELT(result, 6, Rf_ScalarInteger(adjusted_count));
  SET_VECTOR_ELT(result, 7, Rf_ScalarLogical(diag_rescaled ? TRUE : FALSE));
  Rf_setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(3);
  return result;
}
