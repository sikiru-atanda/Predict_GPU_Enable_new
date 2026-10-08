#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <vector>

static int scalar_int(SEXP x, int fallback, int min_value) {
  if (Rf_length(x) < 1) {
    return fallback;
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, INTSXP));
  int value = INTEGER(coerced)[0];
  UNPROTECT(1);
  if (value == NA_INTEGER || value < min_value) {
    return fallback;
  }
  return value;
}

static double scalar_real_value(SEXP x, double fallback) {
  if (Rf_length(x) < 1) {
    return fallback;
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, REALSXP));
  double value = REAL(coerced)[0];
  UNPROTECT(1);
  return R_FINITE(value) || std::isinf(value) ? value : fallback;
}

static std::vector<int> normalized_row_index(SEXP row_index_sexp, int n) {
  std::vector<int> rows;
  if (Rf_isNull(row_index_sexp)) {
    rows.reserve(static_cast<std::size_t>(n));
    for (int i = 1; i <= n; ++i) {
      rows.push_back(i);
    }
    return rows;
  }

  SEXP coerced = PROTECT(Rf_coerceVector(row_index_sexp, INTSXP));
  int m = Rf_length(coerced);
  rows.reserve(static_cast<std::size_t>(m));
  for (int i = 0; i < m; ++i) {
    int value = INTEGER(coerced)[i];
    if (value != NA_INTEGER && value >= 1 && value <= n) {
      rows.push_back(value);
    }
  }
  UNPROTECT(1);

  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

extern "C" SEXP predictpror_kernel_duplicate_pairs(SEXP x_sexp,
                                                   SEXP threshold_sexp,
                                                   SEXP diag_epsilon_sexp,
                                                   SEXP block_size_sexp,
                                                   SEXP max_pairs_sexp,
                                                   SEXP row_index_sexp) {
  if (!Rf_isReal(x_sexp) || !Rf_isMatrix(x_sexp)) {
    Rf_error("x must be a numeric matrix");
  }

  SEXP dim = Rf_getAttrib(x_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (n != p) {
    Rf_error("x must be square");
  }

  double threshold = scalar_real_value(threshold_sexp, 0.95);
  double diag_epsilon = scalar_real_value(diag_epsilon_sexp, 1e-15);
  if (!R_FINITE(diag_epsilon) || diag_epsilon < 0.0) {
    diag_epsilon = 1e-15;
  }
  int block_size = scalar_int(block_size_sexp, 1024, 1);
  double max_pairs_d = scalar_real_value(max_pairs_sexp, R_PosInf);
  bool finite_max_pairs = R_FINITE(max_pairs_d);
  int max_pairs = 0;
  if (finite_max_pairs) {
    max_pairs_d = std::max(0.0, std::min(max_pairs_d, static_cast<double>(INT_MAX)));
    max_pairs = static_cast<int>(max_pairs_d);
  }

  std::vector<int> row_index = normalized_row_index(row_index_sexp, n);
  if (n < 2 || row_index.size() < 2 || threshold >= 1.0 ||
      (finite_max_pairs && max_pairs <= 0)) {
    SEXP row = PROTECT(Rf_allocVector(INTSXP, 0));
    SEXP col = PROTECT(Rf_allocVector(INTSXP, 0));
    SEXP corr = PROTECT(Rf_allocVector(REALSXP, 0));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("Row"));
    SET_STRING_ELT(names, 1, Rf_mkChar("Col"));
    SET_STRING_ELT(names, 2, Rf_mkChar("Corr"));
    SET_VECTOR_ELT(out, 0, row);
    SET_VECTOR_ELT(out, 1, col);
    SET_VECTOR_ELT(out, 2, corr);
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(5);
    return out;
  }

  const double* x = REAL(x_sexp);
  std::vector<double> inv_diag_scale(static_cast<std::size_t>(n), NA_REAL);
  std::vector<unsigned char> valid(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    double diag_value = x[i + i * n];
    double positive = (R_FINITE(diag_value) && diag_value > 0.0) ? diag_value : 0.0;
    double scale = std::sqrt(positive + diag_epsilon);
    if (R_FINITE(scale) && scale > 0.0) {
      inv_diag_scale[static_cast<std::size_t>(i)] = 1.0 / scale;
      valid[static_cast<std::size_t>(i)] = 1;
    }
  }

  std::vector<int> hit_row;
  std::vector<int> hit_col;
  std::vector<double> hit_corr;
  if (finite_max_pairs) {
    hit_row.reserve(static_cast<std::size_t>(max_pairs));
    hit_col.reserve(static_cast<std::size_t>(max_pairs));
    hit_corr.reserve(static_cast<std::size_t>(max_pairs));
  }

  int selected_n = static_cast<int>(row_index.size());
  bool stop = false;
  for (int row_start = 0; row_start < selected_n - 1 && !stop; row_start += block_size) {
    int row_end = std::min(selected_n, row_start + block_size);
    for (int col_start = row_start; col_start < selected_n && !stop; col_start += block_size) {
      int col_end = std::min(selected_n, col_start + block_size);
      for (int col_pos = col_start; col_pos < col_end && !stop; ++col_pos) {
        int col_1 = row_index[static_cast<std::size_t>(col_pos)];
        int col_0 = col_1 - 1;
        if (!valid[static_cast<std::size_t>(col_0)]) {
          continue;
        }
        double col_scale = inv_diag_scale[static_cast<std::size_t>(col_0)];
        int row_limit = (row_start == col_start) ? std::min(row_end, col_pos) : row_end;
        for (int row_pos = row_start; row_pos < row_limit; ++row_pos) {
          int row_1 = row_index[static_cast<std::size_t>(row_pos)];
          int row_0 = row_1 - 1;
          if (!valid[static_cast<std::size_t>(row_0)]) {
            continue;
          }
          double value = x[row_0 + col_0 * n] *
            inv_diag_scale[static_cast<std::size_t>(row_0)] *
            col_scale;
          if (R_FINITE(value) && value > threshold) {
            hit_row.push_back(row_1);
            hit_col.push_back(col_1);
            hit_corr.push_back(value);
            if (finite_max_pairs && static_cast<int>(hit_row.size()) >= max_pairs) {
              stop = true;
              break;
            }
          }
        }
      }
      R_CheckUserInterrupt();
    }
  }

  R_xlen_t out_n = static_cast<R_xlen_t>(hit_row.size());
  SEXP row = PROTECT(Rf_allocVector(INTSXP, out_n));
  SEXP col = PROTECT(Rf_allocVector(INTSXP, out_n));
  SEXP corr = PROTECT(Rf_allocVector(REALSXP, out_n));
  for (R_xlen_t i = 0; i < out_n; ++i) {
    INTEGER(row)[i] = hit_row[static_cast<std::size_t>(i)];
    INTEGER(col)[i] = hit_col[static_cast<std::size_t>(i)];
    REAL(corr)[i] = hit_corr[static_cast<std::size_t>(i)];
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("Row"));
  SET_STRING_ELT(names, 1, Rf_mkChar("Col"));
  SET_STRING_ELT(names, 2, Rf_mkChar("Corr"));
  SET_VECTOR_ELT(out, 0, row);
  SET_VECTOR_ELT(out, 1, col);
  SET_VECTOR_ELT(out, 2, corr);
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(5);
  return out;
}
