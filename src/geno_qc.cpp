#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>
#include <vector>

static void validate_numeric_matrix(SEXP x, const char *name) {
  if (!Rf_isMatrix(x) || (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP)) {
    Rf_error("%s must be a numeric or integer matrix", name);
  }
}

static double numeric_matrix_value(SEXP x, R_xlen_t idx) {
  if (TYPEOF(x) == REALSXP) {
    return REAL(x)[idx];
  }
  int value = INTEGER(x)[idx];
  return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
}

static bool is_missing_numeric(double value) {
  return ISNAN(value);
}

static int scalar_int_value(SEXP x, int fallback, int min_value) {
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

static double round_like_r(double value) {
  if (!R_FINITE(value)) {
    return value;
  }
  return std::nearbyint(value);
}

static void set_matrix_dimnames(SEXP out, SEXP source) {
  SEXP dimnames = Rf_getAttrib(source, R_DimNamesSymbol);
  if (!Rf_isNull(dimnames)) {
    Rf_setAttrib(out, R_DimNamesSymbol, dimnames);
  }
}

static SEXP matrix_dimnames_element(SEXP source, int index) {
  SEXP dimnames = Rf_getAttrib(source, R_DimNamesSymbol);
  if (Rf_isNull(dimnames) || TYPEOF(dimnames) != VECSXP || Rf_length(dimnames) <= index) {
    return R_NilValue;
  }
  return VECTOR_ELT(dimnames, index);
}

extern "C" SEXP predictpror_geno_qc_metrics(SEXP geno_sexp, SEXP ploidy_sexp) {
  validate_numeric_matrix(geno_sexp, "geno_data");
  const int ploidy = scalar_int_value(ploidy_sexp, 2, 1);

  SEXP dim = Rf_getAttrib(geno_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(out_names, 0, Rf_mkChar("monomorphic"));
  SET_STRING_ELT(out_names, 1, Rf_mkChar("marker_missing_rate"));
  SET_STRING_ELT(out_names, 2, Rf_mkChar("individual_missing_rate"));
  SET_STRING_ELT(out_names, 3, Rf_mkChar("maf"));
  SET_STRING_ELT(out_names, 4, Rf_mkChar("heterozygosity"));
  Rf_setAttrib(out, R_NamesSymbol, out_names);

  SEXP marker_missing = PROTECT(Rf_allocVector(REALSXP, m));
  SEXP individual_missing = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP maf = PROTECT(Rf_allocVector(REALSXP, m));
  SEXP heterozygosity = PROTECT(Rf_allocVector(REALSXP, m));
  SEXP row_names = matrix_dimnames_element(geno_sexp, 0);
  SEXP col_names = matrix_dimnames_element(geno_sexp, 1);
  if (!Rf_isNull(col_names)) {
    Rf_setAttrib(marker_missing, R_NamesSymbol, col_names);
    Rf_setAttrib(maf, R_NamesSymbol, col_names);
    Rf_setAttrib(heterozygosity, R_NamesSymbol, col_names);
  }
  if (!Rf_isNull(row_names)) {
    Rf_setAttrib(individual_missing, R_NamesSymbol, row_names);
  }

  std::vector<int> monomorphic;
  monomorphic.reserve(static_cast<std::size_t>(m));

  for (int row = 0; row < n; ++row) {
    REAL(individual_missing)[row] = 0.0;
  }

  for (int col = 0; col < m; ++col) {
    R_xlen_t offset = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n);
    int missing = 0;
    int observed = 0;
    int heterozygous = 0;
    bool have_first = false;
    bool has_multiple_observed_values = false;
    double first_value = NA_REAL;
    double sum = 0.0;

    for (int row = 0; row < n; ++row) {
      double value = numeric_matrix_value(geno_sexp, offset + row);
      if (is_missing_numeric(value)) {
        ++missing;
        REAL(individual_missing)[row] += 1.0;
        continue;
      }
      ++observed;
      sum += value;
      if (value > 0.0 && value < static_cast<double>(ploidy)) {
        ++heterozygous;
      }
      if (!have_first) {
        first_value = value;
        have_first = true;
      } else if (value != first_value) {
        has_multiple_observed_values = true;
      }
    }

    REAL(marker_missing)[col] = n > 0 ? static_cast<double>(missing) / static_cast<double>(n) : NA_REAL;
    REAL(heterozygosity)[col] = n > 0 ? static_cast<double>(heterozygous) / static_cast<double>(n) : NA_REAL;
    if (observed > 0) {
      double phat = (sum / static_cast<double>(observed)) / static_cast<double>(ploidy);
      REAL(maf)[col] = phat < 0.5 ? phat : 1.0 - phat;
    } else {
      REAL(maf)[col] = NA_REAL;
    }
    if (!has_multiple_observed_values) {
      monomorphic.push_back(col + 1);
    }
  }

  for (int row = 0; row < n; ++row) {
    REAL(individual_missing)[row] = m > 0 ? REAL(individual_missing)[row] / static_cast<double>(m) : NA_REAL;
  }

  SEXP mono = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(monomorphic.size())));
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(monomorphic.size()); ++i) {
    INTEGER(mono)[i] = monomorphic[static_cast<std::size_t>(i)];
  }

  SET_VECTOR_ELT(out, 0, mono);
  SET_VECTOR_ELT(out, 1, marker_missing);
  SET_VECTOR_ELT(out, 2, individual_missing);
  SET_VECTOR_ELT(out, 3, maf);
  SET_VECTOR_ELT(out, 4, heterozygosity);

  UNPROTECT(7);
  return out;
}

extern "C" SEXP predictpror_geno_impute_summary(SEXP geno_sexp, SEXP method_sexp) {
  validate_numeric_matrix(geno_sexp, "geno_data");

  SEXP method_char = PROTECT(Rf_coerceVector(method_sexp, STRSXP));
  const char *method = CHAR(STRING_ELT(method_char, 0));
  bool use_mean = std::strcmp(method, "mean") == 0;
  bool use_median = std::strcmp(method, "median") == 0;
  if (!use_mean && !use_median) {
    Rf_error("method must be 'mean' or 'median'");
  }

  SEXP dim = Rf_getAttrib(geno_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, m));
  set_matrix_dimnames(out, geno_sexp);

  double *out_ptr = REAL(out);
  R_xlen_t total = static_cast<R_xlen_t>(n) * static_cast<R_xlen_t>(m);
  for (R_xlen_t idx = 0; idx < total; ++idx) {
    out_ptr[idx] = numeric_matrix_value(geno_sexp, idx);
  }

  std::vector<double> observed_values;
  observed_values.reserve(static_cast<std::size_t>(n));

  for (int col = 0; col < m; ++col) {
    R_xlen_t offset = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n);
    double replacement = NA_REAL;

    if (use_mean) {
      double sum = 0.0;
      int count = 0;
      for (int row = 0; row < n; ++row) {
        double value = out_ptr[offset + row];
        if (!is_missing_numeric(value)) {
          sum += value;
          ++count;
        }
      }
      replacement = count > 0 ? round_like_r(sum / static_cast<double>(count)) : R_NaN;
    } else {
      observed_values.clear();
      for (int row = 0; row < n; ++row) {
        double value = out_ptr[offset + row];
        if (!is_missing_numeric(value)) {
          observed_values.push_back(value);
        }
      }
      if (!observed_values.empty()) {
        std::sort(observed_values.begin(), observed_values.end());
        std::size_t mid = observed_values.size() / 2;
        double median = observed_values.size() % 2 == 1
          ? observed_values[mid]
          : (observed_values[mid - 1] + observed_values[mid]) / 2.0;
        replacement = round_like_r(median);
      }
    }

    for (int row = 0; row < n; ++row) {
      R_xlen_t idx = offset + row;
      if (is_missing_numeric(out_ptr[idx])) {
        out_ptr[idx] = replacement;
      }
    }
    if ((col & 127) == 0) {
      R_CheckUserInterrupt();
    }
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP predictpror_geno_impute_knn(SEXP geno_sexp,
                                             SEXP k_sexp,
                                             SEXP max_features_sexp,
                                             SEXP max_donors_sexp) {
  validate_numeric_matrix(geno_sexp, "geno_data");

  SEXP dim = Rf_getAttrib(geno_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];
  int k = scalar_int_value(k_sexp, 5, 1);
  int max_features = scalar_int_value(max_features_sexp, 2048, 1);
  int max_donors = scalar_int_value(max_donors_sexp, 512, 1);

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, m));
  set_matrix_dimnames(out, geno_sexp);
  double *out_ptr = REAL(out);

  R_xlen_t total = static_cast<R_xlen_t>(n) * static_cast<R_xlen_t>(m);
  std::vector<double> source(static_cast<std::size_t>(total));
  bool has_missing = false;
  for (R_xlen_t idx = 0; idx < total; ++idx) {
    double value = numeric_matrix_value(geno_sexp, idx);
    source[static_cast<std::size_t>(idx)] = value;
    out_ptr[idx] = value;
    if (is_missing_numeric(value)) {
      has_missing = true;
    }
  }

  if (!has_missing || n <= 1 || m < 1) {
    UNPROTECT(1);
    return out;
  }

  std::vector<double> column_fallback(static_cast<std::size_t>(m), NA_REAL);
  for (int col = 0; col < m; ++col) {
    R_xlen_t offset = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n);
    double sum = 0.0;
    int count = 0;
    for (int row = 0; row < n; ++row) {
      double value = source[static_cast<std::size_t>(offset + row)];
      if (!is_missing_numeric(value)) {
        sum += value;
        ++count;
      }
    }
    if (count > 0) {
      column_fallback[static_cast<std::size_t>(col)] = sum / static_cast<double>(count);
    }
  }

  std::vector<double> row_mean(static_cast<std::size_t>(n), 0.0);
  std::vector<int> row_order(static_cast<std::size_t>(n));
  std::vector<int> row_rank(static_cast<std::size_t>(n), 0);
  for (int row = 0; row < n; ++row) {
    double sum = 0.0;
    int count = 0;
    for (int col = 0; col < m; ++col) {
      R_xlen_t idx = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + row;
      double value = source[static_cast<std::size_t>(idx)];
      if (!is_missing_numeric(value)) {
        sum += value;
        ++count;
      }
    }
    row_mean[static_cast<std::size_t>(row)] = count > 0 ? sum / static_cast<double>(count) : 0.0;
    row_order[static_cast<std::size_t>(row)] = row;
  }
  std::sort(row_order.begin(), row_order.end(),
    [&row_mean](int a, int b) {
      double mean_a = row_mean[static_cast<std::size_t>(a)];
      double mean_b = row_mean[static_cast<std::size_t>(b)];
      if (mean_a == mean_b) {
        return a < b;
      }
      return mean_a < mean_b;
    });
  for (int pos = 0; pos < n; ++pos) {
    row_rank[static_cast<std::size_t>(row_order[static_cast<std::size_t>(pos)])] = pos;
  }

  std::vector<int> observed_cols;
  std::vector<int> missing_cols;
  std::vector<int> distance_cols;
  std::vector<int> candidate_rows;
  std::vector<std::pair<double, int> > distances;
  observed_cols.reserve(static_cast<std::size_t>(m));
  missing_cols.reserve(static_cast<std::size_t>(m));
  distance_cols.reserve(static_cast<std::size_t>(std::min(m, max_features)));
  candidate_rows.reserve(static_cast<std::size_t>(std::min(std::max(n - 1, 0), max_donors)));
  distances.reserve(static_cast<std::size_t>(n > 1 ? n - 1 : 0));

  for (int row = 0; row < n; ++row) {
    observed_cols.clear();
    missing_cols.clear();
    for (int col = 0; col < m; ++col) {
      R_xlen_t idx = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + row;
      if (is_missing_numeric(source[static_cast<std::size_t>(idx)])) {
        missing_cols.push_back(col);
      } else {
        observed_cols.push_back(col);
      }
    }

    if (missing_cols.empty()) {
      continue;
    }
    if (observed_cols.empty()) {
      for (std::size_t mi = 0; mi < missing_cols.size(); ++mi) {
        int col = missing_cols[mi];
        out_ptr[static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + row] =
          column_fallback[static_cast<std::size_t>(col)];
      }
      continue;
    }

    distance_cols.clear();
    if (static_cast<int>(observed_cols.size()) <= max_features) {
      distance_cols.insert(distance_cols.end(), observed_cols.begin(), observed_cols.end());
    } else {
      double step = static_cast<double>(observed_cols.size()) / static_cast<double>(max_features);
      int previous = -1;
      for (int i = 0; i < max_features; ++i) {
        int pos = static_cast<int>(std::floor(static_cast<double>(i) * step));
        if (pos >= static_cast<int>(observed_cols.size())) {
          pos = static_cast<int>(observed_cols.size()) - 1;
        }
        int col = observed_cols[static_cast<std::size_t>(pos)];
        if (col != previous) {
          distance_cols.push_back(col);
          previous = col;
        }
      }
    }

    distances.clear();
    candidate_rows.clear();
    if (n - 1 <= max_donors) {
      for (int candidate = 0; candidate < n; ++candidate) {
        if (candidate != row) {
          candidate_rows.push_back(candidate);
        }
      }
    } else {
      int take = std::min(max_donors, n - 1);
      int rank = row_rank[static_cast<std::size_t>(row)];
      int start = rank - take / 2;
      if (start < 0) {
        start = 0;
      }
      if (start + take > n) {
        start = n - take;
      }
      int end = std::min(n, start + take);
      for (int pos = start; pos < end; ++pos) {
        int candidate = row_order[static_cast<std::size_t>(pos)];
        if (candidate != row) {
          candidate_rows.push_back(candidate);
        }
      }
      int left = start - 1;
      int right = end;
      while (static_cast<int>(candidate_rows.size()) < take && (left >= 0 || right < n)) {
        if (left >= 0) {
          int candidate = row_order[static_cast<std::size_t>(left--)];
          if (candidate != row) {
            candidate_rows.push_back(candidate);
          }
        }
        if (static_cast<int>(candidate_rows.size()) >= take) {
          break;
        }
        if (right < n) {
          int candidate = row_order[static_cast<std::size_t>(right++)];
          if (candidate != row) {
            candidate_rows.push_back(candidate);
          }
        }
      }
    }

    for (std::size_t ci = 0; ci < candidate_rows.size(); ++ci) {
      int candidate = candidate_rows[ci];
      double sumsq = 0.0;
      int overlap = 0;
      for (std::size_t fi = 0; fi < distance_cols.size(); ++fi) {
        int col = distance_cols[fi];
        R_xlen_t target_idx = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + row;
        R_xlen_t candidate_idx = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + candidate;
        double candidate_value = source[static_cast<std::size_t>(candidate_idx)];
        if (is_missing_numeric(candidate_value)) {
          continue;
        }
        double diff = source[static_cast<std::size_t>(target_idx)] - candidate_value;
        sumsq += diff * diff;
        ++overlap;
      }
      if (overlap > 0) {
        distances.push_back(std::make_pair(sumsq / static_cast<double>(overlap), candidate));
      }
    }

    if (!distances.empty()) {
      std::sort(distances.begin(), distances.end(),
        [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
          if (a.first == b.first) {
            return a.second < b.second;
          }
          return a.first < b.first;
        });
    }

    for (std::size_t mi = 0; mi < missing_cols.size(); ++mi) {
      int col = missing_cols[mi];
      double sum = 0.0;
      int count = 0;
      for (std::size_t di = 0; di < distances.size() && count < k; ++di) {
        int candidate = distances[di].second;
        R_xlen_t candidate_idx = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + candidate;
        double value = source[static_cast<std::size_t>(candidate_idx)];
        if (!is_missing_numeric(value)) {
          sum += value;
          ++count;
        }
      }
      double replacement = count > 0
        ? sum / static_cast<double>(count)
        : column_fallback[static_cast<std::size_t>(col)];
      out_ptr[static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n) + row] = replacement;
    }

    if ((row & 15) == 0) {
      R_CheckUserInterrupt();
    }
  }

  UNPROTECT(1);
  return out;
}
