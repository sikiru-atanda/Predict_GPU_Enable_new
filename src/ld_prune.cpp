#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <vector>

static int scalar_int_ld(SEXP x, int fallback, int min_value) {
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

static double scalar_real_ld(SEXP x, double fallback, double min_value) {
  if (Rf_length(x) < 1) {
    return fallback;
  }
  SEXP coerced = PROTECT(Rf_coerceVector(x, REALSXP));
  double value = REAL(coerced)[0];
  UNPROTECT(1);
  if (!R_FINITE(value) || value < min_value) {
    return fallback;
  }
  return value;
}

static double matrix_value_ld(SEXP x, R_xlen_t idx) {
  switch (TYPEOF(x)) {
    case REALSXP:
      return REAL(x)[idx];
    case INTSXP: {
      int value = INTEGER(x)[idx];
      return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
    }
    default:
      Rf_error("geno_data must be a numeric or integer matrix");
  }
}

class DisjointSet {
 public:
  explicit DisjointSet(int n) : parent_(n), rank_(n, 0) {
    for (int i = 0; i < n; ++i) {
      parent_[static_cast<std::size_t>(i)] = i;
    }
  }

  int find(int x) {
    int px = parent_[static_cast<std::size_t>(x)];
    if (px != x) {
      parent_[static_cast<std::size_t>(x)] = find(px);
    }
    return parent_[static_cast<std::size_t>(x)];
  }

  void unite(int a, int b) {
    int ra = find(a);
    int rb = find(b);
    if (ra == rb) {
      return;
    }
    unsigned char rank_a = rank_[static_cast<std::size_t>(ra)];
    unsigned char rank_b = rank_[static_cast<std::size_t>(rb)];
    if (rank_a < rank_b) {
      std::swap(ra, rb);
    }
    parent_[static_cast<std::size_t>(rb)] = ra;
    if (rank_a == rank_b) {
      ++rank_[static_cast<std::size_t>(ra)];
    }
  }

 private:
  std::vector<int> parent_;
  std::vector<unsigned char> rank_;
};

static double ld_r2_pair(SEXP geno, int n, int col_a, int col_b, int ploidy) {
  double sum_a = 0.0;
  double sum_b = 0.0;
  double sum_ab = 0.0;
  double sum_sq_a = 0.0;
  double sum_sq_b = 0.0;
  double sum_product = 0.0;
  int count = 0;
  R_xlen_t offset_a = static_cast<R_xlen_t>(col_a) * static_cast<R_xlen_t>(n);
  R_xlen_t offset_b = static_cast<R_xlen_t>(col_b) * static_cast<R_xlen_t>(n);

  for (int row = 0; row < n; ++row) {
    double ga = matrix_value_ld(geno, offset_a + row);
    double gb = matrix_value_ld(geno, offset_b + row);
    if (ISNAN(ga) || ISNAN(gb)) {
      continue;
    }
    sum_a += ga;
    sum_b += gb;
    sum_sq_a += ga * ga;
    sum_sq_b += gb * gb;
    sum_product += ga * gb;
    sum_ab += (ga == 2.0 && gb == 2.0) ? 1.0 : 0.0;
    sum_ab += (ga == 2.0 && gb == 1.0) ? 0.5 : 0.0;
    sum_ab += (ga == 1.0 && gb == 2.0) ? 0.5 : 0.0;
    sum_ab += (ga == 1.0 && gb == 1.0) ? 0.25 : 0.0;
    ++count;
  }

  if (count == 0) {
    return NA_REAL;
  }
  double denom_count = static_cast<double>(count);
  if (ploidy != 2) {
    double centered_cross = sum_product - (sum_a * sum_b / denom_count);
    double centered_a = sum_sq_a - (sum_a * sum_a / denom_count);
    double centered_b = sum_sq_b - (sum_b * sum_b / denom_count);
    double dosage_denom = centered_a * centered_b;
    if (dosage_denom <= 0.0 || !R_FINITE(dosage_denom)) {
      return NA_REAL;
    }
    double r2 = (centered_cross * centered_cross) / dosage_denom;
    return R_FINITE(r2) ? r2 : NA_REAL;
  }
  double p_a = (sum_a / denom_count) / static_cast<double>(ploidy);
  double p_b = (sum_b / denom_count) / static_cast<double>(ploidy);
  double p_ab = sum_ab / denom_count;
  double denom = p_a * (1.0 - p_a) * p_b * (1.0 - p_b);
  if (denom == 0.0 || !R_FINITE(denom)) {
    return NA_REAL;
  }
  double d = p_ab - p_a * p_b;
  double r2 = (d * d) / denom;
  return R_FINITE(r2) ? r2 : NA_REAL;
}

extern "C" SEXP predictpror_ld_prune_graph(SEXP geno_sexp,
                                            SEXP window_sexp,
                                            SEXP r2_threshold_sexp,
                                            SEXP maf_thresh_sexp,
                                            SEXP ploidy_sexp) {
  if (!Rf_isMatrix(geno_sexp) ||
      (TYPEOF(geno_sexp) != REALSXP && TYPEOF(geno_sexp) != INTSXP)) {
    Rf_error("geno_data must be a numeric or integer matrix");
  }

  SEXP dim = Rf_getAttrib(geno_sexp, R_DimSymbol);
  int n = INTEGER(dim)[0];
  int m = INTEGER(dim)[1];
  if (n < 1 || m < 1) {
    return Rf_allocVector(INTSXP, 0);
  }

  int window = scalar_int_ld(window_sexp, 100, 1);
  double r2_threshold = scalar_real_ld(r2_threshold_sexp, 0.9, 0.0);
  double maf_thresh = scalar_real_ld(maf_thresh_sexp, 0.01, 0.0);
  int ploidy = scalar_int_ld(ploidy_sexp, 2, 1);
  if (r2_threshold > 1.0) {
    r2_threshold = 1.0;
  }
  if (maf_thresh > 0.5) {
    maf_thresh = 0.5;
  }

  std::vector<double> maf(static_cast<std::size_t>(m), NA_REAL);
  std::vector<unsigned char> pass_maf(static_cast<std::size_t>(m), 0);
  int pass_count = 0;

  for (int col = 0; col < m; ++col) {
    R_xlen_t offset = static_cast<R_xlen_t>(col) * static_cast<R_xlen_t>(n);
    double sum = 0.0;
    int count = 0;
    for (int row = 0; row < n; ++row) {
      double value = matrix_value_ld(geno_sexp, offset + row);
      if (ISNAN(value)) {
        continue;
      }
      sum += value;
      ++count;
    }
    if (count > 0) {
      double allele_freq = (sum / static_cast<double>(count)) /
        static_cast<double>(ploidy);
      double marker_maf = std::min(allele_freq, 1.0 - allele_freq);
      maf[static_cast<std::size_t>(col)] = marker_maf;
      if (R_FINITE(marker_maf) && marker_maf >= maf_thresh) {
        pass_maf[static_cast<std::size_t>(col)] = 1;
        ++pass_count;
      }
    }
  }

  if (pass_count == 0) {
    return Rf_allocVector(INTSXP, 0);
  }
  if (m == 1) {
    SEXP out_one = PROTECT(Rf_allocVector(INTSXP, pass_count));
    if (pass_count == 1) {
      INTEGER(out_one)[0] = 1;
    }
    UNPROTECT(1);
    return out_one;
  }

  DisjointSet dsu(m);
  for (int i = 0; i < m - 1; ++i) {
    if (!pass_maf[static_cast<std::size_t>(i)]) {
      continue;
    }
    int j_end = std::min(m - 1, i + window);
    for (int j = i + 1; j <= j_end; ++j) {
      if (!pass_maf[static_cast<std::size_t>(j)]) {
        continue;
      }
      double r2 = ld_r2_pair(geno_sexp, n, i, j, ploidy);
      if (!ISNAN(r2) && R_FINITE(r2) && r2 > r2_threshold) {
        dsu.unite(i, j);
      }
    }
    if ((i & 127) == 0) {
      R_CheckUserInterrupt();
    }
  }

  std::vector<int> best(static_cast<std::size_t>(m), -1);
  for (int col = 0; col < m; ++col) {
    if (!pass_maf[static_cast<std::size_t>(col)]) {
      continue;
    }
    int root = dsu.find(col);
    int current = best[static_cast<std::size_t>(root)];
    if (current < 0 ||
        maf[static_cast<std::size_t>(col)] > maf[static_cast<std::size_t>(current)] ||
        (maf[static_cast<std::size_t>(col)] == maf[static_cast<std::size_t>(current)] &&
         col < current)) {
      best[static_cast<std::size_t>(root)] = col;
    }
  }

  std::vector<int> keep;
  keep.reserve(static_cast<std::size_t>(pass_count));
  for (int col = 0; col < m; ++col) {
    int root = dsu.find(col);
    if (pass_maf[static_cast<std::size_t>(col)] &&
        best[static_cast<std::size_t>(root)] == col) {
      keep.push_back(col + 1);
    }
  }

  SEXP out = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(keep.size())));
  for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(keep.size()); ++i) {
    INTEGER(out)[i] = keep[static_cast<std::size_t>(i)];
  }
  UNPROTECT(1);
  return out;
}
