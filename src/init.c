#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP predictpror_kernel_spd_repair(SEXP, SEXP, SEXP);
extern SEXP predictpror_kernel_duplicate_pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predictpror_ld_prune_graph(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predictpror_kernel_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predictpror_grm_dense(SEXP, SEXP, SEXP, SEXP);
extern SEXP predictpror_geno_qc_metrics(SEXP, SEXP);
extern SEXP predictpror_geno_impute_summary(SEXP, SEXP);
extern SEXP predictpror_geno_impute_knn(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"predictpror_kernel_spd_repair", (DL_FUNC) &predictpror_kernel_spd_repair, 3},
  {"predictpror_kernel_duplicate_pairs", (DL_FUNC) &predictpror_kernel_duplicate_pairs, 6},
  {"predictpror_ld_prune_graph", (DL_FUNC) &predictpror_ld_prune_graph, 5},
  {"predictpror_kernel_dense", (DL_FUNC) &predictpror_kernel_dense, 8},
  {"predictpror_grm_dense", (DL_FUNC) &predictpror_grm_dense, 4},
  {"predictpror_geno_qc_metrics", (DL_FUNC) &predictpror_geno_qc_metrics, 2},
  {"predictpror_geno_impute_summary", (DL_FUNC) &predictpror_geno_impute_summary, 2},
  {"predictpror_geno_impute_knn", (DL_FUNC) &predictpror_geno_impute_knn, 4},
  {NULL, NULL, 0}
};

void R_init_PredictProR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
