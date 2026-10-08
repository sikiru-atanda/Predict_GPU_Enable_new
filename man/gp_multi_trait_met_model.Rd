\name{gp_multi_trait_met_model}
\alias{gp_multi_trait_met_model}
\title{Fit a Multi-Trait Multi-Environment Gaussian-Process Model}
\usage{
gp_multi_trait_met_model(
  pheno_data,
  gmatrix = NULL,
  response,
  gen_name,
  heter_groups,
  test_set = NULL,
  env_similarity = NULL,
  env_ids = NULL,
  env_covariates = NULL,
  reaction_norm_feature_qc = TRUE,
  kenv_kernel = "matern32",
  kenv_bandwidth = 1,
  kenv_kernel_kwargs = NULL,
  geno_ids = NULL,
  kernel_list = NULL,
  kernel_weights = NULL,
  force_prediction_se = FALSE,
  varcomp_mode = "reml",
  gp_factor_cache = NULL,
  gp_auto_policy = TRUE,
  tune_gp = "auto",
  prediction_blend = "auto",
  gp_uncertainty_calibration = "auto",
  tuning_validation_fraction = 0.25,
  tuning_max_calibration_rows = NULL,
  tuning_max_validation_rows = NULL,
  return_se = TRUE,
  prediction_output = "test_only",
  return_trait_correlations = FALSE,
  fixed_effects = NULL,
  trait_structure = "unstructured",
  trait_fa_rank = NULL,
  gxe_trait_structure = NULL,
  gxe_trait_fa_rank = NULL,
  mom_pcg_tol = 1e-04,
  mom_pcg_max_iter = 300L,
  dense_reml_max_train = 3000L,
  reml_max_iter = 30L,
  reml_tol = 1e-05,
  seed = 12345L,
  python_bin = NULL,
  project_root = NULL
)
}
\arguments{
\item{pheno_data}{Data frame with genotype, environment, and trait columns.}
\item{gmatrix}{Optional genomic relationship matrix. For scalable large-data
fits, it may be omitted when \code{gp_factor_cache} and matching
\code{geno_ids} are supplied.}
\item{response}{Character vector of at least two trait columns.}
\item{gen_name}{Genotype ID column name.}
\item{heter_groups}{Environment column name.}
\item{test_set}{Optional genotype IDs to predict.}
\item{env_similarity}{Optional environment relationship matrix. A
covariate-derived environment kernel should be supplied here when available; it
replaces any need for an additional environment FA or unstructured environment
covariance choice.}
\item{env_ids}{Optional environment IDs for \code{env_similarity}.}
\item{env_covariates}{Optional environment covariates for backend
environment-kernel construction. When supplied, they define the environment
kernel and replace any need for an additional environment FA or unstructured
environment covariance choice. If both \code{env_similarity} and
\code{env_covariates} are omitted, an identity environment matrix is used. Pass
exactly one of \code{env_similarity} or \code{env_covariates}.}
\item{reaction_norm_feature_qc}{Logical; when \code{TRUE}, use backend feature
QC before converting environment covariates to an environment kernel.}
\item{kenv_kernel}{Environment-covariate kernel family used when
\code{env_covariates} is supplied.}
\item{kenv_bandwidth}{Numeric bandwidth used by the environment kernel.}
\item{kenv_kernel_kwargs}{Optional named list of backend kernel arguments.}
\item{geno_ids}{Optional genotype IDs matching the rows and columns of
\code{gmatrix}, or the rows represented by \code{gp_factor_cache} when
\code{gmatrix} is omitted.}
\item{kernel_list}{Optional named list of additional relationship/kernel matrices. When supplied, multi-trait MET GP fits pass all matrices as backend \code{geno_kernels}.}
\item{kernel_weights}{Optional non-negative weights for the aligned genomic
kernels. Named weights are reordered to the aligned kernel names. The default
assigns weight one to every kernel, so prediction uses their weighted sum. The
backend estimates shared trait and GxE covariance coefficients for the
combined kernel; exported per-kernel components are fixed-weight
decompositions, not independent variance estimates.}
\item{force_prediction_se}{Logical. Prediction standard errors and per-cell trait covariance cost one linear solve per test cell per trait, so they are capped: past the solve budget the fit returns predictions without them and warns rather than running for hours. Set `TRUE` to compute them anyway.}

\item{varcomp_mode}{Variance-component estimation mode passed to the MT-MET GP
backend. The default \code{"reml"} reports REML-derived variance components;
\code{"mom"} remains available for scalable factor-cache fits.}
\item{gp_factor_cache}{Optional Python GRM factor-cache object or serializable
factor-cache specification. If omitted, an in-memory factor cache is built from
\code{gmatrix}.}
\item{gp_auto_policy}{Logical; when \code{TRUE}, record and apply conservative
GP policy metadata, including GPU/CPU routing and large-MET environment-kernel
guidance. Environment covariates, when supplied, remain the environment kernel
and do not trigger a second FA or unstructured environment covariance.}
\item{tune_gp}{Logical or \code{"auto"}; when enabled, run a bounded
genotype-holdout validation pass for MT-MET RMSE fallback/blending and SE
calibration. \code{"auto"} enables it only when enough training rows,
genotypes, and environments are available.}
\item{prediction_blend}{One of \code{"auto"}, \code{"none"}, or
\code{"validation"}. With successful validation, \code{"auto"} applies only
conservative validation-selected baseline blends.}
\item{gp_uncertainty_calibration}{One of \code{"auto"}, \code{"none"}, or
\code{"validation"}. With successful validation and requested SE output, scales
prediction SE/variance columns using validation residuals.}
\item{tuning_validation_fraction}{Fraction of training genotypes held out for
internal validation.}
\item{tuning_max_calibration_rows}{Optional row cap for internal calibration fits.}
\item{tuning_max_validation_rows}{Optional row cap for internal validation predictions.}
\item{return_se}{Logical; request prediction standard errors.}
\item{prediction_output}{One of \code{"test_only"} or \code{"all"}. The latter
queries every genotype-environment-trait cell from the fitted posterior while
the model continues to condition only on training rows.}
\item{return_trait_correlations}{Logical; return genetic, GxE, and residual trait correlation matrices.}
\item{fixed_effects}{Optional fixed-effect columns passed to the backend.}
\item{trait_structure}{Cross-trait covariance structure. Defaults to
\code{"unstructured"}. Use \code{"fa"} only as an explicit, benchmarked
opt-in. This is not an environment covariance structure. Dense MET REML
currently supports only \code{"unstructured"}.}
\item{trait_fa_rank}{Optional factor-analytic rank when \code{trait_structure = "fa"}.}
\item{gxe_trait_structure}{Cross-trait covariance structure for the GxE
component. When \code{NULL}, it follows \code{trait_structure}. This is not an
environment covariance structure. Dense MET REML currently supports only
\code{"unstructured"}.}
\item{gxe_trait_fa_rank}{Optional factor-analytic rank when \code{gxe_trait_structure = "fa"}.}
\item{mom_pcg_tol}{PCG tolerance for the method-of-moments operator path.}
\item{mom_pcg_max_iter}{Maximum PCG iterations.}
\item{dense_reml_max_train}{Maximum training rows for dense REML variance estimation.}
\item{reml_max_iter}{Maximum REML iterations.}
\item{reml_tol}{REML tolerance.}
\item{seed}{Integer random seed.}
\item{python_bin}{Optional Python executable for the GP runtime.}
\item{project_root}{Optional path to the vendored GP Python module directory.}
}
\value{
A list with \code{predictions}, \code{result}, \code{fit}, and \code{info}.
Response-scale covariance matrices are returned even when correlations are not
requested. Per-kernel genetic/GxE decompositions are withheld if the backend
explicitly reports non-convergence.
}
\description{
R wrapper around the multi-trait multi-environment GP backend.
}
