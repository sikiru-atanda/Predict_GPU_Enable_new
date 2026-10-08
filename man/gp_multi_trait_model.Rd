\name{gp_multi_trait_model}
\alias{gp_multi_trait_model}
\title{Fit a Multi-Trait Gaussian-Process Model}
\usage{
gp_multi_trait_model(
  pheno_data,
  gmatrix = NULL,
  response,
  gen_name,
  test_set = NULL,
  geno_ids = NULL,
  kernel_list = NULL,
  kernel_weights = NULL,
  estimate_kernel_weights = FALSE,
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
  residual_structure = "unstructured",
  residual_fa_rank = NULL,
  max_iter = 100L,
  tol_loglik = 1e-04,
  seed = 12345L,
  python_bin = NULL,
  project_root = NULL
)
}
\arguments{
\item{pheno_data}{Data frame with one genotype row and one column per trait.}
\item{gmatrix}{Optional genomic relationship matrix. Required for dense REML
fits. For scalable \code{varcomp_mode = "mom"} fits, it may be omitted when
\code{gp_factor_cache} and matching \code{geno_ids} are supplied.}
\item{response}{Character vector of at least two trait columns.}
\item{gen_name}{Genotype ID column name.}
\item{test_set}{Optional genotype IDs to predict.}
\item{geno_ids}{Optional genotype IDs matching the rows and columns of
\code{gmatrix}, or the rows represented by \code{gp_factor_cache} when
\code{gmatrix} is omitted.}
\item{kernel_list}{Optional named list of additional relationship/kernel matrices. When supplied, multi-trait GP fits pass all matrices as backend \code{geno_kernels}.}
\item{kernel_weights}{Optional non-negative weights for the aligned genomic
kernels. Named weights are reordered to the aligned kernel names. The default
assigns weight one to every kernel. Prediction uses the weighted kernel sum.
The backend estimates one shared trait covariance for that combined kernel, so
the exported per-kernel covariance matrices are explicitly labeled
fixed-weight decompositions rather than independent estimates.}
\item{force_prediction_se}{Logical. Prediction standard errors and per-cell trait covariance cost one linear solve per test cell per trait, so they are capped: past the solve budget the fit returns predictions without them and warns rather than running for hours. Set `TRUE` to compute them anyway.}

\item{estimate_kernel_weights}{Logical. When `TRUE`, fit the relative kernel mixture weights by REML profile likelihood instead of holding them fixed. The first kernel is held at 1 because the overall genetic scale is carried by the trait covariance, so the fitted values are ratios. Requires `varcomp_mode = "reml"` and at least two kernels; otherwise the request is ignored with a warning.}

\item{varcomp_mode}{Variance component mode. Use \code{"reml"} for prediction SE output.}
\item{gp_factor_cache}{Optional Python GRM factor-cache object or serializable
factor-cache specification for \code{varcomp_mode = "mom"}.}
\item{gp_auto_policy}{Logical; when \code{TRUE}, record and apply the
conservative GP device route selected from available GPU/CPU resources.}
\item{tune_gp}{Logical or \code{"auto"}; when enabled, run a bounded
genotype-holdout validation pass for multi-trait RMSE fallback/blending and SE
calibration. \code{"auto"} enables it only when enough training rows and
genotypes are available.}
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
\item{return_se}{Logical; request prediction standard errors where supported.}
\item{prediction_output}{One of \code{"test_only"} or \code{"all"}. The latter
returns posterior predictions for both training and test rows while the fitted
model continues to condition only on the training observations.}
\item{return_trait_correlations}{Logical; return genetic and residual trait correlation matrices.}
\item{fixed_effects}{Optional fixed-effect columns passed to the backend.}
\item{trait_structure}{Trait covariance structure. Defaults to
\code{"unstructured"}. Use \code{"fa"} only as an explicit, benchmarked
opt-in.}
\item{trait_fa_rank}{Optional factor-analytic rank when \code{trait_structure = "fa"}.}
\item{residual_structure}{Residual covariance structure.}
\item{residual_fa_rank}{Optional factor-analytic rank when \code{residual_structure = "fa"}.}
\item{max_iter}{Maximum REML iterations.}
\item{tol_loglik}{Log-likelihood tolerance.}
\item{seed}{Integer random seed.}
\item{python_bin}{Optional Python executable for the GP runtime.}
\item{project_root}{Optional path to the vendored GP Python module directory.}
}
\value{
A list with \code{predictions}, \code{result}, \code{fit}, and \code{info}.
The result includes per-kernel variance and covariance tables. Variance output
is withheld when the variance-component fit explicitly reports
non-convergence.
}
\description{
R wrapper around the single-environment multi-trait GP backend.
}
