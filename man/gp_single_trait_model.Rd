\name{gp_single_trait_model}
\alias{gp_single_trait_model}
\title{Fit a Single-Trait Gaussian-Process Model}
\usage{
gp_single_trait_model(
  pheno_data,
  gmatrix = NULL,
  response,
  gen_name,
  heter_groups = NULL,
  test_set = NULL,
  test_set_source = NULL,
  env_similarity = NULL,
  env_ids = NULL,
  env_covariates = NULL,
  reaction_norm_feature_qc = TRUE,
  kenv_kernel = "matern32",
  kenv_bandwidth = 1,
  kenv_kernel_kwargs = NULL,
  tune_gp = FALSE,
  tuning_strategy = "auto",
  tuning_krr_lams = c(0.03, 0.1, 0.3),
  tuning_kenv_kernels = c("matern32", "matern52", "rbf"),
  tuning_weight_grid = NULL,
  tuning_objective = "rmse",
  mean_adjustment = "auto",
  prediction_blend = "auto",
  env_location = NULL,
  env_year = NULL,
  tuning_validation_fraction = 0.25,
  tuning_max_calibration_rows = NULL,
  tuning_max_validation_rows = NULL,
  tuning_max_genotypes = NULL,
  geno_ids = NULL,
  kernel_list = NULL,
  include_components = NULL,
  fixed_effects = NULL,
  observation_weights = NULL,
  model_name = "KRR",
  gp_backend = "auto",
  gp_output_level = "predict_only",
  gp_return_se = FALSE,
  gp_uncertainty_calibration = "auto",
  gp_full_vc = FALSE,
  gp_factor_cache = NULL,
  gp_tiered_dispatch = FALSE,
  gp_tiered_verbose = FALSE,
  gp_auto_policy = TRUE,
  gp_dtype = "float64",
  large_n_threshold = 50000L,
  operator_tol = 1e-05,
  operator_max_iter = 500L,
  operator_dtype_compute = "float32",
  gp_se_hutchinson_probes = 32L,
  gp_varcomp_mode = "reml",
  gp_fa_rank = 1L,
  gp_prediction_output = "test_only",
  krr_lam = 0.001,
  krr_lams = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10),
  lam_select = "gcv",
  gp_iters = NULL,
  gp_lr = NULL,
  w_g = NULL,
  w_ge = NULL,
  w_e = 0,
  gp_learn_scales = NULL,
  seed = 12345L,
  python_bin = NULL,
  project_root = NULL,
  gp_theta0_warm = NULL,
  gp_engine = "auto"
)
}
\arguments{
\item{pheno_data}{Data frame containing genotype IDs, the response, and optionally an environment column.}
\item{gmatrix}{Genomic relationship matrix. May be \code{NULL} for
operator/factor-cache fits when \code{gp_factor_cache} and \code{geno_ids} are
supplied.}
\item{response}{Response column name.}
\item{gen_name}{Genotype ID column name.}
\item{heter_groups}{Optional environment column name.}
\item{test_set}{Optional genotype IDs to predict. If omitted, rows with missing response values are treated as prediction rows.}
\item{test_set_source}{Optional internal provenance for \code{test_set}. When set to \code{"inferred_missing_response"}, the backend uses missing response rows as cell-level prediction rows rather than holding out every record for those genotype IDs.}
\item{env_similarity}{Optional environment relationship matrix for
multi-environment reaction-norm prediction. For large MET with more than five
environments, supply an environment kernel from covariates or a defensible
similarity matrix instead of relying on an identity environment kernel.}
\item{env_ids}{Optional environment IDs for \code{env_similarity}.}
\item{env_covariates}{Optional environment covariates for backend
environment-kernel construction. When supplied, the backend QC-aligns them to
the phenotype environments and converts them to \code{env_similarity}. Pass
exactly one of \code{env_similarity} or \code{env_covariates}.}
\item{reaction_norm_feature_qc}{Logical; when \code{TRUE}, use backend feature
QC before converting environment covariates to an environment kernel.}
\item{kenv_kernel}{Environment-covariate kernel family used when
\code{env_covariates} is supplied.}
\item{kenv_bandwidth}{Numeric bandwidth used by the environment kernel.}
\item{kenv_kernel_kwargs}{Optional named list of backend kernel arguments.}
\item{tune_gp}{Logical or \code{"auto"}; when \code{TRUE}, run historical
tuning for MET GP fits. The tuning split uses only finite training rows and
masks validation responses before candidate fits. With \code{"auto"},
PredictPro enables tuning only for multi-environment fits with enough training
rows and environments for a defensible internal validation. Environment
covariates are required for prospective/new-environment prediction, but not
for known-environment genotype holdout.}
\item{tuning_strategy}{Tuning split strategy. \code{"auto"} uses
blocked-environment validation for prospective/new-environment prediction and
blocked-genotype validation for known-environment genotype prediction.}
\item{tuning_krr_lams}{Positive ridge/lambda values to evaluate during
blocked tuning.}
\item{tuning_kenv_kernels}{Environment-covariate kernel families to evaluate
during blocked tuning.}
\item{tuning_weight_grid}{Optional data frame with columns \code{w_g},
\code{w_ge}, and \code{w_e} for component-weight tuning.}
\item{tuning_objective}{Validation objective used to select the final tuning
candidate. \code{"auto"} chooses a conservative objective from the prediction
task. \code{"rmse"} preserves the historical raw-RMSE behavior.
\code{"centered_rmse"} removes validation environment means before scoring.
\code{"within_env_spearman"} optimizes weighted within-environment rank
correlation. \code{"centered_cor"} maximizes centered Pearson correlation.}
\item{mean_adjustment}{One of \code{"auto"}, \code{"none"}, or
\code{"location_offset"}. With tuning, \code{"auto"} evaluates both no
adjustment and a training-location mean offset. Without tuning, \code{"auto"}
leaves the response unchanged.}
\item{prediction_blend}{One of \code{"auto"}, \code{"none"}, or
\code{"validation"}. With successful MET tuning, \code{"auto"} only applies
conservative baseline-dominant environment/location mean blends for
prospective/new-environment prediction when they clearly improve the internal
validation objective. Known-environment genotype prediction keeps the GP signal
by default. \code{"validation"} follows the validation-selected blend directly.}
\item{env_location}{Optional environment-location annotation. May be a
phenotype column name, named vector keyed by environment, or data frame with
environment IDs and locations. If omitted, locations are inferred from
environment names by removing a year token.}
\item{env_year}{Optional environment-year annotation. May be a phenotype
column name, named vector keyed by environment, or data frame with environment
IDs and years. If omitted, years are parsed from environment names when
possible.}
\item{tuning_validation_fraction}{Fraction of training genotypes/environments
targeted for blocked validation.}
\item{tuning_max_calibration_rows}{Optional maximum calibration rows for
internal GP tuning. When \code{NULL}, dense tuning keeps the legacy full
blocked split and factor-cache tuning uses a bounded production-safe subset.}
\item{tuning_max_validation_rows}{Optional maximum validation rows for
internal GP tuning.}
\item{tuning_max_genotypes}{Optional maximum unique genotypes used by the
dense tuning kernel materialized from a factor cache.}
\item{geno_ids}{Optional genotype IDs matching the rows and columns of \code{gmatrix}.}
\item{kernel_list}{Optional named list of additional relationship/kernel matrices. When supplied, GP fits pass all matrices as backend \code{geno_kernels}.}
\item{include_components}{Optional GP components, for example \code{c("g", "ge", "e")}.}
\item{fixed_effects}{Optional fixed-effect columns passed to the backend. When \code{NULL}, multi-environment fits use the environment column as a fixed effect if all environments are represented in the training rows.}
\item{observation_weights}{Optional Stage 2 precision weights. Supply a
positive numeric vector, a one-column table, or a column name from
\code{pheno_data}. The GP residual diagonal is fixed to
\code{1 / observation_weights}.}
\item{model_name}{GP backend model label. One of \code{"KRR"}, \code{"GP"},
or \code{"GP_FA"}. \code{"LowRankGP"} is not a distinct single-trait
estimator and is rejected; use \code{"KRR"}. It remains available through
the purpose-specific joint multi-trait route.}
\item{gp_backend}{Python GP backend selection.}
\item{gp_output_level}{One of \code{"predict_only"}, \code{"predict_with_se"}, or \code{"full_vc"}.}
\item{gp_return_se}{Logical; request prediction standard errors and prediction error variances.}
\item{gp_uncertainty_calibration}{One of \code{"auto"}, \code{"none"}, or
\code{"validation"}. With successful GP tuning and requested SE/PEV output,
\code{"auto"} applies validation-derived calibration where available.}
\item{gp_full_vc}{Logical; request full variance-component output where supported.}
\item{gp_factor_cache}{Optional low-rank genomic factor cache for large
operator-backend fits. When supplied with \code{geno_ids}, \code{gmatrix} may
be \code{NULL}; the backend uses the factor cache instead of a dense genomic
relationship matrix.}
\item{gp_tiered_dispatch}{Logical; route MET predictions by environment
connectivity. This is intended for large prospective MET exercises where
all-new test environments should be served by the covariate-derived
environment kernel without fitting an unstructured environment covariance.}
\item{gp_tiered_verbose}{Logical; request verbose tiered-dispatch backend messages.}
\item{gp_auto_policy}{Logical; when \code{TRUE}, let PredictPro choose
conservative backend, dtype, device route, tiered dispatch, and optional
\code{tune_gp = "auto"} behavior from the data shape. Set to \code{FALSE} for
fully manual settings.}
\item{gp_dtype}{Floating-point dtype passed to the Python backend.}
\item{large_n_threshold}{Training-row threshold used by the automatic
operator backend policy.}
\item{operator_tol}{PCG/operator solve tolerance for large operator fits.}
\item{operator_max_iter}{Maximum PCG/operator iterations for large operator fits.}
\item{operator_dtype_compute}{Compute dtype for operator-backend matrix-vector
work.}
\item{gp_se_hutchinson_probes}{Integer number of Hutchinson probes used for
scalable operator/tiered prediction SEs. Dense paths ignore this argument.}
\item{gp_varcomp_mode}{Variance-component estimation mode passed to the GP
backend. Gaussian-process fits use the requested mode; dense
\code{model_name = "KRR"} fits instead report the variance ratio selected by
\code{lam_select}, because a separate REML post-fit would not describe the KRR
predictor. Use \code{"mom"} only for the scalable factor-cache path.}
\item{gp_fa_rank}{Factor-analytic rank for compatible GP backends.}
\item{gp_prediction_output}{Prediction output scope passed to the backend.}
\item{krr_lam}{Numeric KRR/GP ridge parameter used when \code{lam_select = "fixed"}.}
\item{krr_lams}{Numeric candidate residual-to-genetic variance ratios used by
the KRR backend when \code{lam_select} requests automatic selection. The
default spans 0.001 through 10.}
\item{lam_select}{KRR lambda selection mode. The default \code{"gcv"} uses
deterministic exact generalized cross-validation for dense fits up to 2,048
training rows; \code{"fixed"} uses \code{krr_lam}.}
\item{gp_iters}{Optional positive integer optimizer iterations for the exact
GP backend. \code{NULL} preserves the Python backend default.}
\item{gp_lr}{Optional positive learning rate for the exact GP backend.
\code{NULL} preserves the Python backend default.}

\item{gp_learn_scales}{Optional logical override for backend scale learning.
When \code{NULL}, scales are learned only for full variance-component output.}
\item{gp_theta0_warm}{Optional numeric starting values for the GP
variance-component optimizer. Intended for compatible warm-started fits.}
\item{gp_engine}{Character; which REML engine drives GP variance-component
estimation. \code{"auto"} (default) picks the fastest available (currently
resolves to \code{"dense_v"}); \code{"dense_v"} keeps the dense V-formulation;
\code{"mme"} uses the sparse Mixed Model Equations engine (fastest at
n_geno = 100 for MET but regresses at n_geno >= 300 on GBLUP because G^-1
is dense -- see Phase 3.14 NEWS); \code{"eigen"} (Phase 3.15a, experimental)
uses an eigen-projected REML formulation that beats ASReml at large n on
the standalone validation but production integration is incomplete
(Phase 3.15b NEWS) so \code{"eigen"} currently falls back to dense_v.
Requires \code{PREDICTPRO_GP_DISABLE_ENV_MAIN=1} for \code{"mme"} and
\code{"eigen"} on MET. Honors \code{PREDICTPRO_GP_ENGINE} env var.}
\item{w_g}{Optional genomic-kernel component weights.}
\item{w_ge}{Optional genotype-by-environment component weights.}
\item{w_e}{Optional environment component weight.}
\item{seed}{Integer random seed.}
\item{python_bin}{Optional Python executable for the GP runtime.}
\item{project_root}{Optional path to the vendored GP Python module directory.}
}
\value{
A list with \code{predictions}, \code{result}, \code{fit}, and \code{info}.
}
\description{
R wrapper around the vendored GP backend for single-trait, single- or multi-environment Gaussian-process prediction.
}
