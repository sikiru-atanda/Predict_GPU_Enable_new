# Gaussian Process Mixed-Model Framework

This project exposes three public model routes:

- `gp_exact`: exact GP posterior prediction with optional GP posterior SEs.
- `krr_exact`: kernel-ridge / GP-equivalent posterior prediction for fast dense problems.
- `gp_icm_fa`: exact GP posterior prediction with factor-analytic environment covariance.

## Output Levels

Use `output_level` to control work:

- `predict_only`: return predictions only.
- `predict_with_se`: return predictions plus GP posterior prediction SEs.
- `full_vc`: estimate covariance parameters by REML, then compute GP predictions and SEs using that same covariance.

Prediction means and SEs come from Gaussian conditioning. With `full_vc`, the
optional AI refinement precedes conditioning: estimates, parameter SEs and
predictions describe one fitted covariance. Without `full_vc`, no AI reporting
pass is run. KRR retains its separate GCV/fixed-penalty contract.

For `gp_exact` and `gp_icm_fa`, the prediction covariance is itself fitted to a
profiled REML objective when covariance parameters are trainable. Deterministic
Cholesky fits use a bounded Adam warm-up followed by strong-Wolfe L-BFGS;
large systems using stochastic likelihood calculations retain Adam. The
iteration budget is shared between stages. Optimizer evaluation counts and
final gradient diagnostics are returned separately.
After a requested AI refinement, `initial_prediction_reml_*` identifies the
initial optimizer, and `prediction_covariance_source="reported_reml_fit"`
identifies the final covariance. `prediction_reml_*` diagnostics
describe the covariance used for prediction; an iteration limit does not imply
convergence. Multi-trait dense fitting uses a full-Hessian Newton REML method
(the internal `_fit_mt_ai_reml` name is historical). KRR/GCV and operator MoM
routes use their stated estimators rather than REML.

Variance inference returns an explicit `inference_status`. Nonconvergence or
singular information makes parameter SEs unavailable; adding ridge information
or censoring a large SE because its Wald z is small is not allowed. Prediction
SEs remain plug-in estimates conditional on the fitted covariance parameters.
Genetic PEV is separate from uncertainty of the total fixed-plus-genetic mean.
Reliability uses the row's matching genetic prior variance, including the kernel
diagonal, rather than observed phenotypic variance. Fixed-scale GP routes
without an estimated genetic scale leave reliability unavailable with the
reason `genetic variance not estimated on this route`; their SEs are
conditional on the fixed covariance. KRR and operator variance-ratio routes
profile `y'P_lambda y/(n-p)` and scale SEs, PEVs and priors consistently.

REML CV can set `fit_reml=True` while returning `predict_only` output: it
then fits the same estimator as the full REML refit, including heterogeneous
FA residuals. This is an explicit estimator choice, separate from choosing
which output tables to return. FA kernel weights sum to one so loadings and
uniquenesses own the fitted covariance amplitude.

MT-MET REML permits 300 iterations, checks a gradient scaled per residual
degree of freedom, and checks the PSD-projected covariance score at a
boundary. A confirmed boundary estimate is distinguished from a fit that
has not converged; ordinary interior Wald SEs are unavailable at a boundary.
CV0 SEs exclude the new environment's unobserved deviation.

The multi-trait REML routes return `result.varcomp` and
`result.correlation_parameters`, with delta-method SEs using the joint parameter
covariance. These SEs condition on kernel weights and training standardization;
they do not include estimated kernel-weight uncertainty. Boundary covariance
fits explicitly withhold ordinary interior Wald SEs. MET REML uses a common
training scale per trait across environments to preserve the shared covariance
family on the response scale. MoM retains its environment/trait scaling.
For FA-derived variance and heritability, free loadings still contribute
uncertainty when a uniqueness is fixed at zero. These delta-method SEs are
explicitly conditional on the selected active set and fixed kernel weights;
they are not unconstrained Wald SEs for the boundary parameter itself.

The matched exact-GP inference path currently uses `dense_v`. A request for
`mme` or `eigen` warns and falls back until those reporting engines preserve the
same model contract. This does not add dense REML work to operator/MoM prediction.

## Scalable Backend

For large `n`, pass `grm_factor_cache` and use `backend="auto"` or `backend="operator"`. The operator backend uses low-rank GRM roots and supports:

- point predictions;
- diagonal prediction SEs via Hutchinson probing when requested;
- no full prediction covariance.

The canonical import path is `mixed_model_gpu_large_scale_backend.py`; the dated backend implementation remains as the implementation file.

## Multi-Trait MET Status

`fit_multi_trait_gp_met` is a MET-specific GP/operator path, not a full `gp_icm_fa` replacement. It supports scalable multi-trait, multi-environment posterior mean prediction with `varcomp_mode="mom"`, optional exact posterior SE columns via `return_se=True`, optional FA projection for trait covariance matrices, and multi-root GRM caches through `kernel_weights`.

For small/medium problems, `varcomp_mode="reml"` runs a bounded dense REML variance-estimation pass and then returns to the GP/operator prediction path. It refuses large training sets through `dense_reml_max_train` instead of silently doing cubic work.

User-supplied `env_similarity` should be passed as a labeled `pandas.DataFrame` or with `env_ids` so the environment kernel is aligned correctly.

The MET path requires one pre-aggregated row per `(gid, env, trait)` cell. Prediction rows may have missing `y`; training rows must have finite `y`.

Multi-trait fixed effects retain numeric covariates as slopes and use training
levels for categorical contrasts. Aliased columns are removed; unestimable
prediction rows are rejected. MET REML retains environment-trait mean columns
even after centering. For a cell with no training observations (including CV0),
its mean target is the pooled training trait mean, with the corresponding
sample-count weighted GLS correction; it is not a newly estimated environment
effect. Prediction fixed effects use GLS under the fitted covariance in both
REML and MoM routes.

In all prediction tables, `PEV` equals the squared latent prediction SE.
`Reliability_variance_input` and `cov_traits` describe genetic prediction-error
variance (genetic plus GxE for MET), accounting for estimated fixed effects.
`SE_latent` and `Prediction_Var_latent` concern the total fixed-plus-genetic
predictor. `SE`/`Prediction_Var_observed` additionally include independent
observation noise. Reliability is clipped `1 - Reliability_variance_input /
Reliability_reference_variance`, using the row genetic prior. Its basis starts
with `genetic_prior_variance` or `genetic_plus_gxe_prior_variance`. Public R
`Genetic_variance` also means this row prior, not a kernel-independent scalar.
Unavailable genetic uncertainty yields NA reliability with a reason. All these quantities
condition on fitted variance parameters; variance-component estimation
uncertainty is not included. Operator accuracy depends on PCG convergence;
`pcg_all_converged` and related diagnostics cover every prediction/uncertainty
right-hand side, not only the final residual solve.

Single-trait dense and operator GLS means use unpenalized fixed effects, so
changing a covariate's units does not change its fitted mean. PCG verifies
convergence with the actual equation residual; a drifting recurrence residual
triggers a restart within the same iteration budget. A precision or iteration
limit is reported as non-convergence, without loosening the requested tolerance.

The streaming pivoted-Cholesky cache builder copies each working column before
residualization, preserving the input kernel and supporting read-only memmaps.
Raw factor caches store `Phi` in row-major order with `K ~= Phi @ Phi.T`.
For R's column-major binary writer, transpose `Phi` before writing it. Kernel
approximation error and the chosen kernel/noise parameters still affect
prediction accuracy independently of the linear solve.

R batched CV, refits and tiered calls use the same GCV candidate grid by default;
an explicit fixed lambda is preserved. Selection uses training rows only.
Operator selection is bounded by `operator_gcv_max_rows` (or uses exact bounded
feature space when feasible). The selected lambda, method, selection size and
grid-edge flag are reported; grid-edge selections warn. This prediction route
does not estimate REML variance components.

FA rank guards require `n*(k+1)-k*(k-1)/2 <= n*(n+1)/2` for `n` groups and
rank `k`. Thus FA1 is refused at two groups and FA2 at three or four groups.
Equality is saturated (FA3 at six traits), not a guarantee of identification
on a dataset. This applies to genetic/residual trait FA, environment FA, and
MT-MET MoM genetic/GxE trait FA.

FA loadings have a lower-triangular leading block; rank and sample support must
still permit identification. A random covariance supported entirely by the
fixed-effect space is fixed at zero, labelled, and removed from prediction
priors without dropping the other component estimates. Unavailable inference
always carries a status and warning. Integer fixed-effect columns with at
least three values warn because their default interpretation is a slope.

The R result surfaces `covariance_parameters`, `correlation_parameters`, and
`inference_summary`, including boundary conditioning. Exact-GP AI inference
honors positive `PREDICTPRO_GP_MAX_ITER`, `PREDICTPRO_GP_TOL_LOGLIK`, and
`PREDICTPRO_GP_TOL_THETA` overrides. `PREDICTPRO_GP_DISABLE_ENV_MAIN` is no
longer needed; component masks and identifiability checks govern that term.

## Prediction-only and variance-estimation controls

Output selection and variance estimation are separate controls. Dense
`gp_exact`/`gp_icm_fa` with `output_level="predict_only"` skip the AI reporting
pass and prediction SEs, but still optimize their prediction covariance/noise.
They must not be described as wholly REML-free. Explicit single-trait
`backend="operator"` with a factor cache avoids REML; multi-trait and MT-MET
use `varcomp_mode="mom", return_se=False` for their non-REML prediction route.
MoM estimates internal covariance coefficients by moments. None of these flags
make GP prediction independent of a kernel and noise specification.

Preserve the distinction between fast prediction and optional statistical
reporting. Disabling reports must not silently trigger an AI reporting pass,
and comparisons must name the actual estimator and include factor preparation
when it is part of the timed workflow.

## Prediction penalty selection and dense routing

Operator predictions now honor `lam` / `krr_lam`. With `lam_select="gcv"`,
the penalty is selected using training responses only, without REML. Exact
GCV uses a fixed-effect residual projection on at most
`operator_gcv_max_rows=512` observations. Genotype-only kernels with at most
512 training genotypes and four million feature entries instead use exact
feature-space GCV on all training observations. Other larger problems use
an explicitly labelled, seed-controlled training subset. Diagnostics retain
the candidates, scores, sample size and selected penalty. Supplied residual
covariances take precedence; additional random terms currently require a
fixed penalty. For a fixed operator penalty, set `lam_select="fixed"` and
`krr_lam` explicitly (the public default is 0.001).

The dense exact GP route uses deterministic Cholesky through
`dense_cholesky_max_size=2048`, allowing its bounded Adam/L-BFGS schedule to
finish moderate fits. Larger fits retain the iterative/stochastic path and
must still pass their reported convergence checks. The operator prediction
route is unaffected by this dense threshold.

## Exact multi-trait training solves

Single-environment REML uses the training-kernel eigenbasis for complete trait
grids. For unique observed cells with at most 25% missing, it uses an exact
missing-cell precision correction when its workspace estimate fits
`PREDICTPRO_GP_MISSING_MAX_MB` (default 512 MiB). The estimate accounts for
stored precision columns, missing-block factors, kernel/eigen storage and
block-solve intermediates. It replaces the fixed missing-count and eight-million-
entry limits. Precision columns are built in blocks of 32 rather than allocating
the full indicator matrix. This is a workspace allowance, not a hard process
memory cap; autograd transforms, the full GRM and query workspaces also consume
memory. Only the missing block is factored in
addition to the trait-sized eigenbasis blocks. Larger gaps or duplicate cells
use dense observed-data Cholesky. These thresholds choose exact algorithms;
they do not change the fitted covariance or add an approximation. Native
`training_solver` identifies `balanced_kronecker`, `missing_cell_correction`,
or `dense_observed`. Whole-line holdouts do not count as gaps in the training
grid. Training observation order and nonconsecutive genotype indices are
preserved. Fitting, covariance inference and batched prediction share the solve.

Dense unstructured fallback takes exact average-information steps using
first derivatives of the small trait covariances and applications of V_i to Py.
It shares the AI calculation used for final covariance SEs. Objective, jitter,
parameter coordinates, step regularization, line search and score stopping rules
are unchanged. Structured/FA routes retain their Newton policy; MET already
uses first-order L-BFGS. The internal REML fitter's `optimizer_step_basis`
records which step matrix it used.

This reduces covariance algebra costs, not the number of covariance parameters
or a guaranteed number of optimizer iterations. Accurate converged capacity and ASReml comparisons require separate
validation. The non-REML MoM/operator prediction choices remain available.

The shared-full-data exact-GP fast-CV shortcut has been removed. Enabling the
legacy `PREDICTPRO_GP_EXACT_FAST_CV` control, including `force`, warns and uses
independent training-fold refits. KRR's eligible block-CV routes are separate.
Operator fits ignore dense GP `learn_scales`; their own regularization and
response-scale estimation are unchanged.

## Core Dependencies

Runtime: `numpy`, `scipy`, `pandas`, `torch`, `gpytorch`, `linear_operator`, `zarr`.

Tests/benchmarks additionally use `pytest` and optional data readers such as `pyreadr`.
