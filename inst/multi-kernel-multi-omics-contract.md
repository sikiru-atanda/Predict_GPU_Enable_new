# Multi-kernel and multi-omics execution contract

PredictProR retains every explicitly supplied genomic or omics source. The
representation used by a model depends on what that model estimates.

Repeated source names and names that become identical after sanitization are
disambiguated without merging or dropping their matrices. ML/DL feature columns
also receive unique names when omics display labels coincide. This input
preservation applies to Gaussian regression and supported binary, multiclass,
and ordinal classification routes; each model's response-family restrictions
still apply.

## Covariance and kernel models

- ASReml GBLUP fits one independent relationship term per kernel. For joint
  multi-trait models the total genetic covariance is the sum of the fitted
  per-kernel trait covariance matrices. For grouped multi-trait
  multi-environment models, each kernel has its own factor-analytic genetic
  term and the public total is verified with `asreml::vpredict()`.
- GBLUP_BRR and RKHS retain all kernels in separate BGLR ETA terms. Single-
  trait, joint multi-trait, multi-trait multi-environment, and heterogeneous-
  residual MET output exposes the independently fitted covariance contribution
  of every kernel and its posterior SD; the reported total is the sum used by
  prediction.
- Gaussian-process models retain the aligned named kernel bank and use the
  route's additive multi-kernel parameterization and recorded kernel weights
  for both CV predictions and direct/final predictions. Exact single-trait GP
  with full variance-component output estimates and exports one REML scale per
  kernel. Kernel ridge and the joint multi-trait/MET GP backends use a weighted
  combined kernel. For kernel ridge, `kernel_variance_components` marks each
  separate kernel variance as not identifiable (`Components = NA`), while
  `Fixed_weight_allocation` retains the arithmetic split for inspection.
  Equal allocations can result from equal fixed weights and do not imply equal
  fitted genetic variances. Joint multi-trait/MET per-kernel covariance tables
  are labeled fixed-weight decompositions of shared fitted covariance, never
  independent kernel estimates. Joint MT-MET covariance is always
  converted back to response scale, even when correlation matrices were not
  requested.
- Hybrid ASReml, Bayesian, and GP routes retain their structurally distinct
  female-GCA, male-GCA, SCA, and, where requested, environment/reaction-norm
  kernels. Those components are not collapsed into a marker feature table.

For covariance models, fitted genetic components are biological variance
components only when the model converges and its route-specific variance
stability checks pass. Duplicate or nearly collinear kernels can make separate
components non-identifiable; PredictProR must fail or withhold the variance
output rather than relabel prediction dispersion as genetic variance.

## Machine-learning and deep-learning models

- Raw `geno_data`, `omic1_data`, `omic2_data`, and `omic3_data` are aligned by
  sample ID and concatenated as named feature blocks.
- Each explicitly supplied PSD kernel is eigen-decomposed independently. The
  retained eigenfeatures are named by source and concatenated with the raw
  feature blocks.
- The same feature-bank contract is used by single-trait, automatic
  multi-trait, MET ML/DL, and hybrid ML/DL routes in both direct prediction and
  cross-validation.

ML/DL regression and classification do not estimate random-effect covariance
components. Their output therefore reports predictive probability or response
dispersion, predictive error, and predictive reliability diagnostics. Genetic
variance, residual variance, and heritability remain explicitly
non-identifiable.

## Marker-effect Bayesian models

BayesA, BayesB, BayesC, Bayesian LASSO, and marker BRR require raw genomic or
omics feature columns because their priors are defined on marker effects.
Kernel-only inputs must use GBLUP_BRR or RKHS instead. This is a model
definition boundary, not an input-routing limitation.

For one-trait, one-environment analysis, each supplied raw block (`geno_data`,
`omic1_data`, `omic2_data`, and `omic3_data`) is fitted as a separate BGLR ETA
term under the selected Bayes Alphabet prior. Predictions and posterior target
draws sum every term. Gaussian output reports each block's posterior genetic
variance and recomputes total genetic variance from the summed genetic effect,
which retains posterior cross-block dependence. Binary and ordinal output uses
the same blockwise posterior marker-effect draws on the latent probit-liability
scale, with residual liability variance fixed to 1 for identification.

## Reuse and provenance

Automatic multi-model and automatic multi-trait execution prepares phenotype,
feature, and kernel inputs once. Cross-validation folds and the selected final
refit reuse that immutable payload. Public model parameters record the number,
names, and combination strategy for kernel-based fits; grouped ASReml output
also retains separate per-kernel variance and covariance evidence.

Every multi-trait route reports the kernels it actually fitted, so a
multi-kernel run cannot be read as a single-kernel run. This holds for
single-environment and MET panels, for cross-validation and true prediction,
and for the ASReml, Gaussian-process, and Bayesian multi-trait engines. Two
tables carry it:

- `kernel_configuration` names each fitted kernel and its role.
- `model_parameters` carries `multi_kernel_count`, `multi_kernel_names`, and
  `multi_kernel_strategy`.

`multi_kernel_strategy` distinguishes what the engine actually estimated, because
joint multi-trait engines do not treat multiple kernels the same way:

| Route | Strategy | Per-kernel variance separately identified |
| --- | --- | --- |
| Multi-trait ASReml (single-env, MET) | `independent_<cov>_trait_covariance_per_kernel_summed`, `independent_FA_covariance_per_kernel_summed` | Yes - one `vm(GID,Kinv)` term per kernel |
| Multi-trait GBLUP_BRR / RKHS | `independent_bayesian_eta_term_per_kernel` | Yes - one BGLR ETA term per kernel, `<kernel>_G` and `<kernel>_GxE` |
| Joint multi-trait GP, `gp_estimate_kernel_weights = FALSE` | `fixed_weight_combined_kernel` | No - the bank is summed into one `K_geno` using fixed `kernel_weights` |
| Joint multi-trait GP, `gp_estimate_kernel_weights = TRUE` | `reml_estimated_weight_combined_kernel` | Ratios only - see below |
| Joint multi-trait MET GP | `fixed_weight_combined_kernel` | No - method-of-moments backend, no likelihood to profile |

With fixed weights, every supplied kernel still contributes to joint GP
predictions and the weights are recorded, but they are inputs echoed back rather
than estimates: equal weights produce equal contributions and must not be read
as equal fitted genetic variances.

### Fitting the joint GP kernel mixture

`gp_estimate_kernel_weights = TRUE` fits the mixture instead of assuming it. The
model is

    V = Sigma_G (x) sum_k w_k K_k  +  Sigma_eps (x) I

and `w[1]` is held at 1: the overall genetic scale is already carried by
`Sigma_G`, so only the ratios `w[k]/w[1]` are identified. Freeing all of them
would leave the likelihood unchanged along `Sigma_G -> c Sigma_G, w -> w/c`.

The ratios are chosen by maximising the REML *profile* likelihood - each
candidate mixture is scored by a full inner AI-REML fit of `(Sigma_G,
Sigma_eps)`, searched coordinate-wise with a bounded line search over `log w`.
The search is derivative-free on purpose: the combined kernel's
eigendecomposition depends on the weights, and differentiating through `eigh` is
ill-conditioned when that kernel has near-degenerate eigenvalues.

Requirements, all enforced: `gp_varcomp_mode = "reml"`, at least two kernels, and
a single-environment panel. The multi-environment joint GP backend estimates its
variance components by method of moments and has no likelihood to profile, so it
keeps fixed weights and warns rather than reporting an unfitted mixture as a
fitted one.

Read the result with its diagnostic. `multi_kernel_strategy` says whether the
mixture was fitted, and the fit also reports the REML log-likelihood gained over
the starting weights. That gain is the honest measure of how much the data said:
two kernels describing similar relatedness trade off along a very flat ridge, so
a near-zero gain means the mixture is not identified by this dataset, not that
equal weighting was confirmed. A ratio driven to the search bound is reported as
such rather than as an interior optimum.

Use multi-trait ASReml GBLUP or GBLUP_BRR/RKHS when separately estimated
per-kernel variance *components* (not just relative weights) are required.

Cross-validation exports these as `cv_results_processed_kernel_configuration.csv`
and `cv_results_processed_model_parameters.csv`; true prediction exports
`kernel_configuration.csv` and `model_parameters.csv`.

A `gmatrix` built from `geno_data` by `gmatrix_method` and a separately supplied
`gkernel` are distinct model inputs and are both retained; supplying `gkernel`
never replaces a computed marker-derived GRM.

Response-derived feature selection in single-environment multi-trait ASReml CV
rebuilds one genomic kernel per fold from the selected markers and does not
carry the supplied kernels into the folds. That run is reported as a single
rebuilt kernel, with the supplied kernel count recorded in
`multi_kernel_strategy`.
