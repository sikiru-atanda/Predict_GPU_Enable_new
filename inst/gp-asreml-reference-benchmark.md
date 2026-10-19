# GP and ASReml-R reference benchmark

`tools/benchmark_gp_asreml_reference.R` is a licensed, local statistical
reference benchmark for the PredictProR Gaussian-process models. It covers:

1. Single-trait, single-environment GBLUP.
2. Single-trait MET FA1 with one genomic kernel.
3. Single-trait MET FA1 with genomic and omic kernels.
4. Multi-trait, single-environment unstructured GBLUP.
5. Multi-trait, multi-environment main-genetic plus GxE GBLUP.

Run it from the package root with the R installation that can check out the
ASReml-R license:

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools\benchmark_gp_asreml_reference.R
```

On macOS or Linux, including Ubuntu, use the licensed R installation on
`PATH`:

```sh
Rscript tools/benchmark_gp_asreml_reference.R
```

Set `PREDICTPRO_GP_PYTHON` when the GP environment is outside the project. If
it is unset, the runner checks the Windows and POSIX `.venv-gp` layouts and
then `python3`/`python` on `PATH`.

Set `PREDICTPRO_ASREML_LIB` when the licensed ASReml-R package is in a library
other than the active R library. This is especially important when multiple
ASReml-R versions are installed and the license server supports only one of
them. The value is an R library directory containing the `asreml` package, not
the package directory itself. The runner records the selected package version
and library path in `benchmark_bundle.rds`.

ASReml license-server configuration stays external to PredictProR. For example,
define `vsni_LICENSE` in the process environment or the user's `~/.Renviron` as
in the license administrator's instructions. Do not store license settings in
package source. The same environment-variable convention works on Windows,
macOS, Linux, and Ubuntu.

Optional arguments are `--n-genotypes=36`, `--n-holdout=6`, `--seed=20260813`,
and `--out-dir=tmp/gp_asreml_reference_benchmark`. The command writes prediction,
accuracy, timing, covariance, convergence, and engine-parity tables plus an RDS
bundle. Generated results stay under `tmp/` and are not tracked by Git.

Additional controls are `--scenarios=<comma-separated names>`, `--gp-iters=600`
(exact GP only; FA retains its backend budget of 300), `--mt-iters=100`, and
`--met-iters=150`. Native calls, fixtures and prediction-reference JSON files
are saved for audit. `--reuse-gp-results=<directory>` reuses the saved GP calls
and their original timings from an otherwise identical, complete five-scenario
run; it is intended for reference-only repairs, not additional timing replicates.
The bundle identifies this reuse. Keep the dataset, controls and GP source
unchanged when using it.

`accuracy_metrics.csv` retains one overall RMSE, MAE, Pearson correlation, and
Spearman correlation for each scenario and engine. Additional accuracy tables
are emitted only where their strata are meaningful:

- `accuracy_by_environment.csv` covers single-trait MET.
- `accuracy_by_trait.csv` covers multi-trait single-environment and multi-trait
  MET.
- `accuracy_by_environment_trait.csv` covers multi-trait MET at the individual
  environment-trait cell level.
- `accuracy_metrics_stratified.csv` combines overall and stratified results with
  a stable `metric_scope`, `Env`, and `Trait` schema.

For multi-trait MET, the trait table pools across environments and the
environment-trait table is cell-specific. Metrics are not pooled across traits
within environment because trait-specific units make that RMSE scientifically
ambiguous. The legacy overall table is retained for continuity, but trait and
environment-trait rows should be the primary accuracy interpretation when
traits have different scales. The `n` column always records the number of
finite observed-predicted pairs used in a metric row.

## Interpretation boundaries

The benchmark uses common simulated truth and the same genotype holdout for both
engines. It is a correctness and contract benchmark, not a claim that one
engine is faster in production. Fits are sequential, Python worker reuse is
disabled, and timings include each public model call, ASReml kernel preparation,
bounded refinement and prediction. Python startup/serialization is included;
R startup and license checkout are excluded. CPU execution is explicit. These
are sequential end-to-end measurements, not isolated solver or warm-worker times.

The single-trait and multi-trait single-environment models have directly
matched covariance structures. With `full_vc`, the dense GP FA path fits
environment-specific residual variances and uses them for prediction, matching
the single-kernel FA reference model. The
ASReml multi-kernel FA model estimates a separate FA covariance for each kernel,
whereas the GP FA model uses one environmental FA structure with kernel-specific
scalar weights. The simulation lies in the shared submodel, but individual FA
parameters are not expected to match one for one.

For the multi-trait MET comparison, both engines contain a genomic main effect,
an identity-environment GxE effect, and unrestricted positive-definite trait
covariances. PredictProR uses unstructured covariance matrices; ASReml uses the
equivalent heterogeneous-variance plus correlation (`corgh`) parameterization
to keep optimizer updates positive definite. PredictProR REML internally uses
one training scale per trait across environments, so prediction agreement and response-scale
covariance diagnostics are the primary comparison. Treat raw optimizer
parameters as engine-specific.

Every prediction table is normalized to:

`scenario`, `dataset`, `engine`, `model_family`, `GID`, `Env`, `Trait`,
`Observed_value`, `Predicted_value`, `Standard_error`, `PEV`,
`Train_Test_Label`.

In this benchmark's normalized table, `PEV` means the variance of the **full
latent prediction**, including estimated fixed effects, and equals the reported
standard error squared, as in native GP outputs. Genetic-only PEV is in
`Reliability_variance_input`, with its row prior in
`Reliability_reference_variance`. Missing or ambiguously labelled uncertainty stays `NA`. Native outputs
are preserved separately. ASReml predictions freeze covariance parameters; the
single-trait reference uses the sigma parameterization (`idv(units)`) to preserve
the fitted residual scale. MET uses `id(Env)` in the GxE term, avoiding an extra
unidentified variance multiplier and retaining GxE in the classify prediction.

Convergence comes from the actual GP prediction covariance fit, not success of
the separate AI reporting pass. Unavailable covariance reports remain unavailable;
finite covariance matrices must pass symmetry/PSD checks. Scenario failures are
saved and remaining scenarios continue. A returned prediction is not evidence
that its REML endpoint converged. Single-trait AI reporting estimates are labelled
separately from the covariance used for prediction.

`tools/benchmark_gp_asreml_fixed_prediction.py <result directories...>` replays
production GP prediction routines using the saved ASReml covariance estimates.
It checks prediction algebra separately from variance fitting. Multi-kernel FA
is excluded because its covariance families differ. The single-kernel FA replay
uses the equivalent fixed environmental covariance through the exact GP routine.

`tools/benchmark_gp_prediction_routes.R <result directories...>` compares
prediction-only modes on the saved holdouts. It labels the dense GP default,
single-trait memory-mapped operator, and multi-trait MoM routes separately.
`PREDICTPRO_BENCH_ROUTE` can select one route. Dense GP `predict_only` disables
the AI variance report and prediction SEs but still optimizes covariance/noise;
it is not equivalent to a no-REML request. The operator and MoM routes avoid REML.
MoM still obtains internal moment estimates needed for prediction. Memory-mapped
factor construction is included in the operator timing; this small benchmark
uses an exact factor and does not certify low-rank approximation accuracy.
The runner writes the factor in the cache's required row-major order and checks
that its Gram matrix reproduces the input kernel before timing a prediction.
Earlier operator results that used a transposed factor are invalid and are
identified as superseded in the GP statistical-repair report.

ASReml's prediction/update behavior is documented in its
[prediction reference](https://asreml.kb.vsni.co.uk/knowledge-base/predict-asreml/)
and [options reference](https://asreml.kb.vsni.co.uk/knowledge-base/asreml-options/).
