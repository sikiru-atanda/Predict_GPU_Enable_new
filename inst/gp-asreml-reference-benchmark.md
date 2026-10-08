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
disabled, and timings include each public model call.

The single-trait, single-environment; single-kernel FA1 MET; and multi-trait
single-environment models have directly matched covariance structures. The
ASReml multi-kernel FA model estimates a separate FA covariance for each kernel,
whereas the GP FA model uses one environmental FA structure with kernel-specific
scalar weights. The simulation lies in the shared submodel, but individual FA
parameters are not expected to match one for one.

For the multi-trait MET comparison, both engines contain a genomic main effect,
an identity-environment GxE effect, and unrestricted positive-definite trait
covariances. PredictProR uses unstructured covariance matrices; ASReml uses the
equivalent heterogeneous-variance plus correlation (`corgh`) parameterization
to keep optimizer updates positive definite. PredictProR internally standardizes
each environment-trait cell, so prediction agreement and response-scale
covariance diagnostics are the primary comparison. Treat raw optimizer
parameters as engine-specific.

Every prediction table is normalized to:

`scenario`, `dataset`, `engine`, `model_family`, `GID`, `Env`, `Trait`,
`Observed_value`, `Predicted_value`, `Standard_error`, `PEV`,
`Train_Test_Label`.

Under the package output convention, `PEV` is the square of the reported
prediction standard error. A missing uncertainty value remains `NA`; it is not
silently replaced with zero. Covariance estimates must be symmetric and
positive semidefinite within numerical tolerance, and failures stop the run.
