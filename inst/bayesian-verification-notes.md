# Bayesian Verification Notes

These notes record the current Bayesian verification workflow for this project.

## Verified on this machine

- Date verified: 2026-08-15
- Project root: the package source directory
- R version used successfully: `R 4.3.3`
- Working Rscript: `C:/Program Files/R/R-4.3.3/bin/Rscript.exe`

## Supported Bayesian scope

- Single-environment only:
  - `BRR`
  - `BayesA`
  - `BayesB`
  - `BayesC`
  - `BL`
- Multi-environment:
  - `RKHS`
  - `GBLUP_BRR`
- Joint multi-trait and joint multi-trait multi-environment:
  - `RKHS`
  - `GBLUP_BRR`
- Supported input decompositions for the kernel routes:
  - one genomic/relationship kernel
  - multiple named kernels
  - genomic plus one or more omics kernels

Raw marker-regression BayesA, BayesB, BayesC, BL, and BRR are deliberately
single-environment workflows. MET and joint multi-trait Bayesian workflows use
RKHS or marker-derived `GBLUP_BRR` kernels.

## What the verifier checks

`tools/verify_bayesian_models.R` has two default layers of checks and an
optional exhaustive variance-extraction layer:

1. Direct BGLR structure verification
- Confirms the actual returned object fields for:
  - `BRR`
  - `BayesA`
  - `BayesB`
  - `BayesC`
  - `BL`
  - `RKHS`
- Confirms expected on-disk files emitted by each model type.

2. Package extraction verification
- Single-environment true prediction for:
  - `BRR`
  - `BayesA`
  - `BayesB`
  - `BayesC`
  - `BL`
- Single-environment CV prediction path for the same models
- Multi-environment true prediction for `RKHS`
- Multi-environment CV path for `RKHS`
- Multi-environment prep and CV path for `GBLUP_BRR`

3. Direct posterior variance-extraction audit
- Run with `--full-variance-audit` to execute
  `tools/audit_bglr_variance_extraction.R` after the default verifier.
- Fits the complete supported cross-product of marker priors, marker
  multi-omics, RKHS/`GBLUP_BRR`, single/multiple kernels, kernel multi-omics,
  MET, joint multi-trait, and joint multi-trait MET.
- Compares 141 package genetic, GxE, and residual variance estimates with
  direct calculations from the same BGLR posterior chains or fitted covariance
  objects. Every comparison must be finite and agree within `1e-8`.
- The 2026-08-15 release run completed with `141 passed; 0 failed`.

## Variance extraction definitions

- Marker and marker-multi-omics component variance at draw `s` is
  `var(X_j %*% beta_j[s, ])` across training observations.
- Marker/BRR total multi-component variance is
  `var(sum_j X_j %*% beta_j[s, ])`, retaining posterior covariance between
  fitted components; it is not assumed to equal the sum of component means.
- RKHS component variance is the saved BGLR `varU` draw multiplied by the mean
  fitted kernel diagonal over training observations. Independent RKHS kernel
  component draws are summed for the total.
- MET environment-specific genetic variance is the corresponding diagonal of
  the fitted marginal genetic covariance after kernel scaling.
- Joint multi-trait genetic and GxE variance is the appropriate diagonal of
  the sum of scaled BGLR ETA covariance matrices (`Cov$Omega`).
- Residual variance is the posterior mean of BGLR's saved `varE` chain for
  univariate fits, or the appropriate diagonal of `resCov$R` for Multitrait
  fits. Empirical `var(y - yhat)` is not used as residual variance.
- Variance-component standard errors are posterior standard deviations of the
  matching saved component draws.

## Current uncertainty conventions

- Marker-regression models:
  - EBV uncertainty comes from posterior draws of `Xb`
  - prediction intervals are based on posterior predictive summaries
- RKHS / `GBLUP_BRR`:
  - EBV uncertainty comes from BGLR latent effect outputs `u` and `SD.u`
  - response prediction uncertainty remains tied to fitted-response summaries
- Variance component uncertainty is currently reported on posterior SD scale

## Current MET RKHS default

- For multi-environment `RKHS` with random structure like `~ GID + GID:Env`,
  PredictProR now compiles the Bayesian ETA as:
  - `GID` main effect: `BRR`
  - `GID:Env` interaction: `RKHS`
- This is intentional.
- The older MET `RKHS` setup used the same GRM-derived RKHS backbone for both
  `GID` and `GID:Env`, which made the main and interaction terms too redundant
  and depressed the fitted marginal genetic variance and reliability.
- The mixed `BRR` + `RKHS` default avoids that duplication while preserving a
  nonlinear kernel interaction term for `GID:Env`.
- Gaussian `PEV`, `Genetic_variance`, and `Reliability` for the MET Bayesian
  outputs are now computed from that mixed ETA structure and verified against
  the corresponding direct `BGLR` fit.
- On the wheat validation run that motivated this change, plain `RKHS + GRM`
  moved from the older low-reliability regime to materially higher values:
  - within-environment mean reliability around `0.52`
  - across-environment mean reliability around `0.66`
- Those values are still model- and dataset-specific, but they provide a
  practical sanity check that the MET `RKHS` path is no longer being depressed
  by redundant main-effect and interaction kernels built from the same GRM
  backbone.

## Interpreting low BayesB reliability

- In PredictProR, Bayesian `Genetic_variance`, `PEV`, and `Reliability` are
  target-specific.
- `Reliability` is computed as:
  - `1 - PEV / Genetic_variance`
- A low or zero `BayesB` reliability does not automatically mean the predicted
  values are wrong.
- It means the posterior target-specific prediction error variance is close to,
  or larger than, the corresponding target-specific genetic variance for that
  fitted model and dataset.
- This can happen when shrinkage/sparsity is strong, the testing target is weakly
  informed, or the estimated target-level genetic variance is small relative to
  posterior uncertainty.
- Users should interpret low reliability as a warning about uncertainty for that
  prediction target, not as automatic evidence that the model fit failed.

## Run the synthetic verifier

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools/verify_bayesian_models.R
```

Expected final line:

```text
verify_bayesian_models: ok
```

Run the complete release audit with:

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools/verify_bayesian_models.R --full-variance-audit
```

Expected final lines include:

```text
AUDIT SUMMARY: 141 passed; 0 failed; 141 comparisons.
audit_bglr_variance_extraction: ok
verify_case: full-variance-audit ok
verify_bayesian_models: ok
```

The audit intentionally uses short deterministic synthetic chains to verify
extraction identity, not inferential convergence. Scientific analyses still
require adequate production MCMC length, trace/effective-sample-size review,
and preferably replicated chains.

## Suggested real-data spot checks

The synthetic verifier is the main regression guard. For project data, use the wheat dataset:

- `<data dir>/WheatPhenoGenoOLD.Rdata`

Recommended manual spot checks:

- single-environment `BayesA` or `BRR` true prediction on `F5I`
- multi-environment `RKHS` true prediction on `B2IR`, `F5I`, `B5I`
- multi-environment `GBLUP_BRR` CV prep and CV prediction path

## Local R packages added for verification

Installed into project library `.r-lib` during verification:

- `brio`
- `waldo`
- `labeling`

These were needed so test and plotting-related verification paths could run cleanly.
