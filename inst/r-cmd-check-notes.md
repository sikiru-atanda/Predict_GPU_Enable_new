# R CMD check status

Current status on the develop branch:

```
R version 4.3.3 (2024-02-29 ucrt), gcc 12.3.0
* using options '--no-examples --no-tests --no-manual'
Status: 1 NOTE   (0 ERRORs, 0 WARNINGs)
```

The package **builds, installs and loads** cleanly:

- `checking whether package 'PredictProR' can be installed ... OK` — C++ layer
  (`grm_dense`, `kernel_dense`, `kernel_spd_repair`, `kernel_duplicate_scan`,
  `ld_prune`, `geno_qc`) compiles via Rtools43 / gcc 12.2.0 and is registered
  through `src/init.c` with `useDynLib(.registration=TRUE)`.
- `checking compiled code ... OK`.
- `checking package dependencies ... OK`.
- All namespace / load / unload / off-search-path checks pass.

The full suite (`tests/`, examples, PDF manual) is skipped here for a bounded
check; the long-path full check is the next step when needed.

## NOTEs

1. **Dependencies in R code.** The note printed alongside the asreml license
   check-out message; typically a soft note about declared imports vs `:::`
   usage or namespace declarations. No functional impact.

## Resolved

- **Hidden files** (was a NOTE; now OK) — top-level hidden files and folders
  (`^\.[^/]+(/|$)`) are excluded in `.Rbuildignore` and in
  `tools/build_clean_tarball.R`'s `hard_ignores`.
- **C++ specification / C++11** (was a NOTE; now OK) — dropped `CXX_STD = CXX11`
  from `src/Makevars` and `src/Makevars.win`. R 4.3.3 defaults to C++17 (a
  superset for this layer's matrix / LAPACK operations); compiled-code check
  still OK.
- **`Rd \usage` undocumented arguments** (was a WARNING; now OK) — filled
  `@param` descriptions across five batches: small utilities, the `AI_*`
  classical-ML wrappers, the preprocessing cluster, the Bayesian / DL / misc
  helpers, and finally the `model_execute` dispatcher (~250 args, grouped by
  family). Removed two stray `@param core` entries (not in the signature).
- **`Rd contents` `\item` entries with no description** (was a WARNING; now OK)
  — resolved as a by-product of the `@param` fill, since the empty `\item`
  entries came from the same skeleton stubs.

## Run the bounded check locally

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools/build_clean_tarball.R
$env:R_LIBS = '<project>/.r-lib;<your R user library>'
& 'C:\Program Files\R\R-4.3.3\bin\R.exe' CMD check --no-manual --no-tests --no-examples --output=tmp/clean-dist tmp/clean-dist/PredictProR_0.20.7.tar.gz
```

The full check (with tests/examples) is `R CMD check PredictProR_0.20.7.tar.gz`;
expect a longer run and additional skips for the asreml / Python-backed model
runtime-smoke tests.
