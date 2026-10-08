# ASReml-R Local Setup Notes

These notes record the working ASReml-R configuration verified for this project on this machine.

## Verified working setup

- Date verified: 2026-04-22
- Project root: the package source directory
- R version used successfully: `R 4.3.3`
- Working Rscript: `C:/Program Files/R/R-4.3.3/bin/Rscript.exe`
- Working ASReml-R package location: `<your R user library>/asreml`
- Verified ASReml-R version: `4.2.0.279`

## License configuration

The machine already had a valid ASReml-R license configuration:

- `.Renviron` entry: `vsni_LICENSE=5053@horus-hub.bigdata.ag.ndsu.edu`
- License token file found at:
  `<path to your VSNi license>/asreml_r.json`

If a future session cannot load `asreml`, check the `.Renviron` value first and then confirm the license token still exists.

## Quick verification command

Use this command from PowerShell to verify that the correct R and `asreml` package are available:

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' -e "cat(R.version.string, '\n'); cat(.libPaths(), sep='\n'); cat('\n'); library(asreml); cat(as.character(packageVersion('asreml')), '\n')"
```

Expected output should include:

- `R version 4.3.3`
- a user library path under `<your R user library>`
- `Loading ASReml-R version 4.2`

## Project library note

This project also uses a local package library under:

- `<project>/.r-lib`

Most local verification scripts prepend that path with:

```r
.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = TRUE), .libPaths()))
```

That keeps project dependencies stable while still allowing `asreml` to load from the user library.

## Reproducible wheat example

The wheat data file used during verification was:

- `<data dir>/WheatPhenoGenoOLD.Rdata`

To run the reusable example script:

```powershell
& 'C:\Program Files\R\R-4.3.3\bin\Rscript.exe' tools/run_wheat_asreml_examples.R '<data dir>/WheatPhenoGenoOLD.Rdata'
```

This runs:

- `CV2` with `grm`
- `CV2` with `grm + COP`
- true prediction with `grm`
- true prediction with `grm + COP`

Results are written as `.rds` files under:

- `tools/outputs/wheat_asreml_examples`

## Main code-path fixes applied

The verified fixes are in package source files used by `model_execute()`:

- `R/asreml_utilis_new.R`
- `R/random_terms_fit_new.R`

What changed:

- inverse kernels are now built correctly and returned as inverse triplets with `INVERSE = TRUE`
- the multi-kernel MET builder keeps the requested covariance structure on the first kernel and uses `idv(Env)` for additional kernels to avoid the previous two-kernel `corgh(Env)` failure

## Regression coverage

Regression tests added:

- `tests/testthat/test-asreml-multi-kernel-builder.R`

That test covers:

- inverse-triplet construction for `compute_inverse_and_sparse(..., inverse = TRUE)`
- multi-kernel MET random-term construction for `asreml`
