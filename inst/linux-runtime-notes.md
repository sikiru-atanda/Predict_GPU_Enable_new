# Linux runtime notes

Validated on Ubuntu 22.04 (WSL2) with R 4.4.1, GCC 11, ASReml-R and a conda
Python 3.11 env, against PredictProR 0.20.170/0.20.171.

## System libraries

The package compiles with the standard toolchain (`r-base-dev`, gcc/g++,
gfortran, LAPACK/BLAS). One dependency needs an extra system library:

```sh
sudo apt-get install libglpk40   # needed by igraph (libglpk.so.40, libltdl.so.7)
```

Without it `igraph` installs but fails to load
(`libglpk.so.40: cannot open shared object file`), so PredictProR cannot load.
Prebuilt Ubuntu binaries from Posit Package Manager
(`https://packagemanager.posit.co/cran/__linux__/jammy/latest`) install the R
dependencies in a few minutes.

## Python

Python models run as subprocesses. A bare system `python3` usually has no
numpy/scikit-learn; PredictProR no longer auto-selects such an interpreter for
ML models and instead asks for `PREDICTPRO_PYTHON`. Point the env vars at a
conda/venv Python before running models:

```sh
export PREDICTPRO_PYTHON=$HOME/miniconda3/envs/<env>/bin/python      # ML
export PREDICTPRO_DL_PYTHON=$PREDICTPRO_PYTHON                       # DL (torch)
export PREDICTPRO_GP_PYTHON=$PREDICTPRO_PYTHON                       # GP (gpytorch, linear_operator, zarr)
```

## Output files

Exported file names are case-exact and identical on every OS
(`Predicted_Value.csv`, `variance_components.csv`, ...). Linux file systems are
case-sensitive, so scripts must use these exact names.

## Same results on Windows and Linux?

Twenty small workflows (`tarball_user_validation/scripts/regression_check_*.R`)
were run with the same seed on both systems:

| Engine | Windows vs Linux |
|---|---|
| BGLR Bayes (BayesB, GBLUP_BRR, RKHS, MET, hybrid, binary), ASReml (single, MET corgh/fa1, hybrid, CV), Bayes CV and feature-selection CV | bit-identical |
| RandomForest (sklearn) | bit-identical with the same scikit-learn/numpy versions; differs (cor 0.998) across scikit-learn versions |
| XGBoost | bit-identical with `subsample = 1` and `colsample_bytree = 1`. With row/column subsampling it differs (cor 0.995) even with the same XGBoost version: its sampling uses C++ standard-library distributions, which differ between MSVC and GCC builds |
| Deep learning (torch) | not bit-reproducible across torch versions or CPU vs GPU (cor 0.94 for a 10-epoch toy net); use `dl_n_seeds` ensembles for stable rankings |

Each system is reproducible run to run with a fixed `random_state`; the
differences above are only between systems or library versions.
