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

Python models run as subprocesses. PredictProR installs the Python packages
the requested models need on first use (numpy, pandas, scikit-learn,
xgboost/catboost/lightgbm, PyTorch, GPyTorch, ...), once per interpreter and
before any parallel worker starts. PyTorch comes from the PyTorch wheel index:
the newest CUDA build the NVIDIA driver supports, else the CPU build.

A bare system `python3` on Debian/Ubuntu cannot take installs: it is not
writable, and on recent releases it is "externally managed" (PEP 668) and
refuses pip. PredictProR then creates its own virtualenv from it under
`~/.cache/R/PredictProR/python-envs` (or `PREDICTPRO_PYTHON_ENV_HOME`). No sudo
is needed: without the `python3-venv` package the virtualenv is created without
pip, and PredictProR adds pip from its official zipapp
(`https://bootstrap.pypa.io/pip/pip.pyz`; `PREDICTPRO_PIP_BOOTSTRAP_URL` sets
another URL). The same happens for an interpreter you cannot write to.
Checked on Ubuntu 22.04 with only the system `python3` 3.10 (no pip, no
`python3-venv`) and nothing configured.

A virtualenv's `bin/python` is a symlink to the base Python. PredictProR keeps
the path as given (before 0.20.196 it resolved the link and ran the base
Python, outside the virtualenv).

Without these variables, PredictProR checks the interpreters it finds
(virtualenvs, conda environments, system Pythons) to choose one per model
family. Each check starts the interpreter and imports torch, so the results
are kept in `~/.cache/R/PredictProR/python-discovery.rds` (on every OS the
folder is `tools::R_user_dir("PredictProR", "cache")`). A result is used
again only while the interpreter and its site-packages folders are
unchanged (installing, upgrading or removing a package counts as a change),
and for at most 30 days; the conda environment list is kept until conda
creates or removes an environment, for at most 7 days. The same file keeps
the memory of an idle parallel worker, which the parallel policy measures
by starting one: until PredictProR or one of its imports is reinstalled, R
or the ML interpreter changes, for at most 30 days. Set
`PREDICTPRO_PROBE_CACHE=false` to check and measure every time.

To use a specific interpreter, point the env vars at it:

```sh
export PREDICTPRO_PYTHON=$HOME/miniconda3/envs/<env>/bin/python      # ML
export PREDICTPRO_DL_PYTHON=$PREDICTPRO_PYTHON                       # DL (torch)
export PREDICTPRO_GP_PYTHON=$PREDICTPRO_PYTHON                       # GP (gpytorch, linear_operator, zarr)
```

To prepare a server before running jobs, run
`PredictProR::ensure_python_dependencies()` once. On a machine without internet
access, point pip at a mirror (`PIP_INDEX_URL`, and `PYTORCH_INDEX_URL` for
PyTorch wheels) or install the packages yourself and set
`PREDICTPRO_AUTO_INSTALL=false`.

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
