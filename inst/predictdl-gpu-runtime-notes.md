PredictProR deep-learning runtime notes

Preferred reticulate Python
- Preferred env name: `predictdl`
- Windows path pattern:
  `~/.virtualenvs/predictdl/Scripts/python.exe`
- macOS/Linux path pattern:
  `~/.virtualenvs/predictdl/bin/python`

Current state
- The package now prefers the `predictdl` interpreter automatically on load when `RETICULATE_PYTHON` is not already set.
- Runtime selection is OS-aware and checks both:
  standard virtualenv layouts on Windows, macOS, and Linux
  conda environments with the matching name `predictdl`
- When that environment is CUDA-enabled, the package can use `device = "cuda"`.
- End-to-end deep-learning model-family verification should be run against the active `predictdl` environment on each machine.

If a future session binds reticulate to the wrong Python
- Set `RETICULATE_PYTHON` to:
  the Python executable inside the `predictdl` env for that OS
- Restart R and load `PredictProR` again.

Verification helper
- Run:
  `Rscript tools/verify_all_deep_models_runtime.R`
- For a CUDA sweep:
  set `PREDICTPRO_DL_DEVICE=cuda` before running the script.
