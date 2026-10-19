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

GPU and CPU
- Every deep-learning model trains on the GPU when PyTorch sees one: with
  `device = NULL` (the default) a CUDA GPU is used, else Apple MPS, else the
  CPU. `device = "cuda"` (or `"gpu"`), `"cuda:<index>"`, `"mps"` or `"cpu"`
  names one; any other value stops before the data are read. PredictProR
  installs the CUDA build of PyTorch when an NVIDIA driver is present.
- A fit that cannot use its GPU continues on the CPU, and PredictProR says so:
  - `device = "cuda"` was given, but PyTorch sees no CUDA GPU;
  - moving the model to the GPU failed;
  - training or prediction on the GPU hit an error (for example CUDA out of
    memory); training goes on from where it stopped, with the same optimizer
    state.
  If `torch.compile` fails on the GPU, the model trains uncompiled on the same
  GPU (also reported). A machine with an NVIDIA GPU that its PyTorch cannot
  use (a CPU-only build, or a driver older than the build) is reported too:
  uninstall torch from that Python and PredictProR installs the CUDA build on
  the next run.
- Each `model_execute()` run with such a fallback ends with one warning, and
  `Run_metadata` records, for the fits of the run (parallel workers
  included), `dl_device`, the fits per device (for example `cuda: 36 fits`),
  and `dl_device_notes`, `none` or each fallback with its model and number of
  fits.
- On a machine without a GPU, two models are very slow: TabAttention (SAINT)
  and TabTransformer (FT-Transformer). They turn every marker into a token and
  attend over all pairs of markers, so their cost grows with the number of
  markers and a CPU cannot keep up. The other models are as fast or faster on
  a CPU on small panels: a GPU pays a start-up cost that a small network does
  not earn back.
- Measured on one machine, one true-prediction call as in scenario 01 (80
  lines, 400 markers, `epochs = "auto"`, 3 bootstrap fits). CPU: Intel Core
  i7-12700H, 16 threads; GPU: NVIDIA RTX A2000 8 GB laptop; PyTorch 2.14.1,
  CPU and CUDA 13.0 builds. Seconds per call:

  | Model             | GPU | CPU                |
  |-------------------|----:|-------------------:|
  | TabAttention      |  59 | > 900 (stopped)    |
  | TabTransformer    |  74 | > 900 (stopped)    |
  | Conv1DNet         |  49 | 123                |
  | TabNet            |  57 | 26                 |
  | GPNet             |  42 | 23                 |
  | MixtureOfExperts  |  37 | 21                 |
  | ResNet            |  32 | 19                 |
  | DenseNeuralNet    |  25 | 17                 |
  | FactorNet         |  27 | 14                 |
  | CrossNet          |  27 | 16                 |
  | DenseAttentionNet |  25 | 15                 |
  | LightTreeNet      |  23 | 14                 |

- Without a GPU, leave TabAttention and TabTransformer out of `GS_model`, or
  fit them on fewer markers (LD pruning, feature selection). The scenario suite
  has 17 calls of these two models: about 23 minutes in total on the laptop
  GPU, hours on a CPU.

If a future session binds reticulate to the wrong Python
- Set `RETICULATE_PYTHON` to:
  the Python executable inside the `predictdl` env for that OS
- Restart R and load `PredictProR` again.

Verification helper
- Run:
  `Rscript tools/verify_all_deep_models_runtime.R`
- For a CUDA sweep:
  set `PREDICTPRO_DL_DEVICE=cuda` before running the script.
