"""Canonical import shim for the large-scale operator backend.

The implementation file keeps its historical dated name, but runtime code
imports ``mixed_model_gpu_large_scale_backend``. Keeping this small wrapper in
the repo makes the import work outside pytest/benchmark scripts without
manual ``sys.modules`` aliasing.
"""
from __future__ import annotations

from gp_large_backend import *  # noqa: F401,F403
from gp_large_backend import _scatter_add_cells, _unique_inverse  # noqa: F401
