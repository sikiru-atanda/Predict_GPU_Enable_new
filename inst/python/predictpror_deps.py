"""Python packages for PredictProR's ML, DL and GP bridges, installed on first use.

PredictProR installs the packages its models need into the interpreter that
runs them. This module is the one place that does it:

* the check is import-free (``importlib.util.find_spec``), so a bridge whose
  packages are present starts without delay;
* an install holds a lock on the environment (an OS file lock, released when
  the process ends, even if it crashes), so parallel workers never run pip on
  the same environment at once; the missing set is re-checked once the lock
  is held, so a worker that waited finds the packages and installs nothing;
* everything missing goes into one pip call, and the import caches are
  cleared afterwards, so the packages work in the process that installed them;
* PyTorch comes from the PyTorch wheel index: the newest CUDA build the NVIDIA
  driver supports (then older ones), else the CPU build. Packages that need
  torch are installed after it, with the installed torch and numpy pinned,
  so pip neither swaps the CUDA build for a CPU one nor moves numpy;
* a failure raises DependencyError with the reason, pip's output and what to do.

R runs ``python predictpror_deps.py ensure --groups ...`` once in the main
session before it starts parallel workers. Each bridge calls ensure_groups()
again when it starts, which only costs the find_spec checks.

Environment variables: PREDICTPRO_AUTO_INSTALL=false turns installing off
(missing packages then raise DependencyError); PYTORCH_INDEX_URL forces the
PyTorch wheel index; PREDICTPRO_TORCH_DEVICE=cpu installs the CPU build;
PREDICTPRO_PIP_LOCK_TIMEOUT (seconds, default 3600) bounds the wait for
another process's install; PREDICTPRO_PIP_BOOTSTRAP_URL replaces the URL of
pip's zipapp, used when an interpreter has no pip and no ensurepip. pip's own
settings (PIP_INDEX_URL, proxies, ...) apply as usual.
"""
from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import os
import platform
import re
import subprocess
import sys
import sysconfig
import tempfile
import time
from contextlib import contextmanager
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# group -> [(import name, pip requirement)]
GROUPS: Dict[str, List[Tuple[str, str]]] = {
    "numpy": [("numpy", "numpy>=1.24")],
    "tabular": [("numpy", "numpy>=1.24"), ("pandas", "pandas>=1.5")],
    "sklearn": [("scipy", "scipy"), ("sklearn", "scikit-learn")],
    "xgboost": [("xgboost", "xgboost")],
    "catboost": [("catboost", "catboost")],
    "lightgbm": [("lightgbm", "lightgbm")],
    "torch": [("torch", "torch")],
    "gpytorch": [("gpytorch", "gpytorch"), ("linear_operator", "linear-operator")],
    "gp": [("scipy", "scipy"), ("gpytorch", "gpytorch"), ("linear_operator", "linear-operator"),
           ("zarr", "zarr")],
}

# top-level module -> group that provides it (for ensure_module)
MODULE_GROUP: Dict[str, str] = {
    "numpy": "numpy", "pandas": "tabular", "scipy": "sklearn", "sklearn": "sklearn",
    "xgboost": "xgboost", "catboost": "catboost", "lightgbm": "lightgbm", "torch": "torch",
    "gpytorch": "gpytorch", "linear_operator": "gpytorch", "zarr": "gp",
}

# pip requirements that depend on torch: installed after torch, with torch pinned
_NEEDS_TORCH = {"gpytorch", "linear-operator"}

TORCH_INDEX = "https://download.pytorch.org/whl/"
# CUDA wheel tags and the CUDA version each needs from the driver, newest first
_CUDA_TAGS: List[Tuple[str, Tuple[int, int]]] = [
    ("cu130", (13, 0)), ("cu129", (12, 9)), ("cu128", (12, 8)), ("cu126", (12, 6)),
    ("cu124", (12, 4)), ("cu121", (12, 1)), ("cu118", (11, 8)),
]
_MIN_PYTHON = (3, 10)


class DependencyError(RuntimeError):
    """A required package is missing and could not be installed."""

    def __init__(self, reason: str, message: str):
        super().__init__(message)
        self.reason = reason


def _say(msg: str) -> None:
    print(f"PredictProR: {msg}", file=sys.stderr, flush=True)


def _truthy(value: Optional[str], default: bool) -> bool:
    if value is None or not value.strip():
        return default
    return value.strip().lower() not in ("0", "false", "no", "off")


def auto_install_enabled() -> bool:
    return _truthy(os.environ.get("PREDICTPRO_AUTO_INSTALL"), True)


def _requirements(groups: Iterable[str]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    seen = set()
    for group in groups:
        if group not in GROUPS:
            raise ValueError(f"Unknown dependency group '{group}'. Known: {', '.join(sorted(GROUPS))}")
        for module, spec in GROUPS[group]:
            if module not in seen:
                seen.add(module)
                out.append((module, spec))
    return out


def _present(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def _torch_probe() -> Tuple[bool, Dict[str, object]]:
    """Import torch in a child process (a broken torch cannot be unloaded here)."""
    code = (
        "import json, torch\n"
        "torch.tensor([1.0]).sum().item()\n"
        "print(json.dumps({'version': torch.__version__, 'cuda': torch.version.cuda,"
        " 'cuda_available': bool(torch.cuda.is_available())}))\n"
    )
    try:
        proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                              timeout=600, env=_child_env())
    except Exception as exc:  # pragma: no cover - interpreter could not start
        return False, {"error": repr(exc)}
    if proc.returncode == 0:
        try:
            return True, json.loads(proc.stdout.strip().splitlines()[-1])
        except Exception:
            return True, {}
    tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-5:]
    return False, {"error": "\n".join(tail)}


def missing(groups: Iterable[str], probe_torch: bool = False) -> List[Tuple[str, str]]:
    """Requirements of ``groups`` that are not installed (torch also when it cannot import)."""
    out = []
    for module, spec in _requirements(groups):
        if not _present(module):
            out.append((module, spec))
        elif module == "torch" and probe_torch and not _torch_probe()[0]:
            out.append((module, spec))
    return out


# ---------------------------------------------------------------------------
# Can this interpreter take installs?
# ---------------------------------------------------------------------------
def _dir_writable(path: str) -> bool:
    # Not tempfile: on Windows it takes "access denied" for a name clash while
    # os.access() says writable (os.access ignores ACLs), and retries ~2**31 times.
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        return False
    for attempt in range(5):
        probe = os.path.join(path, f".predictpror-write-test-{os.getpid()}-{time.time_ns()}-{attempt}")
        try:
            fd = os.open(probe, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        except OSError:
            return False
        os.close(fd)
        try:
            os.unlink(probe)
        except OSError:
            pass
        return True
    return False


def _externally_managed() -> bool:
    # PEP 668: a distribution's system Python refuses pip outside a virtualenv
    if sys.prefix != getattr(sys, "base_prefix", sys.prefix):
        return False
    if _truthy(os.environ.get("PIP_BREAK_SYSTEM_PACKAGES"), False):
        return False
    stdlib = sysconfig.get_path("stdlib") or ""
    return bool(stdlib) and os.path.exists(os.path.join(stdlib, "EXTERNALLY-MANAGED"))


def install_problem() -> Optional[str]:
    """Why packages cannot be installed into this interpreter, or None."""
    if sys.version_info[:2] < _MIN_PYTHON:
        return "python_too_old"
    if _externally_managed():
        return "externally_managed"
    purelib = sysconfig.get_paths()["purelib"]
    if not _dir_writable(purelib) or not _dir_writable(sys.prefix):
        return "not_writable"
    return None


def _child_env() -> Dict[str, str]:
    env = dict(os.environ)
    # install into (and check against) this environment only, never user site
    env["PYTHONNOUSERSITE"] = "1"
    env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
    env.setdefault("PIP_NO_INPUT", "1")
    return env


_PIP_ZIPAPP_URL = "https://bootstrap.pypa.io/pip/pip.pyz"


def _pip_works() -> bool:
    probe = subprocess.run([sys.executable, "-m", "pip", "--version"], capture_output=True,
                           text=True, env=_child_env())
    return probe.returncode == 0


def _pip_from_zipapp() -> str:
    # pip's official single-file build: it runs on any supported Python and
    # installs pip into the interpreter that runs it. Needed where ensurepip
    # is left out (Debian/Ubuntu without the python3-venv package).
    import urllib.request
    url = os.environ.get("PREDICTPRO_PIP_BOOTSTRAP_URL", _PIP_ZIPAPP_URL)
    target = os.path.join(sys.prefix, f".predictpror-pip-{os.getpid()}.pyz")
    try:
        with urllib.request.urlopen(url, timeout=120) as resp, open(target, "wb") as out:  # honours HTTPS_PROXY
            out.write(resp.read())
    except Exception as exc:
        return f"could not download {url}: {exc}\n"
    try:
        run = subprocess.run([sys.executable, target, "install", "pip"], capture_output=True,
                             text=True, env=_child_env())
        return (run.stdout or "") + (run.stderr or "")
    finally:
        try:
            os.unlink(target)
        except OSError:
            pass


def _ensure_pip() -> None:
    if _pip_works():
        return
    _say("pip is missing; bootstrapping it with ensurepip ...")
    boot = subprocess.run([sys.executable, "-m", "ensurepip", "--upgrade", "--default-pip"],
                          capture_output=True, text=True, env=_child_env())
    output = (boot.stdout or "") + (boot.stderr or "")
    if not _pip_works():
        _say("ensurepip is not available; installing pip from its zipapp ...")
        output += _pip_from_zipapp()
    if not _pip_works():
        raise DependencyError("pip_missing", _message("pip_missing", [], output))


# ---------------------------------------------------------------------------
# Lock on the environment
# ---------------------------------------------------------------------------
def _lock(fh) -> None:
    if os.name == "nt":
        import msvcrt
        fh.seek(0)
        msvcrt.locking(fh.fileno(), msvcrt.LK_NBLCK, 1)
    else:
        import fcntl
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def _unlock(fh) -> None:
    if os.name == "nt":
        import msvcrt
        fh.seek(0)
        msvcrt.locking(fh.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        import fcntl
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


@contextmanager
def environment_lock(timeout: Optional[float] = None):
    """Hold the install lock of this environment (released if the process dies)."""
    if timeout is None:
        try:
            timeout = float(os.environ.get("PREDICTPRO_PIP_LOCK_TIMEOUT", "3600"))
        except ValueError:
            timeout = 3600.0
    path = os.path.join(sys.prefix, ".predictpror-install.lock")
    fh = open(path, "a+")
    deadline = time.monotonic() + timeout
    waited = False
    try:
        while True:
            try:
                _lock(fh)
                break
            except OSError:
                if not waited:
                    _say("another PredictProR process is installing packages into this "
                         "environment; waiting for it ...")
                    waited = True
                if time.monotonic() > deadline:
                    raise DependencyError("lock_timeout", (
                        f"Waited {int(timeout)} s for another process installing packages into "
                        f"{sys.prefix}. If no such process is running, delete {path}; "
                        "PREDICTPRO_PIP_LOCK_TIMEOUT sets the wait."))
                time.sleep(2.0)
        yield
    finally:
        try:
            _unlock(fh)
        except OSError:
            pass
        fh.close()


# ---------------------------------------------------------------------------
# pip
# ---------------------------------------------------------------------------
def _pip(args: Sequence[str], verbose: bool) -> Tuple[int, List[str]]:
    cmd = [sys.executable, "-m", "pip", *args]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            errors="replace", env=_child_env())
    lines: List[str] = []
    assert proc.stdout is not None
    for raw in proc.stdout:
        line = raw.rstrip("\r\n")
        lines.append(line)
        if verbose and line.strip():
            print(f"  {line}", file=sys.stderr, flush=True)
    return proc.wait(), lines


def _pip_install(specs: Sequence[str], verbose: bool, index_url: Optional[str] = None,
                 constraints: Sequence[str] = ()) -> Tuple[int, List[str]]:
    args = ["install", "--progress-bar", "off", "--upgrade-strategy", "only-if-needed"]
    if index_url:
        args += ["--index-url", index_url]
    cfile = None
    if constraints:
        fd, cfile = tempfile.mkstemp(prefix="predictpror-constraints-", suffix=".txt")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write("\n".join(constraints) + "\n")
        args += ["--constraint", cfile]
    try:
        return _pip(args + list(specs), verbose)
    finally:
        if cfile:
            try:
                os.unlink(cfile)
            except OSError:
                pass


def _installed_version(dist: str) -> Optional[str]:
    try:
        from importlib import metadata
        return metadata.version(dist)
    except Exception:
        return None


_NETWORK_SIGNS = (
    "Could not fetch URL", "NewConnectionError", "Failed to establish a new connection",
    "ProxyError", "Max retries exceeded", "getaddrinfo failed", "Name or service not known",
    "Temporary failure in name resolution", "SSLError", "ConnectTimeoutError", "ReadTimeoutError",
    "Network is unreachable", "No connection could be made",
)


def _index_reachable(url: str, timeout: float = 15.0) -> bool:
    import urllib.request
    try:
        with urllib.request.urlopen(url, timeout=timeout):  # honours HTTP(S)_PROXY
            return True
    except Exception:
        return False


def _classify_pip_failure(lines: Sequence[str], index_url: Optional[str] = None) -> str:
    text = "\n".join(lines)
    if any(sign in text for sign in _NETWORK_SIGNS):
        return "network"
    if "ResolutionImpossible" in text or "conflicting dependencies" in text:
        return "conflict"
    if "No matching distribution found" in text or "Could not find a version" in text:
        # "(from versions: none)" is also all pip says when the index cannot
        # be reached (no retries configured), so look at the index itself
        if "from versions: none" in text and "PIP_NO_INDEX" not in os.environ:
            m = re.search(r"Looking in indexes:\s*(\S+)", text)
            url = index_url or (m.group(1).rstrip(",") if m else "https://pypi.org/simple/")
            if not _index_reachable(url):
                return "network"
        return "no_distribution"
    if "Permission denied" in text or "Access is denied" in text or "WinError 5" in text:
        return "not_writable"
    return "pip_failed"


_ADVICE = {
    "network": ("pip could not reach the package index. Check the internet connection or the "
                "proxy (HTTPS_PROXY), or point pip to a reachable mirror with PIP_INDEX_URL "
                "(PYTORCH_INDEX_URL for PyTorch wheels)."),
    "conflict": ("the new packages need versions that conflict with packages already in this "
                 "environment. Use a separate environment for PredictProR "
                 "(setup_predictdl_env(), setup_predictgp_env() or PREDICTPRO_*_PYTHON)."),
    "no_distribution": (f"no wheel was found for this platform/Python ({platform.system()} "
                        f"{platform.machine()}, Python {sys.version_info[0]}.{sys.version_info[1]}). "
                        "Use a supported Python (3.10-3.13) or set PIP_INDEX_URL/PYTORCH_INDEX_URL."),
    "not_writable": ("this environment cannot be written to. From R, PredictProR then creates its "
                     "own environment; otherwise set PREDICTPRO_*_PYTHON to a Python you can "
                     "install into."),
    "externally_managed": ("this is a system Python managed by the operating system (PEP 668), "
                           "which refuses pip. From R, PredictProR then creates its own "
                           "environment; otherwise set PREDICTPRO_*_PYTHON to a virtualenv/conda Python."),
    "python_too_old": (f"Python {_MIN_PYTHON[0]}.{_MIN_PYTHON[1]} or newer is required; set "
                       "PREDICTPRO_*_PYTHON to a newer interpreter."),
    "pip_missing": ("this Python has no pip, ensurepip could not add it and pip's zipapp could not "
                    "be downloaded (bootstrap.pypa.io; PREDICTPRO_PIP_BOOTSTRAP_URL sets another "
                    "URL). Install pip for it (e.g. the python3-pip / python3-venv system package) "
                    "or use another Python."),
    "auto_install_disabled": ("PREDICTPRO_AUTO_INSTALL is false. Install the packages yourself, "
                              "or unset it so PredictProR installs them."),
    "not_importable": "pip reported success but the packages still cannot be imported.",
    "pip_failed": "pip failed; its output is below.",
}


def _message(reason: str, specs: Sequence[str], output: str = "") -> str:
    head = (f"PredictProR needs Python packages that are not installed in {sys.executable}"
            + (f": {', '.join(specs)}" if specs else "") + ".")
    msg = f"{head}\nThey could not be installed: {_ADVICE.get(reason, _ADVICE['pip_failed'])}"
    tail = [ln for ln in output.splitlines() if ln.strip()][-40:]
    if tail:
        msg += "\n--- pip output (last lines) ---\n" + "\n".join(tail)
    return msg


# ---------------------------------------------------------------------------
# PyTorch
# ---------------------------------------------------------------------------
def _driver_cuda_version() -> Optional[Tuple[int, int]]:
    try:
        proc = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=30)
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    m = re.search(r"CUDA Version:\s*(\d+)\.(\d+)", proc.stdout or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def torch_candidates(prefer_gpu: bool = True, cuda: str = "auto",
                     index_url: Optional[str] = None) -> List[Tuple[str, Optional[str]]]:
    """(label, wheel index) to try in order; None is pip's default index."""
    forced = (index_url or os.environ.get("PYTORCH_INDEX_URL") or "").strip()
    if forced:
        return [("forced", forced)]
    if platform.system() == "Darwin":
        return [("default", None)]  # the PyPI wheels carry MPS support
    cpu = [("cpu", TORCH_INDEX + "cpu")]
    cuda = (cuda or "auto").strip().lower()
    if (not prefer_gpu or cuda == "cpu"
            or os.environ.get("PREDICTPRO_TORCH_DEVICE", "").strip().lower() == "cpu"):
        return cpu
    if cuda.startswith("cu"):
        return [(cuda, TORCH_INDEX + cuda)] + cpu
    driver = _driver_cuda_version()
    if driver is None:
        return cpu
    tags = [tag for tag, need in _CUDA_TAGS if need <= driver][:3]
    return [(tag, TORCH_INDEX + tag) for tag in tags] + cpu


def install_torch(verbose: bool = True, prefer_gpu: bool = True, cuda: str = "auto",
                  torch_version: Optional[str] = None,
                  index_url: Optional[str] = None) -> Dict[str, object]:
    """Install a working torch (caller holds the environment lock)."""
    spec = f"torch=={torch_version}" if torch_version else "torch"
    if _present("torch"):
        _say("torch is installed but cannot be imported; reinstalling it ...")
        _pip(["uninstall", "-y", "torch"], verbose)
    log: List[str] = []
    reason = "pip_failed"
    for label, url in torch_candidates(prefer_gpu, cuda, index_url):
        where = url or "the default index"
        _say(f"installing PyTorch ({label} build from {where}); this is a large download ...")
        rc, lines = _pip_install([spec], verbose, index_url=url)
        log += [f"[{label}] {ln}" for ln in lines]
        if rc != 0:
            reason = _classify_pip_failure(lines, url)
            if reason in ("network", "not_writable"):
                break  # another wheel index will not help
            continue
        importlib.invalidate_caches()
        ok, info = _torch_probe()
        if ok:
            info = dict(info)
            info.update({"build": label, "index": where})
            return info
        log.append(f"[{label}] installed but import failed: {info.get('error', '')}")
        reason = "not_importable"
        _pip(["uninstall", "-y", "torch"], verbose)
    raise DependencyError(reason, _message(reason, [spec], "\n".join(log)))


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------
def ensure_groups(groups: Iterable[str], verbose: bool = True, prefer_gpu: bool = True,
                  cuda: str = "auto", torch_version: Optional[str] = None,
                  index_url: Optional[str] = None, probe_torch: bool = False) -> Dict[str, object]:
    """Install whatever ``groups`` need that is missing; return what was done."""
    groups = list(groups)
    report: Dict[str, object] = {"python": sys.executable, "installed": [], "torch": None}
    todo = missing(groups, probe_torch=probe_torch)
    if not todo:
        return report
    specs = [spec for _, spec in todo]
    if not auto_install_enabled():
        raise DependencyError("auto_install_disabled", _message("auto_install_disabled", specs))
    problem = install_problem()
    if problem:
        raise DependencyError(problem, _message(problem, specs))
    with environment_lock():
        # another process may have installed them while this one waited
        todo = missing(groups, probe_torch=probe_torch)
        if not todo:
            return report
        _ensure_pip()
        specs = [spec for _, spec in todo]
        _say(f"installing {', '.join(specs)} into {sys.executable} ...")
        torch_needed = any(module == "torch" for module, _ in todo)
        plain = [spec for module, spec in todo if module != "torch" and spec not in _NEEDS_TORCH]
        after_torch = [spec for module, spec in todo if spec in _NEEDS_TORCH]
        if plain:
            # keep the numpy that compiled packages here were built against
            pins = [f"numpy=={v}" for v in [_installed_version("numpy")] if v]
            rc, lines = _pip_install(plain, verbose, constraints=pins)
            if rc != 0:
                reason = _classify_pip_failure(lines)
                raise DependencyError(reason, _message(reason, plain, "\n".join(lines)))
        if torch_needed:
            report["torch"] = install_torch(verbose, prefer_gpu, cuda, torch_version, index_url)
        if after_torch:
            pins = [f"{dist}=={v}" for dist in ("torch", "numpy")
                    for v in [_installed_version(dist)] if v]
            rc, lines = _pip_install(after_torch, verbose, constraints=pins)
            if rc != 0:
                reason = _classify_pip_failure(lines)
                raise DependencyError(reason, _message(reason, after_torch, "\n".join(lines)))
        importlib.invalidate_caches()
        purelib = sysconfig.get_paths()["purelib"]
        if purelib not in sys.path:
            sys.path.append(purelib)
        still = [spec for module, spec in todo if not _present(module)]
        if still:
            raise DependencyError("not_importable", _message("not_importable", still))
        report["installed"] = specs
        _say("installed " + ", ".join(specs) + ".")
    return report


def ensure_module(name: str, verbose: bool = True):
    """Import ``name``, first installing the group that provides it if it is missing."""
    top = name.split(".", 1)[0]
    if not _present(top) and top in MODULE_GROUP:
        ensure_groups([MODULE_GROUP[top]], verbose=verbose)
    return importlib.import_module(name)


def _check(groups: Sequence[str]) -> Dict[str, object]:
    todo = missing(groups, probe_torch=True)
    return {
        "python": sys.executable,
        "python_version": platform.python_version(),
        "missing": [spec for _, spec in todo],
        "missing_modules": [module for module, _ in todo],
        "install_problem": install_problem() if todo else None,
        "auto_install": auto_install_enabled(),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="predictpror_deps")
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name in ("check", "ensure"):
        p = sub.add_parser(name)
        p.add_argument("--groups", required=True, help="comma-separated dependency groups")
        p.add_argument("--json-out", default=None)
        p.add_argument("--prefer-gpu", action=argparse.BooleanOptionalAction, default=True)
        p.add_argument("--cuda", default="auto")
        p.add_argument("--torch-version", default=None)
        p.add_argument("--index-url", default=None)
        p.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    groups = [g.strip() for g in args.groups.split(",") if g.strip()]
    rc = 0
    try:
        if args.cmd == "check":
            out: Dict[str, object] = _check(groups)
        else:
            out = ensure_groups(groups, verbose=not args.quiet, prefer_gpu=args.prefer_gpu,
                                cuda=args.cuda, torch_version=args.torch_version,
                                index_url=args.index_url, probe_torch=True)
        out["status"] = "ok"
    except DependencyError as exc:
        out = {"status": "error", "reason": exc.reason, "message": str(exc), "python": sys.executable}
        print(str(exc), file=sys.stderr, flush=True)
        rc = 3
    except ValueError as exc:
        out = {"status": "error", "reason": "bad_request", "message": str(exc), "python": sys.executable}
        print(str(exc), file=sys.stderr, flush=True)
        rc = 2
    text = json.dumps(out)
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            fh.write(text)
    else:
        print(text)
    return rc


if __name__ == "__main__":
    sys.exit(main())
