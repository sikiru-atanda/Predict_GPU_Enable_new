from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import numpy as np
except Exception:
    np = None

try:
    import pandas as pd
except Exception:
    pd = None


def _module_dir() -> Path:
    return Path(__file__).resolve().parent


def _load_dl_models():
    mod_path = _module_dir() / "dl_models.py"
    if not mod_path.exists():
        raise FileNotFoundError(f"dl_models.py not found: {mod_path}")
    spec = importlib.util.spec_from_file_location("predictpror_dl_models", mod_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _require_tabular_runtime():
    global np, pd
    if np is None:
        import numpy as _np
        np = _np
    if pd is None:
        import pandas as _pd
        pd = _pd


def _read_matrix_csv(path: str) -> np.ndarray:
    _require_tabular_runtime()
    return pd.read_csv(path).to_numpy(dtype=float)


def _read_target_csv(path: str, task: str):
    _require_tabular_runtime()
    vals = pd.read_csv(path).to_numpy()
    if vals.ndim == 2 and vals.shape[1] == 1:
        vals = vals[:, 0]
    if task in ("gaussian", "regression", "binary", "multitask_regression"):
        return vals.astype(float)
    if vals.ndim == 2:
        vals = vals[:, 0]
    return vals.astype(str)


def _read_json(path: Optional[str], default: Any):
    if path is None:
        return default
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _json_default(value: Any):
    _require_tabular_runtime()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(obj: Dict[str, Any], path: str):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default), encoding="utf-8")


def _read_class_levels(path: Optional[str]):
    levels = _read_json(path, [])
    if levels is None:
        return []
    return [str(x) for x in levels]


def _apply_class_levels(result: Dict[str, Any], task: str, class_levels):
    _require_tabular_runtime()
    if task not in ("binary", "multiclass") or not class_levels:
        return result
    prob = result.get("probabilities")
    if prob is not None:
        prob_arr = np.asarray(prob)
        if prob_arr.ndim == 1:
            prob_arr = prob_arr.reshape(-1, 1)
        result_classes = [str(x) for x in (result.get("classes") or [])]
        if len(class_levels) == prob_arr.shape[1]:
            result["classes"] = [str(x) for x in class_levels]
            if task == "multiclass":
                pred = np.asarray(result.get("predictions", []))
                try:
                    pred_idx = pred.astype(int).reshape(-1)
                    if pred_idx.size == prob_arr.shape[0] and np.all((0 <= pred_idx) & (pred_idx < len(class_levels))):
                        result["predictions"] = np.asarray([class_levels[i] for i in pred_idx])
                except (TypeError, ValueError):
                    pass
        elif result_classes and len(result_classes) == prob_arr.shape[1]:
            full_prob = np.zeros((prob_arr.shape[0], len(class_levels)), dtype=float)
            for j, cls in enumerate(result_classes):
                target_idx = None
                if cls in class_levels:
                    target_idx = class_levels.index(cls)
                else:
                    try:
                        encoded_idx = int(cls)
                        if 0 <= encoded_idx < len(class_levels):
                            target_idx = encoded_idx
                    except (TypeError, ValueError):
                        target_idx = None
                if target_idx is not None:
                    full_prob[:, target_idx] = prob_arr[:, j].astype(float)
            result["probabilities"] = full_prob
            result["classes"] = [str(x) for x in class_levels]
            if task == "multiclass":
                pred_idx = np.argmax(full_prob, axis=1).astype(int)
                result["predictions"] = np.asarray([class_levels[i] for i in pred_idx])
    elif task == "binary" and len(class_levels) == 2:
        result["classes"] = [str(x) for x in class_levels]
    return result


def _write_predictions(result: Dict[str, Any], task: str, out_dir: Path):
    _require_tabular_runtime()
    preds = np.asarray(result["predictions"])
    if task in ("gaussian", "regression", "binary", "multitask_regression"):
        pred_arr = preds.astype(float)
        if pred_arr.ndim == 1:
            pred_df = pd.DataFrame({"prediction": pred_arr})
        else:
            pred_df = pd.DataFrame(
                pred_arr,
                columns=[f"prediction_{i + 1}" for i in range(pred_arr.shape[1])],
            )
    else:
        pred_df = pd.DataFrame({"prediction": preds.astype(str)})
    pred_df.to_csv(out_dir / "predictions.csv", index=False)

    prob = result.get("probabilities")
    if prob is not None:
        prob_arr = np.asarray(prob, dtype=float)
        if prob_arr.ndim == 1:
            prob_arr = prob_arr.reshape(-1, 1)
        prob_df = pd.DataFrame(prob_arr)
        classes = result.get("classes")
        if classes is not None and len(classes) == prob_df.shape[1]:
            prob_df.columns = [str(x) for x in classes]
        prob_df.to_csv(out_dir / "probabilities.csv", index=False)

    meta = {
        "classes": [str(x) for x in (result.get("classes") or [])],
        "task": task,
        "history": result.get("history"),
        "python": sys.executable,
    }
    _write_json(meta, str(out_dir / "meta.json"))


def _fit_predict(args) -> int:
    params = _read_json(args.params_json, {})
    task = str(args.task).strip().lower()
    class_levels = _read_class_levels(args.class_levels_json)
    dl_models = _load_dl_models()
    result = _fit_predict_payload_with_module(
        dl_models=dl_models,
        model_type=str(args.model_type).lower(),
        X_train=_read_matrix_csv(args.x_train_csv),
        y_train=_read_target_csv(args.y_train_csv, task),
        X_test=_read_matrix_csv(args.x_test_csv),
        params=params,
        task=task,
    )
    result = _apply_class_levels(result, task, class_levels)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_predictions(result, task, out_dir)
    return 0


def _fit_predict_payload_with_module(
    dl_models,
    model_type: str,
    X_train,
    y_train,
    X_test,
    params: Optional[Dict[str, Any]] = None,
    task: str = "gaussian",
):
    return dl_models.fit_predict_payload(
        model_type=str(model_type).lower(),
        X_train=np.asarray(X_train),
        y_train=np.asarray(y_train),
        X_test=np.asarray(X_test),
        params=params or {},
        task=str(task).strip().lower(),
    )


def _batch_fit_predict(args) -> int:
    manifest = _read_json(args.manifest_json, {})
    jobs = manifest.get("jobs", manifest if isinstance(manifest, list) else [])
    if not isinstance(jobs, list):
        raise ValueError("Batch manifest must contain a list under 'jobs'.")

    dl_models = _load_dl_models()
    completed = []
    for i, job in enumerate(jobs):
        task = str(job["task"]).strip().lower()
        params = _read_json(job.get("params_json"), {})
        class_levels = _read_class_levels(job.get("class_levels_json"))
        result = _fit_predict_payload_with_module(
            dl_models=dl_models,
            model_type=str(job["model_type"]).strip().lower(),
            X_train=_read_matrix_csv(job["x_train_csv"]),
            y_train=_read_target_csv(job["y_train_csv"], task),
            X_test=_read_matrix_csv(job["x_test_csv"]),
            params=params,
            task=task,
        )
        result = _apply_class_levels(result, task, class_levels)
        out_dir = Path(job["out_dir"])
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_predictions(result, task, out_dir)
        completed.append({"id": str(job.get("id", i + 1)), "out_dir": str(out_dir)})

    summary_path = manifest.get("summary_json") if isinstance(manifest, dict) else None
    if summary_path:
        _write_json({"jobs": completed, "python": sys.executable}, summary_path)
    return 0


def _bootstrap_predict(args) -> int:
    params = _read_json(args.params_json, {})
    training_seeds = _read_json(args.training_seeds_json, None)
    task = str(args.task).strip().lower()
    class_levels = _read_class_levels(args.class_levels_json)
    dl_models = _load_dl_models()
    payload = dl_models.bootstrap_fit_predict_payload(
        model_type=str(args.model_type).lower(),
        X_train=_read_matrix_csv(args.x_train_csv),
        y_train=_read_target_csv(args.y_train_csv, task),
        X_pred=_read_matrix_csv(args.x_pred_csv),
        n_bootstrap=int(args.n_bootstrap),
        seed=int(args.seed),
        params=params,
        task=task,
        training_seeds=training_seeds,
        seed_aggregation=str(args.seed_aggregation),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(payload["bootstrap"]).to_csv(out_dir / "bootstrap.csv", index=False)
    pd.DataFrame(payload["seed_predictions"]).to_csv(
        out_dir / "seed_predictions.csv", index=False
    )
    pd.DataFrame({
        "seed_variance": payload["seed_variance"],
        "seed_range": payload["seed_range"],
    }).to_csv(out_dir / "seed_variability.csv", index=False)
    meta = {
        "task": payload["task"],
        "n_bootstrap": int(payload["n_bootstrap"]),
        "n_prediction_rows": int(payload["n_prediction_rows"]),
        "classes": class_levels or payload["classes"],
        "python": sys.executable,
        "training_seeds": payload["training_seeds"],
        "seed_aggregation": payload["seed_aggregation"],
        "n_model_fits": int(payload["n_model_fits"]),
    }
    _write_json(meta, str(out_dir / "meta.json"))
    return 0


def _prewarm(args) -> int:
    _load_dl_models()
    out_dir = Path(args.out_dir) if args.out_dir else None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_json({"python": sys.executable, "status": "ok"}, str(out_dir / "meta.json"))
    return 0


def _print_setup_info(info: Dict[str, Any]) -> int:
    info = dict(info or {})
    info["python"] = sys.executable
    for key, value in info.items():
        if isinstance(value, (list, tuple)):
            value = ",".join(str(x) for x in value)
        print(f"{key}={value}")
    return 0


def _setup_deps(args) -> int:
    mod = _load_dl_models()
    info = mod.setup_deps(
        prefer_gpu=bool(args.prefer_gpu),
        cuda=args.cuda,
        torch_version=args.torch_version,
        numpy_spec="numpy>=1.24",
        extra_packages=("torchvision", "torchaudio"),
        index_url=args.index_url,
    )
    return _print_setup_info(info)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_setup = sub.add_parser("setup-deps")
    p_setup.add_argument("--prefer-gpu", action=argparse.BooleanOptionalAction, default=True)
    p_setup.add_argument("--cuda", default="auto")
    p_setup.add_argument("--torch-version", default=None)
    p_setup.add_argument("--index-url", default=None)

    p_fit = sub.add_parser("fit-predict")
    p_fit.add_argument("--model-type", required=True)
    p_fit.add_argument("--task", required=True)
    p_fit.add_argument("--x-train-csv", required=True)
    p_fit.add_argument("--y-train-csv", required=True)
    p_fit.add_argument("--x-test-csv", required=True)
    p_fit.add_argument("--params-json", default=None)
    p_fit.add_argument("--class-levels-json", default=None)
    p_fit.add_argument("--out-dir", required=True)

    p_batch = sub.add_parser("batch-fit-predict")
    p_batch.add_argument("--manifest-json", required=True)

    p_boot = sub.add_parser("bootstrap-fit-predict")
    p_boot.add_argument("--model-type", required=True)
    p_boot.add_argument("--task", required=True)
    p_boot.add_argument("--x-train-csv", required=True)
    p_boot.add_argument("--y-train-csv", required=True)
    p_boot.add_argument("--x-pred-csv", required=True)
    p_boot.add_argument("--params-json", default=None)
    p_boot.add_argument("--class-levels-json", default=None)
    p_boot.add_argument("--n-bootstrap", type=int, required=True)
    p_boot.add_argument("--seed", type=int, default=123)
    p_boot.add_argument("--training-seeds-json", default=None)
    p_boot.add_argument("--seed-aggregation", default="mean")
    p_boot.add_argument("--out-dir", required=True)

    p_prewarm = sub.add_parser("prewarm")
    p_prewarm.add_argument("--out-dir", default=None)

    args = parser.parse_args(argv)
    if args.cmd == "setup-deps":
        return _setup_deps(args)
    if args.cmd == "fit-predict":
        return _fit_predict(args)
    if args.cmd == "batch-fit-predict":
        return _batch_fit_predict(args)
    if args.cmd == "bootstrap-fit-predict":
        return _bootstrap_predict(args)
    if args.cmd == "prewarm":
        return _prewarm(args)
    raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
