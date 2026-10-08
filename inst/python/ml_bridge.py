from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


def _module_dir() -> Path:
    return Path(__file__).resolve().parent


def _load_ml_models():
    mod_path = _module_dir() / "ml_models.py"
    if not mod_path.exists():
        raise FileNotFoundError(f"ml_models.py not found: {mod_path}")
    spec = importlib.util.spec_from_file_location("predictpror_ml_models", mod_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _read_matrix_csv(path: str) -> np.ndarray:
    df = pd.read_csv(path)
    return df.to_numpy(dtype=float)


def _read_target_csv(path: str, task: str):
    df = pd.read_csv(path)
    vals = df.iloc[:, 0].to_numpy()
    if task == "gaussian":
        return vals.astype(float)
    return vals.astype(str)


def _read_json(path: Optional[str], default: Any):
    if path is None:
        return default
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_json(obj: Dict[str, Any], path: str):
    Path(path).write_text(json.dumps(obj, indent=2), encoding="utf-8")


def _write_predictions(result: Dict[str, Any], task: str, out_dir: Path):
    preds = np.asarray(result["predictions"])
    if task == "gaussian":
      pred_df = pd.DataFrame({"prediction": preds.astype(float)})
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
    }
    if result.get("selected_alpha") is not None:
        meta["selected_alpha"] = float(result["selected_alpha"])
    _write_json(meta, str(out_dir / "meta.json"))


def _write_feature_importance(result: Dict[str, Any], out_dir: Path):
    importance = np.asarray(result["importance"], dtype=float).reshape(-1)
    pd.DataFrame({"importance": importance}).to_csv(out_dir / "feature_importance.csv", index=False)


def fit_predict_payload(
    model: str,
    task: str,
    X_train,
    y_train,
    X_test,
    params: Optional[Dict[str, Any]] = None,
    class_levels = None,
):
    ml_models = _load_ml_models()
    return _fit_predict_payload_with_module(
        ml_models=ml_models,
        model=model,
        task=task,
        X_train=X_train,
        y_train=y_train,
        X_test=X_test,
        params=params,
        class_levels=class_levels,
    )


def _fit_predict_payload_with_module(
    ml_models,
    model: str,
    task: str,
    X_train,
    y_train,
    X_test,
    params: Optional[Dict[str, Any]] = None,
    class_levels = None,
):
    return ml_models.fit_predict(
        model=model,
        X_train=np.asarray(X_train),
        y_train=np.asarray(y_train),
        X_test=np.asarray(X_test),
        task=task,
        class_levels=class_levels,
        params=params or {},
    )


def feature_importance_payload(
    model: str,
    task: str,
    X_train,
    y_train,
    params: Optional[Dict[str, Any]] = None,
    class_levels = None,
):
    ml_models = _load_ml_models()
    return ml_models.feature_importance(
        model=model,
        X_train=np.asarray(X_train),
        y_train=np.asarray(y_train),
        task=task,
        class_levels=class_levels,
        params=params or {},
    )


def _bootstrap_row(res, task: str, class_levels):
    """One bootstrap replicate's output row and the class labels it used."""
    if task == "multiclass":
        prob = np.asarray(res["probabilities"], dtype=float)
        if prob.ndim == 1:
            prob = prob.reshape(-1, 1)
        classes = [str(x) for x in (res.get("classes") or class_levels or list(range(prob.shape[1])))]
        return np.asarray(prob.T, dtype=float).reshape(-1), classes
    if task == "binary":
        prob = np.asarray(res.get("probabilities"), dtype=float)
        if prob.ndim == 1:
            return prob.reshape(-1), None
        classes = [str(x) for x in (res.get("classes") or class_levels or list(range(prob.shape[1])))]
        positive_class = str((class_levels or classes or ["0", "1"])[min(1, len(class_levels or classes or ["0", "1"]) - 1)])
        try:
            positive_idx = classes.index(positive_class)
        except ValueError:
            positive_idx = min(1, prob.shape[1] - 1)
        return np.asarray(prob[:, positive_idx], dtype=float).reshape(-1), None
    return np.asarray(res["predictions"], dtype=float).reshape(-1), None


# Process-pool workers receive the training data once (initializer) and then
# only the resample indices of each replicate.
_BOOT_STATE: Dict[str, Any] = {}


def _bootstrap_worker_init(model, task, X_train, y_train, X_pred, params, class_levels):
    _BOOT_STATE.update(dict(
        model=model, task=task, X_train=X_train, y_train=y_train, X_pred=X_pred,
        params=params, class_levels=class_levels, ml_models=_load_ml_models(),
    ))


def _bootstrap_worker_fit(idx):
    s = _BOOT_STATE
    res = _fit_predict_payload_with_module(
        ml_models=s["ml_models"], model=s["model"], task=s["task"],
        X_train=s["X_train"][idx, :], y_train=s["y_train"][idx], X_test=s["X_pred"],
        params=s["params"], class_levels=s["class_levels"],
    )
    return _bootstrap_row(res, s["task"], s["class_levels"])


def _rf_infinitesimal_jackknife(X_train, y_train, X_pred, params, n_draws):
    """Gaussian random-forest prediction uncertainty without refits.

    One forest is fitted on all training rows. With per-tree predictions
    t_b(x) and in-bag counts N_bi, the bias-corrected infinitesimal jackknife
    (Wager, Hastie & Efron 2014, JMLR 15:1625) is
        V_IJ(x)   = sum_i Cov_b(N_bi, t_b(x))^2
        V_IJ-U(x) = V_IJ(x) - n * v_N * Var_b(t_b(x)) / B
    (v_N = mean variance of the in-bag counts; ~1 for the bootstrap).
    The variance is returned as `n_draws` deterministic normal-quantile rows
    centred on the forest prediction (column sd = sqrt(V_IJ-U)), so code that
    reads a bootstrap matrix gets the jackknife variance and a normal interval.
    """
    from sklearn.ensemble._forest import _generate_sample_indices, _get_n_samples_bootstrap
    ml_models = _load_ml_models()
    estimator = ml_models._build_estimator(model="randomforest", task="gaussian",
                                          params=ml_models._clean_params(params), n_classes=None)
    X_train = np.asarray(X_train, dtype=float)
    X_pred = np.asarray(X_pred, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    estimator.fit(X_train, y_train)
    n = X_train.shape[0]
    n_boot = _get_n_samples_bootstrap(n, estimator.max_samples)
    trees = estimator.estimators_
    B = len(trees)
    inbag = np.zeros((n, B))
    for b, tree in enumerate(trees):
        idx = _generate_sample_indices(tree.random_state, n, n_boot)
        inbag[:, b] = np.bincount(idx, minlength=n)
    tree_pred = np.column_stack([t.predict(X_pred) for t in trees])  # (m, B)
    mean_pred = tree_pred.mean(axis=1)
    pred_c = tree_pred - mean_pred[:, None]
    inbag_c = inbag - inbag.mean(axis=1, keepdims=True)
    v_ij = np.sum((inbag_c @ pred_c.T / B) ** 2, axis=0)  # (m,)
    v_n = float(np.mean(inbag.var(axis=1)))
    bias = n * v_n * np.sum(pred_c ** 2, axis=1) / (B ** 2)
    v_u = v_ij - bias
    # The correction can undershoot for small B; keep the uncorrected value as
    # a floor proportion so the variance never collapses to zero.
    v_u = np.where(v_u > 0, v_u, 0.05 * v_ij)
    k = max(2, int(n_draws))
    from statistics import NormalDist
    z = np.array([NormalDist().inv_cdf((i + 0.5) / k) for i in range(k)])
    z = (z - z.mean()) / z.std(ddof=1)
    boot = mean_pred[None, :] + z[:, None] * np.sqrt(v_u)[None, :]
    return boot


def bootstrap_fit_predict_payload(
    model: str,
    task: str,
    X_train,
    y_train,
    X_pred,
    n_bootstrap: int,
    seed: int = 12345,
    params: Optional[Dict[str, Any]] = None,
    class_levels = None,
    n_jobs: int = 1,
):
    X_train = np.asarray(X_train)
    y_train = np.asarray(y_train)
    X_pred = np.asarray(X_pred)
    import os
    # Opt-in only. Simulation (n=150-200, 600-2000 markers, 15-20 training
    # sets each): unbiased on average with 1000 trees but noisier than 30
    # bootstrap forests, with lower 90% interval coverage (0.81 vs 0.89);
    # conservative (+15% to +37% SD) with 300 trees.
    if (model == "randomforest" and task == "gaussian" and
            os.environ.get("PREDICTPRO_RF_UNCERTAINTY", "bootstrap").strip().lower() == "jackknife"):
        boot = _rf_infinitesimal_jackknife(X_train, y_train, X_pred, dict(params or {}), n_bootstrap)
        return {
            "bootstrap": boot,
            "task": task,
            "n_bootstrap": int(boot.shape[0]),
            "n_prediction_rows": int(X_pred.shape[0]),
            "classes": [],
            "uncertainty_method": "random_forest_infinitesimal_jackknife",
        }
    rng = np.random.default_rng(int(seed))
    n = X_train.shape[0]
    classes_out = class_levels or []
    # Draw every replicate's indices up front from the one seeded stream, so
    # results are identical however many processes fit them.
    all_idx = [rng.integers(0, n, size=n, endpoint=False) for _ in range(int(n_bootstrap))]
    n_jobs = max(1, min(int(n_jobs or 1), len(all_idx)))
    params = dict(params or {})
    serial_params = dict(params)

    outputs = None
    if n_jobs > 1:
        import concurrent.futures as cf
        pool_params = dict(params)
        # each replicate fits single-threaded; parallelism is across replicates
        for key in ("n_jobs", "rf_n_jobs", "nthread", "thread_count", "xgb_nthread",
                    "lightgbm_nthread", "catboost_thread_count"):
            if key in pool_params:
                pool_params[key] = 1
        # pool processes would all queue on one GPU; keep replicates on CPU
        if str(pool_params.get("device", "")).lower() in ("cuda", "gpu"):
            pool_params["device"] = "cpu"
        try:
            with cf.ProcessPoolExecutor(
                max_workers=n_jobs,
                initializer=_bootstrap_worker_init,
                initargs=(model, task, X_train, y_train, X_pred, pool_params, class_levels),
            ) as pool:
                outputs = list(pool.map(_bootstrap_worker_fit, all_idx))
        except (OSError, MemoryError, RuntimeError, cf.process.BrokenProcessPool) as exc:
            # Starting pool processes can fail under memory pressure (Windows:
            # "OSError: [Errno 22]" while sending the start-up data). The
            # resample indices are pre-drawn, so in-process refits give the
            # same result; fall back instead of losing the whole prediction.
            sys.stderr.write(f"Bootstrap process pool failed ({exc!r}); running the refits in-process.\n")
            outputs = None
    if outputs is None:
        ml_models = _load_ml_models()
        outputs = []
        for idx in all_idx:
            res = _fit_predict_payload_with_module(
                ml_models=ml_models, model=model, task=task,
                X_train=X_train[idx, :], y_train=y_train[idx], X_test=X_pred,
                params=serial_params, class_levels=class_levels,
            )
            outputs.append(_bootstrap_row(res, task, class_levels))

    rows = [row for row, _ in outputs]
    for _, classes in outputs:
        if classes:
            classes_out = classes
    boot = np.vstack(rows) if rows else np.empty((0, X_pred.shape[0]), dtype=float)
    return {
        "bootstrap": boot,
        "task": task,
        "n_bootstrap": int(n_bootstrap),
        "n_prediction_rows": int(X_pred.shape[0]),
        "classes": classes_out,
    }


def _fit_predict(args) -> int:
    params = _read_json(args.params_json, {})
    class_levels = _read_json(args.class_levels_json, None)
    ml_models = _load_ml_models()
    result = _fit_predict_payload_with_module(
        ml_models=ml_models,
        model=args.model,
        task=args.task,
        X_train=_read_matrix_csv(args.x_train_csv),
        y_train=_read_target_csv(args.y_train_csv, args.task),
        X_test=_read_matrix_csv(args.x_test_csv),
        params=params,
        class_levels=class_levels,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_predictions(result, args.task, out_dir)
    return 0


def _batch_fit_predict(args) -> int:
    manifest = _read_json(args.manifest_json, {})
    jobs = manifest.get("jobs", manifest if isinstance(manifest, list) else [])
    if not isinstance(jobs, list):
        raise ValueError("Batch manifest must contain a list under 'jobs'.")

    ml_models = _load_ml_models()
    completed = []
    for i, job in enumerate(jobs):
        task = str(job["task"]).strip().lower()
        params = _read_json(job.get("params_json"), {})
        class_levels = _read_json(job.get("class_levels_json"), None)
        result = _fit_predict_payload_with_module(
            ml_models=ml_models,
            model=str(job["model"]).strip().lower(),
            task=task,
            X_train=_read_matrix_csv(job["x_train_csv"]),
            y_train=_read_target_csv(job["y_train_csv"], task),
            X_test=_read_matrix_csv(job["x_test_csv"]),
            params=params,
            class_levels=class_levels,
        )
        out_dir = Path(job["out_dir"])
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_predictions(result, task, out_dir)
        completed.append({"id": str(job.get("id", i + 1)), "out_dir": str(out_dir)})

    summary_path = manifest.get("summary_json") if isinstance(manifest, dict) else None
    if summary_path:
        _write_json({"jobs": completed}, summary_path)
    return 0


def _bootstrap_predict(args) -> int:
    params = _read_json(args.params_json, {})
    class_levels = _read_json(args.class_levels_json, None)
    payload = bootstrap_fit_predict_payload(
        model=args.model,
        task=args.task,
        X_train=_read_matrix_csv(args.x_train_csv),
        y_train=_read_target_csv(args.y_train_csv, args.task),
        X_pred=_read_matrix_csv(args.x_pred_csv),
        n_bootstrap=int(args.n_bootstrap),
        seed=int(args.seed),
        params=params,
        class_levels=class_levels,
        n_jobs=int(getattr(args, "n_jobs", 1) or 1),
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(payload["bootstrap"]).to_csv(out_dir / "bootstrap.csv", index=False)
    meta = {
        "task": payload["task"],
        "n_bootstrap": int(payload["n_bootstrap"]),
        "n_prediction_rows": int(payload["n_prediction_rows"]),
        "classes": payload["classes"],
    }
    _write_json(meta, str(out_dir / "meta.json"))
    return 0


def _feature_importance(args) -> int:
    params = _read_json(args.params_json, {})
    class_levels = _read_json(args.class_levels_json, None)
    payload = feature_importance_payload(
        model=args.model,
        task=args.task,
        X_train=_read_matrix_csv(args.x_train_csv),
        y_train=_read_target_csv(args.y_train_csv, args.task),
        params=params,
        class_levels=class_levels,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_feature_importance(payload, out_dir)
    return 0


def _setup_deps(args) -> int:
    ml_models = _load_ml_models()
    info = ml_models.setup_deps(prefer_gpu=bool(args.prefer_gpu))
    for k, v in info.items():
        print(f"{k}={v}")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_setup = sub.add_parser("setup-deps")
    p_setup.add_argument("--prefer-gpu", action="store_true")

    p_fit = sub.add_parser("fit-predict")
    p_fit.add_argument("--model", required=True)
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
    p_boot.add_argument("--model", required=True)
    p_boot.add_argument("--task", required=True)
    p_boot.add_argument("--x-train-csv", required=True)
    p_boot.add_argument("--y-train-csv", required=True)
    p_boot.add_argument("--x-pred-csv", required=True)
    p_boot.add_argument("--params-json", default=None)
    p_boot.add_argument("--class-levels-json", default=None)
    p_boot.add_argument("--n-bootstrap", type=int, required=True)
    p_boot.add_argument("--n-jobs", type=int, default=1)
    p_boot.add_argument("--seed", type=int, default=12345)
    p_boot.add_argument("--out-dir", required=True)

    p_imp = sub.add_parser("feature-importance")
    p_imp.add_argument("--model", required=True)
    p_imp.add_argument("--task", required=True)
    p_imp.add_argument("--x-train-csv", required=True)
    p_imp.add_argument("--y-train-csv", required=True)
    p_imp.add_argument("--params-json", default=None)
    p_imp.add_argument("--class-levels-json", default=None)
    p_imp.add_argument("--out-dir", required=True)

    args = parser.parse_args(argv)
    if args.cmd == "setup-deps":
        return _setup_deps(args)
    if args.cmd == "fit-predict":
        return _fit_predict(args)
    if args.cmd == "batch-fit-predict":
        return _batch_fit_predict(args)
    if args.cmd == "bootstrap-fit-predict":
        return _bootstrap_predict(args)
    if args.cmd == "feature-importance":
        return _feature_importance(args)
    raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
