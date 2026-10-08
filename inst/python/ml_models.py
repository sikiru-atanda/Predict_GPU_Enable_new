from __future__ import annotations

import importlib
import subprocess
import sys
import warnings
from typing import Any, Dict, Optional, Sequence

import numpy as np


def _run(cmd: Sequence[str]):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, text=True)


def _pip_install(pkgs: Sequence[str]):
    cmd = [sys.executable, "-m", "pip", "install", "--upgrade", "--quiet", "--no-cache-dir", *pkgs]
    return _run(cmd)


def _lazy_import(name: str, pip_name: Optional[str] = None):
    try:
        return importlib.import_module(name)
    except Exception:
        _pip_install([pip_name or name])
        return importlib.import_module(name)


def setup_deps(prefer_gpu: bool = True):
    sklearn = _lazy_import("sklearn", "scikit-learn")
    xgboost = importlib.import_module("xgboost") if importlib.util.find_spec("xgboost") else None
    catboost = importlib.import_module("catboost") if importlib.util.find_spec("catboost") else None
    lightgbm = importlib.import_module("lightgbm") if importlib.util.find_spec("lightgbm") else None
    return {
        "sklearn_version": getattr(sklearn, "__version__", None),
        "xgboost_version": getattr(xgboost, "__version__", None),
        "catboost_version": getattr(catboost, "__version__", None),
        "lightgbm_version": getattr(lightgbm, "__version__", None),
        "prefer_gpu": bool(prefer_gpu),
    }


def _as_2d_float(x: Any) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def _predict_quiet(estimator: Any, X_test: np.ndarray):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="X does not have valid feature names, but LGBM",
            category=UserWarning,
        )
        pred = estimator.predict(X_test)
        prob = estimator.predict_proba(X_test) if hasattr(estimator, "predict_proba") else None
    return pred, prob


def _as_target(y: Any, task: str) -> np.ndarray:
    arr = np.asarray(y)
    if task == "gaussian":
        return arr.astype(float)
    return arr.astype(str)


def _normalize_kernel(kernel: Optional[str]) -> str:
    mapping = {
        "gaussian": "rbf",
        "radial": "rbf",
        "linear": "linear",
        "polynomial": "poly",
        "hyperbolic_tangent": "sigmoid",
        "sigmoid": "sigmoid",
    }
    key = (kernel or "rbf").strip().lower()
    return mapping.get(key, key)


def _normalize_svm_type(task: str, svm_type: Optional[str]) -> str:
    key = (svm_type or "").strip().lower()
    if task == "gaussian":
        mapping = {
            "eps-regression": "eps-regression",
            "epsilon-regression": "eps-regression",
            "eps": "eps-regression",
            "nu-regression": "nu-regression",
            "nu": "nu-regression",
        }
        return mapping.get(key, "eps-regression")
    mapping = {
        "c-classification": "c-classification",
        "c": "c-classification",
        "nu-classification": "nu-classification",
        "nu": "nu-classification",
    }
    return mapping.get(key, "c-classification")


def _clean_params(params: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    out = {}
    for key, value in dict(params or {}).items():
        if value is None:
            continue
        if isinstance(value, np.generic):
            value = value.item()
        out[str(key)] = value
    return out


def _svm_gamma_from_params(params: Dict[str, Any]):
    if params.get("gamma_value") is not None:
        return float(params["gamma_value"])
    if params.get("sigma_value") is not None:
        sigma = float(params["sigma_value"])
        if sigma > 0:
            return 1.0 / (2.0 * (sigma ** 2))
    return None


def _build_estimator(model: str,
                     task: str,
                     params: Dict[str, Any],
                     n_classes: Optional[int] = None):
    model = model.lower()
    random_state = params.get("random_state", 42)

    if model == "xgboost":
      xgb = _lazy_import("xgboost")
      common = {
          "random_state": random_state,
          "n_estimators": int(params.get("nrounds", 100)),
          "learning_rate": float(params.get("eta", 0.1)),
          "max_depth": int(params.get("max_depth", 6)),
          "subsample": float(params.get("subsample", 1.0)),
          "colsample_bytree": float(params.get("colsample_bytree", 1.0)),
          "min_child_weight": float(params.get("min_child_weight", 1.0)),
          "gamma": float(params.get("xgb_gamma", params.get("gamma", 0.0))),
          "reg_alpha": float(params.get("xgb_alpha", params.get("alpha", 0.0))),
          "reg_lambda": float(params.get("xgb_lambda", params.get("lambda", 1.0))),
          "booster": params.get("xgb_booster", params.get("booster", "gbtree")),
          "n_jobs": int(params.get("xgb_nthread", params.get("n_jobs", 1))),
          "verbosity": 0,
      }
      device = params.get("device")
      if device:
          common["device"] = device
          common.setdefault("tree_method", "hist")
      if task == "gaussian":
          common["objective"] = "reg:squarederror"
          return xgb.XGBRegressor(**common)
      if task == "binary":
          common["objective"] = "binary:logistic"
          common["eval_metric"] = "logloss"
          return xgb.XGBClassifier(**common)
      common["objective"] = "multi:softprob"
      common["num_class"] = int(n_classes or 2)
      common["eval_metric"] = "mlogloss"
      return xgb.XGBClassifier(**common)

    if model == "randomforest":
        sklearn_ensemble = _lazy_import("sklearn.ensemble")
        common = {
            "n_estimators": int(params.get("ntree", 500)),
            "random_state": random_state,
            "n_jobs": int(params.get("n_jobs", params.get("rf_n_jobs", 1))),
        }
        if params.get("mtry") is not None:
            common["max_features"] = int(params["mtry"])
        elif task == "gaussian":
            # R randomForest regression default mtry = floor(p/3). sklearn's
            # regressor default (all features) is bagged trees, not a random
            # forest, and ~2x slower per tree on marker data. Classification
            # keeps sklearn's sqrt(p), which matches R.
            common["max_features"] = 1.0 / 3.0
        if params.get("maxnodes") is not None:
            common["max_leaf_nodes"] = int(params["maxnodes"])
        if params.get("nodesize") is not None:
            common["min_samples_leaf"] = int(params["nodesize"])
        if task == "gaussian":
            return sklearn_ensemble.RandomForestRegressor(**common)
        return sklearn_ensemble.RandomForestClassifier(**common)

    if model == "catboost":
        cb = _lazy_import("catboost")
        common = {
            "iterations": int(params.get("catboost_iterations", params.get("nrounds", 500))),
            "depth": int(params.get("catboost_depth", params.get("max_depth", 6))),
            "learning_rate": float(params.get("catboost_learning_rate", params.get("eta", 0.03))),
            "l2_leaf_reg": float(params.get("catboost_l2_leaf_reg", 3.0)),
            "random_state": random_state,
            "thread_count": int(params.get("catboost_thread_count", params.get("n_jobs", 1))),
            "verbose": False,
            "allow_writing_files": False,
        }
        if task == "gaussian":
            common["loss_function"] = "RMSE"
            return cb.CatBoostRegressor(**common)
        if task == "binary":
            common["loss_function"] = "Logloss"
            return cb.CatBoostClassifier(**common)
        common["loss_function"] = "MultiClass"
        return cb.CatBoostClassifier(**common)

    if model == "lightgbm":
        lgbm = _lazy_import("lightgbm")
        common = {
            "n_estimators": int(params.get("lightgbm_nrounds", params.get("nrounds", 100))),
            "learning_rate": float(params.get("lightgbm_learning_rate", params.get("eta", 0.05))),
            "num_leaves": int(params.get("lightgbm_num_leaves", 31)),
            "feature_fraction": float(params.get("lightgbm_feature_fraction", 1.0)),
            "bagging_fraction": float(params.get("lightgbm_bagging_fraction", 1.0)),
            "min_child_samples": int(params.get("lightgbm_min_data_in_leaf", 20)),
            "reg_alpha": float(params.get("lightgbm_lambda_l1", 0.0)),
            "reg_lambda": float(params.get("lightgbm_lambda_l2", 0.0)),
            "random_state": random_state,
            "n_jobs": int(params.get("lightgbm_nthread", params.get("n_jobs", 1))),
            "verbosity": -1,
        }
        if task == "gaussian":
            return lgbm.LGBMRegressor(objective="regression", **common)
        if task == "binary":
            return lgbm.LGBMClassifier(objective="binary", **common)
        return lgbm.LGBMClassifier(objective="multiclass", num_class=int(n_classes or 2), **common)

    if model == "svm":
        sklearn_svm = _lazy_import("sklearn.svm")
        kernel = _normalize_kernel(params.get("svm_kernel"))
        gamma = _svm_gamma_from_params(params)
        svm_type = _normalize_svm_type(task, params.get("svm_type"))
        common = {
            "kernel": kernel,
            "C": float(params.get("C_value", 1.0)),
        }
        if kernel == "poly":
            common["degree"] = int(params.get("degree_value", 3))
        if kernel in {"rbf", "poly", "sigmoid"} and gamma is not None:
            common["gamma"] = gamma
        if kernel in {"poly", "sigmoid"}:
            common["coef0"] = float(params.get("offset_value", 0.0))
        if task == "gaussian":
            if svm_type == "nu-regression":
                common["nu"] = float(params.get("nu_value", 0.5))
                return sklearn_svm.NuSVR(**common)
            return sklearn_svm.SVR(**common)
        if svm_type == "nu-classification":
            common.pop("C", None)
            common["nu"] = float(params.get("nu_value", 0.5))
            return sklearn_svm.NuSVC(probability=True, random_state=random_state, **common)
        return sklearn_svm.SVC(probability=True, random_state=random_state, **common)

    if model == "pls":
        sklearn_cross = _lazy_import("sklearn.cross_decomposition")
        sklearn_model_selection = _lazy_import("sklearn.model_selection")
        max_components = int(params.get("pls_max_components", params.get("ncomp", 3)))
        if task != "gaussian":
            raise ValueError("PLS only supports gaussian regression.")
        max_components = max(1, max_components)
        auto_components = bool(params.get("pls_auto_components", False))
        if auto_components:
            return {
                "kind": "pls_auto",
                "max_components": max_components,
                "pls_module": sklearn_cross,
                "model_selection_module": sklearn_model_selection,
                "random_state": random_state,
            }
        return sklearn_cross.PLSRegression(n_components=max_components, scale=False)

    if model == "ridge":
        sklearn_linear = _lazy_import("sklearn.linear_model")
        if task != "gaussian":
            raise ValueError("Ridge only supports gaussian regression.")
        alpha = params.get("ridge_alpha", params.get("lambda_rr"))
        if alpha is not None:
            return sklearn_linear.Ridge(alpha=float(alpha), random_state=random_state)
        alphas = params.get("ridge_alphas")
        if alphas is None:
            alphas = np.logspace(-6, 3, 60)
        return sklearn_linear.RidgeCV(alphas=np.asarray(alphas, dtype=float))

    if model == "lasso":
        sklearn_linear = _lazy_import("sklearn.linear_model")
        if task != "gaussian":
            raise ValueError("Lasso only supports gaussian regression.")
        alpha = params.get("lasso_alpha")
        if alpha is not None:
            return sklearn_linear.Lasso(alpha=float(alpha), random_state=random_state, max_iter=10000)
        alphas = params.get("lasso_alphas")
        if alphas is None:
            alphas = np.logspace(-6, 0, 60)
        return sklearn_linear.LassoCV(alphas=np.asarray(alphas, dtype=float), random_state=random_state, max_iter=10000)

    if model == "knn":
        sklearn_neighbors = _lazy_import("sklearn.neighbors")
        k = int(params.get("k", 5))
        if task == "gaussian":
            return sklearn_neighbors.KNeighborsRegressor(n_neighbors=k)
        return sklearn_neighbors.KNeighborsClassifier(n_neighbors=k)

    raise ValueError(f"Unsupported model: {model}")


class _PrecomputedKernelSVM:
    """SVR/SVC fitted on a BLAS-computed kernel matrix.

    libsvm evaluates kernel entries one at a time on a single thread. With
    1,451 lines x 10,346 markers one RBF SVR fit took 564 s and its prediction
    762 s; the same model on a precomputed kernel took 1.2 s, with predictions
    equal to 4e-14 and the same support vectors. Kernels follow libsvm's
    definitions exactly; gamma=None means sklearn's "scale" rule.
    """

    def __init__(self, inner, kernel: str, gamma=None, degree: int = 3, coef0: float = 0.0):
        self.inner = inner
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0

    def _k(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        if self.kernel == "linear":
            return A @ B.T
        if self.kernel == "rbf":
            sq = (A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2.0 * (A @ B.T)
            return np.exp(-self.gamma_ * np.maximum(sq, 0.0))
        if self.kernel == "poly":
            return (self.gamma_ * (A @ B.T) + self.coef0) ** self.degree
        return np.tanh(self.gamma_ * (A @ B.T) + self.coef0)       # sigmoid

    def fit(self, X, y):
        self.X_fit_ = np.asarray(X, dtype=float)
        if self.gamma is None:
            var = self.X_fit_.var()
            self.gamma_ = 1.0 / (self.X_fit_.shape[1] * var) if var != 0 else 1.0
        else:
            self.gamma_ = float(self.gamma)
        self.inner.set_params(kernel="precomputed")
        self.inner.fit(self._k(self.X_fit_, self.X_fit_), y)
        return self

    def predict(self, X):
        return self.inner.predict(self._k(np.asarray(X, dtype=float), self.X_fit_))

    def predict_proba(self, X):
        return self.inner.predict_proba(self._k(np.asarray(X, dtype=float), self.X_fit_))

    @property
    def classes_(self):
        return self.inner.classes_


def _maybe_precompute_svm(estimator: Any, n_train: int) -> Any:
    """Wrap a libsvm estimator so its kernel is computed with BLAS, unless the
    n x n kernel matrix would be too large (PREDICTPRO_SVM_PRECOMPUTE_MAX_N)."""
    import os
    kernel = getattr(estimator, "kernel", None)
    if kernel not in {"linear", "rbf", "poly", "sigmoid"}:
        return estimator
    try:
        max_n = int(os.environ.get("PREDICTPRO_SVM_PRECOMPUTE_MAX_N", "15000"))
    except ValueError:
        max_n = 15000
    if n_train > max_n:
        return estimator
    gamma = getattr(estimator, "gamma", "scale")
    return _PrecomputedKernelSVM(
        inner=estimator,
        kernel=kernel,
        gamma=None if isinstance(gamma, str) else float(gamma),
        degree=int(getattr(estimator, "degree", 3)),
        coef0=float(getattr(estimator, "coef0", 0.0)),
    )


def _centred_rank(X: np.ndarray, rel_tol: float = 1e-10) -> int:
    """Numerical rank of the column-centred matrix (what PLS decomposes)."""
    Xc = X - X.mean(axis=0, keepdims=True)
    if Xc.size == 0:
        return 0
    sv = np.linalg.svd(Xc, compute_uv=False)
    if not sv.size or sv[0] <= 0:
        return 0
    return int(np.sum(sv > sv[0] * rel_tol))


def _cap_pls_components(estimator: Any, X: np.ndarray) -> Any:
    """Keep PLS components within the rank of this fit's training data.

    Components beyond the rank of the centred X only fit residual noise along
    near-null directions and make the coefficients explode. Bootstrap
    resamples and CV folds have fewer unique rows than the outer training set,
    so an ncomp valid for the full data can exceed a resample's rank.
    """
    n_comp = getattr(estimator, "n_components", None)
    if n_comp is None or type(estimator).__name__ != "PLSRegression":
        return estimator
    cap = max(1, _centred_rank(X))
    if n_comp > cap:
        estimator.set_params(n_components=cap)
    return estimator


def fit_predict(model: str,
                X_train: Any,
                y_train: Any,
                X_test: Any,
                task: str = "gaussian",
                class_levels: Optional[Sequence[str]] = None,
                params: Optional[Dict[str, Any]] = None):
    task = (task or "gaussian").strip().lower()
    params = _clean_params(params)
    X_train = _as_2d_float(X_train)
    X_test = _as_2d_float(X_test)
    y_train = _as_target(y_train, task)
    class_levels = [str(x) for x in (class_levels or np.unique(y_train).tolist())] if task != "gaussian" else None
    n_classes = len(class_levels) if task != "gaussian" else None

    y_fit = y_train
    decode_lookup = None
    if task != "gaussian":
        observed = [str(x) for x in np.unique(y_train).tolist()]
        fit_levels = [x for x in class_levels if x in set(observed)]
        if not fit_levels:
            fit_levels = observed
        if len(fit_levels) == 1:
            cls = fit_levels[0]
            prob = np.zeros((X_test.shape[0], len(class_levels)), dtype=float)
            if cls in class_levels:
                prob[:, class_levels.index(cls)] = 1.0
            return {
                "predictions": [cls] * X_test.shape[0],
                "probabilities": prob.tolist(),
                "classes": class_levels,
            }
        decode_lookup = {idx: label for idx, label in enumerate(fit_levels)}
        encode_lookup = {label: idx for idx, label in decode_lookup.items()}
        y_fit = np.asarray([encode_lookup[str(val)] for val in y_train], dtype=int)
        n_classes = len(fit_levels)

    estimator = _build_estimator(model=model, task=task, params=params, n_classes=n_classes)
    if model.strip().lower() == "svm" and not isinstance(estimator, dict):
        estimator = _maybe_precompute_svm(estimator, X_train.shape[0])
    if isinstance(estimator, dict) and estimator.get("kind") == "pls_auto":
        max_components = min(_centred_rank(X_train), X_train.shape[1], int(estimator["max_components"]))
        max_components = max(1, max_components)
        kfold_splits = min(5, X_train.shape[0])
        kfold_splits = max(2, kfold_splits)
        scorer = "neg_mean_squared_error"
        best_comp = 1
        best_score = -np.inf
        for n_comp in range(1, max_components + 1):
            est = estimator["pls_module"].PLSRegression(n_components=n_comp, scale=False)
            cv = estimator["model_selection_module"].KFold(
                n_splits=kfold_splits,
                shuffle=True,
                random_state=estimator["random_state"],
            )
            score = np.mean(
                estimator["model_selection_module"].cross_val_score(
                    est,
                    X_train,
                    y_train.astype(float),
                    cv=cv,
                    scoring=scorer,
                )
            )
            if score > best_score:
                best_score = score
                best_comp = n_comp
        estimator = estimator["pls_module"].PLSRegression(n_components=best_comp, scale=False)

    estimator = _cap_pls_components(estimator, X_train)
    estimator.fit(X_train, y_fit)

    if task == "gaussian":
        if model == "lightgbm":
            pred, _ = _predict_quiet(estimator, X_test)
        else:
            pred = estimator.predict(X_test)
        out = {"predictions": np.asarray(pred, dtype=float).tolist()}
        # Penalty chosen by an internal CV estimator (LassoCV/RidgeCV), so the
        # caller can reuse it for refits instead of repeating the search.
        if hasattr(estimator, "alpha_"):
            out["selected_alpha"] = float(estimator.alpha_)
        return out

    if model == "lightgbm":
        pred, prob = _predict_quiet(estimator, X_test)
    else:
        pred = estimator.predict(X_test)
        prob = estimator.predict_proba(X_test)
    pred = [decode_lookup.get(int(x), str(x)) for x in np.asarray(pred).ravel()]
    fit_classes = [decode_lookup.get(int(x), str(x)) for x in getattr(estimator, "classes_", range(n_classes))]
    classes = fit_classes
    if class_levels:
        full_prob = np.zeros((X_test.shape[0], len(class_levels)), dtype=float)
        for j, cls in enumerate(fit_classes):
            if cls in class_levels and j < np.asarray(prob).shape[1]:
                full_prob[:, class_levels.index(cls)] = np.asarray(prob, dtype=float)[:, j]
        prob = full_prob
        classes = class_levels
    return {
        "predictions": np.asarray(pred).astype(str).tolist(),
        "probabilities": np.asarray(prob, dtype=float).tolist(),
        "classes": classes,
    }


def feature_importance(model: str,
                       X_train: Any,
                       y_train: Any,
                       task: str = "gaussian",
                       class_levels: Optional[Sequence[str]] = None,
                       params: Optional[Dict[str, Any]] = None):
    task = (task or "gaussian").strip().lower()
    params = _clean_params(params)
    X_train = _as_2d_float(X_train)
    y_train = _as_target(y_train, task)
    class_levels = [str(x) for x in (class_levels or np.unique(y_train).tolist())] if task != "gaussian" else None
    n_classes = len(class_levels) if task != "gaussian" else None

    estimator = _build_estimator(
        model=model,
        task=task,
        params=params,
        n_classes=n_classes,
    )
    if isinstance(estimator, dict):
        raise ValueError(f"Feature importance is not available for auto estimator: {model}")

    estimator = _cap_pls_components(estimator, X_train)
    estimator.fit(X_train, y_train)
    if hasattr(estimator, "feature_importances_"):
        importance = np.asarray(estimator.feature_importances_, dtype=float)
    elif hasattr(estimator, "coef_"):
        coef = np.asarray(estimator.coef_, dtype=float)
        importance = np.mean(np.abs(coef), axis=0) if coef.ndim > 1 else np.abs(coef)
    else:
        raise ValueError(f"Feature importance is not available for model: {model}")

    if importance.shape[0] != X_train.shape[1]:
        raise ValueError(
            f"Feature importance length {importance.shape[0]} does not match "
            f"predictor count {X_train.shape[1]}"
        )
    return {"importance": importance.astype(float).tolist()}
