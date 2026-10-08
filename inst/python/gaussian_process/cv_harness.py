"""Cross-validation harness.

Contract:
  - User provides a phenotype dataframe with three columns:
      gid_col  -- genotype / Name
      env_col  -- environment / location
      y_col    -- single trait
  - Rows where y is NaN are TEST rows (to be predicted).
  - Rows where y is observed are TRAIN rows (used to fit).
  - CV scheme is inferred per test row from the (gid, env) presence in train:
      CV1 : geno tested in training AND env tested in training (masked obs)
      CV2 : geno NOT tested, env tested
      CV0 : geno NOT tested AND env NOT tested   (new x new)
      CV0_env : geno tested, env NOT tested      (tested geno, new env)
  - User picks mode:
      mode="fast" -> output_level="predict_with_se", no AI-REML,
                     learn_scales=False. Point predictions + SE via the GP
                     posterior path only.
      mode="full" -> output_level="full_vc", AI-REML runs on the training
                     rows, returning a varcomp table alongside predictions.
                     REFUSED when any env has zero observed rows (CV0 env-
                     holdout) because FA psi/lambda for that env are
                     unidentifiable and predictions flip sign. Override with
                     allow_cv0_env_full=True if you understand the caveat.
  - backend: "auto" (default), "dense", or "operator"
      "operator" routes to the low-rank large-scale backend (n_obs >= ~50k).
      Requires grm_factor_cache (zarr or memmap dict), mode="fast", and
      method in {gp_exact, krr_exact}. Returns test-only predictions
      (train rows NOT included in predictions_df) and no Prediction_SE.
  - Return: (predictions_df, varcomp_df_or_None, info_dict).
    predictions_df has at minimum Name, Env, Prediction, cv_class. Dense path
    also includes Prediction_SE_* / Prediction_Var_* and train rows (with
    cv_class="train"). Operator path returns test rows only.

Post-fit helpers:
  - per_env_reliability(preds, varcomp) adds per-row `Var_g_env` and
    `Reliability_env = 1 - PEV/Var_g`.
  - across_env_predictions(preds, varcomp) returns a per-genotype summary:
    Prediction_across, PEV_across (independent-env approx), Var_g_across
    (exact (1/E^2) 1' Sigma_g 1 from varcomp), Reliability_across.
"""
from __future__ import annotations
import time
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _classify_cv(gid, env, tested_genos, tested_envs):
    g_in = gid in tested_genos
    e_in = env in tested_envs
    if g_in and e_in:
        return "CV1"
    if (not g_in) and e_in:
        return "CV2"
    if (not g_in) and (not e_in):
        return "CV0"
    return "CV0_env"  # geno tested, env untested — also a new-env scenario


def cv_predict(
    fw,  # loaded framework module
    pheno_df: pd.DataFrame,
    gid_col: str,
    env_col: str,
    y_col: str,
    geno_kernels: Dict[str, np.ndarray],
    geno_ids: List[str],
    *,
    method: str = "gp_exact",
    mode: str = "fast",
    env_similarity: Optional[np.ndarray] = None,
    env_structure: Optional[str] = None,
    fa_rank: int = 1,
    fixed_effects: Optional[List[str]] = None,
    standardize: str = "per_env",
    dtype: str = "float64",
    seed: int = 12345,
    ai_steps: int = 6,
    iters: int = 80,
    verbose: bool = True,
    **extra_kwargs,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Dict]:
    if mode not in ("fast", "full"):
        raise ValueError(f"mode must be 'fast' or 'full', got {mode!r}")

    allow_cv0_env_full = bool(extra_kwargs.pop("allow_cv0_env_full", False))
    backend = str(extra_kwargs.pop("backend", "auto")).lower()
    grm_factor_cache = extra_kwargs.pop("grm_factor_cache", None)
    if backend not in ("auto", "dense", "operator"):
        raise ValueError(f"backend must be 'auto'|'dense'|'operator', got {backend!r}")
    if backend == "operator":
        if mode != "fast":
            raise ValueError(
                "backend='operator' only supports mode='fast' (operator backend "
                "does not compute variance components; use mode='full' with "
                "backend='dense' instead)."
            )
        if method not in ("gp_exact", "krr_exact"):
            raise ValueError(
                f"backend='operator' only supports method in {{'gp_exact','krr_exact'}}; "
                f"got method={method!r}. For gp_icm_fa, use backend='dense'."
            )
        if grm_factor_cache is None:
            raise ValueError(
                "backend='operator' requires grm_factor_cache (dict with "
                "type='zarr'|'memmap' + cache paths). See test_operator_scale.py."
            )

    df = pheno_df.reset_index(drop=True).copy()
    for col in (gid_col, env_col, y_col):
        if col not in df.columns:
            raise KeyError(f"column {col!r} not in pheno_df")

    y_raw = pd.to_numeric(df[y_col], errors="coerce").to_numpy()
    train_mask = np.isfinite(y_raw)
    test_mask = ~train_mask
    train_idx = np.where(train_mask)[0].astype(np.int64)
    test_idx = np.where(test_mask)[0].astype(np.int64)

    if train_idx.size == 0:
        raise ValueError("no training rows (all y are NaN)")
    if test_idx.size == 0:
        raise ValueError("no test rows (all y observed) — nothing to predict")

    tested_genos = set(df.loc[train_mask, gid_col].astype(str).unique().tolist())
    tested_envs = set(df.loc[train_mask, env_col].astype(str).unique().tolist())

    all_envs = set(df[env_col].astype(str).unique().tolist())
    masked_envs = sorted(all_envs - tested_envs)
    auto_downgraded = False
    if masked_envs and mode == "full" and not allow_cv0_env_full:
        warnings.warn(
            f"cv_predict(mode='full') is structurally unidentifiable when env(s) "
            f"have zero observed rows (CV0 env-holdout). Masked env(s): "
            f"{masked_envs}. Auto-downgrading to mode='fast' (GP posterior "
            f"only, no AI-REML). Pass allow_cv0_env_full=True to run "
            f"mode='full' anyway and inspect the (unreliable) VC table for "
            f"the masked envs.",
            RuntimeWarning, stacklevel=2,
        )
        mode = "fast"
        auto_downgraded = True
    if masked_envs and mode == "full" and allow_cv0_env_full:
        warnings.warn(
            f"cv_predict(mode='full', allow_cv0_env_full=True) with masked "
            f"env(s) {masked_envs}: VC for these envs is unidentifiable; "
            f"predictions may be unreliable.",
            RuntimeWarning, stacklevel=2,
        )

    cv_class = np.array(["train"] * len(df), dtype=object)
    for i in test_idx:
        cv_class[i] = _classify_cv(
            str(df.at[i, gid_col]), str(df.at[i, env_col]),
            tested_genos, tested_envs,
        )

    # Replace NaN y with 0 so the framework's numeric path is stable;
    # y on test rows is unused during fitting (train_idx masks them out)
    # and unused at prediction time (posterior mean uses train y only).
    y_filled = np.where(train_mask, y_raw, 0.0)
    df_feed = df.copy()
    df_feed[y_col] = y_filled

    if backend == "operator":
        fit_kwargs = dict(
            pheno_df=df_feed, gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=geno_ids,
            env_similarity=env_similarity,
            include_components=("g", "e"),
            fixed_effects=fixed_effects,
            method=method,
            output_level="predict_only",
            point_predictions_only=True,
            prediction_output="test_only",
            backend="operator",
            grm_factor_cache=grm_factor_cache,
            train_idx=train_idx, test_idx=test_idx,
            dtype=dtype, standardize=standardize,
            seed=int(seed),
        )
        fit_kwargs.update(extra_kwargs)
    else:
        if mode == "fast":
            output_level = "predict_with_se"
            learn_scales = False
        else:
            output_level = "full_vc"
            learn_scales = bool(extra_kwargs.pop("learn_scales", True))
        fit_kwargs = dict(
            pheno_df=df_feed, gid_col=gid_col, env_col=env_col, y_col=y_col,
            geno_kernels=geno_kernels, geno_ids=geno_ids,
            env_similarity=env_similarity,
            include_components=("g", "ge", "e"),
            fixed_effects=fixed_effects,
            method=method,
            env_structure=env_structure,
            fa_rank=int(fa_rank),
            output_level=output_level,
            train_idx=train_idx, test_idx=test_idx,
            prediction_output="all",
            ai_hutch_samples=64, ai_steps=int(ai_steps), iters=int(iters),
            dtype=dtype, standardize=standardize,
            reml_normalize="asreml",
            learn_envdiag_noise=True,
            learn_scales=learn_scales,
            seed=int(seed),
        )
        fit_kwargs.update(extra_kwargs)

    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = fw.fit_mixed_model(**fit_kwargs)
    wall = time.perf_counter() - t0

    preds = res.get("result", {}).get("predictions")
    if preds is None or not hasattr(preds, "shape"):
        raise RuntimeError("engine returned no predictions")
    preds = preds.copy()

    # Normalize id column: operator backend uses "Genotype", dense uses "Name".
    if "Name" not in preds.columns and "Genotype" in preds.columns:
        preds = preds.rename(columns={"Genotype": "Name"})

    # Align cv_class onto preds by (Name, Env) lookup against df
    key_to_class = {}
    for i in range(len(df)):
        key = (str(df.at[i, gid_col]), str(df.at[i, env_col]))
        # if duplicates exist, later train rows override test assignment,
        # which is the right behavior: if (g,e) has any observed rep, the
        # cell is tested, not held out.
        if key not in key_to_class or cv_class[i] == "train":
            key_to_class[key] = cv_class[i]

    if "Name" in preds.columns and "Env" in preds.columns:
        keys = list(zip(preds["Name"].astype(str), preds["Env"].astype(str)))
        preds["cv_class"] = [key_to_class.get(k, "unknown") for k in keys]
    else:
        preds["cv_class"] = "unknown"

    vc = res.get("varcomp") if mode == "full" else None
    summary = res.get("summary") or {}

    # The framework suppresses its own VC-failure warning when callers wrap
    # in warnings.catch_warnings(simplefilter='ignore') (as cv_predict does
    # above for the fit call). Re-emit here outside the suppression so users
    # who request mode='full' but get no varcomp are not left guessing.
    if mode == "full" and vc is None:
        warnings.warn(
            f"cv_predict(mode='full', method={method!r}, env_structure={env_structure!r}) "
            f"received no varcomp from the framework. Reliability and across-env "
            f"helpers that need Sigma_g will degrade (return NaN for Var_g / "
            f"classical-reliability columns).",
            RuntimeWarning, stacklevel=2,
        )

    pred_meta = (res.get("diagnostics") or {}).get("prediction_meta") or {}
    info = {
        "mode": mode,
        "auto_downgraded": bool(auto_downgraded),
        "method": method,
        "backend": pred_meta.get("backend_used", backend),
        "n_train": int(train_idx.size),
        "n_test": int(test_idx.size),
        "cv_counts": {k: int((cv_class == k).sum()) for k in
                      ("train", "CV1", "CV2", "CV0", "CV0_env")},
        "wall_seconds": float(wall),
        "converged": bool(summary.get("converged", True)),
        "loglik": summary.get("loglik"),
        "n_ai_iter": summary.get("n_ai_iter"),
    }

    if verbose:
        print(f"[cv_predict] method={method!r} mode={mode!r} "
              f"n_train={info['n_train']} n_test={info['n_test']} "
              f"cv_counts={info['cv_counts']} wall={wall:.2f}s "
              f"converged={info['converged']} loglik={info['loglik']}")

    return preds, vc, info


def _parse_fa_sigma_g(varcomp: pd.DataFrame,
                      env_names: List[str]) -> Optional[np.ndarray]:
    """Reconstruct the genetic covariance Sigma_g (E x E) from an FA varcomp
    table. Sigma_g = Lambda Lambda' + diag(psi), where psi is the !var entry
    per env and Lambda[e,j] is the !{env}!fa{j} entry.

    Returns None if varcomp has no recognizable FA rows.
    """
    if varcomp is None or "component" not in varcomp.columns:
        return None
    comp = varcomp["component"].astype(str)
    est = pd.to_numeric(varcomp["estimate"], errors="coerce")
    E = len(env_names)

    # Detect fa_rank as the MAX j with any !fa{j} entry (scan without early
    # break so sparse/pruned FA patterns like {fa1, fa3} are reconstructed
    # at their actual rank rather than truncated at the first gap).
    fa_rank = 0
    for j in range(1, 64):
        if comp.str.contains(f"!fa{j}", regex=False).any():
            fa_rank = j
    if fa_rank == 0:
        # No FA loadings detected: fall back to diagonal psi only
        psi = np.zeros(E, dtype=float)
        any_hit = False
        for e, env in enumerate(env_names):
            rows = comp.str.contains(f"!{env}!var", regex=False) & \
                   comp.str.contains("fa(", regex=False)
            if rows.any():
                psi[e] = float(est[rows].iloc[0])
                any_hit = True
        if not any_hit:
            return None
        return np.diag(psi)

    Lam = np.zeros((E, fa_rank), dtype=float)
    psi = np.zeros(E, dtype=float)
    for e, env in enumerate(env_names):
        psi_rows = comp.str.contains(f"!{env}!var", regex=False) & \
                   comp.str.contains("fa(", regex=False)
        if psi_rows.any():
            psi[e] = float(est[psi_rows].iloc[0])
        for j in range(1, fa_rank + 1):
            lam_rows = comp.str.contains(f"!{env}!fa{j}", regex=False)
            if lam_rows.any():
                Lam[e, j - 1] = float(est[lam_rows].iloc[0])
    Sigma = Lam @ Lam.T + np.diag(psi)
    return Sigma


def _prior_latent_var_per_env(preds: pd.DataFrame,
                              env_col: str,
                              cv_class_col: str) -> Dict[str, float]:
    """Estimate the 'no-information' prior latent variance per env.

    For a test row where the genotype has no training connectivity, the
    posterior latent variance is near the prior. We use the per-env maximum
    over CV1/CV2 test rows as the empirical "no-info" anchor.

    CV0 and CV0_env rows are EXCLUDED from the anchor: those rows have no
    connected training data (new env and/or new geno), so their posterior
    can blow up to values far exceeding the prior in the standardized scale,
    which would inflate the denominator and make CV1/CV2 rows look
    artificially reliable. Excluding them preserves the monotonicity
    property that CV1 reliability >= CV2 reliability.

    Fallback ladder if no CV1/CV2 rows exist in an env:
      1) all test rows in that env (max)
      2) omit the env (no entry in the returned dict)

    This sidesteps scale mismatches between varcomp estimates and
    pred_lat.variance (per-env standardization, G diag != 1, ge/env-main
    random contributions) by deriving the denominator from the same scale
    as the numerator.
    """
    if "Prediction_Var_latent" not in preds.columns:
        return {}
    var_col = pd.to_numeric(preds["Prediction_Var_latent"], errors="coerce")
    envs = preds[env_col].astype(str)
    has_cv = cv_class_col in preds.columns
    cv = preds[cv_class_col].astype(str) if has_cv else None
    out: Dict[str, float] = {}
    for env in envs.unique():
        base = (envs == env) & np.isfinite(var_col)
        if has_cv:
            anchor = base & cv.isin(["CV1", "CV2"])
            if not anchor.any():
                # No connected test rows in this env: fall back to any test row
                anchor = base & (cv != "train")
        else:
            anchor = base
        if anchor.any():
            out[env] = float(var_col[anchor].max())
    return out


def per_env_reliability(preds: pd.DataFrame,
                        varcomp: Optional[pd.DataFrame] = None,
                        *,
                        env_col: str = "Env",
                        gid_col: str = "Name",
                        cv_class_col: str = "cv_class") -> pd.DataFrame:
    """Augment `preds` with per-env reliability columns.

    Two reliability flavors are computed when sources are available:

      (A) Classical BLUP: `Reliability_env_classical = 1 - PEV_g / Var_g`,
          where PEV_g = `SE_g_latent**2` (pure genetic BLUP variance, only
          produced by the krr_exact engine) and Var_g = diag(Sigma_g) from
          varcomp. Both on the same inv-standardized scale, so the ratio is
          meaningful.

      (B) Data-driven: `Reliability_env = 1 - Prediction_Var_latent / Var_latent_prior_env`,
          where `Var_latent_prior_env` is the per-env max posterior variance
          over test rows (== prior at a no-information reference row).
          Used for gp_icm_fa / gp_exact engines where the per-component
          breakdown is not exposed: Prediction_Var_latent mixes g/ge/e
          contributions and sits on a different scale from varcomp, so
          Var_g as a denominator produces spurious saturation at 0.

    Both are clipped to [0, 1]. Columns added:
      - Var_g_env                  (from varcomp, NaN if absent)
      - Var_latent_prior_env       (data-driven)
      - Reliability_env            (flavor B, data-driven)
      - Reliability_env_classical  (flavor A, only if SE_g_latent is present)
    """
    out = preds.copy()
    if env_col not in out.columns:
        raise KeyError(f"{env_col!r} not in preds")
    env_names = sorted(out[env_col].astype(str).unique().tolist())

    Sigma_g = _parse_fa_sigma_g(varcomp, env_names) if varcomp is not None else None
    if Sigma_g is not None:
        var_g_by_env = {env: float(Sigma_g[i, i]) for i, env in enumerate(env_names)}
        out["Var_g_env"] = out[env_col].astype(str).map(var_g_by_env)
    else:
        out["Var_g_env"] = np.nan

    prior_by_env = _prior_latent_var_per_env(out, env_col, cv_class_col)
    out["Var_latent_prior_env"] = out[env_col].astype(str).map(prior_by_env)

    # Flavor B: data-driven
    if "Prediction_Var_latent" in out.columns:
        pev_lat = pd.to_numeric(out["Prediction_Var_latent"], errors="coerce")
    elif "Prediction_SE_latent" in out.columns:
        pev_lat = pd.to_numeric(out["Prediction_SE_latent"], errors="coerce") ** 2
    else:
        pev_lat = pd.Series(np.nan, index=out.index)
    denom_b = pd.to_numeric(out["Var_latent_prior_env"], errors="coerce")
    with np.errstate(invalid="ignore", divide="ignore"):
        rel_b = 1.0 - (pev_lat / denom_b)
    rel_b[~np.isfinite(rel_b)] = np.nan
    out["Reliability_env"] = rel_b.clip(lower=0.0, upper=1.0)

    # Flavor A: classical BLUP (krr_exact only, when SE_g_latent present)
    if "SE_g_latent" in out.columns and Sigma_g is not None:
        pev_g = pd.to_numeric(out["SE_g_latent"], errors="coerce") ** 2
        denom_a = pd.to_numeric(out["Var_g_env"], errors="coerce")
        with np.errstate(invalid="ignore", divide="ignore"):
            rel_a = 1.0 - (pev_g / denom_a)
        rel_a[~np.isfinite(rel_a)] = np.nan
        out["Reliability_env_classical"] = rel_a.clip(lower=0.0, upper=1.0)

    return out


def across_env_predictions(preds: pd.DataFrame,
                           varcomp: Optional[pd.DataFrame] = None,
                           *,
                           env_col: str = "Env",
                           gid_col: str = "Name",
                           cv_class_col: str = "cv_class",
                           restrict_classes: Optional[List[str]] = None,
                           weight_col: Optional[str] = None,
                           posterior_corr: str = "independent",
                           ) -> pd.DataFrame:
    """Per-genotype across-environment summary.

    Computes, per genotype, with per-env weights w_e (default w_e = 1):
      - Prediction_across   = sum_e w_e y_e / sum_e w_e
      - PEV_across          ~= sum_e w_e^2 PostVar_e / (sum_e w_e)^2
      - Prediction_SE_across = sqrt(PEV_across)
      - Var_g_across        = w' Sigma_g w / (1' w)^2     (exact, from varcomp)
      - Prior_latent_across = sum_e w_e^2 Prior_latent_e / (sum_e w_e)^2
      - Reliability_across  = 1 - PEV_across / Prior_latent_across (clipped [0,1])

    Weights:
      - If `weight_col` is None (default), all envs for a genotype get equal
        weight 1 (back-compatible; formulas reduce to (1/E^2) sum_e ...).
      - If `weight_col` is a column name in `preds`, per-row weights are
        aggregated to per-env via mean (so replicate rows for the same (gid,
        env) produce ONE per-env weight). Typical uses: stage-1 precisions
        from two-stage analysis, or inverse observation variances.
      - Non-finite / non-positive weights are dropped from aggregation for
        that genotype.

    Notes:
      - Across-env SE uses an independent-env approximation (ignores posterior
        cross-env covariance). `posterior_corr="prior_corr"` uses a
        prior-correlation-scaled approximation instead (see below).
      - `Var_g_across` is the true genetic variance (from varcomp). Reported
        as reference but NOT used as the reliability denominator -- see
        per_env_reliability docstring for the scale-mismatch issue.
      - Reliability_across uses a data-driven prior (per-env max
        Prediction_Var_latent over test rows) for scale-consistency with PEV.
      - Partial-env genotypes: n_env < E_full; Var_g_across computed on the
        restricted (E_g x E_g) sub-block of Sigma_g.
      - `restrict_classes`: e.g. ["CV1","CV2"] -- only those rows are averaged.

    posterior_corr:
      - "independent" (default): PEV_across = (1/E^2) sum_e w_e^2 PostVar_e.
        Assumes zero posterior cross-env covariance. Always a LOWER bound on
        PEV_across when Σ_g has positive off-diagonals, so Reliability_across
        (and Reliability_across_classical) reported here are UPPER bounds.
      - "prior_corr": when SE_g_latent AND Σ_g (from varcomp) are both
        available, approximates posterior cross-env cov by scaling SEs with
        the PRIOR correlation from Σ_g:
           cov_post(g, e1, e2) ≈ SE_e1 * SE_e2 * Σ_g[e1,e2]/sqrt(Σ_g[e1,e1]*Σ_g[e2,e2])
        Exact when no training data shrinks the posterior covariance (CV2
        genotypes); close otherwise (training reduces posterior variance
        diagonals, but the correlation structure is approximately preserved
        when every env has similar training coverage). Emits two extra
        columns: `PEV_g_across_priorcorr` and
        `Reliability_across_classical_priorcorr`. Falls back to NaN when
        either SE_g_latent or Σ_g is missing.
    """
    if gid_col not in preds.columns or env_col not in preds.columns:
        raise KeyError(f"preds must contain {gid_col!r} and {env_col!r}")
    if "Prediction" not in preds.columns:
        raise KeyError("preds must contain 'Prediction' column")
    if posterior_corr not in ("independent", "prior_corr"):
        raise ValueError(
            f"posterior_corr must be 'independent' or 'prior_corr'; got "
            f"{posterior_corr!r}"
        )

    work = preds.copy()
    if restrict_classes is not None and cv_class_col in work.columns:
        work = work[work[cv_class_col].isin(list(restrict_classes))].copy()
    if work.empty:
        return pd.DataFrame(columns=[
            gid_col, "n_env", "Prediction_across", "PEV_across",
            "Prediction_SE_across", "Var_g_across", "Prior_latent_across",
            "Reliability_across",
        ])

    if "Prediction_Var_latent" in work.columns:
        work["_pev"] = pd.to_numeric(work["Prediction_Var_latent"], errors="coerce")
    elif "Prediction_SE_latent" in work.columns:
        work["_pev"] = pd.to_numeric(work["Prediction_SE_latent"], errors="coerce") ** 2
    else:
        work["_pev"] = np.nan

    if "SE_g_latent" in work.columns:
        work["_pev_g"] = pd.to_numeric(work["SE_g_latent"], errors="coerce") ** 2
    else:
        work["_pev_g"] = np.nan

    if weight_col is not None:
        if weight_col not in work.columns:
            raise KeyError(
                f"weight_col={weight_col!r} not found in preds columns"
            )
        work["_w"] = pd.to_numeric(work[weight_col], errors="coerce")
    else:
        work["_w"] = 1.0

    env_names = sorted(preds[env_col].astype(str).unique().tolist())
    Sigma_g = _parse_fa_sigma_g(varcomp, env_names) if varcomp is not None else None
    E_full = len(env_names)
    var_g_full = (float(Sigma_g.sum()) / (E_full * E_full)
                  if Sigma_g is not None else np.nan)
    prior_by_env = _prior_latent_var_per_env(preds, env_col, cv_class_col)

    rows = []
    for gid, sub in work.groupby(gid_col, sort=False):
        # Collapse replicate rows: one (gid, env) cell contributes ONE env,
        # not one-per-row. Reps would otherwise inflate E_g and shrink
        # PEV_across by the rep factor (see correctness review).
        sub[env_col] = sub[env_col].astype(str)
        per_env = sub.groupby(env_col, sort=True, as_index=False).agg(
            Prediction=("Prediction", "mean"),
            _pev=("_pev", "mean"),
            _pev_g=("_pev_g", "mean"),
            _w=("_w", "mean"),
        )
        # Drop envs with non-finite or non-positive weights
        w_raw = per_env["_w"].to_numpy(dtype=float)
        keep = np.isfinite(w_raw) & (w_raw > 0)
        per_env = per_env.loc[keep].reset_index(drop=True)
        sub_envs = per_env[env_col].tolist()
        E_g = len(sub_envs)
        if E_g == 0:
            continue
        w = per_env["_w"].to_numpy(dtype=float)
        w_sum = float(w.sum())
        w_sum_sq = w_sum * w_sum
        pred_across = float((w * per_env["Prediction"].to_numpy()).sum() / w_sum)

        pev_vals = per_env["_pev"].to_numpy(dtype=float)
        if np.all(np.isfinite(pev_vals)):
            pev_across = float(np.sum(w * w * pev_vals)) / w_sum_sq
        else:
            pev_across = np.nan

        pev_g_vals = per_env["_pev_g"].to_numpy(dtype=float)
        if np.all(np.isfinite(pev_g_vals)):
            pev_g_across = float(np.sum(w * w * pev_g_vals)) / w_sum_sq
        else:
            pev_g_across = np.nan

        if Sigma_g is not None:
            idx = [env_names.index(e) for e in sub_envs if e in env_names]
            if len(idx) == E_g:
                sub_sigma = Sigma_g[np.ix_(idx, idx)]
                var_g_across = float(w @ sub_sigma @ w) / w_sum_sq
            else:
                var_g_across = np.nan
        else:
            var_g_across = np.nan

        prior_vals = np.array([prior_by_env.get(e, np.nan) for e in sub_envs], dtype=float)
        if np.all(np.isfinite(prior_vals)) and prior_vals.size > 0:
            prior_latent_across = float(np.sum(w * w * prior_vals)) / w_sum_sq
        else:
            prior_latent_across = np.nan

        if (np.isfinite(prior_latent_across) and prior_latent_across > 0
                and np.isfinite(pev_across)):
            rel = 1.0 - pev_across / prior_latent_across
            rel = max(0.0, min(1.0, rel))
        else:
            rel = np.nan

        if (np.isfinite(pev_g_across) and np.isfinite(var_g_across)
                and var_g_across > 0):
            rel_classical = 1.0 - pev_g_across / var_g_across
            rel_classical = max(0.0, min(1.0, rel_classical))
        else:
            rel_classical = np.nan

        pev_g_pc = np.nan
        rel_classical_pc = np.nan
        if posterior_corr == "prior_corr":
            if (Sigma_g is not None and np.all(np.isfinite(pev_g_vals))
                    and pev_g_vals.size == E_g and E_g >= 1):
                idx_pc = [env_names.index(e) for e in sub_envs if e in env_names]
                if len(idx_pc) == E_g:
                    sub_sigma_pc = Sigma_g[np.ix_(idx_pc, idx_pc)]
                    diag_sig = np.diag(sub_sigma_pc)
                    if np.all(diag_sig > 0):
                        # Prior correlation from Sigma_g
                        inv_sd = 1.0 / np.sqrt(diag_sig)
                        rho = sub_sigma_pc * np.outer(inv_sd, inv_sd)
                        # Posterior cov approx: SE_e1 * SE_e2 * rho[e1,e2]
                        se_g = np.sqrt(np.clip(pev_g_vals, 0.0, None))
                        cov_post = np.outer(se_g, se_g) * rho
                        # Weighted quadratic form
                        pev_g_pc = float(w @ cov_post @ w) / w_sum_sq
                        if np.isfinite(var_g_across) and var_g_across > 0:
                            rel_classical_pc = 1.0 - pev_g_pc / var_g_across
                            rel_classical_pc = max(0.0, min(1.0, rel_classical_pc))

        row_out = {
            gid_col: gid,
            "n_env": int(E_g),
            "Prediction_across": pred_across,
            "PEV_across": pev_across,
            "Prediction_SE_across": float(np.sqrt(pev_across)) if np.isfinite(pev_across) else np.nan,
            "Var_g_across": var_g_across,
            "Prior_latent_across": prior_latent_across,
            "Reliability_across": rel,
            "PEV_g_across": pev_g_across,
            "Reliability_across_classical": rel_classical,
        }
        if posterior_corr == "prior_corr":
            row_out["PEV_g_across_priorcorr"] = pev_g_pc
            row_out["Reliability_across_classical_priorcorr"] = rel_classical_pc
        rows.append(row_out)
    return pd.DataFrame(rows)


def score_against_truth(preds: pd.DataFrame, truth: pd.DataFrame,
                        gid_col="Name", env_col="Env",
                        y_col="y_true") -> pd.DataFrame:
    """Optional helper for validation: join preds with ground-truth y and
    report RMSE, predictive correlation, and 95% CI coverage per cv_class.
    """
    merged = preds.merge(
        truth[[gid_col, env_col, y_col]].rename(columns={y_col: "y_true"}),
        on=[gid_col, env_col], how="inner",
    )
    out_rows = []
    for cls, sub in merged.groupby("cv_class"):
        if cls == "train":
            continue
        yt = sub["y_true"].to_numpy(dtype=float)
        yp = sub["Prediction"].to_numpy(dtype=float)
        se_all = (sub["Prediction_SE_observed"].to_numpy(dtype=float)
                  if "Prediction_SE_observed" in sub.columns else None)
        mask = np.isfinite(yt) & np.isfinite(yp)
        if mask.sum() == 0:
            continue
        yt, yp = yt[mask], yp[mask]
        resid = yt - yp
        rmse = float(np.sqrt(np.mean(resid ** 2)))
        r = float(np.corrcoef(yt, yp)[0, 1]) if yp.size > 1 else np.nan
        row = {"cv_class": cls, "n": int(mask.sum()), "rmse": rmse,
               "corr": r, "mean_abs_bias": float(np.mean(np.abs(resid)))}
        if se_all is not None:
            se = se_all[mask]
            se_ok = np.isfinite(se)
            if se_ok.any():
                lo, hi = yp[se_ok] - 1.96 * se[se_ok], yp[se_ok] + 1.96 * se[se_ok]
                row["cov95"] = float(np.mean((yt[se_ok] >= lo) & (yt[se_ok] <= hi)))
        out_rows.append(row)
    return pd.DataFrame(out_rows)
