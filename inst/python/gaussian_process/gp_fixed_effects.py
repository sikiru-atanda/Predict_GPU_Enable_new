"""Training-defined, estimable fixed-effect designs for multi-trait GP fits."""
import numpy as np
import pandas as pd
import warnings
from scipy.linalg import qr


def warn_integer_slope(series, name):
    if pd.api.types.is_integer_dtype(series.dtype) and series.nunique() >= 3:
        warnings.warn(f'Integer fixed-effect column {name!r} is fitted as a numeric slope; '
                      'convert it to a factor/category if its values identify groups.', RuntimeWarning)


def full_rank_design(X, train_mask, column_names=None):
    """Remove aliases without assigning arbitrary predictions to unestimable rows.

    Column scaling precedes rank detection so a covariate's units do not decide
    whether an intercept is retained. Test rows never determine the fitted rank.
    """
    X = np.asarray(X, dtype=np.float64)
    train = np.asarray(train_mask, dtype=bool)
    names = list(column_names) if column_names is not None else [f'column[{i}]' for i in range(X.shape[1])]
    if not np.isfinite(X).all():
        bad = [names[i] for i in np.flatnonzero(~np.isfinite(X).all(axis=0))]
        raise ValueError(f'GP fixed effects {bad} must be finite and non-missing; supply values or remove these columns')
    if not X.shape[1]:
        return X
    norms = np.linalg.norm(X[train], axis=0)
    scaled = X / np.where(norms > 0, norms, 1.)
    _, R, piv = qr(scaled[train], mode='economic', pivoting=True)
    diagonal = np.abs(np.diag(R))
    tolerance = max(scaled[train].shape) * np.finfo(float).eps * (diagonal.max() if diagonal.size else 0.)
    rank = int(np.sum(diagonal > tolerance))
    if rank >= int(train.sum()):
        raise ValueError(f'GP fixed effects {names} leave no residual degrees of freedom; remove saturated terms or add training observations')
    kept = np.sort(piv[:rank])
    reduced = scaled[:, kept]
    coefficients = np.linalg.lstsq(reduced[train], scaled[train], rcond=None)[0]
    mismatch = ~np.isclose(reduced @ coefficients, scaled, rtol=1e-9, atol=1e-10)
    if mismatch.any():
        bad = [names[i] for i in np.flatnonzero(mismatch.any(axis=0))]
        raise ValueError(f'GP prediction fixed effects {bad} are not estimable from training rows; '
                         'include those levels in training or remove the unsupported fixed-effect term')
    # Scaling is only for rank/estimability detection. Preserve the original
    # covariate units in the fitted coefficients returned to callers.
    return X[:, kept]


def trait_fixed_design(df, fixed_effects, trait_idx, train_mask, T, intercepts):
    """Numeric slopes and treatment-coded factors, expanded by trait."""
    train = np.asarray(train_mask, dtype=bool)
    if not len(fixed_effects):
        # These intercept bases have disjoint nonzero support on training rows
        # by construction. Avoid a large QR and alias-reconstruction matrix for
        # the common MET model with no additional covariates.
        X = np.asarray(intercepts, dtype=np.float64)
        norms = np.linalg.norm(X[train], axis=0)
        if (norms <= 0).any() or X.shape[1] >= int(train.sum()):
            raise ValueError('GP mean design is unidentifiable or leaves no residual degrees of freedom')
        return X
    blocks = [np.asarray(intercepts, dtype=np.float64)]
    names = [f'intercept[{i}]' for i in range(blocks[0].shape[1])]
    for name in fixed_effects:
        if name not in df:
            raise ValueError(f'fixed_effects col {name!r} not in pheno_df')
        series = df[name]
        if series.isna().any():
            raise ValueError(f'fixed effect {name!r} contains missing values')
        if pd.api.types.is_numeric_dtype(series.dtype) and not pd.api.types.is_bool_dtype(series.dtype):
            warn_integer_slope(series, name)
            columns = series.to_numpy(dtype=np.float64)[:, None]
        else:
            values = series.astype(str).to_numpy()
            levels = sorted(set(values[train]))
            if not set(values).issubset(levels):
                raise ValueError(f'fixed effect {name!r} has levels absent from training; predictions are not estimable. Include those levels in training or remove this fixed effect')
            columns = np.column_stack([values == level for level in levels[1:]]) if len(levels) > 1 else np.empty((len(df), 0))
        for column in columns.T:
            blocks.append(column[:, None] * np.eye(T)[trait_idx])
            names.extend([f'{name}:trait[{t}]' for t in range(T)])
    return full_rank_design(np.concatenate(blocks, axis=1), train, names)


def met_mean_design(env_idx, trait_idx, train_mask, T, trait_env_sds, *, compact=False):
    """Environment-trait intercepts, with a pooled mean for unobserved cells.

    For CV0/unobserved environment-trait cells, the existing pooled training
    trait mean is the prediction target. Its GLS correction is the sample-count
    weighted average of observed environment means, transformed to that row's
    response scale. No unidentifiable new environment intercept is fitted.
    """
    train = np.asarray(train_mask, dtype=bool)
    if compact:
        from gp_absorbed_means import CellMeanDesign
        return CellMeanDesign(env_idx, trait_idx, train, T, trait_env_sds)
    cells, counts = np.unique(np.column_stack([env_idx[train], trait_idx[train]]), axis=0, return_counts=True)
    X = np.zeros((len(env_idx), len(cells)))
    for j, (env, trait) in enumerate(cells):
        X[(env_idx == env) & (trait_idx == trait), j] = 1.
    missing = X.sum(axis=1) == 0
    for trait in range(T):
        rows = np.flatnonzero(missing & (trait_idx == trait))
        if not rows.size:
            continue
        cols = np.flatnonzero(cells[:, 1] == trait)
        if not cols.size:
            raise ValueError(f'trait index {trait} has no training observations')
        weights = counts[cols] / counts[cols].sum()
        X[np.ix_(rows, cols)] = (weights * trait_env_sds[cells[cols, 0], trait])[None, :] / trait_env_sds[env_idx[rows], trait, None]
    return X
