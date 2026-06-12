"""Feature-matrix construction from a paralog gene table.

Generic re-implementation of the manuscript's ParalogPrior features (WP7/WP18.5): numeric
predictors + one-hot-encoded categorical predictors + log1p(family size). Uses manual
get_dummies (not patsy) to stay robust under numpy>=2 / pandas StringDtype.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .schema import coerce_binary_labels


def build_feature_matrix(
    df: pd.DataFrame,
    *,
    label_col: str = "label",
    group_col: str = "group",
    numeric_cols=("is_wgd", "family_size", "age"),
    categorical_cols=("dup_mode", "function_class"),
    log_cols=("family_size",),
    feature_cols=None,
):
    """Return (X, y, groups, feature_names).

    - numeric_cols present in df are used as-is (coerced to float; NaN->0).
    - log_cols are added as log1p(col) and the raw col is dropped from the numeric set.
    - categorical_cols present in df are one-hot encoded (drop NaN category).
    - feature_cols, if given, overrides the auto-selected feature set (must already exist).
    """
    parts = []
    names = []
    log_set = set(log_cols)
    if feature_cols is not None:
        X = df[list(feature_cols)].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        return X.values, _y(df, label_col), df[group_col].values, list(feature_cols)
    for c in numeric_cols:
        if c not in df.columns:
            continue
        v = pd.to_numeric(df[c], errors="coerce").fillna(0.0).astype(float)
        if c in log_set:
            parts.append(np.log1p(v.values)); names.append(f"log1p_{c}")
        else:
            parts.append(v.values); names.append(c)
    for c in categorical_cols:
        if c not in df.columns:
            continue
        d = pd.get_dummies(df[c].astype("object"), prefix=c, dummy_na=False)
        for col in d.columns:
            parts.append(d[col].astype(float).values); names.append(col)
    if not parts:
        raise ValueError("no usable feature columns found; supply numeric_cols/categorical_cols "
                         "that exist in the table, or pass feature_cols explicitly")
    X = np.column_stack(parts)
    return X, _y(df, label_col), df[group_col].values, names


def _y(df, label_col):
    return coerce_binary_labels(df, label_col).astype(int)  # strict: blanks/garbage raise, never silently -> 0
