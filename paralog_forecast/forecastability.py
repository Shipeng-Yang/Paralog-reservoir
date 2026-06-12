"""Leave-one-group-out forecastability assay + learning curve + within-group oracle.

This is the manuscript's central instrument: it asks whether a class-level enrichment
(captured by paralog features) converts into *usable, portable* member-level prioritization.
A high AUROC with empty top-k is the signature finding ("enrichment without forecastability").
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

from .features import build_feature_matrix
from .metrics import forecast_metrics, tie_diagnostic
from .schema import coerce_binary_labels


def _new_model(model):
    if model is None:
        return LogisticRegression(max_iter=1000, class_weight="balanced", C=1.0)
    return clone(model)  # fresh estimator per fold; user-supplied models must not carry state across folds


def leave_one_group_out_forecast(
    df: pd.DataFrame, *, label_col="label", group_col="group", gene_col="gene_id",
    model=None, neg_cap=None, rng=None, return_predictions=True, **feat_kwargs,
):
    """Train on all groups but one, predict the held-out group; pool predictions and score.

    Returns (predictions_df, metrics_dict, tie_dict). `neg_cap` optionally subsamples training
    negatives for speed on large genomes (a faster approximation; check that the ranking is stable).
    """
    rng = np.random.default_rng(1) if rng is None else rng
    metric_kwargs = _pop_metric_kwargs(feat_kwargs)
    X, y, groups, names = build_feature_matrix(df, label_col=label_col, group_col=group_col, **feat_kwargs)
    pooled = np.full(len(df), np.nan)
    for held in pd.unique(groups):
        te = groups == held; tr = ~te
        ytr = y[tr]
        if ytr.sum() < 2:
            continue
        tr_idx = np.where(tr)[0]
        if neg_cap is not None:
            pos = tr_idx[y[tr_idx] == 1]; neg = tr_idx[y[tr_idx] == 0]
            if len(neg) > neg_cap:
                neg = rng.choice(neg, neg_cap, replace=False)
            tr_idx = np.concatenate([pos, neg])
        ytr_idx = y[tr_idx]
        if ytr_idx.sum() < 2 or (len(ytr_idx) - ytr_idx.sum()) < 2:
            continue  # need both classes in the training fold
        m = _new_model(model)
        m.fit(X[tr_idx], ytr_idx)
        pooled[te] = m.predict_proba(X[te])[:, 1]
    mask = ~np.isnan(pooled)
    met = forecast_metrics(y[mask], pooled[mask], **metric_kwargs)
    tie = tie_diagnostic(pooled[mask], y_true=y[mask])
    preds = None
    if return_predictions:
        preds = pd.DataFrame({gene_col: df[gene_col].values[mask], group_col: groups[mask],
                              label_col: y[mask], "score": pooled[mask]})
    return preds, met, tie


def learning_curve_loco(
    df: pd.DataFrame, *, positive_subsamples=(10, 25, 50, 100, 250, 500), repeats=20,
    label_col="label", group_col="group", model=None, neg_cap=None, random_state=1, **feat_kwargs,
):
    """Subsample training positives at increasing N (LOCO), to test whether top-k concentration
    improves with more labels (sparsity-limited) or plateaus (signal/feature-limited).

    `subsample_N` is a GLOBAL positive-label budget per repeat: one subset of N positives is drawn
    from all positives, then each LOCO fold trains on that subset minus the held-out group's positives
    (so the effective training-positive count per fold is N minus whatever falls in the held-out group).
    This makes N a single, interpretable label budget rather than a per-fold cap.

    Returns a long DataFrame: subsample_N, repeat, plus the prevalence-robust metrics.
    """
    rng = np.random.default_rng(random_state)
    metric_kwargs = _pop_metric_kwargs(feat_kwargs)
    X, y, groups, names = build_feature_matrix(df, label_col=label_col, group_col=group_col, **feat_kwargs)
    avail = int(y.sum())
    all_pos = np.where(y == 1)[0]
    rows = []
    for N in positive_subsamples:
        if N > avail:
            continue
        for rep in range(repeats):
            global_keep = set(rng.choice(all_pos, N, replace=False).tolist())  # global positive label budget
            pooled = np.full(len(df), np.nan)
            for held in pd.unique(groups):
                te = groups == held; tr = ~te
                tr_pos = np.array([i for i in np.where(tr & (y == 1))[0] if i in global_keep], dtype=int)
                if len(tr_pos) == 0:
                    continue
                ytr = np.zeros(len(df), dtype=int); ytr[tr_pos] = 1
                tr_idx = np.where(tr)[0]
                if neg_cap is not None:
                    neg = tr_idx[y[tr_idx] == 0]
                    if len(neg) > neg_cap:
                        neg = rng.choice(neg, neg_cap, replace=False)
                    tr_idx = np.concatenate([tr_pos, neg])
                m = _new_model(model)
                if ytr[tr_idx].sum() < 2:
                    continue
                m.fit(X[tr_idx], ytr[tr_idx])
                pooled[te] = m.predict_proba(X[te])[:, 1]
            mask = ~np.isnan(pooled)
            if mask.sum() == 0 or y[mask].sum() == 0:
                continue
            met = forecast_metrics(y[mask], pooled[mask], **metric_kwargs)
            met.update(subsample_N=N, available_N=avail, repeat=rep)
            rows.append(met)
    return pd.DataFrame(rows)


def oracle_within_group_cv(
    df: pd.DataFrame, *, group_col="group", label_col="label", n_splits=5,
    min_positives=30, model=None, random_state=1, **feat_kwargs,
):
    """Within-group K-fold CV (oracle: the group is known) for each group with enough positives.

    Separates 'no signal' from 'signal that does not transfer across groups'. Exploratory:
    folds share within-group structure, so treat as an upper bound, not deployable performance.
    """
    rows = []
    metric_kwargs = _pop_metric_kwargs(feat_kwargs)
    for g in pd.unique(df[group_col]):
        d = df[df[group_col] == g]
        npos = int(coerce_binary_labels(d, label_col).sum())  # strict: blanks raise, not silently -> 0
        if npos < min_positives:
            rows.append(dict(group=g, n_pos=npos, status=f"skipped (n_pos<{min_positives})"))
            continue
        X, y, _, _ = build_feature_matrix(d, label_col=label_col, group_col=group_col, **feat_kwargs)
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        pooled = np.full(len(y), np.nan)
        for tri, tei in skf.split(X, y):
            m = _new_model(model); m.fit(X[tri], y[tri])
            pooled[tei] = m.predict_proba(X[tei])[:, 1]
        mask = ~np.isnan(pooled)
        met = forecast_metrics(y[mask], pooled[mask], **metric_kwargs)
        met.update(group=g, status="ok")
        rows.append(met)
    return pd.DataFrame(rows)


def _pop_metric_kwargs(feat_kwargs):
    """Remove and return top_ks/top_fracs so they go to forecast_metrics, not build_feature_matrix."""
    return {k: feat_kwargs.pop(k) for k in ("top_ks", "top_fracs") if k in feat_kwargs}
