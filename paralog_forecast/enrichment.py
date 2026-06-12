"""Cluster-robust duplicate-class enrichment (the manuscript's WGD-vs-other harness).

Re-implements the locked enrichment estimand: among duplicates, does WGD status predict the
trait label, with orthogroup-cluster-robust standard errors? Plus a Fisher dup-vs-singleton
contrast (no covariate, because family size is collinear with singleton status).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact

from .schema import coerce_binary_labels


def cluster_robust_wgd_enrichment(
    df: pd.DataFrame,
    *,
    label_col: str = "label",
    wgd_col: str = "is_wgd",
    family_size_col: str = "family_size",
    cluster_col: str = "orthogroup",
    duplicate_col: str | None = "is_duplicate",
    adjust_log_length: bool = False,
    length_col: str = "gene_length",
):
    """logit(label) ~ is_wgd + log1p(family_size) [+ log1p(length)] among DUPLICATES,
    with cluster-robust SEs by `cluster_col`.

    Returns a dict: odds_ratio, ci_low, ci_high, beta, se, p_value, n_genes, n_positive,
    n_clusters, converged. If `duplicate_col` is present, the model is fit on duplicates only.
    """
    d = df if duplicate_col is None or duplicate_col not in df.columns else df[df[duplicate_col] == 1]
    y = coerce_binary_labels(d, label_col)  # strict: blanks/garbage raise, never silently -> 0
    out = dict(odds_ratio=np.nan, ci_low=np.nan, ci_high=np.nan, beta=np.nan, se=np.nan,
               p_value=np.nan, n_genes=int(len(d)), n_positive=int(y.sum()),
               n_clusters=int(d[cluster_col].nunique()), converged=False)
    if y.sum() < 5 or (len(y) - y.sum()) < 5 or d[wgd_col].nunique() < 2:
        out["note"] = "too few positives/negatives or WGD has no contrast"
        return out
    cols = {"is_wgd": pd.to_numeric(d[wgd_col], errors="coerce").fillna(0).astype(float).values,
            "log_fs": np.log1p(pd.to_numeric(d[family_size_col], errors="coerce").fillna(1).astype(float).values)}
    if adjust_log_length and length_col in d.columns:
        cols["log_len"] = np.log1p(pd.to_numeric(d[length_col], errors="coerce").fillna(0).astype(float).values)
    X = sm.add_constant(pd.DataFrame(cols, index=d.index))
    try:
        res = sm.Logit(y, X).fit(disp=0, maxiter=200, cov_type="cluster",
                                 cov_kwds={"groups": d[cluster_col].values})
        b, se = float(res.params.iloc[1]), float(res.bse.iloc[1])
        out.update(odds_ratio=float(np.exp(b)), ci_low=float(np.exp(b - 1.96 * se)),
                   ci_high=float(np.exp(b + 1.96 * se)), beta=b, se=se,
                   p_value=float(res.pvalues.iloc[1]),
                   converged=bool(res.mle_retvals.get("converged", True)))
    except Exception as e:  # separation / singular design
        out["note"] = f"fit failed: {type(e).__name__}"
    return out


def fisher_dup_vs_singleton(
    df: pd.DataFrame,
    *,
    label_col: str = "label",
    duplicate_col: str = "is_duplicate",
    singleton_col: str = "is_singleton",
):
    """2x2 Fisher exact for label x (duplicate vs singleton). Returns odds_ratio, p_value, table."""
    e = coerce_binary_labels(df, label_col).astype(bool)  # strict labels
    dup = pd.to_numeric(df[duplicate_col], errors="coerce").fillna(0).astype(bool)
    sg = pd.to_numeric(df[singleton_col], errors="coerce").fillna(0).astype(bool)
    a, b = int((e & dup).sum()), int((~e & dup).sum())
    c, d = int((e & sg).sum()), int((~e & sg).sum())
    if (a + b) == 0 or (c + d) == 0:
        return dict(odds_ratio=np.nan, p_value=np.nan, table=[[a, b], [c, d]])
    orr, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return dict(odds_ratio=float(orr), p_value=float(p), table=[[a, b], [c, d]])
