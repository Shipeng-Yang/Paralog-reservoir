"""Prevalence-robust retrieval metrics + tie diagnostic for the forecastability assay.

recall@k is deliberately NOT the headline metric: its denominator grows with the positive
count and it is deflated under extreme class imbalance. The honest metrics are hits@k,
enrichment@k = hits / (k * prevalence), top-fraction enrichment, AUPRC-lift, and median rank
percentile, read alongside a tie diagnostic (discrete features tie large blocks at the top).
"""
from __future__ import annotations
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score


def forecast_metrics(y_true, y_score, *, top_fracs=(0.01, 0.05), top_ks=(10, 25, 50, 100)):
    """Return a flat dict of prevalence-robust ranking metrics."""
    y = np.asarray(y_true); s = np.asarray(y_score, dtype=float)
    n = len(y); npos = int(y.sum()); prev = npos / n if n else np.nan
    order = np.argsort(-s); yo = y[order]
    out = {"n": n, "n_pos": npos, "prevalence": prev}
    out["auroc"] = float(roc_auc_score(y, s)) if 0 < npos < n else np.nan
    out["auprc"] = float(average_precision_score(y, s)) if 0 < npos < n else np.nan
    out["auprc_lift"] = out["auprc"] / prev if (prev and np.isfinite(out["auprc"])) else np.nan
    ranks = np.empty(n); ranks[order] = np.arange(n)
    # rank percentile of positives; 0.0 = top-ranked, ~0.5 = random, 1.0 = bottom (divide by n-1 so the
    # worst possible rank maps to exactly 1.0)
    out["median_rank_pct_pos"] = float(np.median(ranks[y == 1]) / max(n - 1, 1)) if npos else np.nan
    for k in top_ks:
        if k > n:
            out[f"hits@{k}"] = np.nan; out[f"enrichment@{k}"] = np.nan; out[f"recall@{k}"] = np.nan
            continue
        hits = int(yo[:k].sum()); exp = k * prev
        out[f"hits@{k}"] = hits
        out[f"enrichment@{k}"] = (hits / exp) if exp > 0 else np.nan
        out[f"recall@{k}"] = hits / npos if npos else np.nan  # descriptive only
    for frac in top_fracs:
        kk = max(1, int(round(frac * n))); hits = int(yo[:kk].sum()); exp = kk * prev
        tag = f"top_{frac:g}"
        out[f"{tag}_hits"] = hits
        out[f"{tag}_enrichment"] = (hits / exp) if exp > 0 else np.nan
    return out


def tie_diagnostic(y_score, *, y_true=None, top_frac=0.01, round_dp=6):
    """Quantify discrete-feature score ties at the top of the ranking.

    Returns n, n_distinct_scores, the max-score block size, and the top-frac cutoff block size;
    if y_true is supplied, also the number of positives in each block. A large max-score block
    means the literal top-k is an arbitrary tie-break, so top-k counts must be read with this.
    """
    s = np.round(np.asarray(y_score, dtype=float), round_dp)
    n = len(s)
    y = None if y_true is None else np.asarray(y_true)
    mx = s.max() if n else np.nan
    max_mask = (s == mx) if n else np.zeros(0, dtype=bool)
    srt = np.sort(s)[::-1]
    kk = max(1, int(round(top_frac * n))) if n else 1
    cut = srt[kk - 1] if n else np.nan
    cut_mask = (s >= cut - 1e-9) if n else np.zeros(0, dtype=bool)
    out = dict(n=n, n_distinct_scores=int(len(np.unique(s))) if n else 0,
               max_score_block=int(max_mask.sum()), topfrac=top_frac,
               topfrac_cutoff_block=int(cut_mask.sum()))
    if y is not None:
        out["max_score_block_pos"] = int(y[max_mask].sum())
        out["topfrac_cutoff_block_pos"] = int(y[cut_mask].sum())
    return out
