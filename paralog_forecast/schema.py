"""Input-table validation for paralog_forecast.

The reusable assay operates on a single gene-level table (one row per gene). All column
names are caller-supplied so the module never depends on project-specific names.
"""
from __future__ import annotations
import pandas as pd


class SchemaError(ValueError):
    """Raised when the supplied gene table does not meet the required schema."""


def coerce_binary_labels(df: pd.DataFrame, label_col: str = "label"):
    """Strictly coerce a label column to a float 0/1 numpy array, raising SchemaError on any
    missing/non-numeric or non-binary value. Use this in every public entry point so that blanks
    or garbage are NEVER silently treated as negatives (unknown genes must be coded explicit 0)."""
    lab = pd.to_numeric(df[label_col], errors="coerce")
    if lab.isna().any():
        raise SchemaError(f"{label_col} has missing/non-numeric values; labels must be explicit 0/1 "
                          "(unknown genes must be coded 0, not left blank)")
    if not lab.isin((0, 1)).all():
        raise SchemaError(f"{label_col} must be binary 0/1")
    return lab.astype(float).values


def validate_gene_table(
    df: pd.DataFrame,
    *,
    gene_col: str = "gene_id",
    group_col: str = "group",
    label_col: str = "label",
    cluster_col: str = "orthogroup",
    required_feature_cols=("is_wgd", "family_size"),
) -> None:
    """Validate the minimal schema for the forecastability / enrichment assay.

    Required: a unique gene id, a grouping column (crop/species, for leave-one-group-out),
    a binary 0/1 label, a cluster id (e.g. orthogroup) for cluster-robust SEs, and the
    required feature columns. Raises SchemaError with an actionable message on failure.
    """
    missing = [c for c in (gene_col, group_col, label_col, cluster_col, *required_feature_cols)
               if c not in df.columns]
    if missing:
        raise SchemaError(f"missing required columns: {missing}. Present: {list(df.columns)}")
    if df[gene_col].duplicated().any():
        n = int(df[gene_col].duplicated().sum())
        raise SchemaError(f"{gene_col} must be unique ({n} duplicated gene ids)")
    lab = coerce_binary_labels(df, label_col)  # strict 0/1; raises on missing/non-numeric/non-binary
    if int(lab.sum()) < 2:
        raise SchemaError(f"{label_col} has < 2 positives; nothing to learn/test")
    if df[group_col].nunique() < 2:
        raise SchemaError(f"{group_col} must have >= 2 groups for leave-one-group-out")
    if df[cluster_col].isna().any():
        raise SchemaError(f"{cluster_col} (cluster id) has missing values")
