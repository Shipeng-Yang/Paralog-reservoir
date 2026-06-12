"""paralog_forecast — a lightweight, reusable assay for the paper
"Crop trait genes are drawn from a functionally structured but non-forecastable paralog reservoir".

It implements, on any user-supplied gene-level table (one row per gene):
  1. cluster-robust WGD-vs-other duplicate-class enrichment        (enrichment.py)
  2. leave-one-group-out ParalogPrior forecastability              (forecastability.py)
  3. learning curve + within-group oracle                          (forecastability.py)
  4. prevalence-robust top-k / enrichment@k metrics + tie diagnostic (metrics.py)

This is NOT a full software package — it is a few plain modules. See README for the input schema.
"""
from .schema import validate_gene_table, SchemaError, coerce_binary_labels
from .features import build_feature_matrix
from .enrichment import cluster_robust_wgd_enrichment, fisher_dup_vs_singleton
from .metrics import forecast_metrics, tie_diagnostic
from .forecastability import (
    leave_one_group_out_forecast,
    learning_curve_loco,
    oracle_within_group_cv,
)

__all__ = [
    "validate_gene_table", "SchemaError", "coerce_binary_labels", "build_feature_matrix",
    "cluster_robust_wgd_enrichment", "fisher_dup_vs_singleton",
    "forecast_metrics", "tie_diagnostic",
    "leave_one_group_out_forecast", "learning_curve_loco", "oracle_within_group_cv",
]
__version__ = "0.1.0"
