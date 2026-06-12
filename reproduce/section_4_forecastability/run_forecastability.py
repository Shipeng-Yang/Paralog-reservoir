#!/usr/bin/env python
"""Reproduce the manuscript forecastability result (Results Section 4 / WP18.5) using the
reusable paralog_forecast module on the project's frozen gene-level foundation.

Requires the Zenodo frozen_inputs archive (gene_level_foundation.parquet). Resolve its location
via PARALOG_ROOT or config/config.yaml (see reproduce/common/paths.py). The reusable assay is the
same code shipped in paralog_forecast/, so this doubles as a worked example on real data.

Run: PYTHONPATH=$PWD python reproduce/section_4_forecastability/run_forecastability.py
"""
import os, sys
import pandas as pd
sys.path.insert(0, os.getcwd())
from reproduce.common.paths import frozen
import paralog_forecast as pf

f = pd.read_parquet(frozen("gene_level_foundation.parquet"))
f = f[f["include_in_main"] == 1].copy()
# map the project's columns onto the generic assay schema
df = pd.DataFrame({
    "gene_id": f["gene_id"], "group": f["species"],
    "label": (f["cloned_domestication_flag"] == True).astype(int),
    "orthogroup": f["orthogroup_id"], "is_wgd": f["is_wgd"].astype(int),
    "is_duplicate": f["is_duplicate"].astype(int), "is_singleton": f["is_singleton"].astype(int),
    "family_size": f["family_size"], "dup_mode": f["dup_class6"],
    "function_class": f["primary_functional_bin"], "age": f["age_tier"],
})
# restrict to the crops that carry cloned labels (leave-one-crop-out is over those groups)
df = df[df["group"].isin(["At", "Gm", "Os", "Sb", "Sl", "Zm"])]
pf.validate_gene_table(df)

print("== cloned-endpoint enrichment (WGD-vs-other, cluster-robust) ==")
print(pf.cluster_robust_wgd_enrichment(df))
print("\n== leave-one-crop-out forecastability ==")
preds, met, tie = pf.leave_one_group_out_forecast(df, neg_cap=30000)
for k in ("auroc", "auprc_lift", "hits@100", "enrichment@100", "top_0.01_enrichment",
          "top_0.05_enrichment", "median_rank_pct_pos"):
    print(f"  {k} = {met.get(k)}")
print(f"  tie: distinct_scores={tie['n_distinct_scores']} max_score_block={tie['max_score_block']}")
print("\nInterpretation: real aggregate signal (AUROC ~0.7, top-5% enrichment) with empty/weak "
      "operational top-k = enrichment without usable member-level forecastability.")
