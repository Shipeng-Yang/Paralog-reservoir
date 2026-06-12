#!/usr/bin/env python
"""Reproduce per-crop WGD-vs-other enrichment (Results Section 2/3) on the frozen foundation,
using the reusable paralog_forecast.cluster_robust_wgd_enrichment. Requires Zenodo frozen_inputs.

Run: PYTHONPATH=$PWD python reproduce/section_3_enrichment/run_enrichment.py
"""
import os, sys
import pandas as pd
sys.path.insert(0, os.getcwd())
from reproduce.common.paths import frozen
import paralog_forecast as pf

f = pd.read_parquet(frozen("gene_level_foundation.parquet"))
f = f[f["include_in_main"] == 1].copy()
for sp in ["Os", "Zm", "Sl", "Gm", "Sb", "At"]:
    d = f[f["species"] == sp]
    df = pd.DataFrame({
        "label": (d["cloned_domestication_flag"] == True).astype(int),
        "is_wgd": d["is_wgd"].astype(int), "is_duplicate": d["is_duplicate"].astype(int),
        "family_size": d["family_size"], "orthogroup": d["orthogroup_id"],
    })
    r = pf.cluster_robust_wgd_enrichment(df)
    print(f"{sp}: OR={r['odds_ratio']:.2f} CI=[{r['ci_low']:.2f},{r['ci_high']:.2f}] "
          f"p={r['p_value']:.3g} pos={r['n_positive']}")
