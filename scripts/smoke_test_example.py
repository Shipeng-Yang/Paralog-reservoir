#!/usr/bin/env python
"""End-to-end smoke test of paralog_forecast on the synthetic example table.

Exercises every public entry point and prints a compact report. Exits non-zero on failure.
Run: PYTHONPATH=$PWD python scripts/smoke_test_example.py
"""
import os, sys
import pandas as pd

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
import paralog_forecast as pf  # noqa: E402

EX = os.path.join(HERE, "..", "data", "example", "synthetic_gene_table.tsv")
if not os.path.exists(EX):
    print("synthetic example missing; run scripts/make_synthetic_example.py first"); sys.exit(2)
df = pd.read_csv(EX, sep="\t")

print("== validate ==")
pf.validate_gene_table(df)
print(f"  OK: {len(df)} genes, {int(df.label.sum())} positives, {df.group.nunique()} groups")

print("== cluster-robust WGD-vs-other enrichment (among duplicates) ==")
enr = pf.cluster_robust_wgd_enrichment(df)
print(f"  OR={enr['odds_ratio']:.2f} CI=[{enr['ci_low']:.2f},{enr['ci_high']:.2f}] "
      f"p={enr['p_value']:.2e} n={enr['n_genes']} pos={enr['n_positive']}")
assert enr["converged"], "enrichment model did not converge"

print("== dup-vs-singleton Fisher ==")
fis = pf.fisher_dup_vs_singleton(df)
print(f"  OR={fis['odds_ratio']:.2f} p={fis['p_value']:.2e}")

print("== leave-one-group-out forecastability ==")
preds, met, tie = pf.leave_one_group_out_forecast(df)
print(f"  AUROC={met['auroc']:.3f} AUPRC_lift={met['auprc_lift']:.2f} "
      f"hits@100={met['hits@100']} enrichment@100={met['enrichment@100']:.2f} "
      f"top_0.01_enrichment={met['top_0.01_enrichment']:.2f} median_rank_pct={met['median_rank_pct_pos']:.3f}")
print(f"  tie diagnostic: distinct_scores={tie['n_distinct_scores']} max_score_block={tie['max_score_block']}")
assert preds is not None and len(preds) == len(df)

print("== learning curve (sparsity vs signal) ==")
lc = pf.learning_curve_loco(df, positive_subsamples=(25, 50, 100, 200), repeats=5)
summ = lc.groupby("subsample_N")[["auroc", "enrichment@100", "top_0.01_enrichment"]].median()
print(summ.to_string())

print("== within-group oracle ==")
orc = pf.oracle_within_group_cv(df, min_positives=20, n_splits=4)
print(orc[[c for c in ("group", "n_pos", "auroc", "enrichment@100", "status") if c in orc.columns]].to_string(index=False))

print("\nSMOKE TEST PASSED")
