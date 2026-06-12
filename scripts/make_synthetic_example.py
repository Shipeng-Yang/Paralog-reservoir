#!/usr/bin/env python
"""Generate a small SYNTHETIC gene table that exercises paralog_forecast end-to-end.

This contains NO manuscript data and NO real biological results — it is simulated only, with a
fixed seed, so users can run the API in seconds and see the expected output formats. It is built
to reproduce the paper's *qualitative* signature (a real class-level enrichment that does NOT
convert into usable top-k forecastability), without using any real labels.

Run: python scripts/make_synthetic_example.py
Writes: data/example/synthetic_gene_table.tsv
"""
import os
import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(__file__), "..", "data", "example", "synthetic_gene_table.tsv")
rng = np.random.default_rng(20260608)

N_GROUPS = 5
GENES_PER_GROUP = 1500
DUP_MODES = ["singleton", "tandem", "proximal", "transposed", "dispersed", "wgd"]
FUNC_CLASSES = ["TF", "hormone", "kinase", "transporter", "enzyme", "defense", "other"]

rows = []
gid = 0
for g in range(N_GROUPS):
    group = f"crop{g+1}"
    # duplication mode composition (WGD ~ 25-40%, singletons ~10%)
    mode = rng.choice(DUP_MODES, GENES_PER_GROUP, p=[0.10, 0.13, 0.06, 0.22, 0.20, 0.29])
    is_wgd = (mode == "wgd").astype(int)
    is_singleton = (mode == "singleton").astype(int)
    is_duplicate = (1 - is_singleton).astype(int)
    family_size = np.where(is_singleton == 1, 1, rng.integers(2, 12, GENES_PER_GROUP))
    # function class: WGD enriched for TF/hormone; tandem/proximal for enzyme/defense
    func = []
    for m in mode:
        if m == "wgd":
            func.append(rng.choice(FUNC_CLASSES, p=[0.30, 0.18, 0.10, 0.08, 0.10, 0.06, 0.18]))
        elif m in ("tandem", "proximal"):
            func.append(rng.choice(FUNC_CLASSES, p=[0.06, 0.05, 0.08, 0.08, 0.30, 0.25, 0.18]))
        else:
            func.append(rng.choice(FUNC_CLASSES, p=[0.12, 0.10, 0.12, 0.12, 0.18, 0.12, 0.24]))
    func = np.array(func)
    age = rng.integers(1, 5, GENES_PER_GROUP)  # ordinal phylostratum proxy
    og = rng.integers(0, 400, GENES_PER_GROUP) + g * 1000  # orthogroup cluster id
    # latent trait-recruitment propensity: real class-level signal (WGD + TF) but with
    # large group-specific and idiosyncratic noise -> enrichment WITHOUT portable top-k
    lin = (-5.0 + 0.9 * is_wgd + 0.8 * (func == "TF") + 0.5 * (func == "hormone")
           + 0.15 * np.log1p(family_size) + rng.normal(0, 1.4, GENES_PER_GROUP)  # idiosyncratic noise
           + rng.normal(0, 0.6))  # group intercept shift
    p = 1 / (1 + np.exp(-lin))
    label = rng.binomial(1, p)
    for i in range(GENES_PER_GROUP):
        rows.append(dict(gene_id=f"g{gid:06d}", group=group, label=int(label[i]),
                         orthogroup=f"OG{og[i]:05d}", is_wgd=int(is_wgd[i]),
                         is_duplicate=int(is_duplicate[i]), is_singleton=int(is_singleton[i]),
                         family_size=int(family_size[i]), dup_mode=mode[i],
                         function_class=func[i], age=int(age[i])))
        gid += 1

df = pd.DataFrame(rows)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
df.to_csv(OUT, sep="\t", index=False)
print(f"wrote {os.path.relpath(OUT)}: {len(df)} genes, {int(df.label.sum())} positives "
      f"({df.label.mean():.3%} prevalence), {df.group.nunique()} groups")
