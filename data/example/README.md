# Synthetic example

`synthetic_gene_table.tsv` is **simulated** (fixed seed; `scripts/make_synthetic_example.py`).
It contains NO manuscript data and NO real biological results. It exists only to exercise the
`paralog_forecast` API and show the expected output formats. It is built so that the assay shows a
real class-level WGD enrichment (large, significant odds ratio) alongside only **weak, non-portable
top-k retrieval** (high-ish AUROC but no reliable top-k concentration / flat learning curve) — the
paper's qualitative signature — without using any real trait-gene labels. Note the synthetic top-k is
weak/noisy rather than literally empty; the example illustrates the API and the *pattern*, and is not a
substitute for the manuscript's full-data results.
