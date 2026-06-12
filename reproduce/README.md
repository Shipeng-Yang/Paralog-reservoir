# reproduce/ — manuscript reproduction

These scripts regenerate the paper's enrichment and forecastability results **from the frozen
inputs**, which are archived on Zenodo (not in this repo; see top-level README → Data availability).

Setup:
1. `conda env create -f ../environment.yml && conda activate paralog-forecast`
2. Download the Zenodo `frozen_inputs` archive.
3. `export PARALOG_ROOT=/path/to/unpacked/frozen_inputs`  (or set `config/config.yaml`).
4. Run a section, e.g. `PYTHONPATH=$PWD python reproduce/section_4_forecastability/run_forecastability.py`.

The drivers call the same `paralog_forecast/` module that is documented for reuse on your own data,
so they double as worked real-data examples. The full per-figure/per-table scripts used in the paper
accompany the Zenodo archive; these section drivers reproduce the load-bearing numbers.
