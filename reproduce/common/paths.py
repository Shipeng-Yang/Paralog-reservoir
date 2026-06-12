"""Resolve the data root for reproduction scripts (no hard-coded NAS paths).

Priority: (1) PARALOG_ROOT env var, (2) config/config.yaml project_root, else a clear error.
This replaces the original project-local absolute path so a cloned repo can run once the
Zenodo frozen_inputs/ archive is downloaded and project_root points at it.
"""
import os


def get_project_root(config_path="config/config.yaml"):
    env = os.environ.get("PARALOG_ROOT")
    if env:
        return env
    if os.path.exists(config_path):
        try:
            import yaml
            with open(config_path) as fh:
                cfg = yaml.safe_load(fh) or {}
            if cfg.get("project_root"):
                return cfg["project_root"]
        except Exception:
            pass
    raise SystemExit(
        "Could not resolve the data root. Set the PARALOG_ROOT environment variable, e.g.\n"
        "    export PARALOG_ROOT=/path/to/paralog_frozen_inputs\n"
        "or copy config/config.example.yaml to config/config.yaml and set project_root.\n"
        "The frozen inputs themselves are on Zenodo (see README, Data availability)."
    )


def frozen(name, root=None):
    """Path to a frozen input by filename, relative to the resolved data root."""
    return os.path.join(root or get_project_root(), name)
