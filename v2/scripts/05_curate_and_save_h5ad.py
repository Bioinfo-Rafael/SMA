#!/usr/bin/env python3
"""05_curate_and_save_h5ad.py -- apply the human-edited curation yaml and write
curated h5ad files.

Per dataset:
  * rename_obs_columns
  * add_constant_obs
  * cell_type_column -> obs['cell_type']  (else 'unknown')
  * cluster_column   -> obs['cluster_label'] (else skipped)
  * original obs columns are NEVER dropped
  * var_names := gene_symbol (make_unique); gene_symbol_upper (re)built
  * X coerced to CSR sparse
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
from reporting import get_logger  # noqa: E402

log = get_logger("05_curate")


def curate_one(ds_id: str, cfg: dict, overwrite: bool):
    in_path = (ROOT / cfg["input_h5ad"]) if not Path(cfg["input_h5ad"]).is_absolute() \
        else Path(cfg["input_h5ad"])
    out_path = (ROOT / cfg["output_h5ad"]) if not Path(cfg["output_h5ad"]).is_absolute() \
        else Path(cfg["output_h5ad"])

    if not in_path.exists():
        raise RuntimeError(f"input h5ad not found: {in_path}")
    if out_path.exists() and not overwrite:
        log.info("[%s] curated exists, skip (%s)", ds_id, out_path.name)
        return out_path

    adata = ad.read_h5ad(in_path)

    rename = cfg.get("rename_obs_columns") or {}
    if rename:
        adata.obs = adata.obs.rename(columns=rename)
        log.info("[%s] renamed obs columns: %s", ds_id, rename)

    for key, value in (cfg.get("add_constant_obs") or {}).items():
        adata.obs[key] = value

    # cell_type (never auto-decided; only what the human specified)
    ctc = cfg.get("cell_type_column")
    if ctc and ctc in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs[ctc].astype(str)
        log.info("[%s] cell_type <- '%s'", ds_id, ctc)
    else:
        if ctc:
            log.warning("[%s] cell_type_column '%s' not in obs; using 'unknown'", ds_id, ctc)
        adata.obs["cell_type"] = "unknown"

    clc = cfg.get("cluster_column")
    if clc and clc in adata.obs.columns:
        adata.obs["cluster_label"] = adata.obs[clc].astype(str)
        log.info("[%s] cluster_label <- '%s'", ds_id, clc)
    elif clc:
        log.warning("[%s] cluster_column '%s' not in obs; skipping", ds_id, clc)

    # var_names := gene_symbol preferred, made unique
    if "gene_symbol" in adata.var.columns:
        adata.var_names = pd.Index(adata.var["gene_symbol"].astype(str))
        adata.var_names_make_unique()
    adata.var["gene_symbol_upper"] = adata.var_names.astype(str).str.upper()

    au.to_csr(adata)
    adata.uns["curation"] = {"cell_type_column": str(ctc),
                             "cluster_column": str(clc),
                             "keep_original_obs": bool(cfg.get("keep_original_obs", True))}

    au.save_h5ad(adata, out_path, overwrite=overwrite)
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curation", default=str(ROOT / "config" / "curation_template.yaml"))
    ap.add_argument("--datasets", nargs="*", help="dataset_id filter")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    paths = au.project_paths(ROOT)
    au.ensure_dirs(paths)
    with open(args.curation) as fh:
        curation = yaml.safe_load(fh)

    ok, failures = [], []
    for ds_id, cfg in (curation.get("datasets") or {}).items():
        if args.datasets and ds_id not in args.datasets:
            continue
        try:
            curate_one(ds_id, cfg, args.overwrite)
            ok.append(ds_id)
        except Exception as exc:
            log.error("[%s] CURATE FAILED: %s", ds_id, exc)
            failures.append(ds_id)

    log.info("curated ok: %s", ok)
    if failures:
        log.error("curate failures: %s", failures)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
