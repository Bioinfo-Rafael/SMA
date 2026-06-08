#!/usr/bin/env python3
"""06_read_saved_h5ad_template.py -- TEMPLATE for later integration / comparison.

This is a scaffold only. It loads the curated h5ad files, compares their
obs/var schemas, looks at shared genes via gene_symbol_upper, and groups
datasets by data_status. It deliberately does NOT integrate, normalize or
batch-correct anything yet.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
from reporting import get_logger  # noqa: E402

log = get_logger("06_template")


def find_curated(curation: dict):
    out = []
    for ds_id, cfg in (curation.get("datasets") or {}).items():
        p = cfg["output_h5ad"]
        p = (ROOT / p) if not Path(p).is_absolute() else Path(p)
        if p.exists():
            out.append((ds_id, p))
        else:
            log.warning("[%s] curated h5ad not found yet: %s", ds_id, p)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curation", default=str(ROOT / "config" / "curation_template.yaml"))
    args = ap.parse_args()

    with open(args.curation) as fh:
        curation = yaml.safe_load(fh)

    curated = find_curated(curation)
    if not curated:
        log.warning("no curated h5ad found; run 05_curate_and_save_h5ad.py first")
        return 0

    by_status = defaultdict(list)
    gene_sets = {}
    obs_cols = {}
    for ds_id, path in curated:
        a = ad.read_h5ad(path, backed="r")
        status = str(a.uns.get("data_status", "unknown"))
        by_status[status].append(ds_id)
        gene_sets[ds_id] = set(a.var["gene_symbol_upper"].astype(str)) \
            if "gene_symbol_upper" in a.var.columns else set(a.var_names.astype(str).str.upper())
        obs_cols[ds_id] = list(a.obs.columns)
        log.info("[%s] %d cells x %d genes  data_status=%s",
                 ds_id, a.n_obs, a.n_vars, status)

    # ---- schema comparison ----
    log.info("=== obs schema comparison ===")
    all_cols = set().union(*obs_cols.values())
    common_cols = set(all_cols)
    for cols in obs_cols.values():
        common_cols &= set(cols)
    log.info("obs columns common to ALL datasets (%d): %s",
             len(common_cols), sorted(common_cols))
    for ds_id, cols in obs_cols.items():
        extra = sorted(set(cols) - common_cols)
        if extra:
            log.info("[%s] dataset-specific obs columns: %s", ds_id, extra)

    # ---- shared genes ----
    if gene_sets:
        shared = set.intersection(*gene_sets.values())
        log.info("=== genes shared (gene_symbol_upper) across ALL %d datasets: %d ===",
                 len(gene_sets), len(shared))

    # ---- data_status grouping + warning ----
    log.info("=== datasets grouped by data_status ===")
    for status, ids in by_status.items():
        log.info("  %-28s : %s", status, ids)
    raw = by_status.get("raw_counts", [])
    non_raw = [s for s in by_status if s != "raw_counts"]
    if non_raw:
        log.warning("NON-raw data_status present: %s. Do NOT pool these with "
                    "raw_counts (%s) without status-aware handling "
                    "(e.g. TPM and SoupX-corrected are not raw UMI counts).",
                    {s: by_status[s] for s in non_raw}, raw)

    log.info("TEMPLATE ONLY: no integration / normalization / batch correction performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
