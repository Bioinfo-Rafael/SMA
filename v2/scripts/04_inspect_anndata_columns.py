#!/usr/bin/env python3
"""04_inspect_anndata_columns.py -- write manual-inspection reports for every
interim h5ad.

This script deliberately does NOT pick cell-type / cluster columns. It dumps
every obs/var column, dtype, missing count, value counts and numeric summaries
so a human can decide and fill config/curation_template.yaml.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
from reporting import get_logger, inspect_anndata  # noqa: E402

log = get_logger("04_inspect")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--datasets", nargs="*", help="dataset_id filter")
    args = ap.parse_args()

    paths = au.project_paths(ROOT)
    au.ensure_dirs(paths)
    with open(args.manifest) as fh:
        manifest = yaml.safe_load(fh)

    done, missing = [], []
    for ds in manifest["datasets"]:
        if args.datasets and ds["dataset_id"] not in args.datasets:
            continue
        interim = Path(paths["interim"]) / ds["output_h5ad"]
        if not interim.exists():
            log.warning("[%s] interim not found: %s", ds["dataset_id"], interim)
            missing.append(ds["dataset_id"])
            continue
        log.info("[%s] inspecting %s", ds["dataset_id"], interim.name)
        adata = ad.read_h5ad(interim)
        out = inspect_anndata(adata, ds["dataset_id"], paths["reports"])
        log.info("  reports -> %s", out["summary"])
        done.append(ds["dataset_id"])

    log.info("inspected: %s", done)
    if missing:
        log.warning("missing interim h5ad for: %s", missing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
