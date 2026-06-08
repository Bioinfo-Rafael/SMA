#!/usr/bin/env python3
"""00_make_manifest.py -- validate the hand-authored dataset_manifest.yaml and
emit a flat overview table.

The YAML is the source of truth; this step only checks it is well-formed and
produces data/reports/manifest_overview.csv for humans.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
from reporting import get_logger  # noqa: E402

log = get_logger("00_make_manifest")

REQUIRED_KEYS = ["dataset_id", "parent_gse", "source_accession", "loader",
                 "output_h5ad", "files"]
KNOWN_LOADERS = {
    "10x_h5_per_sample", "10x_mtx_per_sample", "combined_umi_tsv_with_metadata",
    "dense_or_text_matrix_bundle", "mtx_or_text_bundle", "dense_gene_by_cell_matrix",
    "processed_count_matrix_with_metadata", "nested_tar_dropseq", "rds_bridge",
}


def load_manifest(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def validate(manifest: dict) -> list[str]:
    errors: list[str] = []
    datasets = manifest.get("datasets", [])
    seen_ids, seen_outputs = set(), set()
    for i, ds in enumerate(datasets):
        tag = ds.get("dataset_id", f"<index {i}>")
        for key in REQUIRED_KEYS:
            if not ds.get(key):
                errors.append(f"{tag}: missing required key '{key}'")
        loader = ds.get("loader")
        if loader and loader not in KNOWN_LOADERS:
            errors.append(f"{tag}: unknown loader '{loader}'")
        if ds.get("dataset_id") in seen_ids:
            errors.append(f"{tag}: duplicate dataset_id")
        seen_ids.add(ds.get("dataset_id"))
        if ds.get("output_h5ad") in seen_outputs:
            errors.append(f"{tag}: duplicate output_h5ad")
        seen_outputs.add(ds.get("output_h5ad"))
        for f in ds.get("files", []):
            if not f.get("name") or not f.get("url"):
                errors.append(f"{tag}: file entry missing name/url: {f}")
    return errors


def overview(manifest: dict) -> pd.DataFrame:
    rows = []
    for ds in manifest.get("datasets", []):
        rows.append({
            "dataset_id": ds.get("dataset_id"),
            "parent_gse": ds.get("parent_gse"),
            "source_accession": ds.get("source_accession"),
            "loader": ds.get("loader"),
            "data_status": ds.get("data_status"),
            "processing_status": ds.get("processing_status"),
            "n_files": len(ds.get("files", [])),
            "files": ";".join(f["name"] for f in ds.get("files", [])),
            "output_h5ad": ds.get("output_h5ad"),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--strict", action="store_true", help="exit non-zero on errors")
    args = ap.parse_args()

    paths = au.project_paths(ROOT)
    au.ensure_dirs(paths)

    manifest = load_manifest(Path(args.manifest))
    n = len(manifest.get("datasets", []))
    log.info("loaded %d datasets from %s", n, args.manifest)

    errors = validate(manifest)
    for e in errors:
        log.error(e)

    df = overview(manifest)
    out = Path(paths["reports"]) / "manifest_overview.csv"
    df.to_csv(out, index=False)
    log.info("wrote overview -> %s", out)
    print(df.to_string(index=False))

    if errors:
        log.error("%d validation error(s)", len(errors))
        if args.strict:
            return 1
    else:
        log.info("manifest OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
