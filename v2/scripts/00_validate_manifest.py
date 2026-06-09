#!/usr/bin/env python3
"""00_validate_manifest.py -- config/dataset_manifest.yaml を検証し、
一覧表を data/reports/manifest_overview.csv に書き出す。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import manifest_utils as mf  # noqa: E402

log = mf.get_logger("00_validate")


def overview(manifest: dict) -> pd.DataFrame:
    """manifest を1行1データセットの表に整形。"""
    rows = []
    for ds in mf.list_datasets(manifest):
        rows.append({
            "dataset_id": ds.get("dataset_id"),
            "source_accession": ds.get("source_accession"),
            "parent_gse": ds.get("parent_gse"),
            "loader_hint": ds.get("loader_hint"),
            "data_status": ds.get("data_status"),
            "n_files": len(ds.get("files", [])),
            "files": ";".join(f["name"] for f in ds.get("files", [])),
            "output": ds.get("output"),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--strict", action="store_true", help="エラー時に非0終了")
    args = ap.parse_args()

    paths = mf.project_paths(ROOT)
    mf.ensure_dirs(paths)
    manifest = mf.load_manifest(Path(args.manifest))
    log.info("%d データセットを読み込み", len(mf.list_datasets(manifest)))

    errors = mf.validate_manifest(manifest)
    for e in errors:
        log.error(e)

    df = overview(manifest)
    out = Path(paths["reports"]) / "manifest_overview.csv"
    df.to_csv(out, index=False)
    log.info("一覧を書き出し -> %s", out)
    print(df.to_string(index=False))

    if errors:
        log.error("検証エラー %d 件", len(errors))
        return 1 if args.strict else 0
    log.info("manifest OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
