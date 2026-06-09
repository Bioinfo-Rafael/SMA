#!/usr/bin/env python3
"""04_overview.py -- プロジェクト全体の俯瞰。

manifest のデータセット一覧、各種パス、raw/extracted のファイル有無をまとめて表示する。
（旧 notebooks/python/00_overview.ipynb をスクリプト化したもの。AnnData は読まない。）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import manifest_utils as mf  # noqa: E402
from archive_utils import find_files  # noqa: E402

log = mf.get_logger("04_overview")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    args = ap.parse_args()

    paths = mf.project_paths(ROOT)
    mf.ensure_dirs(paths)
    manifest = mf.load_manifest(Path(args.manifest))

    # 1) データセット一覧
    print("=== データセット一覧 ===")
    rows = []
    for ds in mf.list_datasets(manifest):
        rows.append({
            "dataset_id": ds["dataset_id"],
            "source_accession": ds["source_accession"],
            "parent_gse": ds.get("parent_gse"),
            "loader_hint": ds["loader_hint"],
            "data_status": ds.get("data_status"),
            "output": ds["output"],
        })
    print(pd.DataFrame(rows).to_string(index=False))

    # 2) パス
    print("\n=== パス ===")
    for key, value in paths.items():
        print(f"  {key:20s} {value}")

    # 3) ファイル有無（raw / extracted）
    print("\n=== ファイル有無 ===")
    rows = []
    for ds in mf.list_datasets(manifest):
        raw_dir = mf.dataset_raw_dir(paths, ds)
        ext_dir = mf.dataset_extracted_dir(paths, ds)
        raw_present = sum((raw_dir / f["name"]).exists() for f in ds["files"])
        rows.append({
            "dataset_id": ds["dataset_id"],
            "raw_present": f"{raw_present}/{len(ds['files'])}",
            "extracted_files": len(find_files(ext_dir)),
        })
    print(pd.DataFrame(rows).to_string(index=False))

    print("\n次は notebooks/python/01_load_and_inspect_each_gse.ipynb へ。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
