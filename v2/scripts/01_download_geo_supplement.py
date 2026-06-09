#!/usr/bin/env python3
"""01_download_geo_supplement.py -- manifest の全ファイルを
data/raw/<source_accession>/ にダウンロードする（レジューム対応・完了済みはスキップ）。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import geo_download as gd  # noqa: E402
import manifest_utils as mf  # noqa: E402

log = mf.get_logger("01_download")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--datasets", nargs="*", help="dataset_id / source_accession で絞り込み")
    ap.add_argument("--force", action="store_true", help="既存でも取り直す")
    args = ap.parse_args()

    paths = mf.project_paths(ROOT)
    mf.ensure_dirs(paths)
    manifest = mf.load_manifest(Path(args.manifest))

    failures = []
    for ds in mf.list_datasets(manifest):
        if args.datasets and ds["dataset_id"] not in args.datasets \
                and ds["source_accession"] not in args.datasets:
            continue
        dest = mf.dataset_raw_dir(paths, ds)
        log.info("=== %s -> %s ===", ds["dataset_id"], dest)
        try:
            gd.download_files(ds["files"], dest, force=args.force)
        except Exception as exc:
            log.error("失敗 %s: %s", ds["dataset_id"], exc)
            failures.append(ds["dataset_id"])

    if failures:
        log.error("ダウンロード失敗: %s", failures)
        return 1
    log.info("全ダウンロード完了")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
