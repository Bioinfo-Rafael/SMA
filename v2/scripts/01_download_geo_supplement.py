#!/usr/bin/env python3
"""01_download_geo_supplement.py -- download every file in the manifest into
data/raw/<source_accession>/.

Downloads resume / skip when already complete. Use --datasets to restrict.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
import geo_download as gd  # noqa: E402
from reporting import get_logger  # noqa: E402

log = get_logger("01_download")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--datasets", nargs="*", help="dataset_id filter")
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    args = ap.parse_args()

    paths = au.project_paths(ROOT)
    au.ensure_dirs(paths)
    with open(args.manifest) as fh:
        manifest = yaml.safe_load(fh)

    failures = []
    for ds in manifest["datasets"]:
        if args.datasets and ds["dataset_id"] not in args.datasets:
            continue
        dest = au.dataset_raw_dir(paths, ds)
        log.info("=== %s -> %s ===", ds["dataset_id"], dest)
        try:
            gd.download_files(ds["files"], dest, force=args.force)
        except Exception as exc:
            log.error("FAILED %s: %s", ds["dataset_id"], exc)
            failures.append(ds["dataset_id"])

    if failures:
        log.error("download failures: %s", failures)
        return 1
    log.info("all downloads complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
