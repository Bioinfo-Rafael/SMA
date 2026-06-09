#!/usr/bin/env python3
"""01_download_geo_supplement.py -- download every manifest file into
data/raw/<source_accession>/. Resumable; complete files are skipped.
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
    ap.add_argument("--datasets", nargs="*", help="dataset_id / source_accession filter")
    ap.add_argument("--force", action="store_true")
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
            log.error("FAILED %s: %s", ds["dataset_id"], exc)
            failures.append(ds["dataset_id"])

    if failures:
        log.error("download failures: %s", failures)
        return 1
    log.info("all downloads complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
