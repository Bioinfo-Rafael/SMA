#!/usr/bin/env python3
"""02_extract_archives.py -- safely unpack manifest archives (archive: true)
from data/raw/<acc>/ into data/extracted/<acc>/.

* path-traversal + symlink hardened
* nested tars (GSE178693) extracted recursively
* non-tar .gz supplementary files are left in place (notebooks read them directly)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import archive_utils as arc  # noqa: E402
import manifest_utils as mf  # noqa: E402

log = mf.get_logger("02_extract")

NESTED_LOADERS = {"nested_tar_dropseq"}


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
        raw_dir = mf.dataset_raw_dir(paths, ds)
        ext_dir = mf.dataset_extracted_dir(paths, ds)
        archives = [f for f in ds["files"] if f.get("archive")]
        if not archives:
            log.info("[%s] no archives to extract", ds["dataset_id"])
            continue
        recursive = ds.get("loader_hint") in NESTED_LOADERS

        for entry in archives:
            src = raw_dir / entry["name"]
            if not src.exists():
                log.error("[%s] archive not downloaded: %s", ds["dataset_id"], src)
                failures.append(ds["dataset_id"])
                continue
            marker = ext_dir / (entry["name"] + ".extracted")
            if marker.exists() and not args.force:
                log.info("[%s] already extracted %s (skip)", ds["dataset_id"], entry["name"])
                continue
            try:
                if recursive:
                    arc.extract_tar_recursive(src, ext_dir)
                else:
                    arc.extract_tar_safe(src, ext_dir)
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.write_text("ok\n")
            except Exception as exc:
                log.error("[%s] extraction FAILED for %s: %s",
                          ds["dataset_id"], entry["name"], exc)
                failures.append(ds["dataset_id"])

    if failures:
        log.error("extraction failures: %s", sorted(set(failures)))
        return 1
    log.info("all extractions complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
