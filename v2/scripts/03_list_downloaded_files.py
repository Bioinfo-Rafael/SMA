#!/usr/bin/env python3
"""03_list_downloaded_files.py -- list what has been downloaded / extracted and
check it against the manifest. Writes data/reports/downloaded_files.txt.

This is the last .py step. AnnData loading and everything after happens in the
notebooks under notebooks/.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import manifest_utils as mf  # noqa: E402
from archive_utils import find_files  # noqa: E402

log = mf.get_logger("03_list")


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    args = ap.parse_args()

    paths = mf.project_paths(ROOT)
    mf.ensure_dirs(paths)
    manifest = mf.load_manifest(Path(args.manifest))

    lines: list[str] = []

    def emit(text: str = "") -> None:
        print(text)
        lines.append(text)

    all_present = True
    for ds in mf.list_datasets(manifest):
        raw_dir = mf.dataset_raw_dir(paths, ds)
        ext_dir = mf.dataset_extracted_dir(paths, ds)
        emit(f"\n=== {ds['dataset_id']}  ({ds['source_accession']}) ===")

        emit("  raw/ expected:")
        for f in ds.get("files", []):
            p = raw_dir / f["name"]
            ok = p.exists()
            optional = f.get("optional", False)
            if not ok and not optional:
                all_present = False
            size = _human(p.stat().st_size) if ok else "-"
            flag = "OK " if ok else ("opt" if optional else "MISSING")
            emit(f"    [{flag}] {f['name']:55s} {size}")

        ext_files = find_files(ext_dir)
        emit(f"  extracted/ : {len(ext_files)} files")
        for p in ext_files[:20]:
            emit(f"      {p.relative_to(ext_dir)}  ({_human(p.stat().st_size)})")
        if len(ext_files) > 20:
            emit(f"      ... (+{len(ext_files) - 20} more)")

        if ds.get("loader_hint") == "R_jupyter_kernel_manual":
            r_dir = Path(paths["intermediate_from_r"]) / ds.get("intermediate_from_r", "")
            r_files = find_files(r_dir)
            emit(f"  intermediate_from_r/ : {len(r_files)} files "
                 f"(produced by the R notebook)")

    out = Path(paths["reports"]) / "downloaded_files.txt"
    out.write_text("\n".join(lines))
    log.info("wrote file list -> %s", out)
    log.info("all required raw files present: %s", all_present)
    return 0 if all_present else 1


if __name__ == "__main__":
    raise SystemExit(main())
