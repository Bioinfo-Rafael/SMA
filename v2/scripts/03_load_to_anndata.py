#!/usr/bin/env python3
"""03_load_to_anndata.py -- dispatch each logical dataset to its loader and
write one interim h5ad per dataset into data/interim_h5ad/.

No QC / filtering / normalization happens here. raw counts -> adata.X;
processed status (TPM / SoupX) is recorded via the manifest data_status.
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import anndata_utils as au  # noqa: E402
import io_10x  # noqa: E402
import io_dense  # noqa: E402
import io_rds_bridge  # noqa: E402
from reporting import get_logger  # noqa: E402

log = get_logger("03_load")

LOADERS = {
    "10x_h5_per_sample": io_10x.load_10x_h5_per_sample,
    "10x_mtx_per_sample": io_10x.load_10x_mtx_per_sample,
    "combined_umi_tsv_with_metadata": io_dense.load_combined_umi_tsv_with_metadata,
    "dense_or_text_matrix_bundle": io_dense.load_dense_or_text_matrix_bundle,
    "mtx_or_text_bundle": io_dense.load_mtx_or_text_bundle,
    "dense_gene_by_cell_matrix": io_dense.load_dense_gene_by_cell_matrix,
    "processed_count_matrix_with_metadata": io_dense.load_processed_count_matrix_with_metadata,
    "nested_tar_dropseq": io_dense.load_nested_tar_dropseq,
    "rds_bridge": io_rds_bridge.load_rds_bridge,
}


def load_one(ds: dict, paths: dict, overwrite: bool):
    out_path = Path(paths["interim"]) / ds["output_h5ad"]
    if out_path.exists() and not overwrite:
        log.info("[%s] interim exists, skip (%s)", ds["dataset_id"], out_path.name)
        return out_path

    loader = LOADERS.get(ds["loader"])
    if loader is None:
        raise RuntimeError(f"no loader registered for '{ds['loader']}'")

    log.info("[%s] loading via %s", ds["dataset_id"], ds["loader"])
    adata = loader(ds, paths, log)
    au.finalize_anndata(adata, ds)
    log.info("[%s] assembled AnnData: %d cells x %d genes (data_status=%s)",
             ds["dataset_id"], adata.n_obs, adata.n_vars, ds.get("data_status"))
    au.save_h5ad(adata, out_path, overwrite=overwrite)
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(ROOT / "config" / "dataset_manifest.yaml"))
    ap.add_argument("--datasets", nargs="*", help="dataset_id filter")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    paths = au.project_paths(ROOT)
    au.ensure_dirs(paths)
    with open(args.manifest) as fh:
        manifest = yaml.safe_load(fh)

    ok, failures = [], []
    for ds in manifest["datasets"]:
        if args.datasets and ds["dataset_id"] not in args.datasets:
            continue
        try:
            load_one(ds, paths, args.overwrite)
            ok.append(ds["dataset_id"])
        except Exception as exc:
            log.error("[%s] LOAD FAILED: %s", ds["dataset_id"], exc)
            log.debug(traceback.format_exc())
            failures.append(ds["dataset_id"])

    log.info("loaded ok: %s", ok)
    if failures:
        log.error("load failures: %s", failures)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
