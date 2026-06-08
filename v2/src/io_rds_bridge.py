"""RDS bridge loader (GSE295514).

Python does not read .rds directly. We shell out to scripts/rds_to_h5ad_bridge.R
which inspects the object class and writes a CellRanger-style MTX triplet plus a
meta.csv (Seurat counts + meta.data, or SingleCellExperiment counts + colData).
This module then assembles the AnnData from those intermediate files.

If R is unavailable or the object is neither Seurat nor SingleCellExperiment,
the R script writes a class report and exits non-zero; we surface that report
and raise a clear error.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

import pandas as pd

import anndata_utils as au
import io_10x

log = logging.getLogger("io_rds_bridge")


def _find_rscript() -> str | None:
    return shutil.which("Rscript")


def _bridge_script(paths: dict) -> Path:
    return Path(paths["root"]) / "scripts" / "rds_to_h5ad_bridge.R"


def run_rds_bridge(ds: dict, paths: dict, logger=log) -> Path:
    raw_dir = au.dataset_raw_dir(paths, ds)
    rds_name = ds.get("rds_file") or ds["files"][0]["name"]
    rds_path = raw_dir / rds_name
    if not rds_path.exists():
        raise RuntimeError(f"RDS file missing for {ds['dataset_id']}: {rds_path}")

    out_dir = au.dataset_extracted_dir(paths, ds) / "rds_bridge"
    out_dir.mkdir(parents=True, exist_ok=True)

    rscript = _find_rscript()
    bridge = _bridge_script(paths)
    if rscript is None:
        raise RuntimeError(
            "Rscript not found on PATH; cannot convert RDS for "
            f"{ds['dataset_id']}. Install R (+ Matrix, optionally Seurat / "
            "SingleCellExperiment) and re-run.")
    if not bridge.exists():
        raise RuntimeError(f"bridge R script missing: {bridge}")

    cmd = [rscript, str(bridge), str(rds_path), str(out_dir)]
    logger.info("[%s] running R bridge: %s", ds["dataset_id"], " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # always persist a class/diagnostic report
    report = Path(paths["reports"]) / f"{ds['dataset_id']}_rds_bridge.log"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        f"cmd: {' '.join(cmd)}\nreturncode: {proc.returncode}\n\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n")
    if proc.returncode != 0:
        raise RuntimeError(
            f"R bridge failed for {ds['dataset_id']} (rc={proc.returncode}); "
            f"see {report}")
    return out_dir


def load_rds_bridge(ds: dict, paths: dict, logger=log):
    out_dir = run_rds_bridge(ds, paths, logger)

    mtx = next((p for p in [out_dir / "matrix.mtx.gz", out_dir / "matrix.mtx"]
                if p.exists()), None)
    barcodes = next((p for p in [out_dir / "barcodes.tsv.gz", out_dir / "barcodes.tsv"]
                     if p.exists()), None)
    features = next((p for p in [out_dir / "features.tsv.gz", out_dir / "features.tsv"]
                     if p.exists()), None)
    if not (mtx and barcodes and features):
        raise RuntimeError(
            f"R bridge did not produce an MTX triplet in {out_dir} "
            f"for {ds['dataset_id']}")

    a = io_10x.read_mtx_triplet(mtx, barcodes, features)
    a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=au.UNKNOWN,
                               sample_label=au.UNKNOWN, gsm_id=au.UNKNOWN,
                               source_file=ds.get("rds_file", "rds"))

    meta_path = next((p for p in [out_dir / "meta.csv", out_dir / "meta.csv.gz"]
                      if p.exists()), None)
    if meta_path is not None:
        meta = pd.read_csv(meta_path, index_col=0).add_prefix("meta_")
        meta.index = meta.index.astype(str)
        joined = a.obs.join(meta, how="left")
        joined = joined.loc[a.obs_names]
        a.obs = joined
        logger.info("[%s] joined RDS meta.data (%d cols)",
                    ds["dataset_id"], meta.shape[1])
        # apply per-sample rules against any obvious condition columns
        _apply_rules_from_meta(a, ds, logger)
    return a


def _apply_rules_from_meta(a, ds: dict, logger):
    """Best-effort: if a metadata column clearly holds Control/rNLS8 labels,
    apply sample_rules to it. Otherwise leave as unknown for human curation."""
    rules = ds.get("sample_rules") or []
    if not rules:
        return
    import re
    candidate_cols = [c for c in a.obs.columns
                      if a.obs[c].astype(str).str.contains(
                          r"(?i)control|ctrl|rnls8|nls", regex=True, na=False).any()]
    for col in candidate_cols:
        for idx, value in a.obs[col].astype(str).items():
            setvals, matched = au.classify_sample(value, rules)
            if matched:
                for k, v in setvals.items():
                    a.obs.at[idx, k] = v
        logger.info("[%s] applied condition rules from meta column '%s'",
                    ds["dataset_id"], col)
        break
