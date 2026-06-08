"""Dense / text-matrix loaders and metadata-joined loaders.

Covered loaders:
  dense_or_text_matrix_bundle          (GSE167198)
  mtx_or_text_bundle                   (GSE167327)
  dense_gene_by_cell_matrix            (GSE167331, TPM)
  combined_umi_tsv_with_metadata       (GSE173524)
  processed_count_matrix_with_metadata (GSE206330, SoupX)
  nested_tar_dropseq                   (GSE178693)

Orientation policy: text matrices are assumed gene x cell and transposed to
cell x gene, but the observed shape and the first rows/cols are always written
to a report stub so a human can confirm.
"""
from __future__ import annotations

import gzip
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

import anndata_utils as au
import io_10x

log = logging.getLogger("io_dense")


def _open_text(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path, "rt")


def _sep_for(path: Path) -> str:
    low = Path(path).name.lower()
    if ".csv" in low:
        return ","
    return "\t"


# --------------------------------------------------------------------------
# generic dense matrix reader
# --------------------------------------------------------------------------
def read_dense_matrix(path: Path, *, assume="gene_by_cell"):
    """Read a dense matrix (rows=genes, cols=cells by assumption) and return
    (X cells x genes csr, gene_names, cell_names, orientation_report)."""
    sep = _sep_for(path)
    df = pd.read_csv(path, sep=sep, index_col=0)
    report = {
        "file": Path(path).name,
        "raw_shape_rows_x_cols": list(df.shape),
        "assumed_orientation": assume,
        "head_index": [str(x) for x in df.index[:5]],
        "head_columns": [str(x) for x in df.columns[:5]],
    }
    # rows=genes, cols=cells -> transpose to cells x genes
    genes = [str(x) for x in df.index]
    cells = [str(x) for x in df.columns]
    X = sp.csr_matrix(df.to_numpy(dtype=np.float32).T)
    return X, genes, cells, report


def _anndata_from_dense(path: Path, ds: dict, *, sample_label, gsm_id,
                        data_status_override=None):
    X, genes, cells, report = read_dense_matrix(path)
    a = ad.AnnData(X=X)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.standardize_var(a, gene_symbol=genes)
    sample_id = au.sanitize_id(sample_label)
    a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=sample_id,
                               sample_label=sample_label, gsm_id=gsm_id,
                               source_file=Path(path).name)
    if data_status_override:
        a.obs["data_status"] = data_status_override
    a.uns[f"orientation_report::{Path(path).name}"] = report
    log.info("  dense %-30s shape(raw)=%s -> cells=%d genes=%d",
             Path(path).name, report["raw_shape_rows_x_cols"], a.n_obs, a.n_vars)
    return a


# --------------------------------------------------------------------------
# dense_or_text_matrix_bundle  (GSE167198)
# --------------------------------------------------------------------------
def _text_matrix_files(extracted: Path):
    out = []
    for path in sorted(extracted.rglob("*")):
        if not path.is_file():
            continue
        low = path.name.lower()
        if low.endswith((".txt.gz", ".txt", ".tsv.gz", ".tsv", ".csv.gz", ".csv",
                         ".dge.txt.gz", ".dge.txt")):
            out.append(path)
    return out


def load_dense_or_text_matrix_bundle(ds: dict, paths: dict, logger=log):
    extracted = au.dataset_extracted_dir(paths, ds)
    files = _text_matrix_files(extracted)
    logger.info("[%s] %d text-matrix files under %s",
                ds["dataset_id"], len(files), extracted)
    if not files:
        raise RuntimeError(f"no text matrices for {ds['dataset_id']} in {extracted}")
    default_status = ds.get("data_status") or "unknown_text_matrix"
    parts = []
    for path in files:
        meta = au.parse_geo_filename(path.name)
        parts.append(_anndata_from_dense(
            path, ds, sample_label=meta["sample_label"], gsm_id=meta["gsm_id"],
            data_status_override=default_status))
    return io_10x._concat(parts, ds)


# --------------------------------------------------------------------------
# mtx_or_text_bundle  (GSE167327)
# --------------------------------------------------------------------------
def load_mtx_or_text_bundle(ds: dict, paths: dict, logger=log):
    extracted = au.dataset_extracted_dir(paths, ds)
    raw_dir = au.dataset_raw_dir(paths, ds)

    # report the small side table if present (do not assume it is a matrix)
    side = ds.get("side_table")
    if side:
        side_path = raw_dir / side
        if side_path.exists():
            _report_side_table(side_path, ds, paths, logger)

    groups = io_10x._group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    if complete:
        logger.info("[%s] using %d MTX triplets", ds["dataset_id"], len(complete))
        return io_10x.load_10x_mtx_per_sample(ds, paths, logger)

    logger.info("[%s] no MTX triplets; falling back to dense text matrices",
                ds["dataset_id"])
    return load_dense_or_text_matrix_bundle(ds, paths, logger)


def _report_side_table(path: Path, ds: dict, paths: dict, logger):
    out = Path(paths["reports"]) / f"{ds['dataset_id']}_side_table_{path.name}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"side table: {path}", ""]
    try:
        with _open_text(path) as fh:
            for i, line in enumerate(fh):
                if i >= 50:
                    lines.append("... (truncated at 50 lines)")
                    break
                lines.append(line.rstrip("\n"))
    except Exception as exc:  # pragma: no cover
        lines.append(f"ERROR reading side table: {exc}")
    out.write_text("\n".join(lines))
    logger.info("[%s] wrote side-table report -> %s", ds["dataset_id"], out)


# --------------------------------------------------------------------------
# dense_gene_by_cell_matrix  (GSE167331, TPM)
# --------------------------------------------------------------------------
def load_dense_gene_by_cell_matrix(ds: dict, paths: dict, logger=log):
    raw_dir = au.dataset_raw_dir(paths, ds)
    fname = ds.get("matrix_file") or ds["files"][0]["name"]
    path = raw_dir / fname
    if not path.exists():
        raise RuntimeError(f"matrix file missing for {ds['dataset_id']}: {path}")
    a = _anndata_from_dense(path, ds, sample_label=ds["source_accession"],
                            gsm_id=au.UNKNOWN,
                            data_status_override=ds.get("data_status"))
    logger.info("[%s] TPM matrix loaded: %d cells x %d genes (data_status=%s)",
                ds["dataset_id"], a.n_obs, a.n_vars, ds.get("data_status"))
    return a


# --------------------------------------------------------------------------
# combined_umi_tsv_with_metadata  (GSE173524)
# --------------------------------------------------------------------------
def _read_table(path: Path, index_col=0):
    return pd.read_csv(path, sep=_sep_for(path), index_col=index_col)


def load_combined_umi_tsv_with_metadata(ds: dict, paths: dict, logger=log):
    raw_dir = au.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or "GSE173524_umi.tsv.gz"
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"UMI table missing for {ds['dataset_id']}: {mpath}")

    logger.info("[%s] reading combined UMI table %s (genes x cells -> transpose)",
                ds["dataset_id"], matrix_file)
    df = _read_table(mpath)                       # genes x cells
    genes = [str(x) for x in df.index]
    cells = [str(x) for x in df.columns]
    X = sp.csr_matrix(df.to_numpy(dtype=np.float32).T)   # cells x genes

    a = ad.AnnData(X=X)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.standardize_var(a, gene_symbol=genes)
    a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=au.UNKNOWN,
                               sample_label=au.UNKNOWN, gsm_id=au.UNKNOWN,
                               source_file=matrix_file)

    # join per-cell + per-sample metadata (kept verbatim, prefixed meta_)
    for entry in ds.get("metadata_files", []) or []:
        mp = raw_dir / entry["name"]
        if not mp.exists():
            logger.warning("  metadata file missing: %s", mp)
            continue
        meta = _read_table(mp)
        meta = meta.add_prefix("meta_")
        if entry.get("role") == "per_cell":
            joined = a.obs.join(meta, how="left")
            joined = joined.loc[a.obs_names]
            a.obs = joined
            logger.info("  joined per-cell metadata %s (%d cols)",
                        entry["name"], meta.shape[1])
        else:
            a.uns[f"sample_metadata::{entry['name']}"] = meta.reset_index().astype(str).to_dict("list")
            logger.info("  stashed per-sample metadata %s into uns", entry["name"])
    return a


# --------------------------------------------------------------------------
# processed_count_matrix_with_metadata  (GSE206330, SoupX)
# --------------------------------------------------------------------------
def load_processed_count_matrix_with_metadata(ds: dict, paths: dict, logger=log):
    raw_dir = au.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or ds["files"][0]["name"]
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"processed matrix missing for {ds['dataset_id']}: {mpath}")

    logger.info("[%s] reading processed matrix %s (data_status=%s)",
                ds["dataset_id"], matrix_file, ds.get("data_status"))
    X, genes, cells, report = read_dense_matrix(mpath)
    a = ad.AnnData(X=X)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.standardize_var(a, gene_symbol=genes)
    a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=au.UNKNOWN,
                               sample_label=au.UNKNOWN, gsm_id=au.UNKNOWN,
                               source_file=matrix_file)
    a.uns[f"orientation_report::{matrix_file}"] = report

    for entry in ds.get("metadata_files", []) or []:
        mp = raw_dir / entry["name"]
        if not mp.exists():
            logger.warning("  metadata file missing: %s", mp)
            continue
        meta = _read_table(mp).add_prefix("meta_")
        joined = a.obs.join(meta, how="left").loc[a.obs_names]
        a.obs = joined
        logger.info("  joined metadata %s (%d cols)", entry["name"], meta.shape[1])
    return a


# --------------------------------------------------------------------------
# nested_tar_dropseq  (GSE178693)
# --------------------------------------------------------------------------
def load_nested_tar_dropseq(ds: dict, paths: dict, logger=log):
    """RAW.tar may contain per-sample tars (already expanded by 02 via
    safe_extract_recursive). Detect MTX triplets first, else dense matrices."""
    extracted = au.dataset_extracted_dir(paths, ds)

    formats = {"mtx_triplets": 0, "text_matrices": 0}
    groups = io_10x._group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    formats["mtx_triplets"] = len(complete)
    text_files = _text_matrix_files(extracted)
    formats["text_matrices"] = len(text_files)

    report = Path(paths["reports"]) / f"{ds['dataset_id']}_formats.txt"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        f"dataset: {ds['dataset_id']}\nextracted: {extracted}\n"
        f"mtx_triplets: {len(complete)}\ntext_matrices: {len(text_files)}\n"
        + "files:\n" + "\n".join(f"  {p.relative_to(extracted)}"
                                 for p in sorted(extracted.rglob('*')) if p.is_file())
    )
    logger.info("[%s] nested formats: %s (report -> %s)",
                ds["dataset_id"], formats, report)

    parts = []
    if complete:
        for key, slot in sorted(complete.items()):
            sample_label = slot["label"]
            a = io_10x.read_mtx_triplet(slot["mtx"], slot["barcodes"], slot["features"])
            a.obs = au.make_sample_obs(a.obs_names, ds,
                                       sample_id=au.sanitize_id(sample_label),
                                       sample_label=sample_label, gsm_id=slot["gsm"],
                                       source_file=Path(slot["mtx"]).name)
            parts.append(a)
    else:
        for path in text_files:
            meta = au.parse_geo_filename(path.name)
            try:
                parts.append(_anndata_from_dense(
                    path, ds, sample_label=meta["sample_label"], gsm_id=meta["gsm_id"],
                    data_status_override=ds.get("data_status")))
            except Exception as exc:
                logger.warning("  could not parse %s as matrix: %s", path.name, exc)

    if not parts:
        raise RuntimeError(
            f"no readable matrices for {ds['dataset_id']}; see {report}")
    return io_10x._concat(parts, ds)
