"""10x loaders: per-sample CellRanger .h5 and per-sample MTX triplets.

Each sample becomes its own AnnData (full obs schema) and the per-sample
objects are concatenated with anndata.concat. raw counts -> adata.X.
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

log = logging.getLogger("io_10x")


def _open_text(path: Path):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _concat(parts: list, ds: dict):
    if not parts:
        raise RuntimeError(f"no samples loaded for {ds['dataset_id']}")
    if len(parts) == 1:
        return parts[0]
    return ad.concat(parts, join="outer", merge="first", uns_merge="first",
                     index_unique=None)


# --------------------------------------------------------------------------
# 10x H5
# --------------------------------------------------------------------------
def read_10x_h5_file(path: Path):
    import scanpy as sc
    adata = sc.read_10x_h5(str(path))
    adata.var_names_make_unique()
    return adata


def load_10x_h5_per_sample(ds: dict, paths: dict, logger=log):
    extracted = au.dataset_extracted_dir(paths, ds)
    h5_files = sorted(p for p in extracted.rglob("*.h5") if p.is_file())
    logger.info("[%s] found %d h5 files under %s", ds["dataset_id"], len(h5_files), extracted)
    if not h5_files:
        raise RuntimeError(f"no .h5 files for {ds['dataset_id']} in {extracted}")

    parts = []
    for h5 in h5_files:
        meta = au.parse_geo_filename(h5.name)
        sample_label = meta["sample_label"]
        sample_id = au.sanitize_id(sample_label)
        a = read_10x_h5_file(h5)

        ensembl = a.var["gene_ids"].to_numpy() if "gene_ids" in a.var else None
        ftype = a.var["feature_types"].to_numpy() if "feature_types" in a.var else None
        au.standardize_var(a, gene_symbol=a.var_names.to_numpy(),
                           gene_id=ensembl, ensembl_id=ensembl, feature_type=ftype)
        a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=sample_id,
                                   sample_label=sample_label, gsm_id=meta["gsm_id"],
                                   source_file=h5.name)
        a.X = sp.csr_matrix(a.X)
        parts.append(a)
        logger.info("  loaded sample %-20s cells=%d", sample_id, a.n_obs)
    return _concat(parts, ds)


# --------------------------------------------------------------------------
# 10x MTX triplets
# --------------------------------------------------------------------------
def _read_features(path: Path):
    rows = []
    with _open_text(path) as fh:
        for line in fh:
            rows.append(line.rstrip("\n").split("\t"))
    df = pd.DataFrame(rows)
    ncol = df.shape[1]
    gene_id = df[0].astype(str).to_numpy() if ncol >= 1 else None
    symbol = df[1].astype(str).to_numpy() if ncol >= 2 else gene_id
    ftype = df[2].astype(str).to_numpy() if ncol >= 3 else None
    return gene_id, symbol, ftype


def read_mtx_triplet(mtx_path: Path, barcodes_path: Path, features_path: Path):
    """Read a CellRanger MTX triplet. MTX is genes x cells -> transpose to
    cells x genes. Returns an AnnData with the required var columns set."""
    matrix = scipy.io.mmread(str(mtx_path)).tocsr()  # genes x cells
    with _open_text(barcodes_path) as fh:
        barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    gene_id, symbol, ftype = _read_features(features_path)

    X = matrix.T.tocsr().astype(np.float32)  # cells x genes
    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index([str(b) for b in barcodes])
    var_names = symbol if symbol is not None else gene_id
    adata.var_names = pd.Index([str(s) for s in var_names])
    adata.var_names_make_unique()
    au.standardize_var(adata, gene_symbol=var_names, gene_id=gene_id,
                       ensembl_id=gene_id, feature_type=ftype)
    return adata


def _group_mtx_triplets(extracted: Path) -> dict:
    """Group barcodes/features/matrix files of each sample by their parsed
    prefix. Returns {prefix: {'mtx':p, 'barcodes':p, 'features':p, ...}}."""
    groups: dict = {}
    for path in sorted(extracted.rglob("*")):
        if not path.is_file():
            continue
        meta = au.parse_geo_filename(path.name)
        if meta["kind"] not in ("mtx", "barcodes", "features"):
            continue
        key = str(path.parent / meta["prefix"])
        slot = groups.setdefault(key, {"prefix": meta["prefix"], "gsm": meta["gsm_id"],
                                       "label": meta["sample_label"], "dir": path.parent})
        slot[meta["kind"]] = path
    return groups


def load_10x_mtx_per_sample(ds: dict, paths: dict, logger=log):
    extracted = au.dataset_extracted_dir(paths, ds)
    groups = _group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    logger.info("[%s] found %d MTX sample groups (%d complete) under %s",
                ds["dataset_id"], len(groups), len(complete), extracted)
    if not complete:
        raise RuntimeError(
            f"no complete MTX triplets for {ds['dataset_id']} in {extracted}; "
            f"groups seen: {[g['prefix'] for g in groups.values()]}")

    parts = []
    for key, slot in sorted(complete.items()):
        sample_label = slot["label"]
        sample_id = au.sanitize_id(sample_label)
        a = read_mtx_triplet(slot["mtx"], slot["barcodes"], slot["features"])
        a.obs = au.make_sample_obs(a.obs_names, ds, sample_id=sample_id,
                                   sample_label=sample_label, gsm_id=slot["gsm"],
                                   source_file=Path(slot["mtx"]).name)
        parts.append(a)
        logger.info("  loaded sample %-24s cells=%d genes=%d",
                    sample_id, a.n_obs, a.n_vars)
    return _concat(parts, ds)
