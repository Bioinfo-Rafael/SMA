"""Dense / text / combined-table / processed / nested loaders for notebooks,
plus the reader for R-exported intermediates (GSE295514).

Orientation policy: text matrices are assumed gene x cell and transposed to
cell x gene; the observed shape and head rows/cols are stashed in
adata.uns['orientation_report::<file>'] for the human to confirm.
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
import manifest_utils as mf
from archive_utils import find_files

log = logging.getLogger("io_dense")


def _open_text(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path, "rt")


def _sep_for(path) -> str:
    return "," if ".csv" in Path(path).name.lower() else "\t"


def _base_obs(barcodes, ds, *, sample_id, sample_label, gsm_id, source_file):
    obs = pd.DataFrame(index=pd.Index([str(b) for b in barcodes]))
    obs["cell_id_original"] = obs.index.astype(str)
    obs["sample_id"] = sample_id
    obs["sample_label"] = sample_label
    obs["gsm_id"] = gsm_id
    obs["source_file"] = source_file
    obs["source_accession"] = ds["source_accession"]
    obs["dataset_id"] = ds["dataset_id"]
    return obs


# --------------------------------------------------------------------------
# generic dense reader
# --------------------------------------------------------------------------
def read_dense_gene_by_cell_matrix(path):
    """Read a dense matrix (rows=genes, cols=cells) -> (X cells x genes csr,
    genes, cells, orientation_report)."""
    sep = _sep_for(path)
    df = pd.read_csv(path, sep=sep, index_col=0)
    report = {"file": Path(path).name, "raw_shape_rows_x_cols": list(df.shape),
              "assumed": "rows=genes, cols=cells",
              "head_index": [str(x) for x in df.index[:5]],
              "head_columns": [str(x) for x in df.columns[:5]]}
    genes = [str(x) for x in df.index]
    cells = [str(x) for x in df.columns]
    X = sp.csr_matrix(df.to_numpy(dtype=np.float32).T)
    return X, genes, cells, report


def _anndata_from_dense(path, ds, *, sample_label, gsm_id):
    X, genes, cells, report = read_dense_gene_by_cell_matrix(path)
    a = ad.AnnData(X=X)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.ensure_standard_var_columns(a, gene_symbol=genes)
    a.obs = _base_obs(a.obs_names, ds, sample_id=au.sanitize_id(sample_label),
                      sample_label=sample_label, gsm_id=gsm_id, source_file=Path(path).name)
    a.uns[f"orientation_report::{Path(path).name}"] = report
    log.info("  dense %-30s raw_shape=%s -> cells=%d genes=%d",
             Path(path).name, report["raw_shape_rows_x_cols"], a.n_obs, a.n_vars)
    return a


# --------------------------------------------------------------------------
# dense_or_text_matrix_bundle  (GSE167198)
# --------------------------------------------------------------------------
def _text_matrix_files(extracted: Path):
    return [p for p in find_files(extracted, ("*",))
            if p.name.lower().endswith((".txt.gz", ".txt", ".tsv.gz", ".tsv",
                                        ".csv.gz", ".csv", ".dge.txt.gz", ".dge.txt"))]


def load_dense_or_text_matrix_bundle(ds, paths, logger=log):
    extracted = mf.dataset_extracted_dir(paths, ds)
    files = _text_matrix_files(extracted)
    logger.info("[%s] %d text-matrix files", ds["dataset_id"], len(files))
    if not files:
        raise RuntimeError(f"no text matrices for {ds['dataset_id']} in {extracted}")
    parts = []
    for path in files:
        meta = au.parse_geo_filename(path.name)
        parts.append(_anndata_from_dense(path, ds, sample_label=meta["sample_label"],
                                         gsm_id=meta["gsm_id"]))
    return io_10x._concat(parts, ds["source_accession"])


# --------------------------------------------------------------------------
# mtx_or_text_bundle  (GSE167327)
# --------------------------------------------------------------------------
def report_side_table(path, ds, paths, logger=log, n_lines: int = 50):
    out = Path(paths["reports"]) / f"{ds['dataset_id']}_side_table_{Path(path).name}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"side table: {path}", ""]
    try:
        with _open_text(path) as fh:
            for i, line in enumerate(fh):
                if i >= n_lines:
                    lines.append(f"... (truncated at {n_lines} lines)")
                    break
                lines.append(line.rstrip("\n"))
    except Exception as exc:  # pragma: no cover
        lines.append(f"ERROR: {exc}")
    out.write_text("\n".join(lines))
    logger.info("[%s] side-table report -> %s", ds["dataset_id"], out)
    return out


def load_mtx_or_text_bundle(ds, paths, logger=log):
    extracted = mf.dataset_extracted_dir(paths, ds)
    side = ds.get("side_table")
    if side:
        side_path = mf.dataset_raw_dir(paths, ds) / side
        if side_path.exists():
            report_side_table(side_path, ds, paths, logger)

    groups = io_10x.group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    if complete:
        logger.info("[%s] %d MTX triplets", ds["dataset_id"], len(complete))
        return io_10x.load_10x_mtx_per_sample(ds, paths, logger)
    logger.info("[%s] no MTX triplets; falling back to dense text", ds["dataset_id"])
    return load_dense_or_text_matrix_bundle(ds, paths, logger)


# --------------------------------------------------------------------------
# dense_gene_by_cell_matrix  (GSE167331, TPM)
# --------------------------------------------------------------------------
def load_dense_gene_by_cell_matrix(ds, paths, logger=log):
    raw_dir = mf.dataset_raw_dir(paths, ds)
    fname = ds.get("matrix_file") or ds["files"][0]["name"]
    path = raw_dir / fname
    if not path.exists():
        raise RuntimeError(f"matrix file missing for {ds['dataset_id']}: {path}")
    a = _anndata_from_dense(path, ds, sample_label=ds["source_accession"], gsm_id=au.UNKNOWN)
    logger.info("[%s] matrix loaded: %d cells x %d genes (data_status hint=%s)",
                ds["dataset_id"], a.n_obs, a.n_vars, ds.get("data_status"))
    return a


# --------------------------------------------------------------------------
# combined_umi_tsv_with_metadata  (GSE173524)
# --------------------------------------------------------------------------
def _read_table(path, index_col=0):
    return pd.read_csv(path, sep=_sep_for(path), index_col=index_col)


def read_combined_umi_with_metadata(ds, paths, logger=log):
    raw_dir = mf.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or "GSE173524_umi.tsv.gz"
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"UMI table missing for {ds['dataset_id']}: {mpath}")

    logger.info("[%s] reading %s (genes x cells -> transpose)", ds["dataset_id"], matrix_file)
    df = _read_table(mpath)                                   # genes x cells
    genes = [str(x) for x in df.index]
    cells = [str(x) for x in df.columns]
    a = ad.AnnData(X=sp.csr_matrix(df.to_numpy(dtype=np.float32).T))
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.ensure_standard_var_columns(a, gene_symbol=genes)
    a.obs = _base_obs(a.obs_names, ds, sample_id=au.UNKNOWN, sample_label=au.UNKNOWN,
                      gsm_id=au.UNKNOWN, source_file=matrix_file)

    for entry in ds.get("metadata_files", []) or []:
        mp = raw_dir / entry["name"]
        if not mp.exists():
            logger.warning("  metadata file missing: %s", mp)
            continue
        meta = _read_table(mp).add_prefix("meta_")
        if entry.get("role") == "per_cell":
            meta.index = meta.index.astype(str)
            a.obs = a.obs.join(meta, how="left").loc[a.obs_names]
            logger.info("  joined per-cell metadata %s (%d cols)", entry["name"], meta.shape[1])
        else:
            a.uns[f"sample_metadata::{entry['name']}"] = meta.reset_index().astype(str).to_dict("list")
            logger.info("  stashed per-sample metadata %s in uns", entry["name"])
    return a


# --------------------------------------------------------------------------
# processed_count_matrix_with_metadata  (GSE206330, SoupX)
# --------------------------------------------------------------------------
def load_processed_count_matrix_with_metadata(ds, paths, logger=log):
    raw_dir = mf.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or ds["files"][0]["name"]
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"processed matrix missing for {ds['dataset_id']}: {mpath}")

    logger.info("[%s] reading processed matrix %s", ds["dataset_id"], matrix_file)
    X, genes, cells, report = read_dense_gene_by_cell_matrix(mpath)
    a = ad.AnnData(X=X)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(genes)
    a.var_names_make_unique()
    au.ensure_standard_var_columns(a, gene_symbol=genes)
    a.obs = _base_obs(a.obs_names, ds, sample_id=au.UNKNOWN, sample_label=au.UNKNOWN,
                      gsm_id=au.UNKNOWN, source_file=matrix_file)
    a.uns[f"orientation_report::{matrix_file}"] = report

    for entry in ds.get("metadata_files", []) or []:
        mp = raw_dir / entry["name"]
        if not mp.exists():
            logger.warning("  metadata file missing: %s", mp)
            continue
        meta = _read_table(mp).add_prefix("meta_")
        meta.index = meta.index.astype(str)
        a.obs = a.obs.join(meta, how="left").loc[a.obs_names]
        logger.info("  joined metadata %s (%d cols)", entry["name"], meta.shape[1])
    return a


# --------------------------------------------------------------------------
# nested_tar_dropseq  (GSE178693)
# --------------------------------------------------------------------------
def load_nested_tar_dropseq(ds, paths, logger=log):
    extracted = mf.dataset_extracted_dir(paths, ds)
    groups = io_10x.group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    text_files = _text_matrix_files(extracted)

    report = Path(paths["reports"]) / f"{ds['dataset_id']}_formats.txt"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(
        f"dataset: {ds['dataset_id']}\nextracted: {extracted}\n"
        f"mtx_triplets: {len(complete)}\ntext_matrices: {len(text_files)}\nfiles:\n"
        + "\n".join(f"  {p.relative_to(extracted)}" for p in find_files(extracted)))
    logger.info("[%s] nested formats: mtx=%d text=%d (report -> %s)",
                ds["dataset_id"], len(complete), len(text_files), report)

    parts = []
    if complete:
        for _, slot in sorted(complete.items()):
            a = io_10x.read_10x_mtx_triplet(slot["mtx"], slot["barcodes"], slot["features"])
            a.obs = _base_obs(a.obs_names, ds, sample_id=au.sanitize_id(slot["label"]),
                              sample_label=slot["label"], gsm_id=slot["gsm"],
                              source_file=Path(slot["mtx"]).name)
            parts.append(a)
    else:
        for path in text_files:
            meta = au.parse_geo_filename(path.name)
            try:
                parts.append(_anndata_from_dense(path, ds, sample_label=meta["sample_label"],
                                                 gsm_id=meta["gsm_id"]))
            except Exception as exc:
                logger.warning("  could not read %s: %s", path.name, exc)
    if not parts:
        raise RuntimeError(f"no readable matrices for {ds['dataset_id']}; see {report}")
    return io_10x._concat(parts, ds["source_accession"])


# --------------------------------------------------------------------------
# R intermediate reader  (GSE295514, after the R notebook)
# --------------------------------------------------------------------------
def read_from_r_intermediate(directory, ds, logger=log):
    """Assemble AnnData from files exported by notebooks/R/01_GSE295514_read_rds.ipynb:
    counts.mtx (genes x cells), barcodes.csv, genes.csv, metadata.csv."""
    directory = Path(directory)
    mtx = next((directory / n for n in ("counts.mtx.gz", "counts.mtx") if (directory / n).exists()), None)
    if mtx is None:
        raise RuntimeError(
            f"R intermediate counts.mtx not found in {directory}; "
            "run notebooks/R/01_GSE295514_read_rds.ipynb first")
    matrix = scipy.io.mmread(str(mtx)).tocsr()                # genes x cells
    genes = pd.read_csv(directory / "genes.csv")
    barcodes = pd.read_csv(directory / "barcodes.csv")
    gene_names = genes.iloc[:, 0].astype(str).to_numpy()
    cell_names = barcodes.iloc[:, 0].astype(str).to_numpy()

    a = ad.AnnData(X=matrix.T.tocsr().astype(np.float32))     # cells x genes
    a.obs_names = pd.Index(cell_names)
    a.var_names = pd.Index(gene_names)
    a.var_names_make_unique()
    au.ensure_standard_var_columns(a, gene_symbol=gene_names)
    a.obs = _base_obs(a.obs_names, ds, sample_id=au.UNKNOWN, sample_label=au.UNKNOWN,
                      gsm_id=au.UNKNOWN, source_file=ds.get("rds_file", "rds"))
    meta_path = directory / "metadata.csv"
    if meta_path.exists():
        meta = pd.read_csv(meta_path, index_col=0).add_prefix("meta_")
        meta.index = meta.index.astype(str)
        a.obs = a.obs.join(meta, how="left").loc[a.obs_names]
        logger.info("[%s] joined R meta.data (%d cols)", ds["dataset_id"], meta.shape[1])
    return a
