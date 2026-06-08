"""Shared AnnData helpers: standard obs/var schema, sample classification,
global obs_names, data-status tagging and safe h5ad saving.

Design rules (from the project spec):
* Never auto-decide cell-type / cluster columns.
* Never drop original obs columns.
* raw counts go into adata.X; TPM / SoupX / processed status is recorded in
  adata.uns['data_status'] and adata.obs['processing_status'].
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("anndata_utils")

UNKNOWN = "unknown"

# Minimum obs columns every interim AnnData must carry.
REQUIRED_OBS_COLS = [
    "cell_id_original", "cell_id", "sample_id", "sample_label", "gsm_id",
    "parent_gse", "source_accession", "dataset_id", "tissue", "region",
    "assay", "technology", "enrichment", "disease_status", "disease_model",
    "genotype", "treatment", "age", "age_month", "sex", "replicate",
    "processing_status", "data_status", "source_file",
]

# Minimum var columns.
REQUIRED_VAR_COLS = [
    "gene_id", "gene_symbol", "gene_symbol_upper", "ensembl_id", "feature_type",
]


# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
def project_root() -> Path:
    """v2/ project root (this file lives in v2/src/)."""
    return Path(__file__).resolve().parent.parent


def project_paths(root: Path | None = None) -> dict:
    root = Path(root) if root is not None else project_root()
    data = root / "data"
    return {
        "root": root,
        "config": root / "config",
        "data": data,
        "raw": data / "raw",
        "extracted": data / "extracted",
        "interim": data / "interim_h5ad",
        "curated": data / "curated_h5ad",
        "reports": data / "reports",
    }


def ensure_dirs(paths: dict) -> None:
    for key, value in paths.items():
        if key in ("root", "config"):
            continue
        Path(value).mkdir(parents=True, exist_ok=True)


def dataset_raw_dir(paths: dict, ds: dict) -> Path:
    return Path(paths["raw"]) / ds["source_accession"]


def dataset_extracted_dir(paths: dict, ds: dict) -> Path:
    return Path(paths["extracted"]) / ds["source_accession"]


# --------------------------------------------------------------------------
# file-name parsing / sample classification
# --------------------------------------------------------------------------
_MATRIX_SUFFIX_RE = re.compile(
    r"[_.]?(filtered_feature_bc_matrix\.h5|raw_feature_bc_matrix\.h5|"
    r"feature_bc_matrix\.h5|matrix\.mtx(\.gz)?|barcodes\.tsv(\.gz)?|"
    r"features\.tsv(\.gz)?|genes\.tsv(\.gz)?)$",
    re.IGNORECASE,
)

_KIND_SUFFIXES = [
    ("filtered_feature_bc_matrix.h5", "h5"),
    ("raw_feature_bc_matrix.h5", "h5"),
    ("feature_bc_matrix.h5", "h5"),
    (".h5", "h5"),
    ("matrix.mtx.gz", "mtx"),
    ("matrix.mtx", "mtx"),
    ("barcodes.tsv.gz", "barcodes"),
    ("barcodes.tsv", "barcodes"),
    ("features.tsv.gz", "features"),
    ("features.tsv", "features"),
    ("genes.tsv.gz", "features"),
    ("genes.tsv", "features"),
]


def parse_geo_filename(fname: str) -> dict:
    """Pull GSM id, matrix kind and a sample 'prefix' (used to group the
    barcodes/features/matrix triplet of one 10x sample)."""
    base = Path(fname).name
    low = base.lower()
    gsm = None
    m = re.match(r"(GSM\d+)", base)
    if m:
        gsm = m.group(1)
    kind = "other"
    for suf, k in _KIND_SUFFIXES:
        if low.endswith(suf):
            kind = k
            break
    prefix = _MATRIX_SUFFIX_RE.sub("", base)
    # sample label = prefix without the leading GSM id token
    sample_label = re.sub(r"^GSM\d+[_.\-]?", "", prefix)
    return {
        "file": base,
        "gsm_id": gsm or UNKNOWN,
        "kind": kind,
        "prefix": prefix,
        "sample_label": sample_label or prefix,
    }


def sanitize_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-") or UNKNOWN


def classify_sample(name: str, rules) -> tuple[dict, list]:
    """Apply every rule whose regex matches `name`; merge their `set` dicts in
    list order (later rule wins on a shared key). Returns (values, matched)."""
    out: dict = {}
    matched: list = []
    for rule in rules or []:
        try:
            if re.search(rule["match"], name):
                out.update(rule.get("set", {}) or {})
                matched.append(rule["match"])
        except re.error as exc:  # pragma: no cover - bad manifest regex
            log.warning("bad regex %r in manifest: %s", rule.get("match"), exc)
    return out, matched


# --------------------------------------------------------------------------
# obs / var construction
# --------------------------------------------------------------------------
def make_sample_obs(index, ds: dict, *, sample_id: str, sample_label: str,
                    gsm_id: str, source_file: str, extra: dict | None = None) -> pd.DataFrame:
    """Build a full-schema obs frame for one sample (constants + per-sample
    rules). `index` are the original cell barcodes / ids."""
    idx = pd.Index([str(x) for x in index])
    obs = pd.DataFrame(index=idx)
    for col in REQUIRED_OBS_COLS:
        obs[col] = UNKNOWN

    constants = dict(ds.get("constants", {}) or {})
    constants.setdefault("parent_gse", ds.get("parent_gse", UNKNOWN))
    constants.setdefault("source_accession", ds.get("source_accession", UNKNOWN))
    constants.setdefault("dataset_id", ds.get("dataset_id", UNKNOWN))
    for key, val in constants.items():
        obs[key] = val

    obs["cell_id_original"] = idx.astype(str)
    obs["sample_id"] = sample_id
    obs["sample_label"] = sample_label
    obs["gsm_id"] = gsm_id
    obs["source_file"] = source_file
    obs["data_status"] = ds.get("data_status", UNKNOWN)
    obs["processing_status"] = ds.get("processing_status", UNKNOWN)

    # per-sample condition rules, matched against label + id + file name
    match_target = f"{sample_label} {sample_id} {source_file}"
    setvals, _ = classify_sample(match_target, ds.get("sample_rules"))
    for key, val in setvals.items():
        obs[key] = val

    if extra:
        for key, val in extra.items():
            obs[key] = val
    return obs


def standardize_var(adata, *, gene_symbol=None, gene_id=None,
                    ensembl_id=None, feature_type=None):
    """Ensure the required var columns exist. Arrays may be passed explicitly;
    missing pieces fall back to var_names / 'unknown'."""
    n = adata.n_vars
    if gene_symbol is None:
        gene_symbol = adata.var_names.astype(str).to_numpy()
    gene_symbol = np.asarray([str(x) for x in gene_symbol], dtype=object)

    adata.var["gene_symbol"] = gene_symbol
    adata.var["gene_symbol_upper"] = np.asarray(
        [s.upper() for s in gene_symbol], dtype=object)
    adata.var["gene_id"] = np.asarray(
        [str(x) for x in (gene_id if gene_id is not None else gene_symbol)],
        dtype=object)
    adata.var["ensembl_id"] = np.asarray(
        [str(x) for x in (ensembl_id if ensembl_id is not None else [UNKNOWN] * n)],
        dtype=object)
    adata.var["feature_type"] = np.asarray(
        [str(x) for x in (feature_type if feature_type is not None else [UNKNOWN] * n)],
        dtype=object)
    return adata


def init_obs_schema(adata):
    for col in REQUIRED_OBS_COLS:
        if col not in adata.obs.columns:
            adata.obs[col] = UNKNOWN
    return adata


def set_global_obs_names(adata, source_accession: str):
    """obs_names := {source_accession}_{sample_id}_{original_barcode}, made
    globally unique. Also writes obs['cell_id']."""
    n = adata.n_obs
    sid = (adata.obs["sample_id"].astype(str)
           if "sample_id" in adata.obs else pd.Series([UNKNOWN] * n))
    orig = (adata.obs["cell_id_original"].astype(str)
            if "cell_id_original" in adata.obs else pd.Series(adata.obs_names))
    names = [f"{source_accession}_{s}_{b}" for s, b in zip(sid, orig)]
    adata.obs_names = pd.Index(names)
    adata.obs_names_make_unique()
    adata.obs["cell_id"] = adata.obs_names.astype(str)
    return adata


def set_data_status(adata, ds: dict):
    data_status = ds.get("data_status", UNKNOWN)
    processing_status = ds.get("processing_status", UNKNOWN)
    adata.obs["data_status"] = data_status
    adata.obs["processing_status"] = processing_status
    adata.uns["data_status"] = data_status
    adata.uns["processing_status"] = processing_status
    adata.uns["dataset_id"] = ds.get("dataset_id", UNKNOWN)
    adata.uns["parent_gse"] = ds.get("parent_gse", UNKNOWN)
    adata.uns["source_accession"] = ds.get("source_accession", UNKNOWN)
    adata.uns["title"] = ds.get("title", UNKNOWN)
    adata.uns["loader"] = ds.get("loader", UNKNOWN)
    return adata


def to_csr(adata):
    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    return adata


def _stringify_object_cols(df: pd.DataFrame) -> None:
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)


def save_h5ad(adata, path, overwrite: bool = False):
    path = Path(path)
    if path.exists() and not overwrite:
        log.info("skip existing %s", path)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    _stringify_object_cols(adata.obs)
    _stringify_object_cols(adata.var)
    adata.write_h5ad(path)
    log.info("wrote %s  (%d cells x %d genes)", path, adata.n_obs, adata.n_vars)
    return path


def finalize_anndata(adata, ds: dict):
    """Apply the shared finishing steps after a loader returns an AnnData with
    X (cells x genes), var (symbols) and full-schema obs."""
    init_obs_schema(adata)
    if not set(REQUIRED_VAR_COLS).issubset(adata.var.columns):
        standardize_var(adata)
    set_data_status(adata, ds)
    set_global_obs_names(adata, ds.get("source_accession", UNKNOWN))
    to_csr(adata)
    return adata
