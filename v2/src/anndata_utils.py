"""ノートブックで共通利用する AnnData ヘルパー。

標準 obs/var スキーマの付与、obs_names の一意化、h5ad の保存/一括ロードなど。
merge は anndata.concat をそのまま使えばよいので、ここでは薄い補助のみ置く
（notebooks/python/04 は ad.concat を直接呼ぶ）。

生カウントは adata.X に入れる。TPM / SoupX / RDS などの処理状態は
obs['data_status'] / uns['data_status'] にノートブック側で明示する。
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("anndata_utils")

UNKNOWN = "unknown"

# 標準 obs スキーマ（notebook 03 で埋める。unknown のままでも可）
REQUIRED_OBS_COLS = [
    "cell_id_original", "cell_id", "sample_id", "sample_label",
    "source_accession", "parent_gse", "dataset_id", "species", "disease_area",
    "tissue", "region", "assay", "technology", "enrichment", "disease_status",
    "disease_model", "genotype", "treatment", "age", "age_month", "sex",
    "replicate", "data_status", "processing_status", "source_file",
]
REQUIRED_VAR_COLS = [
    "gene_id", "gene_symbol", "gene_symbol_upper", "ensembl_id", "feature_type",
]


# --------------------------------------------------------------------------
# ファイル名のパース（GEO ファイル名から機械的に sample id を取り出す）
# --------------------------------------------------------------------------
_MATRIX_SUFFIX_RE = re.compile(
    r"[_.]?(filtered_feature_bc_matrix\.h5|raw_feature_bc_matrix\.h5|"
    r"feature_bc_matrix\.h5|matrix\.mtx(\.gz)?|barcodes\.tsv(\.gz)?|"
    r"features\.tsv(\.gz)?|genes\.tsv(\.gz)?)$",
    re.IGNORECASE,
)
_KIND_SUFFIXES = [
    ("filtered_feature_bc_matrix.h5", "h5"), ("raw_feature_bc_matrix.h5", "h5"),
    ("feature_bc_matrix.h5", "h5"), (".h5", "h5"),
    ("matrix.mtx.gz", "mtx"), ("matrix.mtx", "mtx"),
    ("barcodes.tsv.gz", "barcodes"), ("barcodes.tsv", "barcodes"),
    ("features.tsv.gz", "features"), ("features.tsv", "features"),
    ("genes.tsv.gz", "features"), ("genes.tsv", "features"),
]


def sanitize_id(text: str) -> str:
    """id に使える文字だけ残す。"""
    return re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-") or UNKNOWN


def parse_geo_filename(fname: str) -> dict:
    """GEO のファイル名から GSM id・種別(kind)・サンプル prefix を取り出す。"""
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
    prefix = _MATRIX_SUFFIX_RE.sub("", base)               # 三点セットを束ねる用
    sample_label = re.sub(r"^GSM\d+[_.\-]?", "", prefix) or prefix
    return {"file": base, "gsm_id": gsm or UNKNOWN, "kind": kind,
            "prefix": prefix, "sample_label": sample_label}


# --------------------------------------------------------------------------
# obs / var スキーマ
# --------------------------------------------------------------------------
def ensure_standard_obs_columns(adata, defaults: dict | None = None):
    """欠けている標準 obs 列を追加（unknown もしくは defaults で埋める）。既存列は消さない。"""
    defaults = defaults or {}
    if "cell_id_original" not in adata.obs.columns:
        adata.obs["cell_id_original"] = adata.obs_names.astype(str)
    for col in REQUIRED_OBS_COLS:
        if col not in adata.obs.columns:
            adata.obs[col] = defaults.get(col, UNKNOWN)
    for col, val in defaults.items():
        adata.obs[col] = val
    if (adata.obs["cell_id"] == UNKNOWN).all():
        adata.obs["cell_id"] = adata.obs_names.astype(str)
    return adata


def ensure_standard_var_columns(adata, *, gene_symbol=None, gene_id=None,
                                ensembl_id=None, feature_type=None):
    """gene_id / gene_symbol / gene_symbol_upper / ensembl_id / feature_type を整える。
    与えられなければ var_names / 'unknown' で補完。"""
    n = adata.n_vars
    if gene_symbol is None:
        gene_symbol = (adata.var["gene_symbol"].to_numpy()
                       if "gene_symbol" in adata.var.columns
                       else adata.var_names.astype(str).to_numpy())
    gene_symbol = np.asarray([str(x) for x in gene_symbol], dtype=object)
    adata.var["gene_symbol"] = gene_symbol
    adata.var["gene_symbol_upper"] = np.asarray([s.upper() for s in gene_symbol], dtype=object)
    if gene_id is None:
        gene_id = (adata.var["gene_id"].to_numpy()
                   if "gene_id" in adata.var.columns else gene_symbol)
    adata.var["gene_id"] = np.asarray([str(x) for x in gene_id], dtype=object)
    if ensembl_id is None:
        ensembl_id = (adata.var["ensembl_id"].to_numpy()
                      if "ensembl_id" in adata.var.columns else [UNKNOWN] * n)
    adata.var["ensembl_id"] = np.asarray([str(x) for x in ensembl_id], dtype=object)
    if feature_type is None:
        feature_type = (adata.var["feature_type"].to_numpy()
                        if "feature_type" in adata.var.columns else [UNKNOWN] * n)
    adata.var["feature_type"] = np.asarray([str(x) for x in feature_type], dtype=object)
    return adata


def make_obs_names_unique(adata, prefix: str | None = None):
    """obs_names を全体で一意化する。prefix を渡すと
    {prefix}_{sample_id}_{元バーコード} 形式にしてから一意化する。"""
    if "cell_id_original" not in adata.obs.columns:
        adata.obs["cell_id_original"] = adata.obs_names.astype(str)
    if prefix is not None:
        sid = (adata.obs["sample_id"].astype(str)
               if "sample_id" in adata.obs.columns else pd.Series([UNKNOWN] * adata.n_obs))
        orig = adata.obs["cell_id_original"].astype(str)
        adata.obs_names = pd.Index([f"{prefix}_{s}_{b}" for s, b in zip(sid, orig)])
    adata.obs_names_make_unique()
    adata.obs["cell_id"] = adata.obs_names.astype(str)
    return adata


def to_csr(adata):
    """X を CSR sparse に揃える。"""
    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    return adata


# --------------------------------------------------------------------------
# 保存 / ロード
# --------------------------------------------------------------------------
def _stringify_object_cols(df: pd.DataFrame) -> None:
    # h5ad 保存時に object dtype が原因で落ちないよう str に寄せる
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)


def save_h5ad(adata, path, *, overwrite: bool = True, sparse: bool = True):
    """h5ad を保存（overwrite=False なら既存はスキップ）。"""
    path = Path(path)
    if path.exists() and not overwrite:
        log.info("既存のためスキップ %s", path.name)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    if sparse:
        to_csr(adata)
    _stringify_object_cols(adata.obs)
    _stringify_object_cols(adata.var)
    adata.write_h5ad(path)
    log.info("保存 %s (%d cells x %d genes)", path.name, adata.n_obs, adata.n_vars)
    return path


def load_h5ad_collection(directory, pattern: str = "*.h5ad", backed=None) -> dict:
    """ディレクトリ内の h5ad を {ファイル名stem: AnnData} で一括ロード。"""
    directory = Path(directory)
    out: dict = {}
    for p in sorted(directory.glob(pattern)):
        out[p.stem] = ad.read_h5ad(p, backed=backed)
    return out


# --------------------------------------------------------------------------
# 集計（merge 後の細胞数表など）
# --------------------------------------------------------------------------
def cell_count_table(adata, keys) -> pd.DataFrame:
    """obs の各キー値ごとの細胞数（long 形式）。"""
    rows = []
    for key in keys:
        if key in adata.obs.columns:
            for value, count in adata.obs[key].astype(str).value_counts().items():
                rows.append({"column": key, "value": value, "n_cells": int(count)})
    return pd.DataFrame(rows)
