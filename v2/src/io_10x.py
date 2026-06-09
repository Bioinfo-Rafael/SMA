"""ノートブック用の 10x ローダー：サンプルごとの CellRanger .h5 と MTX 三点セット。

生カウントを X に入れた AnnData を作る（var は gene_symbol/ensembl、obs はファイル名
由来の機械的なサンプル情報のみ）。genotype / treatment / condition などの生物学的な
obs はここでは付けず、notebook 03(curate) で人間が設定する。
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
import manifest_utils as mf
from archive_utils import find_files

log = logging.getLogger("io_10x")


def _open_text(path: Path):
    # .gz ならテキストモードで開く
    path = Path(path)
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path, "rt")


def _concat(parts, source_accession):
    # サンプル AnnData を結合し、obs_names を {acc}_{sample}_{barcode} で一意化
    if not parts:
        raise RuntimeError(f"サンプルを1つも読めませんでした: {source_accession}")
    adata = parts[0] if len(parts) == 1 else ad.concat(
        parts, join="outer", merge="first", uns_merge="first", index_unique=None)
    au.make_obs_names_unique(adata, prefix=source_accession)
    return adata


def _base_obs(barcodes, ds, *, sample_id, sample_label, gsm_id, source_file):
    # ファイル名から取れる機械的な obs だけを入れる
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
# 10x H5
# --------------------------------------------------------------------------
def read_10x_h5_file(path):
    """CellRanger の .h5 を1つ AnnData で読む（var_names = 遺伝子シンボル）。"""
    import scanpy as sc
    adata = sc.read_10x_h5(str(path))
    adata.var_names_make_unique()
    return adata


def load_10x_h5_per_sample(ds, paths, logger=log):
    """extracted 配下の *.h5 をサンプルごとに読み、結合して返す。"""
    extracted = mf.dataset_extracted_dir(paths, ds)
    h5_files = find_files(extracted, ("*.h5",))
    logger.info("[%s] %s 配下に h5 が %d 個", ds["dataset_id"], extracted, len(h5_files))
    if not h5_files:
        raise RuntimeError(f"{ds['dataset_id']} の .h5 が {extracted} にありません")

    parts = []
    for h5 in h5_files:
        meta = au.parse_geo_filename(h5.name)
        a = read_10x_h5_file(h5)
        ensembl = a.var["gene_ids"].to_numpy() if "gene_ids" in a.var else None
        ftype = a.var["feature_types"].to_numpy() if "feature_types" in a.var else None
        au.ensure_standard_var_columns(a, gene_symbol=a.var_names.to_numpy(),
                                       gene_id=ensembl, ensembl_id=ensembl,
                                       feature_type=ftype)
        a.obs = _base_obs(a.obs_names, ds, sample_id=au.sanitize_id(meta["sample_label"]),
                          sample_label=meta["sample_label"], gsm_id=meta["gsm_id"],
                          source_file=h5.name)
        a.X = sp.csr_matrix(a.X)
        parts.append(a)
        logger.info("  サンプル %-22s cells=%d", au.sanitize_id(meta["sample_label"]), a.n_obs)
    return _concat(parts, ds["source_accession"])


# --------------------------------------------------------------------------
# 10x MTX 三点セット
# --------------------------------------------------------------------------
def _read_features(path: Path):
    rows = []
    with _open_text(path) as fh:
        for line in fh:
            rows.append(line.rstrip("\n").split("\t"))
    df = pd.DataFrame(rows)
    ncol = df.shape[1]
    # 0列目=遺伝子id(ensembl)、1列目=シンボル、2列目=feature type（v3）
    gene_id = df[0].astype(str).to_numpy() if ncol >= 1 else None
    symbol = df[1].astype(str).to_numpy() if ncol >= 2 else gene_id
    ftype = df[2].astype(str).to_numpy() if ncol >= 3 else None
    return gene_id, symbol, ftype


def read_10x_mtx_triplet(mtx_path, barcodes_path, features_path):
    """CellRanger の MTX 三点セット（genes x cells）を読み、cells x genes の AnnData に。"""
    matrix = scipy.io.mmread(str(mtx_path)).tocsr()           # genes x cells
    with _open_text(barcodes_path) as fh:
        barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    gene_id, symbol, ftype = _read_features(features_path)

    adata = ad.AnnData(X=matrix.T.tocsr().astype(np.float32))  # cells x genes に転置
    adata.obs_names = pd.Index([str(b) for b in barcodes])
    var_names = symbol if symbol is not None else gene_id
    adata.var_names = pd.Index([str(s) for s in var_names])
    adata.var_names_make_unique()
    au.ensure_standard_var_columns(adata, gene_symbol=var_names, gene_id=gene_id,
                                   ensembl_id=gene_id, feature_type=ftype)
    return adata


def group_mtx_triplets(extracted: Path) -> dict:
    """barcodes/features/matrix を prefix でグループ化して三点セットを束ねる。"""
    groups: dict = {}
    for path in find_files(extracted, ("*",)):
        meta = au.parse_geo_filename(path.name)
        if meta["kind"] not in ("mtx", "barcodes", "features"):
            continue
        key = str(path.parent / meta["prefix"])
        slot = groups.setdefault(key, {"prefix": meta["prefix"], "gsm": meta["gsm_id"],
                                       "label": meta["sample_label"]})
        slot[meta["kind"]] = path
    return groups


def load_10x_mtx_per_sample(ds, paths, logger=log):
    """extracted 配下の MTX 三点セットをサンプルごとに読み、結合して返す。"""
    extracted = mf.dataset_extracted_dir(paths, ds)
    groups = group_mtx_triplets(extracted)
    complete = {k: v for k, v in groups.items()
                if {"mtx", "barcodes", "features"}.issubset(v)}
    logger.info("[%s] MTX グループ %d 個（完全 %d 個） @ %s",
                ds["dataset_id"], len(groups), len(complete), extracted)
    if not complete:
        raise RuntimeError(
            f"{ds['dataset_id']} の完全な MTX 三点セットが {extracted} にありません; "
            f"見つかった prefix: {[g['prefix'] for g in groups.values()]}")

    parts = []
    for _, slot in sorted(complete.items()):
        a = read_10x_mtx_triplet(slot["mtx"], slot["barcodes"], slot["features"])
        a.obs = _base_obs(a.obs_names, ds, sample_id=au.sanitize_id(slot["label"]),
                          sample_label=slot["label"], gsm_id=slot["gsm"],
                          source_file=Path(slot["mtx"]).name)
        parts.append(a)
        logger.info("  サンプル %-24s cells=%d genes=%d",
                    au.sanitize_id(slot["label"]), a.n_obs, a.n_vars)
    return _concat(parts, ds["source_accession"])
