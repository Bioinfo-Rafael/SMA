"""ノートブック用の dense/text/結合テーブル/処理済み/ネスト tar ローダーと、
R が書き出した中間ファイル（GSE295514）の読み込み。

向きの方針：テキスト行列は「行=遺伝子・列=細胞」と仮定して cells x genes に転置する。
ただし観測した shape と先頭行/列を adata.uns['orientation_report::<file>'] に残し、
人間が確認できるようにする。
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
    # .csv はカンマ、それ以外はタブ区切りとみなす
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
# 汎用 dense リーダー
# --------------------------------------------------------------------------
def read_dense_gene_by_cell_matrix(path):
    """dense 行列（行=遺伝子・列=細胞）を読み、(X cells x genes, genes, cells, 向きレポート) を返す。"""
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
    """サンプルごとの dense テキスト行列を全部読み、結合して返す。"""
    extracted = mf.dataset_extracted_dir(paths, ds)
    files = _text_matrix_files(extracted)
    logger.info("[%s] テキスト行列 %d 個", ds["dataset_id"], len(files))
    if not files:
        raise RuntimeError(f"{ds['dataset_id']} のテキスト行列が {extracted} にありません")
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
    """小さい付随テキスト（行列とは限らない）の先頭をレポートに書き出す。"""
    out = Path(paths["reports"]) / f"{ds['dataset_id']}_side_table_{Path(path).name}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"side table: {path}", ""]
    try:
        with _open_text(path) as fh:
            for i, line in enumerate(fh):
                if i >= n_lines:
                    lines.append(f"... ({n_lines} 行で打ち切り)")
                    break
                lines.append(line.rstrip("\n"))
    except Exception as exc:  # pragma: no cover
        lines.append(f"ERROR: {exc}")
    out.write_text("\n".join(lines))
    logger.info("[%s] 付随テーブルのレポート -> %s", ds["dataset_id"], out)
    return out


def load_mtx_or_text_bundle(ds, paths, logger=log):
    """MTX があれば MTX、無ければ dense テキストとして読む。付随テキストはレポート。"""
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
        logger.info("[%s] MTX 三点セット %d 個", ds["dataset_id"], len(complete))
        return io_10x.load_10x_mtx_per_sample(ds, paths, logger)
    logger.info("[%s] MTX 無し。dense テキストで読みます", ds["dataset_id"])
    return load_dense_or_text_matrix_bundle(ds, paths, logger)


# --------------------------------------------------------------------------
# dense_gene_by_cell_matrix  (GSE167331, TPM)
# --------------------------------------------------------------------------
def load_dense_gene_by_cell_matrix(ds, paths, logger=log):
    """1枚の dense 行列（TPM 等）を読む。data_status は manifest 側の宣言に従う。"""
    raw_dir = mf.dataset_raw_dir(paths, ds)
    fname = ds.get("matrix_file") or ds["files"][0]["name"]
    path = raw_dir / fname
    if not path.exists():
        raise RuntimeError(f"{ds['dataset_id']} の行列ファイルがありません: {path}")
    a = _anndata_from_dense(path, ds, sample_label=ds["source_accession"], gsm_id=au.UNKNOWN)
    logger.info("[%s] 行列ロード: %d cells x %d genes (data_status hint=%s)",
                ds["dataset_id"], a.n_obs, a.n_vars, ds.get("data_status"))
    return a


# --------------------------------------------------------------------------
# combined_umi_tsv_with_metadata  (GSE173524)
# --------------------------------------------------------------------------
def _read_table(path, index_col=0):
    return pd.read_csv(path, sep=_sep_for(path), index_col=index_col)


def read_combined_umi_with_metadata(ds, paths, logger=log):
    """結合 UMI テーブル（genes x cells）を転置し、付随メタデータを obs に join。"""
    raw_dir = mf.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or "GSE173524_umi.tsv.gz"
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"{ds['dataset_id']} の UMI テーブルがありません: {mpath}")

    logger.info("[%s] %s を読み込み（genes x cells -> 転置）", ds["dataset_id"], matrix_file)
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
            logger.warning("  メタデータ欠落: %s", mp)
            continue
        meta = _read_table(mp).add_prefix("meta_")            # 元メタは meta_ 接頭辞で保持
        if entry.get("role") == "per_cell":
            meta.index = meta.index.astype(str)
            a.obs = a.obs.join(meta, how="left").loc[a.obs_names]
            logger.info("  細胞メタを join %s (%d 列)", entry["name"], meta.shape[1])
        else:
            a.uns[f"sample_metadata::{entry['name']}"] = meta.reset_index().astype(str).to_dict("list")
            logger.info("  サンプルメタを uns に格納 %s", entry["name"])
    return a


# --------------------------------------------------------------------------
# processed_count_matrix_with_metadata  (GSE206330, SoupX)
# --------------------------------------------------------------------------
def load_processed_count_matrix_with_metadata(ds, paths, logger=log):
    """処理済み行列（SoupX 補正）+ メタデータ CSV を読み、obs に join。"""
    raw_dir = mf.dataset_raw_dir(paths, ds)
    matrix_file = ds.get("matrix_file") or ds["files"][0]["name"]
    mpath = raw_dir / matrix_file
    if not mpath.exists():
        raise RuntimeError(f"{ds['dataset_id']} の処理済み行列がありません: {mpath}")

    logger.info("[%s] 処理済み行列 %s を読み込み", ds["dataset_id"], matrix_file)
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
            logger.warning("  メタデータ欠落: %s", mp)
            continue
        meta = _read_table(mp).add_prefix("meta_")
        meta.index = meta.index.astype(str)
        a.obs = a.obs.join(meta, how="left").loc[a.obs_names]
        logger.info("  メタを join %s (%d 列)", entry["name"], meta.shape[1])
    return a


# --------------------------------------------------------------------------
# nested_tar_dropseq  (GSE178693)
# --------------------------------------------------------------------------
def load_nested_tar_dropseq(ds, paths, logger=log):
    """ネスト展開済みの中身を MTX 優先で探し、無ければ dense テキストで読む。形式はレポート。"""
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
    logger.info("[%s] 形式: mtx=%d text=%d (レポート -> %s)",
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
                logger.warning("  %s を読めませんでした: %s", path.name, exc)
    if not parts:
        raise RuntimeError(f"{ds['dataset_id']} で読める行列がありません; {report} を参照")
    return io_10x._concat(parts, ds["source_accession"])


# --------------------------------------------------------------------------
# R 中間ファイルの読み込み  (GSE295514, R ノートブックの後)
# --------------------------------------------------------------------------
def read_from_r_intermediate(directory, ds, logger=log):
    """notebooks/R/01_GSE295514_read_rds.ipynb が書き出した
    counts.mtx(genes x cells) / barcodes.csv / genes.csv / metadata.csv から AnnData を組む。"""
    directory = Path(directory)
    mtx = next((directory / n for n in ("counts.mtx.gz", "counts.mtx") if (directory / n).exists()), None)
    if mtx is None:
        raise RuntimeError(
            f"R 中間ファイル counts.mtx が {directory} にありません; "
            "先に notebooks/R/01_GSE295514_read_rds.ipynb を実行してください")
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
        logger.info("[%s] R の meta.data を join (%d 列)", ds["dataset_id"], meta.shape[1])
    return a
