"""ノートブック用の表示・診断ヘルパー（notebooks/python/01・03 で使用）。

すべて「ノートブックに print / 図を出して人間が判断する」ための関数。
cell type / cluster / condition 列を自動で決めることはしない。

含むもの:
  - 構造の確認 : summarize_adata / show_obs_columns / show_var_columns /
                 show_obs_value_counts / show_numeric_obs_summary / show_sample_counts
  - 前処理段階の診断 : infer_processing_state / processing_state_table /
                       gene_stats / plot_value_distribution / plot_gene_stats
  - 名寄せ履歴 : CurationLog
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# 表示でデフォルトに使うサンプル/条件系の列
DEFAULT_SAMPLE_KEYS = [
    "sample_id", "source_accession", "disease_status", "treatment",
    "tissue", "enrichment", "data_status",
]


# ==========================================================================
# 構造の確認
# ==========================================================================
def summarize_adata(adata, name: str | None = None) -> dict:
    """形・X の素性・obs/var 列などを簡潔に print し、dict でも返す。"""
    X = adata.X
    issparse = sp.issparse(X)
    info = {
        "name": name, "n_obs": adata.n_obs, "n_vars": adata.n_vars,
        "X_sparse": issparse, "X_dtype": str(getattr(X, "dtype", "NA")),
        "obs_columns": list(adata.obs.columns), "var_columns": list(adata.var.columns),
        "layers": list(adata.layers.keys()), "obsm": list(adata.obsm.keys()),
        "uns_keys": list(adata.uns.keys()),
    }
    try:
        info["X_min"] = float(X.min()) if issparse else float(np.min(X))
        info["X_max"] = float(X.max()) if issparse else float(np.max(X))
        info["X_sum"] = float(X.sum())
    except Exception as exc:  # pragma: no cover
        info["X_stats_error"] = str(exc)

    head = f"=== {name} ===\n" if name else ""
    print(f"{head}{adata.n_obs} cells x {adata.n_vars} genes | "
          f"X sparse={info['X_sparse']} dtype={info['X_dtype']} "
          f"min={info.get('X_min')} max={info.get('X_max')} sum={info.get('X_sum')}")
    print(f"obs columns ({len(info['obs_columns'])}): {info['obs_columns']}")
    print(f"var columns ({len(info['var_columns'])}): {info['var_columns']}")
    print(f"layers={info['layers']} obsm={info['obsm']} uns={info['uns_keys']}")
    return info


def show_obs_columns(adata) -> pd.DataFrame:
    """obs 各列の dtype / 非NULL数 / 欠損数 / ユニーク数 / 例 を表にして返す。"""
    rows = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        rows.append({"column": col, "dtype": str(s.dtype),
                     "non_null": int(s.notna().sum()), "n_missing": int(s.isna().sum()),
                     "n_unique": int(s.nunique(dropna=True)),
                     "example": str(s.iloc[0]) if len(s) else ""})
    return pd.DataFrame(rows)


def show_var_columns(adata) -> pd.DataFrame:
    """var 各列の概要を表にして返す。"""
    rows = []
    for col in adata.var.columns:
        s = adata.var[col]
        rows.append({"column": col, "dtype": str(s.dtype),
                     "non_null": int(s.notna().sum()),
                     "n_unique": int(s.nunique(dropna=True)),
                     "example": str(s.iloc[0]) if len(s) else ""})
    return pd.DataFrame(rows)


def show_obs_value_counts(adata, max_unique: int = 50, top: int = 30) -> dict:
    """低カーディナリティ/カテゴリ/object な obs 列の value_counts を print。"""
    out: dict = {}
    for col in adata.obs.columns:
        s = adata.obs[col]
        is_cat = str(s.dtype) == "category"
        if s.dtype.kind in ("O", "b") or is_cat or (s.nunique(dropna=True) <= max_unique):
            vc = s.astype(str).value_counts().head(top)
            out[col] = vc
            print(f"\n--- {col} ({s.nunique(dropna=True)} unique) ---")
            print(vc.to_string())
    return out


def show_numeric_obs_summary(adata) -> pd.DataFrame:
    """数値 obs 列の describe を print して返す。"""
    num = adata.obs.select_dtypes(include=[np.number])
    if num.shape[1] == 0:
        print("(数値の obs 列はありません)")
        return pd.DataFrame()
    desc = num.describe().T
    print(desc.to_string())
    return desc


def show_sample_counts(adata, keys=None) -> dict:
    """指定キーごとの細胞数を print。"""
    keys = keys or DEFAULT_SAMPLE_KEYS
    out: dict = {}
    for key in keys:
        if key in adata.obs.columns:
            vc = adata.obs[key].astype(str).value_counts()
            out[key] = vc
            print(f"\n=== cells by {key} ===")
            print(vc.to_string())
    return out


# ==========================================================================
# 前処理段階の診断（どこまで preprocess されたデータか推定する）
# ==========================================================================
def _dense_sample(X, max_cells: int = 2000):
    """大きい行列はランダムでなく先頭 max_cells 行だけ密にして統計に使う。"""
    n = X.shape[0]
    Xs = X[:max_cells] if n > max_cells else X
    return Xs.toarray() if sp.issparse(Xs) else np.asarray(Xs)


def infer_processing_state(adata, layer: str | None = None, max_cells: int = 2000) -> dict:
    """X（または layer）の値から「どこまで前処理されたか」を推定する。

    判定は確定ではなく**ヒント**。整数性・最小値・最大値・負値・行和(library size)
    などから raw count / 正規化 / log / TPM らしさを推測して返す。
    """
    X = adata.layers[layer] if layer else adata.X
    arr = _dense_sample(X, max_cells)
    finite = arr[np.isfinite(arr)]
    nonzero = finite[finite != 0]

    is_int = bool(np.allclose(nonzero, np.round(nonzero))) if nonzero.size else True
    has_negative = bool((finite < 0).any())
    vmax = float(finite.max()) if finite.size else 0.0
    vmin = float(finite.min()) if finite.size else 0.0
    libsize = arr.sum(axis=1)
    libsize_cv = float(np.std(libsize) / (np.mean(libsize) + 1e-9))

    # ざっくりした推定ロジック（あくまで目安）
    if has_negative:
        guess = "scaled/z-scored または log後中心化（負値あり）"
    elif is_int:
        guess = "raw / filtered UMI counts（整数）"
    elif vmax <= 20:
        guess = "log変換後の可能性（非整数・最大値が小さい）"
    elif abs(np.median(libsize[libsize > 0]) - 1e6) / 1e6 < 0.2 if (libsize > 0).any() else False:
        guess = "CPM/TPM の可能性（行和が約1e6）"
    else:
        guess = "非整数（正規化済み or TPM/処理済みの可能性）"

    return {
        "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
        "layer": layer or "X",
        "sampled_cells": int(arr.shape[0]),
        "is_integer": is_int, "has_negative": has_negative,
        "min": vmin, "max": vmax,
        "libsize_median": float(np.median(libsize)) if libsize.size else 0.0,
        "libsize_cv": libsize_cv,
        "declared_data_status": str(adata.obs["data_status"].iloc[0])
        if "data_status" in adata.obs.columns and adata.n_obs else "unknown",
        "guess": guess,
    }


def processing_state_table(adatas: dict, layer: str | None = None) -> pd.DataFrame:
    """複数 AnnData の推定結果を1表にまとめる（宣言した data_status と並べて比較）。"""
    rows = []
    for name, a in adatas.items():
        st = infer_processing_state(a, layer=layer)
        st = {"dataset": name, **st}
        rows.append(st)
    return pd.DataFrame(rows)


def gene_stats(adata, max_cells: int = 2000) -> pd.DataFrame:
    """遺伝子ごとの統計量（平均・分散・発現細胞数・dropout率・最大値）。"""
    X = adata.X
    Xs = X[:max_cells] if adata.n_obs > max_cells else X
    if sp.issparse(Xs):
        mean = np.asarray(Xs.mean(axis=0)).ravel()
        sq = np.asarray(Xs.multiply(Xs).mean(axis=0)).ravel()
        var = sq - mean ** 2
        n_cells = np.asarray((Xs > 0).sum(axis=0)).ravel()
        gmax = np.asarray(Xs.max(axis=0).todense()).ravel()
    else:
        Xs = np.asarray(Xs)
        mean = Xs.mean(axis=0)
        var = Xs.var(axis=0)
        n_cells = (Xs > 0).sum(axis=0)
        gmax = Xs.max(axis=0)
    n = Xs.shape[0]
    df = pd.DataFrame({
        "gene": adata.var_names.astype(str),
        "mean": mean, "var": var,
        "n_cells_expressing": n_cells.astype(int),
        "dropout_rate": 1.0 - n_cells / n,
        "max": gmax,
    })
    return df.sort_values("mean", ascending=False).reset_index(drop=True)


def plot_value_distribution(adata, name: str | None = None, layer: str | None = None,
                            max_cells: int = 2000):
    """値そのもの・library size・遺伝子あたり発現細胞数の分布を描画。

    matplotlib の Figure を返す（保存はノートブック側で）。
    """
    import matplotlib.pyplot as plt

    X = adata.layers[layer] if layer else adata.X
    arr = _dense_sample(X, max_cells)
    nonzero = arr[arr != 0]
    libsize = arr.sum(axis=1)
    genes_per_cell = (arr > 0).sum(axis=1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].hist(nonzero, bins=60)
    axes[0].set_title("nonzero values")
    axes[0].set_yscale("log")
    axes[1].hist(libsize, bins=60)
    axes[1].set_title("library size (row sum)")
    axes[2].hist(genes_per_cell, bins=60)
    axes[2].set_title("genes per cell")
    fig.suptitle(name or "")
    fig.tight_layout()
    return fig


def plot_gene_stats(adata, max_cells: int = 2000, name: str | None = None):
    """mean-variance 関係と dropout-mean 関係を描画（過分散/正規化の見当用）。"""
    import matplotlib.pyplot as plt

    gs = gene_stats(adata, max_cells)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].scatter(gs["mean"] + 1e-6, gs["var"] + 1e-6, s=3, alpha=0.3)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("mean")
    axes[0].set_ylabel("variance")
    axes[0].set_title("mean-variance")
    axes[1].scatter(gs["mean"] + 1e-6, gs["dropout_rate"], s=3, alpha=0.3)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("mean")
    axes[1].set_ylabel("dropout rate")
    axes[1].set_title("dropout-mean")
    fig.suptitle(name or "")
    fig.tight_layout()
    return fig


# ==========================================================================
# 名寄せ履歴（curate での列名・値の変更を記録して CSV 出力）
# ==========================================================================
class CurationLog:
    """各GSEで「どの obs 列名・値をどう変えたか」を記録するロガー。

    使い方（notebooks/python/02）:
        log = CurationLog()
        log.rename_obs(adata, {"orig_celltype": "author_cell_type"}, "GSE287569")
        log.set_constant(adata, "GSE287569", tissue="spinal cord", assay="scRNA-seq")
        log.map_values(adata, "GSE287569", "sample_label",
                       {"WT": "control"}, new_column="disease_status")
        log.export(paths["reports"] / "curation_rename_log.csv")
    """

    def __init__(self):
        self.records: list = []

    # 第2引数 `dataset` はログ上の GSE 識別子。obs 列として dataset_id を
    # **kv で渡せるよう、引数名はあえて dataset_id を避けている。
    def rename_obs(self, adata, mapping: dict, dataset: str):
        """obs 列名を mapping に従ってリネームし、履歴に残す。"""
        for old, new in mapping.items():
            if old in adata.obs.columns:
                adata.obs.rename(columns={old: new}, inplace=True)
                kind = "rename_obs"
            else:
                kind = "rename_obs_MISSING"  # 元列が無かった（要確認）
            self.records.append({"dataset_id": dataset, "kind": kind,
                                 "column": "", "original": old, "new": new})
        return adata

    def set_constant(self, adata, dataset: str, **kv):
        """obs に定数列をセットし、履歴に残す。"""
        for key, val in kv.items():
            adata.obs[key] = val
            self.records.append({"dataset_id": dataset, "kind": "set_constant",
                                 "column": key, "original": "", "new": str(val)})
        return adata

    def map_values(self, adata, dataset: str, column: str, mapping: dict,
                   new_column: str | None = None):
        """column の値を mapping で写像（無い値はそのまま）し、新/同列に入れて記録。"""
        new_column = new_column or column
        adata.obs[new_column] = adata.obs[column].astype(str).map(
            lambda x: mapping.get(x, x)).astype(str)
        for old, new in mapping.items():
            self.records.append({"dataset_id": dataset, "kind": "map_value",
                                 "column": f"{column}->{new_column}",
                                 "original": old, "new": new})
        return adata

    def note(self, dataset: str, **kw):
        """任意のメモを履歴に残す。"""
        self.records.append({"dataset_id": dataset, "kind": "note", **kw})

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.records)

    def export(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.to_frame().to_csv(path, index=False)
        print(f"名寄せ履歴を書き出しました -> {path}  ({len(self.records)} 件)")
        return path
