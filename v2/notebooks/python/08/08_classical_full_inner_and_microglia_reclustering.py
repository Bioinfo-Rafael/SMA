#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =====================================================================
# 08 classical full inner clustering and microglia reclustering
# =====================================================================
# 04d で作成した full inner-gene AnnData (merged_qc_original_scale_inner.h5ad) を入力に、
# scVI を使わない古典的 Scanpy workflow（log-normalization -> HVG -> scale -> PCA -> kNN ->
# UMAP -> Leiden）で全細胞クラスタリングを行い、marker gene を確認したうえで人手で
# cell type annotation を行うためのスクリプト。
#
# その後、手動 annotation に基づき microglia-like cluster を抽出し、その細胞だけを
# 再クラスタリングして microglia / DAM / IFN / MHC-II / stress / contamination 等の
# 細かい状態を再確認する。
#
# これは 2-pass 設計のスクリプトであり、1回目 / 2回目は --pass 引数で明示的に指定する。
#   --pass 1 : 全細胞クラスタリング + marker 出力 + 手動 annotation 用 CSV(template) 生成
#   （人間が template を編集し manual_annotation_filled_full_clustering.csv として保存）
#   --pass 2 : filled CSV を読み込み、microglia-like subset を抽出して再クラスタリング
#
# -------------------------------------------------------------------
# 解釈上の重要な注意（コメント）:
# 1. merged_qc_original_scale_inner.h5ad は full inner-gene object だが、.X は original-scale で
#    あり raw_count_like / cpm_tpm_like / log_normalized_like が混在している可能性がある。
# 2. このスクリプトは .X と .var を保持し、クラスタリング・可視化用に
#    layers["logexpr_for_clustering"] を作成する。
# 3. HVG selection は PCA/UMAP/clustering のためだけに使う。
# 4. marker gene 抽出は HVG だけでなく full inner genes 上で行う。
# 5. rank_genes_groups は探索的な cluster marker 検出であり、sample 単位の condition DEG ではない。
# 6. microglia-like subset の選択は、出力された cluster marker と marker plot に基づく
#    意図的な手動選択である。
# 7. 手動選択後、microglia-like subset を再クラスタリングして細かい subpopulation を調べる。
# -------------------------------------------------------------------
#
# 実行（repository root から、または v2/notebooks/python/08 から、あるいは SMA_ROOT 指定で任意の場所から）:
#   python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
#   （template を手動で埋めて filled として保存したのち）
#   python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
#
# project root はスクリプト内で自動探索するため cwd には依存しない。
# v2/notebooks/python/08/v2/results/... のような二重 path は作らない。

import os
import re
import argparse
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # ヘッドレスで figure を保存するため
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

sc.settings.verbosity = 1
warnings.simplefilter("ignore", category=FutureWarning)


# =====================================================================
# 設定
# =====================================================================
INPUT_RELPATH = "v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad"
# results に作る出力フォルダは step 番号を先頭に付ける
OUT_RELPATH = "v2/results/08_classical_full_inner_microglia_reclustering"

LAYER = "logexpr_for_clustering"
STATE_COL = "qc_preprocessing_state"
NORMALIZE_STATES = {"raw_count_like", "cpm_tpm_like"}
ASIS_STATES = {"log_normalized_like"}

RANDOM_STATE = 0
SCALE_MAX_VALUE = 10

N_TOP_GENES = 3000
N_PCS = 30
LEIDEN_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

# dotplot / tracksplot / template で使う代表 resolution（存在するものだけ使う）
SELECTED_FULL_COLS = ["leiden_r0_6", "leiden_r1_0"]
SELECTED_MICRO_COLS = ["microglia_leiden_r0_6", "microglia_leiden_r1_0"]

# microglia-like subset 選択時の fallback keyword（case-insensitive）
MICROGLIA_KEYWORDS = ["microglia", "dam", "myeloid", "macrophage"]
MIN_MICRO_CELLS = 30

# plot に使う metadata 候補（存在するものだけ plot）
META_PLOT_COLS = [
    "source_accession", "dataset_id", "qc_preprocessing_state",
    "condition", "Condition", "disease", "genotype", "treatment",
]
# HVG batch_key 候補（存在し、かつ 2 水準以上あれば使う）
BATCH_KEY_CANDIDATES = ["source_accession", "dataset_id"]

# marker group（mouse 式 Title case。var_names の大文字小文字が違う可能性があるため case-insensitive 解決）
marker_groups_priority1 = {
    "DAM_core": ["Apoe", "Tyrobp", "Trem2", "Gpnmb", "Cst7", "Lpl", "C1qa", "C1qb", "C1qc"],
    "Homeostatic_microglia": ["P2ry12", "Tmem119", "Cx3cr1", "Sall1"],
    "DAM_activated": ["Apoe", "Trem2", "Tyrobp", "Gpnmb", "Lpl", "Cst7", "Cd68"],
    "Complement": ["C1qa", "C1qb", "C1qc"],
}
marker_groups_priority2 = {
    "Microglia_identity": ["Hexb", "Fcrls", "Olfml3", "Gpr34", "P2ry13", "Siglech", "Slc2a5", "Csf1r"],
    "Homeostatic_support": ["Tgfbr1", "Mef2c", "Maf", "Selplg", "Sparc"],
    "DAM_support": ["Itgax", "Axl", "Clec7a", "Cd9", "Cd63", "Spp1", "Lgals3", "Csf1", "C3ar1"],
    "Phagocytic_lysosomal": ["Ctsb", "Ctsd", "Ctsz", "Lamp1", "Fcgr3", "Mertk"],
    "IFN": ["Ifit1", "Ifit2", "Ifit3", "Isg15", "Irf7", "Stat1", "Mx1", "Oasl2", "Usp18", "Cxcl10"],
    "MHC_II": ["H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "Cd74", "B2m", "H2-K1", "H2-D1", "Cd83", "Cd86"],
    "Proliferation": ["Mki67", "Top2a", "Birc5", "Mcm5", "Stmn1"],
    "Stress": ["Fos", "Jun", "Junb", "Egr1", "Dusp1", "Atf3", "Hspa1a", "Hspa1b"],
    "Monocyte_macrophage_contamination": [
        "Ccr2", "Ly6c2", "Lyz2", "S100a8", "S100a9", "Ms4a7", "Fcgr1", "Mrc1", "Lyve1", "Pf4", "Cd163"
    ],
}
MARKER_PRIORITIES = [
    ("priority1", marker_groups_priority1),
    ("priority2", marker_groups_priority2),
]


# =====================================================================
# ユーティリティ関数
# =====================================================================
def log(msg: str):
    print(msg, flush=True)


def warn(msg: str):
    print(f"[warn] {msg}", flush=True)


def find_project_root(input_relpath: str = INPUT_RELPATH) -> Path:
    """input h5ad が存在する親ディレクトリを探して SMA root を返す（cwd 非依存）。

    探索順: 環境変数 SMA_ROOT -> __file__ の親群 -> cwd の親群。
    """
    candidates = []
    env = os.environ.get("SMA_ROOT")
    if env:
        candidates.append(Path(env).expanduser().resolve())
    try:
        here = Path(__file__).resolve()
        candidates.append(here.parent)
        candidates.extend(here.parents)
    except NameError:
        pass
    cwd = Path.cwd().resolve()
    candidates.append(cwd)
    candidates.extend(cwd.parents)

    seen = set()
    checked = []
    for c in candidates:
        if c in seen:
            continue
        seen.add(c)
        checked.append(c)
        if (c / input_relpath).exists():
            return c

    msg = [
        "SMA project root を自動検出できませんでした。",
        f"探したファイル: {input_relpath}",
        "",
        "対処: SMA リポジトリ内で実行するか、環境変数 SMA_ROOT を指定してください。",
        "  export SMA_ROOT=/path/to/SMA",
        "",
        "確認した root 候補:",
    ]
    msg.extend([f"  - {p}" for p in checked[:30]])
    raise FileNotFoundError("\n".join(msg))


def sanitize_filename(s) -> str:
    """ファイル名用に英数以外を _ に置換する。"""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(s)).strip("_")


def res_tag(r) -> str:
    """resolution を列名用の tag にする（0.2->0_2, 1.0->1_0）。"""
    return str(r).replace(".", "_")


def savefig(path):
    """現在の matplotlib figure を保存する。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def build_upper_map(adata) -> dict:
    """var_names を {大文字: 実際の var_name} の dict にする（case-insensitive 解決用）。"""
    m = {}
    for v in adata.var_names:
        m.setdefault(str(v).upper(), str(v))
    return m


def dedup_keep_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def flatten_marker_groups(*group_dicts) -> list:
    """複数の marker group dict から重複なしの gene リストを作る。"""
    genes = []
    for gd in group_dicts:
        for _, gl in gd.items():
            genes.extend(gl)
    return dedup_keep_order(genes)


def resolve_present_genes(upper_map: dict, genes):
    """gene のうち var_names に存在するものを実 var_name のリストにして返す。"""
    out = []
    for g in genes:
        vn = upper_map.get(str(g).upper())
        if vn is not None:
            out.append(vn)
    return dedup_keep_order(out)


def marker_presence_table(adata) -> pd.DataFrame:
    """marker group ごとに var_names 上の存在を case-insensitive に判定する。"""
    upper = build_upper_map(adata)
    rows = []
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            for g in genes:
                vn = upper.get(str(g).upper())
                rows.append({
                    "priority": priority,
                    "group": group,
                    "gene": g,
                    "present": vn is not None,
                    "matched_var_name": vn if vn is not None else "",
                })
    return pd.DataFrame(rows)


def parse_bool(x) -> bool:
    return str(x).strip().lower() in {"true", "1", "yes", "y", "t", "include"}


def make_logexpr_layer(adata):
    """qc_preprocessing_state に応じて per-cell 正規化した log-expression を layer に格納する。

    - raw_count_like / cpm_tpm_like : normalize_total(1e4) -> log1p
    - log_normalized_like           : そのまま
    - その他 / 列が無い               : 警告のうえ normalize_total(1e4)+log1p（列が無い場合は全細胞）
    この layer はクラスタリング・UMAP・marker 可視化・探索的 marker 検出のためのものであり、
    厳密な count データとしては使わない。.X（original-scale）は保持する。
    """
    if STATE_COL not in adata.obs.columns:
        warn(f"{STATE_COL} 列が無いため、全細胞に normalize_total(1e4)+log1p を適用して "
             f"{LAYER} layer を作成します（厳密な count としては使わないこと）。")
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        adata.layers[LAYER] = tmp.X
        return

    states = pd.Series(adata.obs[STATE_COL].astype(str).values, index=adata.obs_names)
    log(f"[logexpr] {STATE_COL}: {dict(states.value_counts())}")
    subs = []
    for st in pd.unique(states.values):
        sub = adata[(states.values == st)].copy()
        if st in NORMALIZE_STATES:
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
        elif st in ASIS_STATES:
            pass  # そのまま使う
        else:
            warn(f"未知の state '{st}' は normalize せずそのまま使います。")
        subs.append(sub)

    merged = ad.concat(subs, axis=0, join="outer", merge="same")
    merged = merged[adata.obs_names].copy()  # 元の細胞順に戻す
    adata.layers[LAYER] = merged.X
    adata.uns["logexpr_layer_note_08"] = (
        "layer logexpr_for_clustering: per-cell normalize_total(1e4)+log1p for "
        "raw_count_like/cpm_tpm_like; log_normalized_like kept as-is; "
        "for clustering/visualization only, NOT strict counts. .X kept original-scale."
    )


def run_leiden(work, resolution, key_added):
    """leiden を実行する（igraph -> default leiden -> louvain の順に fallback。seed 固定）。"""
    try:
        sc.tl.leiden(work, resolution=resolution, key_added=key_added,
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=RANDOM_STATE)
    except Exception as e1:
        warn(f"leiden(flavor=igraph) 失敗 ({e1}); default leiden を試す")
        try:
            sc.tl.leiden(work, resolution=resolution, key_added=key_added,
                         random_state=RANDOM_STATE)
        except Exception as e2:
            warn(f"leiden 失敗 ({e2}); louvain に fallback")
            sc.tl.louvain(work, resolution=resolution, key_added=key_added,
                          random_state=RANDOM_STATE)


def run_classical_pipeline(work, n_top_genes, n_pcs, resolutions, leiden_prefix,
                           batch_key=None, label=""):
    """古典的 Scanpy workflow（HVG -> scale -> PCA -> kNN -> UMAP -> Leiden）。

    work.X は log-expression（logexpr layer の値）である前提。HVG に subset して
    scale/PCA/UMAP/Leiden を行い、その HVG-subset object を返す（呼び出し側が
    obsm/obs を full-gene object に転送する）。
    返り値: (clustered_hvg, leiden_cols)
    """
    nt = int(min(n_top_genes, max(2, work.n_vars - 1)))
    hvg_kwargs = dict(n_top_genes=nt, flavor="seurat")
    if batch_key and batch_key in work.obs.columns and work.obs[batch_key].astype(str).nunique() > 1:
        hvg_kwargs["batch_key"] = batch_key
        log(f"[{label}] HVG batch_key = {batch_key}")
    sc.pp.highly_variable_genes(work, **hvg_kwargs)
    n_hvg = int(work.var["highly_variable"].sum())
    log(f"[{label}] HVG = {n_hvg} / {work.n_vars}")

    sub = work[:, work.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=SCALE_MAX_VALUE)
    n_comps = int(min(n_pcs, sub.n_vars - 1, sub.n_obs - 1))
    n_comps = max(2, n_comps)
    sc.pp.pca(sub, n_comps=n_comps, svd_solver="arpack", random_state=RANDOM_STATE)
    use_pcs = int(min(n_comps, sub.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=use_pcs, random_state=RANDOM_STATE)
    sc.tl.umap(sub, random_state=RANDOM_STATE)

    leiden_cols = []
    for r in resolutions:
        key = f"{leiden_prefix}{res_tag(r)}"
        run_leiden(sub, r, key)
        leiden_cols.append(key)
        log(f"[{label}] {key}: {sub.obs[key].astype(str).nunique()} clusters")
    return sub, leiden_cols


def transfer_clustering(target, source, leiden_cols):
    """HVG-subset の clustering 結果を full-gene object へ転送する（細胞順は同一前提）。"""
    target.obsm["X_pca"] = np.asarray(source.obsm["X_pca"]).copy()
    target.obsm["X_umap"] = np.asarray(source.obsm["X_umap"]).copy()
    for c in leiden_cols:
        target.obs[c] = pd.Categorical(source.obs[c].astype(str).values)
    for k in ("neighbors", "umap", "pca"):
        if k in source.uns:
            target.uns[k] = source.uns[k]
    for k in ("distances", "connectivities"):
        if k in source.obsp:
            target.obsp[k] = source.obsp[k]


def export_rank_genes_groups(adata_markers, cluster_col, outdir, fname):
    """rank_genes_groups（wilcoxon, use_raw=False, full inner genes）を CSV 保存して df を返す。

    探索的 cluster marker 検出。condition DEG ではない。
    """
    if cluster_col not in adata_markers.obs.columns:
        warn(f"{cluster_col} が obs に無いため rank_genes_groups をスキップ")
        return None
    if adata_markers.obs[cluster_col].astype(str).nunique() < 2:
        warn(f"{cluster_col} の cluster 数 < 2 のため rank_genes_groups をスキップ")
        return None
    key = f"rgg_{cluster_col}"
    sc.tl.rank_genes_groups(adata_markers, groupby=cluster_col, method="wilcoxon",
                            use_raw=False, key_added=key)
    df = sc.get.rank_genes_groups_df(adata_markers, group=None, key=key)
    df = df.rename(columns={"group": "cluster", "names": "gene", "scores": "score"})
    cols = ["cluster", "gene", "score", "logfoldchanges", "pvals", "pvals_adj"]
    df = df[[c for c in cols if c in df.columns]]
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / fname, index=False)
    log(f"[markers] saved {outdir.name}/{fname}")
    return df


def export_umap_plots(adata, color_cols, outdir, prefix="umap"):
    """UMAP を指定 obs 列で着色して保存する（存在しない列は skip）。"""
    if "X_umap" not in adata.obsm:
        warn("X_umap が無いため UMAP plot をスキップ")
        return
    outdir.mkdir(parents=True, exist_ok=True)
    for c in dedup_keep_order(color_cols):
        if c not in adata.obs.columns:
            continue
        try:
            on_data = c.startswith("leiden") or c.startswith("microglia_leiden")
            sc.pl.embedding(
                adata, basis="umap", color=c, show=False, title=c,
                legend_loc="on data" if on_data else "right margin",
            )
            savefig(outdir / f"{prefix}_{sanitize_filename(c)}.png")
        except Exception as e:
            warn(f"UMAP plot 失敗 ({c}): {e}")
            plt.close("all")


def export_marker_feature_umaps(adata, outdir):
    """priority1 / priority2 の存在 marker gene を group 単位で feature UMAP にする。"""
    if "X_umap" not in adata.obsm:
        warn("X_umap が無いため marker feature UMAP をスキップ")
        return
    upper = build_upper_map(adata)
    fdir = outdir / "marker_feature"
    fdir.mkdir(parents=True, exist_ok=True)
    for gname, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            present = resolve_present_genes(upper, genes)
            if not present:
                continue
            try:
                sc.pl.embedding(adata, basis="umap", color=present, color_map="viridis",
                                ncols=4, show=False, frameon=False)
                savefig(fdir / f"{gname}_{sanitize_filename(group)}.png")
            except Exception as e:
                warn(f"marker feature UMAP 失敗 ({gname}/{group}): {e}")
                plt.close("all")


def export_dotplots_tracksplots(adata, group_dict, cluster_col, outdir, label):
    """marker group dict について dotplot / tracksplot を cluster 列で保存する。"""
    if cluster_col not in adata.obs.columns:
        return
    upper = build_upper_map(adata)
    gdict = {}
    for group, genes in group_dict.items():
        present = resolve_present_genes(upper, genes)
        if present:
            gdict[group] = present
    if not gdict:
        warn(f"{label}: dotplot/tracksplot に使える gene が無い")
        return
    outdir.mkdir(parents=True, exist_ok=True)
    for fn, plotter, scaled in [("dotplot", sc.pl.dotplot, True),
                                ("tracksplot", sc.pl.tracksplot, False)]:
        try:
            if scaled:
                plotter(adata, gdict, groupby=cluster_col, standard_scale="var", show=False)
            else:
                plotter(adata, gdict, groupby=cluster_col, show=False)
            savefig(outdir / f"{fn}_{label}_by_{sanitize_filename(cluster_col)}.png")
        except Exception as e:
            warn(f"{fn} 失敗 ({label}, {cluster_col}): {e}")
            plt.close("all")


def export_manual_annotation_template(adata, cluster_cols, marker_df_by_col, outpath,
                                      analysis_level, kind):
    """cluster 列 × cluster ごとに 1 行の手動 annotation template CSV を作る。

    kind="full"  -> include_for_microglia_recluster 列を含む
    kind="micro" -> cluster_interpretation 列を含む
    """
    rows = []
    for col in cluster_cols:
        if col not in adata.obs.columns:
            continue
        vc = adata.obs[col].astype(str).value_counts()
        df = marker_df_by_col.get(col) if marker_df_by_col else None
        for cid in sorted(vc.index, key=lambda x: (len(x), x)):
            top = ""
            if df is not None:
                top = ", ".join(
                    df[df["cluster"].astype(str) == cid]["gene"].head(15).astype(str).tolist())
            row = {
                "analysis_level": analysis_level,
                "cluster_col": col,
                "cluster_id": cid,
                "n_cells": int(vc[cid]),
                "top_markers": top,
                "manual_annotation": "",
            }
            if kind == "full":
                row["include_for_microglia_recluster"] = ""
            else:
                row["cluster_interpretation"] = ""
            row["notes"] = ""
            rows.append(row)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(outpath, index=False)
    log(f"[template] {analysis_level}: {len(rows)} rows -> {outpath}")


def load_and_apply_manual_annotation(adata, filled_df):
    """filled CSV の手動 annotation を obs に戻す。

    cluster_col ごとに manual_annotation_<col> / include_for_microglia_recluster_<col> を追加。
    返り値: 適用できた cluster_col のリスト。
    """
    if "cluster_col" not in filled_df.columns:
        warn("filled CSV に cluster_col 列が無いため annotation を適用できません。")
        return []
    df = filled_df[filled_df["cluster_col"].notna()].copy()
    applied = []
    for col in df["cluster_col"].astype(str).unique():
        if col not in adata.obs.columns:
            warn(f"cluster_col '{col}' が obs に無いためスキップ")
            continue
        sub = df[df["cluster_col"].astype(str) == col]
        ann_map = {str(c): ("" if pd.isna(a) else str(a))
                   for c, a in zip(sub["cluster_id"], sub.get("manual_annotation", ""))}
        adata.obs[f"manual_annotation_{col}"] = (
            adata.obs[col].astype(str).map(ann_map).fillna("").astype(str).values)
        if "include_for_microglia_recluster" in sub.columns:
            inc_map = {str(c): parse_bool(v)
                       for c, v in zip(sub["cluster_id"], sub["include_for_microglia_recluster"])}
            adata.obs[f"include_for_microglia_recluster_{col}"] = (
                adata.obs[col].astype(str).map(inc_map).fillna(False).astype(bool).values)
        applied.append(col)
        log(f"[annotation] applied for {col}")
    return applied


def select_microglia_like_cells(adata, applied_cols, keywords):
    """include_for_microglia_recluster==TRUE で選択。空なら manual_annotation の keyword で fallback。

    返り値: (mask: np.ndarray[bool], summary: DataFrame)
    """
    final = pd.Series(False, index=adata.obs_names)
    summary = []
    pat = "|".join([re.escape(k) for k in keywords])
    for col in applied_cols:
        inc_col = f"include_for_microglia_recluster_{col}"
        ann_col = f"manual_annotation_{col}"
        method = "include_column"
        if inc_col in adata.obs.columns:
            mask = adata.obs[inc_col].astype(bool)
        else:
            mask = pd.Series(False, index=adata.obs_names)
        if int(mask.sum()) == 0 and ann_col in adata.obs.columns:
            mask = adata.obs[ann_col].astype(str).str.lower().str.contains(pat, na=False, regex=True)
            method = "keyword_fallback"
        final = final | mask.values
        summary.append({
            "cluster_col": col,
            "selection_method": method,
            "n_selected": int(mask.sum()),
        })
    summary.append({
        "cluster_col": "__union_all_cols__",
        "selection_method": "union",
        "n_selected": int(final.sum()),
    })
    return final.values, pd.DataFrame(summary)


def present_meta_cols(adata, candidates):
    return [c for c in candidates if c in adata.obs.columns]


def present_cluster_cols(adata, preferred, fallback):
    sel = [c for c in preferred if c in adata.obs.columns]
    return sel if sel else [c for c in fallback if c in adata.obs.columns]


def note_overwrite(path: Path):
    if path.exists():
        log(f"[overwrite] 既存ファイルを上書きします: {path}")


# =====================================================================
# パス設定
# =====================================================================
def setup_paths() -> SimpleNamespace:
    """SMA root を自動検出し、step 番号付きの出力フォルダ構成を組み立てる。"""
    root = find_project_root(INPUT_RELPATH)
    out = (root / OUT_RELPATH).resolve()
    P = SimpleNamespace(
        root=root,
        input=(root / INPUT_RELPATH).resolve(),
        out=out,
        reports=out / "01_reports",
        full=out / "02_full_clustering",
        full_plots=out / "02_full_clustering" / "plots",
        full_markers=out / "02_full_clustering" / "marker_genes",
        ann=out / "03_manual_annotation",
        micro=out / "04_microglia_reclustering",
        micro_plots=out / "04_microglia_reclustering" / "plots",
        micro_markers=out / "04_microglia_reclustering" / "marker_genes",
    )
    P.full_clustered = P.full / "full_inner_classical_clustered.h5ad"
    P.template_full = P.ann / "manual_annotation_template_full_clustering.csv"
    P.filled_full = P.ann / "manual_annotation_filled_full_clustering.csv"
    P.full_with_ann = P.ann / "full_inner_with_manual_annotation.h5ad"
    P.micro_like = P.ann / "microglia_like_from_manual_annotation.h5ad"
    P.micro_summary = P.ann / "microglia_selection_summary.csv"
    P.template_micro = P.ann / "manual_annotation_template_microglia_reclustering.csv"
    P.micro_reclustered = P.micro / "microglia_classical_reclustered.h5ad"
    return P


def ensure_dirs(P):
    for d in (P.out, P.reports, P.full, P.full_plots, P.full_markers,
              P.ann, P.micro, P.micro_plots, P.micro_markers):
        d.mkdir(parents=True, exist_ok=True)


def detect_batch_key(adata):
    return next((b for b in BATCH_KEY_CANDIDATES
                 if b in adata.obs.columns and adata.obs[b].astype(str).nunique() > 1), None)


# =====================================================================
# Pass 1: 全細胞 classical clustering + marker + template
# =====================================================================
def run_pass1(P):
    log("\n" + "=" * 70)
    log("PASS 1: 全細胞 classical clustering + marker + 手動 annotation template")
    log("=" * 70)

    # --- Part 1: load and inspect ---
    log("\n[Part 1] load and inspect")
    adata = sc.read_h5ad(P.input)
    log(f"  adata.shape = {adata.shape}")

    report = []
    report.append("08 input QC report")
    report.append("=" * 60)
    report.append(f"input: {P.input}")
    report.append(f"adata.shape: {adata.shape} (cells x genes)")
    report.append("")
    report.append(f"obs columns ({adata.obs.shape[1]}):")
    report.append(", ".join(map(str, adata.obs.columns)))
    report.append("")
    report.append(f"var columns ({adata.var.shape[1]}):")
    report.append(", ".join(map(str, adata.var.columns)))
    report.append("")
    report.append(f"obsm keys: {list(adata.obsm.keys())}")
    report.append(f"uns keys : {list(adata.uns.keys())}")
    report.append("")
    for col in [STATE_COL, "source_accession", "dataset_id"]:
        if col in adata.obs.columns:
            report.append(f"[{col}] value_counts:")
            report.append(adata.obs[col].astype(str).value_counts().to_string())
            report.append("")
        else:
            report.append(f"[{col}] (列なし)")
            report.append("")
    cond_like = [c for c in adata.obs.columns
                 if any(k in c.lower() for k in ["condition", "disease", "genotype", "treatment"])]
    report.append(f"condition-like columns: {cond_like}")
    for c in cond_like:
        report.append(f"[{c}] value_counts:")
        report.append(adata.obs[c].astype(str).value_counts().to_string())
        report.append("")
    mpt = marker_presence_table(adata)
    report.append("marker gene presence (full inner var_names):")
    report.append(f"  present {int(mpt['present'].sum())} / {len(mpt)}")
    g = mpt.groupby("group")["present"].agg(["sum", "count"])
    for grp, rrow in g.iterrows():
        report.append(f"    {grp}: {int(rrow['sum'])}/{int(rrow['count'])}")
    (P.reports / "08_input_qc_report.txt").write_text("\n".join(report), encoding="utf-8")
    log(f"  saved: {P.reports / '08_input_qc_report.txt'}")
    mpt.to_csv(P.reports / "marker_presence_full_inner.csv", index=False)
    log(f"  saved: {P.reports / 'marker_presence_full_inner.csv'}")

    # --- Part 2: logexpr layer (.X は保持) ---
    log("\n[Part 2] build logexpr_for_clustering layer (do not overwrite .X)")
    make_logexpr_layer(adata)

    # --- Part 3: full object classical clustering ---
    log("\n[Part 3] full object classical clustering")
    batch_key = detect_batch_key(adata)
    adata_clust = adata.copy()
    adata_clust.X = adata.layers[LAYER].copy()  # クラスタリングは log-expression で
    clustered_hvg, leiden_cols = run_classical_pipeline(
        adata_clust, N_TOP_GENES, N_PCS, LEIDEN_RESOLUTIONS, "leiden_r",
        batch_key=batch_key, label="full")
    del adata_clust
    transfer_clustering(adata, clustered_hvg, leiden_cols)
    adata.uns["clustering_note_08"] = (
        "classical: logexpr layer -> HVG(seurat) -> scale -> PCA(n_pcs=30) -> kNN -> UMAP -> Leiden; "
        ".X kept original-scale; layer logexpr_for_clustering used for clustering/visualization."
    )
    note_overwrite(P.full_clustered)
    adata.write_h5ad(P.full_clustered)
    log(f"  saved: {P.full_clustered}")

    # --- Part 4: marker gene extraction (full inner genes) ---
    log("\n[Part 4] marker gene extraction (full inner genes, exploratory)")
    afm = adata.copy()
    afm.X = adata.layers[LAYER].copy()
    full_marker_df = {}
    for col in leiden_cols:
        df = export_rank_genes_groups(afm, col, P.full_markers, f"markers_{col}.csv")
        if df is not None:
            full_marker_df[col] = df

    # --- Part 5: full clustering plots ---
    log("\n[Part 5] full clustering plots")
    meta_cols = present_meta_cols(afm, META_PLOT_COLS)
    export_umap_plots(afm, meta_cols + leiden_cols, P.full_plots, prefix="umap")
    export_marker_feature_umaps(afm, P.full_plots)
    sel_full_cols = present_cluster_cols(afm, SELECTED_FULL_COLS, leiden_cols)
    for col in sel_full_cols:
        export_dotplots_tracksplots(afm, marker_groups_priority1, col, P.full_plots, "priority1")
        export_dotplots_tracksplots(afm, marker_groups_priority2, col, P.full_plots, "priority2")

    # --- Part 6: manual annotation template ---
    log("\n[Part 6] manual annotation template (full clustering)")
    note_overwrite(P.template_full)
    export_manual_annotation_template(
        afm, sel_full_cols, full_marker_df, P.template_full,
        analysis_level="full_clustering", kind="full")
    del afm

    log("\n" + "-" * 70)
    log("PASS 1 完了。次の手順:")
    log(f"  1. {P.template_full} を開き、cluster ごとに")
    log("     manual_annotation / include_for_microglia_recluster / notes を埋める")
    log(f"  2. {P.filled_full.name} という名前で保存する")
    log("  3. もう一度 --pass 2 で実行する:")
    log("       python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2")


# =====================================================================
# Pass 2: 手動 annotation 適用 + microglia-like subset 再クラスタリング
# =====================================================================
def run_pass2(P):
    log("\n" + "=" * 70)
    log("PASS 2: 手動 annotation 適用 + microglia-like subset 再クラスタリング")
    log("=" * 70)

    if not P.full_clustered.exists():
        raise FileNotFoundError(
            f"{P.full_clustered} が見つかりません。先に --pass 1 を実行してください:\n"
            "  python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1"
        )
    if not P.filled_full.exists():
        log(
            "Manual annotation file not found. Full clustering is complete. "
            "Fill manual_annotation_template_full_clustering.csv and save it as "
            "manual_annotation_filled_full_clustering.csv, then rerun this script with "
            "--pass 2 to perform microglia-like subset re-clustering."
        )
        log(f"  期待した場所: {P.filled_full}")
        return

    # --- Part 7: apply manual annotation and select microglia-like cells ---
    log("\n[Part 7] apply manual annotation and select microglia-like cells")
    adata = sc.read_h5ad(P.full_clustered)
    if LAYER not in adata.layers:
        warn(f"{LAYER} layer が無いため再作成します。")
        make_logexpr_layer(adata)
    batch_key = detect_batch_key(adata)

    filled_df = pd.read_csv(P.filled_full)
    applied_cols = load_and_apply_manual_annotation(adata, filled_df)
    if not applied_cols:
        warn("適用できる cluster_col がありませんでした。Part 7-10 をスキップします。")
        return

    note_overwrite(P.full_with_ann)
    adata.write_h5ad(P.full_with_ann)
    log(f"  saved: {P.full_with_ann}")

    mask, sel_summary = select_microglia_like_cells(adata, applied_cols, MICROGLIA_KEYWORDS)
    sel_summary.to_csv(P.micro_summary, index=False)
    log(f"  saved: {P.micro_summary}")
    n_micro = int(mask.sum())
    log(f"  microglia-like 選択細胞数 = {n_micro}")

    if n_micro < MIN_MICRO_CELLS:
        warn(f"microglia-like 細胞が少なすぎます ({n_micro} < {MIN_MICRO_CELLS})。"
             " Part 8-10 をスキップします。template / include 列を確認してください。")
        return

    micro = adata[mask].copy()
    note_overwrite(P.micro_like)
    micro.write_h5ad(P.micro_like)
    log(f"  saved: {P.micro_like}")

    # --- Part 8: microglia-like subset re-clustering ---
    log("\n[Part 8] microglia-like subset re-clustering")
    if LAYER not in micro.layers:
        warn(f"{LAYER} layer が無いため microglia subset で再作成します。")
        make_logexpr_layer(micro)
    mc = micro.copy()
    mc.X = micro.layers[LAYER].copy()
    micro_hvg, micro_cols = run_classical_pipeline(
        mc, min(N_TOP_GENES, mc.n_vars), min(N_PCS, max(2, mc.n_obs - 1)),
        LEIDEN_RESOLUTIONS, "microglia_leiden_r", batch_key=batch_key, label="micro")
    del mc
    transfer_clustering(micro, micro_hvg, micro_cols)
    micro.uns["clustering_note_08"] = "microglia reclustering (classical); .X original-scale."
    note_overwrite(P.micro_reclustered)
    micro.write_h5ad(P.micro_reclustered)
    log(f"  saved: {P.micro_reclustered}")

    # --- Part 9: microglia marker genes and plots ---
    log("\n[Part 9] microglia marker genes and plots")
    mfm = micro.copy()
    mfm.X = micro.layers[LAYER].copy()
    micro_marker_df = {}
    for col in micro_cols:
        df = export_rank_genes_groups(mfm, col, P.micro_markers, f"markers_{col}.csv")
        if df is not None:
            micro_marker_df[col] = df
    micro_meta = present_meta_cols(mfm, META_PLOT_COLS)
    export_umap_plots(mfm, micro_cols + micro_meta, P.micro_plots, prefix="umap")
    export_marker_feature_umaps(mfm, P.micro_plots)
    sel_micro_cols = present_cluster_cols(mfm, SELECTED_MICRO_COLS, micro_cols)
    for col in sel_micro_cols:
        export_dotplots_tracksplots(mfm, marker_groups_priority1, col, P.micro_plots, "priority1")
        export_dotplots_tracksplots(mfm, marker_groups_priority2, col, P.micro_plots, "priority2")

    # --- Part 10: manual annotation template for microglia subclusters ---
    log("\n[Part 10] manual annotation template (microglia reclustering)")
    note_overwrite(P.template_micro)
    export_manual_annotation_template(
        mfm, sel_micro_cols, micro_marker_df, P.template_micro,
        analysis_level="microglia_reclustering", kind="micro")
    del mfm

    log("\n" + "-" * 70)
    log("PASS 2 完了。microglia subcluster の結果を確認し、")
    log(f"  {P.template_micro}")
    log("  を開いて microglia subcluster を手動 annotation してください。")


# =====================================================================
# main
# =====================================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="08 classical full inner clustering and microglia reclustering (2-pass).",
        epilog=(
            "例:\n"
            "  python 08_classical_full_inner_and_microglia_reclustering.py --pass 1\n"
            "    -> 全細胞クラスタリング + marker + 手動 annotation template を作成\n"
            "  （template を編集して *_filled_full_clustering.csv として保存）\n"
            "  python 08_classical_full_inner_and_microglia_reclustering.py --pass 2\n"
            "    -> 手動 annotation を適用し microglia-like subset を再クラスタリング"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--pass", dest="run_pass", required=True, choices=["1", "2"],
        help="1=全細胞クラスタリング+template生成 / 2=手動annotation適用+microglia再クラスタリング",
    )
    return ap.parse_args()


def main():
    args = parse_args()
    P = setup_paths()
    ensure_dirs(P)

    log("=" * 70)
    log("08 classical full inner clustering and microglia reclustering")
    log("=" * 70)
    log(f"PROJECT_ROOT : {P.root}")
    log(f"INPUT_H5AD   : {P.input}")
    log(f"OUT_DIR      : {P.out}")
    log(f"PASS         : {args.run_pass}")

    if args.run_pass == "1":
        run_pass1(P)
    else:
        run_pass2(P)

    log("\n" + "=" * 70)
    log("08 完了")
    log("=" * 70)
    log(f"出力先: {P.out}")


if __name__ == "__main__":
    main()
