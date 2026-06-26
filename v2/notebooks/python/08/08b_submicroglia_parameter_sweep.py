#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =====================================================================
# 08b submicroglia parameter sweep
# =====================================================================
# 08 pass2 で得た microglia 再クラスタリング結果
#   v2/results/08_classical_full_inner_microglia_reclustering/
#     04_microglia_reclustering/microglia_classical_reclustered.h5ad
# を読み込み、.obs["microglia_leiden_r1_5"] のうち指定 cluster 群（3, 7-17, 20）だけを抽出して、
# PCA次元 × kNN近傍数 × Leiden resolution の 3x3x3 = 27 通りで再クラスタリング・UMAP・
# dotplot/tracksplot・cluster marker・composition 集計を行い、subcluster 構造の安定性と
# Condition/dataset_id への偏りを確認する。
#
# 重要:
# - h5ad は一切保存しない（選択 subset・各条件の再クラスタリング結果とも in-memory のみ）。
#   出力は plots / dotplots / composition / markers / summary の CSV・PNG のみ。
# - 入力 .X は original-scale の可能性があるため上書きしない。解析用 copy で
#   .X = layers["logexpr_for_clustering"] にしてクラスタリング・plot・marker に使う。
# - rank_genes_groups は探索的 cluster marker 検出であり、sample 単位の condition DEG ではない。
# - HVG は PCA/UMAP/clustering 用のみ。marker 出力は full inner genes 全体で行う。
# - 27 条件をまとめて回すため、各段階を try/except で囲み、失敗条件は summary に
#   status="failed" と error_message を記録して次へ進める。
#
# 実装は 08_classical_full_inner_and_microglia_reclustering.py のヘルパを踏襲（コピー）。
#
# 実行:
#   cd /home/suzuki/Learn/SMA/v2/notebooks/python/08
#   python 08b_submicroglia_parameter_sweep.py

import os
import re
import traceback
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
import scanpy.external as sce

sc.settings.verbosity = 1
warnings.simplefilter("ignore", category=FutureWarning)


# =====================================================================
# 設定（入力 / 出力 / スイープ条件）
# =====================================================================
INPUT_RELPATH = (
    "v2/results/08_classical_full_inner_microglia_reclustering/"
    "04_microglia_reclustering/microglia_classical_reclustered.h5ad"
)
OUT_RELPATH = (
    "v2/results/08_classical_full_inner_microglia_reclustering/"
    "05_selected_microglia_parameter_sweep"
)

SOURCE_CLUSTER_KEY = "microglia_leiden_r1_5"
# 抽出する cluster（3, 7..17, 20）。必ず文字列比較する。
SELECT_CLUSTERS = ["3"] + [str(i) for i in range(7, 18)] + ["20"]

# スイープ条件（3 x 3 x 3 = 27）
PCA_DIMS = [20, 30, 50]
N_NEIGHBORS_LIST = [5, 15, 20]
RESOLUTIONS = [0.5, 1.0, 1.5]

LAYER = "logexpr_for_clustering"
STATE_COL = "qc_preprocessing_state"
NORMALIZE_STATES = {"raw_count_like", "cpm_tpm_like"}
ASIS_STATES = {"log_normalized_like"}

RANDOM_STATE = 0
SCALE_MAX_VALUE = 10
N_TOP_GENES = 3000

BATCH_KEY_CANDIDATES = ["source_accession", "dataset_id"]
CONDITION_CANDIDATES = ["Condition", "condition", "disease", "genotype", "treatment"]

# marker group（mouse 式 Title case。var_names の大文字小文字差は case-insensitive で解決）
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
# full cell-type contamination 確認用
marker_groups_celltype = {
    "Neuron": ["Snap25", "Syt1", "Rbfox1", "Rbfox2", "Rbfox3", "Tubb3", "Map2", "Nefl", "Nefm", "Stmn2", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Slc32a1"],
    "Astrocyte": ["Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Gja1", "Sox9", "S100b", "Gfap", "Agt", "Sparcl1", "Clu"],
    "Reactive_astrocyte": ["Gfap", "Vim", "Serpina3n", "Lcn2", "C3", "Hif3a"],
    "Oligodendrocyte": ["Plp1", "Mbp", "Mog", "Mag", "Mobp", "Cnp", "Mal", "Opalin", "Car2", "Ugt8a", "Ermn"],
    "OPC_oligodendrocyte_precursor": ["Pdgfra", "Cspg4", "Vcan", "Olig1", "Olig2", "Sox10", "Tnr", "Bcas1", "Enpp6", "Tcf7l2"],
    "Endothelial": ["Pecam1", "Cldn5", "Kdr", "Flt1", "Tek", "Ly6c1", "Slco1a4", "Bsg", "Esam", "Vwf"],
    "Pericyte_vascular_mural": ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Notch3", "Vtn", "Des", "Acta2", "Tagln", "Myh11"],
    "Ependymal": ["Foxj1", "Tmem212", "Dnah5", "Dnah12", "Cfap126", "Nnat", "Pifo", "Rsph1"],
    "Schwann_cell_peripheral_nerve_contamination": ["Mpz", "Pmp22", "Mbp", "Plp1", "Sox10", "Ncmap", "Prx"],
    "Meningeal_fibroblast_leptomeningeal": ["Dcn", "Col1a1", "Col1a2", "Lum", "Pdgfra", "Fn1", "Cxcl12", "Pi16"],
    "Monocyte_macrophage_contamination_broad": ["Ccr2", "Ly6c2", "Lyz2", "S100a8", "S100a9", "Ms4a7", "Fcgr1", "Itgam", "Cd14", "Mrc1", "Cd163", "Pf4", "Lyve1", "Folr2"],
    "Tcell_NK": ["Cd3d", "Cd3e", "Trac", "Lck", "Nkg7", "Gzma", "Gzmb", "Klrb1c"],
    "Bcell_plasma": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Jchain", "Mzb1", "Xbp1"],
    "Neutrophil": ["S100a8", "S100a9", "Mpo", "Elane", "Lcn2", "Retnlg", "Cxcr2"],
}
# marker presence 出力で参照する全 group
MARKER_SETS = [
    ("priority1", marker_groups_priority1),
    ("priority2", marker_groups_priority2),
    ("celltype", marker_groups_celltype),
]


# =====================================================================
# ユーティリティ（08 から踏襲）
# =====================================================================
def log(msg: str):
    print(msg, flush=True)


def warn(msg: str):
    print(f"[warn] {msg}", flush=True)


def find_project_root(input_relpath: str = INPUT_RELPATH) -> Path:
    """input h5ad が存在する親ディレクトリを探して SMA root を返す（cwd 非依存）。"""
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
    return re.sub(r"[^0-9A-Za-z]+", "_", str(s)).strip("_")


def res_tag(r) -> str:
    """resolution を列名用の tag にする（0.5->0_5, 1.0->1_0, 1.5->1_5）。"""
    return str(r).replace(".", "_")


def savefig(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def build_upper_map(adata) -> dict:
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


def resolve_present_genes(upper_map: dict, genes):
    """gene のうち var_names に存在するものを実 var_name のリストにする（case-insensitive）。"""
    out = []
    for g in genes:
        vn = upper_map.get(str(g).upper())
        if vn is not None:
            out.append(vn)
    return dedup_keep_order(out)


def first_present(columns, candidates):
    cols = set(map(str, columns))
    for c in candidates:
        if c in cols:
            return c
    return None


def detect_batch_key(adata):
    return next((b for b in BATCH_KEY_CANDIDATES
                 if b in adata.obs.columns and adata.obs[b].astype(str).nunique() > 1), None)


def make_logexpr_layer(adata):
    """qc_preprocessing_state に応じて per-cell 正規化した log-expression を layer に格納する。

    .X（original-scale）は保持する。クラスタリング・可視化・探索的 marker 用であり、
    厳密な count としては使わない。
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
            pass
        else:
            warn(f"未知の state '{st}' は normalize せずそのまま使います。")
        subs.append(sub)
    merged = ad.concat(subs, axis=0, join="outer", merge="same")
    merged = merged[adata.obs_names].copy()
    adata.layers[LAYER] = merged.X


def run_leiden(work, resolution, key_added):
    """leiden を実行（igraph -> default leiden -> louvain の順に fallback。seed 固定）。"""
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


def harmony_integrate_to_obsm(sub, batch_key, label=""):
    """PCA(X_pca) を Harmony で batch 補正し sub.obsm['X_pca_harmony'] (n_obs x n_pcs) を作る。

    scanpy の sce.pp.harmony_integrate を優先。scanpy ラッパは harmonypy の Z_corr を一律 .T
    するため harmonypy 2.x（Z_corr が (n_obs, n_pcs)）では shape 不一致で失敗する。その場合は
    harmonypy を直接呼び、向きを補正して格納する（Harmony 自体は必ず実行）。
    どうしても実行できない場合のみ RuntimeError。
    """
    n = sub.n_obs
    try:
        sce.pp.harmony_integrate(
            sub, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony")
        H = np.asarray(sub.obsm.get("X_pca_harmony"))
        if H.ndim == 2 and H.shape[0] == n:
            return
        warn(f"[{label}] scanpy harmony_integrate の出力 shape={getattr(H, 'shape', None)} が不正。"
             "harmonypy 直接呼び出しに fallback します。")
    except Exception as e:
        warn(f"[{label}] scanpy harmony_integrate が失敗 ({e}); harmonypy 直接呼び出しに fallback します。")

    try:
        import harmonypy
        ho = harmonypy.run_harmony(
            np.asarray(sub.obsm["X_pca"]).astype(np.float64),
            sub.obs, [batch_key], random_state=RANDOM_STATE)
        Z = np.asarray(ho.Z_corr)
        if Z.ndim == 2 and Z.shape[0] != n and Z.shape[1] == n:
            Z = Z.T  # (n_pcs, n_obs) -> (n_obs, n_pcs)
        if Z.ndim != 2 or Z.shape[0] != n:
            raise ValueError(f"unexpected Z_corr shape {np.asarray(ho.Z_corr).shape}")
        sub.obsm["X_pca_harmony"] = np.ascontiguousarray(Z)
    except Exception as e:
        raise RuntimeError(
            f"Harmony integration failed with batch_key={batch_key}. "
            f"Please check that scanpy.external and harmonypy are installed. "
            f"Original error: {e}"
        )


def drop_uns_colors(adata):
    """subset 後に category color 不整合で plot が落ちないよう *_colors を削除する。"""
    for k in list(adata.uns.keys()):
        if str(k).endswith("_colors"):
            del adata.uns[k]


# =====================================================================
# スイープのサブステップ
# =====================================================================
def compute_pca_harmony(selected, pca_requested, batch_key):
    """HVG -> scale -> PCA -> (Harmony) を計算し HVG-subset object を返す。

    返り値: (hvg_sub, actual_n_pcs, n_hvg, harmony_used)
    """
    work = selected.copy()
    work.X = work.layers[LAYER].copy()  # クラスタリングは log-expression で
    nt = int(min(N_TOP_GENES, max(2, work.n_vars - 1)))
    hvg_kwargs = dict(n_top_genes=nt, flavor="seurat")
    if batch_key and batch_key in work.obs.columns and work.obs[batch_key].astype(str).nunique() > 1:
        hvg_kwargs["batch_key"] = batch_key
    sc.pp.highly_variable_genes(work, **hvg_kwargs)
    n_hvg = int(work.var["highly_variable"].sum())

    sub = work[:, work.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=SCALE_MAX_VALUE)
    n_comps = int(max(2, min(pca_requested, sub.n_vars - 1, sub.n_obs - 1)))
    sc.pp.pca(sub, n_comps=n_comps, svd_solver="arpack", random_state=RANDOM_STATE)
    actual_n_pcs = int(sub.obsm["X_pca"].shape[1])

    harmony_used = False
    if batch_key and batch_key in sub.obs.columns and sub.obs[batch_key].astype(str).nunique() > 1:
        harmony_integrate_to_obsm(sub, batch_key, f"pca{pca_requested}")
        harmony_used = "X_pca_harmony" in sub.obsm
    else:
        warn(f"[pca{pca_requested}] batch_key 無し -> Harmony skip（X_pca を使用）")
    log(f"[pca{pca_requested}] n_hvg={n_hvg}, actual_n_pcs={actual_n_pcs}, harmony_used={harmony_used}")
    return sub, actual_n_pcs, n_hvg, harmony_used


def compute_neighbors_umap(hvg_sub, n_neighbors, harmony_used, actual_n_pcs):
    """neighbors + UMAP を計算した HVG-subset copy を返す（resolution 非依存）。"""
    sub = hvg_sub.copy()
    if harmony_used and "X_pca_harmony" in sub.obsm:
        sc.pp.neighbors(sub, n_neighbors=n_neighbors, use_rep="X_pca_harmony",
                        random_state=RANDOM_STATE)
    else:
        sc.pp.neighbors(sub, n_neighbors=n_neighbors, n_pcs=actual_n_pcs,
                        random_state=RANDOM_STATE)
    sc.tl.umap(sub, random_state=RANDOM_STATE)
    return sub


def build_analysis_object(selected, kn_sub, cluster_key, resolution):
    """leiden を実行し、full-gene 解析 object（.X=logexpr）に obsm/obs を転送して返す。"""
    hk = kn_sub.copy()
    run_leiden(hk, resolution, cluster_key)

    ana = selected.copy()  # full inner genes, .X=original-scale, layer 保持
    ana.obsm["X_pca"] = np.asarray(hk.obsm["X_pca"]).copy()
    if "X_pca_harmony" in hk.obsm:
        ana.obsm["X_pca_harmony"] = np.asarray(hk.obsm["X_pca_harmony"]).copy()
    ana.obsm["X_umap"] = np.asarray(hk.obsm["X_umap"]).copy()
    ana.obs[cluster_key] = pd.Categorical(hk.obs[cluster_key].astype(str).values)
    ana.X = ana.layers[LAYER]  # 解析用 copy は logexpr を .X にする（h5ad は保存しない）
    drop_uns_colors(ana)
    return ana


# =====================================================================
# 各条件の出力
# =====================================================================
def export_umaps(ana, cluster_key, plots_dir):
    plots_dir.mkdir(parents=True, exist_ok=True)
    # cluster
    try:
        sc.pl.embedding(ana, basis="umap", color=cluster_key, show=False,
                        legend_loc="on data", title=cluster_key)
        savefig(plots_dir / f"umap_by_{cluster_key}.png")
    except Exception as e:
        warn(f"UMAP(cluster) 失敗 ({cluster_key}): {e}")
        plt.close("all")
    # Condition
    cond = first_present(ana.obs.columns, CONDITION_CANDIDATES)
    if cond is None:
        warn("Condition 系の列が無いため umap_by_Condition をスキップ")
    else:
        try:
            sc.pl.embedding(ana, basis="umap", color=cond, show=False, title=cond)
            savefig(plots_dir / "umap_by_Condition.png")
        except Exception as e:
            warn(f"UMAP(Condition={cond}) 失敗: {e}")
            plt.close("all")
    # dataset_id
    if "dataset_id" not in ana.obs.columns:
        warn("dataset_id 列が無いため umap_by_dataset_id をスキップ")
    else:
        try:
            sc.pl.embedding(ana, basis="umap", color="dataset_id", show=False, title="dataset_id")
            savefig(plots_dir / "umap_by_dataset_id.png")
        except Exception as e:
            warn(f"UMAP(dataset_id) 失敗: {e}")
            plt.close("all")


def build_group_dict(ana, group_dict):
    upper = build_upper_map(ana)
    d = {}
    for grp, genes in group_dict.items():
        present = resolve_present_genes(upper, genes)
        if present:
            d[grp] = present
    return d


def export_dot_tracks(ana, cluster_key, dot_dir):
    dot_dir.mkdir(parents=True, exist_ok=True)
    all_groups = {}
    all_groups.update(marker_groups_priority1)
    all_groups.update(marker_groups_priority2)
    all_groups.update(marker_groups_celltype)
    specs = [
        ("all_markers", all_groups),
        ("priority1", marker_groups_priority1),
        ("priority2", marker_groups_priority2),
    ]
    for label, gd in specs:
        d = build_group_dict(ana, gd)
        if not d:
            warn(f"{label}: dotplot/tracksplot に使える gene が無い")
            continue
        try:
            sc.pl.dotplot(ana, d, groupby=cluster_key, standard_scale="var", show=False)
            savefig(dot_dir / f"dotplot_{label}_by_{cluster_key}.png")
        except Exception as e:
            warn(f"dotplot 失敗 ({label}): {e}")
            plt.close("all")
        try:
            sc.pl.tracksplot(ana, d, groupby=cluster_key, show=False)
            savefig(dot_dir / f"tracksplot_{label}_by_{cluster_key}.png")
        except Exception as e:
            warn(f"tracksplot 失敗 ({label}): {e}")
            plt.close("all")
    # marker presence table（全 marker set, case-insensitive）
    upper = build_upper_map(ana)
    rows = []
    for set_name, gd in MARKER_SETS:
        for grp, genes in gd.items():
            for g in genes:
                vn = upper.get(str(g).upper())
                rows.append({
                    "marker_set": set_name,
                    "group": grp,
                    "gene": g,
                    "present": vn is not None,
                    "matched_var_name": vn if vn is not None else "",
                })
    pd.DataFrame(rows).to_csv(dot_dir / f"marker_presence_by_{cluster_key}.csv", index=False)


def export_markers(ana, cluster_key, mk_dir):
    """探索的 cluster marker（wilcoxon, full inner genes）。condition DEG ではない。"""
    mk_dir.mkdir(parents=True, exist_ok=True)
    if ana.obs[cluster_key].astype(str).nunique() < 2:
        warn(f"{cluster_key} の cluster 数 < 2 のため rank_genes_groups をスキップ")
        return ""
    key = f"rgg_{cluster_key}"
    sc.tl.rank_genes_groups(ana, groupby=cluster_key, method="wilcoxon",
                            use_raw=False, key_added=key)
    df = sc.get.rank_genes_groups_df(ana, group=None, key=key)
    df = df.rename(columns={"group": "cluster", "names": "gene", "scores": "score"})
    cols = ["cluster", "gene", "score", "logfoldchanges", "pvals", "pvals_adj"]
    df = df[[c for c in cols if c in df.columns]]
    path = mk_dir / f"markers_{cluster_key}.csv"
    df.to_csv(path, index=False)
    return str(path)


def save_heatmap(df, path, title):
    """matplotlib のみで heatmap を保存（seaborn 不使用）。x ラベルが長い場合は回転。"""
    if df.shape[0] == 0 or df.shape[1] == 0:
        warn(f"heatmap skip (empty): {path.name}")
        return
    arr = df.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * df.shape[1] + 2),
                                    max(3, 0.4 * df.shape[0] + 2)))
    im = ax.imshow(arr, aspect="auto", cmap="viridis")
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels([str(c) for c in df.columns], rotation=90, fontsize=7)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels([str(i) for i in df.index], fontsize=7)
    ax.set_xlabel(str(df.columns.name) if df.columns.name else "")
    ax.set_ylabel(str(df.index.name) if df.index.name else "")
    ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    savefig(path)


def _crosstab(ana, cluster_key, meta_col):
    ct = pd.crosstab(ana.obs[cluster_key].astype(str), ana.obs[meta_col].astype(str))
    ct.index.name = cluster_key
    ct.columns.name = meta_col
    # cluster id を数値順に近い形で並べる
    try:
        ct = ct.reindex(sorted(ct.index, key=lambda x: (len(x), x)))
    except Exception:
        ct = ct.sort_index()
    return ct


def _fraction_clusterwise(counts):
    return counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)


def _fraction_metawise(counts):
    return counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)


def export_composition(ana, cluster_key, comp_dir):
    comp_dir.mkdir(parents=True, exist_ok=True)

    cond_col = first_present(ana.obs.columns, CONDITION_CANDIDATES)
    if "dataset_id" in ana.obs.columns:
        dataset_col = "dataset_id"
    elif "source_accession" in ana.obs.columns:
        dataset_col = "source_accession"
    else:
        dataset_col = None

    meta_cols = dedup_keep_order([c for c in [cond_col, dataset_col] if c])
    if not meta_cols:
        warn("Condition/dataset 系の列が無いため composition をスキップ（cluster_sizes のみ）")

    for meta in meta_cols:
        counts = _crosstab(ana, cluster_key, meta)
        fcw = _fraction_clusterwise(counts)
        fmw = _fraction_metawise(counts)
        base = f"composition_{cluster_key}_by_{meta}"
        counts.to_csv(comp_dir / f"{base}_counts.csv")
        fcw.to_csv(comp_dir / f"{base}_fraction_clusterwise.csv")
        fmw.to_csv(comp_dir / f"{base}_fraction_metawise.csv")
        save_heatmap(counts, comp_dir / f"{base}_counts_heatmap.png",
                     f"{cluster_key} x {meta} (counts)")
        save_heatmap(fcw, comp_dir / f"{base}_fraction_clusterwise_heatmap.png",
                     f"{cluster_key} x {meta} (cluster-wise fraction)")
        save_heatmap(fmw, comp_dir / f"{base}_fraction_metawise_heatmap.png",
                     f"{cluster_key} x {meta} (meta-wise fraction)")

    # 06 互換ファイル（fraction は cluster-wise）
    if cond_col is not None:
        counts = _crosstab(ana, cluster_key, cond_col)
        counts.to_csv(comp_dir / f"composition_{cluster_key}_by_Condition_counts.csv")
        _fraction_clusterwise(counts).to_csv(
            comp_dir / f"composition_{cluster_key}_by_Condition_fraction.csv")
    if "dataset_id" in ana.obs.columns:
        counts = _crosstab(ana, cluster_key, "dataset_id")
        counts.to_csv(comp_dir / f"composition_{cluster_key}_by_dataset_id_counts.csv")
        _fraction_clusterwise(counts).to_csv(
            comp_dir / f"composition_{cluster_key}_by_dataset_id_fraction.csv")

    # cluster size table
    sizes = ana.obs[cluster_key].astype(str).value_counts()
    sizes = sizes.reindex(sorted(sizes.index, key=lambda x: (len(x), x)))
    sizes.rename_axis(cluster_key).to_frame("n_cells").to_csv(
        comp_dir / f"cluster_sizes_{cluster_key}.csv")


def run_one_condition(selected, kn_sub, cluster_key, resolution, cond_dir):
    """1 条件分の出力（h5ad は保存しない）。返り値: (n_clusters, marker_csv_path)。"""
    ana = build_analysis_object(selected, kn_sub, cluster_key, resolution)
    n_clusters = int(ana.obs[cluster_key].astype(str).nunique())

    export_umaps(ana, cluster_key, cond_dir / "plots")
    export_dot_tracks(ana, cluster_key, cond_dir / "dotplots")
    marker_csv = export_markers(ana, cluster_key, cond_dir / "markers")
    export_composition(ana, cluster_key, cond_dir / "composition")
    return n_clusters, marker_csv


# =====================================================================
# main
# =====================================================================
def main():
    root = find_project_root(INPUT_RELPATH)
    INPUT = (root / INPUT_RELPATH).resolve()
    OUT = (root / OUT_RELPATH).resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    summary_path = OUT / "summary_all_parameter_sets.csv"

    log("=" * 70)
    log("08b submicroglia parameter sweep")
    log("=" * 70)
    log(f"PROJECT_ROOT : {root}")
    log(f"INPUT        : {INPUT}")
    log(f"OUT          : {OUT}")

    # --- load ---
    adata = sc.read_h5ad(INPUT)
    log(f"input shape = {adata.shape}")
    if SOURCE_CLUSTER_KEY not in adata.obs.columns:
        raise KeyError(
            f"入力に cluster 列 '{SOURCE_CLUSTER_KEY}' がありません。"
            f" obs columns: {list(adata.obs.columns)}"
        )

    # --- select clusters (必ず文字列比較) ---
    cl = adata.obs[SOURCE_CLUSTER_KEY].astype(str)
    available = set(cl.unique())
    present_sel = [c for c in SELECT_CLUSTERS if c in available]
    missing_sel = [c for c in SELECT_CLUSTERS if c not in available]
    if missing_sel:
        warn(f"選択 cluster のうち入力に存在しないもの: {missing_sel}")
    mask = cl.isin(SELECT_CLUSTERS).values
    n_sel = int(mask.sum())
    if n_sel == 0:
        raise RuntimeError(
            f"選択 cluster {SELECT_CLUSTERS} に該当する細胞が 0 です。"
            f" 入力の cluster: {sorted(available, key=lambda x: (len(x), x))}"
        )
    selected = adata[mask].copy()
    selected.obs[SOURCE_CLUSTER_KEY] = selected.obs[SOURCE_CLUSTER_KEY].astype(str).astype("category")
    drop_uns_colors(selected)
    log(f"selected: {selected.shape}  (clusters used: {present_sel})")

    if LAYER not in selected.layers:
        warn(f"{LAYER} layer が無いため作成します（.X は保持）。")
        make_logexpr_layer(selected)

    batch_key = detect_batch_key(selected)
    log(f"batch_key = {batch_key}")

    n_cells = int(selected.n_obs)
    n_genes = int(selected.n_vars)

    # --- 27 条件スイープ ---
    summary_rows = []
    for pca in PCA_DIMS:
        hvg_sub = None
        actual_n_pcs = None
        n_hvg = None
        harmony_used = False
        pca_err = None
        try:
            hvg_sub, actual_n_pcs, n_hvg, harmony_used = compute_pca_harmony(selected, pca, batch_key)
        except Exception as e:
            pca_err = f"PCA/Harmony failed: {type(e).__name__}: {e}"
            warn(f"[pca{pca}] {pca_err}")
            traceback.print_exc()

        for k in N_NEIGHBORS_LIST:
            kn_sub = None
            kn_err = pca_err
            if hvg_sub is not None:
                try:
                    kn_sub = compute_neighbors_umap(hvg_sub, k, harmony_used, actual_n_pcs)
                except Exception as e:
                    kn_err = f"neighbors/umap failed: {type(e).__name__}: {e}"
                    warn(f"[pca{pca}/knn{k:02d}] {kn_err}")
                    traceback.print_exc()

            for res in RESOLUTIONS:
                cluster_key = f"submicro_pca{pca}_knn{k:02d}_res{res_tag(res)}"
                cond_dir = OUT / f"pca{pca}" / f"knn{k:02d}" / f"res{res_tag(res)}"
                row = {
                    "pca_requested": pca,
                    "actual_n_pcs": actual_n_pcs,
                    "n_neighbors": k,
                    "resolution": res,
                    "cluster_key": cluster_key,
                    "n_cells": n_cells,
                    "n_genes": n_genes,
                    "n_hvg": n_hvg,
                    "n_clusters": None,
                    "batch_key": batch_key if batch_key else "",
                    "harmony_used": harmony_used,
                    "output_dir": str(cond_dir),
                    "marker_csv_path": "",
                    "status": "ok",
                    "error_message": "",
                }
                if kn_err is not None or kn_sub is None:
                    row["status"] = "failed"
                    row["error_message"] = kn_err or "upstream stage failed"
                    summary_rows.append(row)
                    log(f"[skip] {cluster_key}: {row['error_message']}")
                    continue
                try:
                    log(f"[run] {cluster_key}")
                    n_clusters, marker_csv = run_one_condition(
                        selected, kn_sub, cluster_key, res, cond_dir)
                    row["n_clusters"] = n_clusters
                    row["marker_csv_path"] = marker_csv
                except Exception as e:
                    row["status"] = "failed"
                    row["error_message"] = f"{type(e).__name__}: {e}"
                    warn(f"[fail] {cluster_key}: {row['error_message']}")
                    traceback.print_exc()
                summary_rows.append(row)
                # 逐次保存（途中で落ちても summary が残るように）
                pd.DataFrame(summary_rows).to_csv(summary_path, index=False)

    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    n_ok = sum(1 for r in summary_rows if r["status"] == "ok")
    log("\n" + "=" * 70)
    log(f"08b 完了: {n_ok}/{len(summary_rows)} 条件が ok")
    log(f"summary: {summary_path}")
    log(f"出力先: {OUT}")
    log("=" * 70)


if __name__ == "__main__":
    main()
