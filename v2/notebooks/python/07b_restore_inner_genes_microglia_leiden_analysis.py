# %% [markdown]
# # 07b. inner 遺伝子復元 + microglia_leiden_r05 全cluster解析
#
# 07 は target cluster（"6"）中心の解析だったが、07b では microglia subclustering の
# cluster 列（`microglia_leiden_cluster_r05` / `microglia_leiden_r05`）を **全 cluster 横断**の
# グループキーとして使い、04d 由来の full inner genes（~8863）上で
# marker 可視化（dotplot / tracksplot / UMAP feature）・DEG・pseudo-bulk を作り直す。
# Condition（ALS など）を使った比較も含める。
#
# 前提（07 と同じ）:
# - full inner `.X` は original-scale（raw/cpm/log 混在）。HVG/annotation ソースは 06 の
#   microglia-subclustered AnnData。var_names は大文字化 → marker は case-insensitive 解決。
# - joined は full inner X（microglia 細胞のみ ~66850）＋ microglia 再クラスタリング annotation。
# - 可視化・探索的 DEG は per-cell 正規化した logexpr copy を使い、pseudo-bulk は raw_count_like の
#   original-scale を使う。
#
# 07 本体は変更しない。出力は `v2/results/07b_inner_fullgenes_microglia_leiden_analysis/`。
#
# 実行（SMA リポジトリのルートから、または SMA_ROOT 指定で任意ディレクトリから）:
# ```bash
# python v2/notebooks/python/07b_restore_inner_genes_microglia_leiden_analysis.py
# ```

# %%
import os
import re
import copy as _copy
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

import matplotlib
matplotlib.use("Agg")  # ヘッドレスで figure を保存
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

sc.settings.verbosity = 1
warnings.simplefilter("ignore", category=FutureWarning)


# %%
# =====================================================================
# パス（07 と同じ root 自動検出。cwd 非依存・SMA_ROOT 対応・二重 path 防止）
# =====================================================================
REQUIRED_INPUT_RELPATHS = [
    "v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad",
    "v2/results/microglia_subclustering/adata_microglia_subclustered.h5ad",
]


def candidate_roots():
    roots = []
    env_root = os.environ.get("SMA_ROOT")
    if env_root:
        roots.append(Path(env_root).expanduser().resolve())
    try:
        here = Path(__file__).resolve()
        roots.append(here.parent)
        roots.extend(here.parents)
    except NameError:
        pass
    cwd = Path.cwd().resolve()
    roots.append(cwd)
    roots.extend(cwd.parents)
    out, seen = [], set()
    for r in roots:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out


def looks_like_sma_root(root: Path) -> bool:
    return all((root / rel).exists() for rel in REQUIRED_INPUT_RELPATHS)


def find_project_root() -> Path:
    checked = []
    for root in candidate_roots():
        checked.append(root)
        if looks_like_sma_root(root):
            return root
    msg = [
        "SMA project root を自動検出できませんでした。",
        "",
        "対処: SMA リポジトリ内で実行するか、環境変数 SMA_ROOT を指定してください。",
        "  export SMA_ROOT=/home/suzuki/Learn/SMA",
        "",
        "確認した root 候補:",
    ]
    msg.extend([f"  - {p}" for p in checked[:30]])
    msg.append("")
    msg.append("各候補で必要だった入力ファイル:")
    for rel in REQUIRED_INPUT_RELPATHS:
        msg.append(f"  - {rel}")
    raise FileNotFoundError("\n".join(msg))


PROJECT_ROOT = find_project_root()


def rpath(rel: str) -> Path:
    return (PROJECT_ROOT / rel).resolve()


FULL_INNER_PATH = rpath("v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad")
HVG_RESULT_PATH = rpath("v2/results/microglia_subclustering/adata_microglia_subclustered.h5ad")
OUT_DIR = rpath("v2/results/07b_inner_fullgenes_microglia_leiden_analysis")

JOINED_PATH = OUT_DIR / "inner_fullgenes_with_microglia_leiden.h5ad"
LOGEXPR_JOINED_PATH = OUT_DIR / "inner_fullgenes_logexpr_with_microglia_leiden.h5ad"

FIG_DIR = OUT_DIR / "figures"
FIG_UMAP_DIR = FIG_DIR / "umap_feature"
FIG_DOT_DIR = FIG_DIR / "dotplot"
FIG_TRACKS_DIR = FIG_DIR / "tracksplot"
CLUSTER_SUMMARY_DIR = OUT_DIR / "cluster_summary"
DEG_DIR = OUT_DIR / "deg"
PSEUDOBULK_DIR = OUT_DIR / "pseudobulk"

for d in (OUT_DIR, FIG_DIR, FIG_UMAP_DIR, FIG_DOT_DIR, FIG_TRACKS_DIR,
          CLUSTER_SUMMARY_DIR, DEG_DIR, PSEUDOBULK_DIR):
    d.mkdir(parents=True, exist_ok=True)


# %%
# =====================================================================
# 設定（07 から流用 + 07b 用に cluster 列候補を microglia 優先に）
# =====================================================================
STATE_COL = "qc_preprocessing_state"
NORMALIZE_STATES = {"raw_count_like", "cpm_tpm_like"}
ASIS_STATES = {"log_normalized_like"}

# cluster グループキー候補（ユーザー指定の microglia_leiden_cluster_r05 を最優先で自動検出）
CLUSTER_COL_CANDIDATES = [
    "microglia_leiden_cluster_r05", "hvg_microglia_leiden_cluster_r05",
    "microglia_leiden_r05", "hvg_microglia_leiden_r05",
    "microglia_leiden_r0_5", "hvg_microglia_leiden_r0_5",
    "leiden_harmony_r05", "hvg_leiden_harmony_r05",
    "leiden_before_harmony_r05", "hvg_leiden_before_harmony_r05",
]
UMAP_BASIS_CANDIDATES = [
    "umap_microglia", "X_umap_microglia",
    "umap_after_harmony", "X_umap_after_harmony",
    "umap", "X_umap",
]
CONDITION_CANDIDATES = ["Condition", "condition", "disease", "genotype", "treatment"]
SAMPLE_COL_CANDIDATES = [
    "sample_id", "sample_label", "gsm_id", "donor_id",
    "animal_id", "mouse_id", "orig.ident", "source_file",
]
ALS_TOKENS = ["ALS", "als", "SOD1", "TDP", "disease", "Disease"]
MIN_PB_SAMPLES = 2   # pseudo-bulk Welch の各群最低サンプル数

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

STATE = {"join_key": None, "n_common_cells": None, "cluster_col": None,
         "umap_basis": None, "condition_col": None, "sample_col": None,
         "sample_col_provisional": False, "n_clusters": 0}


# %%
# =====================================================================
# ユーティリティ（07 から逐語）
# =====================================================================
def warn(msg: str):
    print(f"[warn] {msg}")


def san(s) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", str(s)).strip("_")


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


def resolve_genes(upper_map: dict, genes):
    found, missing = {}, []
    for g in genes:
        vn = upper_map.get(str(g).upper())
        if vn is not None:
            found[g] = vn
        else:
            missing.append(g)
    return found, missing


def dedup_keep_order(seq):
    seen, out = set(), []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def resolve_umap_basis(adata, candidates):
    for cand in candidates:
        base = cand[2:] if cand.startswith("X_") else cand
        if ("X_" + base) in adata.obsm or base in adata.obsm:
            return base
    return None


def first_present(columns, candidates):
    cols = set(map(str, columns))
    for c in candidates:
        if c in cols:
            return c
    return None


def obs_columns_overview(adata) -> pd.DataFrame:
    rows = []
    for c in adata.obs.columns:
        s = adata.obs[c]
        try:
            nuniq = int(s.nunique(dropna=True))
        except TypeError:
            nuniq = -1
        try:
            example = ", ".join(map(str, pd.Series(s.dropna().unique()[:5])))
        except Exception:
            example = ""
        rows.append({"column": c, "dtype": str(s.dtype), "n_unique": nuniq,
                     "n_missing": int(s.isna().sum()), "example_values": example})
    return pd.DataFrame(rows)


def marker_presence_table(adata) -> pd.DataFrame:
    upper = build_upper_map(adata)
    rows = []
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            for g in genes:
                vn = upper.get(str(g).upper())
                rows.append({"priority": priority, "marker_group": group, "gene_query": g,
                             "found": vn is not None, "matched_var_name": vn if vn else ""})
    return pd.DataFrame(rows)


def make_logexpr_copy(adata):
    """qc_preprocessing_state に応じて per-cell 正規化した copy を返す（sparse 維持、.X 非破壊）。"""
    out = adata.copy()
    if STATE_COL not in out.obs.columns:
        warn(f"{STATE_COL} 列が無いため、正規化せずそのまま log-expr copy とします。")
        return out
    states = pd.Series(out.obs[STATE_COL].astype(str).values, index=out.obs_names)
    print(f"[logexpr] {STATE_COL}: {dict(states.value_counts())}")
    subs, did_norm = [], False
    for st in list(pd.unique(states.values)):
        sub = out[(states.values == st)].copy()
        if st in NORMALIZE_STATES:
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            did_norm = True
        elif st in ASIS_STATES:
            pass
        else:
            warn(f"未知の state '{st}' はそのままにします。")
        subs.append(sub)
    merged = ad.concat(subs, axis=0, join="outer", merge="same")
    merged = merged[out.obs_names].copy()
    merged.var = out.var.copy()
    merged.uns = _copy.deepcopy(out.uns)
    for k in out.obsm.keys():
        if k not in merged.obsm:
            merged.obsm[k] = np.asarray(out.obsm[k]).copy()
    if did_norm:
        merged.uns["log1p"] = {"base": None}
    merged.uns["logexpr_note_07b"] = (
        "per-cell normalize_total(1e4)+log1p for raw_count_like/cpm_tpm_like; "
        "log_normalized_like kept as-is; exploratory visualization/DEG only.")
    return merged


def detect_als_values(series):
    vals = pd.Series(series.astype(str)).dropna()
    toks = [t.lower() for t in ALS_TOKENS]
    return [u for u in sorted(vals.unique()) if any(t in u.lower() for t in toks)]


# --- pseudo-bulk（07 から逐語） ---
def pseudobulk_sum(adata_raw, label_series):
    labels = pd.Series(label_series).astype(str)
    valid = labels.notna().values & (labels.values != "nan")
    A = adata_raw[valid]
    lab = labels.values[valid]
    cats = pd.Categorical(lab)
    codes = cats.codes
    n = A.n_obs
    ncat = len(cats.categories)
    M = sparse.csr_matrix((np.ones(n, dtype=np.float64), (codes, np.arange(n))), shape=(ncat, n))
    X = A.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    summed = np.asarray((M @ X).todense())
    counts_df = pd.DataFrame(summed.T, index=A.var_names.astype(str), columns=list(cats.categories))
    n_cells = pd.Series(np.asarray(M.sum(axis=1)).ravel(), index=list(cats.categories), name="n_cells")
    return counts_df, n_cells


def pseudobulk_metadata(adata_raw, composite_series, factor_cols, extra_maps=None):
    obs = adata_raw.obs.copy()
    obs["_pb"] = pd.Series(composite_series).astype(str).values
    rows = []
    for pb, sub in obs.groupby("_pb"):
        row = {"pb_sample": pb, "n_cells": int(len(sub))}
        for c in factor_cols:
            if c in sub.columns:
                m = sub[c].astype(str).mode()
                row[c] = m.iat[0] if len(m) else ""
        if extra_maps:
            for col, mp in extra_maps.items():
                row[col] = mp.get(pb, "")
        rows.append(row)
    return pd.DataFrame(rows).set_index("pb_sample")


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    out[ok] = q
    return out


def join_full_and_hvg(full, hvg):
    """07 の join ロジック（cell_uid/obs_names overlap でキー選択、共通細胞にそろえ obs/obsm/uns 移植）。"""
    def key_index(adata, key):
        if key == "obs_names":
            return pd.Index(adata.obs_names.astype(str))
        return pd.Index(adata.obs[key].astype(str))

    obs_overlap = len(pd.Index(full.obs_names.astype(str)).intersection(
        pd.Index(hvg.obs_names.astype(str))))
    cands = []
    if "cell_uid" in full.obs.columns and "cell_uid" in hvg.obs.columns:
        uo = len(pd.Index(full.obs["cell_uid"].astype(str)).intersection(
            pd.Index(hvg.obs["cell_uid"].astype(str))))
        cands.append(("cell_uid", uo))
    cands.append(("obs_names", obs_overlap))
    join_key, join_overlap = max(cands, key=lambda kv: kv[1])
    if join_overlap <= 0:
        raise RuntimeError("full inner と microglia result の間で対応細胞が見つかりません。")
    fk, hk = key_index(full, join_key), key_index(hvg, join_key)
    if not fk.is_unique or not hk.is_unique:
        raise RuntimeError(f"対応キー {join_key} が unique でないため結合できません。")
    common = pd.Index(dedup_keep_order(list(hk[hk.isin(set(fk))])))
    fpos = pd.Series(np.arange(full.n_obs), index=fk).loc[common].to_numpy()
    hpos = pd.Series(np.arange(hvg.n_obs), index=hk).loc[common].to_numpy()

    joined = full[fpos].copy()   # full の .X（8863 genes）を維持
    hview = hvg[hpos]            # view（HVG .X はコピーしない）
    for col in hvg.obs.columns:
        vals = hview.obs[col].values
        if col not in joined.obs.columns:
            joined.obs[col] = vals
        else:
            try:
                same = np.array_equal(pd.Series(joined.obs[col].values).astype(str).to_numpy(),
                                      pd.Series(vals).astype(str).to_numpy())
            except Exception:
                same = False
            if not same:
                joined.obs["hvg_" + col] = vals
    for k in list(hvg.obsm.keys()):
        joined.obsm[k] = np.asarray(hview.obsm[k]).copy()
    for k, v in hvg.uns.items():
        newk = k if k not in joined.uns else ("hvg_" + k)
        try:
            joined.uns[newk] = _copy.deepcopy(v)
        except Exception as e:
            warn(f"uns['{k}'] をコピーできませんでした: {e}")
    joined.uns["restore_note_07b"] = (
        f"X=04d original-scale inner genes; annotation from microglia result; "
        f"join_key={join_key}; n_common_cells={len(common)}")
    return joined, join_key, int(len(common))


# %% [markdown]
# ## 1. 確認パート

# %%
for p in (FULL_INNER_PATH, HVG_RESULT_PATH):
    if not p.exists():
        raise FileNotFoundError(f"入力が見つかりません: {p}")

print("loading full inner ...")
full = sc.read_h5ad(FULL_INNER_PATH)
print("loading microglia result ...")
hvg = sc.read_h5ad(HVG_RESULT_PATH)
print(f"full inner: {full.shape}  microglia: {hvg.shape}")

rep = []
rep.append("07b inspection report")
rep.append("=" * 60)
rep.append(f"FULL_INNER_PATH : {FULL_INNER_PATH}")
rep.append(f"HVG_RESULT_PATH : {HVG_RESULT_PATH}")
rep.append(f"full inner shape: {full.shape}")
rep.append(f"microglia shape : {hvg.shape}")
rep.append("")
rep.append(f"microglia obs columns: {list(hvg.obs.columns)}")
rep.append(f"microglia obsm keys  : {list(hvg.obsm.keys())}")
rep.append("")
for c in CLUSTER_COL_CANDIDATES + CONDITION_CANDIDATES:
    if c in hvg.obs.columns:
        rep.append(f"  candidate present: {c}  (n_unique={hvg.obs[c].nunique()})")
(OUT_DIR / "00_inspection_report.txt").write_text("\n".join(rep), encoding="utf-8")
marker_presence_table(full).to_csv(OUT_DIR / "00_marker_presence_full_inner.csv", index=False)
marker_presence_table(hvg).to_csv(OUT_DIR / "00_marker_presence_hvg.csv", index=False)
print("saved: 00_inspection_report.txt / 00_marker_presence_*.csv")


# %% [markdown]
# ## 2. 結合 + logexpr copy

# %%
joined, join_key, n_common = join_full_and_hvg(full, hvg)
STATE["join_key"] = join_key
STATE["n_common_cells"] = n_common
print(f"[join] key={join_key}, joined shape={joined.shape}")
joined.write_h5ad(JOINED_PATH)
print("saved:", JOINED_PATH)

log_adata = make_logexpr_copy(joined)
log_adata.write_h5ad(LOGEXPR_JOINED_PATH)
print("saved:", LOGEXPR_JOINED_PATH)


# %%
# --- cluster 列 / basis / condition / sample 列の決定 ---
cluster_col = first_present(log_adata.obs.columns, CLUSTER_COL_CANDIDATES)
if cluster_col is None:
    raise RuntimeError(
        f"cluster グループ列が見つかりません（候補: {CLUSTER_COL_CANDIDATES}）。"
        f" obs columns: {list(log_adata.obs.columns)}")
STATE["cluster_col"] = cluster_col
log_adata.obs[cluster_col] = log_adata.obs[cluster_col].astype(str).astype("category")
joined.obs[cluster_col] = joined.obs[cluster_col].astype(str).astype("category")
STATE["n_clusters"] = int(log_adata.obs[cluster_col].nunique())

umap_basis = resolve_umap_basis(log_adata, UMAP_BASIS_CANDIDATES)
STATE["umap_basis"] = umap_basis
condition_col = first_present(log_adata.obs.columns, CONDITION_CANDIDATES)
STATE["condition_col"] = condition_col
print(f"[setup] cluster_col={cluster_col} (n_clusters={STATE['n_clusters']}), "
      f"umap_basis={umap_basis}, condition_col={condition_col}")


# %% [markdown]
# ## 3. marker 可視化（全 cluster、full inner genes）

# %%
upper_log = build_upper_map(log_adata)

# 欠落 marker
missing_rows = []
for priority, groups in MARKER_PRIORITIES:
    for group, genes in groups.items():
        _, missing = resolve_genes(upper_log, genes)
        for g in missing:
            missing_rows.append({"priority": priority, "marker_group": group, "gene_query": g})
pd.DataFrame(missing_rows).to_csv(FIG_DIR / "marker_missing_in_logexpr.csv", index=False)

# UMAP feature plot（group ごと multipanel）
if umap_basis is not None:
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            f, _ = resolve_genes(upper_log, genes)
            vns = dedup_keep_order(list(f.values()))
            if not vns:
                continue
            try:
                sc.pl.embedding(log_adata, basis=umap_basis, color=vns, color_map="viridis",
                                ncols=4, show=False, frameon=False)
                plt.suptitle(f"{priority}: {group}")
                savefig(FIG_UMAP_DIR / f"{priority}_{san(group)}.png")
            except Exception as e:
                warn(f"UMAP feature 失敗 ({priority}/{group}): {e}")
                plt.close("all")
    # cluster で着色した UMAP も
    try:
        sc.pl.embedding(log_adata, basis=umap_basis, color=cluster_col, show=False,
                        legend_loc="on data", title=cluster_col)
        savefig(FIG_UMAP_DIR / f"umap_by_{san(cluster_col)}.png")
    except Exception as e:
        warn(f"UMAP(cluster) 失敗: {e}")
        plt.close("all")
    for meta in [condition_col, "dataset_id", "source_accession"]:
        if meta and meta in log_adata.obs.columns:
            try:
                sc.pl.embedding(log_adata, basis=umap_basis, color=meta, show=False, title=meta)
                savefig(FIG_UMAP_DIR / f"umap_by_{san(meta)}.png")
            except Exception as e:
                warn(f"UMAP({meta}) 失敗: {e}")
                plt.close("all")
else:
    warn("UMAP basis が無いため UMAP feature plot をスキップ")

# dotplot / tracksplot（priority ごと、cluster_col でグループ化）
def group_dict(groups):
    d = {}
    for group, genes in groups.items():
        f, _ = resolve_genes(upper_log, genes)
        vns = dedup_keep_order(list(f.values()))
        if vns:
            d[group] = vns
    return d


for priority, groups in MARKER_PRIORITIES:
    gd = group_dict(groups)
    if not gd:
        continue
    try:
        sc.pl.dotplot(log_adata, gd, groupby=cluster_col, standard_scale="var", show=False)
        savefig(FIG_DOT_DIR / f"{priority}_by_{san(cluster_col)}.png")
    except Exception as e:
        warn(f"dotplot 失敗 ({priority}): {e}")
        plt.close("all")
    try:
        sc.pl.tracksplot(log_adata, gd, groupby=cluster_col, show=False)
        savefig(FIG_TRACKS_DIR / f"{priority}_by_{san(cluster_col)}.png")
    except Exception as e:
        warn(f"tracksplot 失敗 ({priority}): {e}")
        plt.close("all")
print(f"[viz] figures -> {FIG_DIR}")


# %% [markdown]
# ## 4. cluster summary（全 cluster の細胞数・metadata 内訳）

# %%
obs = log_adata.obs
cl = obs[cluster_col].astype(str)
cl.value_counts().rename_axis(cluster_col).to_frame("n_cells").to_csv(
    CLUSTER_SUMMARY_DIR / "cluster_sizes.csv")
for meta in [condition_col, "source_accession", "dataset_id", STATE_COL]:
    if meta and meta in obs.columns:
        ct = pd.crosstab(cl, obs[meta].astype(str))
        ct.to_csv(CLUSTER_SUMMARY_DIR / f"cluster_by_{san(meta)}_counts.csv")
        if meta == condition_col:
            frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
            frac.to_csv(CLUSTER_SUMMARY_DIR / f"cluster_by_{san(meta)}_fraction.csv")
print(f"[cluster_summary] -> {CLUSTER_SUMMARY_DIR}")


# %% [markdown]
# ## 5. DEG（探索的）
# - 各 cluster vs rest の marker（全 cluster）
# - Condition（ALS など）DEG（探索的）
# これらは探索的 cluster/condition marker であり、厳密な condition DEG ではない（pseudo-bulk を優先）。

# %%
# 5-1. 各 cluster vs rest
for method, fname in [("wilcoxon", "markers_by_cluster_wilcoxon.csv"),
                      ("t-test_overestim_var", "markers_by_cluster_ttest_overestim_var.csv")]:
    try:
        key = f"rgg_cluster_{san(method)}"
        sc.tl.rank_genes_groups(log_adata, groupby=cluster_col, method=method,
                                use_raw=False, key_added=key)
        df = sc.get.rank_genes_groups_df(log_adata, group=None, key=key)
        df.to_csv(DEG_DIR / fname, index=False)
        print(f"  saved: {fname}")
    except Exception as e:
        warn(f"cluster marker DEG 失敗 ({method}): {e}")

# 5-2. Condition DEG（探索的）
if condition_col is None:
    warn("Condition 系の列が無いため Condition DEG をスキップ")
else:
    for c in [condition_col, "source_accession", "dataset_id"]:
        if c in obs.columns:
            obs[c].astype(str).value_counts().to_csv(
                DEG_DIR / f"value_counts_{san(c)}.csv", header=["n_cells"])
    cond = log_adata.obs[condition_col].astype(str)
    if cond.nunique() >= 2:
        try:
            log_adata.obs["_cond"] = cond.astype("category")
            sc.tl.rank_genes_groups(log_adata, groupby="_cond", method="wilcoxon",
                                    use_raw=False, key_added="rgg_condition")
            sc.get.rank_genes_groups_df(log_adata, group=None, key="rgg_condition").to_csv(
                DEG_DIR / "condition_deg_wilcoxon.csv", index=False)
            print("  saved: condition_deg_wilcoxon.csv")
        except Exception as e:
            warn(f"Condition DEG 失敗: {e}")
    # ALS vs 非ALS
    als_values = detect_als_values(log_adata.obs[condition_col])
    print(f"  ALS 候補 Condition 値: {als_values}")
    if als_values:
        is_als = log_adata.obs[condition_col].astype(str).isin(als_values).values
        grp = np.where(is_als, "ALS", "nonALS")
        log_adata.obs["_als_group"] = pd.Categorical(grp, categories=["ALS", "nonALS"])
        n_als, n_non = int(is_als.sum()), int((~is_als).sum())
        (DEG_DIR / "als_group_counts.csv").write_text(
            f"ALS,{n_als}\nnonALS,{n_non}\n", encoding="utf-8")
        if n_als >= 20 and n_non >= 20:
            try:
                sc.tl.rank_genes_groups(log_adata, groupby="_als_group", groups=["ALS"],
                                        reference="nonALS", method="wilcoxon",
                                        use_raw=False, key_added="rgg_als")
                sc.get.rank_genes_groups_df(log_adata, group="ALS", key="rgg_als").to_csv(
                    DEG_DIR / "als_vs_rest_wilcoxon.csv", index=False)
                print("  saved: als_vs_rest_wilcoxon.csv")
            except Exception as e:
                warn(f"ALS vs nonALS DEG 失敗: {e}")
        else:
            warn(f"ALS/nonALS 群サイズ不足 (ALS={n_als}, nonALS={n_non})。ALS DEG をスキップ。")


# %% [markdown]
# ## 6. pseudo-bulk（raw_count_like のみ、original-scale）
# 厳密な condition DEG はこちらを優先。sample × cluster と sample × cluster × Condition を集約し、
# edgeR 入力（design ~ cluster + Condition）と、各 cluster 内 ALS vs 非ALS の logCPM Welch を出力。

# %%
pb_info = []
if STATE_COL not in joined.obs.columns:
    warn(f"{STATE_COL} 列が無いため pseudo-bulk をスキップ")
else:
    raw_mask = (joined.obs[STATE_COL].astype(str).values == "raw_count_like")
    n_raw = int(raw_mask.sum())
    pb_info.append(f"raw_count_like cells: {n_raw}")
    print(f"[pseudobulk] raw_count_like 細胞数 = {n_raw}")
    if n_raw == 0:
        warn("raw_count_like 細胞が 0。pseudo-bulk をスキップ。")
    else:
        raw_sub = joined[raw_mask].copy()
        raw_sub.obs[cluster_col] = raw_sub.obs[cluster_col].astype(str)

        # sample 列
        sample_col = first_present(raw_sub.obs.columns, SAMPLE_COL_CANDIDATES)
        if sample_col is None:
            warn("sample_id 候補列が無いため source_accession+dataset_id+Condition で暫定 sample_id を作成"
                 "（統計的に弱いので解釈注意）。")
            parts = [raw_sub.obs[c].astype(str) for c in
                     ["source_accession", "dataset_id", condition_col] if c and c in raw_sub.obs.columns]
            prov = parts[0] if parts else pd.Series(["s0"] * raw_sub.n_obs, index=raw_sub.obs_names)
            for p in parts[1:]:
                prov = prov.str.cat(p, sep="|")
            raw_sub.obs["sample_id_provisional"] = prov.values
            sample_col = "sample_id_provisional"
            STATE["sample_col_provisional"] = True
        STATE["sample_col"] = sample_col
        pb_info.append(f"sample_col: {sample_col} (provisional={STATE['sample_col_provisional']})")
        for c in ["source_accession", "dataset_id", condition_col]:
            if c and c in raw_sub.obs.columns:
                pb_info.append(f"[{c}] {dict(raw_sub.obs[c].astype(str).value_counts())}")

        FACTOR_COLS = [c for c in [condition_col, "source_accession", "dataset_id", STATE_COL] if c]

        # 6-1. sample × cluster
        comp_cluster = (raw_sub.obs[sample_col].astype(str) + "||cl=" + raw_sub.obs[cluster_col].astype(str))
        counts_cluster, _ = pseudobulk_sum(raw_sub, comp_cluster)
        cl_map = dict(zip(comp_cluster.values, raw_sub.obs[cluster_col].astype(str).values))
        sid_map = dict(zip(comp_cluster.values, raw_sub.obs[sample_col].astype(str).values))
        meta_cluster = pseudobulk_metadata(raw_sub, comp_cluster, FACTOR_COLS,
                                           extra_maps={"cluster_label": cl_map, "sample_id": sid_map})
        counts_cluster.round().astype(int).to_csv(PSEUDOBULK_DIR / "pseudobulk_counts_cluster.tsv", sep="\t")
        meta_cluster.to_csv(PSEUDOBULK_DIR / "pseudobulk_metadata_cluster.tsv", sep="\t")
        print(f"[pseudobulk] sample x cluster: {counts_cluster.shape}")

        # 6-2. sample × cluster × Condition（edgeR 用）
        if condition_col and condition_col in raw_sub.obs.columns:
            comp_cc = (raw_sub.obs[sample_col].astype(str)
                       + "||cl=" + raw_sub.obs[cluster_col].astype(str)
                       + "||cond=" + raw_sub.obs[condition_col].astype(str))
            counts_cc, _ = pseudobulk_sum(raw_sub, comp_cc)
            cl_map2 = dict(zip(comp_cc.values, raw_sub.obs[cluster_col].astype(str).values))
            cond_map2 = dict(zip(comp_cc.values, raw_sub.obs[condition_col].astype(str).values))
            sid_map2 = dict(zip(comp_cc.values, raw_sub.obs[sample_col].astype(str).values))
            meta_cc = pseudobulk_metadata(raw_sub, comp_cc, FACTOR_COLS,
                                          extra_maps={"cluster_label": cl_map2,
                                                      "condition_label": cond_map2,
                                                      "sample_id": sid_map2})
            counts_cc_int = counts_cc.round().astype(int)
            counts_cc_int.to_csv(PSEUDOBULK_DIR / "pseudobulk_counts_cluster_condition.tsv", sep="\t")
            meta_cc.to_csv(PSEUDOBULK_DIR / "pseudobulk_metadata_cluster_condition.tsv", sep="\t")
            # logCPM
            lib = counts_cc.sum(axis=0).replace(0, np.nan)
            logcpm = np.log2(counts_cc.divide(lib, axis=1) * 1e6 + 1.0)
            logcpm.to_csv(PSEUDOBULK_DIR / "pseudobulk_logCPM_cluster_condition.tsv", sep="\t")
            print(f"[pseudobulk] sample x cluster x condition: {counts_cc.shape}")

            # edgeR 入力（design ~ cluster + Condition）
            counts_cc_int.to_csv(PSEUDOBULK_DIR / "edgeR_counts.tsv", sep="\t")
            meta_cc.to_csv(PSEUDOBULK_DIR / "edgeR_metadata.tsv", sep="\t")
            design_txt = [
                "# edgeR design info (07b: pseudo-bulk by sample x cluster x Condition)",
                "# counts file   : edgeR_counts.tsv (genes x pseudobulk-samples, raw count sum)",
                "# metadata file : edgeR_metadata.tsv (rows = pseudobulk-samples)",
                f"# cluster col   : 'cluster_label'  (levels = {cluster_col} の各 cluster)",
                f"# condition col : 'condition_label'  (= {condition_col})",
                f"# sample col    : 'sample_id'  (provisional={STATE['sample_col_provisional']})",
                "#",
                "# 推奨 edgeR フロー (R) — 例: cluster ごとに Condition 間 DE、あるいは全体で ~cluster+Condition:",
                "#   library(edgeR)",
                "#   counts <- read.delim('edgeR_counts.tsv', row.names=1, check.names=FALSE)",
                "#   meta   <- read.delim('edgeR_metadata.tsv', row.names=1, check.names=FALSE)",
                "#   meta   <- meta[colnames(counts), ]",
                "#   group  <- factor(paste(meta$cluster_label, meta$condition_label, sep='.'))",
                "#   y <- DGEList(counts=counts, group=group)",
                "#   keep <- filterByExpr(y, group=group); y <- y[keep,, keep.lib.sizes=FALSE]",
                "#   y <- calcNormFactors(y)",
                "#   design <- model.matrix(~ 0 + cluster_label + condition_label, data=meta)",
                "#   y <- estimateDisp(y, design); fit <- glmQLFit(y, design)",
                "#   # contrast で目的の cluster / condition 比較を指定して glmQLFTest",
                "#",
            ]
            if STATE["sample_col_provisional"]:
                design_txt.append(
                    "# [warn] sample_id は暫定。biological replicate として弱いので解釈注意。")
            (PSEUDOBULK_DIR / "edgeR_design_info.txt").write_text("\n".join(design_txt), encoding="utf-8")
            print("  saved: edgeR_counts.tsv / edgeR_metadata.tsv / edgeR_design_info.txt")

            # 6-3. 各 cluster 内 ALS vs 非ALS の logCPM Welch（探索）
            als_values = detect_als_values(raw_sub.obs[condition_col])
            if als_values:
                meta_idx = meta_cc.copy()
                meta_idx["is_als"] = meta_idx["condition_label"].astype(str).isin(als_values)
                deg_rows = []
                for cluster_id in sorted(meta_idx["cluster_label"].astype(str).unique(),
                                         key=lambda x: (len(x), x)):
                    sub_meta = meta_idx[meta_idx["cluster_label"].astype(str) == cluster_id]
                    als_cols = list(sub_meta.index[sub_meta["is_als"].values])
                    non_cols = list(sub_meta.index[~sub_meta["is_als"].values])
                    if len(als_cols) < MIN_PB_SAMPLES or len(non_cols) < MIN_PB_SAMPLES:
                        continue
                    A = logcpm[als_cols].to_numpy()
                    B = logcpm[non_cols].to_numpy()
                    t_stat, p_val = stats.ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
                    t_stat = np.asarray(t_stat, dtype=float)
                    p_val = np.asarray(p_val, dtype=float)
                    log2fc = np.nanmean(A, axis=1) - np.nanmean(B, axis=1)
                    fdr = bh_fdr(p_val)
                    sub_df = pd.DataFrame({
                        "cluster": cluster_id, "gene": counts_cc.index,
                        "log2FC_ALS_vs_nonALS": log2fc, "t_stat": t_stat,
                        "p_value": p_val, "FDR_BH": fdr,
                        "n_samples_ALS": len(als_cols), "n_samples_nonALS": len(non_cols),
                    })
                    deg_rows.append(sub_df)
                if deg_rows:
                    out = pd.concat(deg_rows, ignore_index=True).sort_values(
                        ["cluster", "p_value"], na_position="last")
                    out.to_csv(PSEUDOBULK_DIR / "pseudobulk_DEG_per_cluster_ALS_vs_rest_logCPM_welch.csv",
                               index=False)
                    print("  saved: pseudobulk_DEG_per_cluster_ALS_vs_rest_logCPM_welch.csv "
                          f"({out['cluster'].nunique()} clusters)")
                else:
                    warn("各 cluster で ALS/非ALS の pseudo-bulk sample 数が不足。per-cluster Welch をスキップ。")
        else:
            warn("Condition 列が無いため sample x cluster x Condition / edgeR をスキップ。")

(PSEUDOBULK_DIR / "pseudobulk_inspection.txt").write_text("\n".join(pb_info), encoding="utf-8")


# %% [markdown]
# ## 7. summary

# %%
def list_files(d: Path):
    return sorted(str(p.relative_to(OUT_DIR)) for p in d.rglob("*") if p.is_file()) if d.exists() else []


md = []
md.append("# 07b inner full-gene 解析（microglia_leiden_r05 全cluster軸）— サマリー\n")
md.append("## 入力")
md.append(f"- full inner: `{FULL_INNER_PATH}`  shape={tuple(full.shape)}")
md.append(f"- microglia : `{HVG_RESULT_PATH}`  shape={tuple(hvg.shape)}\n")
md.append("## 結合")
md.append(f"- join key: `{STATE['join_key']}`  / 共通細胞数: {STATE['n_common_cells']}")
md.append(f"- joined（.X=04d original-scale）: `{JOINED_PATH.name}`  shape={tuple(joined.shape)}")
md.append(f"- logexpr（可視化・探索用）: `{LOGEXPR_JOINED_PATH.name}`\n")
md.append("## グループ列 / basis")
md.append(f"- cluster_col: `{STATE['cluster_col']}` (n_clusters={STATE['n_clusters']})")
md.append(f"- umap_basis: `{STATE['umap_basis']}`")
md.append(f"- condition_col: `{STATE['condition_col']}`")
md.append(f"- pseudobulk sample_col: `{STATE['sample_col']}` (provisional={STATE['sample_col_provisional']})\n")
for label, d in [("figures", FIG_DIR), ("cluster_summary", CLUSTER_SUMMARY_DIR),
                 ("deg", DEG_DIR), ("pseudobulk", PSEUDOBULK_DIR)]:
    files = list_files(d)
    md.append(f"## {label}（{len(files)} files）")
    for f in files[:60]:
        md.append(f"- `{f}`")
    md.append("")
md.append("## 注意点")
md.append("- joined `.X` は 04d original-scale（混在）。可視化・探索 DEG は logexpr copy を使用。")
md.append("- cluster marker / Condition DEG（scanpy rank_genes_groups）は **探索的**。")
md.append("- 厳密な condition DEG は raw_count_like の pseudo-bulk（edgeR 入力 / logCPM Welch）を優先。")
if STATE["sample_col_provisional"]:
    md.append("- pseudo-bulk の sample_id は暫定。biological replicate として弱いので解釈注意。")
(OUT_DIR / "README_analysis_summary.md").write_text("\n".join(md), encoding="utf-8")
print("saved:", OUT_DIR / "README_analysis_summary.md")


# %%
print("\n" + "=" * 70)
print("07b 完了")
print("=" * 70)
print("出力先:", OUT_DIR)
print("\n実行コマンド:")
print("    python v2/notebooks/python/07b_restore_inner_genes_microglia_leiden_analysis.py")
