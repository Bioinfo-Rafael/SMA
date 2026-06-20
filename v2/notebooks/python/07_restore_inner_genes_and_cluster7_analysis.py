# %% [markdown]
# # 07. inner 遺伝子復元 + cluster 7 解析
#
# このスクリプトの目的は、HVG 3000 だけになっている探索用 AnnData (05 の出力) に含まれる
# クラスタリング・UMAP・annotation 情報を、04d で作成した full inner gene AnnData
# (inner 遺伝子 ~8863 個) に戻し、その full inner 遺伝子で marker 可視化・DEG・pseudo-bulk
# 解析を行うことである。
#
# 重要な前提（04d / 05 を読んで確認済み）:
# - full inner `.X` は **original-scale**（raw_count_like / cpm_tpm_like / log_normalized_like が混在）。
# - HVG result `.X` は探索用 log-expression（HVG 3000 サブセット）。**この `.X` で上書きしない**。
# - 両者の var_names は **大文字化** されている可能性が高い（例: "APOE"）。よって marker は
#   mouse 式 Title case を基本にしつつ、**case-insensitive matching** で解決する。
# - 05 の post-Harmony UMAP は `obsm["X_umap_after_harmony"]`（basis="umap_after_harmony"）。
# - cell_uid は obs_names と一致している想定だが、両方を確認したうえで対応付けを行う。
#
# 実行（SMA リポジトリのルートから）:
# ```bash
# python v2/notebooks/python/07_restore_inner_genes_and_cluster7_analysis.py
# ```

# %%
# =====================================================================
# import と基本設定
# =====================================================================
import os
import re
import copy as _copy
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

import matplotlib
matplotlib.use("Agg")  # ヘッドレスで figure を保存するため
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

sc.settings.verbosity = 1
warnings.simplefilter("ignore", category=FutureWarning)

# =====================================================================
# パス設定
#   - cwd 決め打ちをやめる。
#   - __file__, cwd, それぞれの親ディレクトリを上にたどり、
#     SMA root を自動検出する。
#   - SMA_ROOT 環境変数があればそれを最優先する。
# =====================================================================

REQUIRED_INPUT_RELPATHS = [
    "v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad",
    "v2/results/check_merged_h5ad/inner_logexpr_hvg_pca_umap_harmony_cluster_annotation_check.h5ad",
]


def candidate_roots():
    """SMA root 候補を列挙する。"""
    roots = []

    # 1. 環境変数を最優先
    env_root = os.environ.get("SMA_ROOT")
    if env_root:
        roots.append(Path(env_root).expanduser().resolve())

    # 2. script 実行時の __file__ から親をたどる
    try:
        here = Path(__file__).resolve()
        roots.append(here.parent)
        roots.extend(here.parents)
    except NameError:
        pass

    # 3. notebook / interactive 実行時の cwd から親をたどる
    cwd = Path.cwd().resolve()
    roots.append(cwd)
    roots.extend(cwd.parents)

    # 重複除去
    out = []
    seen = set()
    for r in roots:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out


def looks_like_sma_root(root: Path) -> bool:
    """必要な入力 h5ad が存在する場所を SMA root とみなす。"""
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
        "以下のどちらかで対処してください。",
        "1. SMAリポジトリのrootで実行する",
        "2. 環境変数 SMA_ROOT を指定する",
        "",
        "例:",
        "  export SMA_ROOT=/home/suzuki/Learn/SMA",
        "  python v2/notebooks/python/07_restore_inner_genes_and_cluster7_analysis.py",
        "",
        "確認したroot候補:",
    ]
    msg.extend([f"  - {p}" for p in checked[:30]])
    msg.append("")
    msg.append("各候補で必要だった入力ファイル:")
    for rel in REQUIRED_INPUT_RELPATHS:
        msg.append(f"  - {rel}")
    raise FileNotFoundError("\n".join(msg))


PROJECT_ROOT = find_project_root()


def rpath(rel: str) -> Path:
    """'v2/...' のような相対パスを SMA root 基準の絶対パスに解決する。"""
    return (PROJECT_ROOT / rel).resolve()


FULL_INNER_PATH = rpath("v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad")
HVG_RESULT_PATH = rpath(
    "v2/results/check_merged_h5ad/inner_logexpr_hvg_pca_umap_harmony_cluster_annotation_check.h5ad"
)
OUT_DIR = rpath("v2/results/full_inner_with_hvg_annotation_analysis")

JOINED_PATH = OUT_DIR / "inner_fullgenes_with_hvg_umap_harmony_cluster_annotation.h5ad"
LOGEXPR_JOINED_PATH = (
    OUT_DIR / "inner_fullgenes_logexpr_with_hvg_umap_harmony_cluster_annotation.h5ad"
)

# 出力サブディレクトリ
FIG_DIR = OUT_DIR / "figures"
FIG_UMAP_DIR = FIG_DIR / "umap_feature"
FIG_DOT_DIR = FIG_DIR / "dotplot"
FIG_TRACKS_DIR = FIG_DIR / "tracksplot"
CLUSTER7_DIR = OUT_DIR / "cluster7_summary"
DEG_DIR = OUT_DIR / "deg"
PSEUDOBULK_DIR = OUT_DIR / "pseudobulk"

# 入力ファイルが存在することを確認してから出力ディレクトリを作る
missing_inputs = [p for p in [FULL_INNER_PATH, HVG_RESULT_PATH] if not p.exists()]
if missing_inputs:
    raise FileNotFoundError(
        "入力ファイルが見つかりません:\n"
        + "\n".join([f"  - {p}" for p in missing_inputs])
        + "\n\nPROJECT_ROOT の推定が間違っている場合は SMA_ROOT を指定してください。"
    )

for d in (
    OUT_DIR, FIG_DIR, FIG_UMAP_DIR, FIG_DOT_DIR, FIG_TRACKS_DIR,
    CLUSTER7_DIR, DEG_DIR, PSEUDOBULK_DIR,
):
    d.mkdir(parents=True, exist_ok=True)

print("PROJECT_ROOT     :", PROJECT_ROOT)
print("FULL_INNER_PATH  :", FULL_INNER_PATH)
print("HVG_RESULT_PATH  :", HVG_RESULT_PATH)
print("OUT_DIR          :", OUT_DIR)

# %%
# =====================================================================
# marker group 定義（mouse 式 Title case）
# =====================================================================
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

# (priority ラベル, group dict) のリスト。CSV や図のループで使う。
MARKER_PRIORITIES = [
    ("priority1", marker_groups_priority1),
    ("priority2", marker_groups_priority2),
]

# 探索順の候補列
UMAP_BASIS_CANDIDATES = [
    "umap_after_harmony", "X_umap_after_harmony", "umap", "X_umap",
]
GROUPBY_CANDIDATES = [
    "leiden_harmony_r05", "hvg_leiden_harmony_r05",
    "leiden_before_harmony_r05", "hvg_leiden_before_harmony_r05",
    "auto_cell_type_marker", "hvg_auto_cell_type_marker",
]
CLUSTER_COL_CANDIDATES = [
    "leiden_harmony_r05", "hvg_leiden_harmony_r05",
    "leiden_before_harmony_r05", "hvg_leiden_before_harmony_r05",
]
TARGET_CLUSTER = "7"
SAMPLE_COL_CANDIDATES = [
    "sample_id", "sample_label", "gsm_id", "donor_id",
    "animal_id", "mouse_id", "orig.ident", "source_file",
]
ALS_TOKENS = ["ALS", "als", "SOD1", "TDP", "disease", "Disease"]
MIN_GROUP_CELLS = 20  # DEG の各群の最低細胞数

# 後段の summary / README 生成のために結果を貯めておく状態 dict
STATE = {
    "join_key": None,
    "n_common_cells": None,
    "cluster_col": None,
    "groupby_col": None,
    "umap_basis": None,
    "n_cluster7": 0,
    "deg_done": {},
    "pseudobulk_done": False,
    "sample_col": None,
    "sample_col_provisional": False,
}

# %%
# =====================================================================
# ユーティリティ関数
# =====================================================================
def warn(msg: str):
    print(f"[warn] {msg}")


def san(s) -> str:
    """ファイル名用に英数以外を _ に置換する。"""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(s)).strip("_")


def savefig(path):
    """matplotlib の現在の figure を保存する（06 の慣習に合わせる）。"""
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


def resolve_genes(upper_map: dict, genes):
    """genes を case-insensitive に var_names へ解決する。

    返り値: (found 順序付き dict {query: var_name}, missing list)
    """
    found = {}
    missing = []
    for g in genes:
        vn = upper_map.get(str(g).upper())
        if vn is not None:
            found[g] = vn
        else:
            missing.append(g)
    return found, missing


def dedup_keep_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def resolve_umap_basis(adata, candidates):
    """obsm から最初に見つかる UMAP basis を返す（sc.pl.embedding 用の basis 文字列）。

    scanpy は basis="foo" のとき obsm["X_foo"] -> obsm["foo"] の順に探すので、
    候補から "X_" を剥がした base を返せば両方の格納形式に対応できる。
    """
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
    """obs 列の概要（dtype, n_unique, 欠損, 例）を DataFrame にする。"""
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
        rows.append({
            "column": c,
            "dtype": str(s.dtype),
            "n_unique": nuniq,
            "n_missing": int(s.isna().sum()),
            "example_values": example,
        })
    return pd.DataFrame(rows)


# %% [markdown]
# ## 1. 確認パート
#
# いきなり解析を進めず、最初に full inner / HVG result の対応関係を確認する。
# 結果はテキストログ (`00_inspection_report.txt`) と CSV で保存する。

# %%
# --- ファイル存在確認 ---
for p in (FULL_INNER_PATH, HVG_RESULT_PATH):
    if not p.exists():
        raise FileNotFoundError(
            f"入力ファイルが見つかりません: {p}\n"
            "04d / 05 を実行して h5ad を生成してから、このスクリプトを実行してください。"
        )

print("loading full inner ...")
full = sc.read_h5ad(FULL_INNER_PATH)
print("loading HVG result ...")
hvg = sc.read_h5ad(HVG_RESULT_PATH)

# 確認レポートの行を貯める
REPORT_LINES = []


def rep(line: str = ""):
    print(line)
    REPORT_LINES.append(line)


rep("=" * 70)
rep("07 inspection report: full inner vs HVG result")
rep("=" * 70)
rep(f"FULL_INNER_PATH : {FULL_INNER_PATH}")
rep(f"HVG_RESULT_PATH : {HVG_RESULT_PATH}")
rep("")

# %%
# --- shape ---
rep("-" * 70)
rep("[shape]")
rep(f"full inner shape : {full.shape}  (cells x genes)")
rep(f"HVG result shape : {hvg.shape}  (cells x genes)")
rep("")

# --- obs_names unique ---
rep("-" * 70)
rep("[obs_names uniqueness]")
full_obs_unique = bool(pd.Index(full.obs_names).is_unique)
hvg_obs_unique = bool(pd.Index(hvg.obs_names).is_unique)
rep(f"full obs_names unique : {full_obs_unique}")
rep(f"HVG  obs_names unique : {hvg_obs_unique}")
rep("")

# --- cell_uid 列の有無と unique ---
rep("-" * 70)
rep("[cell_uid]")
full_has_uid = "cell_uid" in full.obs.columns
hvg_has_uid = "cell_uid" in hvg.obs.columns
rep(f"full has cell_uid : {full_has_uid}")
rep(f"HVG  has cell_uid : {hvg_has_uid}")
if full_has_uid:
    rep(f"full cell_uid unique : {bool(full.obs['cell_uid'].is_unique)}")
if hvg_has_uid:
    rep(f"HVG  cell_uid unique : {bool(hvg.obs['cell_uid'].is_unique)}")
rep("")

# %%
# --- obs_names / cell_uid の一致・overlap ---
rep("-" * 70)
rep("[matching obs_names / cell_uid]")

full_obs = pd.Index(full.obs_names.astype(str))
hvg_obs = pd.Index(hvg.obs_names.astype(str))

obs_names_exact = (len(full_obs) == len(hvg_obs)) and bool((full_obs == hvg_obs).all())
obs_names_overlap = len(full_obs.intersection(hvg_obs))
rep(f"obs_names 完全一致 : {obs_names_exact}")
rep(f"obs_names overlap  : {obs_names_overlap} cells")

if full_has_uid and hvg_has_uid:
    full_uid = pd.Index(full.obs["cell_uid"].astype(str))
    hvg_uid = pd.Index(hvg.obs["cell_uid"].astype(str))
    uid_exact = (len(full_uid) == len(hvg_uid)) and bool((full_uid == hvg_uid).all())
    uid_overlap = len(full_uid.intersection(hvg_uid))
    rep(f"cell_uid 完全一致  : {uid_exact}")
    rep(f"cell_uid overlap   : {uid_overlap} cells")
else:
    uid_exact = None
    uid_overlap = None
    rep("cell_uid 完全一致  : (両方に cell_uid が無いため判定不可)")
rep("")

# %%
# --- obs columns 一覧 -> CSV ---
rep("-" * 70)
rep("[obs columns]")
full_obs_overview = obs_columns_overview(full)
hvg_obs_overview = obs_columns_overview(hvg)
full_obs_overview.to_csv(OUT_DIR / "00_obs_columns_full.csv", index=False)
hvg_obs_overview.to_csv(OUT_DIR / "00_obs_columns_hvg.csv", index=False)
rep(f"full obs columns ({full.obs.shape[1]}): {list(full.obs.columns)}")
rep(f"HVG  obs columns ({hvg.obs.shape[1]}): {list(hvg.obs.columns)}")
rep("saved: 00_obs_columns_full.csv / 00_obs_columns_hvg.csv")
rep("")

# --- HVG side obsm / uns ---
rep("-" * 70)
rep("[HVG obsm / uns keys]")
rep(f"HVG obsm keys : {list(hvg.obsm.keys())}")
rep(f"HVG uns  keys : {list(hvg.uns.keys())}")
rep(f"full obsm keys: {list(full.obsm.keys())}")
rep(f"full uns  keys: {list(full.uns.keys())}")
rep("")

# %%
# --- HVG 側の候補列の存在確認 ---
rep("-" * 70)
rep("[HVG candidate columns]")
hvg_candidate_cols = [
    "leiden_harmony_r05", "leiden_before_harmony_r05", "leiden_before_harmony_r10",
    "auto_cell_type_marker", "cell_type", "Condition",
    "source_accession", "dataset_id", "qc_preprocessing_state",
]
for c in hvg_candidate_cols:
    present = c in hvg.obs.columns
    extra = ""
    if present:
        try:
            nun = hvg.obs[c].nunique(dropna=True)
            extra = f"  (n_unique={nun})"
        except Exception:
            pass
    rep(f"  {c:32s}: {'present' if present else 'MISSING'}{extra}")
rep("")

# %%
# --- marker presence（full / HVG）-> CSV ---
rep("-" * 70)
rep("[marker presence]")


def marker_presence_table(adata) -> pd.DataFrame:
    upper = build_upper_map(adata)
    rows = []
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            for g in genes:
                vn = upper.get(str(g).upper())
                rows.append({
                    "priority": priority,
                    "marker_group": group,
                    "gene_query": g,
                    "found": vn is not None,
                    "matched_var_name": vn if vn is not None else "",
                })
    return pd.DataFrame(rows)


marker_full = marker_presence_table(full)
marker_hvg = marker_presence_table(hvg)
marker_full.to_csv(OUT_DIR / "00_marker_presence_full_inner.csv", index=False)
marker_hvg.to_csv(OUT_DIR / "00_marker_presence_hvg.csv", index=False)

for label, tbl in [("full inner", marker_full), ("HVG", marker_hvg)]:
    n_total = len(tbl)
    n_found = int(tbl["found"].sum())
    rep(f"{label}: marker {n_found}/{n_total} 個が var_names に存在")
    # group ごとの found 数
    g = tbl.groupby("marker_group")["found"].agg(["sum", "count"])
    for grp, row in g.iterrows():
        rep(f"    {grp:34s}: {int(row['sum'])}/{int(row['count'])}")
rep("saved: 00_marker_presence_full_inner.csv / 00_marker_presence_hvg.csv")
rep("")

# %%
# --- 確認レポートを書き出す ---
(OUT_DIR / "00_inspection_report.txt").write_text("\n".join(REPORT_LINES), encoding="utf-8")
print("saved:", OUT_DIR / "00_inspection_report.txt")

# %% [markdown]
# ## 2. 結合パート
#
# cell_uid または obs_names が対応していることを確認したうえで、full inner AnnData に
# HVG result の情報（obs / obsm / 必要な uns）を移植する。
#
# - `.X` は **full inner 側を維持**（HVG の `.X` で上書きしない）。
# - full inner 側の遺伝子（~8863 個）を保持する。
# - obs 列名が衝突して値が異なる場合は `hvg_` prefix を付ける。
# - obsm は全てコピー、uns は衝突する場合 `hvg_` prefix を付ける。
# - obsm を正しく整列させるため、両者の共通細胞（overlap）にそろえる。

# %%
# --- 対応キーの決定（overlap が大きい方を採用） ---
candidates_for_key = []
if full_has_uid and hvg_has_uid:
    candidates_for_key.append(("cell_uid", uid_overlap if uid_overlap is not None else 0))
candidates_for_key.append(("obs_names", obs_names_overlap))

# overlap が最大のキーを選ぶ（同点なら cell_uid 優先 = リスト先頭優先）
join_key, join_overlap = max(candidates_for_key, key=lambda kv: kv[1])
if join_overlap <= 0:
    raise RuntimeError(
        "full inner と HVG result の間で対応する細胞が見つかりません "
        "(obs_names / cell_uid いずれも overlap=0)。入力ファイルを確認してください。"
    )
STATE["join_key"] = join_key
print(f"[join] 対応キー = {join_key} (overlap={join_overlap} cells)")


def key_index(adata, key) -> pd.Index:
    if key == "obs_names":
        return pd.Index(adata.obs_names.astype(str))
    return pd.Index(adata.obs[key].astype(str))


fk = key_index(full, join_key)
hk = key_index(hvg, join_key)
if not fk.is_unique or not hk.is_unique:
    raise RuntimeError(
        f"対応キー {join_key} が unique でないため安全に結合できません。"
    )

# HVG 側の順序を保ったまま共通細胞を取る
common = hk[hk.isin(set(fk))]
common = pd.Index(dedup_keep_order(list(common)))
STATE["n_common_cells"] = int(len(common))
print(f"[join] 共通細胞数 = {len(common)}  "
      f"(full={full.n_obs} のうち {full.n_obs - len(common)} 細胞は HVG 側に対応無し)")

# 共通細胞に対する位置 index
fpos = pd.Series(np.arange(full.n_obs), index=fk).loc[common].to_numpy()
hpos = pd.Series(np.arange(hvg.n_obs), index=hk).loc[common].to_numpy()

# %%
# --- full inner を共通細胞にそろえて joined を作る（.X = full inner を維持） ---
joined = full[fpos].copy()  # full の .X(8863 genes) を保持
hview = hvg[hpos]           # view（HVG の .X はコピーしない）

# obs 列の移植（衝突かつ値が違えば hvg_ prefix）
n_added, n_prefixed, n_skipped = 0, 0, 0
for col in hvg.obs.columns:
    vals = hview.obs[col].values
    if col not in joined.obs.columns:
        joined.obs[col] = vals
        n_added += 1
    else:
        same = False
        try:
            a = pd.Series(joined.obs[col].values).astype(str).to_numpy()
            b = pd.Series(vals).astype(str).to_numpy()
            same = np.array_equal(a, b)
        except Exception:
            same = False
        if same:
            n_skipped += 1  # 値が同じなので何もしない
        else:
            joined.obs["hvg_" + col] = vals
            n_prefixed += 1
print(f"[join] obs: added={n_added}, hvg_prefixed={n_prefixed}, skipped(identical)={n_skipped}")

# obsm は全てコピー
for k in list(hvg.obsm.keys()):
    joined.obsm[k] = np.asarray(hview.obsm[k]).copy()
print(f"[join] obsm copied: {list(joined.obsm.keys())}")

# uns は衝突する場合 hvg_ prefix
n_uns = 0
for k, v in hvg.uns.items():
    newk = k if k not in joined.uns else ("hvg_" + k)
    try:
        joined.uns[newk] = _copy.deepcopy(v)
        n_uns += 1
    except Exception as e:  # deepcopy できない特殊 object は警告のみ
        warn(f"uns['{k}'] をコピーできませんでした: {e}")
print(f"[join] uns copied: {n_uns} keys")

# 由来を記録
joined.uns["restore_note_07"] = (
    "X is 04d original-scale inner genes; obs/obsm/uns annotation restored from 05 HVG result; "
    f"join_key={join_key}; n_common_cells={len(common)}"
)

# %%
# --- joined を保存 ---
print(f"[join] joined shape = {joined.shape}  (期待: cells x ~8863 genes)")
joined.write_h5ad(JOINED_PATH)
print("saved:", JOINED_PATH)

# %% [markdown]
# ## 3. 可視化用 log-expression copy の作成
#
# JOINED の `.X` は 04d 由来の original-scale で、raw_count_like / cpm_tpm_like /
# log_normalized_like が混在している。可視化のため per-cell に正規化した copy を作る。
#
# - `raw_count_like` / `cpm_tpm_like` -> `normalize_total(target_sum=1e4)` -> `log1p`
# - `log_normalized_like` -> そのまま
# - その他の state -> 警告のみ（そのまま）

# %%
STATE_COL = "qc_preprocessing_state"
NORMALIZE_STATES = {"raw_count_like", "cpm_tpm_like"}
ASIS_STATES = {"log_normalized_like"}


def make_logexpr_copy(adata):
    """qc_preprocessing_state に応じて per-cell 正規化した copy を返す（sparse 維持）。"""
    out = adata.copy()
    if STATE_COL not in out.obs.columns:
        warn(f"{STATE_COL} 列が無いため、正規化せずそのまま log-expr copy とします。")
        return out

    states = pd.Series(out.obs[STATE_COL].astype(str).values, index=out.obs_names)
    present = list(pd.unique(states))
    print(f"[logexpr] qc_preprocessing_state: {dict(states.value_counts())}")

    subs = []
    did_normalize = False
    for st in present:
        mask = (states.values == st)
        sub = out[mask].copy()
        if st in NORMALIZE_STATES:
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)
            did_normalize = True
        elif st in ASIS_STATES:
            pass  # そのまま
        else:
            warn(f"未知の state '{st}' はそのままにします（正規化しません）。")
        subs.append(sub)

    merged = ad.concat(subs, axis=0, join="outer", merge="same")
    # 元の細胞順に並べ直す
    merged = merged[out.obs_names].copy()
    # concat で落ちる var / uns / obsm を復元
    merged.var = out.var.copy()
    merged.uns = _copy.deepcopy(out.uns)
    for k in out.obsm.keys():
        if k not in merged.obsm:
            merged.obsm[k] = np.asarray(out.obsm[k]).copy()
    if did_normalize:
        merged.uns["log1p"] = {"base": None}
    merged.uns["logexpr_note_07"] = (
        "per-cell normalize_total(1e4)+log1p for raw_count_like/cpm_tpm_like; "
        "log_normalized_like kept as-is; exploratory visualization only."
    )
    return merged


log_adata = make_logexpr_copy(joined)
log_adata.write_h5ad(LOGEXPR_JOINED_PATH)
print("saved:", LOGEXPR_JOINED_PATH)

# %% [markdown]
# ## 4. marker 可視化
#
# LOGEXPR_JOINED を使って UMAP feature plot / dotplot / tracksplot を作る。
# 存在しない gene は落とさず、存在する gene だけで描画し、欠落 gene 一覧を CSV に保存する。

# %%
# --- UMAP basis と marker の解決 ---
umap_basis = resolve_umap_basis(log_adata, UMAP_BASIS_CANDIDATES)
STATE["umap_basis"] = umap_basis
if umap_basis is None:
    warn(f"UMAP basis が見つかりません（候補: {UMAP_BASIS_CANDIDATES}）。UMAP feature plot をスキップします。")
else:
    print(f"[viz] UMAP basis = {umap_basis} (obsm keys: {list(log_adata.obsm.keys())})")

upper_log = build_upper_map(log_adata)

# 欠落 gene 一覧を貯める
missing_rows = []
for priority, groups in MARKER_PRIORITIES:
    for group, genes in groups.items():
        found, missing = resolve_genes(upper_log, genes)
        for g in missing:
            missing_rows.append({"priority": priority, "marker_group": group, "gene_query": g})
missing_df = pd.DataFrame(missing_rows)
missing_df.to_csv(FIG_DIR / "marker_missing_in_logexpr.csv", index=False)
print(f"[viz] logexpr で欠落している marker: {len(missing_df)} 個 -> figures/marker_missing_in_logexpr.csv")

# %%
# --- 4-1. UMAP feature plot ---
def umap_feature_plot(genes_resolved, out_path, title=None):
    """found な var_name のリストを multipanel embedding にして保存する。"""
    var_names = dedup_keep_order(list(genes_resolved.values()))
    if not var_names or umap_basis is None:
        return False
    try:
        sc.pl.embedding(
            log_adata, basis=umap_basis, color=var_names, color_map="viridis",
            ncols=4, show=False, frameon=False,
        )
        if title:
            plt.suptitle(title)
        savefig(out_path)
        return True
    except Exception as e:
        warn(f"UMAP feature plot 失敗 ({out_path.name}): {e}")
        plt.close("all")
        return False


if umap_basis is not None:
    # priority ごとの全 gene（dedup）
    for priority, groups in MARKER_PRIORITIES:
        all_found = {}
        for group, genes in groups.items():
            f, _ = resolve_genes(upper_log, genes)
            all_found.update(f)
        umap_feature_plot(all_found, FIG_UMAP_DIR / f"{priority}_all_genes.png",
                          title=f"{priority} markers")
    # marker group ごと
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            f, _ = resolve_genes(upper_log, genes)
            umap_feature_plot(f, FIG_UMAP_DIR / f"{priority}_{san(group)}.png",
                              title=f"{priority}: {group}")
    # 個別 gene plot（thorough 用、サブフォルダ）
    indiv_dir = FIG_UMAP_DIR / "individual_genes"
    indiv_dir.mkdir(parents=True, exist_ok=True)
    done_genes = set()
    for priority, groups in MARKER_PRIORITIES:
        for group, genes in groups.items():
            f, _ = resolve_genes(upper_log, genes)
            for q, vn in f.items():
                if vn in done_genes:
                    continue
                done_genes.add(vn)
                try:
                    sc.pl.embedding(log_adata, basis=umap_basis, color=vn,
                                    color_map="viridis", show=False, frameon=False)
                    savefig(indiv_dir / f"{san(q)}.png")
                except Exception as e:
                    warn(f"個別 UMAP feature plot 失敗 ({q}): {e}")
                    plt.close("all")
    print(f"[viz] UMAP feature plots -> {FIG_UMAP_DIR}")

# %%
# --- groupby の決定（dotplot / tracksplot 用） ---
groupby_col = first_present(log_adata.obs.columns, GROUPBY_CANDIDATES)
STATE["groupby_col"] = groupby_col
if groupby_col is None:
    warn(f"dotplot/tracksplot 用 groupby 列が見つかりません（候補: {GROUPBY_CANDIDATES}）。dotplot/tracksplot をスキップします。")
else:
    print(f"[viz] groupby = {groupby_col}")
    # category 化（順序の安定化）
    log_adata.obs[groupby_col] = log_adata.obs[groupby_col].astype(str).astype("category")


def resolved_group_dict(groups):
    """{group: [found var_names]} を作る（空の group は除外）。"""
    d = {}
    for group, genes in groups.items():
        f, _ = resolve_genes(upper_log, genes)
        vns = dedup_keep_order(list(f.values()))
        if vns:
            d[group] = vns
    return d

# %%
# --- 4-2. Dotplot ---
if groupby_col is not None:
    for priority, groups in MARKER_PRIORITIES:
        gdict = resolved_group_dict(groups)
        if not gdict:
            warn(f"{priority}: dotplot に使える gene がありません。")
            continue
        try:
            sc.pl.dotplot(log_adata, gdict, groupby=groupby_col,
                          standard_scale="var", show=False)
            savefig(FIG_DOT_DIR / f"{priority}_by_{san(groupby_col)}.png")
        except Exception as e:
            warn(f"dotplot 失敗 ({priority}): {e}")
            plt.close("all")
        # marker group ごとの dotplot
        for group, vns in gdict.items():
            try:
                sc.pl.dotplot(log_adata, vns, groupby=groupby_col,
                              standard_scale="var", show=False)
                savefig(FIG_DOT_DIR / f"{priority}_{san(group)}_by_{san(groupby_col)}.png")
            except Exception as e:
                warn(f"dotplot 失敗 ({priority}/{group}): {e}")
                plt.close("all")
    print(f"[viz] dotplots -> {FIG_DOT_DIR}")

# %%
# --- 4-3. Tracksplot ---
if groupby_col is not None:
    for priority, groups in MARKER_PRIORITIES:
        gdict = resolved_group_dict(groups)
        if not gdict:
            continue
        try:
            sc.pl.tracksplot(log_adata, gdict, groupby=groupby_col, show=False)
            savefig(FIG_TRACKS_DIR / f"{priority}_by_{san(groupby_col)}.png")
        except Exception as e:
            warn(f"tracksplot 失敗 ({priority}): {e}")
            plt.close("all")
        for group, vns in gdict.items():
            try:
                sc.pl.tracksplot(log_adata, vns, groupby=groupby_col, show=False)
                savefig(FIG_TRACKS_DIR / f"{priority}_{san(group)}_by_{san(groupby_col)}.png")
            except Exception as e:
                warn(f"tracksplot 失敗 ({priority}/{group}): {e}")
                plt.close("all")
    print(f"[viz] tracksplots -> {FIG_TRACKS_DIR}")

# %% [markdown]
# ## 5. cluster 7 の抽出と確認

# %%
cluster_col = first_present(log_adata.obs.columns, CLUSTER_COL_CANDIDATES)
STATE["cluster_col"] = cluster_col

cluster7_ok = False
if cluster_col is None:
    warn(f"cluster 列が見つかりません（候補: {CLUSTER_COL_CANDIDATES}）。"
         "セクション 5-7 をスキップします。")
else:
    print(f"[cluster7] cluster 列 = {cluster_col}")
    cluster_values = log_adata.obs[cluster_col].astype(str)
    avail = sorted(cluster_values.unique(), key=lambda x: (len(x), x))
    if TARGET_CLUSTER not in set(cluster_values):
        warn(f"cluster '{TARGET_CLUSTER}' が {cluster_col} に存在しません。"
             f" 存在する cluster: {avail}")
        warn("セクション 5-7（cluster 7 依存の解析）をスキップします。")
    else:
        cluster7_ok = True
        is_c7 = (cluster_values.values == TARGET_CLUSTER)
        n_c7 = int(is_c7.sum())
        STATE["n_cluster7"] = n_c7
        print(f"[cluster7] cluster {TARGET_CLUSTER} の細胞数 = {n_c7}")

        # 各種内訳を集計して保存
        summary_lines = [f"cluster column: {cluster_col}",
                         f"target cluster: {TARGET_CLUSTER}",
                         f"n_cells: {n_c7}", ""]
        c7_obs = log_adata.obs.loc[is_c7]
        for col in ["source_accession", "dataset_id", "Condition",
                    "qc_preprocessing_state", "auto_cell_type_marker"]:
            if col in c7_obs.columns:
                vc = c7_obs[col].astype(str).value_counts()
                vc.to_csv(CLUSTER7_DIR / f"cluster7_counts_by_{san(col)}.csv",
                          header=["n_cells"])
                summary_lines.append(f"[{col}]")
                summary_lines.extend(f"  {k}: {v}" for k, v in vc.items())
                summary_lines.append("")
            else:
                warn(f"cluster7 集計: 列 '{col}' が無いためスキップ。")
        (CLUSTER7_DIR / "cluster7_summary.txt").write_text(
            "\n".join(summary_lines), encoding="utf-8")
        print(f"[cluster7] 集計 -> {CLUSTER7_DIR}")

# %% [markdown]
# ## 6. DEG 解析（探索的）
#
# 注意: ここでの scanpy `rank_genes_groups` による DEG は **探索的解析** であり、
# 厳密な condition DEG ではない。厳密な condition DEG は pseudo-bulk（セクション 7）で行う。

# %%
def export_rank_genes(adata, group, key, out_csv):
    """rank_genes_groups の結果 (uns[key]) を tidy DataFrame で保存する。"""
    df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
    df.to_csv(out_csv, index=False)
    return df


# --- 6-1. cluster 7 vs rest ---
if cluster7_ok:
    print("[deg] 6-1. cluster 7 vs rest (探索的)")
    # cluster 列を category 化
    log_adata.obs[cluster_col] = log_adata.obs[cluster_col].astype(str).astype("category")
    for method, fname in [
        ("wilcoxon", "cluster7_vs_rest_wilcoxon.csv"),
        ("t-test_overestim_var", "cluster7_vs_rest_ttest_overestim_var.csv"),
    ]:
        try:
            key = f"rgg_c7_vs_rest_{san(method)}"
            sc.tl.rank_genes_groups(
                log_adata, groupby=cluster_col, groups=[TARGET_CLUSTER],
                reference="rest", method=method, use_raw=False, key_added=key,
            )
            export_rank_genes(log_adata, TARGET_CLUSTER, key, DEG_DIR / fname)
            STATE["deg_done"][fname] = True
            print(f"  saved: {fname}")
        except Exception as e:
            warn(f"cluster7 vs rest DEG 失敗 ({method}): {e}")

# %%
# --- 6-2. cluster 7 vs ALS reference ---
def detect_als_values(series):
    vals = pd.Series(series.astype(str)).dropna()
    uniq = sorted(vals.unique())
    toks = [t.lower() for t in ALS_TOKENS]
    return [u for u in uniq if any(t in u.lower() for t in toks)]


def make_deg_group_series(obs, cluster_col, condition_col, als_values):
    """cluster7 / ALS_reference / other の group ラベルを返す。"""
    grp = pd.Series("other", index=obs.index, dtype=object)
    is_c7 = obs[cluster_col].astype(str).values == TARGET_CLUSTER
    if condition_col in obs.columns and als_values:
        is_als = obs[condition_col].astype(str).isin(als_values).values
        grp[is_als & ~is_c7] = "ALS_reference"
    grp[is_c7] = "cluster7"  # cluster7 を最優先
    return grp


als_values = []
deg2_ok = False
if cluster7_ok:
    print("[deg] 6-2. cluster 7 vs ALS reference (探索的)")
    # まず value_counts を出す
    vc_lines = []
    for col in ["Condition", "source_accession", "dataset_id"]:
        if col in log_adata.obs.columns:
            vc = log_adata.obs[col].astype(str).value_counts()
            vc.to_csv(DEG_DIR / f"value_counts_{san(col)}.csv", header=["n_cells"])
            vc_lines.append(f"[{col}]")
            vc_lines.extend(f"  {k}: {v}" for k, v in vc.items())
            vc_lines.append("")
            print(f"  value_counts({col}) -> {len(vc)} 値")
        else:
            warn(f"列 '{col}' が無いため value_counts をスキップ。")

    cond_col = "Condition" if "Condition" in log_adata.obs.columns else None
    if cond_col is None:
        warn("Condition 列が無いため cluster7 vs ALS DEG をスキップします。")
    else:
        als_values = detect_als_values(log_adata.obs[cond_col])
        print(f"  ALS 候補 Condition 値: {als_values}")
        vc_lines.append(f"ALS candidate Condition values: {als_values}")

        grp = make_deg_group_series(log_adata.obs, cluster_col, cond_col, als_values)
        log_adata.obs["deg_group_cluster7_vs_ALS"] = pd.Categorical(
            grp, categories=["cluster7", "ALS_reference", "other"])

        gcounts = grp.value_counts()
        gcounts.to_csv(DEG_DIR / "cluster7_vs_ALS_group_counts.csv", header=["n_cells"])
        n_c7 = int(gcounts.get("cluster7", 0))
        n_ref = int(gcounts.get("ALS_reference", 0))
        print(f"  group counts: cluster7={n_c7}, ALS_reference={n_ref}, "
              f"other={int(gcounts.get('other', 0))}")
        vc_lines.append(f"group counts: {dict(gcounts)}")
        (DEG_DIR / "cluster7_vs_ALS_value_counts.txt").write_text(
            "\n".join(vc_lines), encoding="utf-8")

        if n_ref < MIN_GROUP_CELLS or n_c7 < MIN_GROUP_CELLS:
            warn(f"群サイズが小さすぎます (cluster7={n_c7}, ALS_reference={n_ref}, "
                 f"最低={MIN_GROUP_CELLS})。cluster7 vs ALS の rank test をスキップします。")
        else:
            deg2_ok = True
            for method, fname in [
                ("wilcoxon", "cluster7_vs_ALS_reference_wilcoxon.csv"),
                ("t-test_overestim_var", "cluster7_vs_ALS_reference_ttest_overestim_var.csv"),
            ]:
                try:
                    key = f"rgg_c7_vs_als_{san(method)}"
                    sc.tl.rank_genes_groups(
                        log_adata, groupby="deg_group_cluster7_vs_ALS",
                        groups=["cluster7"], reference="ALS_reference",
                        method=method, use_raw=False, key_added=key,
                    )
                    export_rank_genes(log_adata, "cluster7", key, DEG_DIR / fname)
                    STATE["deg_done"][fname] = True
                    print(f"  saved: {fname}")
                except Exception as e:
                    warn(f"cluster7 vs ALS DEG 失敗 ({method}): {e}")

# %% [markdown]
# ## 7. pseudo-bulk 解析
#
# pseudo-bulk は **raw_count_like dataset だけ** を対象にする。
# 理由: 04d full inner `.X` は original-scale だが raw/cpm/log が混在しており、
# negative binomial 系の pseudo-bulk には raw count が必要なため。
# ここでは JOINED（original-scale）の raw_count_like 細胞のみを使う。

# %%
pb_ready = False
raw_sub = None
sample_col = None
if cluster7_ok and STATE_COL in joined.obs.columns:
    raw_mask = (joined.obs[STATE_COL].astype(str).values == "raw_count_like")
    n_raw = int(raw_mask.sum())
    print(f"[pseudobulk] raw_count_like 細胞数 = {n_raw}")
    if n_raw == 0:
        warn("raw_count_like 細胞が 0 です。pseudo-bulk をスキップします。")
    else:
        raw_sub = joined[raw_mask].copy()  # original-scale raw counts（sparse 維持）

        # 確認: raw_count_like 内の構成
        pb_info = [f"raw_count_like cells: {n_raw}"]
        for col in ["source_accession", "dataset_id", "Condition"]:
            if col in raw_sub.obs.columns:
                vc = raw_sub.obs[col].astype(str).value_counts()
                pb_info.append(f"[{col}] {dict(vc)}")
        # cluster7 の raw 内訳
        if cluster_col in raw_sub.obs.columns:
            n_c7_raw = int((raw_sub.obs[cluster_col].astype(str).values == TARGET_CLUSTER).sum())
            pb_info.append(f"cluster7 cells in raw_count_like: {n_c7_raw}")

        # sample_id 候補列
        sample_col = first_present(raw_sub.obs.columns, SAMPLE_COL_CANDIDATES)
        if sample_col is None:
            # 暫定 sample_id を作る（弱い）
            warn("sample_id 候補列が無いため、source_accession + dataset_id + Condition で "
                 "暫定 sample_id を作ります。これは統計的に弱いので解釈に注意してください。")
            parts = []
            for col in ["source_accession", "dataset_id", "Condition"]:
                if col in raw_sub.obs.columns:
                    parts.append(raw_sub.obs[col].astype(str))
            if parts:
                prov = parts[0].astype(str)
                for p in parts[1:]:
                    prov = prov.str.cat(p.astype(str), sep="|")
            else:
                prov = pd.Series(["sample0"] * raw_sub.n_obs, index=raw_sub.obs_names)
            raw_sub.obs["sample_id_provisional"] = prov.values
            sample_col = "sample_id_provisional"
            STATE["sample_col_provisional"] = True
        STATE["sample_col"] = sample_col
        pb_info.append(f"sample_id column used: {sample_col} "
                       f"(provisional={STATE['sample_col_provisional']})")
        print(f"[pseudobulk] sample_id 列 = {sample_col} "
              f"(provisional={STATE['sample_col_provisional']})")

        # deg_group を raw_sub にも付与（log_adata と同じロジック）
        cond_col_raw = "Condition" if "Condition" in raw_sub.obs.columns else None
        if cond_col_raw is not None:
            als_vals_raw = detect_als_values(raw_sub.obs[cond_col_raw]) if not als_values else als_values
            raw_sub.obs["deg_group_cluster7_vs_ALS"] = make_deg_group_series(
                raw_sub.obs, cluster_col, cond_col_raw, als_vals_raw).values
            gc = raw_sub.obs["deg_group_cluster7_vs_ALS"].value_counts()
            # sample 数
            n_samp_c7 = raw_sub.obs.loc[
                raw_sub.obs["deg_group_cluster7_vs_ALS"] == "cluster7", sample_col].nunique()
            n_samp_ref = raw_sub.obs.loc[
                raw_sub.obs["deg_group_cluster7_vs_ALS"] == "ALS_reference", sample_col].nunique()
            pb_info.append(f"deg_group counts (raw): {dict(gc)}")
            pb_info.append(f"sample counts: cluster7={n_samp_c7}, ALS_reference={n_samp_ref}")
            print(f"[pseudobulk] sample 数: cluster7={n_samp_c7}, ALS_reference={n_samp_ref}")

        (PSEUDOBULK_DIR / "pseudobulk_inspection.txt").write_text(
            "\n".join(pb_info), encoding="utf-8")
        pb_ready = True
elif cluster7_ok:
    warn(f"{STATE_COL} 列が無いため pseudo-bulk をスキップします。")

# %%
# --- 7-1. pseudo-bulk count matrix 作成（sparse 維持の one-hot 集約） ---
def pseudobulk_sum(adata_raw, label_series):
    """label ごとに raw count を sum する。

    返り値: (counts_df: genes x pb_samples, n_cells: Series[pb_sample]).
    one-hot 行列との積で sparse を保ったまま集約する。
    """
    labels = pd.Series(label_series).astype(str)
    valid = labels.notna().values & (labels.values != "nan")
    A = adata_raw[valid]
    lab = labels.values[valid]
    cats = pd.Categorical(lab)
    codes = cats.codes
    n = A.n_obs
    ncat = len(cats.categories)
    M = sparse.csr_matrix(
        (np.ones(n, dtype=np.float64), (codes, np.arange(n))), shape=(ncat, n))
    X = A.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    summed = M @ X  # (ncat x genes)
    summed = np.asarray(summed.todense())
    counts_df = pd.DataFrame(
        summed.T, index=A.var_names.astype(str), columns=list(cats.categories))
    n_cells = pd.Series(np.asarray(M.sum(axis=1)).ravel(),
                        index=list(cats.categories), name="n_cells")
    return counts_df, n_cells


def pseudobulk_metadata(adata_raw, composite_series, factor_cols, label_col, label_value_map):
    """pb_sample ごとの metadata（n_cells と各 factor の最頻値）を作る。"""
    obs = adata_raw.obs.copy()
    obs["_pb"] = pd.Series(composite_series).astype(str).values
    rows = []
    for pb, sub in obs.groupby("_pb"):
        row = {"pb_sample": pb, "n_cells": int(len(sub))}
        for c in factor_cols:
            if c in sub.columns:
                m = sub[c].astype(str).mode()
                row[c] = m.iat[0] if len(m) else ""
        if label_value_map is not None:
            row[label_col] = label_value_map.get(pb, "")
        rows.append(row)
    return pd.DataFrame(rows).set_index("pb_sample")


FACTOR_COLS = ["Condition", "source_accession", "dataset_id", STATE_COL]
# notebook でセル単体を再実行しても NameError にならないよう初期化
counts_grp = None
meta_grp = None
counts_grp_int = None

if pb_ready:
    # (a) sample_id x cluster_label
    if cluster_col in raw_sub.obs.columns:
        comp_cluster = (raw_sub.obs[sample_col].astype(str)
                        + "||cl=" + raw_sub.obs[cluster_col].astype(str))
        counts_cluster, ncell_cluster = pseudobulk_sum(raw_sub, comp_cluster)
        # pb_sample -> cluster ラベルの対応
        cl_map = dict(zip(comp_cluster.values, raw_sub.obs[cluster_col].astype(str).values))
        sid_map = dict(zip(comp_cluster.values, raw_sub.obs[sample_col].astype(str).values))
        meta_cluster = pseudobulk_metadata(
            raw_sub, comp_cluster, FACTOR_COLS, "cluster_label", cl_map)
        meta_cluster["sample_id"] = [sid_map.get(i, "") for i in meta_cluster.index]
        meta_cluster["cluster_label"] = [cl_map.get(i, "") for i in meta_cluster.index]

        counts_cluster.round().astype(int).to_csv(
            PSEUDOBULK_DIR / "pseudobulk_counts_cluster.tsv", sep="\t")
        meta_cluster.to_csv(PSEUDOBULK_DIR / "pseudobulk_metadata_cluster.tsv", sep="\t")
        print(f"[pseudobulk] sample x cluster: {counts_cluster.shape} -> pseudobulk_counts_cluster.tsv")

    # (b) sample_id x deg_group（cluster7 vs ALS_reference に限定）
    if "deg_group_cluster7_vs_ALS" in raw_sub.obs.columns:
        dg = raw_sub.obs["deg_group_cluster7_vs_ALS"].astype(str)
        keep = dg.isin(["cluster7", "ALS_reference"]).values
        sub2 = raw_sub[keep].copy()
        if sub2.n_obs > 0:
            comp_grp = (sub2.obs[sample_col].astype(str)
                        + "||grp=" + sub2.obs["deg_group_cluster7_vs_ALS"].astype(str))
            counts_grp, ncell_grp = pseudobulk_sum(sub2, comp_grp)
            grp_map = dict(zip(comp_grp.values,
                               sub2.obs["deg_group_cluster7_vs_ALS"].astype(str).values))
            sid_map2 = dict(zip(comp_grp.values, sub2.obs[sample_col].astype(str).values))
            meta_grp = pseudobulk_metadata(sub2, comp_grp, FACTOR_COLS, "group", grp_map)
            meta_grp["sample_id"] = [sid_map2.get(i, "") for i in meta_grp.index]
            meta_grp["group"] = [grp_map.get(i, "") for i in meta_grp.index]

            counts_grp_int = counts_grp.round().astype(int)
            counts_grp_int.to_csv(
                PSEUDOBULK_DIR / "pseudobulk_counts_cluster7_vs_ALS.tsv", sep="\t")
            meta_grp.to_csv(
                PSEUDOBULK_DIR / "pseudobulk_metadata_cluster7_vs_ALS.tsv", sep="\t")
            print(f"[pseudobulk] sample x deg_group: {counts_grp.shape} "
                  "-> pseudobulk_counts_cluster7_vs_ALS.tsv")
        else:
            warn("cluster7 / ALS_reference に該当する raw_count_like 細胞がありません。")
            counts_grp = None
            meta_grp = None
    else:
        counts_grp = None
        meta_grp = None

# %%
# --- 7-2. pseudo-bulk DEG（logCPM + Welch t-test, Python のみ） ---
def bh_fdr(pvals):
    """Benjamini-Hochberg FDR。NaN は NaN のまま返す。"""
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


if pb_ready and counts_grp is not None:
    print("[pseudobulk] 7-2. logCPM + Welch t-test (cluster7 vs ALS_reference)")
    counts = counts_grp.copy()  # genes x pb_samples
    # 各 pb_sample が cluster7 か ALS_reference か
    grp_of_sample = meta_grp["group"].reindex(counts.columns).astype(str)
    c7_cols = list(grp_of_sample.index[grp_of_sample.values == "cluster7"])
    ref_cols = list(grp_of_sample.index[grp_of_sample.values == "ALS_reference"])
    print(f"  pseudobulk samples: cluster7={len(c7_cols)}, ALS_reference={len(ref_cols)}")

    # CPM -> logCPM
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    logcpm = np.log2(cpm + 1.0)
    logcpm.to_csv(PSEUDOBULK_DIR / "pseudobulk_logCPM_cluster7_vs_ALS.tsv", sep="\t")

    if len(c7_cols) >= 2 and len(ref_cols) >= 2:
        A = logcpm[c7_cols].to_numpy()
        B = logcpm[ref_cols].to_numpy()
        t_stat, p_val = stats.ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
        t_stat = np.asarray(t_stat, dtype=float)
        p_val = np.asarray(p_val, dtype=float)
        mean_c7 = np.nanmean(A, axis=1)
        mean_ref = np.nanmean(B, axis=1)
        log2fc = mean_c7 - mean_ref  # logCPM はすでに log2 スケール
        fdr = bh_fdr(p_val)
        res = pd.DataFrame({
            "gene": counts.index,
            "log2FC_cluster7_vs_ALS": log2fc,
            "mean_logCPM_cluster7": mean_c7,
            "mean_logCPM_ALS_reference": mean_ref,
            "t_stat": t_stat,
            "p_value": p_val,
            "FDR_BH": fdr,
            "n_samples_cluster7": len(c7_cols),
            "n_samples_ALS_reference": len(ref_cols),
        }).sort_values("p_value", na_position="last")
        res.to_csv(
            PSEUDOBULK_DIR / "pseudobulk_DEG_cluster7_vs_ALS_logCPM_welch.csv", index=False)
        STATE["pseudobulk_done"] = True
        print("  saved: pseudobulk_DEG_cluster7_vs_ALS_logCPM_welch.csv")
    else:
        warn(f"pseudo-bulk sample 数が不足 (cluster7={len(c7_cols)}, "
             f"ALS_reference={len(ref_cols)})。Welch t-test をスキップします "
             "(logCPM 行列のみ保存)。各群 sample >=2 が必要です。")

    # --- R / edgeR 用入力ファイル ---
    if counts_grp_int is None:
        counts_grp_int = counts_grp.round().astype(int)
    counts_grp_int.to_csv(PSEUDOBULK_DIR / "edgeR_counts.tsv", sep="\t")
    edger_meta = meta_grp.copy()
    edger_meta.to_csv(PSEUDOBULK_DIR / "edgeR_metadata.tsv", sep="\t")
    design_txt = [
        "# edgeR design info (cluster7 vs ALS_reference pseudo-bulk)",
        f"# counts file   : edgeR_counts.tsv (genes x pseudobulk-samples, raw count sum)",
        f"# metadata file : edgeR_metadata.tsv (rows = pseudobulk-samples)",
        f"# group column  : 'group'  (levels: cluster7, ALS_reference)",
        f"# sample column : 'sample_id'  (provisional={STATE['sample_col_provisional']})",
        "#",
        "# 推奨 edgeR フロー (R):",
        "#   library(edgeR)",
        "#   counts <- read.delim('edgeR_counts.tsv', row.names=1, check.names=FALSE)",
        "#   meta   <- read.delim('edgeR_metadata.tsv', row.names=1, check.names=FALSE)",
        "#   meta   <- meta[colnames(counts), ]",
        "#   group  <- factor(meta$group, levels=c('ALS_reference','cluster7'))",
        "#   y <- DGEList(counts=counts, group=group)",
        "#   keep <- filterByExpr(y, group=group); y <- y[keep,, keep.lib.sizes=FALSE]",
        "#   y <- calcNormFactors(y)",
        "#   design <- model.matrix(~group)",
        "#   y <- estimateDisp(y, design)",
        "#   fit <- glmQLFit(y, design); qlf <- glmQLFTest(fit, coef=2)",
        "#   topTags(qlf, n=Inf)",
        "#",
    ]
    if STATE["sample_col_provisional"]:
        design_txt.append(
            "# [warn] sample_id は暫定（source_accession+dataset_id+Condition）。"
            "biological replicate として弱いので結果の解釈に注意。")
    (PSEUDOBULK_DIR / "edgeR_design_info.txt").write_text(
        "\n".join(design_txt), encoding="utf-8")
    print("  saved: edgeR_counts.tsv / edgeR_metadata.tsv / edgeR_design_info.txt")

# %% [markdown]
# ## 8. 出力まとめ（README_analysis_summary.md）

# %%
def list_files(d: Path):
    if not d.exists():
        return []
    return sorted(str(p.relative_to(OUT_DIR)) for p in d.rglob("*") if p.is_file())


md = []
md.append("# 07. inner 遺伝子復元 + cluster 7 解析 — サマリー\n")
md.append("## 入力 h5ad")
md.append(f"- full inner: `{FULL_INNER_PATH}`  shape={tuple(full.shape)}")
md.append(f"- HVG result: `{HVG_RESULT_PATH}`  shape={tuple(hvg.shape)}\n")

md.append("## full inner と HVG result の対応確認")
md.append(f"- 対応キー: `{STATE['join_key']}`")
md.append(f"- obs_names 完全一致: {obs_names_exact}")
md.append(f"- obs_names overlap: {obs_names_overlap}")
if uid_exact is not None:
    md.append(f"- cell_uid 完全一致: {uid_exact} / overlap: {uid_overlap}")
md.append(f"- 共通細胞数（joined の細胞数）: {STATE['n_common_cells']}\n")

md.append("## 作成した h5ad")
md.append(f"- joined（04d original-scale `.X` + 05 annotation）: "
          f"`{JOINED_PATH.relative_to(OUT_DIR)}`  shape={tuple(joined.shape)}")
md.append(f"- logexpr（可視化・探索用、per-cell normalize+log1p）: "
          f"`{LOGEXPR_JOINED_PATH.relative_to(OUT_DIR)}`  shape={tuple(log_adata.shape)}\n")

md.append("## 使用した列 / basis")
md.append(f"- cluster 列: `{STATE['cluster_col']}`")
md.append(f"- dotplot/tracksplot groupby: `{STATE['groupby_col']}`")
md.append(f"- UMAP basis: `{STATE['umap_basis']}`\n")

md.append("## cluster 7")
if cluster7_ok:
    md.append(f"- 細胞数: {STATE['n_cluster7']}")
    for col in ["Condition", "source_accession"]:
        f = CLUSTER7_DIR / f"cluster7_counts_by_{san(col)}.csv"
        if f.exists():
            s = pd.read_csv(f, index_col=0)
            md.append(f"- {col} 内訳:")
            for k, v in s.iloc[:, 0].items():
                md.append(f"    - {k}: {v}")
    md.append("- 詳細: `cluster7_summary/`")
else:
    md.append("- cluster 7 は検出できませんでした（cluster 列が無い、または cluster '7' が存在しない）。")
md.append("")

md.append("## marker gene の存在 / 欠落")
md.append(f"- full inner: {int(marker_full['found'].sum())}/{len(marker_full)} 個存在 "
          "（`00_marker_presence_full_inner.csv`）")
md.append(f"- HVG result: {int(marker_hvg['found'].sum())}/{len(marker_hvg)} 個存在 "
          "（`00_marker_presence_hvg.csv`）")
md.append(f"- logexpr で欠落: {len(missing_df)} 個（`figures/marker_missing_in_logexpr.csv`）\n")

md.append("## 図（UMAP feature / dotplot / tracksplot）")
for label, d in [("UMAP feature", FIG_UMAP_DIR), ("dotplot", FIG_DOT_DIR), ("tracksplot", FIG_TRACKS_DIR)]:
    files = list_files(d)
    md.append(f"- {label}: {len(files)} 枚")
    for f in files[:40]:
        md.append(f"    - `{f}`")
md.append("")

md.append("## DEG 結果ファイル")
deg_files = list_files(DEG_DIR)
if deg_files:
    for f in deg_files:
        md.append(f"- `{f}`")
else:
    md.append("- （なし）")
md.append("")

md.append("## pseudo-bulk 結果ファイル")
pb_files = list_files(PSEUDOBULK_DIR)
if pb_files:
    for f in pb_files:
        md.append(f"- `{f}`")
else:
    md.append("- （なし）")
md.append("")

md.append("## 注意点")
md.append("- joined h5ad の `.X` は **04d original-scale**（raw/cpm/log 混在）。定量比較には不適。")
md.append("- logexpr h5ad は **可視化・探索用**（per-cell normalize_total(1e4)+log1p、log_normalized_like はそのまま）。")
md.append("- scanpy `rank_genes_groups` による DEG（6-1 / 6-2）は **探索的**。cluster7 vs ALS は厳密な condition DEG ではない。")
md.append("- 厳密な condition DEG は **raw_count_like の pseudo-bulk**（7）を優先（edgeR 入力を用意）。")
if STATE["sample_col_provisional"]:
    md.append("- pseudo-bulk の sample_id は **暫定**（source_accession+dataset_id+Condition）。"
              "biological replicate として弱いため解釈に注意。")
elif STATE["sample_col"]:
    md.append(f"- pseudo-bulk の sample_id 列: `{STATE['sample_col']}`")

(OUT_DIR / "README_analysis_summary.md").write_text("\n".join(md), encoding="utf-8")
print("saved:", OUT_DIR / "README_analysis_summary.md")

# %%
# =====================================================================
# 完了メッセージと実行コマンド
# =====================================================================
print("\n" + "=" * 70)
print("07 解析 完了")
print("=" * 70)
print("出力先:", OUT_DIR)
print("\n実行コマンド（SMA リポジトリのルートから）:")
print("    python v2/notebooks/python/07_restore_inner_genes_and_cluster7_analysis.py")
