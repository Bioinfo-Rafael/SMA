# SMA / ALS scRNA-seq preprocessing (v2)

**download スクリプト + 対話的ノートブック** 構成。

`.py` スクリプトは GEO **Supplementary files** の取得・整理（download → extract →
list → overview）まで。それ以降（AnnData 化・obs/var 確認・metadata 整形・前処理段階の
診断・per-GSE h5ad 保存・merge）はすべて **Jupyter ノートブック**で行い、cell type /
cluster / condition 列は人間が判断する（自動判定しない）。

> コード内のコメント・docstring は日本語。関数定義は `src/`、実行単位はノートブック。

各ディレクトリの詳細はそれぞれの README を参照：
[`scripts/`](scripts/README.md) / [`notebooks/`](notebooks/README.md) / [`data/`](data/README.md)。

## 実行手順（この順番）

```bash
# --- 0. 環境（初回のみ）---
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name sma-v2 --display-name "Python (SMA v2)"
# （v1 の venv を流用するなら PYTHON=../v1/.venv/bin/python ./run.sh）

# --- 毎回の起動 ---
cd v2
source .venv/bin/activate          # スクリプト実行用
# または run.sh で:  PYTHON=.venv/bin/python ./run.sh
jupyter lab                        # ノートブックはカーネル「Python (SMA v2)」を選択

# --- 1. ダウンロード/展開（.py スクリプト）---
python scripts/00_validate_manifest.py        # manifest 検証
python scripts/01_download_geo_supplement.py  # data/raw/ へDL（レジューム）
python scripts/02_extract_archives.py         # data/extracted/ へ展開
python scripts/03_list_downloaded_files.py    # 欠落チェック
python scripts/04_overview.py                 # 俯瞰
#   まとめて：  ./run.sh

# --- 2. （GSE295514 のみ）R ノートブックで RDS を中間ファイル化 ---
#   notebooks/R/01_GSE295514_read_rds.ipynb を R kernel で実行
#   -> data/intermediate_from_r/GSE295514/ に counts.mtx 等を出力

# --- 3. AnnData 以降（Jupyter ノートブックを上から順に）---
jupyter lab
#   notebooks/python/01_load_and_inspect_each_gse.ipynb     # 読み込み＋確認
#   notebooks/python/02_curate_each_gse_and_save_h5ad.ipynb # 整形＋名寄せ履歴＋保存
#   notebooks/python/03_inspect_preprocessing_state.ipynb   # 前処理段階の診断
#   notebooks/python/04a_qc_raw_count_like.ipynb            # raw_count_like の QC（細胞/遺伝子選択のみ）
#   notebooks/python/04b_qc_cpm_tpm_like.ipynb              # cpm_tpm_like の QC
#   notebooks/python/04c_qc_log_normalized_like.ipynb       # log_normalized_like の QC
#   notebooks/python/04d_merge_qc_original_scale.ipynb      # original-scale merge ＋ integrated QC
#   notebooks/python/05_check_merged_h5ad.ipynb             # merged の検収＋探索的可視化
```

要約：**scripts 00→04** →（GSE295514 は **R/01**）→ **python 01→02→03→04a→04b→04c→04d→05**。

## どこで何をするか

| 段階 | 場所 |
|---|---|
| manifest 検証 / download / extract / list / overview | `.py`（`scripts/`） |
| ファイル → AnnData 化 ＋ その場で確認 | notebook `python/01` |
| metadata 整形（GSEごと・名寄せ履歴を記録/CSV出力） | notebook `python/02` |
| per-GSE curated h5ad 保存 | notebook `python/02` |
| 前処理段階の診断（統計量・分布で推定。揃えない） | notebook `python/03` |
| preprocessing state 別 QC（細胞/遺伝子選択のみ・`.X` は不変） | notebook `python/04a` `04b` `04c` |
| QC 済み細胞・遺伝子の original-scale merge ＋ integrated QC | notebook `python/04d` |
| merged h5ad の検収＋探索的可視化（gene set 確認・PCA/UMAP/Harmony/Leiden・provisional annotation） | notebook `python/05` |
| RDS を開いて中間ファイル出力 | **R** notebook `R/01_GSE295514_read_rds.ipynb` |

## データ上の注意

* **GSE242942** → scRNA-seq SubSeries **`GSE242939`** のみ使用。bulk の **`GSE242940` は不使用**。
* **GSE167332** → 3 SubSeries を**別々の h5ad** として保存：`GSE167198`（Drop-seq 全脊髄）、
  `GSE167327`（CD45 enriched, inDrop）、`GSE167331`（FACS microglia, SmartSeq2）。
* **GSE167331** は **TPM** → `data_status = processed_TPM`。生カウントと混ぜない。
* **GSE206330** は **SoupX 補正済み** → `data_status = processed_SoupX_corrected`。生カウントと別扱い。
* **GSE295514** は **RDS**：R kernel で開いて中間ファイルを出力 →
  `data_status = RDS_converted_unknown_or_counts`。
* **GSE173524** は生の `GSE173524_umi.tsv.gz` を使用（`*.sctransform.*` は不使用）。
* QC は **preprocessing state 別**（`04a`/`04b`/`04c`）。細胞・遺伝子を選択するだけで `.X` の値は変えない
  （raw count はそのまま、CPM/TPM・log normalized もそのまま保存）。
* merge（`04d`）は `data/merged_h5ad/` に **original-scale** で保存（正規化しない）。`.X` は
  raw_count_like / cpm_tpm_like / log_normalized_like が **混在**するので、merged 全体に同じ意味の
  `total_counts` / `pct_mt` を当てはめない。**outer**（遺伝子の和集合）と **inner**（共通遺伝子）の2種類。

## 構成

```
v2/
├── config/dataset_manifest.yaml      # 真実の情報源（GSE/files/URL/loader_hint/metadata）
├── scripts/                          # download/extract まで（.py）
│   ├── 00_validate_manifest.py
│   ├── 01_download_geo_supplement.py
│   ├── 02_extract_archives.py        # 安全展開（パストラバーサル対策）+ ネスト tar
│   ├── 03_list_downloaded_files.py
│   └── 04_overview.py                # 俯瞰（旧 00_overview ノートブックをスクリプト化）
├── src/                              # ノートブックから呼ぶ関数群
│   ├── geo_download.py               # レジューム付きダウンロード
│   ├── archive_utils.py              # 安全な tar 展開 + find_files
│   ├── manifest_utils.py             # manifest 読込/検証・パス・ロガー
│   ├── io_10x.py                     # 10x .h5 / MTX 三点セット ローダー
│   ├── io_dense.py                   # dense/結合/処理済み/ネスト + R 中間ファイル読込
│   ├── anndata_utils.py              # obs/var スキーマ・obs_names・保存/一括ロード
│   └── notebook_report_utils.py      # summarize/show_* + 前処理診断 + CurationLog
├── notebooks/
│   ├── python/
│   │   ├── 01_load_and_inspect_each_gse.ipynb
│   │   ├── 02_curate_each_gse_and_save_h5ad.ipynb
│   │   ├── 03_inspect_preprocessing_state.ipynb
│   │   ├── 04a_qc_raw_count_like.ipynb
│   │   ├── 04b_qc_cpm_tpm_like.ipynb
│   │   ├── 04c_qc_log_normalized_like.ipynb
│   │   ├── 04d_merge_qc_original_scale.ipynb
│   │   └── 05_check_merged_h5ad.ipynb
│   └── R/
│       └── 01_GSE295514_read_rds.ipynb
├── data/                             # git 管理外。初回実行で作成
│   ├── raw/<acc>/                    # ダウンロードした supplementary
│   ├── extracted/<acc>/              # 展開済み
│   ├── intermediate_from_r/<acc>/    # R ノートブックの出力（counts.mtx 等）
│   ├── interim_h5ad/                 # 生 AnnData（notebook 01 で保存）
│   ├── curated_h5ad/                 # notebook 02 の出力（+ sidecar CSV）
│   ├── qc_h5ad/<state>/              # notebook 04a/04b/04c の出力（stage1_flagged / stage1_filtered）
│   ├── merged_h5ad/                  # notebook 04d の出力（original-scale outer / inner）
│   └── reports/                      # manifest 一覧・ファイル一覧・名寄せ履歴 等
├── results/                          # git 管理外。ノートブックの出力（QC summary・図・検収結果）
│   ├── preprocessing_state/          # notebook 03
│   ├── qc_original_scale_pipeline/   # notebook 04a–04d
│   └── check_merged_h5ad/            # notebook 05
├── requirements.txt
└── run.sh                            # download/extract/list/overview のみ
```

## R ノートブック（GSE295514 RDS）

`notebooks/R/01_GSE295514_read_rds.ipynb` は **R Jupyter kernel** 上で動く。RDS を読み、
class/assays/`meta.data` を確認し、`counts.mtx` / `metadata.csv` / `genes.csv` /
`barcodes.csv` を `data/intermediate_from_r/GSE295514/` に書き出す。Python 側は
`io_dense.read_from_r_intermediate` でこれを AnnData 化する。

R 側で必要になりうるもの（R 内で各自インストール。自動 install はしない）：

```r
install.packages("IRkernel"); IRkernel::installspec()   # Jupyter 用 R kernel
install.packages("Matrix")
install.packages("Seurat")              # SeuratObject も入る
# BiocManager::install("SingleCellExperiment")
# 任意: zellkonverter, Matrix.utils
```

## スコープ

**download / AnnData 化 / 手動確認 / 前処理段階の診断 / 整形 / per-GSE h5ad 保存 /
preprocessing state 別 QC（細胞・遺伝子選択のみ・`.X` は不変）/ original-scale merge ＋ integrated QC /
merged h5ad の検収** まで。

`05` では検収の一環として、探索用 log-expression copy（copy のみ・original-scale は上書きしない）に
対する PCA / UMAP / Harmony / Leiden clustering と、PanglaoDB ベースの **provisional / exploratory** な
自動アノテーションも行う。ただしこれは sanity check であり、最終的な cell type annotation・DE 解析・
scVI・本番の batch correction は **まだ含めない**。
