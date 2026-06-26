# SMA v2 Python notebooks

このディレクトリには、SMA / ALS 関連の scRNA-seq / snRNA-seq データを読み込み、
Condition の名寄せ、前処理状態の確認、preprocessing state 別 QC、original-scale merge、
merged h5ad の検収・探索的可視化、microglia サブクラスタリング、および
full inner 遺伝子での cluster 解析（marker 可視化・DEG・pseudo-bulk）までを行う
Python notebook が入っている。

各ノートブックは番号順に依存している。基本的には下記の実行順序どおりに上から実行する。

## Recommended execution order

1. `01*.ipynb`
2. `02_curate_each_gse_and_save_h5ad.ipynb`
3. `03_inspect_preprocessing_state.ipynb`
4. `04a_qc_raw_count_like.ipynb`
5. `04b_qc_cpm_tpm_like.ipynb`
6. `04c_qc_log_normalized_like.ipynb`
7. `04d_merge_qc_original_scale.ipynb`
8. `05_check_merged_h5ad.ipynb`
9. `06_microglia_subclustering_annotation.ipynb`
10. `07_restore_inner_genes_and_cluster7_analysis.ipynb`（`.py` 版も同梱）
11. `08/08_classical_full_inner_and_microglia_reclustering.py`（`--pass 1` → 手動annotation → `--pass 2`）

## Notebooks

### 01*.ipynb

Raw / downloaded files を読み込み、各 GSE ごとに interim h5ad を作成する。

Output:
- `v2/data/interim_h5ad/*.h5ad`

### 02_curate_each_gse_and_save_h5ad.ipynb

各 GSE の `obs` を確認しながら、Condition のみを手入力辞書で名寄せする。

この段階では、cell type, genotype, treatment, disease_status などの詳細な標準化は行わない。

curated h5ad の `obs` は最小限にし、元の `obs` 全カラムは sidecar CSV に保存する。

Output:
- `v2/data/curated_h5ad/*.h5ad`
- `v2/data/curated_h5ad/original_obs_metadata_by_cell.csv`
- `v2/data/curated_h5ad/curation_summary.csv`

### 03_inspect_preprocessing_state.ipynb

各 curated h5ad の `.X` がどの前処理状態かを診断する。

推定する preprocessing state:
- `raw_count_like`
- `cpm_tpm_like`
- `log_normalized_like`
- `scaled_zscore_like`
- `unknown`

Output:
- `v2/results/preprocessing_state/preprocessing_state_summary.csv`

### 04a_qc_raw_count_like.ipynb

`raw_count_like` と判定されたデータに対して QC を行う。

QC 条件:
- `n_genes_by_counts >= 200`
- `n_genes_by_counts <= 7000`
- `total_counts > 100`
- `pct_counts_mt < 20`
- `sc.pp.filter_genes(min_cells=3)`

重要:
- `.X` の値は変更しない。
- raw count は raw count のまま保存する。
- stage1 filtered h5ad を保存する。

Output:
- `v2/data/qc_h5ad/raw_count_like/*.stage1_flagged.h5ad`
- `v2/data/qc_h5ad/raw_count_like/*.stage1_filtered.h5ad`
- `v2/results/qc_original_scale_pipeline/04a_raw_count_like/`

### 04b_qc_cpm_tpm_like.ipynb

`cpm_tpm_like` と判定されたデータに対して QC を行う。

QC 条件:
- `n_genes_by_counts >= 200`
- `n_genes_by_counts <= 7000`
- `sc.pp.filter_genes(min_cells=3)`

重要:
- CPM/TPM-like では `total_counts` は raw UMI 数ではない。
- `pct_counts_mt` も raw UMI count 由来の mitochondrial fraction ではない。
- そのため、stage1 QC では `total_counts` や `pct_counts_mt` を除外条件に使わない。
- `.X` の値は CPM/TPM-like のまま保存する。

Output:
- `v2/data/qc_h5ad/cpm_tpm_like/*.stage1_flagged.h5ad`
- `v2/data/qc_h5ad/cpm_tpm_like/*.stage1_filtered.h5ad`
- `v2/results/qc_original_scale_pipeline/04b_cpm_tpm_like/`

### 04c_qc_log_normalized_like.ipynb

`log_normalized_like` と判定されたデータに対して QC を行う。

QC 条件:
- `n_genes_by_counts >= 200`
- `n_genes_by_counts <= 7000`
- `sc.pp.filter_genes(min_cells=3)`

重要:
- log normalized-like では `total_counts` は raw UMI 数ではない。
- `pct_counts_mt` も raw UMI count 由来の mitochondrial fraction ではない。
- そのため、stage1 QC では `total_counts` や `pct_counts_mt` を除外条件に使わない。
- `.X` の値は log normalized のまま保存する。

Output:
- `v2/data/qc_h5ad/log_normalized_like/*.stage1_flagged.h5ad`
- `v2/data/qc_h5ad/log_normalized_like/*.stage1_filtered.h5ad`
- `v2/results/qc_original_scale_pipeline/04c_log_normalized_like/`

### 04d_merge_qc_original_scale.ipynb

04a〜04c で作成した `*.stage1_filtered.h5ad` を読み込み、original-scale のまま outer / inner merge を作る。
その後 merge 後の QC 確認指標を計算・可視化し、手動で追加 QC 条件（integrated QC）を決めてから、
`qc_pass_final == True` の細胞だけに絞った最終 h5ad を保存する。

重要:
- この merge では `normalize_total`, `log1p`, `scale`, `regress_out` を行わない。
- `.X` の値は各 dataset の元 scale を保持する（raw count-like / CPM/TPM-like / log normalized-like が混在）。
- したがって、merged h5ad 全体に対して `total_counts` や `matrix_sum` を同じ意味で解釈してはいけない。
- `matrix_sum` は混在スケールでは意味がないため計算しない。
- 確認指標 `n_nonzero_genes` / `n_nonzero_mt_genes` は現在の `.X` の検出遺伝子数。
- `pct_mt_for_qc` は `qc_preprocessing_state` 別に計算する（log_normalized_like は `np.expm1` で戻してから割合を計算）。
- 追加 QC は `INTEGRATED_QC`（初期値はすべて `None`＝無効）で手動指定し、
  `qc_pass_final = qc_pass_stage1 & qc_pass_integrated`。
- 末尾に、探索用に前処理を揃える補助関数の例をコメントアウトで置く（標準では実行しない）。

Output:
- `v2/data/merged_h5ad/merged_qc_original_scale_outer.h5ad`  ← `qc_pass_final==True` の最終版
- `v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad`  ← 同上（共通遺伝子）
- `v2/data/merged_h5ad/merged_qc_original_scale_{outer,inner}.flagged.h5ad`（QC フラグ付き・未 filter、参照用・任意）
- `v2/results/qc_original_scale_pipeline/04d_merged_original_scale/`（QC 指標 CSV・summary・plot）

### 05_check_merged_h5ad.ipynb

04d で作成した merged h5ad の **検収（acceptance check）＋ 軽い探索的可視化**を行う。本解析ではない。
クラスタリングと自動アノテーションは sanity check / provisional annotation として扱う。

主な内容:
- 基本チェック（shape・obs/var の一意性・outer/inner の細胞一致）。
- `cell_uid` をキーに sidecar から known annotation を戻し、`sex` / `cell_type` を統合（元列も保持）。
- **outer / inner gene set の妥当性確認**（遺伝子存在は outer の 0 埋めではなく各
  `*.stage1_filtered.h5ad` の `var_names` で判定）、HVG 重なり・マーカー保持・MT/ribo 保持・
  gene symbol 診断・inner 採否の判断材料。
- `make_logexpr_copy_from_original_scale()` で **探索用 log-expression copy** を作成
  （original-scale は上書きしない。copy のみ）。
- PCA / UMAP（before / after Harmony）、Leiden clustering（`leidenalg` 無ければ Louvain、両方不可なら skip）。
- **PanglaoDB 参照マーカー**（decoupler / OmniPath。取得できない細胞種のみ fallback）に基づく
  provisional / exploratory な自動アノテーションと、既知 annotation との対応確認。

重要:
- `merged_qc_original_scale_{outer,inner}.h5ad` を上書きしない。
- original-scale h5ad に直接 `normalize_total` / `log1p` / `scale` をかけて上書きしない（必ず copy に対して実行）。
- Harmony / leidenalg / decoupler が無くてもノートブック全体が落ちないようにしている。
- 自動アノテーションは provisional であり、本番 cell type annotation ではない。

Output:
- `v2/results/check_merged_h5ad/`（`gene_set/` `pca_umap/` `clustering/` `annotation/` `plots/`）
- 探索用 AnnData `inner_logexpr_hvg_pca_umap_harmony_cluster_annotation_check.h5ad`
  （HVG 3000・PCA/UMAP/Harmony/Leiden・provisional annotation 済み。06 が入力として読む）。

### 06_microglia_subclustering_annotation.ipynb

05 で保存した探索用 AnnData（`inner_logexpr_hvg_pca_umap_harmony_cluster_annotation_check.h5ad`、
HVG 3000）を読み込み、**microglia 細胞だけを subset して再クラスタリング・サブタイプ自動注釈**を行う。
これも provisional / exploratory な解析であり、最終的な cell type annotation ではない。

主な内容:
- microglia 細胞を subset（66850 cells）し、microglia 専用に Harmony / 近傍グラフ / UMAP
  （`X_umap_microglia`）を作り直し、`microglia_leiden_r05` で再クラスタリング。
- Homeostatic / DAM_activated / Complement の score（`score_*`）と、cluster ごとの score 行列から
  `microglia_subtype_auto`（サブタイプ自動注釈）を付与。
- marker の dotplot / matrixplot / stacked violin、`microglia_leiden_r05` と `microglia_subtype_auto`
  の UMAP、Condition / dataset_id 別 composition、cluster 別・subtype 別 DEG（wilcoxon）。

重要:
- `.X` は 05 由来の探索用 log-expression（HVG 3000）。本解析用の original-scale 行列ではない。
- 自動サブタイプ注釈は provisional であり、本番アノテーションではない。

Output:
- `v2/results/microglia_subclustering/adata_microglia_subclustered.h5ad`
  （66850 microglia cells × 3000 HVG。`obs` に `microglia_leiden_r05` / `microglia_subtype_auto` /
  `score_*`、`obsm` に `X_umap_microglia` / `X_umap_after_harmony` / `X_pca_harmony` 等）
- `marker_{dotplot,matrixplot,stacked_violin}_by_{microglia_leiden_r05,microglia_subtype_auto}.png`
- `umap_microglia_*.png` / `umap_score_*.png`、`composition_*_by_{Condition,dataset_id}.{csv,png}`
- `deg_by_{cluster,subtype}_wilcoxon.csv` / `deg_by_{cluster,subtype}_*.png`
- `marker_availability.csv` / `cluster_subtype_score_matrix.csv` /
  `microglia_subtype_annotation_by_cluster.csv` / `markers/<category>/`

### 07_restore_inner_genes_and_cluster7_analysis.ipynb

06 の microglia-subclustered AnnData（HVG 3000）に載っている **microglia 再クラスタリング /
UMAP / annotation を、04d の full inner gene AnnData（共通遺伝子 ~8863 個）に戻し**、
full inner 遺伝子で marker 可視化・DEG・pseudo-bulk 解析を行う。`.py` 版も同梱（同一内容）。

主な内容:
- `cell_uid` / `obs_names` の overlap で full inner を microglia 細胞（~66850）に整列し、
  06 の `obs` / `obsm` / `uns` を移植する（列名衝突かつ値が異なる場合は `hvg_` prefix）。
  発現行列 `.X` は **04d original-scale を保持**（HVG / logexpr の `.X` で上書きしない）。
- 可視化用に、`qc_preprocessing_state` 別に `normalize_total(1e4)`+`log1p`（log_normalized_like は
  そのまま）した **log-expression copy** を別途作成する。
- cluster 列は `microglia_leiden_r05`、UMAP basis は `X_umap_microglia` を最優先で自動検出。
- **ターゲット cluster = `microglia_leiden_r05` の leiden ラベル "6"**
  （元指示の「cluster7」がこの leiden "6" に対応。出力ファイル名は `cluster7_*` を維持）。
- DEG（探索的）: cluster7(=leiden6) vs rest、cluster7 vs ALS reference（`rank_genes_groups`）。
- pseudo-bulk: `raw_count_like` 細胞のみで sample×cluster / sample×deg_group の count を集約し、
  logCPM + Welch t-test（BH-FDR）と **edgeR 入力ファイル**を出力。

重要:
- joined h5ad の `.X` は **04d original-scale**（raw / cpm / log 混在）。定量比較には不適。
- logexpr h5ad は **可視化・探索用**（exploratory）。
- scanpy `rank_genes_groups` による DEG は探索的。**厳密な condition DEG は raw_count_like の
  pseudo-bulk を優先**（edgeR 入力を用意）。
- 存在しない列 / gene / cluster には警告を出して落ちないようにしている。出力ディレクトリは自動作成。
- SMA root は `__file__` / cwd / 環境変数 `SMA_ROOT` から自動検出する。

実行（SMA リポジトリのルートから、または `SMA_ROOT` 指定で任意ディレクトリから）:

```bash
python v2/notebooks/python/07_restore_inner_genes_and_cluster7_analysis.py
```

Output（`v2/results/full_inner_with_hvg_annotation_analysis/`）:
- `00_inspection_report.txt` / `00_obs_columns_{full,hvg}.csv` / `00_marker_presence_{full_inner,hvg}.csv`
- `inner_fullgenes_with_hvg_umap_harmony_cluster_annotation.h5ad`（joined。`.X` = 04d original-scale）
- `inner_fullgenes_logexpr_with_hvg_umap_harmony_cluster_annotation.h5ad`（logexpr。可視化・探索用）
- `figures/{umap_feature,dotplot,tracksplot}/`
- `cluster7_summary/` / `deg/` / `pseudobulk/`（`edgeR_counts.tsv` / `edgeR_metadata.tsv` /
  `edgeR_design_info.txt` を含む）
- `README_analysis_summary.md`（実行時に自動生成）

### 08/08_classical_full_inner_and_microglia_reclustering.py

04d の full inner-gene AnnData を入力に、**scVI を使わない古典的 Scanpy workflow**
（logexpr layer → HVG → scale → PCA → Harmony → kNN → UMAP → Leiden、Leiden resolution=1.5）で全細胞クラスタリングを行い、
marker を確認して人手で cell type annotation したのち、microglia-like cluster を抽出して
再クラスタリングする。スクリプトと専用 README は `08/` フォルダにまとめてある
（`08/08_classical_full_inner_and_microglia_reclustering.py` /
`08/08_README_classical_full_inner_and_microglia_reclustering.md`）。

2-pass 設計で、1回目 / 2回目は **`--pass` 引数で明示的に指定**する。

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
# 03_manual_annotation/manual_annotation_template_full_clustering.csv を手動で埋め、
# manual_annotation_filled_full_clustering.csv として保存してから
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

重要:
- 入力 `.X` は original-scale 混在のため、`.X`/`.var` を保持して `logexpr_for_clustering` layer を作り、
  クラスタリング・可視化に使う（HVG は PCA/UMAP/clustering 用のみ。marker は full inner genes）。
- `rank_genes_groups` は探索的 cluster marker であり condition DEG ではない。
- microglia-like の選択は `include_for_microglia_recluster` 列（空なら microglia/DAM/myeloid/macrophage
  の keyword fallback）に基づく手動選択。

Output（`v2/results/08_classical_full_inner_microglia_reclustering/`、サブフォルダも番号付き）:
- `01_reports/`（`08_input_qc_report.txt` / `marker_presence_full_inner.csv`）
- `02_full_clustering/`（`full_inner_classical_clustered.h5ad` / `marker_genes/` / `plots/`）
- `03_manual_annotation/`（template・filled・`full_inner_with_manual_annotation.h5ad` /
  `microglia_like_from_manual_annotation.h5ad` / `microglia_selection_summary.csv` / microglia template）
- `04_microglia_reclustering/`（`microglia_classical_reclustered.h5ad` / `marker_genes/` / `plots/`）

## Important notes

### Original-scale merge

`04d_merge_qc_original_scale.ipynb` が作る merged h5ad は、前処理を統一した解析用行列ではない。

これは、QC 後の細胞・遺伝子を集めた original-scale merge である。

そのため、後から用途に応じて以下のように subset して使う。

```python
raw_only = merged[merged.obs["qc_preprocessing_state"] == "raw_count_like"].copy()
```

全データを探索的に可視化したい場合は、04d の末尾にある補助関数を使って、別コピーとして
log-expression 空間に変換する。

### Sidecar metadata

02 では、元の `adata.obs` 全カラムを以下に保存する。

```text
v2/data/curated_h5ad/original_obs_metadata_by_cell.csv
```

curated h5ad や merged h5ad では `obs` は最小限にしている。
細胞アノテーション、sex、subcluster、元論文由来 metadata などが必要になった場合は、
`cell_uid` をキーに sidecar CSV から戻す。

### Acceptance check & provisional annotation (05)

`05_check_merged_h5ad.ipynb` は merged h5ad の検収と探索的可視化のためのもので、本解析ではない。
PCA / UMAP / Harmony / Leiden や PanglaoDB ベースの自動アノテーションは sanity check / provisional
であり、最終的な cell type annotation・DE・scVI などは含めない。original-scale の merged h5ad は
上書きせず、前処理を揃える場合は必ず copy に対して行う。

### Microglia subclustering & inner-gene cluster analysis (06–07)

`06` は 05 の探索用 AnnData（HVG 3000）から microglia を subset して再クラスタリング・
サブタイプ自動注釈を行う。`07` はその microglia 注釈を 04d の full inner genes（~8863）に戻し、
marker 可視化・DEG・pseudo-bulk を行う。

- 06 の自動サブタイプ注釈・07 の `rank_genes_groups` DEG はいずれも **探索的（provisional）**。
- 07 の joined h5ad の `.X` は 04d original-scale（混在スケール）。可視化・探索には logexpr copy を使う。
- **厳密な condition DEG は `raw_count_like` の pseudo-bulk を優先**する（07 が edgeR 入力を出力）。
- 07 の解析対象「cluster7」は `microglia_leiden_r05` の **leiden ラベル "6"**（出力名は `cluster7_*` を維持）。
- pseudo-bulk の sample_id 列が暫定（`source_accession`+`dataset_id`+`Condition`）の場合は
  biological replicate として弱いため、結果の解釈に注意する（07 がログ・README で警告する）。

### Deprecated old merge notebook

以前存在した `04_merge_curated_h5ad.ipynb` は削除済み。

この旧 04 は、QC 前の curated h5ad を直接 merge する旧フロー用だった。
現在の正式フローでは、preprocessing state 別 QC 後に `04d_merge_qc_original_scale.ipynb` で merge する。
