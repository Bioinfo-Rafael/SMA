# SMA v2 Python notebooks

このディレクトリには、SMA / ALS 関連の scRNA-seq / snRNA-seq データを読み込み、
Condition の名寄せ、前処理状態の確認、preprocessing state 別 QC、original-scale merge、
merged h5ad の検収・探索的可視化を行うための Python notebook が入っている。

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
- 探索用 AnnData の保存は任意（保存する場合も original-scale とは別ファイル名にする）。

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

### Deprecated old merge notebook

以前存在した `04_merge_curated_h5ad.ipynb` は削除済み。

この旧 04 は、QC 前の curated h5ad を直接 merge する旧フロー用だった。
現在の正式フローでは、preprocessing state 別 QC 後に `04d_merge_qc_original_scale.ipynb` で merge する。
