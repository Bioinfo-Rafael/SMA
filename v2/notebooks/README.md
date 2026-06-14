# notebooks/

AnnData 化以降の**対話的**な処理。スクリプト（`../scripts/`）でダウンロード/展開を
済ませた後、ここを上から順に実行する。cell type / cluster / condition 列は
原則 **自動判定せず**、obs/var を見ながら人間が決める（`05` の PanglaoDB ベース自動
アノテーションは provisional な sanity check で、最終判断は人間が行う）。関数定義は `../src/`。

> 各ノートブック冒頭のセルが `config/dataset_manifest.yaml` を持つ v2 ルートを
> 自動で探して `src/` を import パスに通すので、起動ディレクトリは問わない。

## python/

| ノートブック | 役割 | 入力 → 出力 |
|---|---|---|
| `01_load_and_inspect_each_gse.ipynb` | GSEごとに AnnData 化し、その場で obs/var/値を確認（旧 load+inspect を統合） | `data/extracted`・`data/raw` → `data/interim_h5ad/` |
| `02_curate_each_gse_and_save_h5ad.ipynb` | GSEごとに obs/var を整形。`CurationLog` で名寄せ履歴を記録しCSV出力。元 obs は残す | `data/interim_h5ad/` → `data/curated_h5ad/`、`data/reports/curation_rename_log.csv` |
| `03_inspect_preprocessing_state.ipynb` | 遺伝子統計量・値分布で「どこまで前処理されたか」を推定する**診断のみ**（揃えない・正規化しない） | `data/curated_h5ad/`（無ければ interim） → `results/preprocessing_state/` |
| `04a_qc_raw_count_like.ipynb` | `raw_count_like` の stage1 QC。細胞/遺伝子を選択するだけで `.X` は不変（n_genes・total_counts・pct_mt で除外） | `data/curated_h5ad/`＋`03` の summary → `data/qc_h5ad/raw_count_like/`・`results/` |
| `04b_qc_cpm_tpm_like.ipynb` | `cpm_tpm_like` の stage1 QC（n_genes のみで選択。total_counts/pct_mt は除外条件に使わない）。`.X` 不変 | 同上 → `data/qc_h5ad/cpm_tpm_like/`・`results/` |
| `04c_qc_log_normalized_like.ipynb` | `log_normalized_like` の stage1 QC（n_genes のみで選択）。`.X` 不変 | 同上 → `data/qc_h5ad/log_normalized_like/`・`results/` |
| `04d_merge_qc_original_scale.ipynb` | `*.stage1_filtered.h5ad` を **original-scale のまま** outer/inner merge ＋ integrated QC（正規化しない・`.X` は混在スケール） | `data/qc_h5ad/*/` → `data/merged_h5ad/merged_qc_original_scale_{outer,inner}.h5ad`・`results/` |
| `05_check_merged_h5ad.ipynb` | merged の検収＋探索的可視化（gene set 確認・log-expr copy・PCA/UMAP/Harmony/Leiden・PanglaoDB provisional annotation）。本解析ではない | `data/merged_h5ad/` → `results/check_merged_h5ad/` |

## R/

| ノートブック | 役割 | 出力 |
|---|---|---|
| `01_GSE295514_read_rds.ipynb` | RDS を R kernel で開き class/assays/meta.data を確認し、Python 用中間ファイルを書き出す | `data/intermediate_from_r/GSE295514/`（counts.mtx, metadata.csv, genes.csv, barcodes.csv） |

GSE295514 は **この R ノートを先に実行**してから `python/01` の GSE295514 セルで読む。

## 実行順

`R/01`（GSE295514 のみ）→ `python/01 → 02 → 03 → 04a → 04b → 04c → 04d → 05`。詳細は `../README.md`。
