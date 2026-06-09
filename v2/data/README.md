# data/

パイプラインの入出力データ置き場。**中身（ファイル本体）は git 管理外**
（リポジトリ直下の `.gitignore` が `data/` を無視）。各サブディレクトリは
スクリプト/ノートブックの初回実行時に自動作成される。この README だけは追跡する。

## サブディレクトリ

| ディレクトリ | 中身 | 生成する処理 |
|---|---|---|
| `raw/<acc>/` | GEO からダウンロードした supplementary ファイル（tar/gz/tsv/rds 等） | `scripts/01_download_geo_supplement.py` |
| `extracted/<acc>/` | tar を展開した中身（10x の mtx/tsv、h5、dense txt 等） | `scripts/02_extract_archives.py` |
| `intermediate_from_r/<acc>/` | R ノートブックが書き出す中間ファイル（counts.mtx, metadata.csv, genes.csv, barcodes.csv） | `notebooks/R/01_GSE295514_read_rds.ipynb` |
| `interim_h5ad/` | GSEごとの「生」AnnData（整形前）。`<dataset_id>.h5ad` | `notebooks/python/01` |
| `curated_h5ad/` | obs/var を整形した GSEごとの AnnData。`<dataset_id>.h5ad` | `notebooks/python/02` |
| `merged_h5ad/` | merge 済み AnnData（`merged_raw_or_filtered_count_*` と `merged_all_status_aware_*`） | `notebooks/python/04` |
| `reports/` | 一覧・診断・履歴などのテキスト/CSV | 各スクリプト・ノートブック |

## `reports/` の主なファイル

| ファイル | 内容 | 出力元 |
|---|---|---|
| `manifest_overview.csv` | データセット一覧表 | `scripts/00` |
| `downloaded_files.txt` | DL/展開済みファイルと欠落 | `scripts/03` |
| `curation_rename_log.csv` | 名寄せ履歴（どの列名・値をどう変えたか） | `notebooks/python/02` |
| `<dataset_id>_formats.txt` / `*_side_table_*.txt` | ネスト tar の形式・付随テーブルの中身 | `notebooks/python/01` のローダー |

`<acc>` は source_accession（例 `GSE208629`）。
