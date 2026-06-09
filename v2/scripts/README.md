# scripts/

GEO **Supplementary files** を取得・整理する `.py` スクリプト群。
ここで行うのは「ダウンロード → 展開 → 一覧 → 俯瞰」まで。**AnnData 化以降は
`notebooks/` で行う**（このディレクトリには置かない）。

すべて `src/` の関数を呼ぶ薄いエントリポイント。実行は単体でも `../run.sh` でも可。

## ファイル

| ファイル | 役割 | 主な出力 |
|---|---|---|
| `00_validate_manifest.py` | `config/dataset_manifest.yaml` の検証（必須キー・loader_hint・重複） | `data/reports/manifest_overview.csv` |
| `01_download_geo_supplement.py` | manifest の全 URL を `data/raw/<acc>/` にDL（レジューム・完了済みskip） | `data/raw/<acc>/*` |
| `02_extract_archives.py` | `archive:true` の tar を `data/extracted/<acc>/` へ安全展開（ネスト対応） | `data/extracted/<acc>/*` |
| `03_list_downloaded_files.py` | DL/展開済みを一覧し manifest と突合（欠落検出） | `data/reports/downloaded_files.txt` |
| `04_overview.py` | データセット一覧・各パス・ファイル有無の俯瞰（旧 00_overview ノート） | 標準出力 |

## よく使うオプション

```bash
# 1データセットだけ処理（dataset_id か source_accession）
python scripts/01_download_geo_supplement.py --datasets GSE208629
python scripts/02_extract_archives.py       --datasets GSE208629

# 取り直し / 展開やり直し
python scripts/01_download_geo_supplement.py --force
python scripts/02_extract_archives.py       --force
```

実行順や全体像は `../README.md` を参照。
