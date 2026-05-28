# SMA / ALS scRNA-seq integration pipeline

6つのマウス脊髄 scRNA-seq / snRNA-seq データセット (SMA, ALS, SOD1G93A, etc.) を
ダウンロード → 同一スキーマに整形 → scVI で統合 → クラスタリング → アノテーションする
End-to-End パイプライン。

## データセット (GEO)

| GSE | 内容 | 形式 |
|-----|------|------|
| [GSE287569](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287569) | SOD1G93A spinal cord scRNA-seq (12 samples) | 10x H5 |
| [GSE173524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173524) | SOD1G93A spinal cord snRNA-seq | 統合UMI tsv + 細胞メタtsv (cell-type付) |
| [GSE167332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167332) | spinal cord scRNA-seq (plate-based) | dense matrix + MTX |
| [GSE219201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219201) | Phd1-KO spinal cord snRNA-seq | 10x MTX |
| [GSE242942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242942) | SOD1G93A + PF-04457845 scRNA-seq | 10x MTX |
| [GSE208629](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208629) | SMA spinal cord scRNA-seq | 10x MTX |

## 構成

```
SMA/
├── unified_preprocessing.ipynb   # ① Colab: GEO から DL → 統一フォーマットの h5ad
├── data/                         # h5ad 入力 (.gitignore)
├── scripts/                      # ② ローカル解析パイプライン (1 step = 1 file)
│   ├── config.py                 #   ハイパラ集約 (env var で上書き可)
│   ├── common.py                 #   共通ユーティリティ + マーカー辞書
│   ├── 01_concat.py              #   全 GSE をロード → var_names UPPER → concat
│   ├── 02_qc.py                  #   セル/遺伝子フィルタ + HVG (per-batch)
│   ├── 03_scvi.py                #   scVI 学習 (batch_key=gse_id)
│   ├── 04_umap.py                #   neighbors + UMAP (scVI latent)
│   ├── 05_cluster.py             #   多解像度 Leiden
│   ├── 06_annotate.py            #   マーカースコアでクラスタアノテーション
│   └── 07_plot.py                #   UMAP/dotplot/composition 出力
├── run.sh                        # オーケストレータ
├── results/                      # 中間 h5ad / モデル / 図 (.gitignore)
└── requirements.txt
```

## ① 前処理 (Colab) — `unified_preprocessing.ipynb`

GEO 各 GSE をダウンロード → 異種フォーマット (10x H5 / MTX / dense / 統合 TSV)
を `AnnData` にロード → obs/var を統一スキーマに正規化 → `{GSE}.h5ad` で保存。

統一スキーマ:

- `obs`: `cell_id` / `gse_id` / `sample_id` / `condition` / `cell_type`
  (無ければ `'unknown'`)。元の obs カラムも保持。
- `var_names`: 遺伝子シンボル優先 (無ければ ENSEMBL ID)。 `gene_symbol`,
  `gene_symbol_upper`, `ensembl_id` を別カラムで保持。
- `X`: 生 UMI カウント (CSR sparse)。
- セルバーコードは `{GSM}_{barcode}` でグローバルユニーク化。

末尾でクロスデータセットの遺伝子オーバーラップ表 (大文字統一比較) を出力。

## ② ローカル解析パイプライン

### セットアップ

```bash
python3.14 -m venv .venv
.venv/bin/pip install -r requirements.txt
mkdir -p data results
# data/ に GSE*.h5ad を配置
```

### 実行

```bash
./run.sh                          # 全ステップ
./run.sh 03 04 05                 # 一部だけ再実行
MAX_EPOCHS=200 ./run.sh 03        # ハイパラ上書き
RESOLUTIONS=0.2,0.5,1.0 ./run.sh 05 06 07
```

各ステップは前ステップの h5ad を読み込み、自身の h5ad を書き出すので、途中で
パラメータを変えて部分実行できる。

### ステップ詳細

| step | 入力 | 出力 | 役割 |
|------|------|------|------|
| `01_concat`   | `data/GSE*.h5ad`     | `results/s01_concat.h5ad` | 全 GSE をロード, `var_names` を UPPER に統一して inner-join concat。生カウントを `layers['counts']` に保持。 |
| `02_qc`       | `s01_concat.h5ad`    | `results/s02_qc.h5ad`     | per-cell filter (`min_genes`, `pct_counts_mt`), per-gene filter, per-batch HVG (`seurat_v3`, `batch_key=gse_id`)。GSE242942 の空ドロップレット ~1.87M も除去。 |
| `03_scvi`     | `s02_qc.h5ad`        | `results/s03_scvi.h5ad`<br>`results/scvi_model/` | HVG 部分集合で scVI 学習 (`batch_key=gse_id`)。MPS / CUDA / CPU を自動選択。`obsm['X_scVI']` を書き戻し。 |
| `04_umap`     | `s03_scvi.h5ad`      | `results/s04_umap.h5ad`   | scVI 潜在空間で kNN + UMAP。 |
| `05_cluster`  | `s04_umap.h5ad`      | `results/s05_cluster.h5ad`| Leiden を複数解像度 (`RESOLUTIONS`) で実行し、`obs['leiden_r{r}']` に保存。 |
| `06_annotate` | `s05_cluster.h5ad`   | `results/s06_annot.h5ad`  | normalize_total + log → マーカー gene-score → クラスタの平均スコア最大のラベルを採用 (`min_score` 未満は `Unknown`)。各解像度に対し `obs['annotation_r{r}']` を生成。 |
| `07_plot`     | `s06_annot.h5ad`     | `results/figures/*.png`<br>`results/summary.json` | UMAP (GSE色分け / 各解像度の Leiden / annotation), マーカー dotplot, GSE×annotation composition stacked bar, JSON サマリ。 |

### マーカー辞書 (mouse spinal cord / CNS)

`scripts/common.py` の `MARKERS` を編集することで追加・差し替え可能 (全 UPPER)。
デフォルトは Motor neuron / Inhib neuron / Excit neuron / Neuron general /
Astrocyte / Microglia / Oligodendrocyte / OPC / Ependymal / Endothelial /
Pericyte / T cell / B cell / Schwann。

### 主要ハイパラ (env var で上書き)

| 変数 | 既定 | 効くステップ |
|------|------|--------------|
| `MIN_GENES`     | 500  | 02 |
| `MIN_CELLS`     | 10   | 02 |
| `MAX_PCT_MT`    | 20.0 | 02 |
| `N_HVG`         | 3000 | 02 |
| `BATCH_KEY`     | gse_id | 02 / 03 |
| `N_LATENT`      | 30   | 03 |
| `N_LAYERS`      | 2    | 03 |
| `MAX_EPOCHS`    | 100  | 03 |
| `LR`            | 1e-3 | 03 |
| `N_NEIGHBORS`   | 15   | 04 |
| `MIN_DIST`      | 0.3  | 04 |
| `RESOLUTIONS`   | 0.3,0.5,0.8,1.2 | 05 / 06 / 07 |
| `MIN_SCORE`     | 0.05 | 06 |

## 出力

- `results/s0X_*.h5ad` ─ 各ステップの状態 (容量重視ならステップ間は削除可)
- `results/scvi_model/` ─ 学習済み scVI モデル
- `results/figures/*.png` ─ 可視化
  - `umap_by_gse.png` (GSE 色分け)
  - `umap_leiden_r{r}.png` (各解像度のクラスタ)
  - `umap_annotation_r{r}.png` (各解像度のアノテーション)
  - `dotplot__markers_r{mid}.png`
  - `composition_r{mid}.png`
- `results/summary.json` ─ 細胞数 / クラスタ数 / アノテーション内訳
