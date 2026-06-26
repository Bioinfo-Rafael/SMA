# 08 classical full inner clustering and microglia reclustering

## 目的

この解析は、`04d_merge_qc_original_scale.ipynb` で作成されたfull inner-gene AnnDataを用いて、古典的なScanpy workflowにより全細胞クラスタリングを行い、marker geneを確認したうえで、人間が手動で細胞種annotationを行うための解析である。

その後、手動annotationに基づいてmicroglia-like clusterを抽出し、その細胞集団だけを再クラスタリングして、microglia / DAM / IFN応答 / MHC-II / stress / contaminationなどの細かい状態を再度確認する。

この08ではscVIは使わない。PCA後にHarmony（batch補正）を入れた古典的Scanpy workflow（HVG -> scale -> PCA -> Harmony -> kNN -> UMAP -> Leiden）のみを実行する。Leiden resolutionは1.5のみを使用する。

1回目（全細胞クラスタリング）と2回目（microglia再クラスタリング）は、`--pass` 引数で明示的に指定する。

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

## ファイル配置

スクリプトとこのREADMEは `08` フォルダにまとめている。

```text
v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py
v2/notebooks/python/08/08_README_classical_full_inner_and_microglia_reclustering.md
```

## 入力ファイル

```text
v2/data/merged_h5ad/merged_qc_original_scale_inner.h5ad
```

このファイルはfull inner-gene objectであり、共通遺伝子を保持している。一方で、`.X` はoriginal-scaleであり、以下の前処理状態が混在している可能性がある。

```text
raw_count_like
cpm_tpm_like
log_normalized_like
```

そのため、`.X` をそのままPCA/UMAP/clusteringには使わない。script内で `logexpr_for_clustering` layerを作成し、それを使ってHVG selection, PCA, kNN, UMAP, Leiden clusteringを行う。`.X`（original-scale）と `.var` は保持する。

`logexpr_for_clustering` layer の作り方:

- `raw_count_like` / `cpm_tpm_like` : `sc.pp.normalize_total(target_sum=1e4)` -> `sc.pp.log1p`
- `log_normalized_like` : そのまま
- `qc_preprocessing_state` 列が無い場合 : 警告のうえ全細胞に `normalize_total + log1p`

この layer はクラスタリング・UMAP・marker 可視化・探索的 marker 検出のためだけに使い、厳密な count データとしては使わない。

## 出力ディレクトリ

出力フォルダは step 番号（08）を先頭に付ける。

```text
v2/results/08_classical_full_inner_microglia_reclustering/
```

サブディレクトリも先頭にナンバリングする。

```text
01_reports/
02_full_clustering/
03_manual_annotation/
04_microglia_reclustering/
```

主な出力は以下。

```text
01_reports/08_input_qc_report.txt
01_reports/marker_presence_full_inner.csv
02_full_clustering/full_inner_classical_clustered.h5ad
02_full_clustering/marker_genes/
02_full_clustering/plots/
03_manual_annotation/manual_annotation_template_full_clustering.csv
03_manual_annotation/full_inner_with_manual_annotation.h5ad
03_manual_annotation/microglia_like_from_manual_annotation.h5ad
03_manual_annotation/microglia_selection_summary.csv
03_manual_annotation/manual_annotation_template_microglia_reclustering.csv
04_microglia_reclustering/microglia_classical_reclustered.h5ad
04_microglia_reclustering/marker_genes/
04_microglia_reclustering/plots/
```

## 実行方法

1回目 / 2回目は `--pass` で明示的に指定する（`--pass` は必須引数）。

repository rootから実行する場合：

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

`v2/notebooks/python/08` に移動してから実行する場合：

```bash
cd v2/notebooks/python/08
python 08_classical_full_inner_and_microglia_reclustering.py --pass 1
python 08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

script内でproject rootを自動探索するため、`Path.cwd()` に依存しない。`__file__` / cwd の親ディレクトリ / 環境変数 `SMA_ROOT` をたどり、入力ファイルが存在する場所をrootとみなす。任意のディレクトリから実行する場合:

```bash
export SMA_ROOT=/path/to/SMA
python /path/to/SMA/v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
```

以下のような誤ったpathを作らないこと。

```text
v2/notebooks/python/08/v2/results/...
```

## 1回目の実行（--pass 1）で行われること

`--pass 1` を指定して実行すると、以下を行う。

1. `merged_qc_original_scale_inner.h5ad` を読み込む
2. 入力AnnDataのshape, obs columns, var columns, metadataの分布を確認する（`01_reports/08_input_qc_report.txt`）
3. marker geneがfull inner objectに存在するか確認する（`01_reports/marker_presence_full_inner.csv`）
4. `logexpr_for_clustering` layerを作成する（`.X` はoriginal-scaleのまま保持）
5. HVG 3000を使ってscale, PCA, Harmony, kNN, UMAP, Leiden clusteringを行う（batch_keyは `source_accession` → `dataset_id` の順。2水準以上ある列でHarmony補正。無い場合はHarmonyをskipして`X_pca`を使う）
6. Leiden resolution 1.5 の cluster を出力する（`leiden_r1_5`）
7. 各clusterのmarker geneを出力する（full inner genesを使用。`02_full_clustering/marker_genes/`）
8. UMAP, dotplot, tracksplot, marker feature plotを出力する（`02_full_clustering/plots/`）
9. 手動annotation用CSVを出力する

`--pass 1` の実行後、以下のファイルが作成される。

```text
v2/results/08_classical_full_inner_microglia_reclustering/03_manual_annotation/manual_annotation_template_full_clustering.csv
```

## 手動annotationのやり方

`manual_annotation_template_full_clustering.csv` を開き、clusterごとにmarker geneとplotを確認して、以下の列を手動で埋める。

```text
manual_annotation
include_for_microglia_recluster
notes
```

例：

```text
manual_annotation: Microglia
include_for_microglia_recluster: TRUE
notes: P2ry12, Cx3cr1, Hexb positive
```

```text
manual_annotation: DAM-like microglia
include_for_microglia_recluster: TRUE
notes: Apoe, Tyrobp, Trem2, Lpl positive
```

```text
manual_annotation: Astrocyte
include_for_microglia_recluster: FALSE
notes: Gfap, Aqp4 positive
```

```text
manual_annotation: Monocyte/macrophage contamination
include_for_microglia_recluster: TRUE or FALSE
notes: Ccr2, Ly6c2, Lyz2, S100a8 positive
```

microglia-likeとして再クラスタリングしたいclusterは、`include_for_microglia_recluster` に `TRUE` と書く（`TRUE` / `true` / `1` / `yes` などをTRUEと解釈する）。`include_for_microglia_recluster` が空の場合は、fallbackとして `manual_annotation` に `microglia` / `DAM` / `myeloid` / `macrophage`（大文字小文字無視）を含むclusterを選択する。

手動annotationが終わったら、ファイル名を以下に変更またはコピーして保存する。

```text
v2/results/08_classical_full_inner_microglia_reclustering/03_manual_annotation/manual_annotation_filled_full_clustering.csv
```

## 2回目の実行（--pass 2）で行われること

手動annotation済みファイルを保存した後、`--pass 2` を指定して同じscriptを実行する。

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

`--pass 2` では、scriptが以下のファイルを読み込む。

```text
02_full_clustering/full_inner_classical_clustered.h5ad   （--pass 1 の出力。無ければエラー）
03_manual_annotation/manual_annotation_filled_full_clustering.csv   （無ければメッセージを出して終了）
```

filled CSVが存在する場合、以下を行う。

1. 手動annotationをAnnDataのobsに戻す（cluster_colごとに `manual_annotation_<col>` / `include_for_microglia_recluster_<col>` を追加）
2. `include_for_microglia_recluster == TRUE` のclusterを抽出する（空ならkeyword fallback）
3. microglia-like subset AnnDataを保存する
4. microglia-like subsetだけで再度HVG, scale, PCA, Harmony, kNN, UMAP, Leiden clusteringを行う
5. microglia subclusterごとのmarker geneを出力する
6. microglia subset用のUMAP, dotplot, tracksplot, marker feature plotを出力する
7. microglia subcluster用の手動annotation templateを出力する

出力される主なファイル：

```text
03_manual_annotation/full_inner_with_manual_annotation.h5ad
03_manual_annotation/microglia_like_from_manual_annotation.h5ad
03_manual_annotation/microglia_selection_summary.csv
04_microglia_reclustering/microglia_classical_reclustered.h5ad
03_manual_annotation/manual_annotation_template_microglia_reclustering.csv
```

## microglia再クラスタリング後の手動annotation

microglia subsetの再クラスタリング後、以下のファイルを開く。

```text
v2/results/08_classical_full_inner_microglia_reclustering/03_manual_annotation/manual_annotation_template_microglia_reclustering.csv
```

microglia subclusterごとに、以下を確認して手動annotationする。

* Homeostatic microglia
* DAM-like / activated microglia
* IFN-response microglia
* MHC-II high microglia
* Proliferating microglia
* Stress-response cluster
* Monocyte/macrophage contamination
* Low-quality / ambiguous cluster

特に以下のmarker groupを確認する。

```text
Homeostatic:
P2ry12, Tmem119, Cx3cr1, Sall1

DAM / activated:
Apoe, Tyrobp, Trem2, Gpnmb, Lpl, Cst7, Cd68

Complement:
C1qa, C1qb, C1qc

IFN:
Ifit1, Ifit2, Ifit3, Isg15, Irf7, Stat1, Mx1, Oasl2, Usp18, Cxcl10

MHC-II:
H2-Aa, H2-Ab1, H2-Eb1, Cd74, B2m

Contamination:
Ccr2, Ly6c2, Lyz2, S100a8, S100a9, Mrc1, Lyve1, Pf4, Cd163
```

## 注意点

この08解析は探索的解析である。

Leiden clusterのmarker geneは、`scanpy.tl.rank_genes_groups` による探索的cluster markerであり、厳密なcondition DEGやdisease DEGではない。

疾患条件間のDEGを行う場合は、sample単位でraw countを集約したpseudobulk解析を別途行う必要がある。

また、HVGはPCA/UMAP/clusteringのためだけに使う。marker gene出力やmarker可視化では、HVGだけでなくfull inner genesを使う。

`merged_qc_original_scale_inner.h5ad` の `.X` はoriginal-scaleであり、raw count-like, CPM/TPM-like, log-normalized-likeが混ざっている可能性がある。そのため、`.X` をそのままPCA/UMAP/clusteringに使ってはいけない。

このscriptでは、`.X` を保持したまま、`logexpr_for_clustering` layerを作成し、それをクラスタリングと可視化に使う。

## 推奨される作業手順

まず `--pass 1` を実行する。

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 1
```

次に、以下を確認する。

```text
02_full_clustering/plots/
02_full_clustering/marker_genes/
03_manual_annotation/manual_annotation_template_full_clustering.csv
```

marker geneとplotを見ながら、`manual_annotation_template_full_clustering.csv` を手動で編集する。

編集後、以下の名前で保存する。

```text
03_manual_annotation/manual_annotation_filled_full_clustering.csv
```

その後、同じscriptを `--pass 2` で実行する。

```bash
python v2/notebooks/python/08/08_classical_full_inner_and_microglia_reclustering.py --pass 2
```

microglia-like subsetの再クラスタリング結果を確認する。

```text
04_microglia_reclustering/plots/
04_microglia_reclustering/marker_genes/
03_manual_annotation/manual_annotation_template_microglia_reclustering.csv
```

最後に、microglia subclusterについて再度手動annotationする。

---

# 08b submicroglia parameter sweep

## 目的

08 pass2 で得た microglia 再クラスタリング結果のうち、特定の cluster 群（`microglia_leiden_r1_5` の
3, 7-17, 20）だけをさらに抽出し、**PCA次元 × kNN近傍数 × Leiden resolution の 3×3×3 = 27 通り**で
再クラスタリング・可視化・composition 集計を行う探索スクリプト。

意図は、これらの subcluster 構造が PCA 次元数・kNN 近傍数・Leiden resolution に対して
どれくらい安定か、また Condition / dataset_id に偏った cluster が出ていないかを確認すること。

scriptは `08b_submicroglia_parameter_sweep.py`（`08/` フォルダ内）。08 本体は変更しない。

## 入力

```text
v2/results/08_classical_full_inner_microglia_reclustering/04_microglia_reclustering/microglia_classical_reclustered.h5ad
```

この AnnData は `.obs["microglia_leiden_r1_5"]`、`.layers["logexpr_for_clustering"]`、
full inner genes の var_names を持つ想定。`logexpr_for_clustering` layer を再クラスタリング・plot・
marker 表示に使う。`.X`（original-scale）は上書きしない（解析用 copy で `.X = layers[...]`）。

## 抽出する cluster

```python
SELECT_CLUSTERS = ["3"] + [str(i) for i in range(7, 18)] + ["20"]   # 3, 7..17, 20
```

`microglia_leiden_r1_5` は文字列として比較する。存在しない選択 cluster は警告のうえスキップ、
0 細胞なら停止。

## パラメータ条件（27通り）

```python
PCA_DIMS         = [20, 30, 50]
N_NEIGHBORS_LIST = [5, 15, 20]
RESOLUTIONS      = [0.5, 1.0, 1.5]
```

PCA 次元が細胞数・遺伝子数より大きい場合は自動的に `min(n_pcs, n_obs-1, n_hvg-1)` に縮める
（ディレクトリ名は指定値のまま。実際に使った次元は summary の `actual_n_pcs` に記録）。

処理フロー（各条件）:

```text
logexpr_for_clustering を .X に入れる
  → HVG (min(3000, n_vars-1)) → scale → PCA → Harmony → kNN → UMAP → Leiden
```

Harmony の batch_key は 08 と同じ方針（`source_accession` → `dataset_id` の順、2 水準以上で採用。
どちらも無ければ Harmony を skip し `X_pca` で neighbors）。Harmony は `scanpy.external.pp.harmony_integrate`
を優先し、必要なら harmonypy を直接呼んで shape 差を吸収する（08 の `harmony_integrate_to_obsm()` と同等）。

cluster column 名は条件が分かる形式:

```text
submicro_pca{PCA}_knn{K:02d}_res{res_tag}      例: submicro_pca20_knn05_res0_5
```

## 出力（h5ad は保存しない）

出力 root:

```text
v2/results/08_classical_full_inner_microglia_reclustering/05_selected_microglia_parameter_sweep/
```

**h5ad は一切保存しない**（選択 subset・各条件の再クラスタリング結果とも in-memory のみ）。
出力は CSV・PNG のみ。ディレクトリ構成:

```text
05_selected_microglia_parameter_sweep/
  summary_all_parameter_sets.csv
  pca20/
    knn05/
      res0_5/
        plots/         # umap_by_{cluster_key}.png / umap_by_Condition.png / umap_by_dataset_id.png
        dotplots/      # dotplot/tracksplot (all_markers, priority1, priority2) + marker_presence_by_*.csv
        markers/       # markers_{cluster_key}.csv（探索的 cluster marker。condition DEG ではない）
        composition/   # counts / cluster-wise / meta-wise fraction の CSV と heatmap、cluster_sizes_*.csv
      res1_0/
      res1_5/
    knn15/ ...
    knn20/ ...
  pca30/ ...
  pca50/ ...
```

各条件で得られるもの:

1. UMAP（new cluster / Condition / dataset_id）
2. dotplot / tracksplot（all_markers = priority1+priority2+celltype, priority1, priority2）
3. marker presence CSV（gene の case-insensitive 解決結果）
4. cluster marker CSV（`rank_genes_groups`, wilcoxon, full inner genes。探索的）
5. composition（Condition / dataset_id 別の counts・cluster-wise fraction・meta-wise fraction、
   06 互換の `by_Condition` / `by_dataset_id` の counts+fraction、各 heatmap、cluster size table）

`summary_all_parameter_sets.csv` に全 27 条件の
`pca_requested, actual_n_pcs, n_neighbors, resolution, cluster_key, n_cells, n_genes, n_hvg,
n_clusters, batch_key, harmony_used, output_dir, marker_csv_path, status, error_message` を記録する。
1 条件が失敗しても全体は止めず、その条件を `status="failed"` として記録し次へ進む。

## 実行方法

08 pass2（`microglia_classical_reclustered.h5ad`）が作成済みであることが前提。

```bash
cd /home/suzuki/Learn/SMA/v2/notebooks/python/08
python 08b_submicroglia_parameter_sweep.py
```

repository root から実行してもよい（project root は自動検出。`SMA_ROOT` 環境変数も使える）:

```bash
python v2/notebooks/python/08/08b_submicroglia_parameter_sweep.py
```

実行後、まず `summary_all_parameter_sets.csv` で各条件の cluster 数・harmony_used・status を確認し、
`pca*/knn*/res*/plots/` の UMAP と `composition/` の heatmap を見て、安定した PCA/kNN/resolution の
組み合わせや、Condition/dataset_id に偏った cluster の有無を判断する。

## 注意点

- この 08b も探索的解析であり、`rank_genes_groups` の marker は探索的 cluster marker。厳密な
  condition / disease DEG ではない（pseudobulk を別途行うこと）。
- HVG は PCA/UMAP/clustering 用のみ。marker 出力・可視化は full inner genes 全体で行う。
- h5ad を保存しないため、再クラスタリング結果を再利用したい場合は条件を絞って再実行する。
- 27 条件すべてを回すため計算時間がかかる。途中失敗は summary に記録され、残りは継続する。
