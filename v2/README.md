# SMA / ALS scRNA-seq preprocessing (v2)

**download scripts + interactive notebooks.**

The `.py` scripts only fetch and organise GEO **Supplementary files** (download
→ extract → list). Everything else — reading data into AnnData, inspecting
obs/var, curating metadata, saving per-GSE h5ad, and merging — happens in
**Jupyter notebooks**, so a human stays in the loop and decides cell-type /
cluster / condition columns. Nothing is auto-detected.

## What runs where

| Stage | Where |
|---|---|
| validate manifest, download, extract, list files | `.py` scripts (`scripts/`) |
| load files → AnnData | notebook `python/01` |
| inspect obs/var/layers/uns | notebook `python/02` |
| curate metadata (manual) | notebook `python/03` |
| save per-GSE curated h5ad | notebook `python/03` |
| merge curated h5ad → merged h5ad | notebook `python/04` |
| check merged h5ad | notebook `python/05` |
| open & export the RDS dataset | **R** notebook `R/01_GSE295514_read_rds.ipynb` |

> Helper *functions* live in `src/`; the *execution unit* is the notebook.
> There is **no** `03_load`/`04_inspect`/`05_curate` script and `run.sh` stops
> after the download/extract/list steps.

## Key data caveats

* **GSE242942** → use only the scRNA-seq SubSeries **`GSE242939`**; the bulk
  RNA-seq SubSeries **`GSE242940` is not used**.
* **GSE167332** → 3 SubSeries are **all included but saved as separate h5ad**:
  `GSE167198` (Drop-seq whole cord), `GSE167327` (CD45-enriched, inDrop),
  `GSE167331` (FACS microglia, SmartSeq2).
* **GSE167331** is a **TPM** matrix → `data_status = processed_TPM`; never pooled
  with raw counts.
* **GSE206330** is **SoupX corrected** processed data →
  `data_status = processed_SoupX_corrected`; kept separate from raw counts.
* **GSE295514** is an **RDS** object: open it in the R Jupyter kernel, inspect,
  and export intermediates that the Python notebook reads
  (`data_status = RDS_converted_unknown_or_counts`).
* **GSE173524** uses the raw `GSE173524_umi.tsv.gz`, not the `*.sctransform.*`
  (normalised) version.
* Merged h5ad files are written to `data/merged_h5ad/`. Merges are
  **status-aware**: one merge of only `raw_or_filtered_count` datasets, and one
  of everything with `data_status` preserved in obs.

## Layout

```
v2/
├── config/dataset_manifest.yaml      # SOURCE OF TRUTH (GSE/files/URLs/loader_hint/metadata)
├── scripts/
│   ├── 00_validate_manifest.py
│   ├── 01_download_geo_supplement.py
│   ├── 02_extract_archives.py        # safe (path-traversal) + nested tars
│   └── 03_list_downloaded_files.py
├── src/                              # functions used BY the notebooks
│   ├── geo_download.py               # resumable download
│   ├── archive_utils.py              # safe tar extraction + find_files
│   ├── manifest_utils.py             # load/validate manifest, paths, logging
│   ├── io_10x.py                     # read_10x_h5_file / read_10x_mtx_triplet / loaders
│   ├── io_dense.py                   # dense/combined/processed/nested + R-intermediate reader
│   ├── anndata_utils.py              # obs/var schema, obs_names, save/load, merge
│   └── notebook_report_utils.py      # summarize_adata / show_* helpers
├── notebooks/
│   ├── python/
│   │   ├── 00_overview.ipynb
│   │   ├── 01_load_each_gse_to_anndata.ipynb
│   │   ├── 02_inspect_each_gse_anndata.ipynb
│   │   ├── 03_curate_each_gse_and_save_h5ad.ipynb
│   │   ├── 04_merge_curated_h5ad.ipynb
│   │   └── 05_check_merged_h5ad.ipynb
│   └── R/
│       └── 01_GSE295514_read_rds.ipynb
├── data/                             # git-ignored; created on first run
│   ├── raw/<acc>/                    # downloaded supplementary files
│   ├── extracted/<acc>/              # unpacked archives
│   ├── intermediate_from_r/<acc>/    # R notebook exports (counts.mtx, metadata.csv, ...)
│   ├── interim_h5ad/                 # optional raw AnnData (notebook 01)
│   ├── curated_h5ad/                 # notebook 03 output
│   ├── merged_h5ad/                  # notebook 04 output
│   └── reports/                      # manifest overview, file lists, side tables
├── requirements.txt
└── run.sh                            # download/extract/list ONLY
```

## Usage

### 1. Python environment

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name sma-v2   # optional: named kernel
# (or reuse the v1 env:  PYTHON=../v1/.venv/bin/python ./run.sh)
```

### 2. Download + extract (scripts)

```bash
python scripts/00_validate_manifest.py
python scripts/01_download_geo_supplement.py        # resumable; --datasets GSE208629 to limit
python scripts/02_extract_archives.py
python scripts/03_list_downloaded_files.py
# or simply:  ./run.sh
```

### 3. AnnData onwards (notebooks)

Open JupyterLab and work through `notebooks/python/` in order:
`00_overview` → `01_load…` → `02_inspect…` → `03_curate…` → `04_merge…` →
`05_check…`. For **GSE295514**, first run the **R** notebook
`notebooks/R/01_GSE295514_read_rds.ipynb`, then load it in `python/01`.

```bash
jupyter lab
```

## R notebook (GSE295514 RDS)

`notebooks/R/01_GSE295514_read_rds.ipynb` runs on an **R Jupyter kernel**. It
reads the RDS, inspects class/assays/`meta.data`, and exports
`counts.mtx` / `metadata.csv` / `genes.csv` / `barcodes.csv` to
`data/intermediate_from_r/GSE295514/`, which `io_dense.read_from_r_intermediate`
turns into AnnData.

R packages that may be needed (install yourself in R; not auto-installed):

```r
install.packages("IRkernel"); IRkernel::installspec()   # R kernel for Jupyter
install.packages("Matrix")
# plus, depending on the object:
install.packages("Seurat")              # also pulls SeuratObject
# or (Bioconductor):
# BiocManager::install("SingleCellExperiment")
# optional bridges: zellkonverter, Matrix.utils
```

## Scope

This stage covers **download, AnnData creation, manual inspection, curation,
per-GSE h5ad saving, and merged h5ad saving** only. Analysis, QC, normalization,
clustering, UMAP, scVI and batch correction are intentionally **not** included
yet.
