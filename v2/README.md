# SMA / ALS scRNA-seq preprocessing (v2)

GEO **Supplementary file** based download → AnnData → manual inspection →
config-driven curation → per-dataset `.h5ad`. Built to be re-run, resumable,
and to keep humans in control of cell-type / cluster decisions.

## What this pipeline does (and does not) do

* It uses **GEO Supplementary files only**. It does **not** use SRA FASTQ.
* It does **not** scrape GEO pages. [`config/dataset_manifest.yaml`](config/dataset_manifest.yaml)
  is the single source of truth for every GSE / SubSeries / file name / URL.
* It does **not** normalize, filter, or QC. Raw counts go into `adata.X`.
* It does **not** auto-decide cell-type or cluster columns. You read the
  inspection reports and fill [`config/curation_template.yaml`](config/curation_template.yaml).
* One **logical dataset → one `.h5ad`**. Material with different properties
  (whole tissue, FACS/CD45 enriched, SmartSeq2 TPM, SoupX processed) is **not**
  merged into a single object.

### Data that is NOT raw UMI counts (handled, but flagged)

| Dataset | Nature | Recorded as |
|---|---|---|
| `GSE167331` (FACS microglia, SmartSeq2) | **TPM** matrix | `uns['data_status'] = processed_TPM` |
| `GSE206330` (cortical glia) | **SoupX corrected** processed data | `uns['data_status'] = processed_SoupX_corrected` |
| `GSE295514` (rNLS8 TDP-43) | **RDS** object (Seurat/SCE) | converted via R bridge |

These must not be pooled with raw-UMI datasets without status-aware handling.

### SuperSeries handling

* **GSE167332** → split into **3 separate** logical datasets and saved separately:
  `GSE167198` (Drop-seq whole cord), `GSE167327` (CD45-enriched, inDrop),
  `GSE167331` (FACS microglia, SmartSeq2 **TPM**).
* **GSE242942** → use **only** the scRNA-seq SubSeries **`GSE242939`**.
  The bulk RNA-seq SubSeries **`GSE242940` is not used**.
* **GSE173524** → uses the raw `GSE173524_umi.tsv.gz`, **not** the
  `*.sctransform.tsv.gz` (normalized) version.

## Layout

```
v2/
├── config/
│   ├── dataset_manifest.yaml     # SOURCE OF TRUTH (GSE, files, URLs, loaders, sample rules)
│   └── curation_template.yaml    # human-edited; cell_type/cluster columns start null
├── scripts/
│   ├── 00_make_manifest.py       # validate manifest + write overview
│   ├── 01_download_geo_supplement.py
│   ├── 02_extract_archives.py    # safe (path-traversal hardened) + nested tars
│   ├── 03_load_to_anndata.py     # dispatch loaders -> interim h5ad
│   ├── 04_inspect_anndata_columns.py  # manual-inspection reports (no auto column choice)
│   ├── 05_curate_and_save_h5ad.py     # apply curation yaml -> curated h5ad
│   ├── 06_read_saved_h5ad_template.py # integration TEMPLATE (not run automatically)
│   └── rds_to_h5ad_bridge.R       # RDS (Seurat/SCE) -> MTX triplet + meta.csv
├── src/
│   ├── geo_download.py            # resumable download + safe tar extraction
│   ├── io_10x.py                  # 10x .h5 and MTX-triplet loaders
│   ├── io_dense.py                # dense/text/combined/processed/nested loaders
│   ├── io_rds_bridge.py           # drives the R bridge, assembles AnnData
│   ├── anndata_utils.py           # obs/var schema, sample rules, save helpers
│   └── reporting.py               # logging + inspection report writer
├── data/                          # created on first run (git-ignored)
│   ├── raw/<acc>/                 # downloaded supplementary files
│   ├── extracted/<acc>/           # unpacked archives
│   ├── interim_h5ad/              # 03 output (one per dataset)
│   ├── curated_h5ad/              # 05 output
│   └── reports/                   # 04 inspection reports
├── requirements.txt
└── run.sh
```

> `data/` is git-ignored (large). The subdirectories are created automatically
> on first run.

## Usage

```bash
# 0. environment (Python 3.10+)
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
# (or reuse the v1 env:  PYTHON=../v1/.venv/bin/python ./run.sh)

# 1. validate the manifest
python scripts/00_make_manifest.py

# 2. download supplementary files (resumable; skips complete files)
python scripts/01_download_geo_supplement.py
#    one dataset only:
python scripts/01_download_geo_supplement.py --datasets GSE208629_sc_spinalcord_sma

# 3. extract archives (safe; nested tars handled)
python scripts/02_extract_archives.py

# 4. build interim h5ad (one per dataset)
python scripts/03_load_to_anndata.py

# 5. write inspection reports, then READ them
python scripts/04_inspect_anndata_columns.py
#    -> data/reports/<dataset_id>_summary.txt
#       data/reports/<dataset_id>_obs_columns.csv
#       data/reports/<dataset_id>_obs_value_counts/<col>.csv
#       data/reports/<dataset_id>_var_columns.csv

# 6. edit config/curation_template.yaml (set cell_type_column / cluster_column
#    based on what you saw), then:
python scripts/05_curate_and_save_h5ad.py

# 7. later integration work starts from the template (not run automatically)
python scripts/06_read_saved_h5ad_template.py
```

`./run.sh` runs steps 00–05 in order. `./run.sh 03 04` runs a subset.

### The RDS dataset (GSE295514)

`GSE295514_ALS_mouse_brain.rds` is converted by `scripts/rds_to_h5ad_bridge.R`,
which is invoked automatically by `io_rds_bridge.py` during step 03. It needs
**R** on `PATH` (`Rscript`) plus the `Matrix` package, and `Seurat`/
`SeuratObject` (for Seurat objects) or `SingleCellExperiment` (for SCE). If the
object is neither, the bridge writes `data/reports/GSE295514_*_rds_bridge.log`
and the dataset is skipped with a clear error (the rest of the pipeline
continues).

## Metadata conventions

Every interim h5ad carries a standard `obs` schema (sample/condition/genotype/
treatment/enrichment/tissue/technology/processing_status/data_status/… — see
`src/anndata_utils.REQUIRED_OBS_COLS`) and `var` schema (`gene_id`,
`gene_symbol`, `gene_symbol_upper`, `ensembl_id`, `feature_type`). `obs_names`
are globally unique: `{source_accession}_{sample_id}_{original_barcode}`.

Condition fields that cannot be derived from file names (e.g. for the
metadata-driven `GSE173524`, `GSE206330`, `GSE295514`) are left as `unknown`
on purpose — the merged original metadata is preserved so you can finalize them
in the curation yaml after reading the inspection reports.
