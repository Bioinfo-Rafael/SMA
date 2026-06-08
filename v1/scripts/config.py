"""Central configuration. Edit values here; every step reads from this module.

Override at run-time with environment variables, e.g.:
    MAX_EPOCHS=50 python scripts/03_scvi.py
"""
import os
from pathlib import Path

ROOT      = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / 'data'
RES_DIR   = ROOT / 'results'
FIG_DIR   = RES_DIR / 'figures'
MODEL_DIR = RES_DIR / 'scvi_model'

# ---- QC ----
MIN_GENES_PER_CELL = int(os.getenv('MIN_GENES', 500))  # 10x raw matrices need this to drop empty droplets
MIN_CELLS_PER_GENE = int(os.getenv('MIN_CELLS', 10))
MAX_PCT_MT         = float(os.getenv('MAX_PCT_MT', 20.0))

# ---- HVG ----
N_HVG       = int(os.getenv('N_HVG', 3000))
HVG_FLAVOR  = os.getenv('HVG_FLAVOR', 'seurat_v3')
BATCH_KEY   = os.getenv('BATCH_KEY', 'gse_id')

# ---- scVI ----
N_LATENT    = int(os.getenv('N_LATENT', 30))
N_LAYERS    = int(os.getenv('N_LAYERS', 2))
MAX_EPOCHS  = int(os.getenv('MAX_EPOCHS', 100))
LR          = float(os.getenv('LR', 1e-3))
TRAIN_SIZE  = float(os.getenv('TRAIN_SIZE', 0.9))

# ---- Neighbors / UMAP ----
N_NEIGHBORS = int(os.getenv('N_NEIGHBORS', 15))
MIN_DIST    = float(os.getenv('MIN_DIST', 0.3))

# ---- Multi-resolution Leiden ----
RESOLUTIONS = [float(x) for x in os.getenv('RESOLUTIONS', '0.3,0.5,0.8,1.2').split(',')]

# Stage filenames (under RES_DIR) used as the chain between steps
STAGES = {
    'concat'  : RES_DIR / 's01_concat.h5ad',
    'qc'      : RES_DIR / 's02_qc.h5ad',
    'scvi'    : RES_DIR / 's03_scvi.h5ad',
    'umap'    : RES_DIR / 's04_umap.h5ad',
    'cluster' : RES_DIR / 's05_cluster.h5ad',
    'annot'   : RES_DIR / 's06_annot.h5ad',
}
