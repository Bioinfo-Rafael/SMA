"""Shared helpers: accelerator detection, marker dictionary, IO."""
import sys
from pathlib import Path
import anndata as ad
import scanpy as sc
import torch

sys.path.insert(0, str(Path(__file__).parent))
import config


def get_accelerator():
    if torch.backends.mps.is_available():
        return 'mps'
    if torch.cuda.is_available():
        return 'gpu'
    return 'cpu'


def setup_figdir():
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(config.FIG_DIR)
    sc.settings.set_figure_params(dpi=80, dpi_save=200, frameon=False)


def read_stage(name):
    p = config.STAGES[name]
    print(f'  read  <- {p}')
    return ad.read_h5ad(p)


def write_stage(adata, name):
    config.RES_DIR.mkdir(parents=True, exist_ok=True)
    p = config.STAGES[name]
    adata.write_h5ad(p, compression='gzip')
    print(f'  write -> {p}  ({p.stat().st_size/1e6:.1f} MB)')


# Mouse spinal-cord / CNS markers (gene symbols UPPER-cased to match the integrated var_names).
MARKERS = {
    'Motor_neuron'    : ['CHAT','MNX1','ISL1','LHX3','SLC18A3'],
    'Inhib_neuron'    : ['GAD1','GAD2','SLC6A5'],
    'Excit_neuron'    : ['SLC17A6','SLC17A7'],
    'Neuron_general'  : ['SNAP25','SYT1','STMN2','RBFOX3'],
    'Astrocyte'       : ['GFAP','AQP4','SLC1A3','SLC1A2','S100B','ALDH1L1'],
    'Microglia'       : ['P2RY12','TMEM119','CX3CR1','CSF1R','C1QA','C1QB'],
    'Oligodendrocyte' : ['MBP','MOG','MAG','PLP1','CLDN11'],
    'OPC'             : ['PDGFRA','CSPG4','SOX10','OLIG1','OLIG2'],
    'Ependymal'       : ['FOXJ1','CD24A'],
    'Endothelial'     : ['CLDN5','PECAM1','CDH5','FLT1'],
    'Pericyte'        : ['PDGFRB','MCAM','RGS5','VTN'],
    'T_cell'          : ['CD3E','CD3D','CD8A','CD4'],
    'B_cell'          : ['CD19','CD79A'],
    'Schwann'         : ['MPZ'],
}
