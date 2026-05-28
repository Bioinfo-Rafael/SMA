"""Step 02 — per-cell / per-gene QC + HVG selection.

Drops empty droplets (notably GSE242942 ships the raw matrix with ~1.8M barcodes)
and selects HVGs per batch so dominant studies do not capture the entire feature set.
"""
import sys
from pathlib import Path
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, write_stage


def main():
    adata = read_stage('concat')

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=config.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS_PER_GENE)
    adata = adata[adata.obs['pct_counts_mt'] < config.MAX_PCT_MT].copy()
    print(f'  QC: cells {before:,} -> {adata.n_obs:,}  ({adata.n_obs/before*100:.1f}% kept)')
    print(f'      genes after filter: {adata.n_vars:,}')
    print('      per GSE after QC:', adata.obs['gse_id'].value_counts().to_dict())

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.N_HVG,
        flavor=config.HVG_FLAVOR,
        layer='counts',
        batch_key=config.BATCH_KEY,
        subset=False,
    )
    print(f'  HVG (per batch={config.BATCH_KEY}, n_top={config.N_HVG}) selected: '
          f'{int(adata.var.highly_variable.sum())}')

    write_stage(adata, 'qc')


if __name__ == '__main__':
    main()
