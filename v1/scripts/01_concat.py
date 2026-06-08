"""Step 01 — load every data/GSE*.h5ad, normalise gene symbols, concatenate."""
import sys
from pathlib import Path
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import write_stage


def main():
    files = sorted(config.DATA_DIR.glob('GSE*.h5ad'))
    if not files:
        sys.exit(f'no h5ad in {config.DATA_DIR}')

    adatas = {}
    for f in files:
        a = sc.read_h5ad(f)
        # Use upper-cased gene symbol as var_names so cross-study case drift cannot split genes.
        if 'gene_symbol_upper' in a.var.columns:
            a.var_names = pd.Index(a.var['gene_symbol_upper'].astype(str))
        else:
            a.var_names = pd.Index(pd.Series(a.var_names.astype(str)).str.upper().values)
        a.var_names_make_unique()
        a.obs['gse_id'] = f.stem
        print(f'  {f.stem}: cells={a.n_obs:>9,}  genes={a.n_vars:>6,}')
        adatas[f.stem] = a

    adata = ad.concat(adatas, join='inner', merge='same', index_unique=None)
    adata.obs['gse_id'] = adata.obs['gse_id'].astype('category')
    # Keep an explicit copy of raw counts in a layer (X will be normalised later).
    adata.layers['counts'] = adata.X.copy()
    print(f'\n  concat   cells={adata.n_obs:>9,}  genes={adata.n_vars:>6,}')
    print('  per GSE :', adata.obs['gse_id'].value_counts().to_dict())
    write_stage(adata, 'concat')


if __name__ == '__main__':
    main()
