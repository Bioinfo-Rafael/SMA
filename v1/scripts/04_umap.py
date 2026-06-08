"""Step 04 — kNN graph + UMAP on the scVI latent space."""
import sys
from pathlib import Path
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, write_stage


def main():
    adata = read_stage('scvi')
    print(f'  neighbors  n_neighbors={config.N_NEIGHBORS}  use_rep=X_scVI')
    sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=config.N_NEIGHBORS)
    print(f'  umap       min_dist={config.MIN_DIST}')
    sc.tl.umap(adata, min_dist=config.MIN_DIST)
    write_stage(adata, 'umap')


if __name__ == '__main__':
    main()
