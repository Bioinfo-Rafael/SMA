"""Step 05 — Leiden at multiple resolutions; one column per resolution."""
import sys
from pathlib import Path
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, write_stage


def main():
    adata = read_stage('umap')
    for r in config.RESOLUTIONS:
        key = f'leiden_r{r}'
        sc.tl.leiden(adata, resolution=r, key_added=key,
                     flavor='igraph', directed=False, n_iterations=2)
        print(f'  {key}: {adata.obs[key].nunique()} clusters')
    write_stage(adata, 'cluster')


if __name__ == '__main__':
    main()
