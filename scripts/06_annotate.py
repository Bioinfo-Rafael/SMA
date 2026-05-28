"""Step 06 — score canonical marker sets, assign each cluster to the highest-scoring set.

Produces one annotation column per resolution: `annotation_r{resolution}`.
Clusters whose top marker score is below MIN_SCORE are labelled 'Unknown'.
"""
import sys
from pathlib import Path
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, write_stage, MARKERS

MIN_SCORE = float(__import__('os').getenv('MIN_SCORE', '0.05'))


def annotate_clusters(adata, cluster_key, prefix='score_', min_score=MIN_SCORE):
    score_cols = [c for c in adata.obs.columns if c.startswith(prefix)]
    means = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    best     = means.idxmax(axis=1).str.replace(prefix, '', regex=False)
    best_val = means.max(axis=1)
    best[best_val < min_score] = 'Unknown'
    return adata.obs[cluster_key].map(best).astype('category')


def main():
    adata = read_stage('cluster')

    # Normalise + log for marker scoring (raw counts already kept in layer 'counts').
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    markers = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
    markers = {k: v for k, v in markers.items() if v}
    print('  marker sets retained:', {k: len(v) for k, v in markers.items()})

    for name, genes in markers.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f'score_{name}')

    for r in config.RESOLUTIONS:
        adata.obs[f'annotation_r{r}'] = annotate_clusters(adata, f'leiden_r{r}')
        counts = pd.Series(adata.obs[f'annotation_r{r}']).value_counts()
        print(f'  r={r}:', counts.to_dict())

    write_stage(adata, 'annot')


if __name__ == '__main__':
    main()
