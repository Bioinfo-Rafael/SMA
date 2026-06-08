"""Step 07 — produce all figures: UMAP-by-GSE, per-resolution UMAP-by-Leiden / by-annotation,
marker dot-plot, GSE-vs-annotation composition stacked bar."""
import sys
import json
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, setup_figdir, MARKERS


def main():
    setup_figdir()
    adata = read_stage('annot')

    # 1) UMAP coloured by GSE (study-of-origin)
    sc.pl.umap(adata, color='gse_id', save='_by_gse.png', show=False,
               legend_loc='right margin', title='UMAP by GSE')
    plt.close('all')

    # 2) Per resolution: Leiden cluster id and marker-derived annotation
    for r in config.RESOLUTIONS:
        sc.pl.umap(adata, color=f'leiden_r{r}', save=f'_leiden_r{r}.png', show=False,
                   legend_loc='on data', legend_fontsize=6,
                   title=f'Leiden r={r}  ({adata.obs[f"leiden_r{r}"].nunique()} clusters)')
        sc.pl.umap(adata, color=f'annotation_r{r}', save=f'_annotation_r{r}.png', show=False,
                   legend_loc='right margin', title=f'Annotation r={r}')
        plt.close('all')

    # 3) Marker dot-plot + composition at the middle resolution
    mid = config.RESOLUTIONS[len(config.RESOLUTIONS) // 2]
    markers = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
    markers = {k: v for k, v in markers.items() if v}
    sc.pl.dotplot(adata, markers, groupby=f'annotation_r{mid}',
                  save=f'_markers_r{mid}.png', show=False, standard_scale='var')
    plt.close('all')

    comp = (pd.crosstab(adata.obs['gse_id'], adata.obs[f'annotation_r{mid}'])
              .apply(lambda r: r / r.sum() * 100, axis=1))
    ax = comp.plot(kind='bar', stacked=True, figsize=(9, 4), width=0.85)
    ax.set_ylabel('% of cells')
    ax.set_title(f'Composition per GSE  (annotation r={mid})')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.savefig(config.FIG_DIR / f'composition_r{mid}.png', dpi=200)
    plt.close('all')

    summary = {
        'n_cells'           : int(adata.n_obs),
        'n_genes'           : int(adata.n_vars),
        'gse_counts'        : adata.obs['gse_id'].value_counts().to_dict(),
        'resolutions'       : config.RESOLUTIONS,
        'clusters_per_res'  : {f'leiden_r{r}': int(adata.obs[f'leiden_r{r}'].nunique())
                                for r in config.RESOLUTIONS},
        'annotation_counts' : {f'r{r}': pd.Series(adata.obs[f'annotation_r{r}']).value_counts().to_dict()
                                for r in config.RESOLUTIONS},
    }
    (config.RES_DIR / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print('  summary written ->', config.RES_DIR / 'summary.json')
    print(f'  figures        ->', config.FIG_DIR)


if __name__ == '__main__':
    main()
