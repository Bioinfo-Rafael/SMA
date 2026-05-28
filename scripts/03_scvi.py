"""Step 03 — train scVI on HVGs with batch_key=gse_id; write the latent into obsm."""
import sys
from pathlib import Path
import scvi

sys.path.insert(0, str(Path(__file__).parent))
import config
from common import read_stage, write_stage, get_accelerator


def main():
    adata = read_stage('qc')
    accel = get_accelerator()
    print(f'  accelerator = {accel}  | n_cells={adata.n_obs:,}  n_hvg={int(adata.var.highly_variable.sum())}')

    adata_hvg = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, layer='counts', batch_key=config.BATCH_KEY)
    model = scvi.model.SCVI(
        adata_hvg,
        n_layers=config.N_LAYERS,
        n_latent=config.N_LATENT,
        gene_likelihood='nb',
    )
    model.view_anndata_setup()
    model.train(
        max_epochs=config.MAX_EPOCHS,
        early_stopping=True,
        accelerator=accel,
        devices='auto',
        train_size=config.TRAIN_SIZE,
        plan_kwargs={'lr': config.LR},
    )

    config.MODEL_DIR.mkdir(parents=True, exist_ok=True)
    model.save(str(config.MODEL_DIR), overwrite=True)
    print(f'  saved model -> {config.MODEL_DIR}')

    adata.obsm['X_scVI'] = model.get_latent_representation()
    print(f'  latent shape = {adata.obsm["X_scVI"].shape}')
    write_stage(adata, 'scvi')


if __name__ == '__main__':
    main()
