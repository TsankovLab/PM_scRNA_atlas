""" -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : integrate.py

* Description : scANVI integration - main

_._._._._._._._._._._._._._._._._._._._._. """

from __imports__ import *
from __config__ import *

# create project directory
os.makedirs('scRNA/main/Plots', exist_ok=True)

# setup parameters
batch = 'sampleID'
n_features = 5000

CELLTYPE_COL = 'celltype_simplified2'

# Integrate srt_tumor
# setup paths
main_root = os.path.join('scRNA', 'main')
ScanpyConfig.figdir = Path(os.path.join(main_root, 'Plots'))

scanvi_model_path = os.path.join('PM_scRNA_atlas', 'data', 'scanvi_models', 'main_model')

# read obj.h5ad
odata = ad.read_h5ad(os.path.join('scRNA', 'adata_tumor.h5ad'))

# create raw counts object
adata = ad.AnnData(odata.raw.X, obs=odata.obs, var=odata.raw.var)
adata.layers['counts'] = adata.X.copy()

sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=n_features,
        layer="counts",
        batch_key=batch,
        subset=True,
    )

# load SCANVI model
scanvi_model = scvi.model.SCANVI.load(scanvi_model_path, adata=adata, use_gpu=USE_GPU)
scanvi_model.to_device(DEVICE)

SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

SCANVI_MDE_KEY = "X_scANVI_MDE"
adata.obsm[SCANVI_MDE_KEY] = mde(adata.obsm[SCANVI_LATENT_KEY])

sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata)

sc.pl.umap(adata, color='sampleID', palette=PATIENT_PALETTE, show=False, save=f'_FIGURE_1C_tumor_sampleID_umap.png')
sc.pl.umap(adata, color=CELLTYPE_COL, palette=TUMOR_PALETTE, show=False, save=f'_FIGURE_1C_tumor_celltype_umap.png')
