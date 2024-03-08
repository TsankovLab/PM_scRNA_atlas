library(Seurat)
library(SeuratDisk)

tumor <- readRDS('scRNA/srt_tumor.rds')
SaveH5Seurat(tumor, filename = "scRNA/adata_tumor.h5Seurat", overwrite=T)
Convert("scRNA/adata_tumor.h5Seurat", dest = "h5ad", overwrite=T)

pbmc <- readRDS('scRNA/srt_pbmc.rds')
SaveH5Seurat(pbmc, filename = "scRNA/adata_pbmc.h5Seurat", overwrite=T)
Convert("scRNA/adata_pbmc.h5Seurat", dest = "h5ad", overwrite=T)

# delete tmp files
system('rm scRNA/*.h5Seurat')
