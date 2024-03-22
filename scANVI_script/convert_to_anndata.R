library(Seurat)
library(SeuratDisk)

tumor <- readRDS('scRNA/GSE190597_srt_tumor.rds')
tumor@meta.data$sampleID = as.character(tumor@meta.data$sampleID)
SaveH5Seurat(tumor, filename = "scRNA/adata_tumor.h5Seurat", overwrite=T)
Convert("scRNA/adata_tumor.h5Seurat", dest = "h5ad", overwrite=T)

pbmc <- readRDS('scRNA/GSE190597_srt_pbmc.rds')
pbmc@meta.data$sampleID = as.character(pbmc@meta.data$sampleID)
SaveH5Seurat(pbmc, filename = "scRNA/adata_pbmc.h5Seurat", overwrite=T)
Convert("scRNA/adata_pbmc.h5Seurat", dest = "h5ad", overwrite=T)

# delete tmp files
system('rm scRNA/*.h5Seurat')
