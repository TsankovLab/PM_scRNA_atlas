# Stromal compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/stromal/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))

setwd (projdir)
source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Load scS-score
scs_sample_avg = read.csv ('../../PM_scRNA_atlas/data/scs_score_per_sample.csv', row.names=1)

# Load Seurat object
srt_tumor = readRDS ('../srt_tumor.rds')
srt = srt_tumor[, srt_tumor$celltype_simplified %in% c('Endothelial','Fibroblasts','SmoothMuscle')]

batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraphKnn = paste0 (paste(batch,collapse='_'),'_harmony_knn')
reductionGraphSnn = paste0 (paste(batch,collapse='_'),'_harmony_snn')

# Compute variable features using scran pkg
nfeat=3000
sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
rowData=rownames(srt)) 
sce = modelGeneVar(sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
vf = getTopHVGs(sce, n=nfeat)
VariableFeatures (srt) = vf

# Process merged data
srt = NormalizeData (object = srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt = ScaleData (srt, features = VariableFeatures (object=srt))
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
  
# Run Harmony
srt = srt %>% 
RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
RunUMAP (reduction = reductionSave, dims = 1:15, reduction.name = reductionName, reduction.key=reductionKey)

# Run denovo clustering on non-adjusted reductions
srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=c(reductionGraphKnn,reductionGraphSnn))


### FIGURE 3A / S3A - celltype umap ####
pdf ('Plots/celltypes_umap.pdf', 5,width = 6)
DimPlot (srt, group.by = 'celltype', reduction = reductionName, cols=palette_stroma)
DimPlot (srt, group.by = 'sampleID2', reduction = reductionName) + 
scale_color_manual (values=palette_sample)
dev.off()

### Make stacked barplot of sample usage of cNMF to put aside of heatmap ####
cp = cellComp (
srt, 
metaGroups = c('celltype','sampleID2'), 
plot_as = 'bar',
pal = palette_sample,
ptable_factor = 1,
prop=T) +
theme_minimal()
#cp$data$cNMF__r_max = factor (cp$data$cNMF__r_max, levels = names (cnmf_spectra_filtered[row_order(hm)])) 
pdf (paste0(projdir, 'Plots/FIGURE_3B_sample_abundance_celltype_stacked_barplot.pdf'))
cp
dev.off()

### FIGURE 3B - celltype dotplot ####
top_markers = c('PECAM1','CLDN5','GJA4','PLVAP','COL4A1','COL4A2','ACKR1','VWF','SELE','TFF3','CCL21','PDPN','COL1A1','COL1A2','DCN','ACTA2','MYL9','MYH11')
srt$celltype  = factor (srt$celltype, levels = c('Artery','PLVAP','Vein','LEC','Fibroblasts','SmoothMuscle'))
dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = factor (top_markers, levels = top_markers),
  x = 'celltype', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  swap_axes =T,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression2) + gtheme_italic
pdf ('Plots/FIGURE_3B_top_markers_celltype_expression.pdf', height=2.8, width=6.5)
dp
dev.off()  

### subset endothelial ####
srt_endo = srt[,srt$celltype %in% c('Artery','PLVAP','Vein')]
cnmf_spectra_unique_comb = as.list (read_excel( "../../data/cnmf_per_compartment.xlsx", sheet = "Ems_20"))

srt_endo = ModScoreCor (
        seurat_obj = srt_endo, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Ems', outdir = paste0(projdir,'Plots/'))

### FIGURE S3C - pairwise CN correlations across sample and cells ####
# Run correlation across samples 
ccomp_df = srt_endo@meta.data[,names(cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_endo@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
ccomp_df = cbind (ccomp_df, sarc = scs_sample_avg[rownames(ccomp_df),, drop=F])
ccomp_df = na.omit (ccomp_df)
cor_mat = cor (ccomp_df, method = 'spearman')
scs_score_hm = cor_mat[,ncol(cor_mat)]
scs_score_hm = scs_score_hm[-length(scs_score_hm)]
cor_mat = cor_mat[-ncol(cor_mat), -ncol(cor_mat)]

# cluster using Heatmap
hm = draw (Heatmap (cor_mat, 
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  col = palette_module_correlation,
  clustering_distance_columns = 'euclidean'))

# Run correlation across cells 
ccomp_df = srt_endo@meta.data[,names (cnmf_spectra_unique_comb)]
cor_mat_cells = cor (ccomp_df, method = 'spearman')
cor_mat_cells = cor_mat_cells[row_order(hm), row_order(hm)]
cor_mat_combined = cor_mat[row_order(hm), row_order(hm)]
cor_mat_combined[upper.tri(cor_mat_combined, diag = TRUE)] = cor_mat_cells [upper.tri(cor_mat_cells, diag = TRUE)]

# Add track showing corelation of each module to sSC-score
ha = HeatmapAnnotation (S_score = scs_score_hm[row_order(hm)], col = list (S_score = palette_sample2_fun), which='row')
hm_combined = draw (Heatmap (cor_mat_combined, 
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  clustering_distance_rows='pearson' ,
  clustering_distance_columns = 'pearson', 
  left_annotation = ha,
  col = palette_module_correlation,
  border=T,

  ))
pdf (paste0('Plots/FIGURE_S3C_CN_pairwise_across_cells_and_samples_heatmap.pdf'),5.5,4.5)
hm_combined
dev.off()


# FIGURE S3B - Make dotplot of top markers for each nmf ####
cnmf_spectra_ordered = cnmf_spectra_unique_comb[row_order(hm)]
marker_genes = unlist (lapply (cnmf_spectra_ordered, function(x) head (x, 5)))
srt_endo$Ems_r_max = factor (srt_endo$Ems_r_max, levels = rev(names(cnmf_spectra_unique_comb)[row_order(hm)]))
dotp = geneDot (
seurat_obj=srt_endo,
#gene = top_tfs2, 
gene = factor (marker_genes, levels = marker_genes),
x = 'Ems_r_max', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
swap_axes=T,
plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_S3B_cNMF_markers_expression.pdf', width=8, height =3)
dotp
dev.off()

srt_endo$celltype = factor (srt_endo$celltype, levels  = c('Artery','PLVAP','Vein'))


# FIGURE 3H - VEGFA receptors ####
dotp = geneDot (
seurat_obj=srt_endo,
#gene = top_tfs2, 
gene = c('KDR','FLT4'),
x = 'celltype', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
swap_axes=T,
plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_4H_VEGFA_receptors.pdf', width=2.5, height =2)
dotp
dev.off()



############################################################################
# Use PLVAP+ markers to compare adult vs fetal lung endothelial dataset ####
############################################################################

# Load He dataset ####
Convert("fetal_lung_endothelial_h5seurat")
srt_fetal = LoadH5Seurat("fetal_lung_endothelial_h5seurat")
srt_fetal$celltype = srt_fetal$new_celltype
head (srt_fetal)

# Load Travaglini dataset ####
srt_adult = readRDS("travaglini_lung_atlas.rds")
srt_adult = UpdateSeuratObject (srt_adult)
srt_adult = srt_adult[, srt_adult$free_annotation %in% c('Artery','Bronchial Vessel 1','Bronchial Vessel 2','Capillary','Capillary Aerocyte','Capillary Intermediate 1','Capillary Intermediate 2','Vein','Vascular Smooth Muscle')]
srt_adult$sampleID = srt_adult$sample
srt_adult$status = 'adult'

srt_adult$celltype = srt_adult$free_annotation
srt_adult$project = 'adult_lung'
srt_fetal$project = 'fetal_lung'

# Merge datasets
intersect_features = intersect (rownames (srt_fetal), rownames (srt_adult))
srt_fetal_adult = merge (srt_fetal[intersect_features,srt_fetal$dissection == 'Distal'], srt_adult[intersect_features,srt_adult$location == 'distal'])
srt_fetal_adult = srt_fetal_adult[,!grepl ('erythr',srt_fetal_adult$celltype)]
srt_fetal_adult = srt_fetal_adult[,!grepl ('reticu',srt_fetal_adult$celltype)]
low_celltype_filter = names(table (srt_fetal_adult$celltype)[table (srt_fetal_adult$celltype) > 200])
srt_fetal_adult = srt_fetal_adult[,srt_fetal_adult$celltype %in% low_celltype_filter]

# FIGURE 3E - make boxplot of fetal module score ####

# Find DEG in PLVAP EC 
Idents(srt_endo) = 'celltype'
degClusters = FindAllMarkers (srt_endo, max.cells.per.ident = 1000, min.pct = .1, logfc.threshold = .25, verbose = T)

fetal_endo_PM_deg_genes = degClusters[fetal_endo_MPM_deg$cluster == 'PLVAP',]
fetal_endo_PM_deg_genes = fetal_endo_PM_deg_genes[fetal_endo_PM_deg_genes$avg_log2FC > 0,]
fetal_endo_PM_deg_genes = fetal_endo_PM_deg_genes[fetal_endo_PM_deg_genes$p_val_adj < 0.05, ]
fetal_endo_PM_deg_genes = head (fetal_endo_PM_deg_genes[order (fetal_endo_PM_deg_genes$p_val_adj), 'gene'],15)

srt_fetal_adult  = ModScoreCor (
        seurat_obj = srt_fetal_adult, 
        geneset_list = list(endo_fetal = fetal_endo_PM_deg_genes), 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'fetal', outdir = paste0(projdir,'Plots/'))

ccomp = srt_fetal_adult@meta.data
box = ggplot(ccomp, aes (x= celltype, y= endo_fetal)) + 
  ggtitle (paste('Fetal vs adult cohort')) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
    plot.title = element_text(size = 8)) + 
    #geom_jitter(alpha=0.5, size =5, position = position_jitter(width = 0.4)) +
  geom_violin (width=1, aes(fill = project),alpha = 0.7, lwd=.2) +
  geom_boxplot (aes (fill = project), color='black', outlier.colour="black", outlier.shape=16,
           outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2,width=.3) +
  scale_fill_manual (values = pallette_fetal_vs_adult) +
  gtheme
fetal_score = box$data %>% group_by (celltype) %>% summarise_at(vars(endo_fetal),
                ~summary(.)[[str_to_title('median')]])
box$data$celltype = factor (box$data$celltype, levels = arrange (fetal_score, -endo_fetal)$celltype)  

pdf ('Plots/FIGURE_3E_fetal_marker_genes_score_fetal_vs_krasnow_vlnplot.pdf', height=3.6, width=5.5)
box
dev.off()

#### FIGURE 3C - correlation heatmap with Travaglini ####
srt_endo$status = 'tumor'
srt_endo$location = 'pleura'
srt_endo$free_annotation = srt_endo$celltype
shared_feat = rownames(srt_endo)[rownames(srt_endo) %in% rownames(srt_adult)]

srt_adult_PM = merge (srt_endo[shared_feat,], srt_adult[shared_feat,])
srt_adult_PM = NormalizeData (object = srt_adult_PM, normalization.method = "LogNormalize", scale.factor = 10000)

srt_adult_PM$celltype = gsub (' ','_',srt_adult_PM$celltype)

nfeat = 2000
srt_adult_PM_distal = srt_adult_PM[,srt_adult_PM$location %in% c('pleura','distal')]

# Run seurat integration ####
DefaultAssay(srt_adult_PM_distal) = 'RNA'
nfeat_int = 3000
k.weight = 30
metaGroupNames = c('location')

srtL = SplitObject (srt_adult_PM_distal, split.by = metaGroupNames[1])

# normalize and identify variable features for each dataset independently
srtL = lapply(X = srtL, FUN = function(x) {
    DefaultAssay(x) = 'RNA'
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = "vst", nfeat = nfeat_int)
})

# select features that are repeatedly variable across datasets for integration
features = SelectIntegrationFeatures (object.list = srtL, nfeatures = nfeat_int)

anchors = FindIntegrationAnchors (object.list = srtL, anchor.features = features)
srt_adult_PM_distal = IntegrateData (anchorset = anchors, k.weight=k.weight)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(srt_adult_PM_distal) = "RNA"

library (scran)
sce = SingleCellExperiment (list(counts=srt_adult_PM_distal@assays$RNA@counts, logcounts = srt_adult_PM_distal@assays$RNA@data),
rowData=rownames(srt_adult_PM_distal)) 
sce = modelGeneVar(sce)

# remove batchy genes
batchy_genes = c('RPL','RPS','MT-','MALAT')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
vf = getTopHVGs(sce, n=nfeat)
    
gcdata.avg = AverageExpression (srt_adult_PM_distal, group.by = 'celltype',return.seurat = T) # ,assays = 'integrated')
VariableFeatures(gcdata.avg) = vf
DefaultAssay (gcdata.avg) = 'integrated'
gcdata.avg <- ScaleData(gcdata.avg)
#gcdata.avg = gcdata.avg[,rownames(gcdata.avg@meta.data) != 'VascularSmoothMuscle']

mtx = as.matrix (GetAssayData(gcdata.avg, assay = 'integrated', slot='scale.data'))
library (circlize)

corr = cor(mtx,method = 'spearman')

srt_adult_PM_distal$project[is.na(srt_adult_PM_distal$project)] = 'PM'
ha = HeatmapAnnotation (' ' = srt_adult_PM_distal$project[match (rownames(corr),srt_adult_PM_distal$celltype)], col=list(' ' = pallette_fetal_vs_adult))
ha2 = rowAnnotation (' ' = srt_adult_PM_distal$project[match (rownames(corr),srt_adult_PM_distal$celltype)], col=list(' ' = pallette_fetal_vs_adult))
rownames (corr) = gsub ('_',' ', rownames(corr))
rownames (corr) = gsub ('PLVAP','PLVAP+ EC', rownames(corr))
colnames (corr) = gsub ('_',' ', colnames(corr))
colnames (corr) = gsub ('PLVAP','PLVAP+ EC', colnames(corr))
pdf ('Plots/FIGURE_3C_cor_mat_travaglini_heatmap_top20.pdf',width = 7,4.6)
Heatmap(corr, top_annotation = ha,
  col = palette_module_correlation_fun, right_annotation = ha2,
  column_names_rot = 70,
  border=T)
dev.off()



### subset fibroblasts ####
library (readxl)
library (circlize)

srt_fib = srt[,srt$celltype %in% c('Fibroblasts')]
cnmf_spectra_unique_comb = as.list (read_excel( "../../cnmf_per_compartment.xlsx", sheet = "Fms_20"))

srt_fib = ModScoreCor (
        seurat_obj = srt_fib, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Fms', outdir = paste0(projdir,'Plots/'))


### FIGURE S3E - pairwise CN correlations across sample and cells ####
# Run correlation across samples 
scs_sample_avg = read.csv ('../../PM_scRNA_atlas/data/scs_score_per_sample.csv', row.names=1)
ccomp_df = srt_fib@meta.data[,names(cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_fib@meta.data[,'sampleID4',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
ccomp_df = cbind (ccomp_df, scs = scs_sample_avg[rownames(ccomp_df),])
cor_mat = cor (ccomp_df, method = 'spearman')
scs_sample_avg = cor_mat[,'scs']
scs_sample_avg = scs_sample_avg[-length(scs_sample_avg)]
cor_mat = cor_mat[-ncol(cor_mat), -ncol(cor_mat)]

# cluster using Heatmap
hm = draw (Heatmap (cor_mat, 
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean'))

# Run correlation across cells 
ccomp_df = srt_fib@meta.data[,names (cnmf_spectra_unique_comb)]
cor_mat_cells = cor (ccomp_df, method = 'spearman')
cor_mat_cells = cor_mat_cells[row_order(hm), row_order(hm)]
cor_mat_combined = cor_mat[row_order(hm), row_order(hm)]
cor_mat_combined[upper.tri(cor_mat_combined, diag = TRUE)] = cor_mat_cells [upper.tri(cor_mat_cells, diag = TRUE)]

col_fun = colorRamp2(c(-1, 0, 1), rev(c(palette_sample[3], palette_sample[length(palette_sample)/2], palette_sample[length(palette_sample)-3])))
palette_module_correlation = paletteer::paletteer_c("pals::kovesi.diverging_bwr_40_95_c42",100)

# Add track showing corelation of each module to scS-score
ha = HeatmapAnnotation (scS_score = scs_sample_avg[row_order(hm)], col = list (scS_score = col_fun), which='row')
hm_combined = draw (Heatmap (cor_mat_combined, 
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  clustering_distance_rows='pearson' ,
  clustering_distance_columns = 'pearson', 
  left_annotation = ha,
  col = palette_module_correlation,
  border=T))

pdf (paste0('Plots/FIGURE_S3E_Fms_pairwise_across_cells_and_samples_heatmap.pdf'),5.5,4.5)
hm_combined
dev.off()


### FIGURE S3D - Make dotplot of top markers for each nmf ####
cnmf_spectra_unique_comb_ordered = cnmf_spectra_unique_comb[row_order(hm)]
marker_genes = unlist (lapply (cnmf_spectra_unique_comb_ordered, function(x) head (x, 5)))
marker_genes[marker_genes == 'DUSP6'] = 'IGFBP2'
marker_genes[marker_genes == 'C9orf16'] = 'IGFBP6'
marker_genes[marker_genes == 'SPTBN1'] = 'MFAP5'
#marker_genes[marker_genes == 'COL8A1'] = 'CXCL14'
marker_genes[marker_genes == 'MALAT1'] = 'COL6A2'
#marker_genes[marker_genes == 'PNISR'] = 'MMP2'

srt_fib$Fms_r_max = factor (srt_fib$Fms_r_max, levels = rev(names(cnmf_spectra_unique_comb_ordered)))

dotp = geneDot (
seurat_obj=srt_fib,
#gene = top_tfs2, 
gene = factor (marker_genes, levels = marker_genes),
x = srt_fib$sampleID4, # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
y = 'Fms_r_max',
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
plotcol = gene_expression_palette)

pdf ('Plots/FIGURE_S3D_cNMF_markers_expression.pdf', width=7, height = 3)
dotp
dev.off()


# Load Travaglini dataset ####
srt_fib$status = 'tumor'
srt_fib$location = 'pleura'
srt_adult = readRDS("travaglini_lung_atlas.rds")
srt_adult = UpdateSeuratObject (srt_adult)
srt_adult = srt_adult[, srt_adult$free_annotation %in% c('Alveolar Fibroblast','Adventitial Fibroblast','Fibromyocyte','Myofibroblast','Pericyte','Lipofibroblast')]
srt_adult$sampleID4 = srt_adult$sample
srt_adult$status = 'adult'
srt_adult$celltype = srt_adult$free_annotation

srt_fib$celltype = srt_fib$Fms_r_max
shared_feat = rownames(srt_fib)[rownames(srt_fib) %in% rownames(srt_adult)]

srt_merged = merge (srt_fib[shared_feat,], srt_adult[shared_feat,])
srt_merged = NormalizeData (object = srt_merged, normalization.method = "LogNormalize", scale.factor = 10000)  

srt_merged$celltype = gsub (' ','_',srt_merged$celltype)

nfeat = 2000
srt_merged_distal = srt_merged[,srt_merged$location %in% c('pleura','distal')]

# Run seurat integration ####
DefaultAssay(srt_merged_distal) = 'RNA'
nfeat_int = 3000
k.weight = 30

srtL = SplitObject (srt_merged_distal, split.by = 'location')
# normalize and identify variable features for each dataset independently
srtL = lapply(X = srtL, FUN = function(x) {
    DefaultAssay(x) = 'RNA'
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = "vst", nfeat = nfeat_int)
})

# select features that are repeatedly variable across datasets for integration
features = SelectIntegrationFeatures (object.list = srtL, nfeatures = nfeat_int)

anchors = FindIntegrationAnchors (object.list = srtL, anchor.features = features)
srt_merged_distal = IntegrateData (anchorset = anchors, k.weight=k.weight)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(srt_merged_distal) = "RNA"

library (scran)
set.seed(1234)
sce = SingleCellExperiment (list(counts=srt_merged_distal@assays$RNA@counts, logcounts = srt_merged_distal@assays$RNA@data),
rowData=rownames(srt_merged_distal)) 
sce = modelGeneVar(sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-','MALAT')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
vf = getTopHVGs(sce, n=nfeat)
    
gcdata.avg = AverageExpression (srt_merged_distal, group.by = 'celltype',return.seurat = T) # ,assays = 'integrated')
VariableFeatures(gcdata.avg) = vf
DefaultAssay (gcdata.avg) = 'integrated'
gcdata.avg <- ScaleData(gcdata.avg)
#gcdata.avg = gcdata.avg[,rownames(gcdata.avg@meta.data) != 'VascularSmoothMuscle']

mtx<-as.matrix(GetAssayData(gcdata.avg, assay = 'integrated', slot='scale.data'))
corr<-cor(mtx,method = 'spearman')
col_fun = colorRamp2(c(-1, 0, 1), c(palette_module_correlation[1], palette_module_correlation[length(palette_module_correlation)/2], palette_module_correlation[length(palette_module_correlation)]))
pdf ('Plots/FIGURE_S3F_cor_mat_travaglini_fibroblasts_heatmap_top20.pdf')
Heatmap(corr,col = col_fun)
dev.off()



### FIGURE 3F - SCENIC only on MPM endothelial ####
# Import results
auc_mtx <- read.csv('../../PM_scRNA_atlas/data/SCENIC_Endothelial_auc_mtx.csv')
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = paste0('SCENIC_',colnames(auc_mtx))
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

bc_scenic = rownames (auc_mtx)
bc_scenic = gsub ('\\.','-', bc_scenic)

celltypes_scenic = srt_endo$celltype

# Generate heatmap of TFs ####
motifs_table = read.csv ('../../PM_scRNA_atlas/data/SCENIC_Endothelial_motifs.csv', skip=2) 

pValThreshold = 1e-2
auc_mtx_avg = aggregate (auc_mtx, by=list(unname(celltypes_scenic)), mean)
rownames (auc_mtx_avg) = auc_mtx_avg[,1]
auc_mtx_avg = auc_mtx_avg[,-1]
colnames (auc_mtx_avg) = sub ('SCENIC_','', colnames(auc_mtx_avg))
filter_TF = unique (motifs_table$TF[motifs_table$X.2 < pValThreshold])
auc_mtx_avg = auc_mtx_avg[,colnames(auc_mtx_avg) %in% filter_TF]
proj_celltypes = rownames(auc_mtx_avg)

hm = Heatmap (scale(auc_mtx_avg),
 col = palette_scenic, column_names_rot = 45,
# right_annotation = ha,
 column_names_gp = gpar(fontsize = 7),
 row_names_gp = gpar(fontsize = 8))
pdf ('Plots/FIGURE_S3G_auc_scores_celltypes_heatmap.pdf', width=7, height=2)
hm
dev.off()

# Make module score of ETS1 and MEF2C regulons ####
tfs = c('MEF2C','ETS1')
motifs_table_flt = lapply (tfs, function(x) unlist (sapply (motifs_table[motifs_table$TF == x, 'X.6'], function(y) unlist(strsplit (y, ','))),use.names=F))
motifs_table_flt = lapply (motifs_table_flt, function(x) x[grep ('\\(',x)])
motifs_table_flt = lapply (motifs_table_flt, function(x) gsub (' \\(','',x))
motifs_table_flt = lapply (motifs_table_flt, function(x) gsub ("'",'',x))
names (motifs_table_flt) = tfs

intersect_features = intersect (rownames (srt_fetal_adult), rownames (srt_endo))
srt_fetal_adult_mpm = merge (srt_fetal_adult[intersect_features,], srt_endo[intersect_features,])
srt_fetal_adult_mpm = NormalizeData (object = srt_fetal_adult_mpm, normalization.method = "LogNormalize", scale.factor = 10000)
srt_fetal_adult_mpm$project[is.na(srt_fetal_adult_mpm$project)] = 'PM'
srt_fetal_adult_mpm$celltype_location = paste0(srt_fetal_adult_mpm$celltype,'_',srt_fetal_adult_mpm$project)

srt_fetal_adult_mpm = ModScoreCor (
        seurat_obj = srt_fetal_adult_mpm, 
        geneset_list = motifs_table_flt, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'fetal', outdir = paste0(projdir,'Plots/'))


ccomp = srt_fetal_adult_mpm@meta.data
table (ccomp$celltype_location, useNA = 'always')
ccomp = ccomp[ccomp$celltype_location %in% names(table (ccomp$celltype_location))[table (ccomp$celltype_location) > 200],]
box = lapply (tfs, function(x) ggplot(ccomp, aes_string (x= 'celltype_location', y= x)) + 
  ggtitle (x) +
  geom_violin (width=1, aes(fill = project),alpha = 0.7, lwd=.2) +
  geom_boxplot (aes (fill = project), color='black', outlier.colour="black", outlier.shape=16,
           outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2,width=.3) +
  scale_fill_manual (values = pallette_fetal_vs_adult) +
  gtheme) 
box = lapply (seq_along (tfs), function(x) 
    {
     ct_order = box[[x]]$data %>% group_by (celltype_location) %>% summarise_at(vars(tfs[x]),
                ~summary(.)[[str_to_title('median')]])
      box[[x]]$data$celltype_location = factor (box[[x]]$data$celltype_location, levels = ct_order$celltype_location[order(-ct_order[,2])])
      labels = levels (box[[x]]$data$celltype_location)
      labels = gsub ('_', ' ', labels)
      labels = gsub (' PM','', labels)
      labels = gsub (' fetal lung','', labels)
      labels = gsub (' adult lung','', labels)
      labels = gsub ('PLVAP','PLVAP+ EC', labels)
      box[[x]] + scale_x_discrete(labels= labels)
    })


pdf ('Plots/FIGURE_3F_SCENIC_ETS1_MEF2C_regulons_modulescores_vlnplots.pdf',5,4)
box
dev.off()
