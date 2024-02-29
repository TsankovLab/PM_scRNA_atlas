use UGER
conda activate scrnatools 
R

set.seed(1234)
options(warn = 1)

projdir = 'scRNA/pbmc/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Load scS-score
scs_sample_avg = read.csv ('../../PM_scRNA_atlas/data/scs_score_per_sample.csv', row.names=1)

# Load Seurat object

srt = readRDS ('../srt_pbmc.rds')

# Normalize, scale and compute UMAP after correcting for sampleID batch with harmony
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


### FIGURE 1C ####
pq1 = DimPlot(srt, reduction = reductionName, group.by = "celltype_simplified", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt), 'cells')) +
scale_color_manual(values = palette_celltype_simplified)
pq2 = DimPlot(srt, reduction = reductionName, group.by = "sampleID", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt), 'cells')) +
scale_color_manual(values = palette_sample)

png (paste0("Plots/FIGURE_1C_celltypes_samples_umap.png"), width = 2800, height = 1000, res=300)
pq1 + pq2
dev.off()

### FIGURE S1A ####
ccomp_df = as.data.frame (table(srt$sampleID))
srt$nFeature_RNAL = log2 (srt$nFeature_RNA)
srt$nCount_RNAL = log2 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=F, group.by = 'sampleID',pt.size = 0, ncol = 3) 
vln_p = lapply (vln_p, function(x) x +
scale_fill_manual(values = palette_sample))

png (paste0("Plots/FIGURE_S1A_QC_nFeat_nCount_m.percent_vlnPlot.png"), 3300, 1000, res=300)
print (vln_p[[1]] | vln_p[[2]] | vln_p[[3]])
dev.off()

### FIGURE S1D ####
RNA.markers = c('CD3E','CD8A','MKI67','NKG7','CD79A','JCHAIN','S100A8','FCGR3A','HLA-DQA1','LILRA4')
DefaultAssay (srt) = 'RNA'
fp = FeaturePlot(srt, features = RNA.markers,  reduction = reductionName, cols = palette_feature_RNA, ncol = 4) & 
theme(plot.title = element_text(size = 10, face='italic')) & NoLegend() & NoAxes()

png (paste0("Plots/FIGURE_S1D_feature_plots_RNA.png"), width = 1800, height = 1000, res=300)
fp
dev.off()

RNA.markers=c('CD3D','CD3E','CD8A','CD8B','CCL5','ITGB1','IL7R',
              'FOXP3','CTLA4','LEF1','CCR7','TCF7','SELL','LTB',
               'MKI67','CD27','TRDV2',
              'NCR3','KLRB1','GATA3','XCL1','XCL2','NKG7','GNLY','PRF1','GZMH','FGFBP2','KLRC1','CMC1',
              'CD79A','MS4A1','CD37','TCL1A','IL4R','FCER2',
              'JCHAIN','IGHA1','IGKC',
              'S100A8','S100A9','FCN1','AIF1','LST1','FCGR3A','IFITM3',
              'HLA-DQA1','HLA-DQB1','CLEC9A',
               'FCER1A','CD1C','CLEC10A','AXL','SIGLEC6','IRF7','LILRA4','PLD4',
              'SOX4','CDK6','HBA1','HBA2','ALAS2','PPBP',
               'NRGN','CAVIN2')
label_order = c('CD4 CTL','Treg','CD4 Naive','CD8 Naive','CD4 TCM','CD8 TCM',
  'CD4 TEM','CD8 TEM','CD4 Proliferating','CD8 Proliferating','dnT','gdT','MAIT',
  'ILC','NK','NK_CD56bright','NK Proliferating','B naive','B intermediate','B memory',
  'CD14 Mono','CD16 Mono','cDC1','cDC2','ASDC','pDC','HSPC','Eryth','Platelet',
  'Plasmablast','Doublet')
srt$predicted.celltype.l2 = factor (srt$predicted.celltype.l2, levels = rev (label_order))

### FIGURE S1E ####
dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = factor (RNA.markers, levels = RNA.markers),
  x = 'predicted.celltype.l2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  swap_axes = T,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression2) +
  gtheme
  
png(paste0("Plots/FIGURE_S1E_DotPlot_topRNAmarkers_predicted.celltype.l2.png"), width = 4500, height = 2200, res = 300)
dp 
dev.off()

### FIGURE S1F ####
DefaultAssay(srt)='ADT'
srt = NormalizeData(object = srt, assay = "ADT", normalization.method = "CLR")
srt = ScaleData(srt,assay = "ADT")

ADT.markers = c('AC-CD3','AC-CD4','AC-CD8','AC-CD19','AC-CD33','AC-CD14','AC-CD16','AC-CD56')
DefaultAssay (srt) = 'ADT'
fp = FeaturePlot(srt, features = ADT.markers,  reduction = reductionName, cols = palette_feature_protein, ncol = 4) & 
theme(plot.title = element_text(size = 10)) & NoLegend() & NoAxes()

png (paste0("Plots/FIGURE_S1F_feature_plots_ADT.png"), width = 1800, height = 1000, res=300)
fp
dev.off()

ADT.markers=c('AC-CD3','AC-CD4','AC-CD8','AC-CD5','AC-CD2','AC-CD57','AC-CD25','AC-CD27',
              'AC-CD7','AC-CD62L','AC-CD26','AC-CD103',
              'AC-CD161','AC-CD127','AC-GPR56','AC-CD56','AC-CD94','AC-IgD','AC-CD21',
              'AC-CD19','AC-CD22','AC-CD1c','AC-CD73','AC-CD38','AC-CD33',
              'AC-CD14','AC-CD16','AC-CD54','AC-CD141','AC-CD11c','AC-CD123',
              'AC-CD303','AC-CD62P','AC-CD47','AC-CD71','AC-CD41','AC-CD42b')

### FIGURE S1G ####
dp = geneDot (
  seurat_obj = srt,
  assay= 'ADT',
  #gene = top_tfs2, 
  gene = factor (ADT.markers, levels = ADT.markers),
  x = 'predicted.celltype.l2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  x_name ='samples',
  y_name = 'celltype',
  swap_axes=T,
  plotcol = palette_protein_expression) +
  gtheme
  
png(paste0("Plots/FIGURE_S1G_DotPlot_topADTmarkers_predicted.celltype.l2.png"), width = 3500, height = 2200, res = 300)
dp 
dev.off()

### UMAP of cell types ####
srt$celltype = srt$predicted.celltype.l2
srt$celltype[grep ('CD4', srt$celltype)] = 'CD4'
srt$celltype[grep ('CD8', srt$celltype)] = 'CD8'
srt$celltype[srt$celltype %in% c('B memory','B naive','CD14 Mono','Eryth','dnT','Platelet','gdT')] = NA
srt = srt[,!is.na (srt$celltype)]
p = DimPlot (srt, group.by= 'celltype', cols = pallette_pbmc_celltype) + theme_void()
p[[1]]$layers[[1]]$aes_params$alpha =  .8

p2 = DimPlot (srt, group.by= 'sampleID', cols = palette_sample) + theme_void() + scale_fill_manual(values = palette_sample)
p[[1]]$layers[[1]]$aes_params$alpha =  .8

png (paste0(projdir, 'Plots/FIGURE_5L_celltypes_umap.png'), width = 1400, height=1000, res = 300)
p
dev.off()

png (paste0(projdir, 'Plots/FIGURE_S5F_sample_umap.png'), width = 1400, height=1000, res = 300)
p2
dev.off()

DefaultAssay (srt) = 'RNA'
markers  = c('CD3D','CD4','CD8A','FOXP3','NKG7','XCL1','MKI67')
fp = FeaturePlot(srt, features = markers,  reduction = reductionName, cols = palette_feature_RNA, ncol = length(markers)) & 
theme(plot.title = element_text(size = 10, face='italic')) & NoLegend() & NoAxes()

DefaultAssay (srt) = 'ADT'
markers_adt = paste0('AC-',c('CD3','CD4','CD8','CD25','CD56','CD127','CD161'))
fp2 = FeaturePlot(srt, features = markers_adt,  reduction = reductionName, cols = palette_feature_protein, ncol = length(markers_adt)) & 
theme(plot.title = element_text(size = 10)) & NoLegend() & NoAxes()

png (paste0(projdir, 'Plots/FIGURE_5M_markers_featureplots.png'), height=1300, width=2900, res=300)
wrap_plots (fp, fp2, nrow=2)
#pc / pc2
dev.off()
  
### Subset to T-NK cells
srt_tnk = srt[,srt$predicted.celltype.l1 %in% c('CD4 T','CD8 T','NK')]

# Normalize, scale and compute UMAP after correcting for sampleID batch with harmony ####
batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraphKnn = paste0 (paste(batch,collapse='_'),'_harmony_knn')
reductionGraphSnn = paste0 (paste(batch,collapse='_'),'_harmony_snn')

# Compute variable features using scran pkg
nfeat=3000
DefaultAssay(srt_tnk) = 'RNA'
sce = SingleCellExperiment (list(counts=srt_tnk@assays$RNA@counts, logcounts = srt_tnk@assays$RNA@data),
rowData=rownames(srt_tnk)) 
sce = modelGeneVar(sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
vf = getTopHVGs(sce, n=nfeat)
VariableFeatures (srt_tnk) = vf

# Process merged data
srt_tnk = NormalizeData (object = srt_tnk, normalization.method = "LogNormalize", scale.factor = 10000)
srt_tnk = ScaleData (srt_tnk, features = VariableFeatures (object=srt_tnk))
srt_tnk = RunPCA (srt_tnk, features = VariableFeatures (object = srt_tnk), npcs = ifelse(ncol(srt_tnk) <= 30,ncol(srt_tnk)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
  
# Run Harmony
srt_tnk = srt_tnk %>% 
RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
RunUMAP (reduction = reductionSave, dims = 1:15, reduction.name = reductionName, reduction.key=reductionKey)

# Run denovo clustering on non-adjusted reductions
srt_tnk = FindNeighbors (object = srt_tnk, reduction = reductionSave, dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=c(reductionGraphKnn,reductionGraphSnn))

srt_tnk$celltype = srt_tnk$predicted.celltype.l2
srt_tnk$celltype[grep ('CD4', srt_tnk$celltype)] = 'CD4'
srt_tnk$celltype[grep ('CD8', srt_tnk$celltype)] = 'CD8'
srt_tnk$celltype[srt_tnk$celltype %in% c('B memory','B naive','CD14 Mono','Eryth','dnT','Platelet','gdT')] = NA

### FIGURE 5E S5F ####
pq1 = DimPlot(srt_tnk, reduction = reductionName, group.by = "celltype", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt_tnk), 'cells')) +
scale_color_manual(values = pallette_pbmc_celltype)
pq2 = DimPlot(srt_tnk, reduction = reductionName, group.by = "sampleID", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt_tnk), 'cells')) +
scale_color_manual(values = palette_sample)

png (paste0("Plots/FIGURE_5L_celltypes_samples_umap.png"), width = 2800, height = 1000, res=300)
pq1 + pq2
dev.off()


#####################
### TCR ANALYSIS ####
#####################
#### Load PBMC TCR-sequencing from biopsy samples ####
meta_pbmc = srt@meta.data[srt$batch == 'batch2',]
meta_pbmc$sampleID = as.character (meta_pbmc$sampleID)

# Unzip TCR contig files before this
data.path_pbmc1 = 'GSM7114958_PBMC_Processed_TCR-2-1_TCR_all_contig_annotations.csv'
data.path_pbmc2 = 'GSM7114959_PBMC_Processed_TCR-2-2_TCR_all_contig_annotations.csv'
vdj.dirs = c(data.path_pbmc1, data.path_pbmc2)

# Add hashing pools to seurat metadata ####
meta_pbmc$hash_pool = gsub ('_batch2','', meta_pbmc$pool_batch)

# Check barcodes per pool are unique ####
lapply (split (meta_pbmc$barcode,meta_pbmc$hash_pool), function(x) length(x) == length(unique(x)))

# Split TCR barcodes by pool map barcodes and split by sample ####
tcrL_pbmc = lapply (vdj.dirs, function(x) read.csv(x))
names (tcrL_pbmc) = unique (meta_pbmc$hash_pool)
tcrL_pbmc = lapply (names(tcrL_pbmc), function(x)
  {
  tcr_hashed = tcrL_pbmc[[x]][tcrL_pbmc[[x]]$barcode %in% meta_pbmc$barcode[meta_pbmc$hash_pool == x],]
  tcr_hashed$sampleID = meta_pbmc$sampleID[match (tcr_hashed$barcode, meta_pbmc$barcode)]
  tcr_hashed
  })
tcrL_pbmc = do.call (rbind, tcrL_pbmc)
tcrL_pbmc = split (tcrL_pbmc, tcrL_pbmc$sampleID)
names (tcrL_pbmc) = paste0(names(tcrL_pbmc) ,'_pbmc')

# Load PBMC TCR-sequencing from resection samples ####
data.path1 = '30-472307723_2021-01-05-2-1_TCR_all_contig_annotations.csv'
data.path2 = '30-472307723_2021-01-05-2-2_TCR_all_contig_annotations.csv'
vdj.dirs = c(data.path1, data.path2)
meta_pbmc_batch1 = srt@meta.data[srt$batch == 'batch1',]
meta_pbmc_batch1$sampleID = as.character (meta_pbmc_batch1$sampleID)

# Add hashing pools to seurat metadata ####
meta_pbmc_batch1$hash_pool = sapply (rownames(meta_pbmc_batch1), function(x) unlist(strsplit (x, '\\-'))[2])
meta_pbmc_batch1$hash_pool = gsub ('_batch1','',meta_pbmc_batch1$pool_batch)

# Strip barcodes to match TCR data ####
#meta_pbmc_batch1$barcode = sapply (rownames (meta_pbmc_batch1), function(x) unlist (strsplit (x, '_'))[1])
head (meta_pbmc_batch1$barcode)
tail (meta_pbmc_batch1$barcode)

# Check barcodes per pool are unique ####
lapply (split (meta_pbmc_batch1$barcode,meta_pbmc_batch1$hash_pool), function(x) length(x) == length(unique(x)))

tcrL_pbmc_batch1 = lapply (vdj.dirs, function(x) read.csv(x))
names (tcrL_pbmc_batch1) = unique (meta_pbmc_batch1$hash_pool)
tcrL_pbmc_batch1 = lapply (names(tcrL_pbmc_batch1), function(x)
  {
  tcr_hashed = tcrL_pbmc_batch1[[x]][tcrL_pbmc_batch1[[x]]$barcode %in% meta_pbmc_batch1$barcode[meta_pbmc_batch1$hash_pool == x],]
  tcr_hashed$sampleID = meta_pbmc_batch1$sampleID[match (tcr_hashed$barcode, meta_pbmc_batch1$barcode)]
  tcr_hashed
  })
tcr_df_batch1 = do.call (rbind, tcrL_pbmc_batch1)
tcr_df_batch1 = split (tcr_df_batch1, tcr_df_batch1$sampleID)
names (tcr_df_batch1) = paste0(names(tcr_df_batch1) ,'_pbmc')

contig_list = c(tcr_df_batch1, tcrL_pbmc)
combinedTCR = combineTCR (contig_list, 
                samples = names (contig_list), 
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

pdf ('Plots/clonotype_overlap.pdf',width=5,3.5)
clonalOverlap(combinedTCR, 
              cloneCall = "strict", 
              method = "morisita")
dev.off()

# Check redundancy contig list
contig_list_unique = unlist(lapply (contig_list, function(x) x$barcode[!duplicated(x$barcode)]))
length (contig_list_unique) ; length (unique(contig_list_unique))

unique_barcode = !duplicated(paste0(srt$sampleID, '_pbmc_', srt$barcode))
srt_unique = srt[,unique_barcode]
srt_unique = RenameCells(
    srt_unique,
    new.names = paste0(srt_unique$sampleID, '_pbmc_', srt_unique$barcode))


# Integrate TCR info with scRNA barcodes ####
srt_tcr = combineExpression (combinedTCR, srt_unique,
  proportion = F,
#  cloneSizes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = Inf))
  cloneSize = c(NonExpanded = 1, Small = 5, Large = Inf))
  #)

# Subset for cells with clonotype info and TNK ####
clonesize_name = setNames (c('NonExpanded','Small','Large','None'), c('NonExpanded (0 < X <= 1)','Small (1 < X <= 5)','Large (5 < X <= Inf)', 'None ( < X <= 0)'))
srt_tcr = srt_tcr[, !is.na (srt_tcr$cloneSize)]
srt_tcr$cloneSize = unname(clonesize_name[as.character(srt_tcr$cloneSize)])
srt_tcr_tnk = srt_tcr[,srt_tcr$predicted.celltype.l1 %in% c('CD4 T','CD8 T','NK')]


# FIGURE S5H - Show expanded clonotypes are found in CD8 cells ####  
metaGroupName = 'cloneSize'
cc_box = cellComp (
  seurat_obj = srt_tcr_tnk, 
  metaGroups = c('predicted.celltype.l2',metaGroupName),
  plot_as = 'bar',
  prop = FALSE,
  pal = palette_clonotype,
  subset_prop = 'Large',  
  #pal = viridis::mako(,
  facet_ncol = 15
  ) + gtheme
cc_box$data$predicted.celltype.l2 = factor (cc_box$data$predicted.celltype.l2, levels = as.character(cc_box$data$predicted.celltype.l2[order(-cc_box$data$Freq)]))

png(paste0('Plots/FIGURE_S5H_Expanded_vs_nonExpanded_barplot.png'),width=1200,height=800, res=300)
cc_box
dev.off()

### Subset to CD8 ####
srt_tcr_cd8 = srt_tcr_tnk[, srt_tcr_tnk$predicted.celltype.l1 %in% c('CD8 T')]
srt_tcr_cd8 = srt_tcr_cd8[, !srt_tcr_cd8$sampleID %in% c('P3')] # remove low counts samples

# Import T cell cNMFs
cnmf_t = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Tms_20"))

DefaultAssay (srt_tcr_cd8) = 'RNA'
srt_tcr_cd8 = ModScoreCor (
        seurat_obj = srt_tcr_cd8, 
        geneset_list = cnmf_t, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Tms_', outdir = paste0(projdir,'Plots/'))

# FIGURE 6E - Show expanded clonotypes have higher exhaustion ####  
ccomp_df = as.data.frame (srt_tcr_cd8@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm5, by=as.list(srt_tcr_cd8@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_no_text + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_6N_expanded_exhaustion.png',height = 1000,width = 600, res=300)
box2
dev.off()

### Check with MHC II module ####
ccomp_df = as.data.frame (srt_tcr_cd8@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm4, by=as.list(srt_tcr_cd8@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_no_text + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_S6J_expanded_MHCII.png',height = 1000,width = 600, res=300)
box2
dev.off()


### Check with Cytotoxic module ####
ccomp_df = as.data.frame (srt_tcr_cd8@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm2, by=as.list(srt_tcr_cd8@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_no_text + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_S6J_expanded_cytotoxic.png',height = 1000,width = 600, res=300)
box2
dev.off()














