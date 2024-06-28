# Myeloid compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/myeloid/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

# Load scripts
source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Load scS-score
scs_sample_avg = read.csv ('../../PM_scRNA_atlas/data/scs_score_per_sample.csv', row.names=1)

# Load Seurat object
srt_tumor = readRDS ('../GSE190597_srt_tumor.rds')
srt = srt_tumor[, srt_tumor$celltype_simplified2 %in% c('Myeloid','pDC')]

batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraphKnn = paste0 (paste(batch,collapse='_'),'_harmony_knn')
reductionGraphSnn = paste0 (paste(batch,collapse='_'),'_harmony_snn')

# Compute variable features using scran pkg
nfeat=2000
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


### FIGURE 4A / S4A - celltype umap ####
dp = DimPlot (srt, group.by = 'celltype', reduction = reductionName, cols=palette_myeloid)
dp1 = DimPlot (srt, group.by = 'sampleID', reduction = reductionName) + 
    scale_color_manual (values=palette_sample)
pdf ('Plots/FIGURE_4A_S4A_celltypes_umap.pdf', 5,width = 6)
print (dp)
print (dp1)
dev.off()



#### Plot markers ####
markers = c('CD16','VCAN','APOE','IL1B', 'LILRA4','CCR7','CD1C','CLEC9A','LST1','TPSB2','VSIR','MKI67')

umap_df = data.frame (srt[[reductionName]]@cell.embeddings)
markers_found = rownames(srt)[rownames(srt) %in% markers]
umap_df = cbind (umap_df, t(srt@assays$RNA@data[markers_found,]))
  umap_p1 = lapply (markers_found, function(x) ggplot(data = umap_df) + 
  geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
  #scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,x])), max (abs(umap_df[,x])))) +
  ggtitle (x) + 
  #paletteer::scale_color_paletteer_c("ggthemes::Green", limits=c(0, max (abs(umap_df[,x])))) +
  scale_colour_gradientn (colours = c('lightgrey','darkred'),limits=c(0, max (abs(umap_df[,x])))) +
  #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
  theme_void() + NoLegend())

png ('Plots/FIGURE_4B_celltypes_markers_fplots2.png',width = 2000,1600,res=300)
print (wrap_plots (umap_p1))
dev.off()


# FIGURE S4B - dotplot of Myeloid subtypes markers ####
top_markers = c('CLEC9A','IDO1','DNASE1L3','CLNK','TACSTD2','HLA-DPB1','HLA-DPA1','HLA-DQA1','HLA-DRA','FTL','TPSB2','TPSAB1','CPA3','MS4A2','GATA2','FCN1','LYZ','VCAN','S100A6','C1QB','LST1','APOE','APOC1','C1QC','C1QA','LAMP3','CCL19','LAD1','NCCRP1','CCR7','GZMB','JCHAIN','CLIC3','LILRA4','ITM2C')
srt$celltype = factor (srt$celltype, levels = c('cDC1','cDC2','Mast','Mono_CD14','Mono_CD16','TAMs','mregDC','pDC'))
dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = factor (top_markers, levels = top_markers),
  x = 'celltype', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = NULL,
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  swap_axes = T,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression2) + gtheme_italic
    

pdf ('Plots/FIGURE S4B_top_markers_celltype_expression.pdf', height=2.8, width=8.5)
print (dp)
dev.off()  


# Show barplots of NMF cell abundance
p = cellComp (
  seurat_obj = srt, 
  metaGroups = c('celltype','sampleID'),
  plot_as = 'bar',
  ptable_factor = c(1),
  prop = TRUE,
  pal = palette_sample
  #facet_ncol = 15
  ) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
png ('Plots/FIGURE_S4B_composition_celltype_barplot.png',1200,900, res=300)
print (p)
dev.off()

# Run SCENIC plots #### 
motif_window = 'tss500bp'

# Generate heatmap of TFs ####
motifs_table = read.csv ('../../PM_scRNA_atlas/data/SCENIC_myeloid_motifs.csv', skip=2) 
auc_mtx = read.csv ('../../PM_scRNA_atlas/data/SCENIC_myeloid_auc_mtx.csv')
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = paste0('SCENIC_',colnames(auc_mtx))
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

pValThreshold = 1e-15
auc_mtx_avg = aggregate (auc_mtx, by=as.list(srt@meta.data[,'celltype',drop=F]), mean)
rownames (auc_mtx_avg) = auc_mtx_avg[,1]
auc_mtx_avg = auc_mtx_avg[,-1]
colnames (auc_mtx_avg) = sub ('SCENIC_','', colnames(auc_mtx_avg))
#palette_scenic = rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(50))
filter_TF = unique (motifs_table$TF[motifs_table$X.2 < pValThreshold])
hm = Heatmap (scale (auc_mtx_avg[,filter_TF]),
 col = palette_scenic, column_names_rot = 45,
 column_names_gp = gpar(fontsize = 5),
 row_names_gp = gpar(fontsize = 8))
pdf (paste0('Plots/FIGURE_S4F_auc_scores_celltypes_heatmap.pdf'), width=7, height=2.2)
print (hm)
dev.off()
  

### FIGURE S4C - TF expression across celltypes ####
TF_G = read.table('../../PM_scRNA_atlas/data/HumanTFs.txt', sep="\t", header = TRUE, row.names=1)

celltype_order = c('Mono_CD14','Mono_CD16','TAMs','cDC1','cDC2','mregDC','pDC','Mast')
Idents (srt) = 'celltype'
levels (srt) = celltype_order
degClusters = FindAllMarkers (srt[rownames(srt)%in% TF_G$HGNC.symbol,], max.cells.per.ident = 1000, group.by = 'celltype',min.pct = .1, verbose = T)
degClusters = degClusters[order (-degClusters$avg_log2FC),]
degClusters_flt = degClusters %>% group_by (cluster) %>% slice_head(n=2)
cluster.averages = AverageExpression(object = srt[degClusters_flt$gene,], group.by = 'celltype', return.seurat = F)
cluster.averages = as.data.frame (cluster.averages[[1]])
cluster.averages = cluster.averages[, celltype_order]

hm = Heatmap (
  t(scale (t(cluster.averages))), 
  cluster_rows=F,
  row_names_side = 'left',
  cluster_columns=F, 
  col = palette_gene_expression_fun(cluster.averages),
  column_names_rot = 45,
  row_names_gp = gpar (fontsize = 6),
  column_names_gp = gpar (fontsize = 6))
pdf(paste0("Plots/FIGURE_S4C_TF_celltype_heatmap.pdf"),width =2.2,height=2, useDingbats=T)
hm
dev.off()



### subset tams ####
srt_tam = srt[, srt$celltype == 'TAMs']
cnmf_spectra_unique_comb = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Mms_20"))

srt_tam = ModScoreCor (
        seurat_obj = srt_tam, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Mms', outdir = paste0(projdir,'Plots/'))

### FIGURE S4C - pairwise spearman correlation across cells ####
# Generate heatmap of nmf modules correlation across cells
# Add track showing corelation of each module to sarcomatoid module

# Load scS-score
scs_sample_avg = read.csv ('../../PM_scRNA_atlas/data/scs_score_per_sample.csv', row.names=1)

ccomp_df = srt_tam@meta.data[,names (cnmf_spectra_unique_comb)]
cor_mat = cor (ccomp_df, method = 'spearman')

color_heatmap_palette <- colorRamp2(c(min(cor_mat), 0, max(cor_mat)), c(palette_module_correlation[1], palette_module_correlation[length(palette_module_correlation)/2], palette_module_correlation[length(palette_module_correlation)]))
# triangle heatmap
hm = draw (Heatmap (cor_mat, 
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=color_heatmap_palette, 
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_S4C_cnmf_cor_heatmap_triangle_cells.pdf'),5.5,4.5)
print (hm)
dev.off()


### FIGURE S4D - Make dotplot of top markers for each nmf ####
cnmf_spectra_unique_comb_ordered = cnmf_spectra_unique_comb[row_order(hm)]
marker_genes = unlist (lapply (cnmf_spectra_unique_comb_ordered, function(x) head (x, 5)))
srt_tam$Mms_r_max = factor (srt_tam$Mms_r_max, levels = rev(names(cnmf_spectra_unique_comb_ordered)))
dotp = geneDot (
srt_tam,
#gene = top_tfs2, 
gene = factor (marker_genes, levels = marker_genes),
x = 'Mms_r_max', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
swap_axes = T,
x_name ='samples',
y_name = 'celltype',
plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_S4D_cNMF_markers_expression.pdf', width=11, height=4)
print (dotp)
dev.off()


### FIGURE 4C - pairwise spearman correlation across samples ####
# Generate heatmap of nmf modules correlation across samples
# Add track showing corelation of each module to sarcomatoid module
ccomp_df = srt_tam@meta.data[,names(cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_tam@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
ccomp_df = cbind (ccomp_df, scScore = scs_sample_avg[rownames(ccomp_df),])
ccomp_df = na.omit (ccomp_df)
cor_mat = cor (ccomp_df, method = 'spearman')

scScore_cor = cor_mat[,'scScore']
cor_mat = cor_mat[-ncol(cor_mat), -ncol(cor_mat)]
scScore_cor = scScore_cor[-length(scScore_cor)]

# triangle heatmap
ha = HeatmapAnnotation(scScore = scScore_cor, col = list (scScore = palette_sample2_fun), which='row')
hm = draw (Heatmap (cor_mat, 
  left_annotation = ha,
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation, 
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_4C_modules_cnmf_cor_heatmap_triangle.pdf'),5.5,4.5)
print (hm)
dev.off()


# FIGURE 4E - Test CXCLs expression in TAMs ####
ext_markers = c(
  'CXCL9','CXCL10','CXCL11')

ext_avg = AverageExpression (srt_tam, features = ext_markers, group.by = c('sampleID','SE_group'))
ext_avg = log2 (as.data.frame (t(ext_avg[[1]]))+1)
ext_avg$SE_group = sapply (rownames(ext_avg), function(x) unlist(strsplit (x,'_'))[2])
ext_avg$SE_group = factor (ext_avg$SE_group, levels = c('S-High','E-High'))
ext_avg = gather (ext_avg, gene, expression, 1:(ncol (ext_avg) - 1))
ext_avg$gene = factor (ext_avg$gene, levels = ext_markers)  
box = ggplot (ext_avg, aes_string (x= 'SE_group', y= 'expression')) +
geom_point(position=position_jitter(width=0.1), alpha=1, color="grey", size=0.7) +
geom_boxplot (aes_string(fill='SE_group'),alpha = 0.7, lwd=.2, outlier.shape = NA) +
scale_fill_manual (values = palette_SE_group) +
facet_wrap (~gene, drop=TRUE, scales = 'free', ncol=length (ext_markers)) + gtheme_no_text

stat.test = box$data %>%
group_by (gene) %>%
t_test(reformulate ('SE_group', 'expression')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'SE_group', step.increase=.4)
box = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = 0, hide.ns = T,
label = "p.adj.signif") + NoLegend()

pdf ('Plots/FIGURE_4E_CXCLs_markers.pdf',2.5,width = 3)
print (box)
dev.off()


# FIGURE 4H - SCENIC regulons heatmap ####
# Run SCENIC plots #### 
motif_window = 'tss500bp'

# Generate heatmap of TFs ####
motifs_table = read.csv ('../../PM_scRNA_atlas/data/SCENIC_myeloid_motifs.csv', skip=2) 
auc_mtx = read.csv ('../../PM_scRNA_atlas/data/SCENIC_myeloid_auc_mtx.csv')
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = paste0('SCENIC_',colnames(auc_mtx))
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

pValThreshold = 1e-15
auc_mtx_avg = aggregate (auc_mtx[match (colnames(srt_tam), rownames(auc_mtx)),], by=as.list(srt_tam@meta.data[,'Mms_r_max',drop=F]), mean)
rownames (auc_mtx_avg) = auc_mtx_avg[,1]
auc_mtx_avg = auc_mtx_avg[,-1]
colnames (auc_mtx_avg) = sub ('SCENIC_','', colnames(auc_mtx_avg))
#palette_scenic = rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(50))
filter_TF = unique (motifs_table$TF[motifs_table$X.2 < pValThreshold])
hm = Heatmap (scale (auc_mtx_avg[,filter_TF]),
 col = palette_scenic, column_names_rot = 45,
 column_names_gp = gpar(fontsize = 5),
 row_names_gp = gpar(fontsize = 8))
pdf (paste0('Plots/FIGURE_4H_auc_scores_Mms_heatmap.pdf'), width=7, height=2.2)
print (hm)
dev.off()



### FIGURE 4F - dotplot CXCLs and receptor across celltypes ####
srt_tumor$celltype_simplified3 = as.character(srt_tumor$celltype_simplified)
srt_tumor$celltype_simplified3[srt_tumor$celltype %in% c('Mono_CD16','Mono_CD14','TAMs')] = srt_tumor$celltype[srt_tumor$celltype %in% c('Mono_CD16','Mono_CD14','TAMs')]
#srt$celltype_simplified3[srt$celltype %in% c('pDC')] = 'DCs'
srt_tumor$celltype_simplified3[srt_tumor$celltype_simplified3 %in% c('T_cells')] = srt_tumor$celltype[srt_tumor$celltype_simplified3 %in% c('T_cells')]
srt_tumor$celltype_simplified3 = factor (srt_tumor$celltype_simplified3, levels = rev(c('Malignant','Mesothelium','Glia','Alveolar','Fibroblasts','Endothelial','SmoothMuscle','Mono_CD16','Mono_CD14','TAMs','Mast','DCs','CD8','CD4','TFH','Tregs','NK','B_cells','Plasma','pDCs')))
cxcl_markers = c('CXCL9','CXCL10','CXCL11','CXCR3')

dp = geneDot (
  seurat_obj = srt_tumor,
  gene = cxcl_markers,
  x = 'celltype_simplified3',
  y = 'SE_group',
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  swap_axes = T,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_4F_cxcls_and_receptor_expression2.pdf', height=4.5, width=4)
print (dp)
dev.off()  

### FIGURE 4C - VISTA expression ####
ext_markers='VSIR'
low_celltypes = c('Alveolar','LEC','Mesothelium','Glia')
srt_tumor2 = srt_tumor[,!srt_tumor$celltype %in% low_celltypes]
srt_tumor2 = srt_tumor2[,srt_tumor2$sampleID != 'P10']
srt_tumor2$celltype = gsub ('_','',srt_tumor2$celltype)
ext_avg = AverageExpression (srt_tumor2, features = ext_markers, group.by = c('sampleID','celltype'))[[1]]
rownames (ext_avg) = ext_markers
ext_avg = log2(as.data.frame (t(ext_avg))+1)
ext_avg$group = rownames (ext_avg)

ext_avg = gather (ext_avg, gene, avg_expression, 1:(ncol(ext_avg)-1))
ext_avg$sample = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[1])
ext_avg$celltype = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[2])
ext_avg$SE_group = setNames(srt_tumor2$SE_group, srt_tumor2$sampleID)[ext_avg$sample]
ext_avg$SE_group = factor (ext_avg$SE_group, levels = c('S-High','E-High'))
celltype_order = ext_avg
celltype_order =  do.call (rbind,lapply (split (celltype_order, celltype_order$celltype), function(x) data.frame(mean(x$avg_expression), x$celltype[1])))
celltype_order = celltype_order[,2][order (-celltype_order[,1])]
ext_avg$celltype = factor (ext_avg$celltype, levels = celltype_order)
boxl = list()
for (x in ext_markers)
  {
  ext_avg_sub = ext_avg[ext_avg$gene == x,]
  box = ggplot (ext_avg_sub, aes_string (x= 'celltype', y= 'avg_expression')) +
  geom_point(position=position_jitter(width=0.1), alpha=1, color="grey", size=0.7) +
  geom_boxplot (aes_string(fill='SE_group'),alpha = 0.7, lwd=.2, outlier.shape = NA) +
  ggtitle (x) +
  scale_fill_manual (values = palette_SE_group) +
  gtheme
  
  stat.test = box$data %>% 
  group_by (celltype) %>%
  t_test (reformulate ('SE_group', 'avg_expression')) %>%
  adjust_pvalue (method = "fdr") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = 'celltype', step.increase=.01)
  boxl[[x]] = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = T,
  label = "p.adj.signif") + NoLegend()
  }

pdf ('Plots/FIGURE_4C_VSIR_per_celltype.pdf',width = 5, height=4)
print (boxl)
dev.off()

