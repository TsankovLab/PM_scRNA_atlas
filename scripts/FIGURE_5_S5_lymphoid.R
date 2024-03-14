# TNK compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/tnk/'
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
srt = srt_tumor[, srt_tumor$celltype_simplified %in% c('T_cells','NK')]

# load palettes
# palettes = readRDS (paste0('../../palettes.rds'))
# for (i in seq_along (palettes)) assign (names(palettes)[i], palettes[[i]])

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

# FIGURE 5A / S5A - celltype UMAP ####
dp = DimPlot (srt, group.by = 'celltype', 
  reduction=reductionName) + theme_void() + scale_color_manual (values = palette_tnk_cells)
dp1 = DimPlot (srt, group.by = 'sampleID', 
  reduction=reductionName) + theme_void() + scale_color_manual (values = palette_sample)

pdf ('Plots/FIGURE_5A_S5A_celltypes_umap.pdf',5,3)
print (dp)
print (dp1)
dev.off()


# FIGURE 5B - dotplot of T cell subtype markers ####
top_markers = c('IL7R','SELL','CCR7','FOXP3','IL2RA','TNFRSF18','CXCL13','IL21','TOX2','CD8A','CD8B','CCL5','CD3D','KLRC1','XCL1','FCER1G','GNLY','FGFBP2','SPON2')
srt$celltype = factor (srt$celltype, levels = c('CD4', 'Tregs','TFH','CD8','NKlike_Tcells','KLRC1_NK','FGFBP2_NK'))

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
  x_name ='samples',
  swap_axes = T,
  y_name = 'celltype',
  plotcol = palette_gene_expression2) +
    gtheme_italic


pdf ('Plots/FIGURE_5B_top_markers_celltype_expression.pdf', height=2.8, width=5.5)
print (dp)
dev.off()  

### FIGURE 5B - Make stacked barplot of sample abundance per celltype ####
cp = cellComp (
srt, 
metaGroups = c('celltype','sampleID'), 
plot_as = 'bar',
pal = palette_sample,
ptable_factor = 1,
prop=T) +
theme_minimal()
#cp$data$cNMF__r_max = factor (cp$data$cNMF__r_max, levels = names (cnmf_spectra_filtered[row_order(hm)])) 
pdf (paste0('Plots/FIGURE_5B_sample_abundance_celltype_stacked_barplot.pdf'))
print (cp)
dev.off()

# Import cNMF modules from T cells and compute scores across TNK cells for each 
cnmf_spectra_unique_comb = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Tms_20"))

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Tms', outdir = paste0('Plots/'))

# Subset seurat object only for T cells populations
srt_t = srt[, srt$celltype %in% c('CD4','Tregs','TFH','CD8')]

### FIGURE S5C - pairwise spearman correlation across cells ####
ccomp_df = srt_t@meta.data[,names (cnmf_spectra_unique_comb)]
cor_mat = cor (ccomp_df, method = 'spearman')

# triangle heatmap
hm = draw (Heatmap (cor_mat, 
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_S5C_cnmf_cor_heatmap_triangle_cells.pdf'),5.5,4.5)
print (hm)
dev.off()

### FIGURE S5B - Make dotplot of top markers for each nmf ####
cnmf_spectra_unique_comb_ordered = cnmf_spectra_unique_comb[row_order(hm)]
marker_genes = unlist (lapply (cnmf_spectra_unique_comb_ordered, function(x) head (x, 5)))
srt_t$Tms_r_max = factor (srt_t$Tms_r_max, levels = rev(names(cnmf_spectra_unique_comb_ordered)))
dotp = geneDot (
srt_t,
#gene = top_tfs2, 
gene = factor (marker_genes, levels = marker_genes),
x = 'Tms_r_max', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
y = NULL,
swap_axes = T,
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_S5B_cNMF_markers_expression.pdf', width=11, height=4)
print (dotp)
dev.off()


### FIGURE 5C - pairwise spearman correlation across samples ####
ccomp_df = srt_t@meta.data[,names(cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt_t@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
ccomp_df = cbind (ccomp_df, scScore = scs_sample_avg[rownames(ccomp_df),])
ccomp_df = na.omit (ccomp_df)
cor_mat = cor (ccomp_df, method = 'spearman')
#col_fun = colorRamp2(c(1, 0, -1), c(palette_sample[3], palette_sample[length(palette_sample)/2], palette_sample[length(palette_sample)-3]))
scScore_cor = cor_mat[,'scScore']
cor_mat = cor_mat[-ncol(cor_mat), -ncol(cor_mat)]
#palette_module_correlation = paletteer::paletteer_c("pals::kovesi.diverging_bwr_40_95_c42",100)
scScore_cor = scScore_cor[-length(scScore_cor)]

# triangle heatmap
ha = HeatmapAnnotation(scScore = scScore_cor, col = list (scScore = palette_sample2_fun), which='row')
hm = draw (Heatmap (cor_mat, 
  left_annotation = ha,
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation, 
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_5C_modules_cnmf_cor_heatmap_triangle.pdf'),5.5,4.5)
print (hm)
dev.off()


# FIGURE 5E - Test individually each exhaustion marker ####
ext_markers = c(
  'TIGIT','PDCD1','CTLA4','LAG3','HAVCR2')

ext_avg = AverageExpression (srt_t, features = ext_markers, group.by = c('sampleID','SE_group'))
ext_avg = log2(as.data.frame (t(ext_avg[[1]]))+1)
ext_avg$SE_group = sapply (rownames(ext_avg), function(x) unlist(strsplit (x,'_'))[2])
ext_avg$SE_group = factor (ext_avg$SE_group, levels = c('S-High','E-High'))
ext_avg = gather (ext_avg, gene, expression, 1:(ncol (ext_avg) - 1))
ext_avg$gene = factor (ext_avg$gene, levels = ext_markers)  
box = ggplot (ext_avg, aes_string (x= 'SE_group', y= 'expression')) +
#geom_violin (trim=TRUE, aes_string (fill = 'SE_group'),alpha = 0.7) +
geom_point(position=position_jitter(width=0.1), alpha=1, color="grey", size=0.7) +
geom_boxplot (aes_string(fill='SE_group'),alpha = 0.7, lwd=.2, outlier.shape = NA) +
gtheme_no_text +
scale_fill_manual (values = palette_SE_group) +
facet_wrap (~gene, drop=TRUE, scales = 'free_x', ncol=length (ext_markers))

stat.test = box$data %>%
group_by (gene) %>%
t_test(reformulate ('SE_group', 'expression')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'SE_group', step.increase=.4)
box = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = -1.2, hide.ns = T,
label = "p.adj.signif") + NoLegend()

pdf ('Plots/FIGURE_5E_exhaustion_markers.pdf',2.5,width = 2.8)
print (box)
dev.off()

# FIGURE S5E - Repeat per subtype ####
ext_avg = AverageExpression (srt_t, features = ext_markers, group.by = c('sampleID','celltype'))
ext_avg = log2(as.data.frame (t(ext_avg[[1]]))+1)
ext_avg$group = rownames (ext_avg)
ext_avg = gather (ext_avg, gene, avg_expression, 1:(ncol(ext_avg)-1))
ext_avg$sample = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[1])
ext_avg$celltype = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[2])
ext_avg$SE_group = setNames(srt_t$SE_group, srt_t$sampleID)[ext_avg$sample]
ext_avg$SE_group = factor (ext_avg$SE_group, levels = c('S-High','E-High'))
ext_avg = ext_avg[!is.na(ext_avg$SE_group),]

boxl = list()
for (x in ext_markers)
  {
  ext_avg_sub = ext_avg[ext_avg$gene == x,]
  box = ggplot (ext_avg_sub, aes_string (x= 'SE_group', y= 'avg_expression')) +
  geom_point(position=position_jitter(width=0.1), alpha=1, color="grey", size=0.7) +
  geom_boxplot (aes_string(fill = 'SE_group'), alpha = 0.7, lwd=.2, outlier.shape = NA) +
  gtheme_no_text +
  ggtitle (x) +
  scale_fill_manual (values = palette_SE_group) +
  facet_wrap (~celltype, drop=TRUE, scales = 'free_x', ncol=4)
  
  stat.test = box$data %>% 
  group_by (celltype) %>%
  t_test (reformulate ('SE_group', 'avg_expression')) %>%
  adjust_pvalue (method = "fdr") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = 'SE_group', step.increase=.01)
  boxl[[x]] = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = T,
  label = "p.adj.signif") + NoLegend()
  }

pdf ('Plots/FIGURE_S5E_exhaustion_markers_per_celltype.pdf',width = 12, height=4)
print (wrap_plots (boxl, ncol=5))
dev.off()

# FIGURE 5F - Repeat per subtype across samples to generate heatmap ####
ext_avg = AverageExpression (srt_t, features = ext_markers, group.by = c('SE_group', 'celltype'))
ext_avg = log2(as.data.frame (t(ext_avg[[1]]))+1)
#ext_avg$group = rownames (ext_avg)
# ext_avg = gather (ext_avg, gene, avg_expression, 1:(ncol(ext_avg)-1))
# ext_avg$sample = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[1])
ext_avg = t (ext_avg)
col_split2 = sapply (colnames(ext_avg), function(x) unlist(strsplit (x, '_'))[2])
col_split2 = col_split2[order(col_split2)]
ext_avg = ext_avg[,names(col_split2)]
ext_avgs = scale (t(ext_avg))
hm = draw (Heatmap (t(ext_avgs), 
  cluster_columns=F,
  column_split = col_split2,
  clustering_distance_rows='pearson',
  clustering_distance_columns = 'pearson', 
  col=palette_gene_expression_fun(ext_avgs), 
  heatmap_legend_param= list (direction = 'horizontal'),
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=T))
pdf ('Plots/FIGURE_5F_exhaustion_markers_per_celltype_heatmap.pdf',width = 4, height=2.8)
print (hm)
dev.off()
  
# FIGURE 5E - Show distribution relevant Tms ####
ccomp_df = srt@meta.data
dp = lapply (c('Tm2','Tm5','Tm7'), function(x) ggplot(ccomp_df, aes_string(x = x, group = 'SE_group')) +
  geom_density(aes(fill = SE_group),linewidth=.1, alpha = 0.6) + scale_fill_manual (values = palette_SE_group) + gtheme + NoLegend())
pdf ('Plots/FIGURE_5E_cell_distributions_Tm2_Tm5_Tm7_density.pdf',width = 4, height=1.6)
print (wrap_plots (dp))
dev.off()


# FIGURE 5G - MHCII per subtype ####
ext_markers = c('HLA-DPB1','HLA-DPA1', 'HLA-DQA1','HLA-DRB5','GZMK')
ext_avg = AverageExpression (srt_t, features = ext_markers, group.by = c('sampleID','SE_group'))
ext_avg = log2(as.data.frame (t(ext_avg[[1]]))+1)
ext_avg$SE_group = sapply (rownames(ext_avg), function(x) unlist(strsplit (x,'_'))[2])
ext_avg$SE_group = factor (ext_avg$SE_group, levels = c('S-High','E-High'))
ext_avg = gather (ext_avg, gene, expression, 1:(ncol (ext_avg) - 1))
ext_avg$gene = factor (ext_avg$gene, levels = ext_markers)  
box = ggplot (ext_avg, aes_string (x= 'SE_group', y= 'expression')) +
#geom_violin (trim=TRUE, aes_string (fill = 'SE_group'),alpha = 0.7) +
geom_point(position=position_jitter(width=0.1), alpha=1, color="grey", size=0.7) +
geom_boxplot (aes_string(fill='SE_group'),alpha = 0.7, lwd=.2, outlier.shape = NA) +
gtheme_no_text +
scale_fill_manual (values = palette_SE_group) +
facet_wrap (~gene, drop=TRUE, scales = 'free_x', ncol=length (ext_markers))

stat.test = box$data %>%
group_by (gene) %>%
t_test(reformulate ('SE_group', 'expression')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'SE_group', step.increase=.4)
box = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = -1, hide.ns = T,
label = "p.adj.signif") + NoLegend()

pdf ('Plots/FIGURE_5G_MHCII_markers.pdf',2.5,width = 2.8)
print (box)
dev.off()

# Subset to B-cells and recompute UMAP
srt_b = srt_tumor[, srt_tumor$celltype_simplified %in% c('B_cells','Plasma')]

batch = 'sampleID'
reductionSave = paste0(paste(batch,collapse='_'),'_harmony')
reductionKey = paste0(paste(batch,collapse='_'),'harmonyUMAP_')
reductionName = paste0 (paste(batch,collapse='_'),'_harmony_umap')
reductionGraphKnn = paste0 (paste(batch,collapse='_'),'_harmony_knn')
reductionGraphSnn = paste0 (paste(batch,collapse='_'),'_harmony_snn')

# Compute variable features using scran pkg
nfeat=3000
sce = SingleCellExperiment (list(counts=srt_b@assays$RNA@counts, logcounts = srt_b@assays$RNA@data),
rowData=rownames(srt_b)) 
sce = modelGeneVar(sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
vf = getTopHVGs(sce, n=nfeat)
VariableFeatures (srt_b) = vf

# Process merged data
srt_b = NormalizeData (object = srt_b, normalization.method = "LogNormalize", scale.factor = 10000)
srt_b = ScaleData (srt_b, features = VariableFeatures (object=srt_b))
srt_b = RunPCA (srt_b, features = VariableFeatures (object = srt_b), npcs = ifelse(ncol(srt_b) <= 30,ncol(srt_b)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
  
# Run Harmony
srt_b = srt_b %>% 
RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
RunUMAP (reduction = reductionSave, dims = 1:15, reduction.name = reductionName, reduction.key=reductionKey)

# Run denovo clustering on non-adjusted reductions
srt_b = FindNeighbors (object = srt_b, reduction = reductionSave, dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T, graph.name=c(reductionGraphKnn,reductionGraphSnn))

# FIGURE S5B -- markers dotplot ####
top_markers = c('CD79B','CD79A','CD19','CD37','RGS13','LRMP','MEF2B','LMO2','HMCES','IGHG3','IGHGP','IGHG4')

dp = geneDot (
  seurat_obj = srt_b,
  gene = factor (top_markers, levels = top_markers),
  x = 'celltype', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  swap_axes=T,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression2) +
  gtheme_italic +
    theme(legend.position="bottom")


metaGroupNames = c('celltype','sampleID')
ccc_bar1 = cellComp (
  seurat_obj = srt_b, 
  metaGroups = metaGroupNames,
  plot_as = 'bar',
  pal = palette_sample,
  prop = TRUE,
  ptable_factor = c(1),
  #subset_prop = 'cycling',
  facet_ncol = 6,
  facet_scales = 'free'
  ) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ccc_bar1$layers[[1]]$aes_params$alpha =   1


pdf ('Plots/FIGURE_5H_top_markers_celltype_expression_legend.pdf', height=4, width=7)
print (dp)
dev.off()
 
pdf ('Plots/FIGURE_5H_top_markers_celltype_expression.pdf', height=2.7, width=4)
print (dp)
dev.off()

pdf(paste0 ('Plots/FIGURE_5H_cell_composition_barplot.pdf'), width=2.5, height=4)
print (ccc_bar1)
dev.off()

### FIGURE 5I - DimPlot ####
dp = DimPlot (srt_b, group.by = 'celltype', cols = palette_b_cells, 
  reduction=reductionName) + theme_void()
pdf ('Plots/FIGURE_5I_celltypes_umap.pdf',width= 4,3)
print (dp)
dev.off()

### FIGURE 5I - cellcycle distribution ####
#theme_set(theme_classic())
cc.genes <- readLines('../../PM_scRNA_atlas/data/regev_lab_cell_cycle_genes.txt')
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
srt_b = CellCycleScoring (object = srt_b, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
srt_b$cc = srt_b$S.Score + srt_b$G2M.Score

ccomp_df = srt_b@meta.data[, c('cc','celltype')]
rp = ggplot(ccomp_df, aes_string(x = 'cc', y = 'celltype')) +
  geom_density_ridges(aes(fill = celltype), alpha = 0.7, color = NA) +
  scale_fill_manual(values = palette_b_cells) +
  gtheme

pdf ('Plots/FIGURE_5I_cellcycle_GCBcells.pdf', height=2.5, width=3)
print (wrap_plots (rp))
dev.off()


### FIGURE 5H - cellcycle distribution ####
gctfh_markers = c('TOX2','CXCR5','PDCD1','BCL6','TCF7')
srt_tfh = srt[,srt$celltype == 'TFH']
dotp3 = geneDot (
seurat_obj = srt_tfh,
gene = gctfh_markers,
x = 'sampleID', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
swap_axes=T,
#scale.data=T,
x_name ='samples',
y_name = 'celltype',
plotcol = palette_gene_expression2) + gtheme_italic

pdf (paste0 ('Plots/FIGURE_5H_gctfh_markers_TFH_dotplot.pdf'), height=3.5, width=2.7)
print (dotp3)
dev.off()




#####################
### TCR analysis ####
#####################
# Import TCR of tumors
data.path_P7_tumor = '../../PM_scRNA_atlas/data/P7_filtered_contig_annotations.csv'
data.path_P6_tumor = '../../PM_scRNA_atlas/data/P6_filtered_contig_annotations.csv'
data.path_P2_tumor = '../../PM_scRNA_atlas/data/P2_filtered_contig_annotations.csv'
data.path_P9_tumor = '../../PM_scRNA_atlas/data/P9_filtered_contig_annotations.csv'
data.path_P11_tumor = '../../PM_scRNA_atlas/data/P11_filtered_contig_annotations.csv'
data.path_P12_tumor = '../../PM_scRNA_atlas/data/P12_filtered_contig_annotations.csv'
data.path_P13_tumor = '../../PM_scRNA_atlas/data/P13_filtered_contig_annotations.csv'


vdj.dirs = c(
  data.path_P7_tumor, 
  data.path_P6_tumor, 
  data.path_P2_tumor, 
  data.path_P9_tumor,
  data.path_P11_tumor,
  data.path_P12_tumor,
  data.path_P13_tumor)
tcrL_tumor = lapply (vdj.dirs, function(x) read.csv(x))
names (tcrL_tumor) = c('P7','P6','P2','P9','P11','P12','P13')

contig_list = c(tcrL_tumor)
combinedTCR <- combineTCR (contig_list, 
                samples = names (contig_list), 
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)
pdf (paste0('Plots/clonotype_overlap.pdf'),width=5,3.5)
print (clonalOverlap(combinedTCR, 
              cloneCall = "strict", 
              method = "morisita"))
dev.off()

pdf (paste0('Plots/occupied_repertoire.pdf'),width=4,height=4)
print (clonalQuant(combinedTCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE))
dev.off()


srt$barcode = sapply (colnames(srt), function(x) unlist(strsplit (x, '_'))[2])
srt = RenameCells(
   srt,
   new.names = paste0(srt$sampleID,'_', srt$barcode))

srt = combineExpression (combinedTCR, srt,
  proportion = F,
#  cloneSizes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = Inf))
  cloneSize = c(NonExpanded = 1, Small = 5, Large = Inf))
  #)
srt_tcr = srt[,srt$sampleID %in% c('P7','P6','P9','P11','P12','P13')]
srt_tcr = srt_tcr[, !is.na (srt_tcr$cloneSize)]
clonesize_name = setNames (c('NonExpanded','Small','Large','None'), c('NonExpanded (0 < X <= 1)','Small (1 < X <= 5)','Large (5 < X <= Inf)', 'None ( < X <= 0)'))
srt_tcr$cloneSize = unname(clonesize_name[as.character(srt_tcr$cloneSize)])

# FIGURE 5I - Show expanded clonotypes are found in CD8 cells ####  
dp = DimPlot (srt_tcr, reduction = reductionName, group.by='cloneSize') + 
scale_color_manual (values = palette_clonotype) + gtheme

pdf(paste0('Plots/Expanded_vs_nonExpanded_umap.pdf'), height=3,width=5)
print (dp)
dev.off()

metaGroupName = 'cloneSize'
cc_box = cellComp (
  seurat_obj = srt_tcr, 
  metaGroups = c('celltype',metaGroupName),
  plot_as = 'bar',
  prop = FALSE,
  pal = palette_clonotype,
  subset_prop = 'Large',  
  facet_ncol = 15
  ) + gtheme
cc_box$data$celltype = factor (cc_box$data$celltype, levels = as.character(cc_box$data$celltype[order(-cc_box$data$Freq)]))

pdf(paste0('Plots/FIGURE_5A_Expanded_vs_nonExpanded_barplot.pdf'),4,3)
print (cc_box)
dev.off()


# FIGURE 6B - Show expanded clonotypes have higher exhaustion ####  
# Add modulesscore of exhaustion markers
# Add exhaustion modules
srt_tcr = srt_tcr[,srt_tcr$celltype == 'CD8']
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm5, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA) +
gtheme_no_text

pdf ('Plots/FIGURE_S5J_expanded_exhaustion.pdf', height = 3.3,width = 2)
print (box2)
dev.off()


### Check with MHC II module ####
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm4, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point (position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA) +
gtheme_no_text

pdf ('Plots/FIGURE_S5J_expanded_MHCII.pdf',height = 3.3,width = 2)
print (box2)
dev.off()


### Check with Cytotoxic module ####
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm2, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA) +
gtheme_no_text


pdf ('Plots/FIGURE_S5J_expanded_cytox.pdf',height = 3.3,width = 2)
print (box2)
dev.off()




# FIGURE S5 - proportion of CD8 expanded and exhausted in SE groups ####
metaGroupName = 'cloneSize'
cc_box = cellComp (
  seurat_obj = srt_tcr, 
  metaGroups = c('sampleID',metaGroupName, 'SE_group'),
  plot_as = 'box',
  ptable_factor = c(1),
  prop = T,
  pal = palette_SE_group,
  subset_prop = 'Large',
  facet_ncol = 15
  ) + gtheme_no_text

cc_box$data$SE_group = factor (cc_box$data$SE_group, levels = unique(cc_box$data$SE_group))
cc_box$data$cloneSize = factor (cc_box$data$cloneSize, levels = rev (c('Large','Small','NonExpanded')))
stat.test = cc_box$data %>%
t_test(reformulate ('SE_group', 'Freq')) %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'SE_group', step.increase=.4)
cc_box = cc_box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = 0, hide.ns = T,
label = "p")

pdf(paste0('Plots/FIGURE_S5K_Exhausted_Expanded_vs_nonExpanded_SE_group_boxplot2.pdf'),width=2.2,2.5)
print (cc_box)
dev.off()



##########################################
### Compare tumor and pbmc clonotypes ####
##########################################
# Load seurat objects from tumor and PBMC ####
# PBMC (select only PRE-polyICLC)
srt_pbmc = readRDS ('../srt_pbmc.rds')
srt_pbmc = srt_pbmc[,srt_pbmc$batch == 'batch2']

# Tumor
srt_tumor = srt[,srt$sampleID %in% c('P7','P2','P9')]

# Harmonize barcode suffixes and merge seurat object ####
# PBMC 
#sample_match = setNames (c('p8','p4','p9'), c('P2','P7','P9'))
#srt_pbmc$sampleID4 = sample_match[srt_pbmc$sampleID3]
meta_pbmc = srt_pbmc@meta.data
srt_pbmc = srt_pbmc[,!duplicated (paste0(srt_pbmc$sampleID, '_',srt_pbmc$barcode))]
srt_pbmc = RenameCells(
    srt_pbmc,
    new.names = paste0(srt_pbmc$sampleID, '_pbmc_',srt_pbmc$barcode))

# Tumor 
srt_tumor$barcode = sapply (colnames (srt_tumor), function(x) unlist (strsplit (x, '_'))[2])
srt_tumor$sampleID = as.character(srt_tumor$sampleID)
srt_tumor = srt_tumor[,!duplicated (paste0(srt_tumor$sampleID, '_tumor_',srt_tumor$barcode))]
srt_tumor = RenameCells(
    srt_tumor,
    new.names = paste0(srt_tumor$sampleID, '_tumor_',srt_tumor$barcode))
# Merge
srt_merged = merge (srt_tumor, srt_pbmc)
srt_merged$site = ifelse (grepl ('pbmc',colnames (srt_merged)), 'pbmc','tumor')

# Import TCR_contigs of PBMCs ####
data.path_pbmc1 = '../../PM_scRNA_atlas/data/PBMC_P2_P7_P9_Processed_TCR_1-1_filtered_contig_annotations.csv'
data.path_pbmc2 = '../../PM_scRNA_atlas/data/PBMC_P2_P7_P9_Processed_TCR_1-2_filtered_contig_annotations.csv'
vdj.dirs = c(data.path_pbmc1, data.path_pbmc2)

# Add hashing pools to seurat metadata ####
#meta_pbmc$hash_pool = sapply (rownames(meta_pbmc), function(x) unlist(strsplit (x, '\\-'))[2])
#meta_pbmc$hash_pool = sapply (meta_pbmc$pool_batch, function(x) unlist(strsplit (x, '\\-'))[2])

# Strip barcodes to match TCR data ####
#meta_pbmc$barcode = sapply (rownames (meta_pbmc), function(x) unlist (strsplit (x, '_'))[3])
head (meta_pbmc$barcode)
tail (meta_pbmc$barcode)

# Split TCR barcodes by pool map barcodes and split by sample ####
tcrL_pbmc = lapply (vdj.dirs, function(x) read.csv(x))
#names (tcrL_pbmc) = unique (meta_pbmc$hash_pool)
names (tcrL_pbmc) = c('1_1_2_batch2','1_2_2_batch2')
tcrL_pbmc = lapply (names(tcrL_pbmc), function(x)
  {
  tcr_hashed = tcrL_pbmc[[x]][tcrL_pbmc[[x]]$barcode %in% meta_pbmc$barcode[meta_pbmc$pool_batch == x],]
  tcr_hashed$sampleID = meta_pbmc$sampleID[match (tcr_hashed$barcode, meta_pbmc$barcode)]
  tcr_hashed
  })
tcrL_pbmc = do.call (rbind, tcrL_pbmc)
tcrL_pbmc = split (tcrL_pbmc, tcrL_pbmc$sampleID)
tcrL_pbmc = tcrL_pbmc[sapply (tcrL_pbmc, nrow) > 1]
names (tcrL_pbmc) = paste0(names(tcrL_pbmc) ,'_pbmc')

# Import TCR of tumors ####
data.path_P7_tumor = '../../PM_scRNA_atlas/data/P7_filtered_contig_annotations.csv'
data.path_P2_tumor = '../../PM_scRNA_atlas/data/P2_filtered_contig_annotations.csv'
data.path_P9_tumor = '../../PM_scRNA_atlas/data/P9_filtered_contig_annotations.csv'

vdj.dirs = c (data.path_P7_tumor, data.path_P2_tumor, data.path_P9_tumor)
tcrL_tumor = lapply (vdj.dirs, function(x) read.csv(x))
names (tcrL_tumor) = c('P7_tumor','P2_tumor','P9_tumor')

contig_list = c(tcrL_pbmc, tcrL_tumor)
#contig_list = contig_list[sapply (contig_list, nrow) >0]
combinedTCR <- combineTCR (contig_list, 
                samples = names (contig_list), 
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)
pdf ('Plots/FIGURE_S6B_clonotype_overlap2.pdf',width=5,3.5)
print (clonalOverlap(combinedTCR, 
              cloneCall = "strict", 
              method = "morisita"))
dev.off()

# Compare clonotypes
pdf (paste0('Plots/FIGURE_S6C_shared_clones.pdf'),width=4.5,3)
cs1 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P7_pbmc", 
              y.axis = "P7_tumor",
              dot.size = "total",
              graph = "proportion")
cs2 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P2_pbmc", 
              y.axis = "P2_tumor",
              dot.size = "total",
              graph = "proportion")
cs3 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P9_pbmc", 
              y.axis = "P9_tumor",
              dot.size = "total",
              graph = "proportion")
print (cs1)
print (cs2)
print (cs3)
dev.off()

srt_merged = combineExpression (combinedTCR, 
                         srt_merged, 
                         cloneCall="strict", 
                         group.by = "sample", 
                         proportion = FALSE,
                         #cloneSize = c(NonExpanded = 1, Expanded = Inf)
                         cloneSize = c(Single = 1, NonExpanded = 5, Expanded = Inf)
                         )

srt_merged = srt_merged[, !is.na(srt_merged$cloneSize)]


cnmf_spectra_unique_comb = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Tms_20"))

srt_merged = ModScoreCor (
        seurat_obj = srt_merged, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Tms', outdir = paste0('Plots/'))


### Find celltypes of shared clones and their exhaustion score ####
cs1 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P7_pbmc", 
              y.axis = "P7_tumor",
              dot.size = "total",
              exportTable = T,
              graph = "proportion")
rownames (cs1) = cs1$Var1
meta_P7 = srt_merged@meta.data[grep ('^P7', colnames(srt_merged)),]
colnames (cs1) = c('CTstrict2','pbmc','tumor','class','sum','pbmc_fraction','tumor_fraction')
meta_P7 = cbind (meta_P7, cs1[meta_P7$CTstrict,])

cs2 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P2_pbmc", 
              y.axis = "P2_tumor",
              dot.size = "total",
              exportTable = T,
              graph = "proportion")
rownames (cs2) = cs2$Var1
meta_P2 = srt_merged@meta.data[grep ('^P2', colnames(srt_merged)),]
colnames (cs2) = c('CTstrict2','pbmc','tumor','class','sum','pbmc_fraction','tumor_fraction')
meta_P2 = cbind (meta_P2, cs2[meta_P2$CTstrict,])

cs3 = clonalScatter(combinedTCR, 
              cloneCall ="strict", 
              x.axis = "P9_pbmc", 
              y.axis = "P9_tumor",
              exportTable = T,
              dot.size = "total",
              graph = "proportion")
if ("size" %in% colnames(cs3)){
  cs3$size <- NULL
}
rownames (cs3) = cs3$Var1
meta_P9 = srt_merged@meta.data[grep ('^P9', colnames(srt_merged)),]
colnames (cs3) = c('CTstrict2','pbmc','tumor','class','sum','pbmc_fraction','tumor_fraction')
meta_P9 = cbind (meta_P9, cs3[meta_P9$CTstrict,])

meta_P7$'NA' <- NULL
meta_P2$'NA' <- NULL

### Find expanded clones in both tumor and blood across samples ####
cs_df = do.call (rbind, list (meta_P7, meta_P2, meta_P9))
cs_df = cs_df[,!is.na(colnames(cs_df))]
cs_df = cs_df[,colnames (cs3)]
srt_merged@meta.data = cbind (srt_merged@meta.data, cs_df[match(colnames(srt_merged), rownames(cs_df)),])
srt_merged$site = ifelse (grepl('tumor',colnames(srt_merged)), 'tumor','pbmc')
dual = sapply (unique(srt_merged$CTstrict), function(x) all(c('pbmc','tumor') %in% srt_merged$site[srt_merged$CTstrict == x]))
names(dual)[dual]
ccomp_df = srt_merged@meta.data[srt_merged$CTstrict %in% names(dual)[dual], ]
table (ccomp_df$site)



# Compute mean exhaustion for each clonotype ####
exhaustion_module = 'Tm5' # define exhaustion module
clones_l = split(srt_merged@meta.data[,exhaustion_module], paste0(srt_merged$CTstrict, srt_merged$sampleID, srt_merged$site))

clones = unlist (lapply (clones_l, function(x) mean (x)))
ext_cut = summary (clones)['Mean']

ccomp_df = split (ccomp_df, paste0(ccomp_df$CTstrict, ccomp_df$sampleID, ccomp_df$site))
ccomp_df = do.call (rbind, unname(lapply (ccomp_df, function(x) data.frame (
  CTstrict = x$CTstrict[1],
  exhaustion = mean(x[,exhaustion_module]), 
  class = x$class[1], 
  sampleID = x$sampleID[1],
  tumor_fraction = x$tumor_fraction[1],
  pbmc_fraction = x$pbmc_fraction[1],
  sum = x$sum[1],
  site = x$site[1]))))

# Plot dual expanded clonotypes ####
ccomp_df2 = ccomp_df[ccomp_df$class == 'dual.expanded',]
ccomp_df2 = ccomp_df2[!is.na(ccomp_df2$class),]

ccomp_df3 = reshape(ccomp_df2, idvar = "CTstrict", timevar = "site", direction = "wide")
ccomp_df3$ex_mean = rowMeans  (ccomp_df3[,c('exhaustion.pbmc','exhaustion.tumor')])
sp = ggplot (ccomp_df3, aes (
  x = exhaustion.pbmc, 
  y = exhaustion.tumor, 
  color = ex_mean, 
  size= sum.tumor#, 
  ))+
geom_point(shape=19) +
gtheme + 
paletteer::scale_color_paletteer_c("ggthemes::Red") +

pdf ('Plots/FIGURE_5O_exhausted_shared_clonotypes2.pdf', height = 2.6,3.5)
print (sp)
dev.off()




### FIGURE 6C ####
gcdata_TIC = readRDS('TICAtlas.rds') # Please download from: https://zenodo.org/records/5186413#%20.YRqbJC1h2v6

table(gcdata_TIC@meta.data$cell_type)
table(gcdata_TIC@meta.data$source)

# Subset for NK cells ####
gcdata_TIC_nk = gcdata_TIC[,gcdata_TIC$cell_type == 'NK']
table (gcdata_TIC_nk$cell_type,gcdata_TIC_nk$source)

# Import NK from MPM
srt_nk = srt[, srt$celltype %in% c('FGFBP2_NK','KLRC1_NK')]

srt_nk_merged = merge (gcdata_TIC_nk, srt_nk)
srt_nk_merged = NormalizeData (object = srt_nk_merged, normalization.method = "LogNormalize", scale.factor = 10000)
srt_nk_merged$source[is.na(srt_nk_merged$source)] = 'PM'
srt_nk_merged$subtype[is.na(srt_nk_merged$subtype)] = 'PM'

gd = geneDot (
srt_nk_merged,
gene = 'KLRC1',
x = 'subtype', # 
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
returnDF = T)

gd$groups = factor (gd$groups, levels = gd$groups[order (-gd$percent)])
col_pal = setNames (c('red',rep('grey',13)), levels(gd$groups))
bp = ggplot (gd, aes (x = groups, y = percent)) +
geom_bar(aes (fill = groups), position="dodge", stat="identity") + 
scale_fill_manual (values = col_pal) + NoLegend() + 
gtheme


png (paste0('Plots/KLRC1_positive.png'), width=900,height=700, res=300)
print (bp)
dev.off()

