# Main compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/main/' # define project directory
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Load Seurat object
srt = readRDS ('../srt_tumor.rds')
srt_pbmc = readRDS ('../srt_pbmc.rds')

### FIGURE 1D ####
top_markers = c('KRT19','CALB2','SLPI','HP','ITLN1','CLDN1','PMP2','VGF','OLIG2','SFTPC','SFTPB','COL1A1','COL3A1','PECAM1','PLVAP','VWF','ACTA2','MYL9','MYH11','LYZ','CD14','C1QA','CD3D','CD3E','CD8A','NKG7','GNLY','GZMA','CD79A','IGHM','CD37','IGLC2','IGHA1','IGLC3','IRF8','IRF4','LILRA4')
srt$celltype_simplified2 = factor (srt$celltype_simplified2, levels = c('Malignant','Mesothelium','Glia','Alveolar','Fibroblasts','Endothelial','SmoothMuscle','Myeloid','T_cells','NK','B_cells','Plasma','pDC'))

dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = factor (top_markers, levels = top_markers),
  x = 'celltype_simplified2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  x_name ='samples',
  y_name = 'celltype',
  swap_axes=T,
  plotcol = palette_gene_expression2) +
  gtheme_italic

pdf ('Plots/FIGURE_1D_top_markers_celltype_expression.pdf', height=4.8, width=8.5)
dp + theme(
    axis.text.x = element_text(face = "italic")
  )
dev.off()  

### FIGURE S1A ####
ccomp_df = as.data.frame (table(srt$sampleID))
srt$nFeature_RNAL = log2 (srt$nFeature_RNA)
srt$nCount_RNAL = log2 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=F, group.by = 'sampleID',pt.size = 0, ncol = 3) 
vln_p = lapply (vln_p, function(x) x +
scale_fill_manual(values = palette_sample))

png ("Plots/FIGURE_S1A_QC_nFeat_nCount_m.percent_vlnPlot.png", 3300, 1000, res=300)
print (vln_p[[1]] | vln_p[[2]] | vln_p[[3]])
dev.off()
 
### FIGURE S1B - TF expression across celltypes ####
celltype_order =c('Malignant','Mesothelium','Fibroblasts','SmoothMuscle','Endothelial','Alveolar','Glia','Myeloid','T_cells','NK','B_cells','pDC','Plasma')
DefaultAssay (srt) = 'RNA'

TF.use <- c(
  'TEAD1','MEIS2',
  'WT1','OSR1',
  'TWIST1','SNAI2',
  'HES4','TBX2',
  'SOX7','SOX18',
  'IRX3','IRX5',
  'ASCL1','OLIG1',    
  'CEBPB','KLF6',
  'BCL11B','TSC22D3',
  'ID2','ETS1',
  'POU2F2','MEF2B',
  'IRF8','IRF4',
  'XBP1','CREB3L2'
  )
cluster.averages = AverageExpression(object = srt[TF.use,], group.by = 'celltype_simplified2', return.seurat = F)
cluster.averages = log2(as.data.frame (cluster.averages[[1]])+1)
cluster.averages = cluster.averages[TF.use, celltype_order]
pdf("Plots/FIGURE_S1C_heatmap.pdf", width =3,height=3, useDingbats=T)
Heatmap (
  t(scale (t(cluster.averages))), 
  cluster_rows=F,
  row_names_side = 'left',
  cluster_columns=F, 
  col = palette_gene_expression_fun(scale (t(cluster.averages))),
  column_names_rot = 45,
  row_names_gp = gpar(fontface = "italic",fontsize = 6),
  column_names_gp = gpar (fontsize = 6))
dev.off()

# FIGURE S1C Composition plots by sampling ####
metaGroupNames = c('celltype_simplified','sampleID')
srt$immune = ifelse (srt$celltype_simplified %in% c('B_cells','T_cells','MonoMac','cDCs','NK','Plasma','Mast','pDC'), 'immune','non-immune')

metaGroupNames = c('sampleID','celltype_simplified2','sampling','immune')
ccc_bar2 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'box',
  pal = palette_sampling,
  prop = TRUE,
  ptable_factor = c(1,4),
  #subset_prop = 'cycling',
  facet_ncol = 6,
  returnDF = T,
  facet_scales = 'free'
  )

### Run Diricthlet regression and extract pvalues ####
drc_res = list()
for (i in levels (ccc_bar2$immune))
  {
  cc_box1 = spread(ccc_bar2[ccc_bar2$immune == i,], key = celltype_simplified2, value = Freq)
  tmp = cc_box1[,4:ncol(cc_box1)]
  tmp[is.na(tmp)] = 0
  AL = DR_data(tmp)
  res = DirichReg (AL ~ as.factor (sampling), cc_box1)
  
  # Code from Smillie UC paper: https://github.com/cssmillie/ulcerative_colitis/blob/40a6846565401804d5e2b08e82b52c06d12f0518/analysis.r#L247
  u = summary (res)
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('sampling', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  pvals = as.data.frame (pvals)
  pvals$sampling = rownames (pvals)
  pvals = gather (pvals, celltype, p, 1:(ncol(pvals)-1))
  pvals$sampling = gsub ('as.factor\\(\\)','',pvals$sampling)
  pvals$group1 = 'biopsy'
  pvals$group2 = pvals$sampling
  pvals$p_adj = p.adjust (pvals$p, method = 'fdr')
  pvals$p_sig = ''
  pvals$p_sig = ifelse (pvals$p_adj <= 0.05, ifelse (pvals$p_adj <= 0.01,ifelse (pvals$p_adj <= 0.001,'***','**'),'*'),'ns')
  drc_res[[i]] = pvals
  }

  drc_res_df = do.call (rbind, drc_res)
  #Get y_positions for pval annotation
  ccc_bar2 = split (ccc_bar2, ccc_bar2$celltype_simplified2)
  y_max = do.call (rbind, lapply(ccc_bar2, function(x) 
    {
    tpm = boxplot(x$Freq ~ x$sampling)$stats
    tpm = max (tpm[5,])
    tmp = data.frame (celltype = x$celltype_simplified2[1], y.position = tpm)
    }))
  drc_res_df = drc_res_df[match (as.character(y_max$celltype), drc_res_df$celltype),]
  drc_res_df$y.position = y_max$y.position    
  drc_res_df$immune = gsub ('\\.\\d','',rownames (drc_res_df))

metaGroupNames = c('sampleID','celltype_simplified2','sampling','immune')
srt$celltype_simplified2 = as.character(srt$celltype_simplified2)

ccc_bar3 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'box',
  pal = palette_sampling,
  prop = TRUE,
  ptable_factor = c(1,4),
  #subset_prop = 'cycling',
  facet_ncol = 6,
  returnDF = F,
  facet_scales = 'free'
  ) + gtheme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))#+ NoLegend() + ggtitle (i)

ccc_bar4 = ccc_bar3 + stat_pvalue_manual(
    drc_res_df, 
    step.increase = 0.05, 
    y.position = 'y.position',
    x = "celltype", hide.ns=T, 
    color='grey20',
    label = "p_sig",
     vjust = 0.5,
    tip.length = 0.02
    )
    
pdf (paste0('Plots/FIGURE_S1C_sampling_celltype_fractions_diritchlet.pdf'),width=6, height=2.5)
ccc_bar4
dev.off()


# inferCNV
infcnv = readRDS ('../../PM_scRNA_atlas/data/PM_infercnv.results.obj.Rds')
plot_cnv (infcnv)




