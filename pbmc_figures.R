use UGER
conda activate scrnatools 
R

set.seed(1234)
options(warn = 1)

projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction/scRNA/pbmc/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)
source ('../../../../scripts/scrna_pipeline/useful_functions.R')
source ('../../../../scripts/scrna_pipeline/load_libraries.R')
source ('../../../../scripts/scrna_pipeline/ggplot_aestetics.R')
source ('../../../../scripts/projects/meso_prj/meso_naive_RNA/MPM_naive_13_pallettes.R')

# Load scS-score
scs_sample_avg = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/malignant_stromal_subset/no_harmony/malignant_subset/no_harmony/scs_score_per_sample.csv', row.names=1)

# Load Seurat object
srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_PBMC_CITEseq2_analysis/_cellranger_filtered_Filter_400_1000_10/sampleID_harmony_cc_nCount_RNA_regressed/srt.rds')

#sample_ID = c('PRJ203-ANWO01_7004-1419'='p4','PRJ203-ANWO01_7008-1420'='p8','PRJ203-ANWO01_7009-1421' = 'p9', p786 = 'p786', p811 = 'p811',p826 = 'p826',p846 ='p846', p848 = 'p848')
sampleID = c(p786 = 'P1',p13 = 'P13',p846 = 'P3', p12 = 'P12', p7 = 'P6', p8 = 'P2', p848 = 'P5', p4 = 'P7', p11 = 'P11', p811 = 'P4', p826 = 'P8', p9 = 'P9')
srt$sampleID = unname (sampleID[as.character (srt$sampleID3)])
srt$sampleID = factor (srt$sampleID, levels = names(palette_sample)[names(palette_sample) %in% unique (srt$sampleID)])


# load palettes
scrna_palettes2 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_PBMC_CITEseq2_analysis/palettes.rds')
scrna_palettes = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/palettes.rds')
for (i in seq_along (scrna_palettes2)) assign (names(scrna_palettes)[i], scrna_palettes[[i]])
for (i in seq_along (scrna_palettes)) assign (names(scrna_palettes)[i], scrna_palettes[[i]])
palette_gene_expression = c("#002F70FF", 'white',"#5F1415FF")
palette_protein_expression = c(low="darkblue",mid= "white",high= "darkgreen") 
palette_feature_RNA = c('lightgrey',"#5F1415FF")
palette_feature_protein = c("lightgrey", "darkgreen")

# set reduction name
reductionName = 'sampleID_harmony_umap'

# Set ggplot theme aestetics
gtheme = theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face='italic'),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

gtheme_prot = theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

### FIGURE S1A ####
ccomp_df = as.data.frame (table(srt$sampleID3))
srt$nFeature_RNAL = log2 (srt$nFeature_RNA)
srt$nCount_RNAL = log2 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=F, group.by = 'sampleID4',pt.size = 0, ncol = 3) 
vln_p = lapply (vln_p, function(x) x +
scale_fill_manual(values = palette_sample))

png (paste0("Plots/FIGURE_S1A_QC_nFeat_nCount_m.percent_vlnPlot.png"), 3300, 1000, res=300)
print (vln_p[[1]] | vln_p[[2]] | vln_p[[3]])
dev.off()

### FIGURE 1C ####
pq1 = DimPlot(srt, reduction = reductionName, group.by = "celltype_simplified", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt), 'cells')) +
scale_color_manual(values = celltype_simplified_palette)
pq2 = DimPlot(srt, reduction = reductionName, group.by = "sampleID3", label = F, label.size = 3 ,repel = TRUE) + ggtitle (paste(ncol(srt), 'cells')) +
scale_color_manual(values = palette_sample)

png (paste0("Plots/FIGURE_1C_celltypes_samples_umap.png"), width = 2800, height = 1000, res=300)
pq1 + pq2
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
  x = 'sampleID3', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = 'predicted.celltype.l2',
  z = NULL, 
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_gene_expression) +
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
  x = 'sampleID3', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = 'predicted.celltype.l2',
  z = NULL, 
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  x_name ='samples',
  y_name = 'celltype',
  plotcol = palette_protein_expression) +
  gtheme_prot
  

png(paste0("Plots/FIGURE_S1G_DotPlot_topADTmarkers_predicted.celltype.l2.png"), width = 3500, height = 2200, res = 300)
dp 
dev.off()

#####################
### TCR ANALYSIS ####
#####################
library (scRepertoire)
#### Load PBMC TCR-sequencing from biopsy samples ####
meta_pbmc = srt@meta.data[srt$batch == 'batch2',]
meta_pbmc$sampleID3 = as.character (meta_pbmc$sampleID3)
data.path_pbmc1 = '/ahg/regevdata/projects/lungCancerBueno/10x/prj203-mtsinai-tsankov-pbmcs/20230511/processed/2023-03-29-1-1/version=1/per_sample_outs/cellranger_multi_run/vdj_t/filtered_contig_annotations.csv'
data.path_pbmc2 = '/ahg/regevdata/projects/lungCancerBueno/10x/prj203-mtsinai-tsankov-pbmcs/20230511/processed/2023-03-29-1-2/version=1/per_sample_outs/cellranger_multi_run/vdj_t/filtered_contig_annotations.csv'
vdj.dirs = c(data.path_pbmc1, data.path_pbmc2)

# Add hashing pools to seurat metadata ####
#meta_pbmc$hash_pool = sapply (meta_pbmc$pool_batch, function(x) unlist(strsplit (x, '\\-'))[2])
meta_pbmc$hash_pool = gsub ('_batch2','',meta_pbmc$pool_batch)

# Check barcodes per pool are unique ####
lapply (split (meta_pbmc$barcode,meta_pbmc$hash_pool), function(x) length(x) == length(unique(x)))

# Split TCR barcodes by pool map barcodes and split by sample ####
tcrL_pbmc = lapply (vdj.dirs, function(x) read.csv(x))
names (tcrL_pbmc) = unique (meta_pbmc$hash_pool)
tcrL_pbmc = lapply (names(tcrL_pbmc), function(x)
  {
  tcr_hashed = tcrL_pbmc[[x]][tcrL_pbmc[[x]]$barcode %in% meta_pbmc$barcode[meta_pbmc$hash_pool == x],]
  tcr_hashed$sampleID4 = meta_pbmc$sampleID3[match (tcr_hashed$barcode, meta_pbmc$barcode)]
  tcr_hashed
  })
tcrL_pbmc = do.call (rbind, tcrL_pbmc)
tcrL_pbmc = split (tcrL_pbmc, tcrL_pbmc$sampleID4)
names (tcrL_pbmc) = paste0(names(tcrL_pbmc) ,'_pbmc')

# Load PBMC TCR-sequencing from resection samples ####
data.path1 = '/ahg/regevdata/projects/lungCancerBueno/10x/immunai_transfer_2020-03-20/processed/30-472307723_2021-01-05-2-1_TCR_all_contig_annotations.csv'
data.path2 = '/ahg/regevdata/projects/lungCancerBueno/10x/immunai_transfer_2020-03-20/processed/30-472307723_2021-01-05-2-2_TCR_all_contig_annotations.csv'
vdj.dirs = c(data.path1, data.path2)
meta_pbmc_batch1 = srt@meta.data[srt$batch == 'batch1',]
meta_pbmc_batch1$sampleID3 = as.character (meta_pbmc_batch1$sampleID3)

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
  tcr_hashed$sampleID4 = meta_pbmc_batch1$sampleID3[match (tcr_hashed$barcode, meta_pbmc_batch1$barcode)]
  tcr_hashed
  })
tcr_df_batch1 = do.call (rbind, tcrL_pbmc_batch1)
tcr_df_batch1 = split (tcr_df_batch1, tcr_df_batch1$sampleID4)
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

unique_barcode = !duplicated(paste0(srt$sampleID3, '_pbmc_', srt$barcode))
srt_unique = srt[,unique_barcode]
srt_unique = RenameCells(
    srt_unique,
    new.names = paste0(srt_unique$sampleID3, '_pbmc_', srt_unique$barcode))

srt_tcr = combineExpression (combinedTCR, srt_unique,
  proportion = F,
#  cloneSizes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = Inf))
  cloneSize = c(Single = 1, NonExpanded = 5, Expanded = Inf))
  #)
clonesize_name = setNames (c('Single','NonExpanded','Expanded','None'), c('Single (0 < X <= 1)','NonExpanded (1 < X <= 5)','Expanded (5 < X <= Inf)', 'None ( < X <= 0)'))
srt_tcr = srt_tcr[, !is.na (srt_tcr$cloneSize)]
srt_tcr$cloneSize = unname(clonesize_name[as.character(srt_tcr$cloneSize)])
# Set palette
#clonotype_levels = c('Small (1e-04 < X <= 0.001)','Medium (0.001 < X <= 0.01)','Large (0.01 < X <= 0.1)')
#clonotype_levels = c('NonExpanded (0 < X <= 1)','Expanded (1 < X <= Inf)', NA)
clonotype_levels = c('Single','NonExpanded','Expanded', 'None', NA)

#palette_clonotype = list(palette_clonotype = setNames (c(as.character (paletteer::paletteer_d("trekcolors::romulan2",9))[c(5,1)],'grey55'),c(clonotype_levels)))

# FIGURE S5I - Show expanded clonotypes are found in CD8 cells ####  
dp = DimPlot (srt_tcr, reduction = reductionName, group.by='cloneSize') + scale_color_manual (values = palette_clonotype)

pdf(paste0('Plots/FIGURE_S5_Expanded_vs_nonExpanded_umap.pdf'), height=3,width=5)
dp
dev.off()

metaGroupName = 'cloneSize'
cc_box = cellComp (
  seurat_obj = srt_tcr, 
  metaGroups = c('predicted.celltype.l2',metaGroupName),
  plot_as = 'bar',
  prop = FALSE,
  pal = palette_clonotype,
  subset_prop = 'Expanded',  
  #pal = viridis::mako(,
  facet_ncol = 15
  ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
cc_box$data$predicted.celltype.l2 = factor (cc_box$data$predicted.celltype.l2, levels = as.character(cc_box$data$predicted.celltype.l2[order(-cc_box$data$Freq)]))

pdf(paste0('Plots/FIGURE_S5I_Expanded_vs_nonExpanded_barplot.pdf'),width=6,3)
cc_box
dev.off()

### Subset to CD8 ####
#projdir_T = paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scRNA_meso/', 'Tcells/')
#dir.create (paste0(projdir_T, 'Plots/'), recursive=T)
srt_tcr = srt_tcr[, srt_tcr$predicted.celltype.l1 %in% c('CD8 T')]
srt_tcr = srt_tcr[, !srt_tcr$sampleID4 %in% c('p846')] # remove low counts samples

cnmf_t = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/TNK_cells_subset/sampleID2_harmony/T_only_subset/sampleID2_harmony/Tm_cnmf_combined.rds')

DefaultAssay (srt_tcr) = 'RNA'
srt_tcr = ModScoreCor (
        seurat_obj = srt_tcr, 
        geneset_list = cnmf_t, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Tms_', outdir = paste0(projdir,'Plots/'))

# FIGURE 6E - Show expanded clonotypes have higher exhaustion ####  
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm5, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID4'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('Single','NonExpanded','Expanded'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID4', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Expanded','Single')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

pdf ('Plots/FIGURE_6E_expanded_exhaustion.pdf',4.3,width = 2)
box2
dev.off()



### Check with MHC II module ####
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm4, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID4'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('Single','NonExpanded','Expanded'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID4', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Expanded','Single')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

box3 = ggpaired (ccomp_df[!ccomp_df$sampleID4 %in% c('p9','p11','p4'),], x = "cloneSize", y = "x", id = 'sampleID4', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box3 = box3 + stat_compare_means (paired = TRUE, comparisons = list(c('Expanded','Single')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box3 = box3 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

box4 = ggpaired (ccomp_df[ccomp_df$sampleID4 %in% c('p9','p11','p4'),], x = "cloneSize", y = "x", id = 'sampleID4', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 

box4 = box4 + stat_compare_means (paired = TRUE, comparisons = list(c('Expanded','Single')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box4 = box4 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

pdf ('Plots/FIGURE_6B_expanded_MHCII.pdf',4.3,width = 2)
box2
box3
box4
dev.off()


### Check with Cytotoxic module ####
ccomp_df = as.data.frame (srt_tcr@meta.data)
ccomp_df = ccomp_df[,!duplicated (colnames(ccomp_df))]
ccomp_df = aggregate (ccomp_df$Tm2, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID4'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('Single','NonExpanded','Expanded'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID4', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Expanded','Single')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

# box2$layers[[1]]$aes_params$alpha =  .7
# box2$layers[[1]]$geom_params$linewdith =  .1
# box2$layers[[1]]$aes_params$lwd =  .1
pdf ('Plots/FIGURE_6B_expanded_cytox.pdf',4.3,width = 2)
box2
dev.off()


# ### muscat DS on genotype ####
srt_orig = srt
srt = srt_tcr
force = FALSE
do.fgsea = TRUE
logfcThreshold = .5
pvalAdjTrheshold = 0.05
ds_method = "DESeq2" #c("edgeR", "DESeq2", "limma-trend", "limma-voom")
metaGroupName1 = 'sampleID4'
metaGroupName2 = 'cloneSize'
metaGroupName3 = 'SE_group2'
#metaGroupName3 = 'treatment'
#muscatIdents = c('CB2_PTENL','GFP')
#muscatIdents = c('CB2_PTENL','CB2_PTENL_C124S')
srt$SE_group2 = ifelse (srt$SE_group == 'S-High','Shigh','Slow')
muscatIdents = c('Shigh','Slow')
pbDS_min_cells = 2
topGenes = 20 
#show_genes = c('Cxcr3','Cxcr4','Ccr3', 'Ccr1', 'Cxcl10','Cxcl11','Ccl12','Ccl9','Ccl5','Ccl8','Ccl6') # check genes
#srt2 = srt
#srt = subset (srt, sampleID2_harmony_snn_res.1 %in% c(0,1,2,3,4,5,6,7))
org='human'
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/DS_muscat.R')

top_pathways = Inf
fgseaRes_flt = fgseaResAll[[1]][[1]]
fgseaRes_flt = split (fgseaRes_flt, fgseaRes_flt$NES > 0)
fgseaRes_flt = lapply (fgseaRes_flt, function(x) head (x[order(x$padj),], top_pathways))
fgseaRes_flt = do.call (rbind, fgseaRes_flt)
fgseaRes_flt$padj_log_signed = -log10 (fgseaRes_flt$padj) * sign (fgseaRes_flt$NES)
fgseaRes_flt = fgseaRes_flt[order (fgseaRes_flt$NES),]
fgseaRes_flt$pathway = factor (fgseaRes_flt$pathway, levels = unique (fgseaRes_flt$pathway))
fgseaRes_flt$direction = as.character(sign(fgseaRes_flt$NES))
fgseaRes_flt$direction = ifelse (fgseaRes_flt$direction == '1', 'S-High', 'E-High')
fgseaRes_flt$pathway = factor (fgseaRes_flt$pathway, levels = fgseaRes_flt$pathway[order (fgseaRes_flt$padj_log_signed)])
fgseaResAll_dp = ggplot (fgseaRes_flt, aes (x = padj_log_signed, y = pathway, fill = direction)) + 
  geom_bar (stat = 'identity') +
  theme_classic () +
  scale_fill_manual(values = palette_SE_group) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf (paste0(projdir_ms,'Plots/fGSEA_annotation_',paste(gmt_annotations,collapse='_'),'_barplots.pdf'), width = 12, height=15)
print(fgseaResAll_dp)
dev.off()

vp = VlnPlot (srt, features = unlist(fgseaRes_flt[fgseaRes_flt$pathway == 'GO_ADAPTIVE_IMMUNE_RESPONSE','leadingEdge']), group.by= 'sampleID4')
pdf ('Plots/check_adaptive_response_genes.pdf',height=10,10)
vp
dev.off()

### Compare DEG with TNK clonotypes ####
tumor_deg = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/TNK_cells_subset/sampleID2_harmony/muscat_sampleID4_celltype_SE_group2_Shigh_vs_Ehigh_method_DESeq2/DEGresults.rds')
tumor_deg = tumor_deg[['pch_filtered']][[1]]
colnames (tumor_deg) = paste0('tumor_',colnames (tumor_deg))
pbmc_deg = tbl_df
int_genes = intersect (tumor_deg$tumor_gene, pbmc_deg$gene)

int_deg = cbind (tumor_deg[match(pbmc_deg$gene[pbmc_deg$gene%in%int_genes], tumor_deg$tumor_gene),], pbmc_deg[pbmc_deg$gene %in% int_genes,])
int_deg = int_deg [int_deg$tumor_baseMean > summary (int_deg$tumor_baseMean)['Mean'],]
cor (int_deg[c('avg_log2FC', 'tumor_logFC')], method = 'spearman')
sp = ggplot (int_deg, aes (x= avg_log2FC, y = tumor_logFC)) + geom_point()

pdf ('Plots/tumor_pbmc_deg_comparison.pdf')
sp
dev.off()



vp = VlnPlot (srt_tcr_cd8_exp[,srt_tcr_cd8_exp$exhausted == 'exhausted'], features = rownames(srt_tcr_cd8_exp)[grep('HLA', rownames(srt_tcr_cd8_exp))], group.by= 'sampleID4')
pdf ('Plots/check_adaptive_response_genes.pdf',height=10,10)
vp
dev.off()

vp = VlnPlot (srt_tcr_cd8_exp[,srt_tcr_cd8_exp$exhausted == 'exhausted'], features = c('CTLA4','HAVCR2','LAG3','TIGIT','PDCD1'), group.by= 'sampleID4')
pdf ('Plots/check_adaptive_response_genes.pdf',height=10,10)
vp
dev.off()



