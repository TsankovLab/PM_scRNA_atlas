set.seed(1234)
options(warn = 1)

projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction/scRNA/pbmc_tnk/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)
source ('../../../../scripts/scrna_pipeline/useful_functions.R')
source ('../../../../scripts/scrna_pipeline/load_libraries.R')
source ('../../../../scripts/scrna_pipeline/ggplot_aestetics.R')
source ('../../../../scripts/projects/meso_prj/meso_naive_RNA/MPM_naive_13_pallettes.R')

# Load scS-score
scs_sample_avg = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/malignant_stromal_subset/no_harmony/malignant_subset/no_harmony/scs_score_per_sample.csv', row.names=1)

# Load Seurat object
srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_PBMC_CITEseq2_analysis/_cellranger_filtered_Filter_400_1000_10/sampleID_harmony_cc_nCount_RNA_regressed/TNK_cells_subset/sampleID_harmony_cc_nCount_RNA_regressed/srt.rds')

# set reduction name
reductionName = 'sampleID_harmony_umap'



#sample_ID = c('PRJ203-ANWO01_7004-1419'='p4','PRJ203-ANWO01_7008-1420'='p8','PRJ203-ANWO01_7009-1421' = 'p9', p786 = 'p786', p811 = 'p811',p826 = 'p826',p846 ='p846', p848 = 'p848')
srt$sampleID = factor (srt$sampleID3, levels = names(palette_sample)[names(palette_sample) %in% unique (srt$sampleID3)])


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
  


#####################
### TCR ANALYSIS ####
#####################
library (scRepertoire)
#### Load PBMC TCR-sequencing from biopsy samples ####
meta_pbmc = srt@meta.data[srt$batch == 'batch2',]
meta_pbmc$sampleID = as.character (meta_pbmc$sampleID)
data.path_pbmc1 = '/ahg/regevdata/projects/lungCancerBueno/10x/prj203-mtsinai-tsankov-pbmcs/20230511/processed/2023-03-29-1-1/version=1/per_sample_outs/cellranger_multi_run/vdj_t/filtered_contig_annotations.csv'
data.path_pbmc2 = '/ahg/regevdata/projects/lungCancerBueno/10x/prj203-mtsinai-tsankov-pbmcs/20230511/processed/2023-03-29-1-2/version=1/per_sample_outs/cellranger_multi_run/vdj_t/filtered_contig_annotations.csv'
vdj.dirs = c(data.path_pbmc1, data.path_pbmc2)

# Add hashing pools to seurat metadata ####
#meta_pbmc$hash_pool = sapply (meta_pbmc$pool_batch, function(x) unlist(strsplit (x, '\\-'))[2])
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
data.path1 = '/ahg/regevdata/projects/lungCancerBueno/10x/immunai_transfer_2020-03-20/processed/30-472307723_2021-01-05-2-1_TCR_all_contig_annotations.csv'
data.path2 = '/ahg/regevdata/projects/lungCancerBueno/10x/immunai_transfer_2020-03-20/processed/30-472307723_2021-01-05-2-2_TCR_all_contig_annotations.csv'
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

srt_tcr = combineExpression (combinedTCR, srt_unique,
  proportion = F,
#  cloneSizes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = Inf))
  cloneSize = c(NonExpanded = 1, Small = 5, Large = Inf))
  #)
clonesize_name = setNames (c('NonExpanded','Small','Large','None'), c('NonExpanded (0 < X <= 1)','Small (1 < X <= 5)','Large (5 < X <= Inf)', 'None ( < X <= 0)'))
srt_tcr = srt_tcr[, !is.na (srt_tcr$cloneSize)]
srt_tcr$cloneSize = unname(clonesize_name[as.character(srt_tcr$cloneSize)])


# FIGURE S5H - Show expanded clonotypes are found in CD8 cells ####  
metaGroupName = 'cloneSize'
cc_box = cellComp (
  seurat_obj = srt_tcr, 
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
srt_tcr = srt_tcr[, srt_tcr$predicted.celltype.l1 %in% c('CD8 T')]
srt_tcr = srt_tcr[, !srt_tcr$sampleID %in% c('P3')] # remove low counts samples

library (readxl)
cnmf_t = as.list (read_excel( "../../cnmf_per_compartment.xlsx", sheet = "Tms_20"))

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
ccomp_df = aggregate (ccomp_df$Tm5, by=as.list(srt_tcr@meta.data[,c('cloneSize','sampleID'),drop=F]), mean)

ccomp_df$cloneSize = factor (ccomp_df$cloneSize, levels = c('NonExpanded','Small','Large'))
box2 = ggpaired (ccomp_df, x = "cloneSize", y = "x", id = 'sampleID', width = .5,
         fill = 'white', color='white', line.color = "gray", line.size = 0.3,
         palette = palette_clonotype) 
box2 = box2 + stat_compare_means (paired = TRUE, comparisons = list(c('Large','NonExpanded')),
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_nt + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_6N_expanded_exhaustion.png',height = 1000,width = 600, res=300)
box2
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
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_nt + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_S6J_expanded_MHCII.png',height = 1000,width = 600, res=300)
box2
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
    tip.length=0.02, method='t.test', label = "p.signif", bracket.nudge.y = -0.3) + gtheme_nt + NoLegend()
box2 = box2 + 
geom_point(position='identity', alpha=.7, color="grey44", size=1.2) +
geom_boxplot (aes_string(fill='cloneSize'),color = 'grey22', width=.5, alpha = 0.7, lwd=.2, outlier.shape = NA)

png ('Plots/FIGURE_S6J_expanded_cytotoxic.png',height = 1000,width = 600, res=300)
box2
dev.off()
