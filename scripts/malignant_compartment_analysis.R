set.seed(1234)

projdir = 'scRNA/malignants/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

scripts_dir = '../../scripts/'
source (paste0(scripts_dir,'palettes.R'))
source (paste0(scripts_dir,'ggplot_aestetics.R'))
source (paste0(scripts_dir,'R_utils.R'))
source (paste0(scripts_dir,'R_libraries.R'))

# Import seurat object
#srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/malignant_stromal_subset/no_harmony/malignant_subset/no_harmony/srt.rds')
srt = readRDS ('../srt.rds')
srt = srt[, srt$celltype == 'Malignant']

#### Import bulk RNA subtype signatures to define malignant score ####

# Bueno ####
top_bueno_genes = readRDS ('../bueno_molecular_subtype_deg.rds')
if (!all (colnames(srt@meta.data) %in% names(top_bueno_genes))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = top_bueno_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Bueno', outdir = paste0(projdir,'Plots/'))
}

# Import Blum S score gene signature ####
top_blum_genes = readRDS ('../blum_se_score_genes.rds')
top_blum_genes = lapply (top_blum_genes, function(x) head (x, 20))
if (!all (colnames(srt@meta.data) %in% names(top_blum_genes))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = top_blum_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Blum', outdir = paste0(projdir,'Plots/'))
}

rerun_cnmf = TRUE
if (run_cnmf)
  {

  #### Run cNMF ####
  cnmf_out = paste0('cNMF_normalized/cNMF_5_30_vf3000')
  library (scran)
  dir.create (paste0(cnmf_out,'/Plots/'), recursive=T)
  sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
  rowData=rownames(srt)) 
  sce = modelGeneVar(sce)
  
  # remove batchy genes
  batchy_genes = c('RPL','RPS','MT-')
  sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
  vf = getTopHVGs(sce, n=3000)
  
  count_mat = t(srt@assays$RNA@counts[vf,])
  norm_mat = t(srt@assays$RNA@data[vf,])
  write.table (count_mat, 'cnmf/counts_nmf_3000.txt', sep='\t', col.names = NA)
  write.table (norm_mat, 'cnmf/norm_nmf_3000.txt', sep='\t', col.names = NA)
  
  system ('cnmf prepare --output-dir ./ --name cnmf -c cnmf/counts_nmf_3000.txt --tpm cnmf/norm_nmf_3000.txt  -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 --n-iter 100 --seed 14 --numgenes 3000')
  system ("qsub -cwd -b y -l h_vmem=2g,h_rt=3:00:00 -o . -e . -N cnmf -t 1-100 'source /broad/software/scripts/useuse; use Anaconda; source activate /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf; cnmf factorize --output-dir . --name cnmf --worker-index $SGE_TASK_ID'")
  system ('cnmf combine --output-dir . --name cnmf')
  system ('cnmf k_selection_plot --output-dir ./ --name cnmf')
  system ('cnmf consensus --output-dir ./ --name cnmf --components 25 --local-density-threshold 0.3 --show-clustering')
  
  # Read in NMF results ####
  cnmf_spectra = read.table (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_25.dt_0_3.consensus.txt'))
  
  # Assign genes uniquely to cNMF modules based on spectra values
  cnmf_spectra = t(cnmf_spectra)
  max_spectra = apply (cnmf_spectra, 1, which.max)
  
  cnmf_spectra_unique = lapply (1:ncol(cnmf_spectra), function(x) 
        {
        tmp = cnmf_spectra[names(max_spectra[max_spectra == x]),x,drop=F]
        tmp = tmp[order(-tmp[,1]),,drop=F]
        head(rownames(tmp),top_nmf_genes)
        })
  names(cnmf_spectra_unique) = paste0('CN',seq_along(cnmf_spectra_unique))
  
  cnmf_spectra_unique_full = lapply (1:ncol(cnmf_spectra), function(x) 
        {
        tmp = cnmf_spectra[names(max_spectra[max_spectra == x]),x,drop=F]
        tmp = tmp[order(-tmp[,1]),,drop=F]
        rownames(tmp)
        #head(rownames(tmp),top_nmf_genes)
        })
  
  names(cnmf_spectra_unique_full) = paste0('CN',seq_along(cnmf_spectra_unique))
  
  write.csv (patchvecs (cnmf_spectra_unique_full), paste0(cnmf_out,'/cnmf_list_',k_selection,'_unique.csv'))
  write.csv (patchvecs (cnmf_spectra_nonunique_full), paste0(cnmf_out,'/cnmf_list_',k_selection,'_nonunique.csv'))
  
  
  srt = ModScoreCor (
          seurat_obj = srt, 
          geneset_list = cnmf_spectra_unique, 
          cor_threshold = NULL, 
          pos_threshold = NULL, # threshold for fetal_pval2
          listName = 'Cms', outdir = paste0(projdir,'Plots/'))
  
    
  # combine redundant cNMFs ####
  combine_nmf = c('CN1','CN2_12','CN3','CN4','CN5','CN6_10','CN7',
    'CN8_9','CN8_9','CN6_10','CN11','CN2_12','CN13','CN14','CN15','CN16','CN17_21',
    'CN18','CN19_24','CN20','CN17_21','CN22','CN23','CN19_24','CN25')
  
  cnmf_spectra_unique_comb = split (cnmf_spectra_unique, combine_nmf)
  cnmf_spectra_unique_comb = lapply (cnmf_spectra_unique_comb, function(x) if (length(x) > 1) lapply (x, function(y) head(y, 20)) else x)
  cnmf_spectra_unique_comb = lapply  (cnmf_spectra_unique_comb, function(x) unlist(x))
  
  cnmf_spectra_unique_full_comb = split (cnmf_spectra_unique_full, combine_nmf)
  cnmf_spectra_unique_full_comb = lapply  (cnmf_spectra_unique_full_comb, function(x) unlist(x))
  
  # change module numbering for paper ####
  old_cnmf_numbers = names (cnmf_spectra_unique_comb)
  new_cnmf_numbers = c(
    CN1 = 'Cm1',
    CN2_12 = 'Cm2',
    CN3 = 'Cm3',
    CN4 = 'Cm4',
    CN5 = 'Cm5',
    CN6_10 = 'Cm6',
    CN7 = 'Cm7',
    CN8_9 = 'Cm8',
    CN11 = 'Cm9',
    CN13 = 'Cm10',
    CN14 = 'Cm11',
    CN15 = 'Cm12',
    CN16 = 'Cm13',
    CN17_21 = 'Cm14',
    CN18 = 'Cm15',
    CN19_24 = 'Cm16',
    CN20 = 'Cm17',
    CN22 = 'Cm18',
    CN23 = 'Cm19',
    CN25 = 'Cm20')
  
  names (cnmf_spectra_unique_comb) = new_cnmf_numbers[names (cnmf_spectra_unique_comb)]
  names (cnmf_spectra_unique_full_comb) = new_cnmf_numbers[names (cnmf_spectra_unique_full_comb)]
  
  saveRDS (cnmf_spectra_unique_comb, 'Cm_cnmf_combined.rds')
  write.csv (patchvecs (lapply (cnmf_spectra_unique_comb, unname)),'Cm_cnmf_combined.csv')
  saveRDS (cnmf_spectra_unique_full_comb, 'Cm_cnmf_combined_full.rds')
  write.csv (patchvecs (lapply (cnmf_spectra_unique_full_comb, unname)),'Cm_cnmf_combined_full.csv')
  } else {
  library (readxl)
  cnmf_spectra_unique_comb = as.list (read_excel( "../../cnmf_per_compartment.xlsx", sheet = "Cms_20"))
  cnmf_spectra_unique_comb = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_13s_analysis/cellbender/_cellranger_raw_Filter_400_1000_25/sampling_harmony/malignant_stromal_subset/no_harmony/malignant_subset/no_harmony/Cm_cnmf_combined.rds')
  cnmf_spectra_unique_comb_full = as.list (read_excel( "../../cnmf_per_compartment.xlsx", sheet = "Cms_full"))
  }

# re-compute nmf score after combining ####
if (!all (colnames(srt@meta.data) %in% names(cnmf_spectra_unique_comb))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))
}

sarc_nmf = 'Cm17'
ccomp_df = srt@meta.data[,names (cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]

### FIGURE 2B - Plot boxplots of sarcomatoid cnmf ####
ccomp_df = srt@meta.data
ccomp_df = aggregate (ccomp_df[,c(sarc_nmf,'S_score_20')], by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]

ccomp_df = srt@meta.data
box = ggplot (ccomp_df, aes_string (x= 'sampleID', y= sarc_nmf)) +
  geom_violin (trim=TRUE, aes_string (fill = 'sampleID'),size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7) +
  geom_boxplot (aes_string(fill = 'sampleID'),
    linewidth = .2,
    width=0.2,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.3, alpha=0.7
     ) +
  gtheme +
  scale_fill_manual (values= palette_sample) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
png(paste0('Plots/FIGURE_2B_',sarc_nmf,'_signatures_boxplot.png'),1200,600, res=300) #width = 10, height = 11,
print (box)
dev.off()


#### FIGURE 2A - Generate Neftel diagram using the four subtypes from bueno ####
max_sarc = pmax(srt$Sarcomatoid, srt$`Biphasic-S`)
max_epit = pmax(srt$Epithelioid, srt$`Biphasic-E`)
srt$y_axis <- log2(abs(max_sarc - max_epit) + 1)
srt$y_axis[max_epit > max_sarc] <- -1 * srt$y_axis[max_epit > max_sarc]

srt$x_axis = 0
srt$x_axis[srt$y_axis > 0] = log2(abs(srt$Sarcomatoid - srt$`Biphasic-S`) + 1)[srt$y_axis > 0]
srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`] <- -1 * srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`]

srt$x_axis[srt$y_axis < 0] = log2(abs(srt$Epithelioid - srt$`Biphasic-E`) + 1)[srt$y_axis < 0]
srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid] <- -1 * srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid]
srt$bueno_color = 0
srt$bueno_color[srt$x_axis < 0 & srt$y_axis < 0] = 'Biphasic-E'
srt$bueno_color[srt$x_axis > 0 & srt$y_axis > 0] = 'Biphasic-S'
srt$bueno_color[srt$x_axis < 0 & srt$y_axis > 0] = 'Sarcomatoid'
srt$bueno_color[srt$x_axis > 0 & srt$y_axis < 0] = 'Epithelioid'

p2l = lapply (levels (srt$sampleID), function(x) 
  {    
  df = srt@meta.data[srt$sampleID == x,]
  tot_cells = nrow(df)    
  ggplot(df, aes(x=x_axis, y=y_axis, fill = bueno_color), color='black') +
  geom_point(alpha=1, shape=21, stroke=.25, linewidth=3,size=1) +
  xlim (c(-2,2)) + ylim (c(-2,2)) +
  geom_vline(xintercept = 0,linetype = 'dashed', size=.1) +
  geom_hline(yintercept = 0,linetype = 'dashed', size=.1) +
  scale_fill_manual(values = palette_bulk) +
  annotate("text", x = -1.1, y = 1.8, label = paste0("Sarco ",round(sum(df$x_axis < 0 & df$y_axis > 0) / tot_cells *100,1),'%') , size=3.5) +
  annotate("text", x = 1.2, y = 1.8, label = paste0("Bi-S ",round(sum(df$x_axis > 0 & df$y_axis > 0)/ tot_cells *100,1),'%'), size=3.5) +
  annotate("text", x = -1.2, y = -1.8, label = paste0("Bi-E ",round(sum(df$x_axis < 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=3.5) +
  annotate("text", x = 1.2, y = -1.8, label = paste0("Epit ",round(sum(df$x_axis > 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=3.5) + 
  xlab('') +
  ylab('') + 
  theme_void() + 
  NoLegend()
  })
png ('Plots/FIGURE_S2A_neftel_diagram_on_malignant_cells_per_sample.png',3200,1300, res=300)
print (wrap_plots (p2l, ncol=6))
dev.off()

# FIGURE S2A - Do it across samples ####
ccomp_df = as.data.frame (srt@meta.data)
tot_cells = nrow(ccomp_df)
p2l = ggplot(ccomp_df, aes(x=x_axis, y=y_axis, fill=bueno_color)) +
  geom_point(alpha=.8, shape=21, size=.8,color = 'black', stroke=0.1) +
  xlim (c(-1.2,1.2)) + ylim (c(-1.5,1.5)) +
  geom_vline(xintercept = 0,linetype = 'dashed', size=0.3, color='blue') +
  geom_hline(yintercept = 0,linetype = 'dashed', size=0.3, color='blue') +
  scale_fill_manual(values = palette_bulk) +
  ggtitle ('Bueno subtypes') + 
  gtheme_no_rot +
  NoLegend()
    
png ('Plots/FIGURE_2A_neftel_diagram_on_malignant_cells.png',1000,1000, res=300)
print (p2l)
dev.off()

# FIGURE S2E - Show correlation of scs with blum score ####
ccomp_df = srt@meta.data
ccomp_df = aggregate (ccomp_df[,c(sarc_nmf,'S_score_20')], by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
sp = ggscatter(ccomp_df, x = sarc_nmf, y = 'S_score_20',
  add.params = list(color = "blue", fill = "lightgray", linewidth=0.3,stroke=0.4, size=0.4),
            add = "reg.line", conf.int = TRUE,# label = rownames(ccomp_df),
            cor.coef = TRUE, cor.method = "spearman", alpha=0.7,
            xlab = sarc_nmf, ylab = "S_score Blum et al") + 
            geom_point (aes (color = rownames(ccomp_df)),alpha=.7) +
            scale_color_manual (values = palette_sample) +
            gtheme
      
pdf (paste0('Plots/FIGURE_S2E_sarcomatoid_vs_S_score.pdf'),height=3,width=5)
print (sp)
dev.off()


### FIGURE 2C - CNMF CORRELATION ACROSS SAMPLES ####
ccomp_df = srt@meta.data[,names (cnmf_spectra_unique_comb)]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1]
cor_mat = cor (ccomp_df, method = 'spearman')
#cor_mat = cor_mat[malig_modules,malig_modules]
cor_p = unlist (lapply (1:ncol(ccomp_df), function(y) cor.test (ccomp_df[,sarc_nmf], ccomp_df[,y], method = 'spearman', alternative='two.sided')$p.value))
names (cor_p) = colnames (ccomp_df)

sarc_nmf = 'Cm17'
scs_score_hm = as.data.frame(cor_mat)[,sarc_nmf]
names (scs_score_hm) = rownames (cor_mat)

# triangle heatmap
ha = HeatmapAnnotation(
  S_score = scs_score_hm, 
  col = list (S_score = palette_sample2_fun), 
  which='row',
  simple_anno_size = unit(.3, "cm"))
hm = draw (Heatmap (cor_mat, 
  left_annotation = ha,
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_2C_modules_cnmf_cor_heatmap_triangle.pdf'),5.5,4.5)
print (hm)
dev.off()

### Make stacked barplot of sample usage of cNMF to put aside of heatmap ####
cp = cellComp (
srt, 
metaGroups = c('Cm_r_max','sampleID'), 
plot_as = 'bar',
pal = palette_sample,
ptable_factor = 1,
prop=T) +
gtheme
cp$data$Cm_r_max = factor (cp$data$Cm_r_max, levels = names (cnmf_spectra_unique_comb[row_order(hm)])) 
pdf (paste0(projdir, 'Plots/FIGURE_2C_sample_abundance_cNMF_stacked_barplot.pdf'))
print (cp)
dev.off()


# FIGURE S2D - Plot cnmf module correlation across cells ####
sarc_nmf = 'Cm17'
ccomp_df = srt@meta.data[,names (cnmf_spectra_unique_comb)]
ccomp_df_cor = cor (ccomp_df, method='spearman')
library (circlize)

# Add track showing corelation of each module to sarcomatoid module
sarc_score_hm = as.data.frame(ccomp_df_cor)[[sarc_nmf]]

# triangle heatmap
ha = HeatmapAnnotation(S_score = sarc_score_hm, col = list (S_score = palette_sample2_fun), which='row')
hmc = draw (Heatmap (ccomp_df_cor,# row_km=15,
  left_annotation = ha,
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
pdf (paste0('Plots/FIGURE_S2D_modules_cnmf_cor_cells_heatmap.pdf'),5.5,4.5)
print (hmc)
dev.off()

# FIGURE S2C - Make dotplot of top markers for each nmf ####
cnmf_spectra_unique_comb_ordered = cnmf_spectra_unique_comb[row_order(hmc)]
marker_genes = unlist (lapply (cnmf_spectra_unique_comb_ordered, function(x) head (x, 5)))
srt$Cm_r_max2 = factor (srt$Cm_r_max, levels = rev(names(cnmf_spectra_unique_comb_ordered)))
dotp = geneDot (
seurat_obj = srt,
#gene = top_tfs2, 
gene = factor (marker_genes, levels = marker_genes),
x = 'Cm_r_max2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
y = NULL,
min_expression = 0,
facet_ncol = 5,
lim_expression = NULL,
scale.data=T,
x_name ='samples',
y_name = 'celltype',
swap_axes = T,
plotcol = palette_gene_expression2) + gtheme_italic

pdf ('Plots/FIGURE_S2C_cNMF_markers_expression.pdf', width=20)
print (dotp)
dev.off()


# FIGURE 2D - Check consistency with bulkRNA of cNMF modules correlation to scS_score ####
study  = c('bueno','tcga')
bulk_stat = readRDS ('bueno_tcga_Cms_stats.rds')

sarc_nmf = 'Cm17'
for (st in study)
  {
  bulk_stat_study = bulk_stat[bulk_stat$dataset == st,]
  # get cNMF correlation to S score across samples
  ccomp_df = srt@meta.data[,names (cnmf_spectra_unique_comb)]
  ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
  rownames(ccomp_df) = ccomp_df[,1]
  ccomp_df = ccomp_df[,-1]
  sarc_cor_cnmf = cor (ccomp_df[, sarc_nmf], ccomp_df, method='spearman')
  sarc_cor_cnmf = as.data.frame (t(cbind (sarc_cor = sarc_cor_cnmf, row.names = names(sarc_cor_cnmf))))
  colnames(sarc_cor_cnmf)[1] = 'sarc_cor'
  sarc_cor_cnmf = sarc_cor_cnmf[order (-sarc_cor_cnmf[,1]), , drop=F]
  sarc_cor_cnmf$cnmf = rownames(sarc_cor_cnmf)
  #sarc_cor_cnmf$cnmf_annotation = cnmf_ann[sarc_cor_cnmf$cnmf]
  sarc_cor_cnmf$bulk_padj = bulk_stat_study$p.adj.signif[match (sarc_cor_cnmf$cnmf, bulk_stat_study$mod)]
  sarc_cor_cnmf$bulk_padj[sarc_cor_cnmf$bulk_padj == 'ns'] = ''
  sarc_cor_cnmf$bulk_direction = bulk_stat_study$statistic[match (sarc_cor_cnmf$cnmf, bulk_stat_study$mod)]
  sarc_cor_cnmf$bulk_direction = ifelse (sign (sarc_cor_cnmf$bulk_direction) < 0, 'E-High','S-High')
  sarc_cor_cnmf$bulk_padj2 = sarc_cor_cnmf$bulk_padj
  sarc_cor_cnmf$bulk_padj[sarc_cor_cnmf$sarc_cor < 0] = '' 
  sarc_cor_cnmf$bulk_padj2[sarc_cor_cnmf$sarc_cor > 0] = '' 
  sarc_cor_cnmf$cnmf = factor (sarc_cor_cnmf$cnmf, levels = unique (sarc_cor_cnmf$cnmf))
  
  bp = ggplot(sarc_cor_cnmf, aes(fill=bulk_direction, y=sarc_cor, x=cnmf)) +  
      geom_bar(position="stack", stat="identity",alpha=.7) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      scale_fill_manual (values=c(palette_SE_group)) +
      geom_text(aes(label=bulk_padj), vjust=0.5, hjust = 0.5) +
      geom_text(aes(label=bulk_padj2), vjust=1, hjust=0.5)
  
  
  pdf (paste0('Plots/FIGURE_2D_bulkRNA_',st,'_on_CNMFs_combined_2.pdf'), width=7, height=4)
  print (bp)
  dev.off()
  }


# FIGURE S2E - Check nmfs overlap with blum gene sets ####
blumL = readRDS ('blum_se_score_genes.rds')

  mat_ov = ovmat (c(cnmf_spectra_unique_comb_full, blumL), 
  compare_lists = list(nmf = names(cnmf_spectra_unique_comb_full), 
    blum = names(blumL)),
    ov_threshold=1,
    palette = as.character(ov_pal))
pdf (paste0('Plots/FIGURE_S2E_blum_nmf_overlap.pdf'), width=2.8, height=4)
print (mat_ov)
dev.off()


#### Import CNV and correlate with cnmf modules scores ####
library (biomaRt)
library (GenomicFeatures)

library(GenomicFeatures)
getChromInfoFromUCSC("hg19") %>%
   head(22)

# generate metacells to run correlation on
library (hdWGCNA)
samples_to_remove = c('P13','P3','P1')
if (!file.exists ('metacells.rds'))
  {
  ### Try running cnmf on metacells
  library (hdWGCNA)
  srt_meta = srt  
    srt_meta = SetupForWGCNA(
    srt,
    gene_select = "custom", # the gene selection approach
    features = rownames (srt),  
    #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "hdWGCNA" # the name of the hdWGCNA experiment
  )
  srt_meta = srt_meta[, !srt_meta$sampleID %in% samples_to_remove]    
  srt_meta = MetacellsByGroups (
    seurat_obj = srt_meta,
    group.by = 'sampleID', # specify the columns in seurat_obj@meta.data to group by
    k = 50, # nearest-neighbors parameter
    max_shared = 30,
    reduction = 'pca',
    ident.group = 'sampleID', # set the Idents of the metacell seurat object
    min_cells = 10
  )
    
  # Number of metacells
  message (paste ('number of metacells:', ncol(GetMetacellObject(srt_meta))))
  # normalize metacell expression matrix:
  srt_meta = NormalizeMetacells (srt_meta)  
  metacells_obj = GetMetacellObject (srt_meta)    
  saveRDS (metacells_obj, 'metacells.rds')
  } else {
  metacells_obj = readRDS ('metacells.rds')  
  }

ext_region = 0
region = list (  
  chr1p = c(chr =1 , start = 1, end = 121884915),
  chr3p = c(chr=3, start = 1, end = 87693207),
  chr4 = c(chr = 4, start = 1, end = 190214555),
  chr13 = c(chr = 13, start=1, end = 115169878),
  chr14 = c(chr = 14, start=1, end = 107349540),
  chr22 = c(chr = 22, start = 1, end = 51304566)
  )
  
region_df = do.call (rbind, region)
# specify the database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.
res = cbind(
  region_df,
  mut_regions = apply(region_df, 1, function(i){
    x <- getBM(attributes=c("external_gene_name"), 
               filters = c("chromosome_name" , "start", "end"), 
               values = list(i[1], i[2], i[3]), 
               mart = ensembl)
    # return genes, comma separated
    paste(x$external_gene_name, collapse = ",")
  }))
  
# Compute score of CNVs  ####
mut_regions = lapply(1:nrow(region_df), function(g) strsplit(res[g ,'mut_regions'], '\\,'))
mut_regions = unlist(mut_regions, recursive=F)    

names(mut_regions) = rownames(region_df)
    metacells_obj = ModScoreCor (
        seurat_obj = metacells_obj, 
        geneset_list = mut_regions, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'mut', outdir = paste0(projdir,'Plots/'))

# Re-compute score of CNMF filtering out genes in CNV ####
cnmf_spectra_unique_comb_cnv = cnmf_spectra_unique_comb  
genes_in_cnv = unique (unlist(strsplit (unlist (res[,4]),',')))
cnmf_spectra_unique_comb_cnv = lapply ( cnmf_spectra_unique_comb_cnv, function(x) x[!x%in%genes_in_cnv])
names (cnmf_spectra_unique_comb_cnv) = paste0('cnv_',names (cnmf_spectra_unique_comb_cnv))    
    metacells_obj = ModScoreCor (
        seurat_obj = metacells_obj, 
        geneset_list = cnmf_spectra_unique_comb_cnv, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'CNV_Cms', outdir = paste0(projdir,'Plots/'))

# Compute correlation CNV scores vs CNMFs per sample ####
ccomp_df_cnmf_l = split (metacells_obj@meta.data[,names(cnmf_spectra_unique_comb_cnv)], metacells_obj$sampleID)
ccomp_df_mut_l = split (metacells_obj@meta.data[,names(mut_regions)], metacells_obj$sampleID)

#high_cnv_samples = c('p12','p7','p8','p848','p4')
high_cnv_samples = samples_to_remove
high_cnv_samples = levels (srt$sampleID)[!levels (srt$sampleID) %in% high_cnv_samples]
ccomp_df_cnmf_l = ccomp_df_cnmf_l[high_cnv_samples]  
ccomp_df_mut_l = ccomp_df_mut_l[high_cnv_samples]

ccomp_df_cor = lapply (names(ccomp_df_cnmf_l), function(x) cor (ccomp_df_cnmf_l[[x]], ccomp_df_mut_l[[x]], method='spearman'))

ccomp_df_cor = lapply (names(ccomp_df_cnmf_l), function(x) cor (ccomp_df_cnmf_l[[x]], ccomp_df_mut_l[[x]], method='spearman'))
ccomp_df_cor = lapply (seq_along(ccomp_df_cor), function(x)
    {
    tmp = ccomp_df_cor[[x]]
    tmp = as.data.frame (tmp)
    tmp$sample = names(ccomp_df_cnmf_l)[x]
    tmp$cNMF = rownames (tmp)
    tmp
    })

ccomp_df_cor_df = as.data.frame (do.call (rbind, ccomp_df_cor))
ccomp_df_cor_df = gather (ccomp_df_cor_df, region, correlation, 1:(ncol(ccomp_df_cor_df) - 2))
ccomp_df_cor_df = aggregate (ccomp_df_cor_df$correlation, by = ccomp_df_cor_df[,c('cNMF','region'), drop=F], 'median')
ccomp_df_cor_df = reshape(ccomp_df_cor_df, idvar = "region", timevar = "cNMF", direction = "wide")
rownames (ccomp_df_cor_df) = ccomp_df_cor_df[,1]
ccomp_df_cor_df = ccomp_df_cor_df[,-1]
colnames (ccomp_df_cor_df) = gsub ('x.cnv_','',colnames (ccomp_df_cor_df))

sarc_nmf = 'Cm17'
sarc_score_hm = ccomp_df_cor_df[,sarc_nmf]
names (sarc_score_hm) = rownames (ccomp_df_cor_df)

ha = HeatmapAnnotation (' ' = sarc_score_hm, col = list (' ' = palette_sample_cnv), which='row',simple_anno_size = unit(.3, "cm"))
hm = Heatmap (
  left_annotation = ha,
  ccomp_df_cor_df, col=palette_cnv_fun, 
  border=T, 
  column_names_rot = 45,
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7))
pdf (paste0('Plots/FIGURE_2F_RNA_CNV_Cm_heatmap.pdf'),width=5,height=1.7)
print (hm)
dev.off()
