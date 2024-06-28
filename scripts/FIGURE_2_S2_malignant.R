# Malignant compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/malignants/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Import seurat object
srt_tumor = readRDS ('../GSE190597_srt_tumor.rds')
srt = srt_tumor[, srt_tumor$celltype == 'Malignant']
srt$sampleID = factor (srt$sampleID, levels = levels (srt$sampleID)[1:12])

#### Import bulk RNA subtype signatures to define malignant score ####

# Bueno ####
top_bueno_genes = readRDS ('../../PM_scRNA_atlas/data/bueno_molecular_subtype_deg.rds')
if (!all (names(top_bueno_genes)  %in% colnames(srt@meta.data))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = top_bueno_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Bueno', outdir = paste0(projdir,'Plots/'))
}

# Import Blum S score gene signature ####
top_blum_genes = readRDS ('../../PM_scRNA_atlas/data/blum_se_score_genes.rds')
top_blum_genes = lapply (top_blum_genes, function(x) head (x, 20))
names (top_blum_genes)= paste0(names(top_blum_genes), '_20')
if (!all (names(top_blum_genes) %in% colnames(srt@meta.data))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = top_blum_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Blum', outdir = paste0(projdir,'Plots/'))
}

# Load cnmfs ####
cnmf_spectra_unique_comb = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Cms_20"))
cnmf_spectra_unique_comb_full = as.list (read_excel( "../../PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Cms_full"))

# re-compute nmf score after combining ####
if (!all (names(cnmf_spectra_unique_comb) %in% colnames(srt@meta.data))) {
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = cnmf_spectra_unique_comb, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
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
  NoLegend()
  
pdf(paste0('Plots/FIGURE_2B_',sarc_nmf,'_signatures_boxplot.pdf'),width=3,2) #width = 10, height = 11,
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
  geom_point(alpha=1, shape=21, stroke=.25, size=1) +
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
pdf ('Plots/FIGURE_S2A_neftel_diagram_on_malignant_cells_per_sample.pdf',10,3)
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
    
pdf ('Plots/FIGURE_2A_neftel_diagram_on_malignant_cells.pdf',3,3)
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
            gtheme_no_rot + NoLegend()
      
pdf (paste0('Plots/FIGURE_S2E_sarcomatoid_vs_S_score.pdf'),height=3,width=3)
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
ptable_factor = c(2),
returnDF= T,
prop=T) 
cp = do.call (rbind, lapply (split (cp, cp$Cm_r_max), function(x) {x$Freq = proportions(x$Freq); x}))
cp = ggplot (cp, aes_string (x= 'Cm_r_max', y= 'Freq')) +
    geom_bar(position="stack", stat="identity", aes_string(fill= 'sampleID')) +              
    scale_fill_manual (values= palette_sample) +
    gtheme
cp$data$Cm_r_max = factor (cp$data$Cm_r_max, levels = names (cnmf_spectra_unique_comb[row_order(hm)])) 
pdf (paste0('Plots/FIGURE_2C_sample_abundance_cNMF_stacked_barplot.pdf'))
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
bulk_stat = readRDS ('../../PM_scRNA_atlas/data/bueno_tcga_Cms_stats.rds')

sarc_nmf = 'Cm17'
for (st in study)
  {
  bulk_stat_study = bulk_stat[bulk_stat$study == st,]
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
  
  
  pdf (paste0('Plots/FIGURE_2D_bulkRNA_',st,'_on_CNMFs_combined_3.pdf'), width=7, height=4)
  print (bp)
  dev.off()
  }


# FIGURE S2E - Check nmfs overlap with blum gene sets ####
blumL = readRDS ('../../PM_scRNA_atlas/data/blum_se_score_genes.rds')

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

ha = HeatmapAnnotation (' ' = sarc_score_hm, col = list (' ' = palette_sample_cnv), which='row',simple_anno_size = unit(.3, "cm"),  heatmap_legend_param = list(
legend_direction = "horizontal"))
hm = Heatmap (
  left_annotation = ha,
  ccomp_df_cor_df, col=palette_cnv_fun, 
  border=T, 
  column_names_rot = 45,
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(
legend_direction = "horizontal"))
pdf (paste0('Plots/FIGURE_2F_RNA_CNV_Cm_heatmap.pdf'),width=5,height=1.7)
print (draw (hm, heatmap_legend_side = "bottom"))
dev.off()




