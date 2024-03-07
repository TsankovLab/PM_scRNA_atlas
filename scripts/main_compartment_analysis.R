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


#### FIGURE 3G - Nichenet analysis in PLVAP+ endothelial vs rest ####
library (nichenetr)
receiver_cells = 'PLVAP+_EC'
srt$celltype_nichenet = as.character (srt$celltype_simplified)
srt$celltype_nichenet[srt$celltype == 'PLVAP'] = 'PLVAP+_EC'

sender_cells = unique(srt$celltype_nichenet)[!unique(srt$celltype_nichenet) %in% receiver_cells]

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

Idents(srt) = srt$celltype_nichenet
DE_table_receiver = FindMarkers (object = srt, ident.1 = 'PLVAP+_EC', ident.2 = 'Endothelial', min.pct = 0.10) %>% rownames_to_column("gene")

## receiver
receiver = "PLVAP+_EC"
Idents(srt) = srt$celltype_nichenet
expressed_genes_receiver = get_expressed_genes (receiver, srt, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
list_expressed_genes_sender = sender_cells %>% unique() %>% lapply(get_expressed_genes, srt, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

DE_table_receiver = FindMarkers (object = srt, ident.1 = 'PLVAP+_EC', ident.2 = 'Endothelial', min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 1e-10 & abs(avg_log2FC) >= 0.5) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities (geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands = ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

#srt$celltype_simplified
gd2 = geneDot (
  seurat_obj = srt,
  gene = factor (rev(best_upstream_ligands), levels = best_upstream_ligands),
  x = 'celltype_simplified2',
  #z = srt@meta.data[,'patient'],
  min_expression = 1,
  facet_ncol = 9,
  swap_axes = F,
  scale.data=T,
  plotcol = palette_gene_expression2) + gtheme_italic_y
pdf ('Plots/FIGURE_3G_nichenet_ligands_sender_cells_dotplot2.pdf', height=3, width=4.5)
gd2 
dev.off()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target = vis_ligand_target[, colnames(vis_ligand_target) %in% head (DE_table_receiver$gene,30)]
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot( "Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1, face= 'italic')) + 
scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))


pdf ('Plots/FIGURE_3G_nichenet_target_genes.pdf', width=5, height=4)
p_ligand_target_network + theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1, face = "italic"),axis.text.y = element_text(face = "italic"))
dev.off()


# ligand activity heatmap
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") 
# ligand expression Seurat dotplot
order_ligands_adapted = str_replace_all (order_ligands, "\\.", "-")
rotated_dotplot = DotPlot(srt[,srt$celltype_nichenet %in% sender_cells] %>% subset(celltype_nichenet %in% sender_cells), features = order_ligands_adapted) + 
coord_flip() + scale_color_gradient2 (palette_gene_expression2) + theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, face='italic'),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  ) +
theme(legend.text = element_text(size = 10), 
  legend.title = element_text(size = 12)) 

figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 45,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  #p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face='italic'), legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(3, 10, 10))

legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 3,
    align = "h", rel_widths = c(1.5, 2, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,10), nrow = 2, align = "hv")
pdf ('Plots/FIGURE_3G_nichenet_combineplot.pdf', width=8, height=5)
combined_plot
dev.off()


### FIGURE 4F - dotplot CXCLs and receptor across celltypes ####
srt$celltype_simplified3 = as.character(srt$celltype_simplified)
srt$celltype_simplified3[srt$celltype %in% c('Mono_CD16','Mono_CD14','TAMs')] = srt$celltype[srt$celltype %in% c('Mono_CD16','Mono_CD14','TAMs')]
#srt$celltype_simplified3[srt$celltype %in% c('pDC')] = 'DCs'
srt$celltype_simplified3[srt$celltype_simplified3 %in% c('T_cells')] = srt$celltype[srt$celltype_simplified3 %in% c('T_cells')]
srt$celltype_simplified3 = factor (srt$celltype_simplified3, levels = rev(c('Malignant','Mesothelium','Glia','Alveolar','Fibroblasts','Endothelial','SmoothMuscle','Mono_CD16','Mono_CD14','TAMs','Mast','DCs','CD8','CD4','TFH','Tregs','NK','B_cells','Plasma','pDCs')))
cxcl_markers = c('CXCL9','CXCL10','CXCL11','CXCR3')
pdf ('Plots/FIGURE_4F_cxcls_and_receptor_expression2.pdf', height=4.5, width=4)
  geneDot (
  seurat_obj = srt,
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
dev.off()  

### FIGURE 4C - VISTA expression ####
ext_markers='VSIR'
low_celltypes = c('Alveolar','LEC','Mesothelium','Glia')
srt2 = srt[,!srt$celltype %in% low_celltypes]
srt2 = srt2[,srt2$sampleID != 'P10']
srt2$celltype = gsub ('_','',srt2$celltype)
ext_avg = AverageExpression (srt2, features = ext_markers, group.by = c('sampleID','celltype'))[[1]]
rownames (ext_avg) = ext_markers
ext_avg = log2(as.data.frame (t(ext_avg))+1)
ext_avg$group = rownames (ext_avg)

ext_avg = gather (ext_avg, gene, avg_expression, 1:(ncol(ext_avg)-1))
ext_avg$sample = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[1])
ext_avg$celltype = sapply (ext_avg$group, function(x) unlist(strsplit(x, '_'))[2])
ext_avg$SE_group = setNames(srt2$SE_group, srt2$sampleID)[ext_avg$sample]
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
boxl
dev.off()

# FIGURE 6B - get NK activators and inhibitors and cross with CPDB database ####
### load cellphoneDB results ####
LR_df = read.csv ('../../PM_scRNA_atlas/data/LR_NK_CPDB.csv')

metaGroupNames = c('celltype_simplified2','sampleID')
NK_inhib_sample = AverageExpression (srt, features = LR_df$nk, group.by = metaGroupNames)
NK_inhib_sample = log2(as.data.frame (t(NK_inhib_sample[[1]]))+1)
NK_inhib_sample = NK_inhib_sample[,LR_df$nk]
NK_inhib_sample = NK_inhib_sample[grep ('^NK_', rownames(NK_inhib_sample)),, drop=F]
NK_inhib_sample$sample = rownames (NK_inhib_sample)
NK_inhib_sample$sample = gsub ('NK_','',NK_inhib_sample$sample)
NK_inhib_sample = gather (NK_inhib_sample,LR, expression, 1:(ncol (NK_inhib_sample) - 1))
nk_gene_order = NK_inhib_sample$LR
NK_inhib_sample = do.call (rbind, lapply (split(NK_inhib_sample, NK_inhib_sample$LR), function(x) data.frame (NKmean = median(x$expression), NKsd = sd(x$expression))))
NK_inhib_sample$rec = rownames (NK_inhib_sample)
NK_inhib_sample$rec = gsub ('\\.\\d','',NK_inhib_sample$rec)

Mal_HLA_sample = AverageExpression (srt, features = LR_df$mal, group.by = metaGroupNames)
Mal_HLA_sample = log2(as.data.frame (t(Mal_HLA_sample[[1]]))+1)
Mal_HLA_sample = Mal_HLA_sample[,LR_df$mal]
Mal_HLA_sample = Mal_HLA_sample[grep ('Malignant', rownames(Mal_HLA_sample)),, drop=F]
Mal_HLA_sample$sample = rownames (Mal_HLA_sample)
Mal_HLA_sample = gather (Mal_HLA_sample,LR, expression, 1:(ncol (Mal_HLA_sample) - 1))
mal_gene_order = Mal_HLA_sample$LR
Mal_HLA_sample = do.call (rbind, lapply (split(Mal_HLA_sample, Mal_HLA_sample$LR), function(x) data.frame (MALmean = median(x$expression), MALsd = sd(x$expression))))
Mal_HLA_sample$lig = rownames (Mal_HLA_sample)
Mal_HLA_sample$lig = gsub ('\\.\\d','',Mal_HLA_sample$lig)

Mal_NK_df = cbind(NK_inhib_sample[unique(nk_gene_order),], Mal_HLA_sample[unique(mal_gene_order),])
Mal_NK_df$label = paste0(Mal_NK_df$rec, ':',Mal_NK_df$lig)
Mal_NK_df$label2 = ifelse (Mal_NK_df$label == 'KLRC1:HLA-E', 'hit','nohit')
sp = ggplot(data = Mal_NK_df, aes(x = NKmean, y = MALmean)) +
  geom_errorbar(
    aes(ymin = MALmean - MALsd, ymax = MALmean + MALsd),
    width = 0.1,color= 'grey', linewidth=0.3,
    position = position_dodge(0.1)
  ) +
  geom_errorbar(
    aes(xmin = NKmean - NKsd, xmax = NKmean + NKsd),
    width = 0.1, color= 'grey',linewidth=0.3,
    position = position_dodge(0.1)
  ) +
  geom_point(aes (color=label2)) +
  labs(title = "NK Inhibitors",
       x = "NK",
       y = "Malignant") + 
  scale_color_manual (values = c(hit = 'red',nohit='grey22')) + 
  geom_text_repel (data = Mal_NK_df, aes(label = label), size=1.4, fontface='italic', segment.size=.1, max.overlaps=10) +
  NoLegend() + gtheme_no_rot


pdf ('Plots/FIGURE_6B_NK_MAL_ligand_receptor_scatterplot.pdf', 3,width=2.7)
sp
dev.off()

