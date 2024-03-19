# FIGURE 6 - NK compartment analysis ####

# Set seeds
set.seed(1234)

# Set option to convert errors to warnings to 1
options(warn = 1)

# Set project directory
projdir = 'scRNA/NK/' # define project directory
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

source ('../../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../../PM_scRNA_atlas/scripts/R_utils.R')
source ('../../PM_scRNA_atlas/scripts/palettes.R')
source ('../../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

# Load Seurat object
srt = readRDS ('../srt_tumor.rds')

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
  geom_text_repel (data = Mal_NK_df, aes(label = label), size=1.8, fontface='italic', segment.size=.1, max.overlaps=10) +
  NoLegend() + gtheme_no_rot


pdf ('Plots/FIGURE_6B_NK_MAL_ligand_receptor_scatterplot.pdf', 3,width=2.7)
print (sp)
dev.off()

### FIGURE 6C - comparison of KLRC1-expressing NK across human cancers ####
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


png (paste0('Plots/FIGURE_6C_KLRC1_positive.png'), width=900,height=700, res=300)
print (bp)
dev.off()


# FIGURE 6E - Make barplots of IFN and degranulation response of NK following NKG2A blockade ####
degran = read.csv ("../../PM_scRNA_atlas/data/NK_Degranulation_NKpercent.csv")
degran$treatment[degran$treatment == 'PBS'] = 'vehicle'
degran$treatment = factor (degran$treatment, levels = c('vehicle','rIFNg'))
degran$group = gsub ('a','A',degran$group)
ifn_prd = read.csv ("../../PM_scRNA_atlas/data/NK_IFN_NKpercent.csv")
ifn_prd$treatment[ifn_prd$treatment == 'PBS'] = 'vehicle'
ifn_prd$treatment = factor (ifn_prd$treatment, levels = c('vehicle','rIFNg'))
ifn_prd$group = gsub ('a','A',ifn_prd$group)

sel_groups = c('NK','NK_NKG2A','NK_MHC1','NK_NKG2A_MHC1')
ifn_prd = ifn_prd[ifn_prd$group %in% sel_groups,]
degran = degran[degran$group %in% sel_groups,]

degran$group = factor (degran$group, levels = sel_groups)
ifn_prd$group = factor (ifn_prd$group, levels = sel_groups)
degran_bar <- degran %>%
  group_by (treatment, group) %>%
  summarize (Mean = mean(percentage),
            SEM = sd(percentage) / sqrt(n()))
ifn_prd$group = factor (ifn_prd$group, levels = sel_groups)
ifn_bar <- ifn_prd %>%
  group_by (treatment, group) %>%
  summarize (Mean = mean(percentage),
            SEM = sd(percentage) / sqrt(n()))


palette_nkg2a = setNames (as.character(paletteer::paletteer_d("colorBlindness::Blue2Gray8Steps",8)[c(2,1,7,8)]), sel_groups)
palette_nkg2a = setNames (c('grey55','brown2','grey22','brown3'),sel_groups)

box2 = ggplot() +
        geom_errorbar (data = degran_bar, aes(x = group, ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.3, color = "black",linewidth=.3) +
        geom_bar (data = degran_bar, aes(x = group, y = Mean, color = group), fill='white',stat = "identity", alpha=1) +
        geom_bar (data = degran_bar, aes(x = group, y = Mean, color = group),fill='grey55', stat = "identity", alpha=0.1) +
        geom_point (data = degran, aes (x = group, y = percentage), position='identity', alpha=.7, color="grey44", size=1.2) +
        scale_fill_manual (values = palette_nkg2a) + 
        scale_color_manual (values = palette_nkg2a) +
        geom_line (data = degran, aes(x = group, y = percentage, group = cell_line), color='grey44',linewidth=.2, alpha=.7) +
        facet_wrap (~treatment) + gtheme #+

stat.test2 <- degran |>
  group_by(treatment) |>
  rstatix::t_test(reformulate ('group', 'percentage'), 
  comparisons = list(c('NK','NK_NKG2A'),
    c('NK_MHC1','NK_NKG2A_MHC1')), paired=T) |>
  rstatix::add_xy_position (x = "group", fun = "max", dodge = 0.8) |> 
  rstatix::adjust_pvalue (method = "none") |>
  rstatix::add_significance()

box2 <- box2 + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
  bracket.nudge.y = -5, hide.ns = TRUE,
  label = "p.adj.signif")


png ('Plots/FIGURE_6E_degran_percentages_barplots.png',800,width = 1200, res=300)
box2
dev.off()



box2 = ggplot() +
        geom_errorbar (data = ifn_bar, aes(x = group, ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.3, color = "black",linewidth=.3) +
        geom_bar (data = ifn_bar, aes(x = group, y = Mean, color = group), fill='white',stat = "identity", alpha=1) +
        geom_bar (data = ifn_bar, aes(x = group, y = Mean, color = group),fill='grey55', stat = "identity", alpha=0.1) +
        geom_point (data = ifn_prd, aes (x = group, y = percentage), position='identity', alpha=.7, color="grey44", size=1.2) +
        scale_fill_manual (values = palette_nkg2a) + 
        scale_color_manual (values = palette_nkg2a) + 
        geom_line (data = ifn_prd, aes(x = group, y = percentage, group = cell_line), color='grey44',linewidth=.2, alpha=.7) +
        facet_wrap (~treatment) + gtheme #+

stat.test2 <- ifn_prd |>
  group_by(treatment) |>
  rstatix::t_test(reformulate ('group', 'percentage'), 
  comparisons = list(c('NK','NK_NKG2A'),
    c('NK_MHC1','NK_NKG2A_MHC1')), paired=T) |>
  rstatix::add_xy_position (x = "group", fun = "max", dodge = 0.8) |> 
  rstatix::adjust_pvalue (method = "none") |>
  rstatix::add_significance()

box2 <- box2 + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = TRUE,
  label = "p.adj.signif")


png ('Plots/FIGURE_6E_ifn_percentages_barplots.png',800,width = 1200, res=300)
box2
dev.off()



### Generate plots of boolean gating ####
bool = read.csv ('../../PM_scRNA_atlas/data/NK_Degranulation_or_IFN_Boolean_percent.csv')


sel_groups = c('NK','NK_NKG2A','NK_MHC1','NK_NKG2A_MHC1')
bool = bool[bool$group %in% sel_groups,]

bool$group = factor (bool$group, levels = sel_groups)
bool_bar <- bool %>%
  group_by (treatment, group) %>%
  summarize (Mean = mean(percentage),
            SEM = sd(percentage) / sqrt(n()))

box2 = ggplot() +
        geom_errorbar (data = bool_bar, aes(x = group, ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.3, color = "black",linewidth=.3) +
        geom_bar (data = bool_bar, aes(x = group, y = Mean, color = group), fill='white',stat = "identity", alpha=1) +
        geom_bar (data = bool_bar, aes(x = group, y = Mean, color = group),fill='grey55', stat = "identity", alpha=0.1) +
        geom_point (data = bool, aes (x = group, y = percentage), position='identity', alpha=.7, color="grey44", size=1.2) +
        scale_fill_manual (values = palette_nkg2a) + 
        scale_color_manual (values = palette_nkg2a) + 
        geom_line (data = bool, aes(x = group, y = percentage, group = cell_line), color='grey44',linewidth=.2, alpha=.7) +
        facet_wrap (~treatment) + gtheme #+

stat.test2 <- bool |>
  group_by(treatment) |>
  rstatix::t_test(reformulate ('group', 'percentage'), 
  comparisons = list(c('NK','NK_NKG2A'),
    c('NK_MHC1','NK_NKG2A_MHC1')), paired=T) |>
  rstatix::add_xy_position (x = "group", fun = "max", dodge = 0.8) |> 
  rstatix::adjust_pvalue (method = "none") |>
  rstatix::add_significance()

box2 <- box2 + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
  bracket.nudge.y = -5, hide.ns = TRUE,
  label = "p.adj.signif")


png ('Plots/FIGURE_S6G_boolean_gating_barplots.png',800,width = 1200, res=300)
box2
dev.off()
