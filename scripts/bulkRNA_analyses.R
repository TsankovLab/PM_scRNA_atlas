# BULK-RNA ANALYSIS 

set.seed(1234)

projdir = 'bulkRNA/'
system (paste('mkdir -p',paste0(projdir,'Plots/')))
setwd (projdir)

source ('../PM_scRNA_atlas/scripts/R_libraries.R')
source ('../PM_scRNA_atlas/scripts/palettes.R')
source ('../PM_scRNA_atlas/scripts/R_utils.R')
source ('../PM_scRNA_atlas/scripts/ggplot_aestetics.R')

blk_meta = readRDS ('../PM_scRNA_atlas/data/bueno_tcga_bulk_covariates.rds')

### FIGURE 5I ####
blk='bueno'
corr_res = ggscatter (
            blk_meta[[blk]][,c('TOX2','GC_B_cells','ConsensusCluster')], 
            x = 'TOX2', 
            y = 'GC_B_cells',
            #palette = bulk_palette 
            shape=16,
            color = 'ConsensusCluster',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste('TOX2', collapse=' '), ylab = paste('GC_B_cells', collapse=' '),
            title = paste0(blk,'(n=',nrow(blk_meta[[blk]]),')'), fullrange = TRUE) + 
            scale_color_manual (values=palette_bulk) + 
            NoLegend() + gtheme
            stat_cor(aes(color = ConsensusCluster,
               label =paste(..rr.label.., ..p.label.., sep = "~`,`~")), # HOW TO CHANGE p.label to show stars???
           label.x = 1)  
  
corr_res = ggMarginal(corr_res, type = "density", groupColour = TRUE, groupFill = TRUE)

pdf (paste0('Plots/FIGURE_5I_TLS_correlation_bueno.pdf'))
print (corr_res)
dev.off()


### FIGURE IE - cell type abundances deconvolved from bulk using scRNA with BayesPrism ####
blk='bueno'
celltype_names = c('B_cells','T_cells','Myeloid','NK','pDC','Malignant','Endothelial','Fibroblasts','Mesothelium','Plasma','SmoothMuscle','Alveolar','Glia')
bueno_bp = gather (blk_meta[[blk]][,c(celltype_names, 'histology')], celltype, fraction, 1:length(celltype_names))
bueno_bp$celltype = factor (bueno_bp$celltype, levels = celltype_names)
bk_fractions = ggplot(bueno_bp, aes (x= histology, y= fraction)) + 
        ggtitle (paste('RNAseq bueno cohort')) + 
        geom_boxplot (aes (fill = histology), outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        facet_wrap (~celltype, drop=TRUE, scales = 'free', ncol=7) +
        gtheme_no_text +
        NoLegend () +
        scale_fill_manual (values = palette_bulk)

### Run Dirichlet regression and extract pvalues ####
# Bueno data
cc_box1 = blk_meta[[blk]][,c(celltype_names,'histology')]
cc_box1 = cc_box1[cc_box1$histology %in% c('Sarcomatoid','Epithelioid'),]
cc_box1$histology = factor (cc_box1$histology, levels = c('Sarcomatoid','Epithelioid'))
AL = DR_data (cc_box1[,1:length (celltype_names)])
res = DirichReg (AL ~ histology, cc_box1)

# Code from Smillie UC paper: https://github.com/cssmillie/ulcerative_colitis/blob/40a6846565401804d5e2b08e82b52c06d12f0518/analysis.r#L247
u = summary (res)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('histology', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
pvals = as.data.frame (pvals)
pvals$histology = rownames (pvals)
pvals = gather (pvals, celltype, p, 1:(ncol(pvals)-1))
pvals$histology = gsub ('as.factor\\(\\)','',pvals$histology)
pvals$group1 = 'Sarcomatoid'
pvals$group2 = pvals$histology
pvals$p.adj = p.adjust (pvals$p, method ='fdr')
pvals$p_sig = ''
pvals$p_sig = ifelse (pvals$p.adj <= 0.05, 
  ifelse (pvals$p.adj <= 0.01,
    ifelse (pvals$p.adj <= 0.001,'***','**'),'*'),'ns')
# Get y_positions for pval annotation
theta_sp1 = split (bueno_bp, bueno_bp$celltype)
y_max = do.call (rbind, lapply(theta_sp1, function(x) 
  {
  tpm = boxplot(x$fraction ~ x$histology)$stats
  tpm = max (tpm[5,])
  tmp = data.frame (celltype = x$celltype[1], y.position = tpm)
  }))
pvals$y.position = y_max$y.position

pdf (paste0('Plots/S1H_bayesPrism_celltype_bueno_fractions.pdf'),width=10, height=4)
print (bk_fractions + stat_pvalue_manual(
    pvals, step.increase = 0.001,
    label = "p_sig", hide.ns=T
    ))
dev.off()

main_celltypes = c('Malignant','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','B_cells','T_cells','NK')
bueno_bp = gather (blk_meta[[blk]][,c(main_celltypes, 'histology')], celltype, fraction, 1:length(main_celltypes))
bueno_bp = bueno_bp[bueno_bp$celltype %in% main_celltypes,]
bueno_bp$celltype = factor (bueno_bp$celltype, levels = main_celltypes)
bk_fractions = ggplot(bueno_bp, aes (x= histology, y= fraction)) + 
        ggtitle (paste('RNAseq bueno cohort')) + 
        geom_boxplot (aes (fill = histology), outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        facet_wrap (~celltype, drop=FALSE, scales = 'free', ncol=4) +
        gtheme_no_text +
        NoLegend () +
        scale_fill_manual (values = palette_bulk)
pvals = pvals[pvals$celltype %in% main_celltypes,]
pvals = pvals[match (main_celltypes, pvals$celltype),]
pvals$celltype = factor (pvals$celltype, levels = main_celltypes)
pdf (paste0('Plots/FIGURE_1E_bayesPrism_celltype_bueno_fractions.pdf'),width=7, height=5)
print (bk_fractions + stat_pvalue_manual(
    pvals, step.increase = 0.00,
    label = "p_sig", hide.ns=T
    ))
dev.off()

main_celltypes = c('Mesothelium','Glia','Alveolar','pDC','Plasma')
bueno_bp = gather (blk_meta[[blk]][,c(main_celltypes, 'histology')], celltype, fraction, 1:length(main_celltypes))
bueno_bp = bueno_bp[bueno_bp$celltype %in% main_celltypes,]
bueno_bp$celltype = factor (bueno_bp$celltype, levels = main_celltypes)
bk_fractions = ggplot(bueno_bp, aes (x= histology, y= fraction)) + 
        ggtitle (paste('RNAseq bueno cohort')) + 
        geom_boxplot (aes (fill = histology), outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        facet_wrap (~celltype, drop=FALSE, scales = 'free', ncol=4) +
        gtheme_no_text +
        NoLegend () +
        scale_fill_manual (values = palette_bulk)
pvals_sub = pvals[pvals$celltype %in% main_celltypes,]
pvals_sub = pvals_sub[match (main_celltypes, pvals_sub$celltype),]
pvals_sub$celltype = factor (pvals_sub$celltype, levels = main_celltypes)
pdf (paste0('Plots/FIGURE_S1H_bayesPrism_celltype_bueno_fractions.pdf'),width=7, height=5)
print (bk_fractions + stat_pvalue_manual(
    pvals_sub, step.increase = 0.00,
    label = "p_sig", hide.ns=T
    ))
dev.off()

# Hmeljak data
blk='tcga'
main_celltypes = c('Malignant','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','B_cells','T_cells','NK','Mesothelium','Glia','Alveolar','pDC','Plasma')
tcga_bp = gather (blk_meta[[blk]][,c(celltype_names, 'histology')], celltype, fraction, 1:length(celltype_names))
tcga_bp$celltype = factor (tcga_bp$celltype, levels = celltype_names)
bk_fractions = ggplot(tcga_bp, aes (x= histology, y= fraction)) + 
        ggtitle (paste('RNAseq Hmeljak cohort')) + 
        geom_boxplot (aes (fill = histology), outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        facet_wrap (~celltype, drop=TRUE, scales = 'free', ncol=7) +
        NoLegend () +
        gtheme_no_text +
        scale_fill_manual (values = palette_bulk)

cc_box1 = blk_meta[[blk]][,c(celltype_names,'histology')]
cc_box1 = cc_box1[cc_box1$histology %in% c('Sarcomatoid','Epithelioid'),]
cc_box1$histology = factor (cc_box1$histology, levels = c('Sarcomatoid','Epithelioid'))
AL = DR_data(cc_box1[,1:length (celltype_names)])
res = DirichReg (AL ~ histology, cc_box1)

# Code from Smillie UC paper: https://github.com/cssmillie/ulcerative_colitis/blob/40a6846565401804d5e2b08e82b52c06d12f0518/analysis.r#L247
u = summary (res)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('histology', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
pvals = as.data.frame (pvals)
pvals$histology = rownames (pvals)
pvals = gather (pvals, celltype, p, 1:(ncol(pvals)-1))
pvals$histology = gsub ('as.factor\\(\\)','',pvals$histology)
pvals$group1 = 'Sarcomatoid'
pvals$group2 = pvals$histology
pvals$p.adj = p.adjust (pvals$p, method ='fdr')
pvals$p_sig = ''
pvals$p_sig = ifelse (pvals$p.adj <= 0.05, 
  ifelse (pvals$p.adj <= 0.01,
    ifelse (pvals$p.adj <= 0.001,'***','**'),'*'),'ns')
# Get y_positions for pval annotation
theta_sp1 = split (tcga_bp, tcga_bp$celltype)
y_max = do.call (rbind, lapply(theta_sp1, function(x) 
  {
  tpm = boxplot(x$fraction ~ x$histology)$stats
  tpm = max (tpm[5,])
  tmp = data.frame (celltype = x$celltype[1], y.position = tpm)
  }))
pvals$y.position = y_max$y.position  

pdf (paste0('Plots/FIGURE_S1I_bayesPrism_celltype_hmeljak_fractions.pdf'),width=10, height=4)
print (bk_fractions + stat_pvalue_manual(
    pvals, step.increase = 0.001,
    label = "p_sig", hide.ns=T
    ))
dev.off()


# FIGURE 2D - Validate Cm modules ####

# For malignant modules ####
blk = c('bueno','tcga')
modules = paste0('Cm', 1:20)

bp_l = lapply (seq_along(blk_meta), function(y) 
  {
  blk_long = gather (blk_meta[[y]][,c(modules, 'histology')],module, avg_expression, 1:(ncol(blk_meta[[y]][,c(modules, 'histology')])-1)) 
  tmp = blk_long %>% group_by (module) %>%
    t_test(reformulate ('histology', 'avg_expression'), comparisons = list (c('Sarcomatoid','Epithelioid'))) %>%
    adjust_pvalue (method = "fdr") %>%
    add_significance ()
  tmp$study = blk[y]
  as.data.frame(tmp)
  })
    
stat_test_df = do.call (rbind, bp_l)
    
stat_test_df = stat_test_df[stat_test_df$group1 == 'Sarcomatoid' & stat_test_df$group2 == 'Epithelioid',]
stat_test_df$mod = stat_test_df[,1]
saveRDS (stat_test_df, '../PM_scRNA_atlas/data/bueno_tcga_Cms_stats.rds')  


bp_l = lapply (seq_along(blk_meta), function(y) 
  {
    blk_long = gather (blk_meta[[y]][,c(modules, 'histology')],module, avg_expression, 1:(ncol(blk_meta[[y]][,c(modules, 'histology')])-1)) 
    tmp = ggplot(blk_long, aes (x= histology, y= avg_expression)) + 
    ggtitle (paste('RNAseq',names(blk_meta)[y],'cohort')) +
    geom_boxplot (aes (fill = histology), 
          outlier.colour="black", outlier.shape=16,
          outlier.size=2, 
          outlier.alpha = 0.2, 
          notch=FALSE,
          alpha = 0.7, 
          lwd=.2) +
    NoLegend () +
    scale_fill_manual (values = palette_bulk) +
    gtheme_no_text +
    facet_wrap (~module)
  
    stat.test = tmp$data %>% group_by (module) %>%
      t_test(reformulate ('histology', 'avg_expression'), comparisons = list (c('Sarcomatoid','Epithelioid'))) %>%
      adjust_pvalue (method = "fdr") %>%
      add_significance ()
    stat.test = stat.test %>% add_xy_position (x = 'histology', step.increase=0.5)
    
    tmp + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = T,
      lwd=.2,
      label = "p.adj.signif") + 
      NoLegend()        
    })
    
pdf ('Plots/Cm_bulk_data.pdf', 8,8)
print (wrap_plots (bp_l[[1]]))
print (wrap_plots (bp_l[[2]]))
dev.off()

### FIGURE 4G ####
blk='bueno'
corr_res = lapply (c('CXCL9','CXCL10','CXCL11'), function(x) ggscatter (
            blk_meta[[blk]][,c(x, 'Tcell_ab','ConsensusCluster')], 
            x = x, 
            y = 'Tcell_ab',
            #palette = bulk_palette 
            shape=16,
            color = 'ConsensusCluster',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(x, collapse=' '), ylab = paste('Tcell_ab', collapse=' '),
            title = paste0(blk,'(n=',nrow(blk_meta[[blk]]),')'), fullrange = TRUE) + 
            #scale_fill_manual (values=bulk_palette) + 
            scale_color_manual (values=palette_bulk) +
            gtheme)
#            facet_wrap (~ConsensusCluster, drop=TRUE, scales = 'free_x', ncol=length (unique(bueno_meta2$ConsensusCluster))) +
            #NoLegend() +
            #stat_cor(aes(color = ConsensusCluster,
            #   label =paste(..rr.label.., ..p.label.., sep = "~`,`~")), # HOW TO CHANGE p.label to show stars???
           #label.x = -4)  
corr_res = lapply (corr_res, function(x) ggMarginal(x, type = "density", groupColour = TRUE, groupFill = TRUE))
  

pdf (paste0('Plots/FIGURE_4G_CXCLs_correlation_bueno.pdf'), width = 11, height = 5)
print (wrap_plots (corr_res))
dev.off()


# FIGURE 5E - Validate Tm modules ####
modules = c('Tm2','Tm5','Tm7')
bp_l = lapply (seq_along(blk_meta), function(y) 
  {
    blk_long = gather (blk_meta[[y]][,c(modules, 'histology')],module, avg_expression, 1:(ncol(blk_meta[[y]][,c(modules, 'histology')])-1)) 
    tmp = ggplot(blk_long, aes (x= histology, y= avg_expression)) + 
    ggtitle (paste('RNAseq',names(blk_meta)[y],'cohort')) +
    geom_boxplot (aes (fill = histology), 
          outlier.colour="black", outlier.shape=16,
          outlier.size=2, 
          outlier.alpha = 0.2, 
          notch=FALSE,
          alpha = 0.7, 
          lwd=.2) +
    NoLegend () +
    scale_fill_manual (values = palette_bulk) +
    gtheme_no_text +
    facet_wrap (~module)
  
    stat.test = tmp$data %>% group_by (module) %>%
      t_test(reformulate ('histology', 'avg_expression'), comparisons = list (c('Sarcomatoid','Epithelioid'))) %>%
      adjust_pvalue (method = "fdr") %>%
      add_significance ()
    stat.test = stat.test %>% add_xy_position (x = 'histology', step.increase=0.5)
    
    tmp + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = T,
      lwd=.2,
      label = "p.adj.signif") + 
      NoLegend()        
    })
    
#stat_testL2[[mod_name]] = stat_testL
pdf (paste0 ('Plots/FIGURE_5E_bueno_tcga_T_modules_boxplots.pdf'), height=2, width=6)
print (wrap_plots (bp_l[[1]],bp_l[[2]]))
dev.off ()
    

### FIGURE 3H - PLVAP expression in bulk ####
# Run genes on bulk datasets
modules = c('PLVAP_corrected')
bp_l = lapply (seq_along(blk_meta), function(y) 
  {
   blk_long = gather (blk_meta[[y]][,c(modules, 'histology')],module, avg_expression, 1:(ncol(blk_meta[[y]][,c(modules, 'histology')])-1)) 
    tmp = ggplot(blk_long, aes (x= histology, y= avg_expression)) + 
    ggtitle (paste('RNAseq',names(blk_meta)[y],'cohort')) +
    geom_boxplot (aes (fill = histology), 
          outlier.colour="black", outlier.shape=16,
          outlier.size=2, 
          outlier.alpha = 0.2, 
          notch=FALSE,
          alpha = 0.7, 
          lwd=.2) +
    NoLegend () +
    scale_fill_manual (values = palette_bulk) +
    gtheme_no_text +
    facet_wrap (~module)
  
    stat.test = tmp$data %>% group_by (module) %>%
      t_test(reformulate ('histology', 'avg_expression'), comparisons = list (c('Sarcomatoid','Epithelioid'))) %>%
      adjust_pvalue (method = "fdr") %>%
      add_significance ()
    stat.test = stat.test %>% add_xy_position (x = 'histology', step.increase=0.5)
    
    tmp + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = T,
      lwd=.2,
      label = "p.adj.signif") + 
      NoLegend()        
    })
    
    #stat_testL[[blk]] = data.frame (stat.test, dataset = blk)
    
#stat_testL2[[mod_name]] = stat_testL
pdf (paste0 ('Plots/FIGURE_3H_bueno_tcga_PLVAP_boxplots.pdf'),height=3,width=4)
print (wrap_plots (bp_l[[1]],bp_l[[2]]))
dev.off ()
    
    

### FIGURE 2G S2J - Correlate chr22 region expression module vs Cm2 / Cm17 ####
corr_res = list()
blk = 'bueno'
mod = c('Cm2','Cm17')

for (i in mod)
  {
  corr_res = ggscatter (
            blk_meta[[blk]], 
            x = i,
            y = 'chr22',
            #palette = bulk_palette 
            shape=16,
            color = 'histology',
            add = "reg.line", conf.int = FALSE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(i, collapse=' '), ylab = paste('chr22', collapse=' ')) +
  #          title = paste0(blk,'(n=',nrow(exp_df),')'), fullrange = TRUE) + 
            scale_fill_manual (values=palette_bulk) + 
            scale_color_manual (values=palette_bulk) +
            NoLegend()  + gtheme_no_rot# + 
            #geom_point (data = exp_df, aes (x = your.gene1, y = your.gene2, color = histology, fill= histology))
            #geom_smooth(method = "lm", color = "black") +
            #stat_cor (label.x = 3)
  corr_res = ggMarginal(corr_res, type = "density", groupColour = TRUE, groupFill = TRUE)            
  pdf (paste0 ('Plots/FIGURE_2G_S2J_chr22_',i,'_correlation_scatterplots2.pdf'),4,4)
  print (corr_res)
  dev.off()
  }


### Survival Analysis ####

# FIGURE S2G-H - Run Cox hazard ratio regression survival analysis on nmf modules ####
# Set variables per dataset to use
#hist_column = c(bueno = 'HistologyReduced', mesomics = 'Type', tcga = 'TUMOR_TYPE')
#hist_column = c(bueno = 'HistologyReduced', tcga = 'TUMOR_TYPE')
hist_column = c(bueno = 'ConsensusCluster', tcga = 'TUMOR_TYPE')
event_column = c(bueno = 'Survival.from.surgery..years.', tcga = 'OS_MONTHS')
death_column = c(bueno = 'DEATH_EVENT', tcga = 'DEATH_EVENT')
low='1st Qu.'
high='3rd Qu.'  

modules = c(paste0('Cm',1:20))
cfit_study = list()
cox_l = list()
sc_p = list()
plot_study= list()
for (study in names(blk_meta))
    {
    cfit = list()
    for (mod in modules)
        {
        meta_surv = blk_meta[[study]]
        meta_surv = meta_surv[!is.na(meta_surv[,event_column[[study]]]),]
        meta_surv$DEATH_EVENT = meta_surv[,death_column[[study]]]    
        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        CI <- round(exp(confint(cfit[[mod]])), 2)
        cox_df = data.frame (
          HR = round(exp(coef(cfit[[mod]])), 2),
          CI = paste0('(',paste (CI, collapse=' - '),')'),
          LL = CI[1],
          UL = CI[2],
          P_value_C = round(summary(cfit[[mod]])$coefficients[, 5],2),
          label = mod
          )
        
        raw.vec=meta_surv[,mod]
        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        meta_surv[,mod] = factor (classified.vec, levels = c('Low','Med','High'))

        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        s = summary (cfit[[mod]])
        cox_df$P_value_S = round(s$logtest[3],2)  #s$logtest[3]
        cox_l[[mod]] = cox_df
        
        sc_p[[mod]] = ggadjustedcurves (cfit[[mod]], 
                data = meta_surv, 
                method = "conditional",
                variable=mod,
                palette = 'aaas',
                size=0.4,
                surv.median.line = 'hv',
                ggtheme = theme_classic()) +
                labs (title = mod,
                subtitle = paste0('log-rank = ',round(s$logtest[3],2)),
                caption = paste("n = ", nrow(meta_surv))) +
                      ylim(0,1)+ geom_hline(yintercept = 0.5,c(0.5,0.5),linetype='dotted', col='grey22') + gtheme_text
        }
    cox_df = do.call (rbind, cox_l)
    cox_df$Index = rownames (cox_df)
    cox_df$Index = factor (cox_df$Index, levels = rev(cox_df$Index))
    cfit_study[[study]] = cox_df
    plot_study[[study]] = sc_p
    }

for (study in names (cfit_study))
  {
  forest <- ggplot(cfit_study[[study]], aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()
      
  tab <- ggplot(cfit_study[[study]], aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    geom_text(aes(x = 3, label = P_value_S), size = 4) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (paste0('Plots/cox_regression_',study,'.pdf'), 5,5)
    print (forest + tab + plot_layout(ncol = 2, widths = c(1, 3)))
    dev.off()
    pdf (paste0('Plots/cox_regression_',study,'_stratified.pdf'), height = 2.8,3)
    print (plot_study[[study]])
    dev.off()
  }

### FIGURE S3I - Run Cox regression model for fetal endothelial (PLVAP) ####
hist_column = c(bueno = 'ConsensusCluster', tcga = 'TUMOR_TYPE')
event_column = c(bueno = 'Survival.from.surgery..years.', tcga = 'OS_MONTHS')
death_column = c(bueno = 'DEATH_EVENT', tcga = 'DEATH_EVENT')
low='1st Qu.'
high='3rd Qu.'  

modules = c('fetal_endo')
cfit_study = list()
cox_l = list()
sc_p = list()
plot_study= list()
for (study in names(blk_meta))
    {
    cfit = list()
    for (mod in modules)
        {
        meta_surv = blk_meta[[study]]
        meta_surv = meta_surv[!is.na(meta_surv[,event_column[[study]]]),]
        meta_surv$DEATH_EVENT = meta_surv[,death_column[[study]]]    
        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        CI <- round(exp(confint(cfit[[mod]])), 2)
        cox_df = data.frame (
          HR = round(exp(coef(cfit[[mod]])), 2),
          CI = paste0('(',paste (CI, collapse=' - '),')'),
          LL = CI[1],
          UL = CI[2],
          P_value_C = round(summary(cfit[[mod]])$coefficients[, 5],2),
          label = mod
          )
        
        raw.vec=meta_surv[,mod]
        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        meta_surv[,mod] = factor (classified.vec, levels = c('Low','Med','High'))

        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        s = summary (cfit[[mod]])
        cox_df$P_value_S = round(s$logtest[3],2)  #s$logtest[3]
        cox_l[[mod]] = cox_df
        
        sc_p[[mod]] = ggadjustedcurves (cfit[[mod]], 
                data = meta_surv, 
                method = "conditional",
                variable=mod,
                palette = 'aaas',
                surv.median.line = 'hv',
                ggtheme = theme_classic()) +
                labs (title = mod,
                subtitle = paste0('log-rank = ',round(s$logtest[3],2)),
                caption = paste("n = ", nrow(meta_surv))) +
                      ylim(0,1)+ geom_hline(yintercept = 0.5,c(0.5,0.5),col='grey')
        }
    cox_df = do.call (rbind, cox_l)
    cox_df$Index = rownames (cox_df)
    cox_df$Index = factor (cox_df$Index, levels = rev(cox_df$Index))
    cfit_study[[study]] = cox_df
    plot_study[[study]] = sc_p
    }

cfit_study = do.call (rbind, cfit_study)
cfit_study$Index = c('Bueno','Hmeljak')
  forest <- ggplot(cfit_study, aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()
      
  tab <- ggplot(cfit_study, aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    geom_text(aes(x = 3, label = P_value_S), size = 4) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (paste0('Plots/FIGURE_S3I_cox_regression_PLVAP.pdf'), 5,1.6)
    print (forest + tab + plot_layout(ncol = 2, widths = c(1, 3)))
    dev.off()
    pdf (paste0('Plots/FIGURE_S3I_regression_PLVAP_stratified.pdf'), 5,5)    
    print (plot_study[[study]])
    dev.off()
  

# Run Kaplan Mayer survival curves ####
studies = c('bueno','tcga')
hist_column = c(bueno = 'ConsensusCluster',  tcga = 'TUMOR_TYPE')
subset_hist = c(bueno = 'Sarcomatoid',tcga = 'Epithelioid Mesothelioma')
subset_hist = c(bueno = 'Epithelioid',tcga = 'Epithelioid Mesothelioma')
event_column = c(bueno = 'Survival.from.surgery..years.', tcga = 'OS_MONTHS')
low='1st Qu.'
high='3rd Qu.'  

modules = c('Cm5','Cm12','NK')
km_p_study = list()
for (study in names(blk_meta))
    {
    km_p = list()
    for (mod in modules)
        {
        if (!is.null(subset_hist[[study]])) meta_surv = blk_meta[[study]][blk_meta[[study]][,hist_column[[study]]] %in% subset_hist[[study]],] else
        meta_surv = blk_meta[[study]]
        raw.vec=meta_surv[,mod]

        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        
        
        meta_surv$classified.vec = factor (classified.vec, levels = c('Low','Med','High'))
        meta_surv = meta_surv[!is.na(meta_surv[,event_column[[study]]]),]
        
        meta_surv$binary_histo = meta_surv[,hist_column[[study]]]

        meta_surv$DEATH_EVENT = meta_surv[,death_column[[study]]]            
        meta_surv$OS_MONTHS = meta_surv[,event_column[[study]]]
        
        model_fit = survfit(Surv(OS_MONTHS, DEATH_EVENT) ~ classified.vec, data = meta_surv)
        model_diff = survdiff(Surv(OS_MONTHS, DEATH_EVENT) ~ classified.vec, data = meta_surv)
        
        km_p[[mod]] = ggsurvplot(
        model_fit,
        data = meta_surv,
        size = 1,                 # change line size
        palette ="lancet",        # custom color palettes
        conf.int = FALSE,          # Add confidence interval
        pval = FALSE,              # Add p-value
        risk.table = FALSE,        # Add risk table
        risk.table.col = "strata",# Risk table color by groups
        legend.labs =levels (meta_surv$classified.vec),    # Change legend labels
        risk.table.height = 0.25, # Useful to change when you have multiple groups
        ylab = 'Survival rate',
        #ggtheme = theme_survminer()
        ggtheme = theme_classic(),#,
        title = mod,
        subtitle = paste('pval',round(model_diff$pvalue,2)),
        caption = paste("n = ", nrow(meta_surv))
        )
        }
    km_p_study[[study]] = km_p
    }    

for (study in studies)
    {
    pdf (paste0('Plots/FIGURE_S2I_6A_Kaplan_Mayer_survival_curves_',study,'subset_',subset_hist[[study]],'.pdf'), 3,4)
    print (km_p_study[[study]])
    dev.off()
    }

# Re-run only in Sarcomatoid
subset_hist = c(bueno = 'Sarcomatoid')
low='1st Qu.'
high='3rd Qu.'  
modules = c('Cm5','Cm12','NK')
km_p_study = list()
for (study in names(blk_meta)[1])
    {
    km_p = list()
    for (mod in modules)
        {
        if (!is.null(subset_hist[[study]])) meta_surv = blk_meta[[study]][blk_meta[[study]][,hist_column[[study]]] %in% subset_hist[[study]],] else
        meta_surv = blk_meta[[study]]
        raw.vec=meta_surv[,mod]

        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        
        
        meta_surv$classified.vec = factor (classified.vec, levels = c('Low','Med','High'))
        meta_surv = meta_surv[!is.na(meta_surv[,event_column[[study]]]),]
        
        meta_surv$binary_histo = meta_surv[,hist_column[[study]]]

        meta_surv$DEATH_EVENT = meta_surv[,death_column[[study]]]            
        meta_surv$OS_MONTHS = meta_surv[,event_column[[study]]]
        
        model_fit = survfit(Surv(OS_MONTHS, DEATH_EVENT) ~ classified.vec, data = meta_surv)
        model_diff = survdiff(Surv(OS_MONTHS, DEATH_EVENT) ~ classified.vec, data = meta_surv)
        
        km_p[[mod]] = ggsurvplot(
        model_fit,
        data = meta_surv,
        size = 1,                 # change line size
        palette ="lancet",        # custom color palettes
        conf.int = FALSE,          # Add confidence interval
        pval = FALSE,              # Add p-value
        risk.table = FALSE,        # Add risk table
        risk.table.col = "strata",# Risk table color by groups
        legend.labs =levels (meta_surv$classified.vec),    # Change legend labels
        risk.table.height = 0.25, # Useful to change when you have multiple groups
        ylab = 'Survival rate',
        #ggtheme = theme_survminer()
        ggtheme = theme_classic(),#,
        title = mod,
        subtitle = paste('pval',round(model_diff$pvalue,2)),
        caption = paste("n = ", nrow(meta_surv))
        )
        }
    km_p_study[[study]] = km_p
    }    

for (study in studies)
    {
    pdf (paste0('Plots/FIGURE_S2I_6A_Kaplan_Mayer_survival_curves_',study,'subset_',subset_hist[[study]],'.pdf'), 3,4)
    print (km_p_study[[study]])
    dev.off()
    }

### FIGURE S6B - Run Cox regression model on BayesPrism predictions ####
celltype_names = c('B_cells','T_cells','Myeloid','NK','pDC','Malignant','Endothelial','Fibroblasts','Mesothelium','Plasma','SmoothMuscle','Alveolar','Glia')
hist_column = c(bueno = 'ConsensusCluster', tcga = 'TUMOR_TYPE')
event_column = c(bueno = 'Survival.from.surgery..years.', tcga = 'OS_MONTHS')
death_column = c(bueno = 'DEATH_EVENT', tcga = 'DEATH_EVENT')
low='1st Qu.'
high='3rd Qu.'  

modules = celltype_names
cfit_study = list()
cox_l = list()
sc_p = list()
plot_study= list()
for (study in names(blk_meta))
    {
    cfit = list()
    for (mod in modules)
        {
        meta_surv = blk_meta[[study]]
        meta_surv[[mod]] = meta_surv[[mod]]+1
        meta_surv = meta_surv[!is.na(meta_surv[,event_column[[study]]]),]
        meta_surv$DEATH_EVENT = meta_surv[,death_column[[study]]]    
        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        CI = round(confint(cfit[[mod]]), 2)
        cox_df = data.frame (
          HR = round(coef(cfit[[mod]]), 2),
          CI = paste0('(',paste (CI, collapse=' - '),')'),
          LL = CI[1],
          UL = CI[2],
          P_value_C = round(summary(cfit[[mod]])$coefficients[, 5],2),
          label = mod
          )
        
        raw.vec=meta_surv[,mod]
        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        meta_surv[,mod] = factor (classified.vec, levels = c('Low','Med','High'))

        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,event_column[[study]]])),
                            DEATH_EVENT) ~', mod, '+ strata (histology)'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        s = summary (cfit[[mod]])
        cox_df$P_value_S = round(s$logtest[3],2)  #s$logtest[3]
        cox_l[[mod]] = cox_df
        
        sc_p[[mod]] = ggadjustedcurves (cfit[[mod]], 
                data = meta_surv, 
                method = "conditional",
                variable=mod,
                palette = 'aaas',
                surv.median.line = 'hv',
                ggtheme = theme_classic()) +
                labs (title = mod,
                subtitle = paste0('log-rank = ',round(s$logtest[3],2)),
                caption = paste("n = ", nrow(meta_surv))) +
                      ylim(0,1)+ geom_hline(yintercept = 0.5,c(0.5,0.5),col='grey')
        }
    cox_df = do.call (rbind, cox_l)
    cox_df$Index = rownames (cox_df)
    cox_df$Index = factor (cox_df$Index, levels = rev(cox_df$Index))
    cfit_study[[study]] = cox_df
    plot_study[[study]] = sc_p
    }

for (study in names (cfit_study))
  {
  forest <- ggplot(cfit_study[[study]], aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()
      
  tab <- ggplot(cfit_study[[study]], aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    geom_text(aes(x = 3, label = P_value_S), size = 4) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (paste0('Plots/FIGURE_S6B_cox_regression_',study,'_BayesPrism.pdf'), width=8,5)
    print (forest + tab + plot_layout(ncol = 2, widths = c(1, 3)))
    dev.off()
    pdf (paste0('Plots/FIGURE_S6B_cox_regression_',study,'_stratified_BayesPrism.pdf'), height = 2.8,3)
    print (plot_study[[study]])
    dev.off()
  }


# FIGURE_2G - Check chr22 from bueno clinical var ####
module = 'Cm17_cnv'
blk = 'bueno'
blk_meta[[blk]]$FISH.chrom22 = as.character (blk_meta[[blk]]$FISH.chrom22)
blk_meta[[blk]]$FISH.chrom22[blk_meta[[blk]]$FISH.chrom22 == ''] = 'normal'
blk_meta[[blk]]$FISH.chrom22 = factor (blk_meta[[blk]]$FISH.chrom22, levels = c('normal','del 22'))
bueno_fish = blk_meta[[blk]][blk_meta[[blk]]$FISH.chrom22 %in% levels (blk_meta[[blk]]$FISH.chrom22),]

bp_l = ggplot(bueno_fish, aes_string (x= 'FISH.chrom22', y= module)) + 
        ggtitle (paste('RNAseq cohort')) +
        geom_boxplot (aes (fill = FISH.chrom22),outlier.colour="black", outlier.shape=16,
                 outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
        ylab (module) +
        ggtitle ('FISH.chrom22') +
        gtheme_no_text +
        scale_fill_manual (values=c(rev(paletteer::paletteer_d("nord::mountain_forms")[c(1,3)])))
        
        stat.test = bp_l$data %>%
              t_test(reformulate ('FISH.chrom22', module)) %>%
              adjust_pvalue (method = "none") %>%
              add_significance ()
        stat.test = stat.test %>% add_xy_position (x = 'FISH.chrom22', step.increase=0)
#        stat.test = stat.test[stat.test$group1 %in% c('Sarcomatoid','Epithelioid') & stat.test$group2 %in% c('Sarcomatoid','Epithelioid'),]
        bp_l = bp_l + stat_pvalue_manual (stat.test, 
          remove.bracket=FALSE,
          bracket.nudge.y = -.4,
          size = 2.5, 
          hide.ns = T,
          label = "p.adj.signif")

pdf (paste0('Plots/FIGURE_2G_bueno_chr22_deletion_',module,'_covar_FISH_chrom22_boxplot.pdf'), width = 2.3,height=2)
print (bp_l)
dev.off()

