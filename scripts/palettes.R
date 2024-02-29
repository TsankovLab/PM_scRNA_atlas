


palette_celltype_simplified = c(
  B_cells = 'magenta2',
  T_cells = 'firebrick1',
  MonoMac = 'cornflowerblue',
  Myeloid = 'cornflowerblue',
  DCs = 'deepskyblue',
  NK = 'gold1',
  Fibroblasts = 'azure4',
  Endothelial = 'brown',
  SmoothMuscle = 'blueviolet',
  Plasma = 'lightsalmon1',
  Alveolar = 'black',
  #Malignant = 'palegreen4',
  Malignant = 'plum4',
  Mast = 'royalblue4',
  Mesothelium = 'olivedrab1',
  Glia = 'lawngreen',
  pDC = 'tomato')



palette_gene_expression = as.character(paletteer::paletteer_c("grDevices::Blue-Red 3", 100))
palette_gene_expression2 = as.character(paletteer::paletteer_c("grDevices::Blue-Red 3", 100))[c(1,50,100)]
palette_gene_expression_fun = function(x) {return(colorRamp2(c(min(x), 0, max(x)), c(palette_gene_expression2)))}
palette_sample = c(
'#7c10a3
#9443b3
#ac6ac2
#c28fd2
#d7b4e1
#ebd9f0
#ffffff
#d6dbe2
#adb9c5
#8697aa
#5f778f
#375974
#003c5b')
palette_sample = c(unlist(strsplit (palette_sample, '\n')))
palette_sample = c(rev(as.character(paletteer::paletteer_c("pals::ocean.curl",12))), 'grey')
palette_sample = setNames (palette_sample , c('P1','P13','P3','P12','P6','P2','P5','P7','P11','P4','P8','P9','P10'))
#palette_sample = setNames (as.character(paletteer::paletteer_c("pals::ocean.dense",13)))
#palette_sample = c(palette_sample, P10 = 'grey')
palette_sample2 = c(palette_sample[3], palette_sample[6], palette_sample[10])
palette_SE_group = c('S-High'= unname(palette_sample[3]), 'E-High' = unname(palette_sample[10]))
palette_sample2_fun = colorRamp2(c(1, 0, -1), c(palette_sample2))
palette_sample_cnv = colorRamp2(c(.5, 0, -.5), c(palette_sample2))
palette_cnv = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",3))
palette_cnv_fun = colorRamp2(c(-.5,0,.5), palette_cnv)

ov_pal = rev(paletteer::paletteer_c("grDevices::Purple-Blue",100))

palette_module_correlation = paletteer::paletteer_c("pals::kovesi.diverging_bwr_40_95_c42",100)
palette_module_correlation_fun = colorRamp2(c(-1,0,1), c(palette_module_correlation[1],'white',palette_module_correlation[100]))

# Set palette
palette_bulk = setNames (as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3,1,7)]), c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic','E_score','S_score'))

palette_t_cells = setNames (as.character (paletteer::paletteer_d("fishualize::Bodianus_rufus",4)), c('CD8','CD4','Tregs','TFH'))
palette_nk_cells = setNames (c("#D3E3CAFF", "#92A587FF", "#2F3525FF"), c('FGFBP2_NK','KLRC1_NK','NKlike_Tcells'))
palette_tnk_cells = c(palette_t_cells, palette_nk_cells)

palette_protein_expression = c(low="darkblue",mid= "white",high= "darkgreen") 
palette_feature_RNA = c('lightgrey',"#5F1415FF")
palette_feature_protein = c("lightgrey", "darkgreen")

palette_b_cells = c(B_cells = 'mediumorchid1', GC_B_cells = 'magenta4', Plasma = 'lightskyblue1')

palette_clonotype = setNames (c(as.character (paletteer::paletteer_d("beyonce::X58"))),c('NonExpanded','Small','Large'))
#palette_clonotype = palette_clonotype[1:4]

pallette_pbmc_celltype = setNames (rev(as.character(paletteer::paletteer_d("khroma::smoothrainbow")[c(1,3,5,7,9,11,13,15)])), c('CD4','CD8','ILC','MAIT','NK','NK Proliferating','NK_CD56bright','Treg'))


palette_endothelial = setNames(as.character(paletteer::paletteer_d("colRoz::shark_bay", 3)), c('Artery','PLVAP+EC','Vein'))
pallette_fetal_vs_adult = setNames (c('#DB3EB1FF', 'grey','purple'), c('fetal_lung','adult_lung','PM'))
palette_scenic = rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(50))

palette_celltypes_normal = setNames(brewer.pal(11,'Paired'),rev(c('AT1','AT2','Basal','Secretory','Ionocytes','Neuroendocrine','Tuft.like','Mesothelium','Fibroblast','SmoothMuscle','Ciliated')))


#for (pal in ls()[grep('palette',ls())]) store_pal (list(pal = get(pal)))

