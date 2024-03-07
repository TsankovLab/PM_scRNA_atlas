# load libraries
message ('Load R packages')
packages = c(
  'Seurat',
  'nichenetr',
  'ggplot2',
  'gplots',
  'ggtree',
  'dplyr',
  'Matrix',
  'ggpubr',
  'biomaRt',
  'patchwork',
  'ComplexHeatmap',
  'RColorBrewer',
  'ggrepel',
  'fgsea',
  'readxl',
  'scran',
  'harmony',
  'org.Hs.eg.db',
  'viridis',
  'clusterProfiler',
  'tidyverse',
  'rstatix',
  'paletteer',
  'scCustomize',
  'circlize',
  'SeuratDisk',
  'ggridges',
  'scRepertoire',
  'survival',
  'survminer',
  'ggExtra',
  'DirichletReg',
  'GenomicFeatures'
)
lapply(packages, require, character.only = TRUE)

