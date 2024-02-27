# load libraries
message ('Load R packages')
packages = c(
  'Seurat',
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
  'tidyr',
  'gdata',
  'GSA',
  'readxl',
  'scran',
  'harmony',
  'org.Mm.eg.db',
  'org.Hs.eg.db',
  'viridis',
  'BiocParallel',
  'clusterProfiler',
  'igraph',
  'tidyverse',
  'rstatix',
  'paletteer',
  'extrafont',
  'scCustomize',
  'circlize',
  'SeuratDisk'
)
lapply(packages, require, character.only = TRUE)

