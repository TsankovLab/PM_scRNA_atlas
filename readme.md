### This repository contains custom code to reproduce main analyses in our manuscript, "Single-cell view of tumor microenvironment gradients in Pleural Mesothelioma"
Scripts are organized per cellular compartment and figures in the paper. To download Seurat object inclusive of all cell annotation and additional metadata, please download it from GEO GSE190597.

## Installation

To install this repository, run the following command:

```bash

git clone https://github.com/TsankovLab/PM_scRNA_atlas.git

```
and download the seurat objects srt_tumor.rds and srt_pbmc.rds available at GSE190597 in same folder

These are the R packages required to reproduce the scRNA and bulkRNA analyses. 
System Requirements:

Operating System: Tested on RedHat 7.7 x86_64.  
R Version: 4.2.3. 
R Packages:  
Seurat (v4.4.0) (https://satijalab.org/seurat/).   
SeuratObject (v5.0.1).  
scran (v1.26.0) (https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html)
circlize (v0.4.15) (https://jokergoo.github.io/circlize_book/book/).  
ComplexHeatmap (v2.14.0) (https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html).  
clusterProfiler (v4.6.0) (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).  
ggpubr (v0.6.0) (https://rdrr.io/cran/ggpubr/).  
biomaRt (v2.54.0) (https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/biomaRt/inst/doc/biomaRt.html).  
paletteer (v1.5.0) (https://emilhvitfeldt.github.io/paletteer/).  
rstatix (v0.7.2) (https://rpkgs.datanovia.com/rstatix/).  
readxl (v1.4.3) (https://readxl.tidyverse.org/).  
harmony (v0.1.0) (https://portals.broadinstitute.org/harmony/articles/quickstart.html).   
hdWGCNA (v0.2.18) (https://smorabit.github.io/hdWGCNA/).  
ggplot2 (v3.4.4) (https://ggplot2.tidyverse.org).  
tidyverse (v2.0.0) (https://www.tidyverse.org/)    
patchwork (v1.1.2) (https://patchwork.data-imaginist.com/).  
scRepertoire (v2.0.0) (https://github.com/ncborcherding/scRepertoire).  
DirichletReg (v0.7.1) (http://cran.nexr.com/web/packages/DirichletReg/vignettes/DirichletReg-vig.pdf)

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Analysis packages
install.packages("Seurat", version = "4.4.0")
remotes::install_github("mojaveazure/seurat-disk", version = "0.0.0.9020")
BiocManager::install("scran", version = "1.26.0")
install.packages("harmony", version = "0.1.0")
BiocManager::install("scRepertoire", version = "2.0.0")
BiocManager::install("clusterProfiler", version = "4.6.0")
BiocManager::install("biomaRt", version = "2.54.0")
install.pakcages ("DirichletReg", version="0.7.1")
install.packages("survival", version = "3.5.5")
install.packages("survminer", version = "0.4.9")
devtools::install_github('smorabit/hdWGCNA', ref='v0.2.18')

# Data manipulation  packages
install.packages("dplyr", version = "1.1.2")
install.packages ("tidyverse", version = "2.0.0")
install.packages ("readxl", version = "v1.4.3")
BiocManager::install("GenomicFeatures", version = "1.50.4")
BiocManager::install("org.Hs.eg.db", version = "3.16.0")

# Visualization packages
install.packages ("circlize", version = "0.4.15")
BiocManager::install("ComplexHeatmap", version = "2.14.0")
BiocManager::install("ggtree", version = "3.6.0")
install.packages("ggrepel", version = "0.9.3")
install.packages("ggpubr", version = "0.6.0")
install.packages("rstatix", version = "0.7.2")

# Palettes
install.packages ("RColorBrewer")
install.packages ("viridis")
install.packages("paletteer", version = "1.5.0")
install.packages("ggplot2", version = "3.4.4")
install.packages ("gplots", version = "3.1.3")
install.packages("ggridges", version = "0.5.4")
install.packages ("ggExtra", version = "0.10.1")
install.packages ("patchwork", version = "1.1.2")

```


## Fetal endothelial analysis
For comparison analysis of PLVAP+ EC with fetal and adult normal lungs please download data from:
https://fetal-lung.cellgeni.sanger.ac.uk/scRNA.html
from study: He, Peng, et al. "A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates." Cell 185.25 (2022): 4841-4860.

and
https://www.synapse.org/#!Synapse:syn21041850/wiki/600865
from: Travaglini, Kyle J., et al. "A molecular cell atlas of the human lung from single-cell RNA sequencing." Nature 587.7835 (2020): 619-625.

## NK cells pan-cancer atlas
For comparison analysis of NK PM cells with NKs from an immune pan-cancer study (FIGURE 6C) please download immune pan-cancer altas from here:
https://zenodo.org/records/5186413#%20.YRqbJC1h2v6 from following study:
Nieto, Paula, et al. "A single-cell tumor immune atlas for precision oncology." Genome research 31.10 (2021): 1913-1926.



