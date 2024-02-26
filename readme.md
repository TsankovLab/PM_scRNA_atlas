### This repository contains custom code to reproduce main analyses in our manuscript, "Single-cell view of tumor microenvironment gradients in Pleural Mesothelioma"
Scripts are organized per cellular compartment and figures in the paper. To download Seurat object inclusive of all cell annotation and additional metadata, please download it from GEO GSE190597.

## Installation

To install this repository, run the following commands:

```bash

git clone https://github.com/TsankovLab/PM_scRNA_atlas.git

```

These are the R packages required to reproduce the scRNA analysis. Additionally, to install cNMF (https://github.com/dylkot/cNMF ) and SCENIC (https://github.com/aertslab/SCENIC) softwares please refer to the corresponding GitHub pages

System Requirements:

Operating System: Tested on RedHat 7.7 x86_64.
R Version: 4.2.3
R Packages:
Seurat (v4.4.0) (https://satijalab.org/seurat/)
SeuratObject (v5.0.1)
circlize (v0.4.15) 
ComplexHeatmap (v2.14.0) (https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
clusterProfiler (v4.6.0) (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
ggpubr (v0.6.0) (https://rdrr.io/cran/ggpubr/)
biomaRt (v2.54.0)
paletteer (v1.5.0)
rstatix (v0.7.2)
readxl (v1.4.3)
harmony (v0.1.0) (https://portals.broadinstitute.org/harmony/articles/quickstart.html)
hdWGCNA (v0.2.18) (https://smorabit.github.io/hdWGCNA/)
ggplot2 (v3.4.4) (https://ggplot2.tidyverse.org)
tidyverse (v2.0.0)


```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

install.packages("Seurat", version = "4.4.0")
install.packages ("SeuratObject", version = "5.0.1")
install.packages ("circlize", version = "0.4.15")
BiocManager::install("ComplexHeatmap", version = "2.14.0")
BiocManager::install("clusterProfiler", version = "4.6.0")
install.packages("ggpubr", version = "0.6.0")
BiocManager::install("biomaRt", version = "2.54.0")
install.packages("paletteer", version = "1.5.0")
install.packages("rstatix", version = "0.7.2")
install.packages ("readxl", version = "v1.4.3")
install.packages("harmony", version = "0.1.0")
devtools::install_github('smorabit/hdWGCNA', ref='v0.2.18')
install.packages("ggplot2", version = "3.4.4")
install.packages ("tidyverse", version = "2.0.0")

```
