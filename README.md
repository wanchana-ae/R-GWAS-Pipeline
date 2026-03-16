# R GWAS Pipeline

This repository contains an R-based workflow for genome-wide association studies (GWAS) using SNP data. It covers VCF processing, phenotype analysis, PCA, phylogenetic tree construction, and GWAS analysis using GAPIT.

---

## Features

- Read and process VCF files using `vcfR`.
- Visualize phenotype distributions using `ggplot2`.
- Perform PCA on genotype data (`SNPRelate`) and plot results.
- Construct phylogenetic trees from SNP data (`ape`, `ggtree`).
- Conduct GWAS using GAPIT (`FarmCPU`, `GLM`, `MLM`, `CMLM`, `SUPER`, `MLMM`, `FarmCPU`, and `Blink`).
- Extract specific chromosomes for targeted analysis.
- Generate Manhattan and QQ plots for GWAS results.

---

## Installation

Install required R packages:

```R
install.packages(c(
  "vcfR", "ggplot2", "dplyr", "geneHapR", "reshape2", "tidyverse",
  "seqinr", "rsvg", "ggrepel", "qqman", "ape", "remotes"
))
remotes::install_github('coolbutuseless/svgparser')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("LEA", "ggtree", "SNPRelate"))
```

## Load libraries in your script

```R
library(vcfR)
library(ggplot2)
library(dplyr)
library(geneHapR)
library(LEA)
library(reshape2)
library(svgparser)
library(readr)
library(seqinr)
library(ggtree)
library(rsvg)
library(SNPRelate)
library(ggrepel)
library(ape)
```

## Load GAPIT functions

```R
source("http://zzlab.net/GAPIT/gapit_functions.txt")
```

## Usage

### 1.Input files

- `snp_set.recode.vcf`: SNP dataset in VCF format
- `trait.txt`: Phenotype data file with two columns (Sample and phenotype)

### 2.Phenotype histogram

```R
ggplot(pheno.dat, aes(x = phenotype)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  labs(title = "Histogram of Phenotype", x = "Phenotype Value", y = "Frequency") +
  theme_minimal()
```

### 3.Phylogenetic tree construction

- Convert VCF to GDS
- Generate IBS matrix
- Cluster and visualize with ggtree
- Export tree in Newick format

### 4.PCA analysis

```R
pca <- snpgdsPCA(genofile, autosome.only=FALSE)
# Visualize with ggplot2
```

### 5.GWAS using GAPIT

- Convert VCF to HapMap format
- Run `GAPIT` with FarmCPU model

```R
Rscript GAPIT.R
or
Rscript workshop_master.R
```

- Extract top SNPs and save CSV
- Generate Manhattan and QQ plots using qqman

### 6.Chromosome-specific analysis

```R
vcf.chr12 <- vcf[vcf@fix[,"CHROM"] == 12,]
write.vcf(vcf.chr12, "tomato_chr12.vcf.gz")
```

## References

- **GAPIT Manual**: [https://zzlab.net/GAPIT/gapit_help_document.pdf](https://zzlab.net/GAPIT/gapit_help_document.pdf)  
- **GAPIT GitHub**: [https://github.com/jiabowang/GAPIT](https://github.com/jiabowang/GAPIT)  
- **R SNPRelate**: [https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)  
- **vcfR Package**: [https://cran.r-project.org/web/packages/vcfR/index.html](https://cran.r-project.org/web/packages/vcfR/index.html)  
- **ggtree Package**: [https://bioconductor.org/packages/release/bioc/html/ggtree.html](https://bioconductor.org/packages/release/bioc/html/ggtree.html)  
- **qqman Package**: [https://cran.r-project.org/web/packages/qqman/index.html](https://cran.r-project.org/web/packages/qqman/index.html)  
- **LEA Package**: [https://bioconductor.org/packages/release/bioc/html/LEA.html](https://bioconductor.org/packages/release/bioc/html/LEA.html)  
- **geneHapR Package**: [https://cran.r-project.org/web/packages/geneHapR/index.html](https://cran.r-project.org/web/packages/geneHapR/index.html)
