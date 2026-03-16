#install.packages('vcfR')
#install.packages('ggplot2')
#install.packages('dplyr')
#install.packages('geneHapR')
#install.packages('reshape2')
#install.packages("tidyverse")
#install.packages("seqinr")
#install.packages("rsvg")
#install.packages("ggrepel")
#install.packages("qqman")
#install.packages("ape")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
#BiocManager::install("ggtree")
#BiocManager::install("SNPRelate")


library(vcfR)
library(ggplot2)
library(dplyr)
library(LEA)
library(reshape2)
library(readr)
library(seqinr)
library(ggtree)
library(rsvg)
library(SNPRelate)
library(ggrepel)
library(ape)

#GAPIT
source("http://zzlab.net/GAPIT/gapit_functions.txt")
setwd("/fs4/Mining_SNP/Borassus_flabellifer_Linn/Germplasm/VCF/diversity/GWAS/")

pheno.file <- "trait.txt"
pheno.dat <- read.table(pheno.file,header = T)
colnames(pheno.dat) <- c("Sample","phenotype")

# read hapmap from TASSEL
hmp <- read.delim("final_analysis_set_converted.hmp.txt",header = F)

################################################################################
### GWAS using GAPIT 
### manual - https://zzlab.net/GAPIT/gapit_help_document.pdf
### github - https://github.com/jiabowang/GAPIT?tab=readme-ov-file#installing-gapit
################################################################################

# run GWAS using GAPIT
myGAPIT <- GAPIT(
  Y=pheno.dat,
  G=hmp,
  model=c("GLM","MLM","FarmCPU","Blink"),
  PCA.total=3,
  Random.model = F
)
