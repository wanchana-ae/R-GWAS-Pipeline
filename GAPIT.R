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
#install.packages("remotes")
#remotes::install_github('coolbutuseless/svgparser')

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
BiocManager::install("ggtree")
#BiocManager::install("SNPRelate")

#Call tools
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

#GAPIT
source("http://zzlab.net/GAPIT/gapit_functions.txt")


#setwd("C:/Users/wanch/Downloads/work")
################################################################################
# input phenotype files
################################################################################
vcf.filename <- "snp_set.recode.vcf"
vcf <- read.vcfR(vcf.filename)

pheno.file <- "trait.txt"
pheno.dat <- read.table(pheno.file,header = T)
colnames(pheno.dat) <- c("Sample","phenotype")

################################################################################
### Phenotype histogram
################################################################################
ggplot(pheno.dat, aes(x = phenotype)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  labs(title = "Histogram of Phenotype",
       x = "Phenotype Value",
       y = "Frequency")+
  theme_minimal()

################################################################################
### Phylogenetic tree
################################################################################
#convert VCF to GDS
snpgdsVCF2GDS(vcf.filename, "genotype.gds",  method="biallelic.only")
genofile <- snpgdsOpen("genotype.gds")
#create Identity-By-State (IBS) matrix from GDS
ibs <- snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))

#Cut dendrogram to create group (clusters)
rvExon <- snpgdsCutTree(ibs)

#Plot dendrogram
plot(rvExon$dendrogram,horiz=T, main ="Genotypes SNP Tree")

#Get dendrogram snpgdsCutTree
tree <- rvExon$dendrogram

#Convert dendrogram -> hclust -> phylo
phylo_tree <- as.phylo(as.hclust(tree))

#Plot
tree2 = ggtree(phylo_tree,
               layout="circular",color='darkgreen', 
               branch.length="branch.length")+
  geom_tiplab(aes(label = label), parse = FALSE)+ 
  ggtitle("Genotypes SNP Tree2")

tree2
ggplot(phylo_tree, aes(x, y)) + geom_tree() + theme_tree()
ggtree(phylo_tree, layout="circular")

## Writing to a newick tree file
write.tree(phy=phylo_tree , file="phylo_tree.newick")

#Plot tree https://itol.embl.de/

################################################################################
### PCA
################################################################################
pca<-snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2)) #check Percent PC
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

gg<-ggplot(tab,aes(x=EV1,y=EV2))+
  geom_point(size=2,shape=16,col="blue")+
  labs(subtitle="", 
       x= paste("PC1:",sprintf("%.2f%%",pc.percent[1])),
       y= paste("PC2:",sprintf("%.2f%%",pc.percent[2])),
       title="Principal Component Analysis (PCA)")
gg
theme_set(theme_bw())
options(ggrepel.max.overlaps = Inf) #For label overlaps
gg + geom_label_repel(aes(label = sample.id),
                      size=3,
                      box.padding   = 0.1, 
                      point.padding = 0.1,
                      segment.color = 'grey50')
## Writing to a CSV file
write.csv(tab,"PCA.csv")

################################################################################
### GWAS using GAPIT 
### manual - https://zzlab.net/GAPIT/gapit_help_document.pdf
### github - https://github.com/jiabowang/GAPIT?tab=readme-ov-file#installing-gapit
################################################################################
# convert vcf to hapmap format
hmp <- vcfR2hapmap(vcf)
head(pheno.dat)
# read hapmap from TASSEL
#hmp <- read.delim("file.hmp.txt", header = F)

# run GWAS using GAPIT
myGAPIT <- GAPIT(
  Y=pheno.dat[,c("Trait","Seed_color_DPC_1_11")],
  G=hmp,
  model="FarmCPU",
  PCA.total=3,
  Random.model = F
)
myGAPIT$GWAS[which(myGAPIT$GWAS$P.value == min(myGAPIT$GWAS$P.value)),]

#sort by P.value
sorted_gwas <- myGAPIT$GWAS[order(myGAPIT$GWAS$P.value), ]

# add -log10(P.value)
sorted_gwas$logP <- -log10(sorted_gwas$P.value)

head(sorted_gwas)

#select 100 line
top100_gwas <- sorted_gwas[1:100, ]

# save CSV
write.csv(top100_gwas, "top100_gwas_logP.csv", row.names = FALSE)

################################################################################
### extract only Chromosome 12
################################################################################
vcf.chr12 <- vcf[which(vcf@fix[,"CHROM"] == 12),]
vcf.chr12
write.vcf(vcf.chr12, "tomato_chr12.vcf.gz")

#https://r-graph-gallery.com/101_Manhattan_plot.html
library(qqman)
gwas_Results <- myGAPIT$GWAS
head(gwas_Results)

#covert to numbers
gwas_Results$Chromosome <- as.numeric(gwas_Results$Chromosome)
gwas_Results$Position <- as.numeric(gwas_Results$Position)

# Make the Manhattan plot
manhattan(gwas_Results, chr="Chromosome", bp="Position", snp="SNP", p="P.value" )

#QQ plot
qq(gwas_Results$P.value)
