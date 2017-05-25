setwd("~/Desktop/CBIO243")

# Function initialization
source("https://bioconductor.org/biocLite.R")
biocLite("MethylMix")
library(MethylMix)

# Initial data from Xena
methyl.xena <- read.table("~/Desktop/CBIO243/Xena/HumanMethylation450", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
cnv_threshold.xena <- read.table("~/Desktop/CBIO243/Xena/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
cnv.xena <- read.table("~/Desktop/CBIO243/Xena/Gistic2_CopyNumber_Gistic2_all_data_by_genes", header=T, stringsAsFactors = F, sep = "\t", row.names=1)
geneexp.xena <- read.table("~/Desktop/CBIO243/Xena/HiSeqV2", header=T, stringsAsFactors = F,  sep = "\t", row.names = 1)
LUNG_clinicalMatrix <- read.table("~/Desktop/CBIO243/Xena/LUNG_clinicalMatrix", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)

# Change all "-" to "." for TCGA patient IDs
rownames(LUNG_clinicalMatrix) <- gsub("-",".",rownames(LUNG_clinicalMatrix)) 

# Separate out methylmix tumor/normal samples
tumor_pt_ids <- LUNG_clinicalMatrix[LUNG_clinicalMatrix$sample_type=="Primary Tumor",]
normal_pt_ids <- LUNG_clinicalMatrix[LUNG_clinicalMatrix$sample_type=="Solid Tissue Normal",]


source('~/Desktop/TCGA-Assembler/Module_A.R')
DownloadBiospecimenClinicalData(cancerType = "LUAD", outputFileName = "Clinicaldata")