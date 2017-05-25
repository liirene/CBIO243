setwd("~/Desktop/243datasets")

# Read in initial data
methyl <- read.table("~/Desktop/243datasets/TCGA_LUNG_hMethyl450-2015-02-24/genomicMatrix", header=T, stringsAsFactors = F, sep = "\t")
cnv_threshold <- read.table("~/Desktop/243datasets/TCGA_LUNG_gistic2thd-2015-02-24/genomicMatrix", header=T, stringsAsFactors = F, sep = "\t")
cnv_threshold <- read.table("~/Desktop/243datasets/TCGA_LUNG_gistic2-2015-02-24/genomicMatrix", header=T, stringsAsFactors = F, sep = "\t")
genexp <- read.table("~/Desktop/243datasets/TCGA_LUNG_exp_HiSeqV2-2015-02-24/genomicMatrix", header=T, stringsAsFactors = F,  sep = "\t")

source('~/Desktop/TCGA-Assembler/Module_A.R')
DownloadBiospecimenClinicalData(cancerType = "LUAD", outputFileName = "Clinicaldata")
DownloadM