setwd("~/Desktop/CBIO243")

########################################################################
###Initializing files

#Packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("FDb.InfiniumMethylation.hg19", "MethylMix"))
library(MethylMix)

cancersite <- "LUAD"
targetDirectory <- "~/Desktop/CBIO243/Xena"

# Initial data from Xena
methyl.xena <- read.table("~/Desktop/CBIO243/Xena/HumanMethylation450", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
cnv_threshold.xena <- read.table("~/Desktop/CBIO243/Xena/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
cnv.xena <- read.table("~/Desktop/CBIO243/Xena/Gistic2_CopyNumber_Gistic2_all_data_by_genes", header=T, stringsAsFactors = F, sep = "\t", row.names=1)
geneexp.xena <- read.table("~/Desktop/CBIO243/Xena/HiSeqV2", header=T, stringsAsFactors = F,  sep = "\t", row.names = 1)
LUNG_clinicalMatrix <- read.table("~/Desktop/CBIO243/Xena/LUNG_clinicalMatrix", header=T, stringsAsFactors = F, sep = "\t", row.names = 1)

# Change all "-" to "." for TCGA patient IDs
rownames(LUNG_clinicalMatrix) <- gsub("-",".",rownames(LUNG_clinicalMatrix)) 

# Separate out methylmix tumor/normal samples
tumor_pt_ids <- rownames(LUNG_clinicalMatrix[LUNG_clinicalMatrix$sample_type=="Primary Tumor",])
normal_pt_ids <- rownames(LUNG_clinicalMatrix[LUNG_clinicalMatrix$sample_type=="Solid Tissue Normal",])

########################################################################

### Use annotation data to convert CpG site to closest gene name (from Wei)

library(FDb.InfiniumMethylation.hg19)

hm450 <- get450k()
probenames <- as.character(row.names(methyl.xena))
probes <- hm450[probenames]
mapping <- getNearestTranscript(probes)
nearestGeneSymbol <- mapping[,4]

# Total gene-matched data
LUNG_meth_genename <- cbind(nearestGeneSymbol, methyl.xena[,2:dim(methyl.xena)[2]])

########################################################################

### Data concordance b/t methylation and RNAseq

pt_all <- rownames(LUNG_clinicalMatrix)
pt_methyl <- colnames(methyl.xena)
pt_geneexp <- colnames(geneexp.xena)

common1 <- intersect(pt_all, pt_methyl)
common_pt <- intersect(pt_geneexp, common1) # Total 477 IDs b/t gene exp, clinical, and methyl in LUNG

common_normal <- intersect (common_pt, normal_pt_ids) # total 21
common_tumor <- intersect (common_pt, tumor_pt_ids) # total 454

# Common genes b/t methylation and gene expression 
common_gene <- intersect(LUNG_meth_genename$nearestGeneSymbol, rownames(geneexp.xena))

# Build tumor gene expression data
geneexp_common_pts <- geneexp.xena[common_gene, common_pt]

# Build tumor methylation data: take patient samples that are common between
# gene/methylation datasets AND from Solid Tumors
methyl_tumor <- cbind(LUNG_meth_genename[,1], 
                      LUNG_meth_genename[na.omit(match(common_tumor, 
                                                       colnames(LUNG_meth_genename)))])
# Cannot have NAs. 

# Build normal methylation data: take patient samples that are common between
# gene/methylation datasets AND from Normal Solid tissue
methyl_normal <- cbind(LUNG_meth_genename[,1], LUNG_meth_genename[,common_normal])

# Delete NA rows in datasets to be tested. Geneexp doesn't need to be cleaned.
methyl_tumor.naclean <- methyl_tumor[complete.cases(methyl_tumor),]
methyl_normal.naclean <- methyl_normal[complete.cases(methyl_normal),]

# Average values for each duplicated Cpg island value 
tumor_levels <- as.factor(methyl_tumor.naclean$LUNG_meth_genename)
methyl_tumor.naclean.agg <- aggregate(x = methyl_tumor.naclean, list(tumor_levels), FUN = "mean")

normal_levels <- as.factor(methyl_normal.naclean$LUNG_meth_genename)
methyl_normal.naclean.agg <- aggregate(x = methyl_normal.naclean, list(normal_levels), FUN = "mean")

# Make into matrices with correct rownames
rownames(methyl_tumor.naclean.agg) <- methyl_tumor.naclean.agg$Group.1
rownames(methyl_normal.naclean.agg) <- methyl_normal.naclean.agg$Group.1

# Remove non-character rows
methyl_tumor.naclean.agg <- data.matrix(methyl_tumor.naclean.agg[,c(-1, -2)])
methyl_normal.naclean.agg <- data.matrix(methyl_normal.naclean.agg[,c(-1, -2)])
geneexp_common_pts <- data.matrix(geneexp_common_pts)

# Methylmix run
methylmixresult <- MethylMix(methyl_tumor.naclean.agg, geneexp_common_pts, methyl_normal.naclean.agg)
drivers <- methylmixresult$MethylationDrivers
write.csv(drivers, "drivers.csv")
