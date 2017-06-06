########################################################################
### Initializing files

#Packages
source("https://bioconductor.org/biocLite.R")
biocLite(c("FDb.InfiniumMethylation.hg19", "MethylMix"))     
library(MethylMix)

cancersite <- "LUAD"
targetDirectory <- paste0("~/Desktop/CBIO243/", cancersite, "/")
setwd(targetDirectory)

# Initial data from Xena
methyl.xena <- read.table(paste0(targetDirectory, "HumanMethylation450"), header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
geneexp.xena <- read.table(paste0(targetDirectory, "HiSeqV2"), header=T, stringsAsFactors = F,  sep = "\t", row.names = 1)
clinicalMatrix <- read.table(paste0(targetDirectory, cancersite, "_clinicalMatrix"), header=T, stringsAsFactors = F, sep = "\t", row.names = 1)

# Change all "-" to "." for TCGA patient IDs
rownames(clinicalMatrix) <- gsub("-",".",rownames(clinicalMatrix)) 

# Separate out methylmix tumor/normal samples
pt_ids_tumor <- rownames(clinicalMatrix[clinicalMatrix$sample_type=="Primary Tumor",])
pt_ids_normal <- rownames(clinicalMatrix[clinicalMatrix$sample_type=="Solid Tissue Normal",])

########################################################################

### Use annotation data to convert CpG site to closest gene name (from Wei)

library(FDb.InfiniumMethylation.hg19)

hm450 <- get450k()
probenames <- as.character(row.names(methyl.xena))
probes <- hm450[probenames]
mapping <- getNearestTranscript(probes)
nearestGeneSymbol <- mapping[,4]

# Total gene-matched data
methyl_genename <- cbind(nearestGeneSymbol, methyl.xena[,2:dim(methyl.xena)[2]])

########################################################################

### Intersecting patients and genes b/t methylation and RNAseq datasets
pt_all <- rownames(clinicalMatrix)
pt_methyl <- colnames(methyl_genename[,-1])
pt_geneexp <- colnames(geneexp.xena)

common1 <- intersect(pt_all, pt_methyl)
common_pt <- intersect(pt_geneexp, common1) # Total 477 IDs b/t gene exp, clinical, and methyl in LUNG

common_normal <- intersect (common_pt, pt_ids_normal) # total 21
common_tumor <- intersect (common_pt, pt_ids_tumor) # total 454

# Common genes b/t methylation and gene expression 
common_gene <- intersect(methyl_genename$nearestGeneSymbol, rownames(geneexp.xena))

########################################################################
### Building and formatting datasets for Methylmix

# Build tumor gene expression data
geneexp_common_pts <- geneexp.xena[common_gene, common_pt]

# Build tumor methylation data: take patient samples that are common between
# gene/methylation datasets AND from Solid Tumors
methyl_tumor <- cbind(methyl_genename[,1], 
                      methyl_genename[na.omit(match(common_tumor, 
                                                       colnames(methyl_genename)))])
# Cannot have NAs. 

# Build normal methylation data: take patient samples that are common between
# gene/methylation datasets AND from Normal Solid tissue
methyl_normal <- cbind(methyl_genename[,1], methyl_genename[,common_normal])

# Delete NA rows in datasets to be tested. Geneexp doesn't need to be cleaned.
methyl_tumor.naclean <- methyl_tumor[complete.cases(methyl_tumor),]
methyl_normal.naclean <- methyl_normal[complete.cases(methyl_normal),]

# Average values for each duplicated Cpg island value 
tumor_levels <- as.factor(methyl_tumor.naclean$methyl_genename)
methyl_tumor.naclean.agg <- aggregate(x = methyl_tumor.naclean, list(tumor_levels), FUN = "mean")

normal_levels <- as.factor(methyl_normal.naclean$methyl_genename)
methyl_normal.naclean.agg <- aggregate(x = methyl_normal.naclean, list(normal_levels), FUN = "mean")

# Make into matrices with correct rownames
rownames(methyl_tumor.naclean.agg) <- methyl_tumor.naclean.agg$Group.1
rownames(methyl_normal.naclean.agg) <- methyl_normal.naclean.agg$Group.1

# Remove non-character rows
methyl_tumor.naclean.agg <- data.matrix(methyl_tumor.naclean.agg[,c(-1, -2)])
methyl_normal.naclean.agg <- data.matrix(methyl_normal.naclean.agg[,c(-1, -2)])
geneexp_common_pts <- data.matrix(geneexp_common_pts)

########################################################################
### Running MethylMix. Note: May take up to 20 minutes to run, depending
### on how many samples you have! Plan accordingly!

# Methylmix run and save results in file
methylmixresult <- MethylMix(methyl_tumor.naclean.agg, geneexp_common_pts, methyl_normal.naclean.agg)
driver_methyl <- methylmixresult$MethylationDrivers
methyl_states <- methylmixresult$MethylationStates

write.csv(drivers, "methyl_drivers.csv")
write.csv(methyl_states, "methylstates.csv")
