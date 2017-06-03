### Linear regression for gene expression and cnvs. 

########################################################################

# Function to extract the overall ANOVA p-value out of a linear model object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

########################################################################
### Build gene expression matrix of only tumor patients that are also present in
### methylation and CNV

methyl.xena <- read.table(paste0(targetDirectory, "HumanMethylation450"), header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
geneexp.xena <- read.table(paste0(targetDirectory, "HiSeqV2"), header=T, stringsAsFactors = F,  sep = "\t", row.names = 1)
clinicalMatrix <- read.table(paste0(targetDirectory, cancersite, "_clinicalMatrix"), header=T, stringsAsFactors = F, sep = "\t", row.names = 1)
cnv <- read.table(paste0(targetDirectory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes"), header=T, stringsAsFactors = F, sep = "\t", row.names = 1)

pt_all <- rownames(clinicalMatrix)
pt_methyl <- colnames(methyl_genename[,-1])
pt_geneexp <- colnames(geneexp.xena)
pt_cnv <- colnames(cnv)

# Include the intersections SOON
common1 <- intersect(pt_all, pt_methyl)
common2 <- intersect(pt_geneexp, common1) 
common3 <- intersect(pt_cnv, common2) # CNV only includes tumor patients

common_normal <- intersect (common2, pt_ids_normal) # total 21 in LUAD
common_tumor <- intersect (common3, pt_ids_tumor) # total 451 in LUAD


# Build tumor methylation data: take patient samples that are common between
# gene/methylation datasets AND from Solid Tumors
methyl_tumor <- cbind(methyl_genename[,1], 
                      methyl_genename[na.omit(match(common_tumor, 
                                                    colnames(methyl_genename)))])

# Delete NA rows in datasets to be tested. Geneexp doesn't need to be cleaned.
methyl_tumor.naclean <- methyl_tumor[complete.cases(methyl_tumor),]

# Average values for each duplicated Cpg island value 
tumor_levels <- methyl_tumor.naclean$`methyl_genename[, 1]`
methyl_tumor.naclean.agg <- aggregate(x = methyl_tumor.naclean[,-1], list(tumor_levels), FUN = "mean")
#20190 entries
methyl_tumor.naclean.agg <- na.omit(methyl_tumor.naclean.agg[match(methyl_tumor.naclean.agg$Group.1, common_gene),])
# 16142

# Make into matrices with correct rownames
rownames(methyl_tumor.naclean.agg) <- methyl_tumor.naclean.agg$Group.1

# Make all into matrices

# Common genes b/t methylation and gene expression 
common4 <- intersect(rownames(geneexp.xena), rownames(cnv))
common_gene <- intersect(common4, rownames(methyl_tumor.naclean.agg))

# Remove non-character cols
methyl_tumor_common_pts <- data.matrix(methyl_tumor.naclean.agg[common_gene,c(-1)])

### Build gene exp and cnv datasets
geneexp_common_pts <- data.matrix(geneexp.xena[common_gene, common_tumor])
cnv_common_pts <- data.matrix(cnv[common_gene, common_tumor])

# LM With a loop. 
lm_results <- list()
lm_values <- data.frame(matrix(data=NA, nrow = nrow(cnv_common_pts), ncol=3))
rownames(lm_values) <- rownames(cnv_common_pts)
colnames(lm_values) <- c("rsquared", "methyl_pval", "cnv_pval") 

for (i in 1:nrow(cnv_common_pts)){
 result <- lm(geneexp_common_pts[i,] ~ 0 
                     + methyl_tumor_common_pts[i,] 
                     + cnv_common_pts[i,], 
                     model = T)
 lm_results[[i]] <- result # for reference
 lm_values$rsquared[i] <- summary(result)$r.squared
 lm_values$methyl_pval[i] <- summary(result)$coefficients[1,4]
 lm_values$cnv_pval[i] <- summary(result)$coefficients[2,4]
}

# Order list and draw out top 200 in R-squared for candidate driver "List 1"
# There is some discordance in the methyl genes and the rest of the genes. Try to find out what is missing... 
lm_values <- lm_values[order(-lm_values$rsquared),]
write.csv(lm_values, file=paste0(cancersite, "_lm_values.csv"))

# Write list of genes with R^2 values above 0.9
Rsquared_filter <- 0.95
driver_linreg <- row.names(na.omit(lm_values[lm_values$rsquared >= Rsquared_filter,]))