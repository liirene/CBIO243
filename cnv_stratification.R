setwd("~/Desktop/CBIO243")
source("https://bioconductor.org/biocLite.R")

fileDirectory <- "~/Desktop/CBIO243/"
cancersite <- "LUAD"
targetDirectory <- paste0(fileDirectory, cancersite, "/")
setwd(targetDirectory)

### Load in files
cnv <- read.table(paste0(targetDirectory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes"),
                       header=T, stringsAsFactors = F, sep = "\t", row.names=1)
clinicalMatrix <- read.table(paste0(targetDirectory, cancersite, "_clinicalMatrix"),
                                  header=T, stringsAsFactors = F, sep = "\t", row.names = 1)

rownames(clinicalMatrix) <- gsub("-",".",rownames(clinicalMatrix)) 

### Separate out methylmix tumor/normal samples
pt_ids_tumor <- rownames(clinicalMatrix[clinicalMatrix$sample_type=="Primary Tumor" | clinicalMatrix$sample_type=="Recurrent Tumor",])
pt_ids_normal <- rownames(clinicalMatrix[clinicalMatrix$sample_type=="Solid Tissue Normal",])

# NOTE: GISTIC dataset only contains tumor samples so no separation is needed.

### Stratify CNVs
cnv.tumor <- cnv
cnv.tumor$sums <- rowSums(cnv.tumor)
cnv.tumor$avg <- rowMeans(cnv.tumor)

# Define amplification threshold as the genes which have >1 CNV and are in the top 10% 
# the CNV sums among pts. Deletion threshold as genes which have <-1 CNV and are in the
# bottom 10% of CNV sums among pts.

x <- ncol(cnv.tumor)
y <- x-1

del_avg <- cnv.tumor[order(cnv.tumor$avg),][1:500,y:x]
amp_avg <- cnv.tumor[order(-cnv.tumor$avg),][1:500,y:x]

del_sum <- cnv.tumor[order(cnv.tumor$sums),][1:500,y:x]
amp_sum <- cnv.tumor[order(-cnv.tumor$sums),][1:500,y:x]

driver_del <- rownames(del_avg[match(rownames(del_avg), rownames(del_sum)),])
driver_amp <- rownames(amp_avg[match(rownames(amp_avg), rownames(amp_sum)),])
  
driver_cnv <- c(driver_amp, driver_del)
write.csv(driver_cnv, "driver_cnv.csv")
# In which the first 500 are amp, last 500 are del
  