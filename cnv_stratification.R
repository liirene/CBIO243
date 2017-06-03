setwd("~/Desktop/CBIO243")
source("https://bioconductor.org/biocLite.R")

fileDirectory <- "~/Desktop/CBIO243/"
cancersite <- "LUAD"
targetDirectory <- paste0(fileDirectory, cancersite, "/")
setwd(targetDirectory)

### Load in files
cnv <- read.table(paste0(targetDirectory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes"),
                       header=T, stringsAsFactors = F, sep = "\t", row.names=1)
cnv_thresholded <- read.table("Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
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

# Frequency thresholds (total 516 patients)
# Right now operating on >150?? 

driver_amp <- c()
driver_del <- c()

dataset <- cnv_thresholded

for (i in 1:nrow(dataset)){  
  x <- count(as.numeric(as.vector(dataset[i,])))
  max_freq <- x[which(x$freq == max(x$freq)),]
  
  # Note: if two frequencies are the same, any() will help test if either of
  # them are above zero. 
  # Currently this is slightly biased in favor of filtering things out if a
  # lack of change is largely represented.
  
  if (any(max_freq$x == 0)){
    next
  }
  else if (any(max_freq$x > 0)){
    driver_amp <- c(driver_amp, rownames(dataset)[i])
  }
  else if (any(max_freq$x < 0)){
    driver_del <- c(driver_del, rownames(dataset)[i])
  }
}

driver_cnv <- c(driver_amp, driver_del)
write.csv(driver_cnv, "driver_cnv.csv")

