### Build list of drivers

fileDirectory <- "~/Desktop/CBIO243/"
cancersite <- "LUAD"
targetDirectory <- paste0(fileDirectory, cancersite, "/")
setwd(targetDirectory)

driver_methyl <- read.table("methyl_drivers.csv", header = T, 
                            stringsAsFactors = F, sep = ",",
                            colClasses = c("NULL", "character"))

driver_methyl <- driver_methyl[[1]]

# All lists are combined in driver_methyl, amp (cnv), del (cnv), and linreg.
# Thus need to find intersect between all. 

a <- intersect(driver_linreg, driver_methyl)
b <- intersect(driver_linreg, driver_amp)
c <- intersect(driver_linreg, driver_del)

driver_final <- c(a, b, c)

# Build driver gene expression matrix 
driver_geneexp <- as.matrix(geneexp_common_pts[driver_final,])

write.csv(driver_geneexp, file = paste0("driver_geneexp_Rsq", Rsquared_filter, "_", cancersite, ".csv"))
