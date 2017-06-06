#### Functional KEGG-based annotation of gene sets. 

source("https://bioconductor.org/biocLite.R")
biocLite(c("pathview", "gageData", "gage"))

## Functionalization. 

pathview_annotation <- function{
   require(pathview)
  require(gageData)
  require(gage)
}
 
