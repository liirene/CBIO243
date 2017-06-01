#### Clustering of gene expression
install.packages("matrixStats", "ggplot2")
library(matrixStats)
library(ggplot2)

# Use common tumor patients gene expression values, generated previously

# Stratify out genes with highest variance 
ngenes <- 5000

# How to get genes with highest variance? 
# Is it based on the comparison to normal? 

geneexp_common_vars <- data.frame(rowVars(geneexp_common_pts), row.names = rownames(geneexp_common_pts))
colnames(geneexp_common_vars) <- "vars"

ordered <- geneexp_common_vars[order(-geneexp_common_vars$vars), , drop = F]
var_genelist <- rownames(geneexp_common_vars)[1:5000]

geneexp_common_pts.vars <- geneexp_common_pts[var_genelist,]
set.seed(40)
ncenters <- 100

fit_km <- kmeans(geneexp_common_pts.vars, centers = ncenters, nstart=20, iter.max = 50)
km_means <- aggregate(geneexp_common_pts.vars, by=list(cluster = fit_km$cluster), mean)

km_sizes <- data.frame(cbind(as.character(1:100), fit_km$size))
colnames(km_sizes) <- c("cluster", "elements")
km_sizes$elements <- as.numeric(as.character(km_sizes$elements))
km_sizes$cluster <- as.character(km_sizes$cluster)

# Filter out modules with small number of elements per cluster? 
nfilters <- 10
large_modules <- km_sizes[which(km_sizes$elements >= nfilters),]

for (i in 1:nrow(large_modules)){
  large_modules$clustername[i] <- paste0("cluster_", large_modules$cluster[i])
}

# Compose list of elements w/in each larger cluster
clustergenes <- vector("list", length(large_modules$cluster))
names(clustergenes) <- large_modules$clustername

km_elements_all <- fit_km$cluster

for (i in 1:length(clustergenes)){
  a <- list(clusternumber = as.numeric(large_modules$cluster[i]),
            elements = c(), drivers = c())
  clustergenes[[i]] <- a
  clustergenes[[i]]$elements <- names(km_elements_all[which(km_elements_all == clustergenes[[i]]$clusternumber)]) # all that match !! a vector of names
}  
