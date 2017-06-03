#### Clustering of gene expression
install.packages("matrixStats", "ggplot2")
library(matrixStats)
library(ggplot2)

# Use common tumor patients gene expression values, generated previously

# Stratify out genes with highest variance 
ngenes <- 5000

# How to get genes with highest variance? 
# Is it based on the comparison to normal? 

geneexp_common_vars <- data.frame(rowVars(geneexp_common_pts), row.names = 
                                    rownames(geneexp_common_pts))
colnames(geneexp_common_vars) <- "vars"

ordered <- geneexp_common_vars[order(-geneexp_common_vars$vars), , drop = F]
var_genelist <- rownames(geneexp_common_vars)[1:5000]

geneexp_common_pts.vars <- geneexp_common_pts[var_genelist,]

# Set variables for clustering
set.seed(40) # For randomization reproducibility
ncenters <- 200 # How many clusters to make

fit_km <- kmeans(geneexp_common_pts.vars, centers = ncenters, nstart=20, iter.max = 50)

# Calculate KM means... maybe delete this later since using PCA instead
km_means <- aggregate(geneexp_common_pts.vars, by=list(cluster = fit_km$cluster), mean)

km_sizes <- data.frame(cbind(as.character(1:100), fit_km$size))
colnames(km_sizes) <- c("cluster", "elements")
km_sizes$elements <- as.numeric(as.character(km_sizes$elements))
km_sizes$cluster <- as.character(km_sizes$cluster)

# Filter out modules with small number of elements per cluster? 
nfilters <- 60 # minimum number of elements that must be present per cluster
large_modules <- km_sizes[which(km_sizes$elements >= nfilters),]

for (i in 1:nrow(large_modules)){
  large_modules$clustername[i] <- paste0("cluster_", large_modules$cluster[i])
}

# Compose list of elements w/in each larger cluster
clustergenes <- vector("list", length(large_modules$cluster))
names(clustergenes) <- large_modules$clustername

km_elements_all <- fit_km$cluster # Which genes are in which cluster... just pulling this out

##########################################################################
#### Variable determination for glmnet####
# Original variables 
n_lamb <- 100
n_alpha <- 1
n_stringency <- 10 # max number of coefficients (drivers) per module
# Where alpha = the elastic net mixing parameter. 0 <= a <= 1. The definition 
# of the penalty: closer to 0, ridge. CLoser to 1, lasso. 
# Q: Which is more stringent? Ridge uses L2 penalty: limits the size of coefficient. 
# Lasso is L1, which imposes sparsity among the coefficients. Thus **elastic net** is
# a mixture of the two systems (via how you set alpha)

# More information: 
# https://stats.stackexchange.com/questions/93181/ridge-lasso-and-elastic-net

for (i in 1:length(clustergenes)){
  a <- list(clusternumber = as.numeric(large_modules$cluster[i]),
            elements = c(), geneexp = data.frame(), pc1 = data.frame(), possible_drivers = c())
  clustergenes[[i]] <- a
  
  genes_included <- names(km_elements_all[which(km_elements_all == clustergenes[[i]]$clusternumber)]) # all that match !! a vector of names
  clustergenes[[i]]$elements <- genes_included
  
  # Principal Component Analysis
  b <- prcomp(t(geneexp_common_pts[genes_included,]), scale=F, center=F) # Not sure if this will change anything?
  c <- as.matrix(b$x[,1])
  clustergenes[[i]]$pc1 <- c # Pull out first PC

  # GLMNET
  test <- glmnet(t(driver_geneexp), t(c), alpha = n_alpha, nlambda = n_lamb) 
  test_all_coefs <- coef(test)
  assoc_drivers <- rownames(test_all_coefs[which(test_all_coefs[,n_stringency]!=0),])
  clustergenes[[i]]$possible_drivers <- assoc_drivers[2:length(assoc_drivers)] # don't include intercept 
    
  # include gene expression of all the genes in each cluster... might not need this? 
  clustergenes[[i]]$geneexp <- geneexp_common_pts[genes_included,]
}  
# any log transformations?

# Where a is the first principal comp currently
# test now contains the glmnet object


# coefs contain all the coefficients for each observation (x) to generate the 
# response (y) coefficients that contribute to the sparse model
# Usually with the sparse model, lambda 1 has the fewest... What degree of 
# Lambda to choose? 

# Each combination of lambda is a different set of coefficients (s0, etc)





