#### Clustering of gene expression
install.packages("matrixStats", "ggplot2")
library(matrixStats)
library(ggplot2)

##########################################################################
#### Function calculations

Calculate_corravgs <- function(matrix){
  # Where matrix is the gene expression values AFTER transposition
  x <- cor(matrix)
  
  # average calculation
  result <- mean(x[upper.tri(x)])
  return(result)
}

##########################################################################
# Use common tumor patients gene expression values, generated previously

# Stratify out genes with highest variance 
ngenes <- 5000

# How to get genes with highest variance? 
# Is it based on the comparison to normal? 

geneexp_common_vars <- data.frame(rowVars(geneexp_common_pts), row.names = 
                                    rownames(geneexp_common_pts))
colnames(geneexp_common_vars) <- "vars"

ordered <- geneexp_common_vars[order(-geneexp_common_vars$vars), , drop = F]
var_genelist <- rownames(ordered)[1:5000]

geneexp_common_pts.vars <- geneexp_common_pts[var_genelist,]

##########################################################################

#### Clustering running

# Set variables for clustering
set.seed(40) # For randomization reproducibility s

ncenters <- 100 # How many clusters to make (change this later)

fit_km <- kmeans(geneexp_common_pts.vars, centers = ncenters, nstart=20, iter.max = 50)

# Calculate KM means... maybe delete this later since using PCA instead
km_means <- aggregate(geneexp_common_pts.vars, by=list(cluster = fit_km$cluster), mean)


##########################################################################

#### DELETE THIS LATER? 
km_sizes <- data.frame(cbind(as.character(1:ncenters), fit_km$size))
colnames(km_sizes) <- c("cluster", "elements")
km_sizes$elements <- as.numeric(as.character(km_sizes$elements))
km_sizes$cluster <- as.character(km_sizes$cluster)

##########################################################################

### Cluster filtering

# Obtain coexpression: take the average of correlation coefficients between gene
# pairs and take the average to get average clusters. 

#### Co-expression between gene clusters calculation ( put them into the km_sizes)

km_coexp <- data.frame(cbind(as.character(1:ncenters), rep(NA, ncenters)))
colnames(km_coexp) <- c("cluster", "avg_correlation")
km_coexp$avg_correlation <- as.numeric(as.character(km_coexp$avg_correlation))
km_coexp$cluster <- as.character(km_coexp$cluster)

km_elements_all <- fit_km$cluster # Which genes are in which cluster... just pulling this out

gene_exps_all <- vector("list", ncenters)
names(gene_exps_all) <- paste("cluster", 1:ncenters)

for (i in 1:ncenters) {
  genes_included <- names(km_elements_all[which(km_elements_all == i)]) 
  
  exp_matrix <- t(geneexp_common_pts[genes_included,])
  
  result <- Calculate_corravgs(exp_matrix)
  
  km_coexp$avg_correlation[i] <- result
}


# Genes that have high coexpression are highly coregulated statistically and if 
# they are also functionally enriched, we know this is potentiall an interesting
# process. 

# Filter out modules with small number of elements per cluster?

coexp_filter <- 0.5
large_modules <- km_coexp[which(km_coexp$avg_correlation >= coexp_filter),]

#nfilters <- 60 # minimum number of elements that must be present per cluster
#large_modules <- km_sizes[which(km_sizes$elements >= nfilters),]

for (i in 1:nrow(large_modules)){
  large_modules$clustername[i] <- paste0("cluster_", large_modules$cluster[i])
}

# Compose list of elements w/in each larger cluster
clustergenes <- vector("list", length(large_modules$cluster))
names(clustergenes) <- large_modules$clustername



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

### Lambda testing (do this within the loop?)

# This is to define stringency for the model

##########################################################################
for (i in 1:length(clustergenes)){
  a <- list(clusternumber = as.numeric(large_modules$cluster[i]),
            elements = c(), geneexp = data.frame(), pc1 = data.frame(), possible_drivers = c())
  
  clustergenes[[i]] <- a
  
  genes_included <- names(km_elements_all[which(km_elements_all == clustergenes[[i]]$clusternumber)]) 
  clustergenes[[i]]$elements <- genes_included
  
  # Principal Component Analysis
  b <- prcomp(t(geneexp_common_pts[genes_included,]), scale=F, center=F) 
  # Not sure if the vars included will change anything?
  # Pull out first PC
  c <- as.matrix(b$x[,1])
  clustergenes[[i]]$pc1 <- c 
  
  # GLMNET
  d <- glmnet(t(driver_geneexp), t(c), alpha = n_alpha, nlambda = n_lamb) 
  test_all_coefs <- coef(d)
  
  assoc_drivers <- rownames(test_all_coefs[which(test_all_coefs[,n_stringency]!=0),])
  clustergenes[[i]]$possible_drivers <- assoc_drivers[2:length(assoc_drivers)] 
  # don't include the intercept column
  
  # include gene expression of all the genes in each cluster
  clustergenes[[i]]$geneexp <- geneexp_common_pts[genes_included,]
}  


# coefs contain all the coefficients for each observation (x) to generate the 
# response (y) coefficients that contribute to the sparse model
# Usually with the sparse model, lambda 1 has the fewest... What degree of 
# Lambda to choose? 

# Each combination of lambda is a different set of coefficients (s0, etc)





