#### Graphing 

install.packages("igraph")
library(igraph)

### Print graphs

nodes <- data.frame(matrix(nrow = 0, ncol = 2))
for (i in 1:length(clustergenes)){
  x <- clustergenes[[i]]$possible_drivers
  n_drivers <- length(x)
  y <- names(clustergenes)[i]
  rep_matrix <- cbind(rep(y, n_drivers), x)
  nodes <- rbind(nodes, rep_matrix)
}

cv.nodes <- data.frame(matrix(NA,0,2))
for (i in 1:length(clustergenes)){
  x <- clustergenes[[i]]$lambda_names
  n_drivers <- length(x)
  y <- names(clustergenes)[i]
  rep_matrix <- cbind(rep(y, n_drivers), x)
  cv.nodes <- rbind(cv.nodes, rep_matrix)
}

all_cluster_elements <- data.frame(matrix(NA,0,2))
for (i in 1:length(clustergenes)){
  x <- clustergenes[[i]]$elements
  n_elements <- length(x)
  y <- names(clustergenes)[i]
  rep_matrix <- cbind(rep(y, n_elements), x)
  all_cluster_elements <- rbind(all_cluster_elements, rep_matrix)
}

colnames(nodes) <- c("module", "driver")
colnames(cv.nodes) <- c("module", "driver")
colnames(all_cluster_elements) <- c("module", "gene")

write.csv(nodes, paste0(cancersite, "_Rsq", Rsquared_filter, "_coexpfilters", coexp_filter,
                        "_clusters", ncenters, "_nodes.csv"))

write.csv(cv.nodes, paste0("CV_", cancersite, "_Rsq", Rsquared_filter, "_coexpfilters", coexp_filter, 
                        "_clusters", ncenters, "_nodes.csv"))

write.csv(all_cluster_elements, paste0(cancersite, "_clusterelements.csv"))
