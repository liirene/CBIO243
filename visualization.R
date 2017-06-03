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

colnames(nodes) <- c("module", "driver")

write.csv(nodes, paste0(cancersite,"_filters", nfilters, 
                        "_coef", n_stringency, 
                        "_clusters", ncenters, "_nodes.csv"))
