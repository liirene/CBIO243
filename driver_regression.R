#### Sparse linear regression with glmnet
install.packages("glmnet")
library(glmnet)

# Building gene expression for all drivers

test <- geneexp_common_pts[driver_final,]

glmnet(test,clustergenes[[1]]$geneexp)
