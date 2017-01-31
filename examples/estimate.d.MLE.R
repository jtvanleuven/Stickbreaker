n.muts <- length(Khan.data[1,])-1
geno.matrix <- Khan.data[,seq(1, n.muts)]
fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))
