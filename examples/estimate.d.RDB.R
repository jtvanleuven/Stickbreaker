n.muts <- length(Khan.data[1,])-1
geno.matrix <- Khan.data[,seq(1, n.muts)]
fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
estimate.d.RDB(geno.matrix, fit.matrix,-100)
