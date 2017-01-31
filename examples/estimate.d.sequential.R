n.muts <- length(Khan.data[1,])-1
geno.matrix <- Khan.data[,seq(1, n.muts)]
fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))
d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix,-100)$d.hat.RDB
estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, c(0.1, 10), 1.1)
