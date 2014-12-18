S.A <- div.abiotic[1:39,]
S.B <- div.abiotic[40:79,]
S.C <- div.abiotic[80:121,]
S.D <- div.abiotic[122:152,]
S.E <- div.abiotic[153:19,]
GenA<- data.frame(Means=colMeans(S.A[, 6:9]), Median=apply(S.A[, 6:9], 2, median), Sd=apply(S.A[, 6:9], 2, sd)) 
GenB<- data.frame(Means=colMeans(S.B[, 6:9]), Median=apply(S.B[, 6:9], 2, median), Sd=apply(S.B[, 6:9], 2, sd)) 
GenC<- data.frame(Means=colMeans(S.C[, 6:9]), Median=apply(S.C[, 6:9], 2, median), Sd=apply(S.C[, 6:9], 2, sd)) 
GenD<- data.frame(Means=colMeans(S.C[, 6:9]), Median=apply(S.C[, 6:9], 2, median), Sd=apply(S.C[, 6:9], 2, sd))
GenE<- data.frame(Means=colMeans(S.E[, 6:9]), Median=apply(S.E[, 6:9], 2, median), Sd=apply(S.E[, 6:9], 2, sd)) #no fruitful results 4 now
UpDepthC <- S.C[which(S.C[, 3]<51),] #chose C as most promising variability due to middle position. Selected only 0-50m samples
SurfC <- S.C[which(S.C[, 3]<1),] #Surface
MiddleC <- UpDepthC[which(UpDepthC[,3] >0),]
LowC <- UpDepthC[which(UpDepthC[, 3]>31),]