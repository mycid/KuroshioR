depth <- as.numeric(All$depth)
origin <- OriginalPhytoplankton.csv #csv for all species and groups
#For all groups and species
DIphyto.or <- cbind(depth, origin)
DIphyto.or$STATION <- as.list(DIphyto.or$STATION)
#For just species
depth <- depth[-c(186:190)]
DIphyto <- cbind(depth, test[ -c(186:190), 1:2], centitest[4:ncol(centitest)])
DIphyto$STATION <- as.list(DIphyto$STATION)

# number of unique stations and corresponding rows
SL <- (unique(DIphyto$STATION))

phytoInt <- matrix(0,length(SL),71)

for(i in 1:length(SL)){
  Station <- toString(SL[i])
  n <- nrow(DIphyto[which(DIphyto$STATION==Station),])
  data <-DIphyto[which(DIphyto$STATION==Station),]
  for (z in 1:(n-1)){
    z1=z+1
    dz <- data[z,1]-data[z1,1] #apply(data[1,z], 2, function(x) x-x
    phyto <- (data[z1, 4:74]+data[z, 4:74])/2
    phytoInt[i,] <- as.matrix(phytoInt[i,] + (-dz*phyto))
}
}
DIp <- (All[which(All[, "depth..m."]==0),])
DIp <- DIp[-40,]
phytoInt <- phytoInt[-7,]
PhyInt<- phytoInt
rownames(DIp)<- NULL
IntDiverse <- Create.Diversity.I(phytoInt, 1, 71)
IntBiotic <- cbind(IntDiverse, DIp)
colnames(phytoInt) <- c(as.list(colnames(test[3:73])))
#------------------------phytoInt-------------------
#Now lets cluster the plankton
phytoInt[phytoInt>0]<- 1 
phytoInt <- phytoInt[-40,]
IntD <- vegdist(phytoInt, method="jaccard")
IntH <- hclust(IntD, method="complete")
intclus <- cutree(IntH, 3)
IntBiotic <- cbind(intclus, IntBiotic[-40,])
summary(lm(intclus~S, data=IntBiotic))
ColorDendrogram(IntH, y = intclus, labels = names(intclus), main = "IntBiotic", branchlength = 2)
NPfactorInColor(IntBiotic, xvar="lon", yvar= "lat", IntBiotic$intclus, xlab="Latitude", ylab="Depth", title="IntBiotic")
#

#Cluster the water
Surf <- cbind(DIp$T.C., DIp$S)
Surf <- Surf[-40,]
Surf <- scale(Surf)
SurfD <- vegdist(Surf, method="euclidean")
Surflist <- as.matrix(SurfD)
SurfH <- hclust(SurfD, method="complete")
Surfclus <- cutree(SurfH, 3)
SurfBiotic <- IntBiotic[-40,]
SurfBiotic <- cbind(Surfclus, SurfBiotic)
ColorDendrogram(SurfH, y = Surfclus, labels = names(Surfclus), main = "Surface", branchlength = 2)
NPfactorInColor(SurfBiotic, xvar="lon", yvar= "lat", SurfBiotic$Surfclus, xlab="Latitude", ylab="Longitude", title="Surface")
#Look at species specific to the cluster
phytoInt2 <- as.data.frame(cbind(intclus, PhyInt))
SplitData(phytoInt2, 1, 'PhyInt')

PhyInt1 <- PhyInt1[,colSums(PhyInt1)!=0 ]                  
PhyInt2 <- PhyInt2[,colSums(PhyInt2)!=0 ] 
PhyInt3 <- PhyInt3[,colSums(PhyInt3)!=0 ]

