depth <- Adiv.abiotic2$depth
DIphyto <- cbind(depth, test2)
DIphyto$STATION <- as.list(DIphyto$STATION)

# number of unique stations and corresponding rows
SL <- (unique(DIphyto$STATION))

phytoInt <- matrix(0,length(SL),72)

for(i in 1:length(SL)){
print(i)
  Station <- toString(SL[i])
  n <- nrow(DIphyto[which(DIphyto$STATION==Station),])
  data <- as.matrix(DIphyto[which(DIphyto$STATION==Station),])
  for (z in 1:(n-1)){
    print(z)
    z1=z+1
    dz <- data[2,(z1)]-data[1,z]
    print(dz)
    phyto <- rowSums(data[z:z1, 4:75])
    phyto+phytop
    phytoInt[i,] <- as.matrix(phytoInt[i,] + (dz*phyto))
  print(z)
}
print(i)
}

#practice version

phytoInt <- matrix(0,length(SL),72)

for(i in 1:length(SL)){
  print(i)
  Station <- toString(SL[i])
  n <- nrow(DIphyto[which(DIphyto$STATION==Station),])
  data <- as.matrix(DIphyto[which(DIphyto$STATION==Station),])}
for (z in 1:(n-1)){
  print(z)
  z1=z+1
  dz <- data[2,(z1)]-data[1,z]
  print(dz)
  phyto <- (data[z1, 3:75]+data[z, 3:75])/2
  phytoInt[i,] <- as.matrix(phytoInt[i,] + (dz*phyto))
  print(z)
}
print(i)
}

#apply(test[SL], 3:75], 1, function(x) sum(x>0))