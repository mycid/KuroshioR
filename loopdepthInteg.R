depth <- as.numeric(Adiv.abiotic2$depth)
DIphyto <- cbind(depth, test2)
DIphyto$STATION <- as.list(DIphyto$STATION)

# number of unique stations and corresponding rows
SL <- (unique(DIphyto$STATION))

phytoInt <- matrix(0,length(SL),72)

for(i in 1:length(SL)){
print(i)
  Station <- toString(SL[i])
  n <- nrow(DIphyto[which(DIphyto$STATION==Station),])
  data <-DIphyto[which(DIphyto$STATION==Station),]
  for (z in 1:(n-1)){
    print(z)
    z1=z+1
    dz <- data[z,1]-data[z1,1] #apply(data[1,z], 2, function(x) x-x2)
    print(dz)
    phyto <- (data[z1, 4:75]+data[z, 4:75])/2
    phytoInt[i,] <- as.matrix(phytoInt[i,] + (-dz*phyto))
  print(z)
}
print(i)
}
DIp <- rbind(Adiv.abiotic2[which(Adiv.abiotic2[,5]=="S"),], Adiv.abiotic2[2,])
DIp <- DIp[,-c(5,6,9:38)]
DIp <- cbind(DIp, phytoInt)
S <- c(Surface$S, Adiv.abiotic2[31,16])
T.C <- c(Surface$T.C., Adiv.abiotic2[31,15])
SpR <- apply(DIp[7:78], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
SpA <- apply(DIp[7:78], 1, function(x) sum(x))
SD <- apply(DIp[7:78], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(DIp[7:78], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness 
DiversityI <- data.frame(Richness=SpR, Cellcount=SpA, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
IntDiverse <- cbind(DIp[,1:6], DiversityI, S, T.C)
