depth <- as.numeric(Adiv.abiotic2$depth)
origin <- read.csv("OriginalPhytoplankton.csv")#csv for all species and groups
#For all groups and species
DIphyto.or <- cbind(depth, origin)
DIphyto.or$STATION <- as.list(DIphyto$STATION)
#For just species
DIphyto <- cbind(depth, test)
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
DIp <- (Adiv.abiotic2[which(Adiv.abiotic2[,5]=="S"),])
DIp <- rbind(DIp[1:6,], Adiv.abiotic2[31,], DIp[7:40,])
rownames(DIp)<- NULL
DIp <- DIp[,-c(5,6,9:38)]
DIp <- cbind(DIp, phytoInt)
sigPoDen <- c(Surface[1:6,18], Adiv.abiotic2[31, 18], Surface[7:40, 18])
S <- c(Surface[1:6,16], Adiv.abiotic2[31, 16], Surface[7:40, 16])
Theta <- c(Surface[1:6,14], Adiv.abiotic2[31, 14], Surface[7:40, 14])
SpR <- apply(DIp[7:77], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
SpA <- apply(DIp[7:77], 1, function(x) sum(x))
SD <- apply(DIp[7:77], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(DIp[7:77], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness 
DiversityI <- data.frame(Richness=SpR, Cellcount=SpA, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
IntDiverse <- cbind(DIp[1:6], DiversityI, S, Theta, sigPoDen)
colnames(phytoInt) <- c(as.list(colnames(test[3:73])))
as.factor(colnames(phytoInt))
#For all species and group
SL <- (unique(DIphyto.or$STATION))

phytoInt.or <- matrix(0,length(SL),79)

for(i in 1:length(SL)){
  Station <- toString(SL[i])
  n <- nrow(DIphyto.or[which(DIphyto.or$STATION==Station),])
  data <-DIphyto.or[which(DIphyto.or$STATION==Station),]
  for (z in 1:(n-1)){
    z1=z+1
    dz <- data[z,1]-data[z1,1] #apply(data[1,z], 2, function(x) x-x2)
    phyto <- (data[z1, 4:82]+data[z, 4:82])/2
    phytoInt.or[i,] <- as.matrix(phytoInt.or[i,] + (-dz*phyto))
  }
}
DIp.or <- (Adiv.abiotic2[which(Adiv.abiotic2[,5]=="S"),])
DIp.or <- rbind(DIp.or[1:6,], Adiv.abiotic2[31,], DIp.or[7:40,])
rownames(DIp.or) <- NULL
DIp.or <- DIp.or[,-c(5,6,9:38)]
DIp.or <- cbind(DIp.or, phytoInt.or)

SpR <- apply(DIp.or[7:84], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
SpA <- apply(DIp.or[7:84], 1, function(x) sum(x))
SD <- apply(DIp.or[7:84], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(DIp.or[7:84], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness 
DiversityI.or <- data.frame(Richness=SpR, Cellcount=SpA, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
IntDiverse.or <- cbind(DIp.or[1:6], DiversityI.or, S, Theta, sigPoDen)
colnames(phytoInt.or) <- c(as.list(colnames(origin[3:81])))
as.factor(colnames(phytoInt.or))
