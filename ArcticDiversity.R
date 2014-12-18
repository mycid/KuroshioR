
Arc <- read.csv("Arctic Ocean data sheet.csv")
Arc[,117] <-NULL
Arc[,3] <-NULL
Arc[,51] <-NULL
ArcSpR <- apply(Arc[, 3:114], 1, function(x) sum(x>0)) #species richness
ArcSD <- apply(Arc[, 3:114], 1, function(x) (sum(x*(x-1))/(sum(x)*(sum(x)-1))))
ArcSimE <- (1/ArcSD)/ArcSpR
ArcSW <- apply(Arc[, 3:114], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
ArcSWD <- colSums (ArcSW, na.rm=T) #Shannon Wiener Diversity Index 
ArcShannonE <- ArcSWD/log(ArcSpR) #Eveness 
ArcDiversityI <- data.frame(Richness=ArcSpR, ShannonWiener=ArcSWD, Simpson=ArcSD, Evenness.SW=ArcShannonE, Evenness.Sim=ArcSimE)
ArcD<- cbind(Arc[,1:2],ArcDiversityI)
#Diatoms
DASpR <- apply(Arc[, 3:50], 1, function(x) sum(x>0)) #species richness
DASD <- apply(Arc[, 3:50], 1, function(x) (sum(x*(x-1))/(sum(x)*(sum(x)-1))))
DASimE <- (1/ArcSD)/DASpR
DASW <- apply(Arc[, 3:50], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
DASWD <- colSums (DASW, na.rm=T) #Shannon Wiener Diversity Index 
DAShannonE <- ArcSWD/log(DASpR) #Eveness 
DADiversityI <- data.frame(Richness=DASpR, ShannonWiener=DASWD, Simpson=DASD, Evenness.SW=DAShannonE, Evenness.Sim=DASimE)
DAD <- cbind(Arc[,1:2],DADiversityI)
#Dino
dinoArcSpR <- apply(Arc[, 51:114], 1, function(x) sum(x>0)) #species richness
dinoArcSD <- apply(Arc[, 51:114], 1, function(x) (sum(x*(x-1))/(sum(x)*(sum(x)-1))))
dinoArcSimE <- (1/dinoArcSD)/dinoArcSpR
dinoArcSW <- apply(Arc[, 51:114], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
dinoArcSWD <- colSums (dinoArcSW, na.rm=T) #Shannon Wiener Diversity Index 
dinoArcShannonE <- dinoArcSWD/log(dinoArcSpR) #Eveness 
dinoArcDiversityI <- data.frame(Richness=dinoArcSpR, ShannonWiener=dinoArcSWD, Simpson=dinoArcSD, Evenness.SW=dinoArcShannonE, Evenness.Sim=dinoArcSimE)
dinoArcD <- cbind(Arc[,1:2],dinoArcDiversityI)
Group <- c("Dino", "Dino","Dino", "Dino", "Dino","Dino", "Dino", "Dino","Dino", "Diatom", "Diatom","Diatom","Diatom", "Diatom","Diatom", "Diatom", "Diatom","Diatom")
ArcCompare <- rbind( dinoArcD, DAD)
ArcCompare <- cbind(Group, ArcCompare)
