#Only diatoms and dinoflagellates for this one!
#Both dinos and diatoms
bigtest <- read.csv("Allbigstuff.csv") #Microscopy counts
KC <- read.csv("Kuro_Phytoplankton_coords.csv")
SpR <- apply(bigtest[, c(3:9)], 1, function(x) sum(x>0)) #species richness
SD <- apply(bigtest[, c(3:9)], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(bigtest[, c(3:9)], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness
GroupsBioM <- apply(bigtest[,c(3:9)], 1, function(x) sum(x))
DiversityI <- data.frame(Richness=SpR, GroupsBiomass=GroupsBioM, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
DiversityIndex <- cbind(KC, DiversityI) 
Al<- read.csv("KuroAlldata.csv")   
All <-Al[complete.cases(Al$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton
the <- All$depth..m.
t <- All$T.C.
Theta <- t-((.1*the)/1000)
T <- Theta
S <- All$S
p0<- 999.842594+6.793952*10^(-2)*T-9.095290*10^(-3)*T^(2)+1.001685*10^(-4)*T^(3)-1.120083*10^(-6)*T^(4)+6.536332*10^(-9)*T^(5)
A <- 8.24493*10^(-1)-4.0899*10^(-3)*T+7.6438*10^(-5)*T^(2)-8.2467*10^(-7)*T^(3)+5.3875*10^(-9)*T^(4)
B <- (-5.72466*10^(-3))+1.0227*10^(-4)*T-1.6546*10^(-6)*T^(2)
C <- 4.8314*10^(-4)
PoDen <- p0+A*S+B*S^(1.5)+C*S^(2)
sigPoDen <- PoDen-1000
Bdiv.abio <- cbind(Adiv.abiotic2[,1:4], DiversityIndex, Theta, All[, 7:9], sigPoDen, All[, 11:30])
Bdiv.abio$theta=NULL
do.call(data.frame,lapply(Bdiv.abio, function(x) replace(x, is.infinite(x),1000)))
View(Bdiv.abio)
#Only dinos

DSpR <- apply(test[, 3:24], 1, function(x) sum(x>0)) #species richness
DSD <- apply(test[, 3:24], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
DSimE <- (1/DSD)/DSpR
DSW <- apply(test[, 3:24], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
DSWD <- colSums (DSW, na.rm=T) #Shannon Wiener Diversity Index 
DShannonE <- DSWD/log(DSpR) #Eveness 
DDiversityI <- data.frame(Richness=DSpR, ShannonWiener=DSWD, Simpson=DSD, Evenness.SW=DShannonE, Evenness.Sim=DSimE)
DinoDiversityIndex <- cbind(KC, DDiversityI) 
Al<- read.csv("KuroAlldata.csv")   
All <-Al[complete.cases(Al$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton
the <- All$depth..m.
t <- All$T.C.
Theta <- t-((.1*the)/1000)
T <- Theta
S <- All$S
p0<- 999.842594+6.793952*10^(-2)*T-9.095290*10^(-3)*T^(2)+1.001685*10^(-4)*T^(3)-1.120083*10^(-6)*T^(4)+6.536332*10^(-9)*T^(5)
A <- 8.24493*10^(-1)-4.0899*10^(-3)*T+7.6438*10^(-5)*T^(2)-8.2467*10^(-7)*T^(3)+5.3875*10^(-9)*T^(4)
B <- (-5.72466*10^(-3))+1.0227*10^(-4)*T-1.6546*10^(-6)*T^(2)
C <- 4.8314*10^(-4)
PoDen <- p0+A*S+B*S^(1.5)+C*S^(2)
sigPoDen <- PoDen-1000
dinodiv.abio <- cbind(DinoDiversityIndex, Theta, All[, 7:9], sigPoDen, All[, 11:30])
dinodiv.abio$theta=NULL
View(dinodiv.abio)
dinodiv.abio2 <- dinodiv.abio[-c(186:190),]
#Only diatoms
DiSpR <- apply(test[, 25:73], 1, function(x) sum(x>0)) #species richness
DiSD <- apply(test[, 25:73], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
DiSimE <- (1/DiSD)/DiSpR
DiSW <- apply(test[, 25:73], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
DiSWD <- colSums (DiSW, na.rm=T) #Shannon Wiener Diversity Index 
DiShannonE <- DiSWD/log(DiSpR) #Eveness 
DiDiversityI <- data.frame(Richness=DiSpR, ShannonWiener=DiSWD, Simpson=DiSD, Evenness.SW=DiShannonE, Evenness.Sim=DiSimE)
DiatomDiversityIndex <- cbind(KC, DiDiversityI) 
Al<- read.csv("KuroAlldata.csv")   
All <-Al[complete.cases(Al$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton
the <- All$depth..m.
t <- All$T.C.
Theta <- t-((.1*the)/1000)
T <- Theta
S <- All$S
p0<- 999.842594+6.793952*10^(-2)*T-9.095290*10^(-3)*T^(2)+1.001685*10^(-4)*T^(3)-1.120083*10^(-6)*T^(4)+6.536332*10^(-9)*T^(5)
A <- 8.24493*10^(-1)-4.0899*10^(-3)*T+7.6438*10^(-5)*T^(2)-8.2467*10^(-7)*T^(3)+5.3875*10^(-9)*T^(4)
B <- (-5.72466*10^(-3))+1.0227*10^(-4)*T-1.6546*10^(-6)*T^(2)
C <- 4.8314*10^(-4)
PoDen <- p0+A*S+B*S^(1.5)+C*S^(2)
sigPoDen <- PoDen-1000
Diatomdiv.abio <- cbind(DiatomDiversityIndex, Theta, All[, 7:9], sigPoDen, All[, 11:30])
Diatomdiv.abio$theta=NULL
View(Diatomdiv.abio)
Diatomdiv.abio2 <- Diatomdiv.abio[-c(186:190),]
DinoPercent <- dinodiv.abio$Richness/Diatomdiv.abio$Richness
div.abio2 <- cbind(DinoPercent, Diatomdiv.abio)
