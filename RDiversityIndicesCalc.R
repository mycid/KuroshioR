test <- read.csv("Kuroshio_Phytoplankton.csv") #Microscopy counts
KC <- read.csv("Kuro_Phytoplankton_coords.csv")
SpR <- apply(test[, 3:73], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
SD <- apply(test[, 3:73], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(test[, 3:73], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness 
DiversityI <- data.frame(Richness=SpR, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
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
div.abiotic <- cbind(DiversityIndex, Theta, All[, 7:9], sigPoDen, All[, 11:30])
div.abiotic$theta=NULL
View(div.abiotic)
div.abiotic2 <- div.abiotic[-c(186:190),]
#all groups 
test2 <- read.csv("OriginalPhytoplankton.csv") #Microscopy counts
KC <- read.csv("Kuro_Phytoplankton_coords.csv")
SpR <- apply(test2[, 3:81], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
SD <- apply(test2[, 3:81], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1)))
SimE <- (1/SD)/SpR
SW <- apply(test2[, 3:81], 1, function(x) (x/sum(x))*(-log(x/sum(x))))
SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
ShannonE <- SWD/log(SpR) #Eveness 
DiversityI <- data.frame(Richness=SpR, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE)
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
Aldiv.abiotic <- cbind(DiversityIndex, Theta, All[, 7:9], sigPoDen, test2$Cellcount, All[, 11:30])
Aldiv.abiotic$Theta=NULL

#div.abiotic2 <- div.abiotic[-c(186:190),]

#dissimilarity Index scripts
d <- betadiver(test[3:73], "w")
range(d - vegdist(test[3:73], binary=TRUE))

