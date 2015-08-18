 #Order R files should be executed. Some files depend on others to run and so there order is particular while others do not. Still this should be the go to.
source("RDiversityIndicesCalc.R")
source("SectionCompCode.R")
source("PlotSalandSimpson.R")
source("PotentialDn.R")
#source("PoDenPlots.R") #had to debug
#source("DinovsDiatomgraphs.R")
source("ColorPlot.R")
source("ArcticDiversity.R")
#
SpA <- apply(test[, 3:73], 1, function(x) sum(x))
bigtest <- read.csv("Allbigstuff.csv")
#TotA <-bigtest$Allbig
Transect <-(read.csv("Untitled3.csv"))
#STstat <- ((div.abio$S+div.abio$Theta)*(div.abio$S*div.abio$Theta)) #stats method suggested by Muareen, stats professor, not sure if correctly applied. 
Transect <- rbind("A", Transect)
Adiv.abiotic2 <- cbind(Transect, SpA, div.abiotic )
#
source("DinoDiatomsDiversity.R")
source("Kuro_Vs_Oya.R")
#
Surface<- Adiv.abiotic2[which(Adiv.abiotic2[, 7]<1),] 
S.A <-Surface[1:7,]
S.B <- Surface[8:15,]
S.C<- Surface[16:24,]
S.D <-Surface[25:32,]
S.E <- Surface[33:40,]
#
A<- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="A"),]
Adist<- cbind(A$longitude, A$latitude)
dist1<- distm(Adist)
distance <- dist1[1,]*-.001
ADfromF<-distance+89.70125
B <- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="B"),]
Bdist<- cbind(B$longitude, B$latitude)
bdist2<- distm(Bdist)
distance1 <- bdist2[1,]*.001
BDfromF <-distance1
C<- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="C"),]
Cdist<- cbind(C$longitude, C$latitude)
dist3 <- distm(Cdist)
distance2 <- dist3[1,]*-.001
CDfromF <-distance2+79.3708
D <-Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="D"),]
Ddist<- cbind(D$longitude, D$latitude)
dist4 <- distm(Ddist)
distance3 <- dist4[1,]*.001
DDfromF <-distance3-23.065398
E<- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="E"),]
Edist<- cbind(E$longitude, E$latitude)
dist5 <- distm(Edist)
distance4 <- dist5[1,]*-.001
EDfromF <-distance4
DfromF <- c(ADfromF, BDfromF, CDfromF, DDfromF,EDfromF)
Adiv.abiotic2<-cbind(DfromF,Adiv.abiotic2)
#
source("loopdepthInteg.R")
source("Euclidean.R")