#Abundance vs. Richness
SpA <- apply(test2[, 3:74], 1, sum)
Transect <- (read.csv("Untitled3.csv"))
Transect <- rbind("A", Transect)
Adiv.abiotic2 <- cbind(Transect, div.abio, SpA)
#Scatterplot iterations
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=S)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Salinity")+stat_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x,2), size=1)
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=sigPoDen)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=ShannonE)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradient()+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=NO3..uM.)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener)+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = Evenness.SW, y = Richness, colour=sigPoDen)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = SpA , y = ShannonWiener)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Shannon vs. Abundance")
ggplot(Adiv.abiotic2, aes(x = SpA , y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Abundance vs. Evenness")                              
ggplot(Adiv.abiotic2, aes(x =Richness,y=Evenness.SW, colour=depth)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$Richness)+scale_color_gradientn(colours=jet.colors(13), na.value="black", space="rgb", guide="colourbar")+ggtitle("Kuroshio Richness vs. Evenness") 


Surface <- Adiv.abiotic2[which(Adiv.abiotic2[, 4]<1),] 
Surface[2] <- NULL
S.A <-Surface[1:7,]
S.B <- Surface[8:15,]
S.C<- Surface[16:24,]
S.D <-Surface[25:32,]
S.E <- Surface[33:39,]



Unimodal <- ggplot(Adiv.abiotic2, aes(x = SpA , y = Richness, color=A)) + geom_point(size=5, alpha=.4, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Richness vs. Abundance")
SurfaceAm+stat_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x,2), size=1)
ggplot(Medium, aes(x = longitude, y = latitude, colour=S)) + geom_point(size=5, alpha=.6, label= "Salinity") + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Salinity at the 26-50m")

#A
xw
#Distance from the front calculated from surface front position across who sample then added back to Adiv.abiotic2
AdivA <-Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="A"),]
AdivB <- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="B"),]
AdivC<- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="C"),]
AdivD <-Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="D"),]
AdivE <- Adiv.abiotic2[which(Adiv.abiotic2[, 1]=="E"),]
distA <- cbind(AdivA$longitude, AdivA$latitude)
distA <- distm(distA)
ADfromF <- (c(distA[1,])*-.001)+80.228
distB <- cbind(AdivB$longitude, AdivB$latitude)
distB <- distm(distB)
BDfromF <- (distB[1,]*.001)-14.32125
distC <- cbind(AdivC$longitude, AdivC$latitude)
distC <- distm(distC)
CDfromF <- (c(distC[1,])*-.001)+79.3708
distD <- cbind(AdivD$longitude, AdivD$latitude)
distD <- distm(distD)
DDfromF <- (c(distD[1,])*.001)-23.065398
distE <- cbind(AdivE$longitude, AdivE$latitude)
distE <- distm(distE)
EDfromF <- (c(distE[1,])*-.001)+32.471895
DfromF <- c(ADfromF,BDfromF,CDfromF,DDfromF,EDfromF)
Adiv.abiotic <-cbind(DfromF, Adiv.abiotic2)
#
AllDist1 <- rbind(S.A , S.B , S.C , S.D , S.E)
Surf <- ggplot(AllDist1, aes(x = DfromF, y =S, color=Theta)) + geom_point(size=5, alpha=.4, label= AllDist1$Richness)+scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Richness and Distance from the Front")
Surf+stat_smooth(aes(group = 1), method = "lm", size=1)
RichAbund<-(lm(Richness~DfromF, data=Adiv.abiotic))
plot(lm(AllDist$DfromF~AllDist$Richness))
 