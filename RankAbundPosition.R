#Abundance vs. Richness
#Scatterplot iterations
ggplot(Adiv.abiotic2, aes(x = SpA , y = Richness, color=A)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Richness vs. Abundance")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=S)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=sigPoDen)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Cell Count, Density, and Richness")
ggplot(Adiv.abiotic2, aes(x = Richness, y = SpA, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Richness, Abundance, and Shannon Wiener")
ggplot(Adiv.abiotic2, aes(x = SpA, y = Richness, colour=ShannonE)) + geom_point(size=6, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Cell Count, Density, and Richness")
ggplot(Adiv.abiotic2, aes(x = Richness, y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), limits=c(15,26), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = Evenness.SW, y = Richness, colour=sigPoDen)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Abundance, Richness, and Theta")
ggplot(Adiv.abiotic2, aes(x = PoDen, y = SpA, colour=Richness)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Cell Count, Density, and Richness")

ggplot(Adiv.abiotic2, aes(x = Richness , y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Shannon vs. Abundance")
ggplot(Adiv.abiotic2, aes(x = SpA , y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Abundance vs. Evenness")                              
ggplot(Adiv.abiotic2, aes(x =Richness,y=Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$Richness)+ggtitle("Kuroshio Richness vs. Evenness") 



S.A <-Surface[1:7,]
S.B <- Surface[8:15,]
S.C<- Surface[16:24,]
S.D <-Surface[25:32,]
S.E <- Surface[33:40,]


#A
Adist<- cbind(S.A$lon, S.A$lat)
dist1 <- distm(Adist)
distance <- dist1[1,]*.001
SalGradientA <- with(S.A, (S[-1] - S[-length(S)])/(distance[-1]-distance[-length(distance)]))
AdjDistA <- distance[-7]+((distance[-1] - distance[-length(distance)])/2)
Front4A <- as.data.frame(cbind(SalGradientA,AdjDistA))
#ggplot(Front4A, aes(x =AdjDistA , y =SalGradientA)) + geom_point(size=5, alpha=.6, label=Front4A$SalGradientA )+ggtitle("Salinity Gradient over Distance T.A")
S.A <- cbind(distance, S.A)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistA, SalGradientA, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect A change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,135), ylim=c(-.15,.15))
par(new=TRUE)
plot(S.A$distance, S.A$S, axes=FALSE, xlab="", ylab="", ylim=c(33.6, 34.3), xlim=c(0,135), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.6,34.3), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist1[1,])*-.001)+89.70125
S.A <- S.A[-1]
S.A <- cbind(DfromF, S.A)
#B
Bdist<- cbind(S.B$lon, S.B$lat)
dist2 <- distm(Bdist)
distance1 <- dist2[1,]*.001
SalGradientB <- with(S.B, (S[-1] - S[-length(S)])/(distance1[-1]-distance1[-length(distance1)]))
AdjDistB <- distance1[-8]+((distance1[-1] - distance1[-length(distance1)])/2)
Front4B <- as.data.frame(cbind(SalGradientB, AdjDistB))
#ggplot(Front4B, aes(x =AdjDistB , y =SalGradientB)) + geom_point(size=5, alpha=.6, label=Front4B$SalGradientB )+ggtitle("Salinity Gradient over Distance T.B")
S.B <- cbind(distance1, S.B)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistB, SalGradientB, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect B change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.10,.10))
par(new=TRUE)
plot(S.B$distance1, S.B$S, axes=FALSE, xlab="", ylab="", ylim=c(33.1, 34.2), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.2), col="blue" ,col.axis="blue",las=1)
legend("topleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (dist2[5,]*.001)
S.B <- S.B[-1]
S.B <- cbind(DfromF, S.B)

#C
Cdist<- cbind(S.C$lon, S.C$lat)
dist3 <- distm(Cdist)
distance3 <- dist3[1,]*.001
SalGradientC <- with(S.C, (S[-1] - S[-length(S)])/(distance3[-1]-distance3[-length(distance3)]))
AdjDistC <- distance3[-9]+((distance3[-1] - distance3[-length(distance3)])/2)
Front4C <- as.data.frame(cbind(SalGradientC,AdjDistC))
#ggplot(Front4C, aes(x =AdjDistC , y =SalGradientC)) + geom_point(size=5, alpha=.6, label=Front4C$SalGradientC )+ggtitle("Salinity Gradient over Distance T.C")
S.C <- cbind(distance3, S.C)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistC, SalGradientC, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect C change in S Gradient/Distance", col.main=51, 
sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,90), ylim=c(-.12,.12))
par(new=TRUE)
plot(S.C$distance3, S.C$S, axes=FALSE, xlab="", ylab="", ylim=c(32.6, 34.5), xlim=c(0,90), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.5), col="blue" ,col.axis="blue",las=1)
legend("topleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.6,
text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist3[1,])*-.001)+79.3708
S.C <- S.C[-1]
S.C <- cbind(DfromF, S.C)
#
#D
Ddist<- cbind(S.D$lon, S.D$lat)
dist4 <- distm(Ddist)
distance4 <- dist4[1,]*.001
SalGradientD <- with(S.D, (S[-1] - S[-length(S)])/(distance4[-1]-distance4[-length(distance4)]))
AdjDistD <- distance4[-8]+((distance4[-1] - distance4[-length(distance4)])/2)
Front4D <- as.data.frame(cbind(SalGradientD,AdjDistD))
#ggplot(Front4D, aes(x =AdjDistD , y =SalGradientD)) + geom_point(size=5, alpha=.6, label=Front4D$SalGradientD )+ggtitle("Salinity Gradient over Distance T.D")
S.D <- cbind(distance4, S.D)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistD, SalGradientD, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect D change in S Gradient/Distance", col.main=51, 
sub="Plot Overlayed by Initial salinity and distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.12,.12))
par(new=TRUE)
plot(S.D$distance4, S.D$S, axes=FALSE, xlab="", ylab="", ylim=c(32.75, 34.5), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(32.75,34.5), col="blue" ,col.axis="blue",las=1)
legend("topright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.9,
text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist4[3,])*.001)-23.065398
S.D <- S.D[-1]
S.D <- cbind(DfromF, S.D)
#

#E
Edist<- cbind(S.E$lon, S.E$lat)
dist5 <- distm(Edist)
distance5 <- dist5[1,]*.001
SalGradientE <- with(S.E, (S[-1] - S[-length(S)])/(distance5[-1]-distance5[-length(distance5)]))
AdjDistE <- distance5[-7]+((distance5[-1] - distance5[-length(distance5)])/2)
Front4E <- as.data.frame(cbind(SalGradientE, AdjDistE))
#ggplot(Front4E, aes(x =AdjDistE , y =SalGradientE)) + geom_point(size=5, alpha=.6, label=Front4E$SalGradientE )+ggtitle("Salinity Gradient over Distance T.E")
S.E <- cbind(distance5, S.E)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistE, SalGradientE, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect E change in S Gradient/Distance", col.main=51, 
sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,55), ylim=c(-.12,.12))
par(new=TRUE)
plot(S.E$distance5, S.E$S, axes=FALSE, xlab="", ylab="", ylim=c(33.4, 34.5), xlim=c(0,55), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33,35), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
text.col=c("green","blue"),pch=c(23, 24), col=c("green","blue"))
DfromF <- (c(dist5[1,])*-.001)
S.E <- S.E[-1]
S.E <- cbind(DfromF, S.E)
#
AllDist <- rbind(S.A, S.B, S.C, S.D, S.E)
#Richness Gradient
RichGradientE <- with(S.E, (Richness[-1] - Richness[-length(Richness)])/(distance5[-1]-distance5[-length(distance5)]))
RichGradientD <- with(S.D, (Richness[-1] - Richness[-length(Richness)])/(distance4[-1]-distance4[-length(distance4)]))
RichGradientC <- with(S.C, (Richness[-1] - Richness[-length(Richness)])/(distance3[-1]-distance3[-length(distance3)]))
RichGradientB <- with(S.B, (Richness[-1] - Richness[-length(Richness)])/(distance1[-1]-distance1[-length(distance1)]))
RichGradientA <- with(S.A, (Richness[-1] - Richness[-length(Richness)])/(distance[-1]-distance[-length(distance)]))
RichGradient<- c(RichGradientA, RichGradientB, RichGradientC, RichGradientD, RichGradientE)
SalGradient<- c(SalGradientA, SalGradientB, SalGradientC, SalGradientD, SalGradientE)
Transect <- c("A",  "A",  "A","A","A",  "A",  "A",  "B",  "B", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "C", "D", "D", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E", "E", "E")
AllDist <- cbind(Transect, AllDist)
ggplot(AllDist, aes(x = DfromF , y = Richness, color=Transect)) + geom_point(size=5, alpha=.6, label= AllDist$ShannonWiener)+ggtitle("Richness and Distance from the Front")
summary(lm(sigPoDen~SpA, data=Adiv.abiotic2))

#Surface DistFrom Front imposed over whole data set!

A<- Adiv.abiotic[which(Adiv.abiotic[,"lon"]==143.5),]
Adist<- cbind(A$lon, A$lat)
dist1<- distm(Adist)
distance <- dist1[1,]*-.001
ADfromF<-distance+89.70125
B <- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==144.0),]
Bdist<- cbind(B$lon, B$lat)
bdist2<- distm(Bdist)
distance1 <- bdist2[1,]*.001
BDfromF <-distance1
C<- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==144.5),]
Cdist<- cbind(C$lon, C$lat)
dist3 <- distm(Cdist)
distance2 <- dist3[1,]*-.001
CDfromF <-distance2+79.3708
D <-Adiv.abiotic[which(Adiv.abiotic[, "lon"]==145.0),]
Ddist<- cbind(D$lon, D$lat)
dist4 <- distm(Ddist)
distance3 <- dist4[1,]*.001
DDfromF <-distance3-23.065398
E<- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==145.5),]
Edist<- cbind(E$lon, E$lat)
dist5 <- distm(Edist)
distance4 <- dist5[1,]*-.001
EDfromF <-distance4
DfromF <- c(ADfromF, BDfromF, CDfromF, DDfromF,EDfromF)
Adiv.abiotic<-cbind(Adiv.abiotic, DfromF)
#Making graphs for this stuff
Adiv.abiotic2$A<- as.factor(Adiv.abiotic2$A)
Ad <- ggplot(Adiv.abiotic2, aes(x = DfromF, y = SpA, colour=Richness)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$Richness) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Cell Count Abundance/Distance from the Front")
Ad+stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)+guides(colour = "colorbar", barwidth = 0.5, nbin = 15, barheight = 10)
summary(lm(DfromF~Gccom, data=Adiv.abiotic))
#guide_colourbar(title = waiver(), title.position = NULL, title.theme = NULL, title.hjust = NULL, title.vjust = NULL, label = TRUE, label.position = NULL, label.theme = NULL, label.hjust = NULL, label.vjust = NULL, barwidth = NULL, barheight = 3, nbin = 20, raster = TRUE, ticks = TRUE, draw.ulim = TRUE, draw.llim = TRUE, direction = NULL, default.unit = "line", reverse = FALSE, order = 0)

