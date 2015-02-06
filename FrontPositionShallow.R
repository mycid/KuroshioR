Shallow <- Adiv.abiotic2[which(Adiv.abiotic2[4]>1),] 
Shallow <-Shallow[which(Shallow[4]<=25),]
Sh.A <-Shallow[which(Shallow[, 1]=="A"),]
Sh.B <- Shallow[which(Shallow[, 1]=="B"),]
Sh.C<- Shallow[which(Shallow[, 1]=="C"),]
Sh.D <-Shallow[which(Shallow[, 1]=="D"),]
Sh.E <- Shallow[which(Shallow[, 1]=="E"),]
#A
Adist<- cbind(Sh.A$longitude, Sh.A$latitude)
dist1 <- distm(Adist)
distance <- dist1[1,]*.001
SalGradientA <- with(Sh.A, (S[-1] - S[-length(S)])/(distance[-1]-distance[-length(distance)]))
AdjDistA <- distance[-6]+((distance[-1] - distance[-length(distance)])/2)
Front4A <- as.data.frame(cbind(SalGradientA,AdjDistA))
#ggplot(Front4A, aes(x =AdjDistA , y =SalGradientA)) + geom_point(size=5, alpha=.6, label=Front4A$SalGradientA )+ggtitle("Salinity Gradient over Distance T.A")
Sh.A <- cbind(distance, Sh.A)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistA, SalGradientA, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect A change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,135), ylim=c(-.15,.15))
par(new=TRUE)
plot(Sh.A$distance, Sh.A$S, axes=FALSE, xlab="", ylab="", ylim=c(33.6, 34.3), xlim=c(0,135), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.6,34.3), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist1[1,])*-.001)+80.22796
Sh.A <- Sh.A[-1]
Sh.A <- cbind(DfromF, Sh.A)
#B
Bdist<- cbind(Sh.B$longitude, Sh.B$latitude)
dist2 <- distm(Bdist)
distance1 <- dist2[1,]*.001
SalGradientB <- with(Sh.B, (S[-1] - S[-length(S)])/(distance1[-1]-distance1[-length(distance1)]))
AdjDistB <- distance1[-8]+((distance1[-1] - distance1[-length(distance1)])/2)
Front4B <- as.data.frame(cbind(SalGradientB, AdjDistB))
#ggplot(Front4B, aes(x =AdjDistB , y =SalGradientB)) + geom_point(size=5, alpha=.6, label=Front4B$SalGradientB )+ggtitle("Salinity Gradient over Distance T.B")
Sh.B <- cbind(distance1, Sh.B)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistB, SalGradientB, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect B change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.10,.15))
par(new=TRUE)
plot(Sh.B$distance1, Sh.B$S, axes=FALSE, xlab="", ylab="", ylim=c(33.1, 34.4), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.2), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (dist2[5,]*.001)-23.82237
Sh.B <- Sh.B[-1]
Sh.B <- cbind(DfromF, Sh.B)

#C
Cdist<- cbind(Sh.C$longitude, Sh.C$latitude)
dist3 <- distm(Cdist)
distance3 <- dist3[1,]*.001 
SalGradientC <- with(Sh.C, (S[-1] - S[-length(S)])/(distance3[-1]-distance3[-length(distance3)]))
AdjDistC <- distance3[-8]+((distance3[-1] - distance3[-length(distance3)])/2)
Front4C <- as.data.frame(cbind(SalGradientC, AdjDistC))
#ggplot(Front4C, aes(x =AdjDistC , y =SalGradientC)) + geom_point(size=5, alpha=.6, label=Front4C$SalGradientC )+ggtitle("Salinity Gradient over Distance T.C")
Sh.C <- cbind(distance3, Sh.C)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistC, SalGradientC, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect C change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,90), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sh.C$distance3, Sh.C$S, axes=FALSE, xlab="", ylab="", ylim=c(32.6, 34.4), xlim=c(0,90), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.5), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.6,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist3[1,])*-.001)+28.44770
Sh.C <- Sh.C[-1]
Sh.C <- cbind(DfromF, Sh.C)
#
#D
Ddist<- cbind(Sh.D$longitude, Sh.D$latitude)
dist4 <- distm(Ddist)
distance4 <- dist4[1,]*.001
SalGradientD <- with(Sh.D, (S[-1] - S[-length(S)])/(distance4[-1]-distance4[-length(distance4)]))
AdjDistD <- distance4[-8]+((distance4[-1] - distance4[-length(distance4)])/2)
Front4D <- as.data.frame(cbind(SalGradientD,AdjDistD))
#ggplot(Front4D, aes(x =AdjDistD , y =SalGradientD)) + geom_point(size=5, alpha=.6, label=Front4D$SalGradientD )+ggtitle("Salinity Gradient over Distance T.D")
Sh.D <- cbind(distance4, Sh.D)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistD, SalGradientD, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect D change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sh.D$distance4, Sh.D$S, axes=FALSE, xlab="", ylab="", ylim=c(32.75, 34.5), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(32.75,34.5), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.6,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist4[3,])*.001)-13.692297
Sh.D <- Sh.D[-1]
Sh.D <- cbind(DfromF, Sh.D)
#

#E
Edist<- cbind(Sh.E$longitude, Sh.E$latitude)
dist5 <- distm(Edist)
distance5 <- dist5[1,]*.001
SalGradientE <- with(Sh.E, (S[-1] - S[-length(S)])/(distance5[-1]-distance5[-length(distance5)]))
AdjDistE <- distance5[-8]+((distance5[-1] - distance5[-length(distance5)])/2)
Front4E <- as.data.frame(cbind(SalGradientE,AdjDistE))
#ggplot(Front4E, aes(x =AdjDistE , y =SalGradientE)) + geom_point(size=5, alpha=.6, label=Front4E$SalGradientE )+ggtitle("Salinity Gradient over Distance T.E")
Sh.E <- cbind(distance5, Sh.E)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistE, SalGradientE, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect E change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,55), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sh.E$distance5, Sh.E$S, axes=FALSE, xlab="", ylab="", ylim=c(33.4, 34.5), xlim=c(0,55), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33,35), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23, 24), col=c("green","blue"))
DfromF <- (c(dist5[1,])*-.001)+57.7637
Sh.E <- Sh.E[-1]
Sh.E <- cbind(DfromF, Sh.E)
AllDist <- rbind( Sh.A , Sh.B , Sh.C , Sh.D , Sh.E)
ggplot(AllDist, aes(x = DfromF, y =Richness, color=S)) + geom_point(size=5, alpha=.4, label= AllDist$ShannonWiener)+scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("Richness and Distance from the Front")





Medium <- Adiv.abiotic2[which(Adiv.abiotic2[4]>25),] 
Medium <-Medium[which(Medium[4]<=50),]
Sm.A <-Medium[which(Medium[, 1]=="A"),]
Sm.B <- Medium[which(Medium[, 1]=="B"),]
Sm.C<- Medium[which(Medium[, 1]=="C"),]
Sm.D <-Medium[which(Medium[, 1]=="D"),]
Sm.E <- Medium[which(Medium[, 1]=="E"),]
#A
Adist<- cbind(Sm.A$longitude, Sm.A$latitude)
dist1 <- distm(Adist)
distance <- dist1[1,]*.001
MultiPav <- function(distance=distance){AverageSal <- }
SalGradientA <- with(Sm.A, (S[-1] - S[-length(S)])/(distance[-1]-distance[-length(distance)]))
AdjDistA <- distance[-10]+((distance[-1] - distance[-length(distance)])/2)
Front4A <- as.data.frame(cbind(SalGradientA,AdjDistA))
#ggplot(Front4A, aes(x =AdjDistA , y =SalGradientA)) + geom_point(size=5, alpha=.6, label=Front4A$SalGradientA )+ggtitle("Salinity Gradient over Distance T.A")
Sm.A <- cbind(distance, Sm.A)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistA, SalGradientA, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect A change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,135), ylim=c(-.15,.15))
par(new=TRUE)
plot(Sm.A$distance, Sm.A$S, axes=FALSE, xlab="", ylab="", ylim=c(32, 34.8), xlim=c(0,135), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(32,34.5), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist1[1,])*-.001)+27.74082
Sm.A <- Sm.A[-1]
Sm.A <- cbind(DfromF, Sm.A)
#B
Bdist<- cbind(Sm.B$longitude, Sm.B$latitude)
dist2 <- distm(Bdist)
distance1 <- dist2[1,]*.001
SalGradientB <- with(Sm.B, (S[-1] - S[-length(S)])/(distance1[-1]-distance1[-length(distance1)]))
AdjDistB <- distance1[-13]+((distance1[-1] - distance1[-length(distance1)])/2)
Front4B <- as.data.frame(cbind(SalGradientB, AdjDistB))
#ggplot(Front4B, aes(x =AdjDistB , y =SalGradientB)) + geom_point(size=5, alpha=.6, label=Front4B$SalGradientB )+ggtitle("Salinity Gradient over Distance T.B")
Sm.B <- cbind(distance1, Sm.B)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistB, SalGradientB, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect B change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.10,.15))
par(new=TRUE)
plot(Sm.B$distance1, Sm.B$S, axes=FALSE, xlab="", ylab="", ylim=c(33.1, 34.8), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.8), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (dist2[5,]*.001)-23.82237
Sm.B <- Sm.B[-1]
Sm.B <- cbind(DfromF, Sm.B)

#C
Cdist<- cbind(Sm.C$longitude, Sm.C$latitude)
dist3 <- distm(Cdist)
distance3 <- dist3[1,]*.001
SalGradientC <- with(Sm.C, (S[-1] - S[-length(S)])/(distance3[-1]-distance3[-length(distance3)]))
AdjDistC <- distance3[-17]+((distance3[-1] - distance3[-length(distance3)])/2)
Front4C <- as.data.frame(cbind(SalGradientC,AdjDistC))
#ggplot(Front4C, aes(x =AdjDistC , y =SalGradientC)) + geom_point(size=5, alpha=.6, label=Front4C$SalGradientC )+ggtitle("Salinity Gradient over Distance T.C")
Sm.C <- cbind(distance3, Sm.C)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistC, SalGradientC, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect C change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,90), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sm.C$distance3, Sm.C$S, axes=FALSE, xlab="", ylab="", ylim=c(32.6, 34.4), xlim=c(0,90), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.1,34.5), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.6,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist3[1,])*-.001)+28.447696
Sm.C <- Sm.C[-1]
Sm.C <- cbind(DfromF, Sm.C)
#
#D
Ddist<- cbind(Sm.D$longitude, Sm.D$latitude)
dist4 <- distm(Ddist)
distance4 <- dist4[1,]*.001
SalGradientD <- with(Sm.D, (S[-1] - S[-length(S)])/(distance4[-1]-distance4[-length(distance4)]))
AdjDistD <- distance4[-7]+((distance4[-1] - distance4[-length(distance4)])/2)
Front4D <- as.data.frame(cbind(SalGradientD,AdjDistD))
#ggplot(Front4D, aes(x =AdjDistD , y =SalGradientD)) + geom_point(size=5, alpha=.6, label=Front4D$SalGradientD )+ggtitle("Salinity Gradient over Distance T.D")
Sm.D <- cbind(distance4, Sm.D)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistD, SalGradientD, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect D change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,70), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sm.D$distance4, Sm.D$S, axes=FALSE, xlab="", ylab="", ylim=c(32.75, 34.5), xlim=c(0,70), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(32.75,34.5), col="blue" ,col.axis="blue",las=1)
legend("bottomright",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.6,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))
DfromF <- (c(dist4[3,])*.001)-13.965030
Sm.D <- Sm.D[-1]
Sm.D <- cbind(DfromF, Sm.D)
#

#E
Edist<- cbind(Sm.E$longitude, Sm.E$latitude)
dist5 <- distm(Edist)
distance5 <- dist5[1,]*.001
SalGradientE <- with(Sm.E, (S[-1] - S[-length(S)])/(distance5[-1]-distance5[-length(distance5)]))
AdjDistE <- distance5[-10]+((distance5[-1] - distance5[-length(distance5)])/2)
Front4E <- as.data.frame(cbind(SalGradientE,AdjDistE))
#ggplot(Front4E, aes(x =AdjDistE , y =SalGradientE)) + geom_point(size=5, alpha=.6, label=Front4E$SalGradientE )+ggtitle("Salinity Gradient over Distance T.E")
Sm.E <- cbind(distance5, Sm.E)
par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(AdjDistE, SalGradientE, type="b", pch=23, col="black",bg="green", cex=1.5, main="Transect E change in S Gradient/Distance", col.main=51, 
     sub="Plot Overlayed by Initial S and Distance", col.sub="blue", xlab="Distance (km)", ylab="Salinity Gradient", xlim=c(0,55), ylim=c(-.12,.12))
par(new=TRUE)
plot(Sm.E$distance5, Sm.E$S, axes=FALSE, xlab="", ylab="", ylim=c(33.4, 34.5), xlim=c(0,55), type="b",lty=2, pch=24, col="black", bg=125)
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33,35), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23, 24), col=c("green","blue"))
DfromF <- (c(dist5[1,])*-.001)-57.763684
Sm.E <- Sm.E[-1]
Sm.E <- cbind(DfromF, Sm.E)
AllDist1 <- rbind(AllDist1, Sm.A , Sm.B , Sm.C , Sm.D , Sm.E)
