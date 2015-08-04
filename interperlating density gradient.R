#Calculating the linear distance for each transect                         
#
# Interpolating without depth, just surface for salinity and density only
distance <- -distance+max(distance)
distance3 <- -distance3+max(distance3)
distance5 <- -distance5+max(distance5)
km <- c(distance, distance1, distance3, distance4, distance5)
length(km)
#S.A Sig
SigGradientA <- abs(with(S.A, (sigPoDen[-1] - sigPoDen[-length(sigPoDen)])/(distance[-1]-distance[-length(distance)])))
AdjDistA <- distance[-7]+((distance[-1] - distance[-length(distance)])/2)
intPoDenA<- abs(interp1(x=distance, y=S.A$sigPoDen, xi=c(AdjDistA), method="linear"))
plot(SalGradientA, interpPoDen)
#S.B Sig
SigGradientB <- abs(with(S.B, (sigPoDen[-1] - sigPoDen[-length(sigPoDen)])/(distance1[-1]-distance1[-length(distance1)])))
AdjDistB <- distance1[-8]+((distance1[-1] - distance1[-length(distance1)])/2)
intPoDenB<- abs(interp1(x=distance1, y=S.B$sigPoDen, xi=c(AdjDistB), method="linear"))
#S.C Sig
SigGradientC <- abs(with(S.C, (sigPoDen[-1] - sigPoDen[-length(sigPoDen)])/(distance3[-1]-distance3[-length(distance3)])))
AdjDistC <- distance3[-9]+((distance3[-1] - distance3[-length(distance3)])/2)
intPoDenC<- abs(interp1(x=distance3, y=S.C$sigPoDen, xi=c(AdjDistC), method="linear"))
#S.D Sig
SigGradientD <- abs(with(S.D, (sigPoDen[-1] - sigPoDen[-length(sigPoDen)])/(distance4[-1]-distance4[-length(distance4)])))
AdjDistD <- distance4[-8]+((distance4[-1] - distance4[-length(distance4)])/2)
intPoDenD<- abs(interp1(x=distance4, y=S.D$sigPoDen, xi=c(AdjDistD), method="linear"))
#S.E sig
SigGradientE <- abs(with(S.E, (sigPoDen[-1] - sigPoDen[-length(sigPoDen)])/(distance5[-1]-distance5[-length(distance5)])))
AdjDistE <- distance5[-8]+((distance5[-1] - distance5[-length(distance5)])/2)
intPoDenE<- abs(interp1(x=distance5, y=S.E$sigPoDen, xi=c(AdjDistE), method="linear"))
#All surface Sig
SigGradient <- c(SigGradientA, SigGradientB, SigGradientC, SigGradientD,SigGradientE)
intPoDen <- c(intPoDenA, intPoDenB, intPoDenC, intPoDenD,intPoDenE)
plot(intPoDen, SigGradient) 
#
#S.A S
SGradientA <- abs(with(S.A, (S[-1] - S[-length(S)])/(distance[-1]-distance[-length(distance)])))
AdjDistA <- distance[-7]+((distance[-1] - distance[-length(distance)])/2)
intSA<- abs(interp1(x=distance, y=S.A$S, xi=c(AdjDistA), method="linear"))
#S.B S
SGradientB <- abs(with(S.B, (S[-1] - S[-length(S)])/(distance1[-1]-distance1[-length(distance1)])))
AdjDistB <- distance1[-8]+((distance1[-1] - distance1[-length(distance1)])/2)
intSB<- abs(interp1(x=distance1, y=S.B$S, xi=c(AdjDistB), method="linear"))
#S.C S
SGradientC <- abs(with(S.C, (S[-1] - S[-length(S)])/(distance3[-1]-distance3[-length(distance3)])))
AdjDistC <- distance3[-9]+((distance3[-1] - distance3[-length(distance3)])/2)
intSC<- abs(interp1(x=distance3, y=S.C$S, xi=c(AdjDistC), method="linear"))
#S.D S
SGradientD <- abs(with(S.D, (S[-1] - S[-length(S)])/(distance4[-1]-distance4[-length(distance4)])))
AdjDistD <- distance4[-8]+((distance4[-1] - distance4[-length(distance4)])/2)
intSD<- abs(interp1(x=distance4, y=S.D$S, xi=c(AdjDistD), method="linear"))
#S.E S
SGradientE <- abs(with(S.E, (S[-1] - S[-length(S)])/(distance5[-1]-distance5[-length(distance5)])))
AdjDistE <- distance5[-8]+((distance5[-1] - distance5[-length(distance5)])/2)
intSE<- abs(interp1(x=distance5, y=S.E$S, xi=c(AdjDistE), method="linear"))
#All S
SGradient <- c(SGradientA, SGradientB, SGradientC, SGradientD,SGradientE)
intS <- c(intSA, intSB, intSC, intSD,intSE)
plot(intS, SGradient) 
plot(Surface$S, Surface$Richness)
#Finding the Salinity gradient at the surface and at the actual stations
#S.A. S
AdjDistA <- distance[-7]+((distance[-1] - distance[-length(distance)])/2)
intSA<- abs(interp1(x=distance, y=S.A$S, xi=c(AdjDistA), method="linear"))
SGradientA <- abs(with(S.A, (intSA[-1] - intSA[-length(intSA)])/(AdjDistA[-1]-AdjDistA[-length(AdjDistA)])))
#S.B S
AdjDistB <- distance1[-8]+((distance1[-1] - distance1[-length(distance1)])/2)
intSB<- abs(interp1(x=distance1, y=S.B$S, xi=c(AdjDistB), method="linear"))
SGradientB <- abs(with(S.B, (intSB[-1] - intSB[-length(intSB)])/(AdjDistB[-1]-AdjDistB[-length(AdjDistB)])))
#S.C S
AdjDistC <- distance3[-9]+((distance3[-1] - distance3[-length(distance3)])/2)
intSC<- abs(interp1(x=distance3, y=S.C$S, xi=c(AdjDistC), method="linear"))
SGradientC <- abs(with(S.C, (intSC[-1] - intSC[-length(intSC)])/(AdjDistC[-1]-AdjDistC[-length(AdjDistC)])))

#S.D S
AdjDistD <- distance4[-8]+((distance4[-1] - distance4[-length(distance4)])/2)
intSD<- abs(interp1(x=distance4, y=S.D$S, xi=c(AdjDistD), method="linear"))
SGradientD <- abs(with(S.D, (intSD[-1] - intSD[-length(intSD)])/(AdjDistD[-1]-AdjDistD[-length(AdjDistD)])))

#S.E S
AdjDistE <- distance5[-8]+((distance5[-1] - distance5[-length(distance5)])/2)
intSE<- abs(interp1(x=distance5, y=S.E$S, xi=c(AdjDistE), method="linear"))
SGradientE <- abs(with(S.E, (intSE[-1] - intSE[-length(intSE)])/(AdjDistE[-1]-AdjDistE[-length(AdjDistE)])))
SGradient <- c(SGradientA, SGradientB, SGradientC, SGradientD,SGradientE)
IntSalDiv <- IntDiverse[-c(1,7,8,9,16,17,25,26,33,34,41),]
plot(SGradient, IntSalDiv$ShannonWiener)
