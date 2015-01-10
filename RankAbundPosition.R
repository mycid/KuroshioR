#Abundance vs. Richness
SpA <- apply(test2[, 3:74], 1, sum) 
Adiv.abiotic <- cbind(SpA, div.abio)
Adiv.abiotic2 <- Adiv.abiotic[-c(186:190)]
ggplot(Adiv.abiotic2, aes(x = SpA , y = Richness)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Richness vs. Abundance")                              
ggplot(Adiv.abiotic2, aes(x = SpA , y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Richness vs. Evenness")                              
ggplot(Adiv.abiotic2, aes(x = Richness , y = Evenness.SW)) + geom_point(size=5, alpha=.6, label= Adiv.abiotic2$SpA)+ggtitle("Kuroshio Richness vs. Evenness") 
#Binned by Salinity

Surface <- div.abio2[which(div.abio2[, 3]<1),] 
S.A <-Surface[1:7,]
S.B <- Surface[8:15,]
S.C<- Surface[16:24,]
S.D <-Surface[25:32,]
S.E <- Surface[33:39,]
#A
Adist<- cbind(S.A$longitude, S.A$latitude)
dist1 <- distm(Adist)
distance <- dist[1,]*.001
S.A <- cbind(distance,S.A)
plot(S.A$distance, S.A$S)
#B
Bdist<- cbind(S.B$longitude, S.B$latitude)
dist2 <- distm(Bdist)
distance1 <- dist2[1,]*.001
S.B <- cbind(distance1,S.B)
plot(S.B$distance1, S.B$S)
#C
Cdist<- cbind(S.C$longitude, S.C$latitude)
dist3 <- distm(Cdist)
distance2 <- dist3[1,]*.001
S.C <- cbind(distance2,S.C)
plot(S.C$distance2, S.C$S)
#D
Ddist<- cbind(S.D$longitude, S.D$latitude)
dist4 <- distm(Ddist)
distance3 <- dist4[1,]*.001
S.D <- cbind(distance3,S.D)
plot(S.D$distance3, S.D$S)
#E
Edist<- cbind(S.E$longitude, S.E$latitude)
dist5 <- distm(Edist)
distance4 <- dist5[1,]*.001
S.E <- cbind(distance4,S.E)
plot(S.E$distance4, S.E$sigma_t)

