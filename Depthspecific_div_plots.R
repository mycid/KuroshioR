Upperdepth <- Adiv.abiotic2[which(Adiv.abiotic2[, 7]<40),]
#Making water mass groups
Edge <- Adiv.abiotic2[which(Adiv.abiotic2[, 17 ]>=34.1),]
Edge <- Edge[which(Edge[,17]<34.3),]
Edge <- cbind("mass"=rep("E", nrow(Edge)), Edge)
Kuroshio <- Adiv.abiotic2[which(Adiv.abiotic2[, 17 ]>=34.0),]
Kuroshio <- cbind("mass"=rep("K", nrow(Kuroshio)), Kuroshio)
Oyashio<- Adiv.abiotic2[which(Adiv.abiotic2[, 17] <= 33.55),]
Oyashio <- cbind("mass"=rep("O", nrow(Oyashio)), Oyashio)
Mixed <- Adiv.abiotic2[which(Adiv.abiotic2[,17]>33.55),]
Mixed <- Mixed[which(Mixed[,17]<34.1),]
Mixed <- cbind("mass"=rep("M", nrow(Mixed)), Mixed)
Adiv.abiotic3<-rbind(Kuroshio, Edge, Mixed,Oyashio)
Adiv.abiotic3$mass<- as.factor(Adiv.abiotic3$mass)
Upperdepth <- Adiv.abiotic3[which(Adiv.abiotic3[, 8]<50),]
#Area of apparent activity, ~34.2S, test set up 
Surface <- Adiv.abiotic2[which(Adiv.abiotic2[, 7]<1),]
Edge <- Adiv.abiotic2[which(Adiv.abiotic2[,17] >=33.5),]
Edge <- Edge[which(Edge[,17]<=34.0),]
Edge <- cbind("action"=rep("Active", nrow(Edge)), Edge)
Oth <- Adiv.abiotic2[which(Adiv.abiotic2[,17]>34.0),]
Othe <- Adiv.abiotic2[which(Adiv.abiotic2[,17]<33.5),]
Other <- rbind(Oth, Othe)
Other <- cbind("action"=rep("Other", nrow(Other)), Other)
Adiv.abiotic4 <- rbind(Edge, Other)
Adiv.abiotic4$action <- as.factor(Adiv.abiotic4$action)
Other1 <- Other[which(Other[, 8]<45),]
Edge1 <- Edge[which(Edge[,8]< 45),]
#
#Stats tests for significance
test1<-aov(Richness~mass, data=Upperdepth)
summary(test1) #only significant at top 40m of the water column
TukeyHSD(test1)
Ttest1 <- rstandard(test1)
plot(resid(test1))
qqnorm(test1)
qqline(test1)
#
#Action site versus everywhere else
t.test(Edge$Richness,Other$Richness)# significant and current S/Theta seems optimal within .05
var.test(Edge$Richness,Other$Richness)# significant and same as above for richness, not significant for Shannon Wiener
#
#this test for upper 40m so as to consider depth
t.test(Edge1$ShannonWiener,Other1$ShannonWiener)#pvalue goes down compared to no depth filter test but still highly significant
var.test(Edge1$Richness,Other1$Richness)# pvalue goes down, somewhat significant



#Interpilating the Potential desnity gradient
distance <- -distance+max(distance)
distance3 <- -distance3+max(distance3)
distance5 <- -distance5+max(distance5)
km <- c(distance, distance1, distance3, distance4, distance5)
length(km)
SalGradientA <- with(S.A, (S[-1] - S[-length(S)])/(distance[-1]-distance[-length(distance)]))
AdjDistA <- distance[-7]+((distance[-1] - distance[-length(distance)])/2)
Front4A <- as.data.frame(cbind(SalGradientA,AdjDistA))


interp(A$)
mean(Adiv.abiotic2$Richness)
plot(Adiv.abiotic2$Richness)



filled.contour(interp(A$latitude, A$depth, A$sigPoDen, duplicate="mean"), xlab="latitude", ylab="depth", color=function(x)rev(rainbow(x)), main="Transect A", key.title=title("Potential Density"))
filled.contour(interp(Middle$S, Middle$sigma_t, Middle$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main="Middle Simpson CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Low$sigPoDen, Low$S, Low$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="SigmaT", color=function(x)rev(rainbow(x)), main="Low Eveness CI Plot", key.title=title("     Shannon #"))
