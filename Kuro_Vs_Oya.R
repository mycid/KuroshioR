KuroshioWaters <- Adiv.abiotic2[which(div.abiotic2$S >=34.1),]
KuroshioWaters <- KuroshioWaters[which(KuroshioWaters$S <=34.4),]
KuroshioWaters <- KuroshioWaters[which(KuroshioWaters$Theta >=21),]
OyashioWaters <- Adiv.abiotic2[which(Adiv.abiotic2$S <=33.7),]
OyashioWaters <- OyashioWaters[which(OyashioWaters$Theta<21),]
#Kuroshio
ShEK <- ggplot(KuroshioWaters, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= KuroshioWaters$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Evenness SW All")
SimEK <- ggplot(KuroshioWaters, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= KuroshioWaters$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Simpson Evenness All")
REK <- ggplot(KuroshioWaters, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= KuroshioWaters$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Richness All")
SimK <- ggplot(KuroshioWaters, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= KuroshioWaters$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Simpson All")
grid.arrange(ShEK, SimEK, SimK, REK, ncol=2)
mean(OyashioWaters$Richness)

#Oyashio
ShEO <- ggplot(OyashioWaters, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= OyashioWaters$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Oyashio Evenness SW All")
SimEO <- ggplot(OyashioWaters, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label=OyashioWaters$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Oyashio Simpson Evenness All")
REO <- ggplot(OyashioWaters, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= OyashioWaters$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Oyashio Richness All")
SimO <- ggplot(OyashioWaters, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= OyashioWaters$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Oyashio Simpson All")
grid.arrange(ShEO, SimEO, SimEO, REO, ncol=2)
#models and tables
Kuro <- KuroshioWaters[,c(6:9,13:15)]
Oya <- OyashioWaters[,c(6:9,13:15)]
Mix <- MixedWaters[,c(6:9,13:15)]
Kplot<- 
Oplot<-
Mplot<-

modelSWK <- lm(ShannonWiener~Theta+S, data=KuroshioWaters)#Not significant with any single or paired independent variables
modelRichK <- lm(Richness~S+Theta, data=KuroshioWaters) #This value is significant at Theta but more significant with S
modelESWK <- lm(Evenness.SW~Theta, data=KuroshioWaters)#This value is significant at Theta but more significant with S
modelESK <- lm(Evenness.Sim~Theta, data=KuroshioWaters)#Is significant at Theta but much more significant when including S
modelSimK <- lm(Simpson~Theta, data=KuroshioWaters)#only signficant as Theta
#Mixed
MixedWater <- div.abiotic2[which(div.abiotic2$S >=33.81),]
MixedWaters <- MixedWater[which(MixedWater$S <34.2),]
ShEM <- ggplot(MixedWaters, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=6, alpha=.6, label= MixedWaters$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Mixed Evenness SW All")
SimEM <- ggplot(MixedWaters, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=6, alpha=.6, label= MixedWaters$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Mixed Simpson Evenness All")
REM <- ggplot(MixedWaters, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=6, alpha=.6, label= MixedWaters$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Mixed Richness All")
SimM <- ggplot(MixedWaters, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=6, alpha=.6, label= MixedWaters$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Mixed Simpson All")
grid.arrange(ShEM, SimEM, SimEM, REM, ncol=2)

