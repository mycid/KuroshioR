which( colnames(div.abiotic2)=="sigPoDen" )
fstQDen <- div.abiotic2[ div.abio2$sigPoDen <= quantile(div.abio2$sigPoDen ,0.25),]
sndQDen <- div.abiotic2[ div.abio2$sigPoDen > quantile(div.abio2$sigPoDen, 0.25),]
sndQDen <- sndQDen[sndQDen$sigPoDen <= median(sndQDen$sigPoDen),]
trdQDen <- div.abio2[ div.abio2$sigPoDen > median(div.abio2$sigPoDen),]
trdQDen <- trdQDen[trdQDen$sigPoDen <= quantile(trdQDen$sigPoDen, 0.75),]
frthQDen <- div.abio2[ div.abio2$sigPoDen > quantile(div.abio2$sigPoDen, 0.75),]
#Shannon Wiener Plots
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
SW1stQ <- ggplot(fstQDen, aes(x = S , y = Theta, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= fstQDen$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), na.value="black", space="rgb", guide="colourbar")+ggtitle("1st Quartile sigmaT")
SW2ndQ <- ggplot(sndQDen, aes(x = S , y = Theta, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= sndQDen$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("2nd Quartile sigmaT")
SW3rdQ <- ggplot(trdQDen, aes(x = S , y = Theta, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= trdQDen$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("3rd Quartile sigmaT")
SW4thQ <- ggplot(frthQDen, aes(x = S , y = Theta, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= frthQDen$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("4th Quartile sigmaT")
SWall <- ggplot(div.abio2, aes(x = S , y = Theta, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= div.abio2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("SigmaT All")
grid.arrange(SW1stQ, SW2ndQ, SW3rdQ, SW4thQ, SWall, ncol=2)
#Richness Plots
R1stQ <- ggplot(fstQDen, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= fstQDen$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("1st Quartile sigmaT")
R2ndQ <- ggplot(sndQDen, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= sndQDen$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("2nd Quartile sigmaT")
R3rdQ <- ggplot(trdQDen, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= trdQDen$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("3rd Quartile sigmaT")
R4thQ <- ggplot(frthQDen, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= frthQDen$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("4th Quartile sigmaT")
Rall <- ggplot(div.abio2, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= div.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("SigmaT All")
grid.arrange(R1stQ, R2ndQ, R3rdQ, R4thQ, Rall, ncol=2)
#Simpson Eveness Plots
SE1stQ <- ggplot(fstQDen, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= fstQDen$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("1st Quartile sigmaT")
SE2ndQ <- ggplot(sndQDen, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= sndQDen$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("2nd Quartile sigmaT")
SE3rdQ <- ggplot(trdQDen, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= trdQDen$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("3rd Quartile sigmaT")
SE4thQ <- ggplot(frthQDen, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= frthQDen$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("4th Quartile sigmaT")
SEall <- ggplot(div.abio2, aes(x = S , y = Theta, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= div.abio2$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("SigmaT All")
grid.arrange(SE1stQ, SE2ndQ, SE3rdQ, SE4thQ, SEall, ncol=2)
#Shannon Eveness Plots
ShE1stQ <- ggplot(fstQDen, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= fstQDen$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("1st Quartile sigmaT")
ShE2ndQ <- ggplot(sndQDen, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= sndQDen$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("2nd Quartile sigmaT")
ShE3rdQ <- ggplot(trdQDen, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= trdQDen$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("3rd Quartile sigmaT")
ShE4thQ <- ggplot(frthQDen, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= frthQDen$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("4th Quartile sigmaT")
ShEall <- ggplot(div.abiotic2, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= div.abiotic2$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("SigmaT All")
grid.arrange(ShE1stQ, ShE2ndQ, ShE3rdQ, ShE4thQ, ShEall, ncol=2)
#Simpson Plots
Sim1stQDen <- ggplot(fstQDen, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= fstQDen$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("1st Quartile sigmaT")
Sim2ndQDen <- ggplot(sndQDen, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= sndQDen$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("2nd Quartile sigmaT")
Sim3rdQDen <- ggplot(trdQDen, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= trdQDen$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("3rd Quartile sigmaT")
Sim4thQDen <- ggplot(frthQDen, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= frthQDen$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("4th Quartile sigmaT")
SimAll <- ggplot(div.abio2, aes(x = S , y = Theta, colour=Simpson)) + geom_point(size=5, alpha=.6, label= div.abio2$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("SigmaT All")
grid.arrange(Sim1stQDen, Sim2ndQDen, Sim3rdQDen, Sim4thQDen, SimAll, ncol=2)
#All data, All Indices
grid.arrange(ShEall, SEall, SWall, SimAll, Rall, ncol=2)

