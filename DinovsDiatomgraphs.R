#Richness and Shannon Evenness for dinos and diatoms in the Arctic and on the Kuroshio 

#Kuroshio Dino vs. Diatoms vs. All
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
DiatomR <- ggplot(Diatomdiv.abio2, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= Diatomdiv.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Diatom Richness")
DinoR <- ggplot(dinodiv.abio2, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= dinodiv.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio dinoflagellate Richness")
BothR <- ggplot(div.abio2, aes(x = S , y = Theta, colour=Richness)) + geom_point(size=5, alpha=.6, label= div.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Both Kuroshio Richness")
DiatomSimE <- ggplot(Diatomdiv.abio2, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= Diatomdiv.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", guide = "colourbar")+ggtitle("Kuroshio Diatom Shannon Evenness")
DinoRSimE <- ggplot(dinodiv.abio2, aes(x = S , y = Theta, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= dinodiv.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "blue", guide = "colourbar")+ggtitle("Kuroshio dinoflagellate Shannon Evenness")
BothRSimE <- ggplot(div.abio2, aes(x = S , y = Theta, colour=DinoandDiatomP)) + geom_point(size=5, alpha=.6, label= div.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "#00007F", guide = "colourbar")+ggtitle("Both Kuroshio Shannon Evenness")
grid.arrange(DiatomR, DinoR, DiatomSimE, S=DinoRSimE, BothR, BothRSimE, ncol=2)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
DiatomR <- ggplot(Diatomdiv.abio2, aes(x = latitude, y = depth, colour=Richness)) + geom_point(size=5, alpha=.6, label= Diatomdiv.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio Diatom Richness")
DinoR <- ggplot(dinodiv.abio2, aes(x = latitude , y = depth, colour=Richness)) + geom_point(size=5, alpha=.6, label= dinodiv.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio dinoflagellate Richness")
BothR <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=Richness)) + geom_point(size=5, alpha=.6, label= div.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Both Kuroshio Richness")
DiatomSimE <- ggplot(Diatomdiv.abio2, aes(x = latitude , y = depth, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= Diatomdiv.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", guide = "colourbar")+ggtitle("Kuroshio Diatom Shannon Evenness")
DinoRSimE <- ggplot(dinodiv.abio2, aes(x = latitude , y = depth, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= dinodiv.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "blue", guide = "colourbar")+ggtitle("Kuroshio dinoflagellate Shannon Evenness")
BothRSimE <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= div.abio2$EvennessSW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "#00007F", guide = "colourbar")+ggtitle("Both Kuroshio Shannon Evenness")
grid.arrange(DiatomR, DinoR, DiatomSimE, S=DinoRSimE, BothRSimE, BothR, ncol=2)

filled.contour(interp(div.abio2$latitude, div.abio2$depth, div.abio2$S, duplicate="median"), color.palette = heat.colors)
#Southern Ocean Data sin abiotic 
DinoArcRich <- plot(dinoArcD$Location, DinoArcD$Richness)
DiatomArcR <- boxplot(DAD$location, DAD$Richness)
AllArcR <- plot(ArcD$location, ArcD$Richness)
DinoArcE <- plot(dinoArcD$location, DinoArcD$Evenness.SW)
DiatomArcE <- plot(DAD$location, DAD$Evenness.SW)
AllArcE <- plot(ArcD$location, ArcD$Evenness.SW)
grid.arrange(DinoArcR,DinoArcE,DiatomArcR,DiatomArcE,AllArcR,AllArcE, ncol=2)

