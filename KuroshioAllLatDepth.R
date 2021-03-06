SWall <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=ShannonWiener)) + geom_point(size=5, alpha=.6, label= div.abio2$ShannonWiener) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio All")
Rall <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=Richness)) + geom_point(size=5, alpha=.6, label= div.abio2$Richness) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio All")
SEall <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=Evenness.Sim)) + geom_point(size=5, alpha=.6, label= div.abio2$Evenness.Sim) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio All")
ShEall <- ggplot(div.abiotic2, aes(x = latitude , y = depth, colour=Evenness.SW)) + geom_point(size=5, alpha=.6, label= div.abiotic2$Evenness.SW) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio All")
SimAll <- ggplot(div.abio2, aes(x = latitude , y = depth, colour=Simpson)) + geom_point(size=5, alpha=.6, label= div.abio2$Simpson) + scale_color_gradientn(colours=jet.colors(7), space ="rgb", na.value = "black", guide = "colourbar")+ggtitle("Kuroshio All")
grid.arrange(ShEall, SEall, SWall, SimAll, Rall, ncol=2)
#
ArcD$Location <- as.factor(ArcD$Location)
SWArc <- ggplot(ArcCompare, aes(x = Location , y = ShannonWiener, color=Group)) + geom_point(size=7, alpha=.7, label= ArcCompare$ShannonWiener)+ggtitle("Southern Ocean All")
RArc <- ggplot(ArcCompare, aes(x = Location , y = Richness, colour=Group)) + geom_point(size=7, alpha=.7, label= ArcCompare$Richness) +ggtitle("Southern Ocean All")
SEArc <- ggplot(ArcCompare, aes(x = Location , y = Evenness.Sim, colour=Group)) + geom_point(size=7, alpha=.7, label= ArcCompare$Evenness.Sim) +ggtitle("Southern Ocean All")
ShEArc <- ggplot(ArcCompare, aes(x = Location , y = Evenness.SW, colour=Group)) + geom_point(size=7, alpha=.7, label= ArcCompare$Evenness.SW) +ggtitle("Southern Ocean All")
SimArc <- ggplot(ArcCompare, aes(x = Location, y =Simpson, colour=Group)) + geom_point(size=7, alpha=.7, label= ArcCompare$Simpson)+ggtitle("Southern Ocean All")
grid.arrange(ShEArc, SEArc, SWArc, SimArc, RArc, ncol=2)
