colors <- colorRampPalette(c("black", "#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "purple", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "#7F0000"))

#plot for Richness, PoDen and Cell Count
ADR <- ggplot(Adiv.abiotic2, aes(x=PoDen, y=SpA, colour=Richness))+ geom_point(size=5, alpha=.6, label=Adiv.abiotic2$ShannonWiener, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), na.value="black", space="rgb", guide="colourbar")
ADR <- ADR+labs(x="Potential Density", y="Total Sample Cell Count")
ADR <- ADR+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
ADR<- ADR+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ADR <- ADR+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
#plot for Unimodal relationship 
ADA <- ggplot(IntDiverse, aes(x =DfromF, y = Cellcount, colour=Richness))+ geom_point(size=5, alpha=.6, label=IntDiverse$Richness, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), na.value="black", space="rgb", guide="colourbar")
ADA <- ADA+labs(x="Total Sample Cell Count", y="Richness")
ADA <- ADA+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
ADA<- ADA+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ADA <- ADA+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
ADA+stat_smooth(method = "lm", size = 1)
#Dino vs. Diatoms
DDP <- ggplot(div.abio2, aes(x = S, y = Theta, colour=DinoPercent))+ geom_point(size=5, alpha=.6, label="Percentage of Diatoms", position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), limits=c(0,1), space="rgb", guide="colourbar")
DDP <- DDP+labs(x="Salinity", y="Theta")
DDP <- DDP+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
DDP<- DDP+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
DDP <- DDP+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
#Distance from the front
DF <- ggplot(Upperdepth, aes(x = DfromF, y =Richness , colour=SpA))+ geom_point(size=5, alpha=.6, label=Upperdepth$Evenness.SW, position=position_dodge(), stat="identity")+ scale_color_gradientn(colours=jet.colors(100), na.value="black", space="rgb", guide="colourbar")
DF<- DF+labs(x="Distance from the Front", y="Total Sample Cell Count")
DF <- DF+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
DF<- DF+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
DF <- DF+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
DF+stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)
#Kuroshio All Zoomed in
AZ <- ggplot(KuroshioWaters, aes(x = S , y = Theta, colour=Richness))+ geom_point(size=7, alpha=.5, label=Adiv.abiotic2$ShannonWiener, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100),limits=c(0,1), na.value="black", space="rgb", guide="colourbar")
AZ<- AZ+labs(x="Salinity", y="Theta")
AZ <- AZ+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
AZ<- AZ+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
AZ <- AZ+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
#Kuroshio All
ALK <- ggplot(Adiv.abiotic2, aes(x =longitude, y = latitude, colour=S))+ geom_point(size=6, alpha=.6, label=IntDiverse$Richness, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(7), limits=c(33.5,34.8), space="rgb", guide="colourbar")
ALK <- ALK+labs(x="Richess", y="Eveness")
ALK <- ALK+theme(legend.title = "Shannon Wiener", element_text(colour="chocolate", size=21, face="bold"))
ALK <- ALK+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ALK <- ALK+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
ALK+stat_smooth(method = "lm", formula = y ~ poly(x, 2) , size = 1)
#Euclidean difference plot
EU <- ggplot(IntDiverse, aes(x =longitude, y = latitude, colour=factor(groups)))+ geom_point(size=7, label=groups, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), limits=c(15,26), space="rgb", guide="colourbar")
EU <- EU+labs(x="Richess", y="Eveness")
#EU <- EU+theme(legend.title = "Shannon Wiener", element_text(colour="chocolate", size=21, face="bold"))
EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
EU <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
EU+stat_smooth(method = "lm", formula = y ~ poly(x, 2) , size = 1)

 