colors <- colorRampPalette(c("black", "#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "purple", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "#7F0000"))

#plot for Richness, PoDen and Cell Count)
ADR <- ggplot(D, aes(x=latitude, y=depth, colour=sigPoDen))+ geom_point(size=7, alpha=.5, label=Adiv.abiotic2$Richness, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), na.value="white", space="rgb", guide="colourbar")
ADR <- ADR+labs(x="latitude", y="depth")
ADR <- ADR+geom_text(aes(label=STATION),hjust=0, vjust=-.5)
ADR <- ADR+ scale_y_reverse( lim=c(150,0))
ADR <- ADR+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
ADR<- ADR+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ADR <- ADR+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
ADR
#plot for Unimodal relationship 
ADA <- ggplot(, aes(x =Evenness.SW, y =SpA, colour=Theta))+ geom_point(size=7, alpha=.6, label=Aldiv.abiotic$Richness, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), na.value="black", space="rgb", guide="colourbar")
ADA <- ADA+labs(x="Evenness.SW", y="Total Sample Cell Count")
ADA <- ADA+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
ADA<- ADA+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ADA <- ADA+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
ADA+stat_smooth(method = "lm")
#Dino vs. Diatoms
DDP <- ggplot(Adiv.abiotic2, aes(x =Richness, y = bigtest$Diatoms, colour=depth))+ geom_point(size=7, alpha=.6, label="Percentage of Diatoms", position=position_dodge(), stat="identity")+ scale_color_gradientn(colours=jet.colors(100), limits=c(0,100), space="rgb", guide="colourbar")
DDP <- DDP+labs(x=Richness, y=cells/mL)
DDP <- DDP+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
DDP<- DDP+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
DDP <- DDP+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
#Distance from the front
DF <- ggplot(Upperdepth, aes(x =SpA, y =Richness, colour=GroupSpeciesRatio))+ geom_point(size=5, alpha=.6, label=SaltyShal$GroupSpeciesRatio, position=position_dodge(), stat="identity")+ scale_color_gradientn(colours=jet.colors(100), na.value="black", limits=c(1,100), space="rgb", guide="colourbar")
DF<- DF+labs(x="Distance from the Front", y="Cell count")
DF <- DF+theme(legend.title = element_text(colour="chocolate", size=20, face="bold"))
DF<- DF+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
DF <- DF+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
DF+stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)
#Kuroshio All Zoomed in
#b <- Adiv.abiotic2[which(test2$Prorocentrum.sp.>0),]
#phyto <-as.numeric(test2[test2$Prorocentrum.sp.>0, which(colnames(test2)=="Prorocentrum.sp.")])
AZ <- ggplot(Adiv.abiotic, aes(x =S, y =Theta, colour=Chlorophyll))+ geom_point(size=9, alpha=.6, label="Cell Count", position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(100), na.value="transparent", space="rgb", name="Richness",guide="colourbar")
#AZ<- AZ+labs(x="latitude", y="depth")+ggtitle("")
AZ<- AZ++geom_text(aes(label=STATION))
AZ <- AZ+theme(panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = 'black'), 
      legend.title = element_text( colour="chocolate", size=25, face="bold"), legend.key.height=unit(3, "cm"), legend.key.width=unit(1, "cm"), 
      legend.background = element_rect(colour = 'black', fill = 'black'), plot.margin = unit(c(1,.5,.5,.5), "cm"), legend.text=element_text(color="white", vjust=.5, size=21) )
 AZ<- AZ+theme(axis.title.x = element_text(color="cadetblue", vjust=-0.35, size=30, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=30, face="bold"))
AZ <- AZ+theme(axis.text.x=element_text(size=25, color="white", vjust=0.5), axis.text.y=element_text(size=25, color="white", vjust=.5))
par(mfrow=c(1,2))
#AZ<- AZ+ scale_y_reverse( lim=c(150,0))
 AZ+stat_smooth(data=Upperdepth$S, method="loess", )
#geom_segment(aes(x = 2, y = 15, xend = 3, yend = 15)
#Kuroshio All
ALK <- ggplot(Adiv.abiotic, aes(x =S, y = Theta, colour=Cl))+ geom_point(size=6, alpha=.6, label=Aldiv.abiotic$Theta, position=position_dodge(), stat="identity", )+ scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
ALK <- ALK+labs(x="Richess", y="Evenness")
ALK <- ALK+theme(legend.title = "Shannon Wiener", element_text(colour="chocolate", size=21, face="bold"))
ALK <- ALK+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
ALK <- ALK+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
ALK+stat_smooth(method = "lm", formula = y ~ poly(x, 2) , size = 1)
#Euclidean difference plot

EU <- ggplot(Adiv.abiotic3, aes(x =S, y =Theta, colour=factor(mass), label=STATION))+ geom_point(size=7, alpha=.9, label=clade, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
EU <- EU+geom_text(aes(label=STATION),hjust=0, vjust=-.5)
EU <- EU+labs(x="Salinity", y="Theta")
#EU<- EU+ scale_y_reverse( lim=c(150,0))
#EU <- EU+theme(legend.title = "Shannon Wiener", element_text(colour="chocolate", size=21, face="bold"))
EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
EU <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
EU
#Salinity profile from front
AZZ <- ggplot(Adiv.abiotic2, aes(x =DfromF, y =S))+geom_smooth(method=loess)
AZZ <- AZZ+theme(panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = 'black'), 
               legend.title = element_text( colour="chocolate", size=25, face="bold"), legend.key.height=unit(4, "cm"), legend.key.width=unit(1, "cm"), 
               legend.background = element_rect(colour = 'black', fill = 'black'), plot.margin = unit(c(1,.5,.5,.5), "cm"), legend.text=element_text(color="white", vjust=.5, size=21) )
AZZ<- AZZ+theme(axis.title.x = element_text(color="cadetblue", vjust=-0.35, size=30, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=30, face="bold"))
AZZ <- AZZ+theme(axis.text.x=element_text(size=25, color="white", vjust=0.5), axis.text.y=element_text(size=25, color="white", vjust=.5))
AZZ
lay_out(list(AZ, 1, 1), list(AZZ, 1, 2))
