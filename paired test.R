Sampl <- cbind(Adiv.abiotic2[1:11], test)
BP <- cbind(KC, test$Rhizosolenia.stolterfothii, test$Aulacoseira.distans, test$Gymnodiniales, test$Thalassiosira.spp., test$Prorocentrum.minimum,test$Aulacoseira.ambigua,test$Prorocentrum.sp.,test$Thalassionema.nitzschioides,test$Oscillatoriaceae,test$Thalassiosiraceae,test$Cryptomonadaceae,test$Peridiniales)
ggplot(Sampl, aes(x = SpA, y = Richness, colour=Coscinodiscus.sp.)) + geom_point(size=5, alpha=.6, label= Sampl$Gymnodinium.spp., min=0) + scale_colour_gradientn(guide = "colourbar", limits = c(0,2500),breaks=c(2500,2000, 1500, 1250, 1000, 800, 600,400, 200, 100, 50, 0),values=c(2500,2000, 1500, 1250, 1000, 800, 600,400, 200, 100, 50, 0),na.value="black", space="rgb", colours=c("white", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F","green", "yellow", "#FF7F00", "red", "purple", "#7F0000"))+ggtitle("Rhizosolenia.stolterfothii Distribution")
ggplot(Sampl, aes(x=DfromF, y=SpA, colour=Achnanthes.sp.)) + geom_point(size=5, alpha=.5, label=Sampl$Achnanthes.sp.) + scale_color_gradientn(colours=jet.colors(7), limits=c(0,500), na.value="black", space="rgb", guide="colourbar")

scatter3D(z =Sampl$depth, x = Sampl$longitude, y = Sampl$latitude, colvar=Sampl$Aulacoseira.ambigua., col=jet.colors(20), pch = 18, cex = 2, theta = -10, phi = 30, ticktype = "detailed",xlab = "longitude", ylab = "latitude", zlab = "depth")

jet.colors <-colorRampPalette(c("white", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F","green", "yellow", "#FF7F00", "red", "purple", "#7F0000"))
c.min<-0
c.max<-3000
c.n<-190
c.brks<-seq(c.min,c.max,length=c.n)
c.ramp<-jet.colors(c.n)
c.ramp[findInterval(c.min,c.brks)] <-NA

stat_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x,2), size=1)
Rhizosolenia stolterfothii
Aulacoseira distanst
Gymnodiniales ##Greatest arround Front with some trailing off into deeper colder waters
Thalassiosira spp.#patchily distributed everywhere
Proro centrum minimum
Aulacoseira ambigua
Prorocentrum sp.#Interesting, looks strongest in the middle
Thalassionema nitzschioides #Greatest arround Front with some trailing off into deeper colder waters
#Oscillatoriaceae #Hugely front associated! Very high presence, linear path a bit?
Thalassiosiraceae# Wide spread, moderately strong tail end presence
#Cryptomonadaceae#Strongly warm water associated but with definate presence in cold waters too
Peridiniales#All over the place, super abundant, near universal presence

