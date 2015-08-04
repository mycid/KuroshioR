#Plotting interesting Kuroshio associated species distance from the Front 
a<- as.data.frame(test$Protoperidinium.pellucidum)
a[a>0]<- 1
a<-cbind(DfromF,a)
a<-a[which(a[,2]==1),]
colnames(a) <- c("DfromFront","species")
#
b<- as.data.frame(test$Rhizosolenia.alata)
b[b>0]<- 2
b<-cbind(DfromF,b)
b<-b[which(b[,2]==2),]
colnames(b) <- c("DfromFront","species")
#
c<- as.data.frame(test$Climacodium.frauenfeldianum) 
c[c>0]<- 3
c<-cbind(DfromF,c)
c<-c[which(c[,2]==3),]
colnames(c) <- c("DfromFront","species")
#
d<- as.data.frame(test$Thalassiothrix.sp.)
d[d>0]<- 4
d<-cbind(DfromF,d)
d<-d[which(d[,2]==4),]
colnames(d) <- c("DfromFront","species")
#e
e<- as.data.frame(test$Cerataulina.pelagica)
e[e>0]<- 5
e<-cbind(DfromF,e)
e<-e[which(e[,2]==5),]
colnames(e) <- c("DfromFront","species")
coolsp<- rbind(a,b,c,d,e)
par(mar=c(5,5,3,3))
#f
f<- as.data.frame(test$Prorocentrum.sp.)
f[f>0]<- 6
f<-cbind(DfromF,f)
f<-f[which(f[,2]==6),]
colnames(f) <- c("DfromFront","species")
#g
g<- as.data.frame(test$Prorocentrum.minimum)
g[g>0]<- 7
g<-cbind(DfromF,g)
g<-g[which(g[,2]==7),]
colnames(g) <- c("DfromFront","species")
#graph
coolsp<- rbind(a,b,c,d,e,f,g)
Factor <- c(as.vector(rep(c("Kuroshio"),each=nrow(rbind(a,b,c,d,e)))), as.vector(rep(c("Mixed"),each=nrow(rbind(f,g)))))
coolsp <- cbind(Factor, coolsp)
par(mar=c(5,5,3,3))
plot(coolsp$DfromFront, coolsp$species, pch=23, col="black",bg="green",cex.lab=1.5, cex=2, xlab="Distance (km) from the Front", ylab="Phytoplankton Species")

PA <- ggplot(coolsp, aes(x =DfromFront, y =species, colour=factor(Factor)))+ geom_point(size=7, alpha=.9, label=clade, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), limits=c(15,26), space="rgb", guide="colourbar")
PA <- PA+labs(x="Distance from the Front", y="Species Presence-Absence")
#EU <- EU+theme(legend.title = "Shannon Wiener", element_text(colour="chocolate", size=21, face="bold"))
PA <- PA+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
PA <- PA+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
PA
