#Making a Euclidean distance plot for IntDiverse
#mydata <- na.omit(phytoInt)
#mydata <- scale(phytoInt)
Euclid <- cbind(test[c(1,3:length(test))])
colnames(Euclid) <- c("STATION", as.list(colnames(origin[3:(length(origin)-1)])))
Euclid[Euclid>0]<- 1 
d <- distance(Euclid[2:length(Euclid)], method = "chord") # distance matrix
fit <- hclust(as.dist(d), method="ward.D") 
#mfit <- Mclust(Euclid)
#dend <- plot(fit, hang=-1)
library(sparcl)
# colors the leaves of a dendrogram
clade = cutree(fit, 3)
ColorDendrogram(fit, y = clade, labels = names(Adiv.abiotic2$), main = "Species only", branchlength = 4) 

# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

#Determine number of clusters
Kmean<- cbind(test[3:length(test)])
wss <- (nrow(Kmean)-1)*sum(apply(Kmean,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(Kmean, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
#
# K-Means Cluster Analysis
fit <- kmeans(Kmean, 14)
# get cluster means 
aggregate(Kmean,by=list(fit$cluster),FUN=mean)
# append cluster assignment
Adiv.abiotic2 <- data.frame(Adiv.abiotic2, fit$cluster)
#for origin and DIp.or 4, maybe 5 clusters
#for test and DIp ~13 clusters
hist(fit$cluster)


phytoIn <- t(phytoInt)


rad <- radfit(phytoInt)

vdist <- vegdist(phytoInt, Binary=TRUE)
summary(rad)

#heat map
dstat <- Adiv.abiotic2[2]
Euclid <- cbind(dstat, Euclid)
AEuclid<-Euclid[which(Euclid[,1]=="A"),]
AEuclid <- AEuclid[3:length(AEuclid)]
Col <- colSums(AEuclid[,1:ncol(AEuclid)])>0
Apop <- AEuclid[,which(Col==TRUE)]
#Apop <- Apop[2:length(Apop)]
d <- distance(Apop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Apop <- as.matrix(Apop)
plot(fit) # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Apop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Apop ) ), labels= colnames( Apop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Apop ) ), labels= (A$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Apop ) ), labels= (A$STATION), las= 1)
#B presence absence image
BEuclid<-Euclid[which(Euclid[,1]=="B"),]
BEuclid <- BEuclid[3:length(BEuclid)]
Col <- colSums(BEuclid[,1:ncol(BEuclid)])>0
Bpop <- BEuclid[,which(Col==TRUE)]
#Apop <- Apop[2:length(Apop)]
d <- distance(Bpop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Bpop <- as.matrix(Bpop)
plot(fit) # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Bpop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Bpop ) ), labels= colnames( Bpop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$STATION), las= 1)
#C Transect presence absence image
CEuclid<-Euclid[which(Euclid[,1]=="C"),]
CEuclid <- CEuclid[3:length(CEuclid)]
Col <- colSums(CEuclid[,1:ncol(CEuclid)])>0
Bpop <- BEuclid[,which(Col==TRUE)]
d <- distance(Bpop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Bpop <- as.matrix(Bpop)
plot(fit) # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Bpop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Bpop ) ), labels= colnames( Bpop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$STATION), las= 1)

# D Transect

DEuclid<-Euclid[which(Euclid[ 2]=="D"),]
Col <- colSums(DEuclid[,8:ncol(DEuclid)])>0
DEuclid<- DEuclid[,8:length(DEuclid)]
Dpop <- DEuclid[, which(Col==TRUE)]
Dpop <- as.matrix(Dpop)
par(mar=c(5,10,4.1,2.1))
image4 <- image(Dpop, xaxt= "n", yaxt= "n", col = grey(seq(0, 1)))
image4<- image4+axis( 2, at=seq(0,1,length.out=ncol( Dpop ) ), labels= colnames( Dpop ), las= 2 )
image4<- image4+axis( 1, at=seq(0,1,length.out=nrow( Dpop ) ), labels= rownames( Dpop ), las= 2)
#E Transect
EEuclid<-Euclid[which(Euclid[ 2]=="E"),]
Col <- colSums(EEuclid[,8:ncol(EEuclid)])>0
EEuclid<- EEuclid[,8:length(EEuclid)]
Epop <- EEuclid[,which(Col==TRUE)]
Epop <- as.matrix(Epop)
par(mar=c(5,10,4.1,2.1))
image5 <- image(Epop, xaxt= "n", yaxt= "n", col = grey(seq(0, 1)))
image5 <- image5+axis( 2, at=seq(0,1,length.out=ncol( Epop ) ), labels= colnames( Epop), las= 2 )
image5<- image5+axis( 1, at=seq(0,1,length.out=nrow( Epop ) ), labels= rownames( Epop), las= 2)

#presence-absence/abundance matrix
prabmatrix <- prabinit(prabmatrix=DIp.or[8:length(DIp.or)], rows.are.species=FALSE)
prabt <- prabtest(prabmatrix)
summary(prabt)
library(xlsx)
write.xlsx(DIp.or, "/home/trevor/Desktop/2009KuroshioData/KuroshioR/DIp.or.xlsx")
getwd()
