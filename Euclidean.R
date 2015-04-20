#Making a Euclidean distance plot for IntDiverse
#mydata <- na.omit(phytoInt)
#mydata <- scale(phytoInt)
Euclid <- cbind(test[c(1,3:length(test))])
colnames(Euclid) <- c("STATION", as.list(colnames(origin[3:(length(origin)-1)])))
Euclid[Euclid>0]<- 1 
Euclid <- cbind(Adiv.abiotic2[c(2,6)], Euclid)
d <- distance(Euclid[2:length(Euclid)], method = "chord") # distance matrix
fit <- hclust(as.dist(d), method="ward.D") 
#mfit <- Mclust(Euclid)
 plot(fit) # display dendogram
groups <- cutree(fit, k=5) 

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
AEuclid<-Euclid[which(Euclid[1 ]=="A"),]
Col <- colSums(AEuclid[,4:ncol(AEuclid)])>0
AEuclid<- AEuclid[4:length(AEuclid)]
Apop <- AEuclid[, which(Col==TRUE)]
Apop <- as.matrix(Apop)
fit <- hclust(as.dist(Apop), method="ward.D") 
plot(fit) # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Apop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Apop ) ), labels= colnames( Apop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Apop ) ), labels= rownames( Apop ), las= 2)
#B presence absence image
 BEuclid<-Euclid[which(Euclid[ 2]=="B"),]
Col <- colSums(BEuclid[,8:ncol(BEuclid)])>0
BEuclid<- BEuclid[8:length(BEuclid)]
Bpop <- BEuclid[, which(Col==TRUE)]
Bpop <- as.matrix(Bpop)
par(mar=c(5,10,4.1,2.1))
image2 <- image(Bpop, xaxt= "n", yaxt= "n", col = grey(seq(0, 1)))
image2 <- image2+axis( 2, at=seq(0,1,length.out=ncol( Bpop ) ), labels= colnames( Bpop ), las= 2 )
image2 <- image2+axis( 1, at=seq(0,1,length.out=nrow( Bpop ) ), labels= rownames( Bpop ), las= 2)
#C Transect presence absence image
CEuclid<-Euclid[which(Euclid[ 2]=="C"),]
Col <- colSums(CEuclid[, 8:ncol(CEuclid)])>0
CEuclid<- CEuclid[, 8:length(CEuclid)]
Cpop <- CEuclid[ , which(Col==TRUE)]
Cpop <- as.matrix(Cpop)
par(mar=c(5,10,4.1,2.1))
image3<- image(Cpop, xaxt= "n", yaxt= "n", col = grey(seq(0, 1)))
image3 <- image3+axis( 2, at=seq(0,1,length.out=ncol( Cpop ) ), labels= colnames( Cpop ), las= 2 )
image3<- image3+axis( 1, at=seq(0,1,length.out=nrow( Cpop ) ), labels= rownames( Cpop ), las= 2)
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
