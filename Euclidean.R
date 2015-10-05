#Making a Euclidean distance plot for IntDiverse
#mydata <- na.omit(phytoInt)
#mydata <- scale(phytoInt)
Euclid <- test[c(1,3:length(test))]
Euclid <-Adiv.abiotic2[1:190,c(15,17, 19)]
#colnames(Euclid) <- c("STATION", as.list(colnames(origin[3:(length(origin)-1)])))
Euclid[Euclid>0]<- 1 
d <- dist(Euclid, method="ward") # distance matrix
fit <- hclust(d)
#mfit <- Mclust(Euclid)
#dend <- plot(fit, hang=-1)
library(sparcl)
# colors the leaves of a dendrogram
clade = cutree(fit, 5)
ColorDendrogram(fit, y = clade, labels = names(clade), main = "Species only", branchlength = 4) 
x <- identify(fit)
identify(fit, function(k) print(table(Euclid[k,2:72])))
nD <- dev.cur() 
identify(fit, function(k) barplot(table(Euclid[,2:72]),col=1:73), DEV.FUN = nD)
#For Coloring just the branch labels and creating arbitrary color labels
groupCodes <- as.character(Adiv.abiotic2$A)
rownames(Euclid) <- make.unique(groupCodes)
colorCodes <- c(A="red", B="green", C="blue", D="yellow", E="purple")
distSamples <-distance(Euclid[2:length(Euclid)], method = "chord") # distance matrix
fit <- hclust(as.dist(d))
dend <- as.dendrogram(fit)
library(dendextend)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
# Plotting the new dendrogram
plot(dend)
par(cex = 1)
plot(dend[[.5]], horiz = TRUE)



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
dstat <- Adiv.abiotic2[3]
Euclid <- cbind(dstat, Euclid)
AEuclid<-Euclid[which(Euclid[,1]=="A"),]
AEuclid <- AEuclid[3:length(AEuclid)]
 
#Apop <- Apop[2:length(Apop)]
d <- distance(Apop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Apop <- as.matrix(Apop)
plot(fit, main="Transect A") # display dendogram
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
plot(fit, main="Transect B") # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Bpop,  xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Bpop ) ), labels= colnames( Bpop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Bpop ) ), labels= (B$STATION), las= 1)
#C Transect presence absence image
CEuclid<-Euclid[which(Euclid[,1]=="C"),]
CEuclid <- CEuclid[3:length(CEuclid)]
Col <- colSums(CEuclid[,1:ncol(CEuclid)])>0
Cpop <- CEuclid[,which(Col==TRUE)]
d <- distance(Cpop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Cpop <- as.matrix(Cpop)
plot(fit, main="Transect C") # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Cpop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Cpop ) ), labels= colnames( Cpop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Cpop ) ), labels= (C$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Cpop ) ), labels= (C$STATION), las= 1)
# D Transect
DEuclid<-Euclid[which(Euclid[,1]=="D"),]
DEuclid <- DEuclid[3:length(DEuclid)]
Col <- colSums(DEuclid[,1:ncol(DEuclid)])>0
Dpop <- DEuclid[,which(Col==TRUE)]
d <- distance(Dpop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Dpop <- as.matrix(Dpop)
plot(fit, main="Transect D") # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Dpop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Dpop ) ), labels= colnames( Dpop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Dpop ) ), labels= (D$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Dpop ) ), labels= (D$STATION), las= 1)

#E Transect
EEuclid<-Euclid[which(Euclid[,1]=="E"),]
EEuclid <- EEuclid[3:length(EEuclid)]
Col <- colSums(EEuclid[,1:ncol(EEuclid)])>0
Epop <- EEuclid[,which(Col==TRUE)]
#Apop <- Apop[2:length(Apop)]
d <- distance(Epop, method = "chord")
fit <- hclust(as.dist(d), method="ward.D")
Epop <- as.matrix(Epop)
plot(fit, main = "Transect E") # display dendogram
par(mar=c(5,10,4.1,2.1))
image1 <-image(Epop, xaxt= "n", yaxt= "n",lwd=5, useRaster = FALSE,  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol( Epop ) ), labels= colnames( Epop ), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow( Epop ) ), labels= (E$depth), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Epop ) ), labels= (E$STATION), las= 1)

#presence-absence/abundance matrix
prabmatrix <- prabinit(prabmatrix=DIp.or[8:length(DIp.or)], rows.are.species=FALSE)
prabt <- prabtest(prabmatrix)
summary(prabt)
library(xlsx)
write.xlsx(DIp.or, "/home/trevor/Desktop/2009KuroshioData/KuroshioR/DIp.or.xlsx")
getwd()
