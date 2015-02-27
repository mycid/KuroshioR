ColSum<- rank(colSums(test[3:74]))
ColSum<- as.data.frame(cbind(ColSum)
ColSum<-order(ColSum$ColSum)
colnames(Ranks) <- c("Species","rank")


test3 <- as.matrix(test2[1:39,])
d <- dist(test3)
hc <- hclust(t(d))
plot(hc)
head(d)

#Tru <- apply(test2[, 3:74], 1, function(x) (x>0))
#Tru <- t(Tru)
#Tru <- cbind(KC,Tru)
#Sample <- paste(Tru[,1], Tru[,2], sep="_")
#Tru <- cbind(Sample, Tru)


