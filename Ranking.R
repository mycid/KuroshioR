ColSum<- rank(colSums(test[3:74]))
ColSum<- as.data.frame(cbind(ColSum)
ColSum<-order(ColSum$ColSum)
colnames(Ranks) <- c("Species","rank")

