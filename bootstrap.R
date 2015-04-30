DIp.or <- read.csv("DIp.or.csv")
Euclid[Euclid>0]<- 1 
prabmatrix <- prabinit(prabmatrix=DIp.or[8:length(DIp.or)], rows.are.species=FALSE)
prabt <- prabtest(prabmatrix, times = 1000)
summary(prabt)