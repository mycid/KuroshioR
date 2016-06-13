#Neutral Test on the linux server machine#
library("untb")
CleanTetame <- read.csv("CleanTetame.csv")
CleanTetame <- CleanTetame[which(CleanTetame[,6]<6),]
centitest <- read.csv("centitest.csv")
centitest <- centitest[-c(97:99),]
CleanTetame1 <- CleanTetame[1:20,]
CleanTetame2 <- CleanTetame[21:40,]
CleanTetame3 <- CleanTetame[41:60,]
CleanTetame4 <- CleanTetame[61:80,]
CleanTetame5 <- CleanTetame[81:100,]
CleanTetame6 <- CleanTetame[101:120,]
CleanTetame7 <- CleanTetame[121:140,]
CleanTetame8 <- CleanTetame[141: 156,]

centitest1 <- centitest[1:20,]
centitest2 <- centitest[21:40,]
centitest3 <- centitest[41:60,]
centitest4 <- centitest[61:80,]
centitest5 <- centitest[81:100,]
centitest6 <- centitest[101:120,]
centitest7 <- centitest[121:140,]
centitest8 <- centitest[141: 156,]
source("RDiversityIndicesCalc.R")

#then run NeutralTest(CleanTetame1, centitest1, rep=1000)
