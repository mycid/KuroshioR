#Neutral Test on the linux server machine#
library("untb")
CleanTetame <- read.csv("CleanTetame.csv")
CleanTetame <- CleanTetame[which(CleanTetame[,6]<15),]
centitest <- read.csv("centitest.csv")
centitest <- centitest[-c(97:99),]
CleanTetame1 <- CleanTetame[1:30,]
CleanTetame2 <- CleanTetame[31:60,]
CleanTetame3 <- CleanTetame[61:90,]
CleanTetame4 <- CleanTetame[91:120,]
CleanTetame5 <- CleanTetame[121:150,]
CleanTetame6 <- CleanTetame[151:182,]
centitest1 <- centitest[1:30,]
centitest2 <- centitest[31:60,]
centitest3 <- centitest[61:90,]
centitest4 <- centitest[91:120,]
centitest5 <- centitest[121:150,]
centitest6 <- centitest[151:182,]
source("RDiversityIndicesCalc.R")

#then run NeutralTest(CleanTetame1, centitest1, rep=1000)
