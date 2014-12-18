#Model Selection
library(leaps)
leaps=regsubsets(ShannonWiener~S+Theta+sigPoDen, data=div.abiotic2, nbest=10)
plot(leaps, scale="adjr2")
plot(leaps, scale="bic")
cor.test(div.abiotic4n$ShannonWiener, div.abiotic4n$S )
cor.test(div.abiotic4n$ShannonWiener, div.abiotic4n$sigPoDen)
cor.test(div.abiotic4n$ShannonWiener, div.abiotic4n$sigma_t)

names(div.abiotic2)
modelSW <- lm(ShannonWiener~sigPoDen, data=div.abio2)#Not significant with any single or paired independent variables
modelRich <- lm(Richness~sigPoDen, data=div.abio2) #This value is significant at Theta but more significant with S
modelESW <- lm(Evenness.SW~Theta, data=div.abio2)#This value is significant at Theta but more significant with S
modelES <- lm(Evenness.Sim~sigPoDen, data=div.abio2)#Is significant at Theta but much more significant when including S
modelSim <- lm(Simpson~sigPoDen, data=div.abio2)#only signficant as Theta

Rsq <- c(summary(modelSim)$r.squared, summary(modelSW)$r.squared, summary(modelRich)$r.squared, summary(modelES)$r.squared, summary(modelESW)$r.squared)
Adj.Rsq <- c(summary(modelSim)$adj.r.squared, summary(modelSW)$adj.r.squared, summary(modelRich)$adj.r.squared, summary(modelES)$adj.r.squared, summary(modelESW)$adj.r.squared)
P.values <- c(summary(modelSim)$coefficients[2,4], summary(modelSW)$coefficients[2,4], summary(modelRich)$coefficients[2,4], summary(modelES)$coefficients[2,4], summary(modelESW)$coefficients[2,4])
z <- coef(modelSim)
y <- coef(modelSW)
x <- coef(modelRich)
w <- coef(modelES)
v <- coef(modelESW)
linModel <- c(as.vector(z,y, x,w, v))
modelresults <- data.frame(Rsq, Adj.Rsq, P.values)