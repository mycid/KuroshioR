<<<<<<< HEAD
#Neutral Parameters calculation function. Data input MUST be a matrix for the function to work. 
#S=Shannon Wiener, J=number of individuals, 
centitest<- as.matrix(centitest[4:74]) #Converts centitest to propper format and gets rid of clusters
DIPhyt<- DIphyto[4:74]
centiStation <- as.matrix(centiStation)
  
  
#THIS FUNCTION IS NO LONGER NEEDED. WORK IS BEING DONE BY PROGRAM TETAME
TetameResults <- read.csv("TetameResults.csv")
#NeutralParams <- function (data) 
#{  
#   S<- apply(data, 1, function(x) sum(x>0))
#   J<- apply(data, 1, function(x) sum(x))
#   H<- apply(data, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
#   H <- colSums (H, na.rm=T) #Shannon Wiener Diversity Index 
#   Neut <- cbind(S, J, H)
#   Neutral <- data.frame(matrix(ncol=2, nrow=nrow(data)))
#   colnames(Neutral) <- c(x="Theta", y="m")
#    for(i in 1:nrow(data)) {
#   d<- c("Theta", "m")
#   Neutral <- data.frame(x=2,y=39)
#   for(i in 1:nrow(data)) {
    
#     l<- logkda.R(data[i,], use.brob=TRUE)  # Use logkda() if pari/gp is available
#     z<- optimal.params(data[i,], log.kda=l)
#    Neutral[i,] <- as.matrix(z)
#   }
#   colnames(Neutral) <- c(x="Theta", y="m")
#   Neutral <- cbind(Neut, Neutral)
#   return(Neutral)
#}}
#Neutral<- NeutralParams(centiStation)#Comand to create the dataframe
#In this case centitest is every sample

#Every sample is is to sparse to create Neutral parameters for
#I need to do every station instead. 
#Can't use itegrated depth because that does not give exact number of individuals
#Station test will create a new dataframe with all depth samples added
#This creates a single station sample that contains all samples added over that station

Stationtest<- function(data) { #Data must have depth in first row, and a column STATION written in all caps
  SL <- (unique(data[,"STATION"]))
  centiStation <- matrix(0,length(SL),ncol(data))
    
    for(i in 1:length(SL)){
    Station <- toString(SL[i])
    datas <- data[which(data$STATION==Station),]
      centi<- colSums(datas[,4:74])
      print(centi)
      centiStation[i,] <- as.matrix(centiStation[i,]+centi)
    }
  return(centiStation)}
centiStation <- Stationtest(DIphyto)#Diphyto is a data frame made for depth integration from the loopdepthinteg.R script
#Build 

 NeutralTest <- function(data, data2, rep) #data & data2 must have the same number of rows for this to work.
  #data & data2 must have only abundances of organisms for each sample, nothing else.
{ startt=Sys.time()
  P <- as.numeric(rep(0, nrow(data)))
  ci <- as.numeric(rep(0, nrow(data)))
  reps <- c(1:rep)
  Hlist <- as.numeric(rep(0, nrow(data)))
  Ho <- apply(data2, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
  Ho <- colSums(Ho, na.rm=T)
  Neutral_MeanH <- as.numeric(rep(0, nrow(data)))
  Station <- data[, "Station"]
  depth <- data[, "depth"]
  
  
  for (i in 1:nrow(data)) {
    S <- data[i, "S"]
    Theta <- data[i, "Theta"]
    m <- data[i, "m"]
    J <- data[i, "J"]  
    Hlist <- as.numeric(rep(0, rep))
    
    for(z in 1:rep) {
   #Nsim <- matrix(nrow=rep, ncol=S)
      repeat{
        Nsim <- rand.neutral(J, Theta)
        s    <- length(Nsim) 
        if (s==S) break
        }
       #Nsim[z,] <- as.numeric(Ns)
    
      H <- (Nsim/sum(Nsim))*(-log(Nsim/sum(Nsim)))
      H<- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
      Hlist[z] <- H
      }
  
    Neutral_MeanH[i] <- mean(Hlist)
    p<- t.test(Hlist, mu=Ho[i])
    ci[i]<- list(p[4])
    P[i] <- as.numeric(p[3]) # Ho: mu=3
    endt=Sys.time()
    diff<- endt-startt
    #print(P[i])
    print(paste("station:", Station[i], "depth:", depth[i], "p-value:", P[i], "CI:", ci[i], "H observed:", Ho[i], "mean H:", Neutral_MeanH[i], "calc.time:", Sys.time()))
    }
 NTest <-  data.frame(Station, depth, P, Neutral_MeanH, Ho) 
 return(NTest) 
}
#H <- Create.Diversity.I(centitest, 4, 82  )
#H <- H[-c(97:99), "ShannonWiener"] #Must remove these rows for now because they had all decimals for their abundances
# Adds Shannon values to the tetame data so a test can be performed
NeutralTest(CleanTetame, centitest, rep=10)
centitest <- centitest[-c(97:99),]

--------------------------------------------------------------
# scripts sent by Chust

mystation <- data.frame(centitest[1,])
mycount <- as.count(mystation)
plot(mycount)
preston(mycount,orig=TRUE)
jj <- c(3.74597,  0.45542)
r <- plot(preston(mycount, n=length(volkov(no.of.ind(mycount), jj, bins=TRUE)), orig=TRUE), ylab="Number of species",cex.lab=1.5, las=2)
points(r,volkov(no.of.ind(mycount), jj, bins=TRUE),type="b", pch=19)
------------------------------------------------------------------
# stnum
stationinfo<- Adiv.abiotic[-c(97, 98,99), c("station", "depth..m.")]
tetame <- cbind(stationinfo, TetameResults)

CleanupTetame <- function(tetame, cutoff){
J <- tetame[,"J"]
S <- tetame[,"S"]
Theta <- as.numeric(rep(0, nrow(tetame)))
m     <- as.numeric(rep(0, nrow(tetame)))
Theta2<- as.numeric(rep(0, nrow(tetame)))
m2     <- c( as.numeric(rep(0, nrow(tetame))))
Station <- tetame[,"station"]
depth <- tetame[, "depth..m."]
  
for (i in 1:nrow(tetame)) {
  theta1 <- tetame[1, "Theta"]
  theta2 <- tetame[i, "Theta2"]
  if (theta1 < theta2 && theta1<cutoff) {
       Theta[i] <- theta1
       m[i] <- tetame[i, "m"]
  } else
        if(theta2< theta1 & theta2<cutoff){
            Theta[i] <- theta2
            m[i] <- tetame[i, "m2"]
        } else
              if(theta1==theta2 & theta1<cutoff) {
                 Theta[i] <- theta1
                 m[i] <- tetame[i, "m"]
              } else
                  if (theta1 & theta2 >cutoff) {
                      Theta[i] <- theta1
                      m[i] <- tetame[i, "m"]
                      Theta2[i] <- theta2
                      m2[i] <- tetame[i, "m2"]}
}
cbind(Station, depth, S, J, Theta, m, Theta2, m2)
}
CleanTetame <- CleanupTetame(tetame, 100)



   