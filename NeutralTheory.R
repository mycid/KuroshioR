<<<<<<< HEAD
#Neutral Parameters calculation function. Data input MUST be a matrix for the function to work. 
#S=Shannon Wiener, J=number of individuals, 
#centitest<- as.matrix(centitest[4:74]) #Converts centitest to propper format and gets rid of clusters
#DIPhyt<- DIphyto[4:74]
#centiStation <- as.matrix(centiStation)
  
  
#THIS FUNCTION IS NO LONGER NEEDED. WORK IS BEING DONE BY PROGRAM TETAME
TetameResults <- read.csv("TetameResults.csv")

#Neutral<- NeutralParams(centiStation)#Comand to create the dataframe
#In this case centitest is every sample

 

 NeutralTest <- function(data, data2, rep, t.lim, Theta_Ewens=FALSE) #data & data2 must have the same number of rows for this to work.
  #data & data2 must have only abundances of organisms for each sample, nothing else.
{
  startt= Sys.time()
  P <- as.numeric(rep(0, nrow(data)))
  ci <- as.numeric(rep(0, nrow(data)))
  reps <- c(1:rep)
  Hlist <- as.numeric(rep(0, nrow(data)))
  Ho <- apply(data2, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
  Ho <- colSums(Ho, na.rm=T)
  Neutral_MeanH <- as.numeric(rep(0, nrow(data)))
  Station <- data[, "Station"]
  depth <- data[, "depth"]
  Precise <- as.vector(rep("NA", nrow(data)))
  
  for (i in 1:nrow(data)) {
    if (Theta_Ewens==TRUE) {
      Theta <- data[i, "Theta_Ewens"]
    } else
    {Theta<- data[i, "Theta"]
    }
    S <- data[i, "S"]
    m <- data[i, "m"]
    J <- data[i, "J"]  
    Hlist <- as.numeric(rep(0, rep))
   
    ptm <- proc.time()
    for(z in 1:rep) {
      elapsed <- ptm-proc.time()
      if (elapsed[3]< -t.lim) {
        exact=1
        repeat{
        Nsim <- rand.neutral(J, Theta)
          s    <- length(Nsim) 
        if (s>S-1 || s<S+1) break
        }
          } else {
               repeat{
                          exact=0
                         Nsim <- rand.neutral(J, Theta)
                          s    <- length(Nsim) 
                          if (s==S) break
               }}
      
      H <- (Nsim/sum(Nsim))*(-log(Nsim/sum(Nsim)))
      H<- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
      Hlist[z] <- H
      }

  if(exact==1) {
    Precise[i] <- "FALSE"
  } else {
    Precise[i] <- "TRUE"
  }
    Neutral_MeanH[i] <- mean(Hlist)
    p<- t.test(Hlist, mu=Ho[i])
    ci[i]<- list(p[4])
    P[i] <- as.numeric(p[3]) # Ho: mu=3
    time= ptm-proc.time()
    diff<- endt-startt
    print(paste("station:", Station[i], "depth:", depth[i], "p-value:", P[i], "CI:", ci[i], "H observed:", Ho[i], "mean H:", Neutral_MeanH[i], "calc.time:", time, "precise?:", Precise[i]))
    }
 NTest <-  data.frame(Station, depth, P, Neutral_MeanH, Ho, Precise) 
 return(NTest) 
}
#H <- Create.Diversity.I(centitest, 4, 82  )
#H <- H[-c(97:99), "ShannonWiener"] #Must remove these rows for now because they had all decimals for their abundances
# Adds Shannon values to the tetame data so a test can be performed
NTResults <- NeutralTest(CleanTetame, centitest, rep=10, 30, Theta_Ewens=FALSE)
#centitest <- centitest[-c(97:99),]

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
Tet <- cbind(stationinfo, TetameResults)


CleanTetame <- CleanupTetame(Tet, 100)
Tetame


   