<<<<<<< HEAD
#Neutral Parameters calculation function. Data input MUST be a matrix for the function to work. 
#S=Shannon Wiener, J=number of individuals, 
centitest<- as.matrix(centitest[4:74]) #Converts centitest to propper format and gets rid of clusters
centiStation <- as.matrix(centiStation)
  
  
=======
>>>>>>> master
NeutralParams <- function (data) 
{  
   S<- apply(data, 1, function(x) sum(x>0))
   J<- apply(data, 1, function(x) sum(x))
   H<- apply(data, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
   H <- colSums (H, na.rm=T) #Shannon Wiener Diversity Index 
   Neut <- cbind(S, J, H)
<<<<<<< HEAD
   Neutral <- data.frame(matrix(ncol=2, nrow=nrow(data)))
   colnames(Neutral) <- c(x="Theta", y="m")
    for(i in 1:nrow(data)) {
=======
   d<- c("Theta", "m")
   Neutral <- data.frame(x=2,y=39)
   for(i in 1:nrow(data)) {
>>>>>>> master
    
     l<- logkda.R(data[i,], use.brob=TRUE)  # Use logkda() if pari/gp is available
     z<- optimal.params(data[i,], log.kda=l)
     Neutral[i,] <- as.matrix(z)
   }
<<<<<<< HEAD
   colnames(Neutral) <- c(x="Theta", y="m")
   Neutral <- cbind(Neut, Neutral)
   return(Neutral)
}
Neutral<- NeutralParams(centiStation)#Comand to create the dataframe
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

NeutralTest <- function(data, rep=1000)
{ 
  P <- as.numeric(rep(0, nrow(data)))
  reps <- c(1:rep)
  Hlist <- as.numeric(rep(0, rep))
  for (i in 1:nrow(data)) {
    S <- data[i, "S"]
    Theta <- data[i, "Theta"]
    m <- data[i, "m"]
    J <- data[i, "J"]
    for(z in reps) {
        repeat{
        Nsim <- rand.neutral(J, Theta, string="spp.")
        s <- length(Nsim) 
        if (s==S) break
        }
    H <- (Nsim/sum(Nsim))*(-log(Nsim/sum(Nsim)))
    H<- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
    Hlist[i] <- H 
    }
     p<- t.test(Hlist, mu=data[i, "H"])
    P[i] <- as.numeric(p[3]) # Ho: mu=3
    print(P[i])
  }
  
  return(P)
 }
NeutralTest(Neutral)
=======
   colnames(Neutral) <- c("Theta", "m")
   Neutral <- cbind(Neut, Neutral)
   return(Neutral)
}
#Build 
NeutralTest <- function(data, reps=1000)
{
  for (i in 1:nrow(data)) {
    S <- data[i, 1]
    Theta <- data[i, "Theta"]
    m <- data[i, "m"]
    J <- data[i, "J"]
    Hdist<- list(replicate(
      repeat{
        NSim <- rand.neutral(J, Theta, string="spp.")
        s <- sum(which(NSim>0))
        if (s!=S) {
      break
        }
    H<- (Nsim/sum(Nsim))*(-log(Nsim/sum(Nsim)))
    H<- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
        }
    ,reps))
    p<- t.test(Hdist,mu=data[i, "H"])
    P[i] <- as.numeric(p[3]) # Ho: mu=3
  }
  return(P)
 }
S <- Neutral[1, "S"]
Theta <- Neutral[1, "Theta"]
m <- Neutral[1, "m"]
J <- Neutral[1, "J"]
>>>>>>> master

  repeat{
    a <- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
    NSim <- rand.neutral(J, Theta, string="spp.")
    H<- (NSim/sum(NSim))*(-log(NSim/sum(NSim)))
    s <- length(NSim)
<<<<<<< HEAD
     if (s==S) {
      break
    }
  }
  
l<- logkda.R(centitest[1,], use.brob=TRUE)  # Use logkda() if pari/gp is available
z<- optimal.params(centitest[1,], log.kda=l)
l1<- logkda.R(centitest[2,], use.brob=TRUE)  # Use logkda() if pari/gp is available
z1<- optimal.params(centitest[2,], log.kda=l1)
data.frame(z,z1)





NeutralTest <- function(Neutral, reps=1000)
{
  for (i in 1:nrow(Neutral)) {
    S <- Neutral[1, "S"]
    Theta <- Neutral[1, "Theta"]
    m <- Neutral[1, "m"]
    J <- Neutral[1, "J"]
   
     Hdist<- list(replicate(
      repeat{
        NSim <- rand.neutral(J, Theta,  string="spp.")
        s <- sum(which(NSim>0))
        if (s!=S) {
          break
        }
        H<- (Nsim/sum(Nsim))*(-log(Nsim/sum(Nsim)))
        H<- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
      }
      ,reps))
    p<- t.test(Hlist,mu=Neutral[1, "H"])
    P[i] <- as.numeric(p[3]) # Ho: mu=3
  }
  return(P)
}
S <- Neutral[1, "S"]
Theta <- Neutral[1, "Theta"]
m <- Neutral[1, "m"]
J <- Neutral[1, "J"]




=======
     if (s!=S) {
      break
    }
  }
>>>>>>> master
  