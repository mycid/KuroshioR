NeutralParams <- function (data) 
{  
   S<- apply(data, 1, function(x) sum(x>0))
   J<- apply(data, 1, function(x) sum(x))
   H<- apply(data, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
   H <- colSums (H, na.rm=T) #Shannon Wiener Diversity Index 
   Neut <- cbind(S, J, H)
   d<- c("Theta", "m")
   Neutral <- data.frame(x=2,y=39)
   for(i in 1:nrow(data)) {
    
     l<- logkda.R(data[i,], use.brob=TRUE)  # Use logkda() if pari/gp is available
     z<- optimal.params(data[i,], log.kda=l)
     Neutral[i,] <- as.matrix(z)
   }
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

  repeat{
    a <- sum(H, na.rm=T) #Shannon Wiener Diversity Index 
    NSim <- rand.neutral(J, Theta, string="spp.")
    H<- (NSim/sum(NSim))*(-log(NSim/sum(NSim)))
    s <- length(NSim)
     if (s!=S) {
      break
    }
  }
  