#Making diversity Indices for any data frame with species and abundance
  #s=the column starting point in the data frame for the function
  #e=the end column in the data frame for the function
  #name=name of the data.frame to work from

## @knitr Create.Diversity.I
Create.Diversity.I <- function(name, s, e) {

  SpR <- apply(name[, s:e], 1, function(x) sum(x>0)) #species richness all phytoplankton not just diatom and dinoflagellate
  Abund <-apply(name[, s:e], 1, function(x) sum(x)) #cell count
  SD <- apply(name[, s:e], 1, function(x) (sum(x*(x-1)))/(sum(x)*(sum(x)-1))) #Simpson Diversity
  SimE <- (1/SD)/SpR #Simpson Evenness
  SW <- apply(name[, s:e], 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
  SWD <- colSums (SW, na.rm=T) #Shannon Wiener Diversity Index 
  ShannonE <- SWD/log(SpR) #Shannon Wiener Evenness 
  return(data.frame(Richness=SpR, ShannonWiener=SWD, Simpson=SD, Evenness.SW=ShannonE, Evenness.Sim=SimE, Cellcount=Abund))
}
  
#Calculating Potential Temperature and Salinity
## @knitr PotentialST
PotentialST <- function(d, t, s) { #t=temperature, s=salinity, d=depth
T <- t-((.1*d)/500) #Potential temperature
p0<- 999.842594+6.793952*10^(-2)*T-9.095290*10^(-3)*T^(2)+1.001685*10^(-4)*T^(3)-1.120083*10^(-6)*T^(4)+6.536332*10^(-9)*T^(5)
A <- 8.24493*10^(-1)-4.0899*10^(-3)*T+7.6438*10^(-5)*T^(2)-8.2467*10^(-7)*T^(3)+5.3875*10^(-9)*T^(4)
B <- (-5.72466*10^(-3))+1.0227*10^(-4)*T-1.6546*10^(-6)*T^(2)
C <- 4.8314*10^(-4)
PoDen <- p0+A*s+B*s^(1.5)+C*s^(2)
#Lines above are the standard, complex formula for ataining the Potential desnity
sigPoDen <- PoDen-1000 #gets sigma
return(data.frame(All[,1:7], Theta=T, S=All[, 9], sigPoDen, All[, 11:31])) #Creates the final product data frame
}
#Group stations by abiotic factor
## @AbioticCluster
AbioticCluster <- function(data, clade) { #Creates a cluster analysis for grouping samples by quantitative, continuous physical features
  #data is the data.frame to be used for clustering
  #clade allows the user to specify the number of clades from the dendragrams to be extracted
  
  data <- scale(data) #Scales the data to a mean of 0, sd=1, absolute value
  d <- dist(data, method="euclidean")  # distance matrix
  par(mar=c(3,4,1,1)+.1) #Creates plot margins
  par(mfrow=c(1,3)) #Creates space for 3 plots
  csin <- hclust(d, method="single") #1st method for clustering
  Gcsin <- cutree(csin, clade) #creates clades
  ColorDendrogram(csin, y = Gcsin, labels = names(Gcsin), main = "Single Method", branchlength = 2) 
  ccom <- hclust(d, method="complete") #2nd method for clustering
  Gccom <<- cutree(ccom, clade) #Creates three clades
  ColorDendrogram(ccom, y = Gccom, labels = names(Gccom), main = "Complete Method", branchlength = 2)
  caver <- hclust(d, method="aver") #third method for clustering
  Gcaver <<- cutree(caver, clade) # Creates clades
  clades <- cbind(Gcsin, Gccom, Gcaver) #Defines clades
  assign("clades", clades, envir = globalenv()) #this brings clades into the global environment
  return(ColorDendrogram(caver, y = Gcaver, labels = names(Gcaver), main = "Aver Method", branchlength = 2)) #returns a plot for caver
  } 
  
  #Uses ggplot to plot x and y with a factor in color

NPfactorInColor <- function(data, xvar="", yvar="", factors="", xlab="", ylab="", title=""){ #Plots the clade data for abiotic
  #data is the data frame, x is the x axis, y is the y axis, factor is the factor to be put in color
  myColors <- brewer.pal(3, "Set1")
  names(myColors) <- levels(factor)
  colScale <- scale_colour_manual(name = "Cluster",values = myColors)
  EU <- ggplot(data, aes_string(x=xvar,y=yvar, colour=aes(factor(factors))),  environment = environment())+geom_point(size=7, alpha=.9, position=position_dodge(), stat="identity")
  EU <- EU+labs(x=xlab, y=ylab, colour="Cluster")+ggtitle(title)+theme_bw(base_size = 12, base_family = "Helvetica")
  #theme(panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = 'blue'),
  #legend.title = element_text( colour="black", size=33, face="bold"), legend.key.height=unit(1, "cm"), legend.key.width=unit(1, "cm"), 
  #legend.background = element_rect(colour = 'white', fill = 'white'), plot.margin = unit(c(1,.5,.5,.5), "cm"), legend.text=element_text(color="black", vjust=.5, size=21) )
  EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
  CF <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
  return(CF)
}

  #Function for splitting data by numbered factors
##@knitr SplitData
  SplitData <- function(data, column, header) { #data=dataframe, column=the the number of the column desired to be used as a numerical factor, header is the title of the data frames without the numbers at the end. MUST BE IN PARATHESES
  libr <- setNames(split(data, data[, column]), paste0(header, unique(data[,column]))) #Seperates the dataframe into multiple dataframes based on the factor
  list2env(libr, globalenv()) #Exports both s and each individual, new data frame to the global environment 
  }
  
#Experimenting with functions for column division

orderlySpliting <- function(x, column, n, header) { 
#Function evenly (as possible)  
#divides data into a set number of groups
#Orders data into a numerical sequence
#Use just the column name for 'column', no quotes
#Returns new dataframes to Global 
#header choses a header name for the new dataframes
  x <- x[order(x[,column]),]
  Depths <- split(x,rep(1:n, ceiling(length(x)/n),length.out = length(x)))
  Depths <- setNames(Depths, paste0(header, unique(rep(1:n))))
  list2env(Depths, globalenv())
}
# ordering a triangular heat map for correlation
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)}
#### An attempted hack of preston lines. Hasn't worked so far.
lines.preston <-
  function(x, xadjust = 0.5, ...)
  {
    oct <- as.numeric(names(x)) - xadjust 
    lines(oct, x, ...)
  }

#Function based off of Ettiene's sampling formula 
#I believe this is correct and it seems to work well though I don't understand why acestrylabelind exists
NeutralCommunityAbund <- function(J, theta, m, other.info=FALSE) {
  specieslabelind<- vector(mode="numeric", length=J)
  ancestrylabelind<- vector(mode="numeric", length=J)
  specieslabelanc <- vector(mode="numeric", length=J)
  abund <-numeric()
  
  if(m < 1){
    I<- (m*(J−1))/(1-m)
  }
  a<-0
  s<-0
  for(i in 2:J) {
    if(m < 1) {
      R1 <- I/(I+i−1)
    }
    else{
      R1<- 1
    }
    x<- runif(1)
    if(x<=R1) {
      a<-a+1
      ancestrylabelind[i] <-a
      R2<- theta/(theta+a-1)
      y <- runif(1)
      if(y<=R2){
        s<-s+1
        specieslabelind[i] <-s
        specieslabelanc[a] <-s
      }
      if(y>R2){
        randA <- sample(1:(a-1),1)
        specieslabelind[i] <- specieslabelanc[randA]
        specieslabelanc[a] <- specieslabelanc[randA] 
      } else{ 
        randj <- sample(1:(i-1),1)
        ancestrylabelind[i] <- ancestrylabelind[randj]
        specieslabelind[i] <- specieslabelind[randj] 
      }
    }}
  if(other.info==TRUE){
    assign("specieslabelanc", specieslabelanc, envir = globalenv())
    assign("ancestrylabelind", ancestrylabelind, envir = globalenv())
    ancS<- c(a,s)
    assign("AncS", ancS, envir=globalenv())
  }
  return(table(specieslabelind))
}


#This code makes a variable length this with as many items as rows of the data frame
#This code is for rarefaction of every sample
#Not needed for my purposes currently
df2numeric <- function(data, incidence=TRUE) {
  names <- c()
  new <- list()
  for(i in 1:nrow(data)) 
  {
    names[i] <- paste("s", i, sep='') 
  }
  if(incidence==FALSE) {
    for(i in 1:nrow(data)) {
      new[[i]] <- as.numeric(as.count(data[i,]))
    }}
  else {
    for(i in 1:nrow(data)) {
      new[[i]] <- c(sum(as.numeric(as.count(data[i,]))),as.numeric(as.count(data[i,])))}
  }
  return(new)
  list2env(names, globalenv())
}
#Function for cleaning up Tetame Results. Picks the most reasonable Tetame value, adds station info. Makes a more orderly dataframe
#Cutoff is the maximum reasonable Theta value expected. Arbitrary aproach.

CleanupTetame <- function(tetame, cutoff){
  tetame <- cbind(c(1:nrow(tetame)), tetame)
  tetame <- tetame[which(tetame$Std_Theta!="inf"),]
  J <- tetame[,"J"]
  S <- tetame[,"S"]
  sampleinfo <- tetame$sampleinfo
  Theta <- as.numeric(rep(0, nrow(tetame)))
  m     <- as.numeric(rep(0, nrow(tetame)))
  Theta2<- as.numeric(rep(0, nrow(tetame)))
  m2     <- c( as.numeric(rep(0, nrow(tetame))))
  Theta_Ewens <- tetame$Theta_Ewens
  Std_Theta <- as.numeric(rep(0, nrow(tetame)))
  
  Station <- tetame[,"station"]
  depth <- tetame[, "depth..m."]
  
  for (i in 1:nrow(tetame)) {
    theta1 <- tetame[i, "Theta"]
    theta2 <- tetame[i, "Theta2"]
    if (theta1 < theta2 && theta1<cutoff) {
      Theta[i] <- theta1
      m[i] <- tetame[i, "m"]
      Std_Theta[i] <- tetame[i, "Std_Theta"]
    } else
      if(theta2< theta1 & theta2<cutoff){
        Theta[i] <- theta2
        m[i] <- tetame[i, "m2"]
        Std_Theta[i] <- tetame[i, "Std_Theta2"]
      } else
        if(theta1==theta2 & theta1<cutoff) {
          Theta[i] <- theta1
          m[i] <- tetame[i, "m"]
          Std_Theta[i] <- tetame[i, "Std_Theta"]
        } else
          if (theta1 & theta2 >cutoff) {
            Theta[i] <- theta1
            m[i] <- tetame[i, "m"]
            Std_Theta[i] <- tetame[i, "Std_Theta"]
            Theta2[i] <- theta2
            m2[i] <- tetame[i, "m2"]}
  }