#Downloading all csv files to the Global Environment
  folder <- getwd() #sets directory as an object/path to folder that holds multiple .csv files
  file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
  # read in each .csv file in file_list and create a data frame with the same name as the .csv file
  for (i in 1:length(file_list)){
    assign(file_list[i], 
         read.csv(paste(folder, '/', file_list[i], sep=''))
  )
  }

#Making diversity Indices for any data frame with species and abundance
  #s=the column starting point in the data frame for the function
  #e=the end column in the data frame for the function
  #name=name of the data.frame to work from
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
PotentialST <- function(d, t, s) { #t=temperature, s=salinity, d=depth
T <- t-((.1*d)/500) #Potential temperature
p0<- 999.842594+6.793952*10^(-2)*T-9.095290*10^(-3)*T^(2)+1.001685*10^(-4)*T^(3)-1.120083*10^(-6)*T^(4)+6.536332*10^(-9)*T^(5)
A <- 8.24493*10^(-1)-4.0899*10^(-3)*T+7.6438*10^(-5)*T^(2)-8.2467*10^(-7)*T^(3)+5.3875*10^(-9)*T^(4)
B <- (-5.72466*10^(-3))+1.0227*10^(-4)*T-1.6546*10^(-6)*T^(2)
C <- 4.8314*10^(-4)
PoDen <- p0+A*s+B*s^(1.5)+C*s^(2)
#Lines above are the standard, complex formula for ataining the Potential desnity
sigPoDen <- PoDen-1000 #gets sigma
return(data.frame(DiversityIndex, All[,1:6], Theta=T, S=All[, 8], sigPoDen, All[, 10:31])) #Creates the final product data frame
}
#Group stations by abiotic factor

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
  #Temporary
  #EU <- ggplot(Adiv.abiotic[which("depth..m."==0)], aes(x =lat, y =depth..m., colour=factor(Gccom)))+ geom_point(size=7, alpha=.9, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
  #EU <- ggplot(Adiv.abiotic[which("depth..m."==0)], aes(x =lat, y =depth..m., colour=factor(Gccom), label=station))+ geom_point(size=7, alpha=.9, label=c, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
  #EU <- EU+geom_text(aes(label=station),hjust=0, vjust=-.5)
  #EU <- EU+labs(x="Salinity", y="Theta")+stat_contour(z=Adiv.abiotic$sigPoDen, binwidth = 2)
  #EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
  #CF <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
  #CF #Temporary

  #Function for splitting data by numbered factors
  SplitData <- function(data, column, header) { #data=dataframe, column=the the number of the column desired to be used as a numerical factor, header is the title of the data frames without the numbers at the end. MUST BE IN PARATHESES
  libr <- setNames(split(data, data[, column]), paste0(header, unique(data[,column]))) #Seperates the dataframe into multiple dataframes based on the factor
  list2env(libr, globalenv()) #Exports both s and each individual, new data frame to the global environment 
  }
 
  #Experimenting with functions for column division

#x=Adiv.abiotic
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
#An attempted hack of logkda.pari(a,...). This function currently not working.

logkda.pari <- function (a, numerical = TRUE)
{
  if ((system("gp --version", intern = FALSE, ignore.stderr = TRUE)) !=
      0) {
    warning("pari/gp not installed: method changed to 'polyn'")
    return(logkda.polyn(a))
  }
  pari_string <- "\n allocatemem(); \n allocatemem();\n allocatemem(); \n allocatemem()
  logKDAvec(abund) =\n {\n
  local(S,J,n,m,k,k0,k1,k2,Told,Tnew,specabund,i,j,Sdiff,cnt1,cnt2);\n
  abund = vecsort(abund);\n S = length(abund);\n J = sum(k =
  1,S,abund[k]);\n maxabund = abund[S];\n Sdiff = 1; for(i =
  2,S,if(abund[i] != abund[i - 1],Sdiff++));\n specabund =
  matrix(2,Sdiff,i,j,0); specabund[1,1] = abund[1]; specabund[2,1] = 1;\n
  cnt1 = 1; cnt2 = 1;\n for(i = 2,S,\n    if(abund[i] != abund[i - 1],\n
  cnt1++;\n        cnt2 = 1;\n        specabund[1,cnt1] =
  abund[i];\n        specabund[2,cnt1] = cnt2\n    ,\n        cnt2++;\n
  specabund[2,cnt1] = cnt2\n    )\n );\n\n polyn = vector(1,i,1);\n i
  = 1;\n if(specabund[1,i] == 1,i++);\n Told = vector(1,i,1);\n for(n =
  2,maxabund,\n    Tnew = vector(n,m,(n > m) * Told[min(n-1,m)] +
  Told[max(1,m - 1)] * (m - 1)/(n - 1) + 0. );\n    if(n ==
  specabund[1,i],\n        for(k0 = 1,specabund[2,i],\n
  lenpolyn2 = length(polyn) + length(Tnew) - 1;\n           polyn =
  vector(lenpolyn2,k1,sum(k2 = max(1,k1 + 1 -
  length(Tnew)),min(length(polyn),k1),polyn[k2] * Tnew[k1 + 1 - k2]));\n\n
  );\n        i++;\n    );\n    Told = vector(n,m,Tnew[m]);\n );\n
  logKDA = log(polyn);\n\n logKDA\n }\n "
  if (isTRUE(as.logical(Sys.info()[1] == "Windows"))) {
    return(.logkda.pari.windows(a, numerical, pari_string))
  }
  else {
    return(.logkda.pari.unix(a, numerical, pari_string))
  }
}

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

NeutralParams <- function (data) 
{  
  S<- apply(data, 1, function(x) sum(x>0))
  J<- apply(data, 1, function(x) sum(x))
  H<- apply(data, 1, function(x) (x/sum(x))*(-log(x/sum(x)))) #Shannon Wiener Diversity precursor
  H <- colSums (H, na.rm=T) #Shannon Wiener Diversity Index 
  Neut <- cbind(S, J, H)
  Neutral <- data.frame(matrix(ncol=2, nrow=nrow(data)))
  colnames(Neutral) <- c(x="Theta", y="m")
  for(i in 1:nrow(data)) {
    d<- c("Theta", "m")
    Neutral <- data.frame(x=2,y=39)
    for(i in 1:nrow(data)) {
      
      l<- logkda.R(data[i,], use.brob=TRUE)  # Use logkda() if pari/gp is available
      z<- optimal.params(data[i,], log.kda=l)
      Neutral[i,] <- as.matrix(z)
    }
    colnames(Neutral) <- c(x="Theta", y="m")
    Neutral <- cbind(Neut, Neutral)
    return(Neutral)
  }}

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
#centiStation <- Stationtest(DIphyto)#Diphyto is a data frame made for depth integration from the loopdepthinteg.R script
#Build

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

