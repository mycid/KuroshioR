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

NPfactorInColor <- function(data, xvar, yvar, factor, xlab="", ylab=""){ #Plots the clade data for abiotic
  #data is the data frame, x is the x axis, y is the y axis, factor is the factor to be put in color
  EU <- ggplot(data, aes_string(x=xvar,y=yvar, colour=factor))+geom_point(size=7, alpha=.9, label=factor, position=position_dodge(), stat="identity")
  EU <- EU+labs(x=xlab, y=ylab)
  EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
  CF <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
  return(CF)
  }
  #Temporary
  EU <- ggplot(Adiv.abiotic, aes(x =S, y =Theta, colour=factor(Gccom), label=station))+ geom_point(size=7, alpha=.9, label=c, position=position_dodge(), stat="identity", )#+scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
  EU <- EU+geom_text(aes(label=station),hjust=0, vjust=-.5)
  EU <- EU+labs(x="Salinity", y="Theta")
  EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
  CF <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
  CF #Temporary

  #Function for splitting data by numbered factors
  SplitData <- function(data, column, header) { #data=dataframe, column=the the number of the column desired to be used as a numerical factor, header is the title of the data frames without the numbers at the end. MUST BE IN PARATHESES
  libr <- setNames(split(data, data[, column]), paste0(header, unique(data[,column]))) #Seperates the dataframe into multiple dataframes based on the factor
  list2env(libr, globalenv()) #Exports both s and each individual, new data frame to the global environment 
  }
  
 
  #Experimenting with functions for column division

x=Adiv.abiotic
orderlySpliting <- function(x, column, n, header) { 
#Function evenly (as possible)  
#divides data into a set number of groups
#Orders data into a numerical sequence
#Use just the column name for 'column', no quotes
#Returns new dataframes to Global 
#header choses a header name for the new dataframes
  arguments <- as.list(match.call()) # function piece comes from http://www.r-bloggers.com/passing-columns-of-a-dataframe-to-a-function-without-quotes/
  column = eval(arguments$column, x) 
  print(column)
  x <- x[,order(column),]
  Depths <- setNames(split(x,rep(1:n, ceiling(length(x)/n),length.out = length(x))), paste0(header, unique(data[,column])))
  list2env(Depths, globalenv())
}
