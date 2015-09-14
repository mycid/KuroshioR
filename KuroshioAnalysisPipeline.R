%Kuroshio Analysis Pipline%
  
setwd("/home/trevor/Desktop/2009KuroshioData/KuroshioR")
folder <- getwd() #sets directory as an object/path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list[i], 
         read.csv(paste(folder, '/', file_list[i], sep=''))
  )} #Downloading all csv files to the Global Environment

#Loading packages
packages <- c("akima", "geosphere", "ggplot2", "gridExtra", "plyr", "prabclus", "sparcl", "vegan")
library("akima", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("geosphere", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("ggplot2", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("plyr", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("gridExtra", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("prabclus", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("sparcl", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("vegan", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("iNEXT", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("rareNMtests", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("RColorBrewer", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
library("grid", lib.loc="/home/trevor/R/i686-pc-linux-gnu-library/3.2")
#Preparing the main dataset
test <- Kuroshio_Phytoplankton.csv #shortens name
DiversityIndex <- Create.Diversity.I(test, 3,73) #Makes a new dataframe with diversity indexes 
#| Can use Original dataset which has all phytoplankton or test dataframe which is just those identified to species/genus

All <-KuroAlldata.csv[complete.cases(KuroAlldata.csv$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton. Also changes the name to Al.

#Calculating Potential Density and Temperature
d <- All$depth..m. #depth from the data set
t <- All$T.C.        #Temperature from the data set
s <- All$S #Salinity from the data set
div.abiotic <- PotentialST(d,t,s ) #function for calculating the potential density and temperature 

#Expanding the data set; more things to analyze
bigtest <- Allbigstuff.csv #changing the datset to an easier name
#|This data set has the generalized large classes of micro plankton found in the samples. It has a total cells per liter for all cells
TotA <-apply(bigtest[, 3:9], 1, function(x) sum(x)) #This is all cells for all microplankton
Transect <-Untitled3.csv #This file created to add the categorial transect data
Transect <- rbind("A", Transect) #A fix for Transect, adds a missing row
Adiv.abiotic <- cbind(Transect, TotA, div.abiotic) #creates a new master dataset with additional information
Adiv.abiotic <- div.abiotic[-c(186:190),] #removes the last 5 rows because they have missing data

#Preparing Chlorophyll data
KuroChl <- na.omit(Kuro_chl_coords.csv) #omits unknown Chlorophyll samples
calibrate <- lm(fluorescence~chlorophyll..ug.l., data=KuroChl) #Creates regression for calibration
summary(calibrate) #summarizes that regression
coeffs = coefficients(calibrate) #Retrieves coeffs
KuroChlNA <- subset(Kuro_chl_coords.csv, is.na(Kuro_chl_coords.csv$chlorophyll..ug.l.)) #subsets dataframe for NA values
Chlorophyll <- coeffs[1] + coeffs[2]*KuroChlNA$fluorescence #Creates estimated chlorophyll for the NA subset dataframe
Kuro_chl_coords.csv$chlorophyll..ug.l.[is.na(Kuro_chl_coords.csv$chlorophyll..ug.l.)] <-Chlorophyll[1:sum(is.na(Kuro_chl_coords.csv$chlorophyll..ug.l.))] 
#Above this takes the estimated Chlorophyll for all the NA's and writes them over the NA's in the original dataframe
Chlorophyll <- Kuro_chl_coords.csv$chlorophyll..ug.l.[c(1:185)] #shortens the Chlorophyll length to the main dataset, creates seperate vector
Fluorescence <- Kuro_chl_coords.csv$fluorescence[c(1:185)]
Adiv.abiotic <- cbind(Chlorophyll, Fluorescence, Adiv.abiotic) #adds Chlorophyll to the dataset
#
#Clustering data based on abiotic factors
Abiotic <- cbind(Adiv.abiotic$S, Adiv.abiotic$Theta ) #Creates the desired dataframe with all the abiotic factors to be clustered
AbioticCluster(Abiotic, 3) #Creates three different method dendrograms and three seperate clade lists based on each of those methods


#names of the cladedata sets it saves are: 
#1. Gccom which is the "complete" method
#2. Gcaver which is the "aver" method
#3. Gcsin which is the "single" method
#Objects will be returned together as "clade"

Adiv.abiotic <- cbind(clades, Adiv.abiotic ) #adds the abiotic clades in a number vector to the main data set
#Time to graph:
NPfactorInColor(Adiv.abiotic, xvar="S", yvar= "Theta", Gccom, xlab="Salinity", ylab="Theta") #Doesnt work just yet
#Problem here is that it needs a FUN function somewhere for reasons I don't quite understand
#First some color
colors <- colorRampPalette(c("black", "#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "purple", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "#7F0000"))
#Next divide up the data into five distinct depth groups of about equal size. 
orderlySpliting(Adiv.abiotic, "depth..m.", 5, "Depth")
#Manual code for subseting data by depth

Surface <- subset(Adiv.abiotic, depth..m. == 0)
Shallow <- subset(Adiv.abiotic, depth..m. > 1 & depth..m. <= 25, select=c(Gccom:Si.))
Medium <- subset(Adiv.abiotic, depth..m. > 25 & depth..m. <= 40, select=c(Gccom:Si.))
MeDeep <- subset(Adiv.abiotic, depth..m. > 40 & depth..m. <= 75 , select=c(Gccom:Si.))
Deep <- subset(Adiv.abiotic, depth..m. > 75)

#Ploting depth layers 
NPfactorInColor(Surface, xvar="lon", yvar= "lat", Surface$Gccom, xlab="Latitude", ylab="Depth", title="Surface")
NPfactorInColor(Shallow, xvar="lon", yvar= "lat", Shallow$Gccom, xlab="Latitude", ylab="Depth",title="1:25m")
NPfactorInColor(Medium, xvar="lon", yvar= "lat", Medium$Gccom, xlab="Latitude", ylab="Depth",title="26:40m")
NPfactorInColor(MeDeep, xvar="lon", yvar= "lat", MeDeep$Gccom, xlab="Latitude", ylab="Depth",title="41:75m")
NPfactorInColor(Deep, xvar="lon", yvar= "lat", Deep$Gccom, xlab="Latitude", ylab="Depth", title=">75m") 
NPfactorInColor(Adiv.abiotic, xvar="lat", yvar= "depth..m.", Adiv.abiotic$Gccom, xlab="Latitude", ylab="Depth", title="All") 
#Preparing a data set for rarefaction
#----------------------------------------------------------------------------------
#Rarefaction scripts 
Kuroshio_Phytoplankton <- Kuroshio_Phytoplankton.csv[-c(186:190), -c(1,2)]
Kuroshio_Phytoplankton <- cbind(clades, Kuroshio_Phytoplankton)
#Creating a list for rarefaction based on the clusters
AbundInf <- as.numeric(apply(test[1:190, 3:73], 2, function(x) sum(x>0)))
SplitData(Kuroshio_Phytoplankton, 2, "dfcluster") #Splits Adiv.abiotic into multiple data frames based on the clusters 
libr <- list(dfcluster1, dfcluster2,dfcluster3) #,dfcluster4,dfcluster5) #Combines into one big list for the purpose of using a function on it to do everything at once

#turns the seperate datasets into incidence data then removes zeros, finally adding them together. There should be a function for this. 
clus1 <- as.numeric(apply(dfcluster1[,4:ncol(dfcluster1)], 2, function(x) sum(x>0)))
clus1 <- c(nrow(dfcluster1), clus1[clus1>0])
clus2 <- as.numeric(apply(dfcluster2[,4:ncol(dfcluster2)], 2, function(x) sum(x>0)))
clus2 <- c(nrow(dfcluster2), clus2[clus2>0])
clus3 <- as.numeric(apply(dfcluster3[,4:ncol(dfcluster3)], 2, function(x) sum(x>0)))
clus3 <- c(nrow(dfcluster3), clus3[clus3>0])
#clus4 <- as.numeric(apply(dfcluster4[,4:ncol(dfcluster4)], 2, function(x) sum(x>0)))
#clus4 <- c(nrow(dfcluster4), clus4[clus4>0])
#clus5 <- as.numeric(apply(dfcluster5[,4:ncol(dfcluster5)], 2, function(x) sum(x>0)))
#clus5 <- c(nrow(dfcluster5), clus5[clus5>0])
Rarefaction <- list("1"=clus1, "2"=clus2, "3"=clus3)# "5"=clus5) 

z <- iNEXT(Rarefaction, q=0, datatype="incidence", se=TRUE) #z is object class iNEXT. It provides all the data shown in the rarefaction graphs
ggiNEXT(z, type = 1, se = TRUE, facet.var = "order", color.var = "order") #Creates rarefaction graphs
#--------------------------------------------------------------------------
#
clus1 <- as.numeric(apply(dfcluster1[,4:ncol(dfcluster1)], 2, function(x) sum(x>0)))
clus2 <- as.numeric(apply(dfcluster2[,4:ncol(dfcluster2)], 2, function(x) sum(x>0)))
clus3 <- as.numeric(apply(dfcluster3[,4:ncol(dfcluster3)], 2, function(x) sum(x>0)))
Incidence <- t(cbind(clus1, clus2, clus3))
Kuroshio_Phytoplankton <- Kuroshio_Phytoplankton.csv[-c(186:190), -c(1,2)]
colnames(Incidence)<-colnames(Kuroshio_Phytoplankton)
Incidence <- as.matrix(Incidence)
image((Incidence[nrow(Incidence):1,] ), axes=FALSE, zlim=c(0,25), useRaster = TRUE, col = rainbow(21))
  col = grey(seq(0, 1)))
image1 <- image1+axis( 2, at=seq(0,1,length.out=ncol(Incidence ) ), labels= colnames( Incidence), las= 2 )
image1 <- image1+axis( 1, at=seq(0,1,length.out=nrow(Incidence ) ), labels= ("go"), las= 2)
image1 <- image1+axis( 3, at=seq(0,1,length.out=nrow( Incidence ) ), labels= ("be"), las= 1)

rc <- rainbow(nrow(Incidence), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
IncidenceMatrix <- heatmap(t(Incidence), Rowv=NA, Colv=NA, scale="row",margins=c(10,5))
