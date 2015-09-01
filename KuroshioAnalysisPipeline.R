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

#Preparing the main dataset
test <- Kuroshio_Phytoplankton.csv #shortens name
DiversityIndex <- Create.Diversity.I(test, 3,73) #Makes a new dataframe with diversity indexes 
#| Can use Original dataset which has all phytoplankton or test dataframe which is just those identified to species/genus

All <-KuroAlldata.csv[complete.cases(KuroAlldata.csv$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton. Also changes the name to Al.
DiversityIndex <- cbind(DiversityIndex, All)#adding all the sample information to the newly made diversity Indices

#Calculating Potential Density and Temperature
d <- All$depth..m. #depth from the data set
t <- All$T.C.        #Temperature from the data set
s <- All$S #Salinity from the data set
div.abiotic <- PotentialST(d,t,s ) #function for calculating the potential density and temperature 
div.abiotic$theta=NULL #gets rid of other theta in the data set

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
Abiotic <- cbind(Adiv.abiotic$S, Adiv.abiotic$Theta) #Creates the desired dataframe with all the abiotic factors to be clustered
AbioticCluster(Abiotic, 5) #Creates three different method dendrograms and three seperate clade lists based on each of those methods
#names of the clade data sets it saves are: 
#1. Gccom which is the "complete" method
#2. Gcaver which is the "aver" method
#3. Gcsin which is the "single" method
#Objects will be returned together as "clade"

Adiv.abiotic <- cbind(clades, Adiv.abiotic ) #adds the abiotic clades in a number vector to the main data set
#Time to graph:
NPfactorInColor(Adiv.abiotic, S, Theta, Gccom, xlab="Salinity", ylab="Theta") #Doesnt work just yet
#Problem here is that it needs a FUN function somewhere for reasons I don't quite understand

#Experimenting with rarefaction
Kuroshio_Phytoplankton <- Kuroshio_Phytoplankton.csv[,-c(1,2)] #Takes off bottle data
#Creates a rarified data set based on number of samples. Interpolated
KuroSampleR <- rarefaction.sample(Kuroshio_Phytoplankton, method = "sample-size", q = 0)
#Creates a rareified data set based on the Coverage. Interpolated
KuroCoverageR <-rarefaction.sample(Kuroshio_Phytoplankton, method = "coverage", q = 0)
#
#Plot the results
par(mfcol=c(1,2))
plot(KuroSampleR[,1], KuroSampleR[,2], lwd=2, xlab="Sampling units", ylab="Hill numbers",
     type="l", main="Sample-based rarefaction")
plot(KuroCoverageR[,1], KuroCoverageR[,2], lwd=2, xlab="Coverage", ylab="Hill numbers",
     type="l", main="Coverage-based rarefaction")
