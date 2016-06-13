#Messing with dataframes to get them to various shapes for analysis

#First working with dataframes to export for using in Tetame
#Will be making centitest, and centiall
#
#Adiv.abiotic: the phsical measurements for every sample
----------------------------------------------------------------
  DiversityIndex <- Create.Diversity.I(test, 3,length(test)) #Makes a new dataframe with diversity indexes 
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
  TotA <-data.frame(apply(bigtest[, 3:9], 1, function(x) sum(x))) #This is all cells for all microplankton
  Adiv.abiotic <- cbind(div.abiotic, All, TotA) #creates a new master dataset with additional information
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
  #Automated version
  Abiotic <- as.data.frame(cbind(Adiv.abiotic$S, Adiv.abiotic$Theta )) #Creates the desired dataframe with all the abiotic factors to be clustered
  AbioticCluster(Abiotic, 3) #Creates three different method dendrograms and three seperate clade lists based on each of those methods
  
#centitest
#Changing cell volumes to centi litters. A,B,C transects are .08, D and E are .1
----------------------------------------------------------------
#centitest for just diatoms and dinoflagellates
cent <- read.csv("test.csv")
cent <- cent[-c(186:190), c(4:length(cent))]
cent <- cbind(Adiv.abiotic$lon, cent)
SplitData(cent, 1, "centi")
centitest <- rbind(centi143.5[2:ncol(centi143.5)]/80, centi144[2:ncol(centi144)]/80, centi144.5[2:ncol(centi144.5)]/80, centi145[2:ncol(centi145)]/100, centi145.5[2:ncol(centi145.5)]/100)
#centiorigin <- centitest[-c(97:99), ]#need to remove rows 97:99 because for some reason these rows do not have number of individuals 
write.csv(centitest, "centitest.csv")
centiDiverse <- Create.Diversity.I(centitest, 1, ncol(centitest))
Adiv.abiotic <- cbind(centiDiverse, Adiv.abiotic)
#
---------------------------------------------------
#centitest for all organisms
phyto <-read.csv("OriginalPhytoplankton.csv")
Cent <- phyto[-c(186:190), c(3:81)]
Cent <- cbind(Adiv.abiotic$lon, Cent)
SplitData(Cent, 1, "centi")
centiorigin <- rbind(centi143.5[2:ncol(centi143.5)]/80, centi144[2:ncol(centi144)]/80, centi144.5[2:ncol(centi144.5)]/80, centi145[2:ncol(centi145)]/100, centi145.5[2:ncol(centi145.5)]/100)
centiorigin <- centitest[-c(97:99), ]#need to remove rows 97:99 because for some reason these rows do not have number of individuals 
write.csv(centiorigin, "centiorigin.csv")
centiODiverse <- Create.Diversity.I(centitest, 1, ncol(centitest))
COTRM <- read.csv("centiOriginTetResultsM.csv")
Adiv.abioCO <- cbind(centiDiverse, Adiv.abioticN[, -c(6:11)])
Adiv.abioCO <- cbind("m"=COTRM$m, "I"=COTRM$I, Adiv.abioticCenti)
#  
---------------------------------------------------
#Creating summed phytoplankton counts for each water mass (total of 3, Kuroshio, Oyashio and Deep)
Adiv.abioticN <- Adiv.abiotic[-c(97:99),]# for use later in analysis
centitest<- cbind(Adiv.abioticN$Gccom, centitest)
SplitData(centitest, 1, "centiclus")
#Kuroshio cluster
Kurorows <- centiclus1[2]
Kuroshio <- centiclus1[,-c(1:2), -80]
KuroSum <- apply(Kuroshio, 2, function(x) sum(x))#sumrows
#Oyashio cluster
Oyarows <- centiclus3[2]
Oyashio<- centiclus3[,-c(1:2), -80]
OyaSum <- apply(Oyashio, 2, function(x) sum(x))#sumrows
#Deep cluster
Deeprows <- centiclus2[2]
Deep<- centiclus2[,-c(1:2), -80]
DeepSum <- apply(Deep, 2, function(x) sum(x))#sumrows
#Putting it all together
phytoclusters<- rbind(KuroSum, OyaSum, DeepSum)
write.csv(phytoclusters, "phytoclusters.csv")
