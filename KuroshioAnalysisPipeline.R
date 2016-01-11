#Kuroshio Analysis Pipline#
folder <- getwd() #sets directory as an object/path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list[i], 
         read.csv(paste(folder, '/', file_list[i], sep=''))
  )} #Downloading all csv files to the Global Environment

#Loading packages
packages <- c("akima","iNEXT", "rareNMtests", "RColorBrewer", "grid", "gplots", "geosphere", "ggplot2", "gridExtra", "plyr", "prabclus", "sparcl", "vegan")
library("akima", lib.loc="~/R/win-library/3.2")
library("colorspace", lib.loc="~/R/win-library/3.2")
library("geosphere", lib.loc="~/R/win-library/3.2")
library("ggplot2", lib.loc="~/R/win-library/3.2")
library("gplots", lib.loc="~/R/win-library/3.2")
library("gridExtra", lib.loc="~/R/win-library/3.2")
library("iNEXT", lib.loc="~/R/win-library/3.2")
library("plyr", lib.loc="~/R/win-library/3.2")
library("prabclus", lib.loc="~/R/win-library/3.2")
library("RColorBrewer", lib.loc="~/R/win-library/3.2")
library("sparcl", lib.loc="~/R/win-library/3.2")
library("vegan", lib.loc="~/R/win-library/3.2")
library("untb", lib.loc="~/R/win-library/3.2")
library("GUILDS", lib.loc="~/R/win-library/3.2")
#Preparing the main dataset
test <- Kuroshio_Phytoplankton.csv #shortens name
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
#Ploting relative to the front
#First adding distance from the front information to Adiv.abiotic
#Surface DistFrom Front imposed over whole data set!

A<- Adiv.abiotic[which(Adiv.abiotic[,"lon"]==143.5),]
Adist<- cbind(A$lon, A$lat)
dist1<- distm(Adist)
distance <- dist1[1,]*-.001
ADfromF<-distance+89.70125
B <- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==144.0),]
Bdist<- cbind(B$lon, B$lat)
bdist2<- distm(Bdist)
distance1 <- bdist2[1,]*.001
BDfromF <-distance1
C<- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==144.5),]
Cdist<- cbind(C$lon, C$lat)
dist3 <- distm(Cdist)
distance2 <- dist3[1,]*-.001
CDfromF <-distance2+79.3708
D <-Adiv.abiotic[which(Adiv.abiotic[, "lon"]==145.0),]
Ddist<- cbind(D$lon, D$lat)
dist4 <- distm(Ddist)
distance3 <- dist4[1,]*.001
DDfromF <-distance3-23.065398
E<- Adiv.abiotic[which(Adiv.abiotic[, "lon"]==145.5),]
Edist<- cbind(E$lon, E$lat)
dist5 <- distm(Edist)
distance4 <- dist5[1,]*-.001
EDfromF <-distance4
DfromF <- c(ADfromF, BDfromF, CDfromF, DDfromF,EDfromF)
Adiv.abiotic<-cbind(Adiv.abiotic, DfromF)
#Time to plot
NPfactorInColor(Adiv.abiotic, xvar="DfromF", yvar= "Cellcount", Adiv.abiotic$Gccom, xlab="DfromF", ylab="Cellcount", title="All") 

##########################################################
#Preparing a data set for rarefaction
#----------------------------------------------------------------------------------
#Rarefaction scripts 
Kuroshio_Phytoplankton <- Kuroshio_Phytoplankton.csv[-c(186:190), -c(1,2)]
Kuroshio_Genus <- test[-c(186:190), -c(1,2)]
Kuroshio_Phytoplankton <- cbind(clades, Kuroshio_Phytoplankton)
Kuroshio_Genus <- cbind(clades, Kuroshio_Genus)
#Creating a list for rarefaction based on the clusters
#AbundInf <- as.numeric(apply(test[1:185, 3:ncol(test)], 2, function(x) sum(x>0)))
Arc2 <- Arc[,-1] #Removes one factor from the arc dataset so that it can be clustered. 
#Can use either Arc2 or Kuroshio_Phytoplankton
SplitData(Arc2, 1, "dfcluster") #Splits Adiv.abiotic into multiple data frames based on the clusters 
SplitData(Kuroshio_Phytoplankton, 2, 'dfcluster')
libr <- list(dfcluster1, dfcluster2,dfcluster3) #,dfcluster4,dfcluster5) #Combines into one big list for the purpose of using a function on it to do everything at once

#turns the seperate datasets into incidence data then removes zeros, finally adding them together. There should be a function for this. 
#for Kuroshio
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
#--------------------------------------
z <- iNEXT(Rarefaction, q=0,  datatype="incidence", se=TRUE) #z is object class iNEXT. It provides all the data shown in the rarefaction graphs
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
#
IncidenceMatrix <- heatmap(t(Incidence), Rowv=NA, Colv=NA, scale="row",margins=c(10,5))
#
#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(" blue", "green", "red", "purple", "magenta"))(n = 499)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0, 1, length=100), 
               seq(1, 5, length=100),
               seq(5,10,length=100),              
               seq(10,15,length=100),
               seq(15,40,length=100))     

g <- heatmap.2(Incidence,
          cellnote = Incidence, # same data set for cell labels
          main = "Abundance",   # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",      # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=jet.colors(100),       # use on color palette defined earlier 
          #breaks=,    # enable color transition at specified limits
          dendrogram="none",   # only draw a row dendrogram
          Colv=NA,              # turn off column clustering
          key=TRUE,             # Turns on the key          # Key size adjustment
          )
###################################
#MEAN CHLOROPHYLL AND CELL COUNTS PER CLUSTER
###################################
SplitData(Adiv.abiotic, 2, "Adiv.cluster") #Splits Adiv.abiotic into multiple data frames based on the clusters 
Adiv.abiotic$Gccom<- as.factor(Adiv.abiotic$Gccom)
Chlorotest <- aov(Chlorophyll~Gccom, data=Adiv.abiotic)
summary(Chlorotest)
plot(Chlorotest)
kruskal.test(Chlorophyll~Gccom, data=Adiv.abiotic)
MeanChlor <- c(mean(Adiv.cluster1$Chlorophyll), mean(Adiv.cluster2$Chlorophyll), mean(Adiv.cluster3$Chlorophyll, na.rm=TRUE))
#
Celltest <- aov(Cellcount~Gccom, data=Adiv.abiotic)
summary(Celltest)
plot(Celltest)
kruskal.test(Richness~S, data=Adiv.abioticC)
MeanCell <- c(mean(Adiv.cluster1$Cellcount), mean(Adiv.cluster2$Cellcount), mean(Adiv.cluster3$Cellcount, na.rm=TRUE))
SDCell <- c(sd(Adiv.cluster1$Cellcount), sd(Adiv.cluster2$Cellcount), sd(Adiv.cluster3$Cellcount, na.rm=TRUE))
summary(aov(Richness~, data=Adiv.abiotic))
#
#Number of Species which appear once or twice for each cluster
Kurolets <- sum(clus1[which(clus1==2 | clus1==1)])
Deeplets <- sum(clus2[which(clus2==2 | clus2==1)])
Oyalets <- sum(clus3[which(clus3==2 | clus3==1)])

####################################################
Euclid <- bigtest[-c(186:190), c(3:length(bigtest))]
#Adiv.abioticC <- Adiv.abiotic[-c(which(rowSums(Euclid!=0)==FALSE))]
#Euclid <- Euclid[rowSums(Euclid)!=0, ] 
Euclid[Euclid>0]<- 1 
Euclid <- cbind(Adiv.abiotic$lon, Euclid)
SplitData(Euclid, 1, "EuclidLon")
SplitData(Adiv.abiotic, 15, "Transect")
dA <- vegdist(EuclidLon143.5, method="jaccard")
dB <- vegdist(EuclidLon144, method="jaccard", na.rm=TRUE)
dC <- vegdist(EuclidLon144.5, method="jaccard", na.rm=TRUE)
dD<- vegdist(EuclidLon145, method="jaccard", na.rm=TRUE)
dE <- vegdist(EuclidLon145.5, method="jaccard", na.rm=TRUE)

#Euclid[Euclid>0]<- 1 
#d <- vegdist(Euclid, na.rm=TRUE) # distance matrix

wardA <- hclust(dA, method="ward.D")
wardB <- hclust(dB, method="ward.D")
wardC <- hclust(dC, method="ward.D")
wardD <- hclust(dD, method="ward.D")
wardE <- hclust(dE, method="ward.D")

WardA <<- cutree(wardA, 3) #Creates three clades
WardB <<- cutree(wardB, 3)
WardC <<- cutree(wardC, 3)
WardD <<- cutree(wardD, 3)
WardE <<- cutree(wardE, 3)

TransectA <- cbind(Transect143.5, WardA)
TransectB <- cbind(Transect144, WardB)
TransectC <- cbind(Transect144.5, WardC)
TransectD <- cbind(Transect145, WardD)
TransectE <- cbind(Transect145.5, WardE)

ColorDendrogram(wardE, y = WardE, labels = names(WardE), main = "Transect E (with Abund)", branchlength = 1.5)
NPfactorInColor(TransectE, xvar="DfromF", yvar= "depth..m.", WardE, xlab="DfromF", ylab="depth..m.", title="Transect E (with Abund)") 
NPfactorInColor(TransectE, xvar="S", yvar= "Theta", WardE,xlab="Sal", ylab="Theta", title="Transect E (with Abund)")
###################################################
#Comparing the distance matrices
dA <- as.matrix(dA)
Alist <- as.list(dA)
Adist <- as.list(dist1)
plot(Alist, Adist)
#
dB <- as.matrix(dB)
Blist <- as.list(db)
Bdist <- as.list(dist2)
#
dC <- as.matrix(dC)
Clist <- as.list(dC)
Cdist <- as.list(dist3)
#
dD <- as.matrix(dD)
Dlist <- as.list(dD)
Ddist <- as.list(dist4)
#
dE <- as.matrix(dE)
Elist <- as.list(dE)
Edist <- as.list(dist5)
#
#All Distance Matrix comparison
Euclid <- test[-c(186:190), c(3:length(test))]
Adiv.abioticC <- Adiv.abiotic[-c(which(rowSums(Euclid!=0)==FALSE)),]
Alldist <- cbind(Adiv.abioticC$lon, Adiv.abioticC$lat)
distAll <- distm(Alldist)
distAll <- as.matrix(distAll)
Alllist <- as.list(distAll)
#This plots the total community dissimilarity and 
Euclid <- Euclid[rowSums(Euclid)!=0, ] 
Euclid[Euclid>0]<- 1 
BioticDist <- vegdist(Euclid, method="jaccard", na.rm=TRUE)
Allbio <- hclust(BioticDist, method="ward.D")
Allclust <- cutree(Allbio, 5)
ColorDendrogram(Allbio, y = Allclust, labels = names(Allclust), main = "Transect C (with Abund)", branchlength = 1.5)
BioticDist <- as.matrix(BioticDist)
AllBiolist <- as.list(BioticDist)
mantel(AllBioDist, distAll, method="pearson", permutations=999)
bandplot(Alllist, AllBiolist)
#trying to figure out rank abundance stuff below. Not working yet. 
rankabuncomp(x,y="",factor,scale="abundance",scaledx=F,type="o",rainbow=T,
             legend=T,xlim=c(1,max1), ylim=c(0,max2), ...)
rankabuncomp(Kuroshio_Phytoplankton[4:74], scaledx="abundance", rainbow = T, legend=T)
#-----------------------------------------------------
#below is a script for creating a species correlation matrx in heat map style
#It covers everything from organizing the dataset to graphing it and adding numbers
bigmatrix <- as.matrix(round(cor(test[,3:length(test)], method="pearson"), 3))
cormat <- reorder_cormat(bigmatrix)
bigmatrix <- get_upper_tri(cormat)
upper_tri <- as.data.frame(bigmatrix)
melted_cormat <- melt(bigmatrix)
melted_cormat <- na.omit(melted_cormat)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)
#
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
#-----------------------------------------------------------

DiaOver <- (bigtest$Diatoms/bigtest$Dinoflagellates)
DiaOver <- DiaOver[-c(186:190)]
Adiv.abiotic<- cbind(Adiv.abiotic, DiaOver="DiaOverDino")
EU <- ggplot(Adiv.abiotic[which("depth..m."==0),], aes(x =lat, y =depth..m., colour=DiaOverDino, label=round(DiaOverDino,2)))+ geom_point(size=7, alpha=.9, label=c, position=position_dodge(), stat="identity", )+scale_color_gradientn(colours=jet.colors(7), space="rgb", guide="colourbar")
EU <- EU+geom_text(aes(label=round(DiaOverDino,2)),hjust=0, vjust=-.5)
EU <- EU+labs(x="Salinity", y="Theta")+stat_contour(z=Adiv.abiotic$sigPoDen, binwidth = 2)
EU <- EU+theme(axis.title.x = element_text(color="cadet blue", vjust=-0.35, size=20, face="bold"), axis.title.y = element_text(color="cadetblue" , vjust=0.35, size=20, face="bold"))
CF <- EU+theme(axis.text.x=element_text(size=20, vjust=0.5), axis.text.y=element_text(size=20, vjust=.05))
CF #Temporary
Adiv.abiotic["depth..m."==0]
##---------------------------------------------------------
species.count(t(test[,c(3:length(test))]))
plot(species.count(t(test[,c(3:length(test))])),type="b")
matplot(species.table(t(test[,c(3:length(test))])),type="l",lty=1)
optimal.params.sloss(t(test[,c(3:73)]))
optimal.params.sloss(t(test[,c(3:73)]), nbres = 10, ci = TRUE, cint = c(0.025, 0.975))
#Examples
a <- untb(start=rep(1,50), prob=0.01, gens=2000, keep=TRUE)
plot(species.count(a),type="c")
matplot(species.table(a),type="l",lty=1)
##
PrestonR <- count(colSums(test[,c(3:length(test))]))
Preston<- preston(PrestonR, 16)*.1
summary.count(PrestonR)
optimal.prob(PrestonR)
plot(PrestonR)
plot.preston(Preston)
#Changing cell volumes to centi litters. A,B,C transects are .08, D and E are .1
centitest <- rbind(EuclidLon143.5[3:ncol(EuclidLon143.5)]/80, EuclidLon144[3:ncol(EuclidLon144)]/80, EuclidLon144.5[3:ncol(EuclidLon144.5)]/80, EuclidLon145[3:ncol(EuclidLon145)]/100, EuclidLon145.5[3:ncol(EuclidLon145.5)]/100)
centitest <- cbind(clades, centitest)
SplitData(centitest, 2, 'dfcluster')

#finding the optimal theta and m coefficients for all three water mass clusters
#Also setting up PrestonPlots for each
#Kuroshio
Kuroshio <- dfcluster1[,-c(1:3)]
Kuroshio <- Kuroshio[,which(colSums(Kuroshio)>0)]
KuroSum <- apply(Kuroshio, 2, function(x) sum(x))#sumrows
PresK <- count(KuroSum)
PrestonK<-preston(PresK)
plot.preston(PresK, main="Kuroshio")
plot.preston(PrestonK, main="Kuroshio")
#Oyashio
Oyashio <- dfcluster3[,-c(1:4)]
Oyashio <- Oyashio[,which(colSums(Oyashio)>0)]
OyaSum <- apply(Oyashio, 2, function(x) sum(x))#sumrows
PresO <- count(OyaSum)
PrestonO<-preston(PresO)
plot.preston(PresO, main="Oyashio")
plot.preston(PrestonO, main="Oyashio")
#Deep
Deep <- dfcluster2[,-c(1:4)]
Deep <- Deep[,which(colSums(Deep)>0)]
DeepSum <- apply(Deep, 2, function(x) sum(x))#sumrows
PresD <- count(DeepSum)
PrestonD<-preston(PresD)
plot.preston(PresD, main="Deep")
plot.preston(PrestonD, main="Deep")
plot()
#Forming a main dataset to estimate parameters for each one of these
####

optimal.theta(PresK, N=100)
optimal.params.sloss(t(Kuroshio), nbres=10)
optimal.params.gst(t(Kuroshio), exact = TRUE, ci = FALSE, cint = c(0.025, 0.975), nbres = 100)
#
l <- logkda.R(phytoInt[2,], use.brob=TRUE)  # Use logkda() if pari/gp is available

T <- optimal.params(phytoInt[1,], log.kda=l)  #compare his answer of 7.047958 and 0.22635923.
T1<- optimal.params(phytoInt[2,], log.kda=l)  #compare his answer of 7.047958 and 0.22635923.


#
KuroSim <- replicate(10, rand.neutral(no.of.ind(PresK),optimal.theta(PresK), string=("spp.")))
#
plot(PresK)
replicate(100, lines(rand.neutral(no.of.ind(PresK),optimal.theta(PresK), string=("spp.")), type="S", col="grey"))
#
Lines <- as.table(preston(rand.neutral(no.of.ind(PresK),optimal.theta(PresK), string=("spp."))))

par(mar=c(5.1,4.1,4.1,5.1)) #this must be before both plots are made
plot(PrestonK)
par(new=TRUE)
plot(Lines, axes=FALSE, xlab="", ylab="", type="b",lty=2, pch=24, col="black", bg=125)
axis(1, at=1:5, lab=c("Mon","Tue","Wed","Thu","Fri"))
mtext("Salinity",side=4,col="red",line=3) 
axis(4, ylim=c(33.6,34.3), col="blue" ,col.axis="blue",las=1)
legend("bottomleft",legend=c("S Gradient/Adjusted Distance","Salinity/Station Distance"), cex=.75,
       text.col=c("green","blue"),pch=c(23,24),col=c("green","blue"))


plot.count(PrestonR, uncertainty = TRUE, expectation = TRUE, theta = 6.810761, n = 100)
rand.neutral(no.of.ind(PrestonR), optimal.theta(PrestonR),string="NULL")
