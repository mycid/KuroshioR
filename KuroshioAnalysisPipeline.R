%Kuroshio Analysis Pipline%
  
setwd("/home/trevor/Desktop/2009KuroshioData/KuroshioR")
folder <- getwd() #sets directory as an object/path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder
# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list[i], 
         read.csv(paste(folder, '/', file_list[i], sep=''))
  )} #Downloading all csv files to the Global Environment

print(file_list) #View current names
names <-list("Arctic","data", "DIp.or", "Al", "KC", "test", "origin", "Untitled3") #New names
#Renaming the csv files
print(file_list) #View current names
names <-list("Arctic","data", "DIp.or", "Al", "KC", "test", "origin", "Untitled3") #New names
for (i in 1:length(file_list)) {
  rename(file_list[i]<-names)}

#Preparing the main dataset
DiversityIndex <- Create.Diversity.I(test, 3,73) #Makes a new dataframe with diversity indexes 
#| Can use Original dataset which has all phytoplankton or test dataframe which is just those identified to species/genus

All <-Al[complete.cases(Al$Diatoms..cells.l.),] #selected diatoms..cells.l. because any rows reading NA would not match Kuroshio phytoplankton
DiversityIndex <- cbind(DiversityIndex, KuroAlldata.csv)#adding all the sample information to the newly made diversity Indices
