# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         24 June 2015
# Purpose:      This script will load the csv files output by the Matlab
#               function "Run_Teleporter_Pepisode_Timing_Onset_Offset.m"


# Script ------------------------------------------------------------------

# data is stored in three csv files:
# observation_characteristics
# NT_data
# FT_data

charFiles <- dir('../Pepisode_Timing_Onset_Offset/csv/', pattern = "*characteristics.csv")
ntFiles <- dir('../Pepisode_Timing_Onset_Offset/csv/', pattern = "*NT_data.csv")
ftFiles <- dir('../Pepisode_Timing_Onset_Offset/csv/', pattern = "*FT_data.csv")

charData <- data.frame()
ntData   <- data.frame()
ftData   <- data.frame()

for (thisFile in 1:length(charFiles)) {
  
  tempData <- read.csv(paste0('../Pepisode_Timing_Onset_Offset/csv/', charFiles[thisFile]), header = TRUE)
  charData <- rbind(charData, tempData)
  
  # add information about subject, session, and electrode to data frame based on
  # file name
  fileNameParts <- strsplit(ntFiles[thisFile], split = "_")
  
  tempData <- read.csv(paste0('../Pepisode_Timing_Onset_Offset/csv/', ntFiles[thisFile]), header = TRUE)
  tempData$Subject   <- fileNameParts[[1]][1]
  tempData$Session   <- fileNameParts[[1]][2]
  tempData$Electrode <- fileNameParts[[1]][5]
  
  ntData   <- rbind(ntData, tempData)
  
  tempData <- read.csv(paste0('../Pepisode_Timing_Onset_Offset/csv/', ftFiles[thisFile]), header = TRUE)
  tempData$Subject   <- fileNameParts[[1]][1]
  tempData$Session   <- fileNameParts[[1]][2]
  tempData$Electrode <- fileNameParts[[1]][5]
  
  ftData   <- rbind(ftData, tempData)
  
}

remove(tempData)

# save concatenated data
dir.create('Rda')
save(file = 'Rda/allRawData.Rda', list = c('charData', 'ntData', 'ftData'))
