# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         19 June 2015
# Purpose:      This script will load the csv files output by the Matlab
#               script "runCalculateBoundaryCrossingLatency.m", which contains
#               the latency to enter the correct arm of the maze for trials in
#               which the patient was not already facing the correct arm. It 
#               will also load the cleaned pepisode data from 
#               Pepisode_SpaceXTimepoint_ANOVA/Rda/allCleanData.Rda. Finally, it
#               loads the list of good epochs (trials) for each electrode/session.

library(R.matlab)
library(dplyr)
library(tidyr)

# Load latency csv files --------------------------------------------------

# latency data is stored in csv files
inputFiles  <- dir('csv/', pattern = "*.csv")
latencyData <- data.frame()

for (thisFile in 1:length(inputFiles)) {
  
  tempData <- read.csv(paste0('csv/', inputFiles[thisFile]), header = TRUE)
  
  # add information about subject, session, and electrode to data frame based on
  # file name
  fileNameParts <- strsplit(inputFiles[thisFile], split = "_")
  
  tempData$Subject   <- fileNameParts[[1]][1]
  tempData$Session   <- fileNameParts[[1]][2]
  latencyData        <- rbind(latencyData, tempData)
  
}

remove(tempData)


# Load clean pepisode data ------------------------------------------------

load('../Pepisode_SpaceXTimepoint_ANOVA/Rda/allCleanData.Rda')

pepisodeData <- cleanData


# Load good epochs data ---------------------------------------------------

# the pepisode data does not contain the true trial numbers -- to retrieve this
# data, we'll load in the goodEpochs data that was output by Matlab during
# pre-processing

# get subject/session/electrode data from the pepisode data
observationInfo <- data.frame(ElectrodeID = cleanData$ElectrodeID) %>%
  separate(ElectrodeID, c('Subject', 'Session', 'Electrode')) %>%
  unique()

# initialize output
goodEpochs <- data.frame()

for (thisObservation in 1:nrow(observationInfo)) {
  
  # load in the data for this electrode/session
  matPath <- paste0('/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/', 
                    observationInfo$Subject[thisObservation],
                    '/Mat Files/',
                    observationInfo$Subject[thisObservation],
                    '_',
                    observationInfo$Session[thisObservation],
                    '_',
                    gsub("[[:digit:]]*$", "", observationInfo$Electrode[thisObservation]),
                    '_noSpikes_noWaves_goodEpochs.mat')
  matData <- readMat(matPath)$goodEpochs
  
  # create data frame with good epoch info
  tempData <- data.frame(ElectrodeID = paste(observationInfo$Subject[thisObservation],
                                             observationInfo$Session[thisObservation],
                                             observationInfo$Electrode[thisObservation],
                                             sep = "_"),
                         RealTrialNumber = matData)
  
  # add it to the output
  goodEpochs <- rbind(goodEpochs, tempData)
}


# Save data ---------------------------------------------------------------

dir.create('Rda')
save(file = 'Rda/allRawData.Rda', list = c('latencyData', 'pepisodeData', 'goodEpochs'))
