# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         26 June 2015
# Purpose:      This script will load the analyzed data from Pepisode_Timing_Onset_Offset
#               as well as Matlab data containing the true trial numbers and the
#               associated conditions (necessary for determining the previous
#               trial's time type).

library(dplyr)
library(R.matlab)


# Functions ---------------------------------------------------------------

getFilePaths <- function(electrodeID) {
  
  infoParts <- strsplit(electrodeID, "_")
  goodEpochsPath <- paste0('/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/',
                           infoParts[[1]][1],
                           '/Mat Files/',
                           infoParts[[1]][1],
                           '_',
                           infoParts[[1]][2],
                           '_',
                           substr(infoParts[[1]][3], 1, 3),
                           '_noSpikes_noWaves_goodEpochs.mat')
  epochInfoPath <- paste0('/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/',
                          infoParts[[1]][1],
                          '/Mat Files/',
                          infoParts[[1]][1],
                          '_',
                          infoParts[[1]][2],
                          '_Epochs_Entry.mat')
  return(list(goodEpochsPath = goodEpochsPath,
              epochInfoPath = epochInfoPath))
}

# Load pepisode timing data -----------------------------------------------

load('../Pepisode_Timing_Onset_Offset/Rda/allAnalyzedData.Rda')
remove(episodeData)

# Load Matlab data --------------------------------------------------------
allRealTrialNumbers <- data.frame()
allTimeTypes <- data.frame()

for (thisElectrode in 1:nlevels(onOffData$ElectrodeID)) {
  
  thePaths <- getFilePaths(levels(onOffData$ElectrodeID)[thisElectrode])
  realTrialNumbers <- readMat(thePaths$goodEpochsPath)$goodEpochs
  timeTypes <- readMat(thePaths$epochInfoPath)$eTime
  allRealTrialNumbers <- rbind(allRealTrialNumbers,
                               data.frame(ElectrodeID = levels(onOffData$ElectrodeID)[thisElectrode],
                                          RealTrialNumber = realTrialNumbers))
  allTimeTypes <- rbind(allTimeTypes,
                        data.frame(ElectrodeID = levels(onOffData$ElectrodeID)[thisElectrode],
                                   TimeType = timeTypes))
  
}
dir.create('Rda')
save(file = 'Rda/allRawData.Rda', list = c("onOffData", "allRealTrialNumbers", "allTimeTypes"))
