# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         19 June 2015
# Purpose:      This script will clean the data loaded by 1-Load-Data.R

library(plyr)
library(dplyr)
library(tidyr)

load('Rda/allRawData.Rda')

# Clean up pepisodeData ------------------------------------------------------

pepisodeData <- pepisodeData %>%
  select(-ObservationID) %>%
  separate(ElectrodeID, c('Subject', 'Session', 'Electrode'), remove = FALSE) %>%
  select(-Electrode)
pepisodeData$Subject <- factor(pepisodeData$Subject)
pepisodeData$Session <- factor(pepisodeData$Session)

# initialize output
cleanPepisodeData <- data.frame()

# fix trial numbers
for (thisElectrode in 1:nlevels(pepisodeData$ElectrodeID)) {
  
  thisData <- pepisodeData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode])
  wrongTrialNumber <- sort(unique(thisData$TrialNumber))
  rightTrialNumber <- goodEpochs %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode]) %>%
    select(RealTrialNumber) %>%
    unlist()
  thisData <- thisData %>%
    mutate(TrialNumber = mapvalues(TrialNumber, wrongTrialNumber, rightTrialNumber))
  
  cleanPepisodeData <- rbind(cleanPepisodeData, thisData)
  
}

# Clean up latencyData ----------------------------------------------------

latencyData$Subject <- factor(latencyData$Subject)
latencyData$Session <- factor(latencyData$Session)


# Join the data frames ----------------------------------------------------

cleanData <- inner_join(cleanPepisodeData, latencyData) 

save(file = 'Rda/allCleanData.Rda', list = 'cleanData')