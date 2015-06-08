# Script Name:  3-Analyze-Data.R
# Author:       Lindsay Vass
# Date:         5 June 2015
# Purpose:      This script will analyze the clean data produced by 
#               "2-Clean-Data.R". It will first perform a within-electrode
#               analysis to identify when episodes occur during NT or FT trials
#               (NT and FT are analyzed separately since trial structure 
#               differs). It will then combine data across electrodes to look
#               at when episodes occur across the population.

library(dplyr)
library(reshape2)


# Within-electrode analysis -----------------------------------------------

episodeData <- data.frame()
onOffData <- data.frame()

for (thisElectrode in 1:nlevels(cleanCharData$ElectrodeID)) {
  
  for (thisFrequency in 1:nlevels(cleanCharData$FrequencyBand)) {
    
    # use joins to get the data for this electrode/frequency
    ntData <- cleanCharData %>%
      filter(ElectrodeID == levels(cleanCharData$ElectrodeID)[thisElectrode] &
             FrequencyBand == levels(cleanCharData$FrequencyBand)[thisFrequency]) %>%
      inner_join(cleanNtData) %>%
      group_by(ObservationType) %>%
      melt(c("ElectrodeID", 
             "ObservationID", 
             "TrialNumber", 
             "TrialTimeType", 
             "Frequency",
             "FrequencyBand", 
             "ObservationType"),
           variable.name = "Time",
           value.name = "Value")
    
    ftData <- cleanCharData %>%
      filter(ElectrodeID == levels(cleanCharData$ElectrodeID)[thisElectrode] &
             FrequencyBand == levels(cleanCharData$FrequencyBand)[thisFrequency]) %>%
      inner_join(cleanFtData) %>%
      group_by(ObservationType) %>%
      melt(c("ElectrodeID", 
             "ObservationID", 
             "TrialNumber", 
             "TrialTimeType", 
             "Frequency",
             "FrequencyBand", 
             "ObservationType"),
           variable.name = "Time",
           value.name = "Value")
    
    # R added X to positive timepoints and X. to negative timepoints, so fix that
    # now that they're no longer variable names
    thisData <- rbind(ntData, ftData)
    thisData$Time <- sub("^X\\.", "-", thisData$Time)
    thisData$Time <- sub("^X", "", thisData$Time)
    thisData$Time <- as.numeric(thisData$Time)
    
    # For episodes, get the mean and SEM for each time point
    tempEpisodeData <- thisData %>%
      filter(ObservationType == "Episode") %>%
      group_by(ElectrodeID, FrequencyBand, TrialTimeType, Time, ObservationType) %>%
      summarise(Mean = mean(Value), SEM = sd(Value) / sqrt(n()))
    
    episodeData <- rbind(episodeData, tempEpisodeData)
    
    # For onsets and offsets, keep only the timepoints with 1
    tempOnOffData <- thisData %>%
      filter(ObservationType != "Episode") %>%
      select(-ObservationID) %>%
      group_by(ElectrodeID, Frequency, TrialTimeType, Time, ObservationType) %>%
      filter(Value == 1) %>%
      arrange(TrialNumber)
    
    onOffData <- rbind(onOffData, tempOnOffData)
  }
  
}

# Save
save(file = 'Rda/allAnalyzedData.Rda', list = c("episodeData", "onOffData"))

