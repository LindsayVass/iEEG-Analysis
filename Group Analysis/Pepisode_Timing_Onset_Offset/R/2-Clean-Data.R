# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         5 June 2015
# Purpose:      This script will clean the data frames produced by "1-Load-Data.R"

library(dplyr)


# Clean charData ----------------------------------------------------------

# make new variable electrodeID to capture the fact that we treat the same
# electrode in different sessions as different observations
cleanCharData <- charData %>%
  select(-c(TrialSpaceType, TrialType:TimePoint)) %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
  select(-c(SubjectID, Teleporter, Electrode)) %>%
  mutate(ElectrodeID = factor(ElectrodeID))

# make ObservationID unique across subjects
cleanCharData <- cleanCharData %>%
  mutate(ObservationID = paste(ElectrodeID, ObservationID, sep = "_"))

# Cut up the frequencies into bands
frequencies    <- unique(cleanCharData$Frequency)
freqBandBreaks <- c(0, 4, 8, 12, 30, 182)
freqBandNames  <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
cleanCharData <- cleanCharData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames)) %>%
  select(ElectrodeID, ObservationID:TrialTimeType, Frequency, FrequencyBand)

# Put our conditions in order
timeOrder <- c('NT','FT')
cleanCharData$TrialTimeType <- factor(cleanCharData$TrialTimeType, levels = timeOrder)

# Update ObservationID for ntData and ftData
cleanNtData <- ntData %>%
  mutate(ObservationID = paste(Subject, Session , Electrode, ObservationID, sep = "_")) %>%
  select(-c(Subject, Session, Electrode))
cleanFtData <- ftData %>%
  mutate(ObservationID = paste(Subject, Session, Electrode, ObservationID, sep = "_")) %>%
  select(-c(Subject, Session, Electrode))

# Save
save(file = 'Rda/allCleanData.Rda', list = c('cleanCharData', 'cleanNtData', 'cleanFtData'))

