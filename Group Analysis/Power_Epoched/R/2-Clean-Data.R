# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         3 June 2015
# Purpose:      This script will clean the data frame produced by "1-Load-Data.R"

library(dplyr)

cleanData <- allData

# Put our conditions in order
timePointOrder <- c('Pre3', 'Pre2', 'Pre1', 'Tele', 'Post1', 'Post2', 'Post3')
cleanData$TimePoint <- factor(cleanData$TimePoint, levels = timePointOrder)

trialTypeOrder <- c('NSNT','NSFT','FSNT','FSFT')
cleanData$TrialType <- factor(cleanData$TrialType, levels = trialTypeOrder)

spaceOrder <- c('NS','FS')
cleanData$TrialSpaceType <- factor(cleanData$TrialSpaceType, levels = spaceOrder)

timeOrder <- c('NT','FT')
cleanData$TrialTimeType <- factor(cleanData$TrialTimeType, levels = timeOrder)

# Make a new variable for electrode ID (subjectID + Teleporter + Electrode)
cleanData <- cleanData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
  mutate(ElectrodeID = factor(ElectrodeID))

# Cut up the frequencies into bands
frequencies    <- unique(cleanData$Frequency)
freqBandBreaks <- c(0, 4, 8, 12, 30, 182)
freqBandNames  <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")

cleanData <- cleanData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames)) %>%
  select(ElectrodeID, TrialNumber, TrialSpaceType, TrialTimeType, TimePoint, FrequencyBand, Power)