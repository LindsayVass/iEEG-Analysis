# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         16 June 2015
# Purpose:      This script will clean the data frames produced by "1-Load-Data.R"

library(dplyr)
load('Rda/allRawData.Rda')


# Functions ---------------------------------------------------------------

cutVariablesByTime <- function(varNames, cutPoints, binNames) {
  
  # R added X to positive timepoints and X. to negative timepoints so remove
  # them before converting to numeric
  origVarNames <- varNames
  varNames <- sub("^X\\.", "-", varNames)
  varNames <- sub("^X", "", varNames)
  varNames <- as.numeric(varNames)
  
  # Cut into time bins
  varNames <- cut(varNames, cutPoints, binNames)
  
  # Return a data frame that contains the original variable names and their
  # respective time bins
  df <- data.frame(Names = origVarNames, Bins = varNames, stringsAsFactors = FALSE)
  
}

getMeanPepisodeWithinBin <- function(varNames, binName, pepisodeData) {
  
  # get the column numbers corresponding to this binName
  thisBinNames <- varNames %>%
    filter(Bins == binName) %>%
    select(Names) 
  colNums <- match(thisBinNames$Names, names(pepisodeData))
  
  # take the mean across the columns indicated by colNums
  tempData <- pepisodeData %>%
    select(colNums) %>%
    rowMeans()
  
  # combine the mean pepisode data with the info in the original data frame
  output <- pepisodeData %>%
    select(ObservationID) %>%
    mutate(Pepisode = tempData,
           TimeBin = binName) 
  
}

# Clean charData ----------------------------------------------------------

# make new variable electrodeID to capture the fact that we treat the same
# electrode in different sessions as different observations
cleanCharData <- charData %>%
  select(ObservationID:TrialTimeType, Frequency) %>%
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
  select(-Frequency)

# Divide time points into pre/tele/post bins
timeBinBreaksNt <- c(-1830, 0, 1831, 3661)
timeBinBreaksFt <- c(-2830, 0, 2831, 5661)
timeBinNames  <- c("Pre", "Tele", "Post")

ntNames <- ntData %>%
  select(starts_with("X")) %>%
  names() %>%
  cutVariablesByTime(timeBinBreaksNt, timeBinNames)
ftNames <- ftData %>%
  select(starts_with("X")) %>%
  names() %>%
  cutVariablesByTime(timeBinBreaksFt, timeBinNames)

# Keep only episode data and make the ObservationID variable
tempNtData <- ntData %>%
  filter(ObservationType == "Episode") %>%
  mutate(ObservationID = paste(Subject, Session, Electrode, ObservationID, sep = "_")) 
tempFtData <- ftData %>%
  filter(ObservationType == "Episode") %>%
  mutate(ObservationID = paste(Subject, Session, Electrode, ObservationID, sep = "_"))

# Get the mean pepisode for each time bin
cleanNtData <- data.frame(ObservationID = NA, Pepisode = NA, TimeBin = NA)
for (thisBin in 1:nlevels(ntNames$Bins)) {
  
  tempData <- getMeanPepisodeWithinBin(ntNames, levels(ntNames$Bins)[thisBin], tempNtData)
  cleanNtData <- rbind(cleanNtData, tempData)
  
}

cleanFtData <- data.frame(ObservationID = NA, Pepisode = NA, TimeBin = NA)
for (thisBin in 1:nlevels(ftNames$Bins)) {
  
  tempData <- getMeanPepisodeWithinBin(ftNames, levels(ftNames$Bins)[thisBin], tempFtData)
  cleanFtData <- rbind(cleanFtData, tempData)
  
}

# Join the data together
cleanNtData <- inner_join(cleanCharData, cleanNtData)
cleanFtData <- inner_join(cleanCharData, cleanFtData)
cleanData <- rbind(cleanNtData, cleanFtData)

# Put our conditions in order
timeOrder  <- c('NT','FT')
spaceOrder <- c('NS', 'FS')
timePointOrder <- c('Pre', 'Tele', 'Post')
cleanData$TrialTimeType  <- factor(cleanData$TrialTimeType, levels = timeOrder)
cleanData$TrialSpaceType <- factor(cleanData$TrialSpaceType, levels = spaceOrder)
cleanData$TimeBin <- factor(cleanData$TimeBin, levels = timePointOrder)


# Save
save(file = 'Rda/allCleanData.Rda', list = 'cleanData')

