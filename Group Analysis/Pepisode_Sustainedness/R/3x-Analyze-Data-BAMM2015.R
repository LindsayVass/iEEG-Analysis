# Date:         19 August 2015
# Purpose:      This script will analyze the data from 2-Clean-Data.R. It will
#               perform the following analyses:
#               
# Compare the duration of oscillatory episodes for teleporter and navigation
# epochs. This analysis will identify episodes that cross a boundary (teleporter
# entry) and will test whether the mean duration of these episodes differs 
# between teleporter epochs and their matched navigation epochs. Specifically,
# it will look at how long the episode lasts after the boundary event rather
# than the full duration of the event. NB There is no true "boundary" event for
# the navigation epochs; we are simply using the same temporal event structure
# as the teleporter epochs.
#
# For the BAMM 2015 talk, this will focus on Delta/Theta only, teleporter entry
# only, and will also plot the means for navigation/teleportation.

library(permute)
library(reshape2)
library(dplyr)
# Functions ---------------------------------------------------------------

filterData <- function(dataFrame, electrode, freqband, timepoint) {
  dataFrame <- dataFrame %>%
    filter(ElectrodeID == electrode,
           FrequencyBand == freqband,
           TimePoint == timepoint)
}

testOscDuration <- function(navDf, teleDf, electrode, freqBand, timeType, boundaryType, varName, minTrials = 5) {
  
  navVarInd  <- which(names(navDf) == varName)
  teleVarInd <- which(names(teleDf) == varName)
  navData <- navDf %>%
    filter(ElectrodeID == electrode & FrequencyBand == freqBand)
  teleData <- teleDf %>%
    filter(ElectrodeID == electrode & FrequencyBand == freqBand)
  if (nrow(navData) < minTrials | nrow(teleData) < minTrials){
    output <- data.frame(ElectrodeID = electrode,
                         FrequencyBand = freqBand,
                         TimeType = timeType,
                         BoundaryType = boundaryType,
                         W = NA,
                         UncorrP = NA)
    inputData <- NA
  } else{
    theTest <- wilcox.test(navData[[navVarInd]], teleData[[teleVarInd]])
    output <- data.frame(ElectrodeID = electrode,
                         FrequencyBand = freqBand,
                         TimeType = timeType,
                         BoundaryType = boundaryType,
                         W = theTest$statistic,
                         UncorrP = theTest$p.value)
    
    navDataLab  <- cbind(navData, data.frame(Condition = "Navigation"))
    teleDataLab <- cbind(teleData, data.frame(Condition = "Teleporter"))
    inputData   <- rbind(navDataLab, teleDataLab)
  }
  
  return(list(trueResult = output, inputData = inputData))
  
}

meanPostEventOscDuration <- function(dataFrame, eventTime, trialTimeType) {
  output <- dataFrame %>%
    filter(Onset < eventTime & Offset > eventTime & TrialTimeType == trialTimeType) %>%
    mutate(PostEventDuration = Offset - eventTime) %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber) %>%
    summarise(MeanPostEventDuration = mean(PostEventDuration))
}


# Run analysis ------------------------------------------------------------

navNtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Navigation")
navFtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Navigation")

teleNtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Teleportation")
teleFtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Teleportation")

allData <- rbind(navNtPostEntryDur, navFtPostEntryDur, teleNtPostEntryDur, teleFtPostEntryDur) %>%
  filter(FrequencyBand == "Delta-Theta") %>%
  group_by(ElectrodeID, TimeType, Condition) %>%
  filter(n() > 5) %>%
  summarise(MeanDuration = mean(MeanPostEventDuration),
            SEM = sd(MeanPostEventDuration) / sqrt(n()))

navData <- allData %>%
  filter(Condition == "Navigation")
teleData <- allData %>%
  filter(Condition == "Teleportation")

validData <- inner_join(navData, teleData, by = c('ElectrodeID', 'TimeType')) %>%
  mutate(Difference = (MeanDuration.y - MeanDuration.x)) %>%
  select(ElectrodeID, TimeType, Difference) %>%
  inner_join(allData)

save(file = 'Rda/BAMM2015.Rda', list = 'validData')