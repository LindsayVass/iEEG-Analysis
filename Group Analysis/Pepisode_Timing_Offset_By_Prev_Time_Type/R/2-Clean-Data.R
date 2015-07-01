# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         26 June 2015
# Purpose:      This script will clean the data loaded in by 1-Load-Data.R

library(dplyr)

# Load data ---------------------------------------------------------------

load('Rda/allRawData.Rda')


# Combine into one data frame ---------------------------------------------

offsetData <- data.frame()

for (thisElectrode in 1:nlevels(onOffData$ElectrodeID)) {
  
  thisOffsetData <- onOffData %>%
    ungroup() %>%
    filter(ElectrodeID == levels(onOffData$ElectrodeID)[thisElectrode] &
             ObservationType == "Offset")
  
  badTrialNumbers <- sort(unique(thisOffsetData$TrialNumber))
  
  realTrialNumbers <- allRealTrialNumbers %>%
    filter(ElectrodeID == levels(onOffData$ElectrodeID)[thisElectrode]) %>%
    mutate(BadTrialNumber = seq(1, length(badTrialNumbers)))
  
  thisTimeType <- allTimeTypes %>%
    filter(ElectrodeID == levels(onOffData$ElectrodeID)[thisElectrode]) 
  thisTimeType <- thisTimeType %>%
    mutate(NextTrialNumber = seq(1, nrow(thisTimeType)) + 1,
           PrevTimeType = ifelse(TimeType == 1, "NT", "FT"))
  if (1 %in% realTrialNumbers$RealTrialNumber) {
    thisTimeType <- rbind(data.frame(ElectrodeID = levels(onOffData$ElectrodeID)[thisElectrode],
                                     TimeType = NA,
                                     NextTrialNumber = 1,
                                     PrevTimeType = NA), thisTimeType)
  }
  thisTimeType <- thisTimeType %>% 
    select(-TimeType)
  
  thisOffsetData <- inner_join(thisOffsetData, realTrialNumbers, by = c("ElectrodeID", "TrialNumber" = "BadTrialNumber")) %>%
    select(ElectrodeID, RealTrialNumber, TrialTimeType, Frequency, Time) %>%
    inner_join(thisTimeType, by = c("ElectrodeID", "RealTrialNumber" = "NextTrialNumber"))
  
  offsetData <- rbind(offsetData, thisOffsetData)
}

offsetData$PrevTimeType <- factor(offsetData$PrevTimeType, levels = c("NT", "FT"))
save(file = 'Rda/allCleanData.Rda', list = 'offsetData')
