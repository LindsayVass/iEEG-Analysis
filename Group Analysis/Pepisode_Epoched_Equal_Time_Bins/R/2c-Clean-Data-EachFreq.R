# Script Name:  2c-Clean-Data-EachFreq.R
# Author:       Lindsay Vass
# Date:         19 November 2015
# Purpose:      This script will clean the data frame produced by "1-Load-Data.R".
#               It differs from 2-Clean-Data in that each frequency is kept 
#               separate (i.e., not grouped into bands.)

library(dplyr)

load('Rda/allRawData.Rda')

cleanData <- allData

# Put our conditions in order
timePointOrder <- c('Pre1', 'Tele', 'Post1')
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

save(file = 'Rda/allCleanData_EachFreq.Rda', list = 'cleanData')