# Script Name:  2-Clean-Data.R
# Author:       Lindsay Vass
# Date:         23 October 2015
# Purpose:      This script will clean the data loaded by 1-Load-Data.R

library(reshape2)

load('Rda/allRawData.Rda')

blankPepisode <- blankData %>%
  select(ElectrodeID, Condition:Pepisode) %>%
  rename(TrialTimeType = EpochLength) %>%
  filter(Frequency < 8)
blankPepisode$Condition <- plyr::mapvalues(blankPepisode$Condition, from = c('FreeExplore', 'Navigation'), to = c('BlankFree', 'BlankNav'))

navTelePepisode <- allPepisode %>%
  ungroup() %>%
  select(ElectrodeID, TrialTimeType, TimePoint, Frequency, TelePepisode, NavPepisode) %>%
  filter(Frequency < 8) 
navTelePepisode <- melt(navTelePepisode, id.vars = c('ElectrodeID', 'TrialTimeType', 'TimePoint', 'Frequency'), measure.vars = c('TelePepisode', 'NavPepisode'), variable.name = 'Condition', value.name = 'Pepisode')
navTelePepisode$TrialTimeType <- plyr::mapvalues(navTelePepisode$TrialTimeType, from = c('NT', 'FT'), to = c('Short', 'Long'))
navTelePepisode$Condition <- plyr::mapvalues(navTelePepisode$Condition, from = c('NavPepisode', 'TelePepisode'), to = c('Navigation', 'Teleportation'))
navTelePepisode <- navTelePepisode %>%
  mutate(Condition = paste(Condition, TimePoint, sep = "_")) %>%
  select(ElectrodeID, Condition, TrialTimeType, Frequency, Pepisode)

allPepisode <- rbind(blankPepisode, navTelePepisode) %>%
  group_by(ElectrodeID, Condition, TrialTimeType, Frequency)

save(file = 'Rda/allCleanData_sepBlank.Rda', list = 'allPepisode')