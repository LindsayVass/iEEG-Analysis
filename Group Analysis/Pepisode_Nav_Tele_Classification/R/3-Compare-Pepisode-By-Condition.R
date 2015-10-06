# Script Name:  3-Compare-Pepisode-By-Condition.R
# Author:       Lindsay Vass
# Date:         2 October 2015
# Purpose:      Since low-frequency episode during teleportation could 
#               distinguish navigation from teleportation, this script will test
#               whether pepisode is overall higher for one condition or the other
#               at the group level.


library(dplyr)
library(coin)
library(reshape2)

load('Rda/allClassificationResults.Rda')
load('../Pepisode_Sustainedness/Rda/allCleanData.Rda')

# Results from all electrodes ---------------------------------------------

teleData <- allPepisode %>%
  filter(TimePoint == "Tele") %>%
  select(ElectrodeID, RealTrialNumber, Frequency, TelePepisode, NavPepisode) %>%
  melt(id.vars = c("ElectrodeID", "RealTrialNumber", "Frequency"), measure.vars = c("TelePepisode", "NavPepisode"), variable.name = "Condition", value.name = "Pepisode")

meanPepisode <- teleData %>%
  ungroup() %>%
  group_by(ElectrodeID, Condition) %>%
  summarise(MeanPepisode = mean(Pepisode))

groupMeanPepisode <- meanPepisode %>%
  ungroup() %>%
  group_by(Condition) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            SEMPepisode = sd(MeanPepisode) / sqrt(n()))

groupWilcox <- wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = meanPepisode, distribution = approximate(B = 10000))


# Results from only signify classification electrodes ---------------------

goodElectrodes <- allMeanClassificationResults %>%
  filter(CorrP < 0.05) %>%
  select(ElectrodeID) %>%
  unique()

goodTeleData <- inner_join(teleData, goodElectrodes)

goodMeanPepisode <- goodTeleData %>%
  ungroup() %>%
  group_by(ElectrodeID, Condition) %>%
  summarise(MeanPepisode = mean(Pepisode))
goodMeanPepisode$ElectrodeID <- factor(goodMeanPepisode$ElectrodeID)

goodGroupMeanPepisode <- goodMeanPepisode %>%
  ungroup() %>%
  group_by(Condition) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            SEMPepisode = sd(MeanPepisode) / sqrt(n()))

goodGroupWilcox <- wilcoxsign_test(MeanPepisode ~ Condition | ElectrodeID, data = goodMeanPepisode, distribution = approximate(B = 10000))
