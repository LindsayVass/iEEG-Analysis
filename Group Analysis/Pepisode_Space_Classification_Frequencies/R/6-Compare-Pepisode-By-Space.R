# Script Name:  6-Compare-Pepisode-By-Space.R
# Author:       Lindsay Vass
# Date:         2 October 2015
# Purpose:      Since low-frequency episode during teleportation could 
#               distinguish short and long teleportation events, this script
#               will test whether pepisode is overall higher for a given
#               spatial condition, within each electrode.

library(dplyr)
library(coin)
library(ggplot2)

load('Rda/allCleanData.Rda')
load('Rda/allClassificationResults.Rda')


# Results from all electrodes ---------------------------------------------

teleData <- cleanData %>%
  filter(TimeBin == "Tele",
         Frequency < 8) %>%
  select(ElectrodeID, TrialNumber, TrialSpaceType, Frequency, Pepisode) %>%
  group_by(ElectrodeID, TrialSpaceType, TrialNumber)

meanPepisode <- teleData %>%
  ungroup() %>%
  group_by(ElectrodeID, TrialSpaceType) %>%
  summarise(MeanPepisode = mean(Pepisode))

groupMeanPepisode <- meanPepisode %>%
  ungroup() %>%
  group_by(TrialSpaceType) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
         SEMPepisode = sd(MeanPepisode) / sqrt(n()))

groupWilcox <- wilcoxsign_test(MeanPepisode ~ TrialSpaceType | ElectrodeID, data = meanPepisode, distribution = approximate(B = 10000))


# Results from only signify classification electrodes ---------------------

goodElectrodes <- allClassificationResults %>%
  filter(CorrP < 0.05) %>%
  select(ElectrodeID) %>%
  unique()

goodTeleData <- inner_join(teleData, goodElectrodes)

goodMeanPepisode <- goodTeleData %>%
  ungroup() %>%
  group_by(ElectrodeID, TrialSpaceType) %>%
  summarise(MeanPepisode = mean(Pepisode))
goodMeanPepisode$ElectrodeID <- factor(goodMeanPepisode$ElectrodeID)

goodGroupMeanPepisode <- goodMeanPepisode %>%
  ungroup() %>%
  group_by(TrialSpaceType) %>%
  summarise(GroupMeanPepisode = mean(MeanPepisode),
            SEMPepisode = sd(MeanPepisode) / sqrt(n()))

goodGroupWilcox <- wilcoxsign_test(MeanPepisode ~ TrialSpaceType | ElectrodeID, data = goodMeanPepisode, distribution = approximate(B = 10000))

# test whether distribution of frequencies differs
groupMeanFreqPepisode <- goodTeleData %>%
  group_by(Frequency, TrialSpaceType) %>%
  summarise(MeanPepisode = mean(Pepisode),
         SEMPepisode = sd(Pepisode) / sqrt(n()))
ksResult <- ks.test(groupMeanFreqPepisode$MeanPepisode[which(groupMeanFreqPepisode$TrialSpaceType == "NS")], 
                    groupMeanFreqPepisode$MeanPepisode[which(groupMeanFreqPepisode$TrialSpaceType == "FS")])
meanHist <- ggplot(groupMeanFreqPepisode, aes(x = Frequency, y = MeanPepisode, ymin = MeanPepisode - SEMPepisode, ymax = MeanPepisode + SEMPepisode, fill = TrialSpaceType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(position = "dodge")
meanHist