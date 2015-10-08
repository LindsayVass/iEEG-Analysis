# Script Name:  4-Select-Best-Trials.R
# Author:       Lindsay Vass
# Date:         8 October 2015
# Purpose:      Select 5 trials from favorite 2 electrodes. This will be output 
#               as a .mat file so the raw data can be extracted in Matlab.

library(R.matlab)
library(dplyr)

load('Rda/allClassificationResults.Rda')
load('../Pepisode_Sustainedness/Rda/allCleanData.Rda')

# number of trials
numTrials <- 5

topElectrodes <- data.frame(ElectrodeID = c("UCDMC13_TeleporterA_LAD1", "UCDMC14_TeleporterB_RHD2"))

# select top electrodes for the "both" model
topClass <- allMeanClassificationResults %>%
  filter(Model == "Both") %>%
  inner_join(topElectrodes)

# select 5 trials with nonzero pepisode from the top electrodes
bestData <- allPepisode %>%
  inner_join(topClass) %>%
  group_by(ElectrodeID, RealTrialNumber, FrequencyBand, TimePoint) %>%
  summarise(MeanTele = mean(TelePepisode), 
            MeanNav  = mean(NavPepisode)) 
bestData <- bestData %>%
  filter(TimePoint == "Tele",
         FrequencyBand == "Delta-Theta", 
         MeanTele != 0,
         MeanNav != 0) %>%
  mutate(MeanPepisode = (MeanTele + MeanNav)/2) %>%
  ungroup() %>%
  group_by(ElectrodeID) %>%
  arrange(desc(MeanPepisode)) %>%
  top_n(numTrials) %>%
  select(-c(MeanTele, MeanNav, MeanPepisode)) %>%
  inner_join(allPepisode) %>%
  select(ElectrodeID, RealTrialNumber, TrialSpaceType, TrialTimeType, Frequency, TelePepisode, NavPepisode) %>%
  filter(Frequency > 7)

# Convert from RealTrialNumber, to TrialNumber, which only indexes valid trials
bestTele <- bestData %>%
  select(ElectrodeID:RealTrialNumber) %>%
  inner_join(teleSustain) %>%
  select(ElectrodeID, RealTrialNumber, TrialNumber) %>%
  unique() %>%
  inner_join(bestData) %>%
  select(-NavPepisode) %>%
  rename(Pepisode = TelePepisode)

bestNav <- bestData %>%
  select(ElectrodeID:RealTrialNumber) %>%
  inner_join(navSustain) %>%
  select(ElectrodeID, RealTrialNumber, TrialNumber) %>%
  unique() %>%
  inner_join(bestData) %>%
  select(-TelePepisode) %>%
  rename(Pepisode = NavPepisode)

#dir.create('mat')
writeMat('mat/bestPepisodeTrials.mat', bestTele = bestTele, bestNav = bestNav)
save(file = 'Rda/bestPepisodeTrials.Rda', list = c('bestTele', 'bestNav', 'topClass'))
