# Script Name:  4b-Select-Single-Trial-Data.R
# Author:       Lindsay Vass
# Date:         3 September
# Purpose:      This script will identify trials that show similar pepisode for
#               navigation and teleportation, so that the raw data can be 
#               extracted in Matlab.


library(R.matlab)
library(reshape2)
library(dplyr)

load('Rda/allCleanData.Rda')
load('Rda/allAnalyzedData.Rda')
load('Rda/allRawData.Rda')
remove(sessionInfo, durationPermResults, durationTrueDataCorrected, durationTrueNSigElectrodesCorrected, permDurationData, navSustain, teleSustain, allNavigationPepisodeData, allTeleporterPepisodeData, allNavigationSustainData, allTeleporterSustainData)

sigElectrodes <- pepisodeTrueDataCorrected %>%
  filter(FrequencyBand == "Delta-Theta",
         TimePoint == "Post1",
         CorrP < 0.05) %>%
  select(-c(statistic, p.value, Iteration, CorrP, TimePoint))

sigElectrodeTrials <- allPepisode %>%
  inner_join(sigElectrodes) %>%
  ungroup() %>%
  arrange(ElectrodeID, RealTrialNumber, Frequency, TimePoint) %>%
  mutate(NavMinusTele = NavPepisode - TelePepisode) %>%
  dcast(ElectrodeID + RealTrialNumber + TrialSpaceType + TrialTimeType + TrialType + Frequency ~ TimePoint) %>%
  filter(Pre1 != 0,
         Tele != 0,
         Post1 != 1) %>%
  arrange(desc(Post1), abs(Pre1), abs(Tele)) # look for large difference at Post and small difference at Pre/Tele
  

freqs <- unique(allPepisode$Frequency)

navTrialNumbers <- realNavigationTrialNumbers %>%
  group_by(Subject, Session, Electrode) %>%
  arrange(Subject, Session, Electrode, Trial) %>%
  rename(RealTrialNumber = Trial) %>%
  mutate(NavTrialNumber = row_number())
teleTrialNumbers <- realTeleporterTrialNumbers %>%
  group_by(Subject, Session, Electrode) %>%
  arrange(Subject, Session, Electrode, Trial) %>%
  rename(RealTrialNumber = Trial) %>%
  mutate(TeleTrialNumber = row_number())


manualBestTrials <- data.frame(ElectrodeID = c("UCDMC15_TeleporterB_LHD3",
                                                 "UCDMC15_TeleporterB_RHD3",
                                                 "UCDMC15_TeleporterB_LAD5",
                                                 "UCDMC15_TeleporterB_RHD4",
                                                 "UCDMC15_TeleporterB_RHD1",
                                                 "UCDMC14_TeleporterA_RHD1",
                                                 "UCDMC15_TeleporterB_RHD4"),
                                 RealTrialNumber = c(57, 27, 13, 39, 12, 2, 60),
                                 Frequency = c(freqs[10],
                                               freqs[10],
                                               freqs[10],
                                               freqs[8],
                                               freqs[7],
                                               freqs[10],
                                               freqs[8])) %>%
  inner_join(allPepisode) %>%
  rowwise() %>%
  mutate(Subject = strsplit(ElectrodeID, '_')[[1]][1],
         Session = strsplit(ElectrodeID, '_')[[1]][2],
         Electrode = substr(strsplit(ElectrodeID, '_')[[1]][3], 1, 3)) %>%
  inner_join(navTrialNumbers) %>%
  inner_join(teleTrialNumbers) %>%
  select(-c(Subject, Session, Electrode))


save(file = 'Rda/bestSingleTrials.Rda', list = 'manualBestTrials')
writeMat('mat/bestSingleTrials.mat', manualBestTrials = manualBestTrials)

