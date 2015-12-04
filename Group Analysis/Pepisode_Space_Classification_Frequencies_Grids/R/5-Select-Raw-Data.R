library(dplyr)
library(R.matlab)

load('Rda/allCleanData.Rda')
load('Rda/allClassificationByRegion.Rda')

parietal <- c('LPA3', 'LPA4', 'LPA5', 'LPA6', 'LPA7')
occipital <- c('LPA8', 'LMO6', 'LMO7', 'LMO8', 'LPO5', 'LPO6', 'LPO7', 'LPO8')

goodOccipElec <- occipitalResults %>%
  filter(Model == "NT",
         CorrP < 0.05)
goodPariElec <- parietalResults %>%
  filter(Model == 'NT',
         CorrP < 0.05)

occipData <- cleanData %>%
  inner_join(goodOccipElec) %>%
  filter(TimeBin == "Tele",
         TrialTimeType == 'NT',
         Frequency > 3.98,
         Frequency < 4) %>%
  select(ElectrodeID, TrialNumber, TrialSpaceType, TrialTimeType, Frequency, Pepisode)

parieData <- cleanData %>%
  inner_join(goodPariElec) %>%
  filter(TimeBin == 'Tele',
         TrialTimeType == 'NT',
         Frequency > 3.98,
         Frequency < 4) %>%
  select(ElectrodeID, TrialNumber, TrialSpaceType, TrialTimeType, Frequency, Pepisode)

writeMat('mat/pepisodeTrials.mat', occipData = occipData, parieData = parieData)
save(file = 'Rda/pepisodeTrials.Rda', list = c('occipData', 'parieData'))
